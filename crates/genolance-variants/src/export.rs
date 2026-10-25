//! VCF export - reconstruct a VCF for a single sample (or region) from
//! the Lance store. Per-alt split rows are grouped back into one VCF
//! record per site using `alt_index` order; the stored `info_raw` and
//! `gt_raw` fields make this lossless for fields we preserve.
//!
//! Supports single-sample export and multi-sample merge export (the
//! `bcftools merge` replacement for samples already in the store).

use std::collections::BTreeMap;
use std::io::{BufWriter, Write};

use arrow_array::{Array, Float32Array, RecordBatch, StringArray, UInt32Array, UInt64Array};
use futures::TryStreamExt;
use genolance_core::{Error, Result, ResultExt};
use lancedb::query::{ExecutableQuery, QueryBase};

use genolance_core::schema::{SAMPLES_TABLE, VARIANTS_TABLE};
use genolance_core::store::Store;

use crate::common::{chrom_predicate, missing_column, normalize_chrom, require_samples};

/// Run `genolance export`.
///
/// Writes a VCF for `sample_name` to `output_path` (or stdout if `None`).
/// Optional region filters mirror `genolance query`.
pub async fn run(
    store_path: &str,
    sample_name: &str,
    chrom: Option<&str>,
    start: Option<u64>,
    end: Option<u64>,
    output_path: Option<&str>,
) -> Result<()> {
    let store = Store::open(store_path).await?;

    let var_tables = store.variants.table_names().execute().await?;
    if !var_tables.iter().any(|n| n == VARIANTS_TABLE) {
        return Err(Error::missing_table(
            VARIANTS_TABLE,
            "ingest a sample VCF first",
        ));
    }

    // 1. Header. The sample must be registered; an unknown name is an error,
    // not an empty VCF.
    require_samples(&store, std::slice::from_ref(&sample_name.to_string())).await?;
    let header_text = load_sample_header(&store, sample_name).await?;

    // 2. Variants for this sample, with optional region filter.
    let variants = store.variants.open_table(VARIANTS_TABLE).execute().await?;
    let mut preds: Vec<String> = Vec::new();
    preds.push(format!("sample_name = '{}'", sql_escape(sample_name)));
    if let Some(c) = chrom {
        preds.push(chrom_predicate(c));
    }
    if let Some(s) = start {
        preds.push(format!("pos >= {s}"));
    }
    if let Some(e) = end {
        preds.push(format!("pos <= {e}"));
    }
    let stream = variants
        .query()
        .only_if(preds.join(" AND "))
        .execute()
        .await?;
    let batches: Vec<RecordBatch> = stream.try_collect().await?;

    // 3. Flatten to owned rows, sort, group, emit.
    let mut rows: Vec<Row> = Vec::new();
    for b in &batches {
        collect_rows(b, &mut rows)?;
    }
    // Stable sort by (chrom, pos, ref, alt_index) - rows of the same site
    // land contiguous in alt_index order for group_consecutive().
    rows.sort_by(|a, b| {
        (chrom_key(&a.chrom), a.pos, &a.ref_allele, a.alt_index).cmp(&(
            chrom_key(&b.chrom),
            b.pos,
            &b.ref_allele,
            b.alt_index,
        ))
    });

    if rows.is_empty() && (chrom.is_some() || start.is_some() || end.is_some()) {
        return Err(Error::region_empty(
            format!("sample '{sample_name}'"),
            region_text(chrom, start, end),
        ));
    }

    // Group rows into source records before any output is created, so a
    // store that cannot be exported faithfully fails without a partial file.
    let records = split_records(&rows)?;

    let writer: Box<dyn Write> = match output_path {
        Some(p) => Box::new(BufWriter::new(
            std::fs::File::create(p).with_context(|| format!("creating {p}"))?,
        )),
        None => Box::new(BufWriter::new(std::io::stdout().lock())),
    };
    let mut w = writer;

    // 4. Emit header; the #CHROM line carries the stored sample name.
    if let Some(text) = header_text {
        write_header_for_sample(&mut w, &text, sample_name)?;
    } else {
        write_minimal_header(&mut w, sample_name)?;
    }

    // 5. Emit data lines, one per source record.
    let mut sites = 0usize;
    for site in &records {
        write_site(&mut w, site, sample_name)?;
        sites += 1;
    }
    w.flush()?;
    eprintln!(
        "[export] wrote {sites} sites ({} rows) for {sample_name}",
        rows.len()
    );
    Ok(())
}

// ---------------------------------------------------------------------------
// Sample header lookup
// ---------------------------------------------------------------------------

async fn load_sample_header(store: &Store, sample_name: &str) -> Result<Option<String>> {
    let samples = store.conn.open_table(SAMPLES_TABLE).execute().await?;
    let stream = samples
        .query()
        .only_if(format!("sample_name = '{}'", sql_escape(sample_name)))
        .execute()
        .await?;
    let batches: Vec<RecordBatch> = stream.try_collect().await?;
    // Use the most recent ingest - in practice there's one row per sample.
    for b in batches.iter().rev() {
        let header = b
            .column_by_name("vcf_header")
            .and_then(|c| c.as_any().downcast_ref::<StringArray>());
        if let Some(h) = header {
            if h.len() > 0 {
                return Ok(Some(h.value(h.len() - 1).to_string()));
            }
        }
    }
    Ok(None)
}

/// Write a stored header with its `#CHROM` line rewritten to name exactly
/// the exported sample.
fn write_header_for_sample<W: Write>(w: &mut W, text: &str, sample_name: &str) -> Result<()> {
    let mut saw_chrom = false;
    for line in text.lines() {
        if line.starts_with("#CHROM") {
            let fixed: Vec<&str> = line.split('\t').take(9).collect();
            writeln!(w, "{}\t{sample_name}", fixed.join("\t"))?;
            saw_chrom = true;
        } else {
            writeln!(w, "{line}")?;
        }
    }
    if !saw_chrom {
        return Err(invalid!(
            "stored header for '{sample_name}' has no #CHROM line"
        ));
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Row extraction + grouping
// ---------------------------------------------------------------------------

struct Row {
    chrom: String,
    pos: u64,
    ids: Option<String>,
    ref_allele: String,
    alt_allele: String,
    alt_index: u32,
    alt_count: u32,
    quality: Option<f32>,
    filter: Option<String>,
    gt_raw: Option<String>,
    read_depth: Option<u32>,
    format_ad: Option<String>,
    format_gq: Option<u32>,
    format_pl: Option<String>,
    format_extra: Option<String>,
    format_key_order: Option<String>,
    info_raw: Option<String>,
}

/// Read the rows of one batch. Every column must be present with its
/// expected type; a store that drifted from the schema is an error.
fn collect_rows(b: &RecordBatch, out: &mut Vec<Row>) -> Result<()> {
    let t = VARIANTS_TABLE;
    let chrom = col_str(b, "chrom").ok_or_else(|| missing_column(t, "chrom"))?;
    let pos = col_u64(b, "pos").ok_or_else(|| missing_column(t, "pos"))?;
    let ids = Some(col_str(b, "ids").ok_or_else(|| missing_column(t, "ids"))?);
    let ref_a = col_str(b, "ref_allele").ok_or_else(|| missing_column(t, "ref_allele"))?;
    let alt_a = col_str(b, "alt_allele").ok_or_else(|| missing_column(t, "alt_allele"))?;
    let alt_index = col_u32(b, "alt_index").ok_or_else(|| missing_column(t, "alt_index"))?;
    let alt_count = col_u32(b, "alt_count").ok_or_else(|| missing_column(t, "alt_count"))?;
    let quality = Some(col_f32(b, "quality").ok_or_else(|| missing_column(t, "quality"))?);
    let filter = Some(col_str(b, "filter").ok_or_else(|| missing_column(t, "filter"))?);
    let gt_raw = Some(col_str(b, "gt_raw").ok_or_else(|| missing_column(t, "gt_raw"))?);
    let rd = Some(col_u32(b, "read_depth").ok_or_else(|| missing_column(t, "read_depth"))?);
    let ad = Some(col_str(b, "format_ad").ok_or_else(|| missing_column(t, "format_ad"))?);
    let gq = Some(col_u32(b, "format_gq").ok_or_else(|| missing_column(t, "format_gq"))?);
    let pl = Some(col_str(b, "format_pl").ok_or_else(|| missing_column(t, "format_pl"))?);
    let extra = Some(col_str(b, "format_extra").ok_or_else(|| missing_column(t, "format_extra"))?);
    let key_order =
        Some(col_str(b, "format_key_order").ok_or_else(|| missing_column(t, "format_key_order"))?);
    let info_raw = Some(col_str(b, "info_raw").ok_or_else(|| missing_column(t, "info_raw"))?);

    for i in 0..b.num_rows() {
        out.push(Row {
            chrom: chrom.value(i).to_string(),
            pos: pos.value(i),
            ids: str_opt(ids, i),
            ref_allele: ref_a.value(i).to_string(),
            alt_allele: alt_a.value(i).to_string(),
            alt_index: alt_index.value(i),
            alt_count: alt_count.value(i),
            quality: f32_opt(quality, i),
            filter: str_opt(filter, i),
            gt_raw: str_opt(gt_raw, i),
            read_depth: u32_opt(rd, i),
            format_ad: str_opt(ad, i),
            format_gq: u32_opt(gq, i),
            format_pl: str_opt(pl, i),
            format_extra: str_opt(extra, i),
            format_key_order: str_opt(key_order, i),
            info_raw: str_opt(info_raw, i),
        });
    }
    Ok(())
}

/// Rebuild source records from the per-ALT rows of one `(chrom, pos, ref)`.
///
/// A record is either one row with `alt_count == 1` (also what a
/// split-biallelic input looks like, several of which can share a position)
/// or `alt_count` rows with `alt_index` 0..alt_count. Anything else cannot be
/// attributed to records without guessing and is an error.
fn records_at_site<'a>(site: &[&'a Row]) -> Result<Vec<Vec<&'a Row>>> {
    let head = site[0];
    let ambiguous = |why: String| -> Error {
        Error::malformed_record(head.chrom.clone(), head.pos, "alt rows", why)
    };
    let mut records: Vec<Vec<&'a Row>> = Vec::new();
    let mut by_count: BTreeMap<u32, Vec<&'a Row>> = BTreeMap::new();
    for &r in site {
        if r.alt_count == 0 {
            return Err(ambiguous("row has alt_count 0".into()));
        }
        if r.alt_count == 1 {
            if r.alt_index != 0 {
                return Err(ambiguous(format!(
                    "alt_index {} with alt_count 1",
                    r.alt_index
                )));
            }
            records.push(vec![r]);
        } else {
            by_count.entry(r.alt_count).or_default().push(r);
        }
    }
    for (k, rows) in by_count {
        if rows.len() != k as usize {
            return Err(ambiguous(format!(
                "{} rows with alt_count {k}; cannot tell which source records they belong to",
                rows.len()
            )));
        }
        let mut sorted = rows;
        sorted.sort_by_key(|r| r.alt_index);
        if sorted
            .iter()
            .enumerate()
            .any(|(i, r)| r.alt_index as usize != i)
        {
            return Err(ambiguous(format!(
                "alt_index values do not cover 0..{k} exactly once"
            )));
        }
        records.push(sorted);
    }
    for rec in &records {
        check_record_consistent(rec).map_err(ambiguous)?;
    }
    Ok(records)
}

/// Rows of one record duplicate its site-level fields; they must agree.
fn check_record_consistent(rec: &[&Row]) -> std::result::Result<(), String> {
    let h = rec[0];
    for r in &rec[1..] {
        if r.gt_raw != h.gt_raw
            || r.info_raw != h.info_raw
            || r.quality != h.quality
            || r.filter != h.filter
            || r.ids != h.ids
            || r.format_ad != h.format_ad
            || r.format_pl != h.format_pl
            || r.read_depth != h.read_depth
            || r.format_gq != h.format_gq
            || r.format_extra != h.format_extra
        {
            return Err(format!(
                "rows for ALT {} and {} disagree on site-level fields",
                h.alt_allele, r.alt_allele
            ));
        }
    }
    Ok(())
}

/// Split sorted rows into source records (see [`records_at_site`]).
fn split_records(rows: &[Row]) -> Result<Vec<Vec<&Row>>> {
    let refs: Vec<&Row> = rows.iter().collect();
    let mut out = Vec::new();
    let mut start = 0usize;
    while start < rows.len() {
        let head = &rows[start];
        let mut end = start + 1;
        while end < rows.len()
            && rows[end].chrom == head.chrom
            && rows[end].pos == head.pos
            && rows[end].ref_allele == head.ref_allele
        {
            end += 1;
        }
        out.extend(records_at_site(&refs[start..end])?);
        start = end;
    }
    Ok(out)
}

// ---------------------------------------------------------------------------
// VCF line emission
// ---------------------------------------------------------------------------

fn write_site<W: Write>(w: &mut W, site: &[&Row], _sample_name: &str) -> Result<()> {
    let head = site[0];

    // ALTs in alt_index order (group_consecutive already sorted).
    let alts: Vec<&str> = site.iter().map(|r| r.alt_allele.as_str()).collect();

    let id_col = head.ids.as_deref().unwrap_or(".");
    let qual_col = head.quality.map(format_qual).unwrap_or_else(|| ".".into());
    let filter_col = head.filter.as_deref().unwrap_or(".");

    // Build a key→value map from stored fields so we can emit in any order.
    let info_col = head.info_raw.as_deref().unwrap_or(".");

    // Map of FORMAT key → value for quick reconstruction.
    let mut kv: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    if let Some(ref gt) = head.gt_raw {
        kv.insert("GT".into(), gt.clone());
    }
    if let Some(ref ad) = head.format_ad {
        kv.insert("AD".into(), ad.clone());
    }
    if let Some(dp) = head.read_depth {
        kv.insert("DP".into(), dp.to_string());
    }
    if let Some(gq) = head.format_gq {
        kv.insert("GQ".into(), gq.to_string());
    }
    if let Some(ref pl) = head.format_pl {
        kv.insert("PL".into(), pl.clone());
    }
    if let Some(extra) = head.format_extra.as_deref() {
        for part in extra.split(';') {
            if part.is_empty() {
                continue;
            }
            match part.split_once('=') {
                Some((k, v)) => {
                    kv.insert(k.to_string(), crate::ingest::decode_extra_value(v));
                }
                None => {
                    kv.insert(part.to_string(), ".".to_string());
                }
            }
        }
    }

    // Emit FORMAT keys in original source order if we have it; else fall back
    // to canonical GT:AD:DP:GQ:PL:<extras> order.
    let fmt_keys: Vec<String> = if let Some(order) = head.format_key_order.as_deref() {
        order.split(':').map(str::to_string).collect()
    } else {
        let mut keys = vec!["GT", "AD", "DP", "GQ", "PL"]
            .into_iter()
            .filter(|k| kv.contains_key(*k))
            .map(str::to_string)
            .collect::<Vec<_>>();
        // append any extras not already listed
        for k in kv.keys() {
            if !["GT", "AD", "DP", "GQ", "PL"].contains(&k.as_str()) {
                keys.push(k.clone());
            }
        }
        keys
    };
    let fmt_vals: Vec<String> = fmt_keys
        .iter()
        .map(|k| kv.get(k).cloned().unwrap_or_else(|| ".".into()))
        .collect();

    let (format_col, sample_col) = if fmt_keys.is_empty() {
        (".".to_string(), ".".to_string())
    } else {
        (fmt_keys.join(":"), fmt_vals.join(":"))
    };

    writeln!(
        w,
        "{chrom}\t{pos}\t{id}\t{ref_a}\t{alt}\t{qual}\t{filter}\t{info}\t{fmt}\t{sample}",
        chrom = head.chrom,
        pos = head.pos,
        id = id_col,
        ref_a = head.ref_allele,
        alt = alts.join(","),
        qual = qual_col,
        filter = filter_col,
        info = info_col,
        fmt = format_col,
        sample = sample_col,
    )?;
    Ok(())
}

fn write_minimal_header<W: Write>(w: &mut W, sample_name: &str) -> Result<()> {
    // Fallback when no header was recorded for this sample.
    writeln!(w, "##fileformat=VCFv4.2")?;
    writeln!(w, "##source=genolance-export")?;
    writeln!(
        w,
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">"
    )?;
    writeln!(
        w,
        "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled likelihoods\">"
    )?;
    writeln!(
        w,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}"
    )?;
    Ok(())
}

fn format_qual(q: f32) -> String {
    if q.fract() == 0.0 {
        format!("{}", q as i64)
    } else {
        format!("{q}")
    }
}

// ---------------------------------------------------------------------------
// Arrow column helpers
// ---------------------------------------------------------------------------

fn col_str<'a>(b: &'a RecordBatch, name: &str) -> Option<&'a StringArray> {
    b.column_by_name(name)
        .and_then(|c| c.as_any().downcast_ref::<StringArray>())
}
fn col_u64<'a>(b: &'a RecordBatch, name: &str) -> Option<&'a UInt64Array> {
    b.column_by_name(name)
        .and_then(|c| c.as_any().downcast_ref::<UInt64Array>())
}
fn col_u32<'a>(b: &'a RecordBatch, name: &str) -> Option<&'a UInt32Array> {
    b.column_by_name(name)
        .and_then(|c| c.as_any().downcast_ref::<UInt32Array>())
}
fn col_f32<'a>(b: &'a RecordBatch, name: &str) -> Option<&'a Float32Array> {
    b.column_by_name(name)
        .and_then(|c| c.as_any().downcast_ref::<Float32Array>())
}

fn str_opt(a: Option<&StringArray>, i: usize) -> Option<String> {
    let arr = a?;
    if arr.is_null(i) {
        None
    } else {
        Some(arr.value(i).to_string())
    }
}
fn u32_opt(a: Option<&UInt32Array>, i: usize) -> Option<u32> {
    let arr = a?;
    if arr.is_null(i) {
        None
    } else {
        Some(arr.value(i))
    }
}
fn f32_opt(a: Option<&Float32Array>, i: usize) -> Option<f32> {
    let arr = a?;
    if arr.is_null(i) {
        None
    } else {
        Some(arr.value(i))
    }
}

fn region_text(chrom: Option<&str>, start: Option<u64>, end: Option<u64>) -> String {
    let mut parts = Vec::new();
    if let Some(c) = chrom {
        parts.push(format!("chrom {c}"));
    }
    if let Some(s) = start {
        parts.push(format!("start {s}"));
    }
    if let Some(e) = end {
        parts.push(format!("end {e}"));
    }
    parts.join(", ")
}

fn sql_escape(s: &str) -> String {
    s.replace('\'', "''")
}

/// Sort key for chromosome names: numeric chroms first (by number), then
/// X, Y, M/MT, then everything else alphabetically. Matches the order
/// most tools and references use.
fn chrom_key(c: &str) -> (u8, u32, String) {
    let bare = c.trim_start_matches("chr");
    if let Ok(n) = bare.parse::<u32>() {
        (0, n, String::new())
    } else {
        match bare {
            "X" => (1, 0, String::new()),
            "Y" => (1, 1, String::new()),
            "M" | "MT" => (1, 2, String::new()),
            other => (2, 0, other.to_string()),
        }
    }
}

// ---------------------------------------------------------------------------
// Multi-sample merge export
// ---------------------------------------------------------------------------

/// Run multi-sample merge export.
///
/// Queries each requested sample, unions ALT alleles per site, renumbers
/// allele indices consistently, and emits one VCF record per site with
/// one FORMAT/SAMPLE column per listed sample (samples without the site
/// get no-calls). GT, AD (Number=R) and PL (Number=G) are remapped to the
/// merged allele order. Field set is normalized to `GT:AD:DP:GQ:PL`;
/// non-standard FORMAT fields (`format_extra`) are dropped because their
/// cardinality can disagree between samples at the same site. INFO is kept
/// only when every sample at the site agrees on it and on the ALT list;
/// otherwise it is `.`. FILTER is the union of the samples' non-PASS filters.
///
/// A sample must have at most one source record per `(chrom, pos, ref)`;
/// several (split-biallelic input) cannot be merged without guessing and
/// are an error.
pub async fn run_merge(
    store_path: &str,
    sample_names: &[String],
    chrom: Option<&str>,
    start: Option<u64>,
    end: Option<u64>,
    output_path: Option<&str>,
) -> Result<()> {
    if sample_names.len() < 2 {
        return Err(invalid!(
            "merging needs at least two samples; got {}",
            sample_names.len()
        ));
    }
    let store = Store::open(store_path).await?;
    let var_tables = store.variants.table_names().execute().await?;
    if !var_tables.iter().any(|n| n == VARIANTS_TABLE) {
        return Err(Error::missing_table(VARIANTS_TABLE, "ingest it first"));
    }
    require_samples(&store, sample_names).await?;
    let variants = store.variants.open_table(VARIANTS_TABLE).execute().await?;

    // Per-sample records, indexed by sample position in sample_names.
    // Key: (chrom, pos, ref) -> the sample's single record (one Row per alt).
    let mut per_sample: Vec<BTreeMap<(String, u64, String), Vec<Row>>> =
        (0..sample_names.len()).map(|_| BTreeMap::new()).collect();

    for (idx, s) in sample_names.iter().enumerate() {
        let mut preds = vec![format!("sample_name = '{}'", sql_escape(s))];
        if let Some(c) = chrom {
            preds.push(chrom_predicate(c));
        }
        if let Some(v) = start {
            preds.push(format!("pos >= {v}"));
        }
        if let Some(v) = end {
            preds.push(format!("pos <= {v}"));
        }
        let batches: Vec<RecordBatch> = variants
            .query()
            .only_if(preds.join(" AND "))
            .execute()
            .await?
            .try_collect()
            .await?;
        let mut rows: Vec<Row> = Vec::new();
        for b in &batches {
            collect_rows(b, &mut rows)?;
        }
        rows.sort_by(|a, b| {
            (&a.chrom, a.pos, &a.ref_allele, a.alt_index).cmp(&(
                &b.chrom,
                b.pos,
                &b.ref_allele,
                b.alt_index,
            ))
        });
        // Group by site, then require exactly one record per site.
        let mut groups: BTreeMap<(String, u64, String), Vec<Row>> = BTreeMap::new();
        for r in rows {
            groups
                .entry((normalize_chrom(&r.chrom), r.pos, r.ref_allele.clone()))
                .or_default()
                .push(r);
        }
        for (head_key, group) in groups {
            let n_records = {
                let refs: Vec<&Row> = group.iter().collect();
                records_at_site(&refs)?.len()
            };
            if n_records != 1 {
                return Err(Error::malformed_record(head_key.0, head_key.1, "merge", format!(
                        "sample '{s}' has {n_records} separate records at this position and reference (split multi-allelic input); they cannot be merged into one column without guessing"
                    )));
            }
            per_sample[idx].insert(head_key, group);
        }
    }

    // Union of site keys across samples, sorted by (chrom_key, pos, ref).
    let mut sites: Vec<(String, u64, String)> =
        per_sample.iter().flat_map(|m| m.keys().cloned()).collect();
    sites.sort_by(|a, b| (chrom_key(&a.0), a.1, &a.2).cmp(&(chrom_key(&b.0), b.1, &b.2)));
    sites.dedup();

    if sites.is_empty() && (chrom.is_some() || start.is_some() || end.is_some()) {
        return Err(Error::region_empty(
            "the listed samples",
            region_text(chrom, start, end),
        ));
    }

    // Build every line first so a bad site fails before any output exists.
    let mut lines: Vec<String> = Vec::new();
    for key in &sites {
        if let Some(line) = merged_site_line(&per_sample, key)? {
            lines.push(line);
        }
    }

    let writer: Box<dyn Write> = match output_path {
        Some(p) => Box::new(BufWriter::new(
            std::fs::File::create(p).with_context(|| format!("creating {p}"))?,
        )),
        None => Box::new(BufWriter::new(std::io::stdout().lock())),
    };
    let mut w = writer;
    write_merge_header(&mut w, &store, sample_names).await?;
    for l in &lines {
        writeln!(w, "{l}")?;
    }
    w.flush()?;
    eprintln!(
        "[export] wrote {} sites for {} samples",
        lines.len(),
        sample_names.len()
    );
    Ok(())
}

/// Emit a merged header: prefer the first sample's stored header (for
/// `##contig`, `##FORMAT`, `##INFO` lines) and swap the `#CHROM` line
/// for our multi-sample variant.
async fn write_merge_header<W: Write>(w: &mut W, store: &Store, samples: &[String]) -> Result<()> {
    let root_tables = store.conn.table_names().execute().await?;
    let mut wrote_meta = false;
    if root_tables.iter().any(|n| n == SAMPLES_TABLE) {
        if let Some(header) = load_sample_header(store, &samples[0]).await? {
            // Strip the stored header's trailing #CHROM line; we write our own.
            for line in header.lines() {
                if line.starts_with("#CHROM") {
                    continue;
                }
                writeln!(w, "{line}")?;
            }
            writeln!(w, "##genolance_merge_samples={}", samples.join(","))?;
            wrote_meta = true;
        }
    }
    if !wrote_meta {
        writeln!(w, "##fileformat=VCFv4.2")?;
        writeln!(w, "##source=genolance-merge-export")?;
        writeln!(
            w,
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        )?;
        writeln!(
            w,
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">"
        )?;
        writeln!(
            w,
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">"
        )?;
        writeln!(
            w,
            "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">"
        )?;
        writeln!(
            w,
            "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled likelihoods\">"
        )?;
    }
    writeln!(
        w,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}",
        samples.join("\t")
    )?;
    Ok(())
}

type SiteKey = (String, u64, String);

fn merged_site_line(
    per_sample: &[BTreeMap<SiteKey, Vec<Row>>],
    key: &(String, u64, String),
) -> Result<Option<String>> {
    let mut present: Vec<Vec<&Row>> = Vec::new();
    for sm in per_sample {
        if let Some(rows) = sm.get(key) {
            let mut sorted: Vec<&Row> = rows.iter().collect();
            sorted.sort_by_key(|r| r.alt_index);
            present.push(sorted);
        }
    }
    // Union of ALTs across samples, in first-seen order.
    let mut merged_alts: Vec<String> = Vec::new();
    for rows in &present {
        for r in rows {
            if !merged_alts.iter().any(|a| a == &r.alt_allele) {
                merged_alts.push(r.alt_allele.clone());
            }
        }
    }
    if merged_alts.is_empty() {
        return Ok(None);
    }

    // Site-level fields: keep a value only where every sample agrees.
    let agree = |f: &dyn Fn(&Row) -> Option<String>| -> Option<String> {
        let first = f(present[0][0]);
        if present.iter().all(|rows| f(rows[0]) == first) {
            first
        } else {
            None
        }
    };
    let same_alts = present
        .iter()
        .all(|rows| rows.iter().map(|r| &r.alt_allele).eq(merged_alts.iter()));
    let rep_info = if same_alts {
        agree(&|r| r.info_raw.clone())
    } else {
        None
    };
    let rep_qual = agree(&|r| r.quality.map(format_qual));
    let mut ids: Vec<String> = Vec::new();
    for rows in &present {
        if let Some(i) = &rows[0].ids {
            for part in i.split(';') {
                if !ids.iter().any(|e| e == part) {
                    ids.push(part.to_string());
                }
            }
        }
    }
    let mut failed: Vec<String> = Vec::new();
    let mut any_pass = false;
    for rows in &present {
        match rows[0].filter.as_deref() {
            Some("PASS") => any_pass = true,
            Some(f) => {
                for part in f.split(';') {
                    if !failed.iter().any(|e| e == part) {
                        failed.push(part.to_string());
                    }
                }
            }
            None => {}
        }
    }
    let rep_filter = if !failed.is_empty() {
        failed.join(";")
    } else if any_pass {
        "PASS".to_string()
    } else {
        ".".to_string()
    };

    let mut sample_cols: Vec<String> = Vec::new();
    for sm in per_sample {
        let col = match sm.get(key) {
            Some(rows) => build_merged_sample_column(rows, &merged_alts, key)?,
            None => missing_sample_column(),
        };
        sample_cols.push(col);
    }

    Ok(Some(format!(
        "{chrom}\t{pos}\t{id}\t{ref_a}\t{alt}\t{qual}\t{filter}\t{info}\tGT:AD:DP:GQ:PL\t{sample}",
        chrom = present[0][0].chrom,
        pos = key.1,
        id = if ids.is_empty() {
            ".".to_string()
        } else {
            ids.join(";")
        },
        ref_a = key.2,
        alt = merged_alts.join(","),
        qual = rep_qual.unwrap_or_else(|| ".".into()),
        filter = rep_filter,
        info = rep_info.unwrap_or_else(|| ".".into()),
        sample = sample_cols.join("\t"),
    )))
}

/// Build a `GT:AD:DP:GQ:PL` column for one sample at one site, remapping
/// the sample's local allele indices into the merged ALT ordering. GT, AD
/// (Number=R) and PL (Number=G) are all re-indexed.
fn build_merged_sample_column(
    rows: &[Row],
    merged_alts: &[String],
    key: &(String, u64, String),
) -> Result<String> {
    let bad = |field: &str, reason: String| -> Error {
        Error::malformed_record(key.0.clone(), key.1, field, reason)
    };
    let mut sorted: Vec<&Row> = rows.iter().collect();
    sorted.sort_by_key(|r| r.alt_index);
    // local allele index (0 = REF) -> merged allele index
    let mut local_to_merged: Vec<usize> = vec![0];
    for r in &sorted {
        let p = merged_alts
            .iter()
            .position(|a| a == &r.alt_allele)
            .ok_or_else(|| {
                bad(
                    "ALT",
                    format!("allele {} missing from merged list", r.alt_allele),
                )
            })?;
        local_to_merged.push(p + 1);
    }
    let identity = local_to_merged.len() == merged_alts.len() + 1
        && local_to_merged.iter().enumerate().all(|(i, m)| i == *m);
    let n_local = local_to_merged.len();
    let n_merged = merged_alts.len() + 1;

    let head = sorted[0];
    let gt = match &head.gt_raw {
        Some(g) => remap_gt(g, &local_to_merged).map_err(|e| bad("FORMAT/GT", e))?,
        None => "./.".to_string(),
    };
    let ad = match &head.format_ad {
        None => ".".to_string(),
        Some(v) if identity => v.clone(),
        Some(v) => {
            let vals: Vec<&str> = v.split(',').collect();
            if vals.len() != n_local {
                return Err(bad(
                    "FORMAT/AD",
                    format!("{} value(s) for {n_local} alleles", vals.len()),
                ));
            }
            let mut out = vec!["."; n_merged];
            for (i, val) in vals.iter().enumerate() {
                out[local_to_merged[i]] = val;
            }
            out.join(",")
        }
    };
    let pl = match &head.format_pl {
        None => ".".to_string(),
        Some(v) if identity => v.clone(),
        Some(v) => {
            let vals: Vec<&str> = v.split(',').collect();
            if vals.len() != n_local * (n_local + 1) / 2 {
                return Err(bad(
                    "FORMAT/PL",
                    format!(
                        "{} value(s) is not the diploid genotype count for {n_local} alleles; cannot re-index",
                        vals.len()
                    ),
                ));
            }
            let mut out = vec!["."; n_merged * (n_merged + 1) / 2];
            let mut k = 0usize;
            for b in 0..n_local {
                for a in 0..=b {
                    let (ma, mb) = {
                        let (x, y) = (local_to_merged[a], local_to_merged[b]);
                        (x.min(y), x.max(y))
                    };
                    out[mb * (mb + 1) / 2 + ma] = vals[k];
                    k += 1;
                }
            }
            out.join(",")
        }
    };
    let dp = head
        .read_depth
        .map(|x| x.to_string())
        .unwrap_or_else(|| ".".to_string());
    let gq = head
        .format_gq
        .map(|x| x.to_string())
        .unwrap_or_else(|| ".".to_string());
    Ok(format!("{gt}:{ad}:{dp}:{gq}:{pl}"))
}

fn missing_sample_column() -> String {
    "./.:.:.:.:.".to_string()
}

/// Remap allele indices in a GT string. Preserves phasing separators.
/// `.` passes through; any other token must be a local allele index.
fn remap_gt(gt: &str, local_to_merged: &[usize]) -> std::result::Result<String, String> {
    let mut out = String::with_capacity(gt.len());
    let mut token = String::new();
    for ch in gt.chars() {
        if ch == '|' || ch == '/' {
            flush_allele(&mut token, &mut out, local_to_merged)?;
            out.push(ch);
        } else {
            token.push(ch);
        }
    }
    flush_allele(&mut token, &mut out, local_to_merged)?;
    Ok(out)
}

fn flush_allele(
    token: &mut String,
    out: &mut String,
    local_to_merged: &[usize],
) -> std::result::Result<(), String> {
    if token.is_empty() {
        return Ok(());
    }
    if token == "." {
        out.push('.');
    } else {
        let n: usize = token
            .parse()
            .map_err(|_| format!("allele '{token}' is not an index"))?;
        let merged = local_to_merged.get(n).ok_or_else(|| {
            format!(
                "allele {n} is beyond the sample's {} ALT(s)",
                local_to_merged.len() - 1
            )
        })?;
        out.push_str(&merged.to_string());
    }
    token.clear();
    Ok(())
}
