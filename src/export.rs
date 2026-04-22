//! VCF export — reconstruct a VCF for a single sample (or region) from
//! the Lance store. Per-alt split rows are grouped back into one VCF
//! record per site using `alt_index` order; the stored `info_raw` and
//! `gt_raw` fields make this lossless for fields we preserve.
//!
//! Supports single-sample export and multi-sample merge export (the
//! `bcftools merge` replacement for samples already in the store).

use std::collections::BTreeMap;
use std::io::{BufWriter, Write};

use anyhow::{anyhow, Context, Result};
use arrow_array::{Array, Float32Array, RecordBatch, StringArray, UInt32Array, UInt64Array};
use futures::TryStreamExt;
use lancedb::query::{ExecutableQuery, QueryBase};

use crate::schema::{SAMPLES_TABLE, VARIANTS_TABLE};
use crate::store::Store;

/// Run `biolance export`.
///
/// Writes a VCF for `sample_name` to `output_path` (or stdout if `None`).
/// Optional region filters mirror `biolance query`.
pub async fn run(
    store_path: &str,
    sample_name: &str,
    chrom: Option<&str>,
    start: Option<u64>,
    end: Option<u64>,
    output_path: Option<&str>,
) -> Result<()> {
    let store = Store::open(store_path).await?;

    let tables = store.conn.table_names().execute().await?;
    if !tables.iter().any(|n| n == VARIANTS_TABLE) {
        return Err(anyhow!(
            "store has no '{VARIANTS_TABLE}' table; ingest a sample VCF first"
        ));
    }

    // 1. Header
    let header_text = if tables.iter().any(|n| n == SAMPLES_TABLE) {
        load_sample_header(&store, sample_name).await?
    } else {
        None
    };

    // 2. Variants for this sample, with optional region filter.
    let variants = store.conn.open_table(VARIANTS_TABLE).execute().await?;
    let mut preds: Vec<String> = Vec::new();
    preds.push(format!("sample_name = '{}'", sql_escape(sample_name)));
    if let Some(c) = chrom {
        preds.push(format!("chrom = '{}'", sql_escape(c)));
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
        collect_rows(b, &mut rows);
    }
    // Stable sort by (chrom, pos, ref, alt_index) — rows of the same site
    // land contiguous in alt_index order for group_consecutive().
    rows.sort_by(|a, b| {
        (chrom_key(&a.chrom), a.pos, &a.ref_allele, a.alt_index).cmp(&(
            chrom_key(&b.chrom),
            b.pos,
            &b.ref_allele,
            b.alt_index,
        ))
    });

    let writer: Box<dyn Write> = match output_path {
        Some(p) => Box::new(BufWriter::new(
            std::fs::File::create(p).with_context(|| format!("creating {p}"))?,
        )),
        None => Box::new(BufWriter::new(std::io::stdout().lock())),
    };
    let mut w = writer;

    // 4. Emit header.
    if let Some(text) = header_text {
        // Header already ends in newline(s) from noodles serializer.
        w.write_all(text.as_bytes())?;
    } else {
        write_minimal_header(&mut w, sample_name)?;
    }

    // 5. Emit data lines, one per site.
    let mut sites = 0usize;
    for site in group_consecutive(&rows) {
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
    // Use the most recent ingest — in practice there's one row per sample.
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

fn collect_rows(b: &RecordBatch, out: &mut Vec<Row>) {
    let chrom = col_str(b, "chrom");
    let pos = col_u64(b, "pos");
    let ids = col_str(b, "ids");
    let ref_a = col_str(b, "ref_allele");
    let alt_a = col_str(b, "alt_allele");
    let alt_index = col_u32(b, "alt_index");
    let quality = col_f32(b, "quality");
    let filter = col_str(b, "filter");
    let gt_raw = col_str(b, "gt_raw");
    let rd = col_u32(b, "read_depth");
    let ad = col_str(b, "format_ad");
    let gq = col_u32(b, "format_gq");
    let pl = col_str(b, "format_pl");
    let extra = col_str(b, "format_extra");
    let key_order = col_str(b, "format_key_order");
    let info_raw = col_str(b, "info_raw");

    let (Some(chrom), Some(pos), Some(ref_a), Some(alt_a), Some(alt_index)) =
        (chrom, pos, ref_a, alt_a, alt_index)
    else {
        return;
    };

    for i in 0..b.num_rows() {
        out.push(Row {
            chrom: chrom.value(i).to_string(),
            pos: pos.value(i),
            ids: str_opt(ids, i),
            ref_allele: ref_a.value(i).to_string(),
            alt_allele: alt_a.value(i).to_string(),
            alt_index: alt_index.value(i),
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
}

/// Yield slices of consecutive rows sharing `(chrom, pos, ref_allele)`.
fn group_consecutive(rows: &[Row]) -> impl Iterator<Item = &[Row]> {
    let mut start = 0usize;
    std::iter::from_fn(move || {
        if start >= rows.len() {
            return None;
        }
        let head = &rows[start];
        let mut end = start + 1;
        while end < rows.len() {
            let r = &rows[end];
            if r.chrom == head.chrom && r.pos == head.pos && r.ref_allele == head.ref_allele {
                end += 1;
            } else {
                break;
            }
        }
        let slice = &rows[start..end];
        start = end;
        Some(slice)
    })
}

// ---------------------------------------------------------------------------
// VCF line emission
// ---------------------------------------------------------------------------

fn write_site<W: Write>(w: &mut W, site: &[Row], _sample_name: &str) -> Result<()> {
    let head = &site[0];

    // ALTs in alt_index order (group_consecutive already sorted).
    let alts: Vec<&str> = site.iter().map(|r| r.alt_allele.as_str()).collect();

    let id_col = head.ids.as_deref().unwrap_or(".");
    let qual_col = head
        .quality
        .map(|q| format_qual(q))
        .unwrap_or_else(|| ".".into());
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
                    kv.insert(k.to_string(), v.to_string());
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
    writeln!(w, "##source=biolance-export")?;
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
/// one FORMAT/SAMPLE column per listed sample (missing samples become
/// `./.`). Field set is normalized to `GT:AD:DP:GQ:PL`; non-standard
/// FORMAT fields (`format_extra`) are dropped during merge because their
/// cardinality can disagree between samples at the same site.
pub async fn run_merge(
    store_path: &str,
    sample_names: &[String],
    chrom: Option<&str>,
    start: Option<u64>,
    end: Option<u64>,
    output_path: Option<&str>,
) -> Result<()> {
    if sample_names.len() < 2 {
        return Err(anyhow!(
            "--merge requires at least two samples; got {}",
            sample_names.len()
        ));
    }
    let store = Store::open(store_path).await?;
    let tables = store.conn.table_names().execute().await?;
    if !tables.iter().any(|n| n == VARIANTS_TABLE) {
        return Err(anyhow!("store has no '{VARIANTS_TABLE}' table"));
    }
    let variants = store.conn.open_table(VARIANTS_TABLE).execute().await?;

    // Per-sample row collections, indexed by sample position in sample_names.
    // Key: (chrom, pos, ref) -> Vec<Row>  (one Row per alt for that sample).
    let mut per_sample: Vec<BTreeMap<(String, u64, String), Vec<Row>>> =
        (0..sample_names.len()).map(|_| BTreeMap::new()).collect();

    for (idx, s) in sample_names.iter().enumerate() {
        let mut preds = vec![format!("sample_name = '{}'", sql_escape(s))];
        if let Some(c) = chrom {
            preds.push(format!("chrom = '{}'", sql_escape(c)));
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
            collect_rows(b, &mut rows);
        }
        for r in rows {
            per_sample[idx]
                .entry((r.chrom.clone(), r.pos, r.ref_allele.clone()))
                .or_default()
                .push(r);
        }
    }

    // Union of site keys across samples, sorted by (chrom_key, pos, ref).
    let mut sites: Vec<(String, u64, String)> =
        per_sample.iter().flat_map(|m| m.keys().cloned()).collect();
    sites.sort_by(|a, b| (chrom_key(&a.0), a.1, &a.2).cmp(&(chrom_key(&b.0), b.1, &b.2)));
    sites.dedup();

    // Writer + header.
    let writer: Box<dyn Write> = match output_path {
        Some(p) => Box::new(BufWriter::new(
            std::fs::File::create(p).with_context(|| format!("creating {p}"))?,
        )),
        None => Box::new(BufWriter::new(std::io::stdout().lock())),
    };
    let mut w = writer;
    write_merge_header(&mut w, &store, sample_names).await?;

    let mut emitted = 0usize;
    for key in &sites {
        write_merged_site(&mut w, &per_sample, sample_names, key)?;
        emitted += 1;
    }
    w.flush()?;
    eprintln!(
        "[export --merge] wrote {emitted} sites for {} samples",
        sample_names.len()
    );
    Ok(())
}

/// Emit a merged header: prefer the first sample's stored header (for
/// `##contig`, `##FORMAT`, `##INFO` lines) and swap the `#CHROM` line
/// for our multi-sample variant.
async fn write_merge_header<W: Write>(w: &mut W, store: &Store, samples: &[String]) -> Result<()> {
    let tables = store.conn.table_names().execute().await?;
    let mut wrote_meta = false;
    if tables.iter().any(|n| n == SAMPLES_TABLE) {
        if let Some(header) = load_sample_header(store, &samples[0]).await? {
            // Strip the stored header's trailing #CHROM line; we write our own.
            for line in header.lines() {
                if line.starts_with("#CHROM") {
                    continue;
                }
                writeln!(w, "{line}")?;
            }
            writeln!(w, "##biolance_merge_samples={}", samples.join(","))?;
            wrote_meta = true;
        }
    }
    if !wrote_meta {
        writeln!(w, "##fileformat=VCFv4.2")?;
        writeln!(w, "##source=biolance-merge-export")?;
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

fn write_merged_site<W: Write>(
    w: &mut W,
    per_sample: &[BTreeMap<(String, u64, String), Vec<Row>>],
    sample_names: &[String],
    key: &(String, u64, String),
) -> Result<()> {
    // Union of ALTs across samples, ordered by (alt_index, first-seen).
    let mut merged_alts: Vec<String> = Vec::new();
    for sm in per_sample {
        if let Some(rows) = sm.get(key) {
            // Sort within-sample by alt_index for a stable order.
            let mut sorted: Vec<&Row> = rows.iter().collect();
            sorted.sort_by_key(|r| r.alt_index);
            for r in sorted {
                if !merged_alts.iter().any(|a| a == &r.alt_allele) {
                    merged_alts.push(r.alt_allele.clone());
                }
            }
        }
    }
    if merged_alts.is_empty() {
        return Ok(());
    }

    // Pick ID, QUAL, FILTER, INFO from the first sample that has it.
    // INFO is per-sample specific in our split storage, but we just take
    // the first as a representative; downstream tools should recompute.
    let mut rep_id = ".".to_string();
    let mut rep_qual = ".".to_string();
    let mut rep_filter = ".".to_string();
    let mut rep_info = ".".to_string();
    for sm in per_sample {
        if let Some(rows) = sm.get(key) {
            if let Some(r) = rows.iter().min_by_key(|r| r.alt_index) {
                if let Some(i) = &r.ids {
                    rep_id = i.clone();
                }
                if let Some(q) = r.quality {
                    rep_qual = format_qual(q);
                }
                if let Some(f) = &r.filter {
                    rep_filter = f.clone();
                }
                if let Some(i) = &r.info_raw {
                    rep_info = i.clone();
                }
                break;
            }
        }
    }

    // Per-sample FORMAT: fixed GT:AD:DP:GQ:PL.
    let fmt_col = "GT:AD:DP:GQ:PL";
    let mut sample_cols: Vec<String> = Vec::new();
    for (idx, sm) in per_sample.iter().enumerate() {
        let col = match sm.get(key) {
            Some(rows) => build_merged_sample_column(rows, &merged_alts),
            None => missing_sample_column(),
        };
        sample_cols.push(col);
        let _ = idx;
    }

    writeln!(
        w,
        "{chrom}\t{pos}\t{id}\t{ref_a}\t{alt}\t{qual}\t{filter}\t{info}\t{fmt}\t{sample}",
        chrom = key.0,
        pos = key.1,
        id = rep_id,
        ref_a = key.2,
        alt = merged_alts.join(","),
        qual = rep_qual,
        filter = rep_filter,
        info = rep_info,
        fmt = fmt_col,
        sample = sample_cols.join("\t"),
    )?;
    let _ = sample_names;
    Ok(())
}

/// Build a `GT:AD:DP:GQ:PL` column for one sample at one site, remapping
/// the sample's local allele indices into the merged ALT ordering.
fn build_merged_sample_column(rows: &[Row], merged_alts: &[String]) -> String {
    // Map local alt_index (1-based in VCF) -> merged allele index.
    // alt_index 0 is always REF, local indices >=1 correspond to the
    // sample's own sorted alt rows.
    let mut sorted: Vec<&Row> = rows.iter().collect();
    sorted.sort_by_key(|r| r.alt_index);
    let mut local_to_merged: Vec<u32> = vec![0]; // index 0 -> REF (0)
    for r in &sorted {
        // Find the ALT in merged_alts; should always be present.
        let merged_idx = merged_alts
            .iter()
            .position(|a| a == &r.alt_allele)
            .map(|p| (p + 1) as u32)
            .unwrap_or(0);
        local_to_merged.push(merged_idx);
    }

    // Remap GT string from the head row (all split rows share the same raw GT).
    let head = &sorted[0];
    let gt = match &head.gt_raw {
        Some(g) => remap_gt(g, &local_to_merged),
        None => "./.".to_string(),
    };
    let ad = head.format_ad.clone().unwrap_or_else(|| ".".to_string());
    let dp = head
        .read_depth
        .map(|x| x.to_string())
        .unwrap_or_else(|| ".".to_string());
    let gq = head
        .format_gq
        .map(|x| x.to_string())
        .unwrap_or_else(|| ".".to_string());
    let pl = head.format_pl.clone().unwrap_or_else(|| ".".to_string());
    format!("{gt}:{ad}:{dp}:{gq}:{pl}")
}

fn missing_sample_column() -> String {
    "./.:.:.:.:.".to_string()
}

/// Remap allele indices in a GT string. Preserves phasing separators.
/// Non-numeric tokens (e.g. ".") pass through unchanged.
fn remap_gt(gt: &str, local_to_merged: &[u32]) -> String {
    let mut out = String::with_capacity(gt.len());
    let mut token = String::new();
    for ch in gt.chars() {
        if ch == '|' || ch == '/' {
            flush_allele(&mut token, &mut out, local_to_merged);
            out.push(ch);
        } else {
            token.push(ch);
        }
    }
    flush_allele(&mut token, &mut out, local_to_merged);
    out
}

fn flush_allele(token: &mut String, out: &mut String, local_to_merged: &[u32]) {
    if token.is_empty() {
        return;
    }
    if let Ok(n) = token.parse::<usize>() {
        let merged = local_to_merged.get(n).copied().unwrap_or(n as u32);
        out.push_str(&merged.to_string());
    } else {
        out.push_str(token);
    }
    token.clear();
}
