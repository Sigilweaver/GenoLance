use std::fs::File;
use std::io::{BufRead, BufReader};
use std::num::NonZero;
use std::path::Path;
use std::sync::Arc;
use std::thread;

use anyhow::{anyhow, Context, Result};
use arrow_array::{Float32Array, RecordBatch, StringArray, UInt32Array, UInt64Array};
use arrow_schema::SchemaRef;
use indicatif::{ProgressBar, ProgressStyle};
use lancedb::connection::Connection;
use lancedb::index::scalar::{BTreeIndexBuilder, BitmapIndexBuilder};
use lancedb::index::Index;
use lancedb::Table;
use noodles::bgzf;
use noodles::vcf::{
    self as vcf,
    variant::record::{
        samples::{series::Value as SampleValue, Sample as SampleTrait},
        AlternateBases, Filters as FiltersTrait, Ids as IdsTrait, Info as InfoTrait,
    },
};
use tokio::sync::mpsc;

use biolance_core::schema::{
    clinvar_schema, samples_schema, variant_schema, CLINVAR_TABLE, SAMPLES_TABLE, VARIANTS_TABLE,
};
use biolance_core::store::Store;

/// Number of rows to buffer before flushing a batch to the Lance table.
const BATCH_SIZE: usize = 10_000;

/// Ingest one or more VCF/BCF files into the BioLance store at `store_path`.
///
/// Files whose filename contains "clinvar" or that have no sample columns
/// are ingested into the `clinvar` annotation table. All others are ingested
/// into the `variants` table, one row per (sample, chrom, pos, ref, alt).
pub async fn run(store_path: &str, files: &[String], sample_override: Option<&str>) -> Result<()> {
    let store = Store::open(store_path).await?;

    for file in files {
        ingest_one(&store, file, sample_override).await?;
    }
    Ok(())
}

async fn ingest_one(store: &Store, path: &str, sample_override: Option<&str>) -> Result<()> {
    let filename = Path::new(path)
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or(path);
    println!("[ingest] {filename}");

    let mut reader = open_vcf_reader(path).with_context(|| format!("opening {path}"))?;
    let header = reader
        .read_header()
        .with_context(|| format!("reading VCF header of {path}"))?;

    let looks_like_clinvar = filename.to_lowercase().contains("clinvar");
    let has_samples = !header.sample_names().is_empty();

    if looks_like_clinvar || !has_samples {
        ingest_clinvar(&store.variants, &mut reader, &header).await
    } else {
        let sample_name = sample_override
            .map(str::to_owned)
            .or_else(|| header.sample_names().iter().next().cloned())
            .or_else(|| {
                Path::new(path)
                    .file_stem()
                    .and_then(|s| s.to_str())
                    .map(String::from)
            })
            .ok_or_else(|| anyhow!("could not determine sample name for {path}"))?;

        register_sample(&store.conn, &sample_name, path, &header).await?;
        ingest_variants(&store.variants, &mut reader, &header, &sample_name).await
    }
}

/// Open a VCF reader. For BGZF-compressed inputs (`.vcf.gz` / `.vcf.bgz`),
/// decompress in parallel via `bgzf::io::MultithreadedReader`; for plain
/// `.vcf` use a standard buffered file reader. The returned reader is
/// erased to `Box<dyn BufRead + Send>` so both paths share the same
/// downstream ingest code.
fn open_vcf_reader(path: &str) -> Result<vcf::io::Reader<Box<dyn BufRead + Send>>> {
    let is_bgzf = path.ends_with(".gz") || path.ends_with(".bgz") || path.ends_with(".bcf");
    let file = File::open(path).with_context(|| format!("open {path}"))?;
    if is_bgzf {
        let workers = thread::available_parallelism()
            .map(|n| n.get().min(4).max(1))
            .unwrap_or(2);
        let mt =
            bgzf::io::MultithreadedReader::with_worker_count(NonZero::new(workers).unwrap(), file);
        Ok(vcf::io::Reader::new(Box::new(mt) as Box<dyn BufRead + Send>))
    } else {
        Ok(vcf::io::Reader::new(
            Box::new(BufReader::new(file)) as Box<dyn BufRead + Send>
        ))
    }
}

/// Open the named table or create an empty one with the given schema.
async fn open_or_create(conn: &Connection, name: &str, schema: SchemaRef) -> Result<Table> {
    let existing = conn.table_names().execute().await?;
    if existing.iter().any(|n| n == name) {
        Ok(conn.open_table(name).execute().await?)
    } else {
        Ok(conn.create_empty_table(name, schema).execute().await?)
    }
}

// ---------------------------------------------------------------------------
// Sample registry
// ---------------------------------------------------------------------------

/// Serialize the parsed header back to VCF text via the noodles writer so
/// `biolance export` can reproduce it byte-for-byte modulo ordering.
fn serialize_header(header: &vcf::Header) -> Result<String> {
    let mut buf: Vec<u8> = Vec::new();
    {
        let mut writer = vcf::io::Writer::new(&mut buf);
        writer.write_header(header)?;
    }
    Ok(String::from_utf8(buf)?)
}

fn extract_reference(header: &vcf::Header) -> Option<String> {
    use noodles::vcf::header::record::value::Collection;
    // `##reference=...` lands in `other_records` under the "reference" key.
    let other = header.other_records();
    let entries = other.get("reference")?;
    match entries {
        Collection::Unstructured(list) => list.first().cloned(),
        Collection::Structured(_) => None,
    }
}

async fn register_sample(
    conn: &Connection,
    sample_name: &str,
    source_path: &str,
    header: &vcf::Header,
) -> Result<()> {
    let table = open_or_create(conn, SAMPLES_TABLE, samples_schema()).await?;

    let header_text = serialize_header(header)?;
    let reference = extract_reference(header);
    let ingested_at = chrono::Utc::now().to_rfc3339();

    let schema = samples_schema();
    let batch = RecordBatch::try_new(
        schema,
        vec![
            Arc::new(StringArray::from(vec![sample_name.to_string()])),
            Arc::new(StringArray::from(vec![Some(source_path.to_string())])),
            Arc::new(StringArray::from(vec![header_text])),
            Arc::new(StringArray::from(vec![Some(ingested_at)])),
            Arc::new(StringArray::from(vec![reference])),
        ],
    )?;
    table.add(vec![batch]).execute().await?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Variant ingest (samples VCF)
// ---------------------------------------------------------------------------

#[derive(Default)]
struct VariantBuffer {
    sample_name: Vec<String>,
    chrom: Vec<String>,
    pos: Vec<u64>,
    ids: Vec<Option<String>>,
    ref_allele: Vec<String>,
    alt_allele: Vec<String>,
    alt_index: Vec<u32>,
    alt_count: Vec<u32>,
    quality: Vec<Option<f32>>,
    filter: Vec<Option<String>>,
    genotype: Vec<Option<String>>,
    gt_raw: Vec<Option<String>>,
    read_depth: Vec<Option<u32>>,
    format_ad: Vec<Option<String>>,
    format_gq: Vec<Option<u32>>,
    format_pl: Vec<Option<String>>,
    allele_freq: Vec<Option<f32>>,
    format_extra: Vec<Option<String>>,
    format_key_order: Vec<Option<String>>,
    info_raw: Vec<Option<String>>,
}

impl VariantBuffer {
    fn len(&self) -> usize {
        self.pos.len()
    }

    fn into_batch(self) -> Result<RecordBatch> {
        let schema = variant_schema();
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(self.sample_name)),
                Arc::new(StringArray::from(self.chrom)),
                Arc::new(UInt64Array::from(self.pos)),
                Arc::new(StringArray::from(self.ids)),
                Arc::new(StringArray::from(self.ref_allele)),
                Arc::new(StringArray::from(self.alt_allele)),
                Arc::new(UInt32Array::from(self.alt_index)),
                Arc::new(UInt32Array::from(self.alt_count)),
                Arc::new(Float32Array::from(self.quality)),
                Arc::new(StringArray::from(self.filter)),
                Arc::new(StringArray::from(self.genotype)),
                Arc::new(StringArray::from(self.gt_raw)),
                Arc::new(UInt32Array::from(self.read_depth)),
                Arc::new(StringArray::from(self.format_ad)),
                Arc::new(UInt32Array::from(self.format_gq)),
                Arc::new(StringArray::from(self.format_pl)),
                Arc::new(Float32Array::from(self.allele_freq)),
                Arc::new(StringArray::from(self.format_extra)),
                Arc::new(StringArray::from(self.format_key_order)),
                Arc::new(StringArray::from(self.info_raw)),
            ],
        )?;
        Ok(batch)
    }
}

async fn ingest_variants<R: std::io::BufRead>(
    conn: &Connection,
    reader: &mut vcf::io::Reader<R>,
    header: &vcf::Header,
    sample_name: &str,
) -> Result<()> {
    let table = open_or_create(conn, VARIANTS_TABLE, variant_schema()).await?;
    let progress = make_progress();
    let mut buffer = VariantBuffer::default();
    let mut total: u64 = 0;

    // Pipeline: parsing/encoding happens on this task while another task
    // awaits `table.add`. A small bounded channel keeps memory bounded
    // while letting one batch encode in parallel with the previous flush.
    let (tx, mut rx) = mpsc::channel::<RecordBatch>(2);
    let writer_table = table.clone();
    let writer = tokio::spawn(async move {
        while let Some(batch) = rx.recv().await {
            writer_table.add(vec![batch]).execute().await?;
        }
        Ok::<_, anyhow::Error>(())
    });

    for result in reader.records() {
        let record = result?;
        append_variant_rows(&record, header, sample_name, &mut buffer)?;

        if buffer.len() >= BATCH_SIZE {
            let n = buffer.len() as u64;
            let batch = std::mem::take(&mut buffer).into_batch()?;
            if tx.send(batch).await.is_err() {
                break;
            }
            total += n;
            progress.set_position(total);
        }
    }

    if buffer.len() > 0 {
        let n = buffer.len() as u64;
        let batch = buffer.into_batch()?;
        tx.send(batch).await.ok();
        total += n;
    }
    drop(tx);
    writer
        .await
        .map_err(|e| anyhow!("writer task panicked: {e}"))??;

    let done_msg = format!("ingested {total} rows for sample {sample_name}");
    progress.finish_with_message(done_msg.clone());
    eprintln!("{done_msg}");

    ensure_variant_indices(&table).await?;
    Ok(())
}

/// Create `(chrom, pos, sample_name)` scalar indices on the variants table
/// if they don't already exist. Region queries and per-sample filters rely
/// on these to avoid full scans once the store grows past a few samples.
async fn ensure_variant_indices(table: &Table) -> Result<()> {
    let existing: Vec<String> = table
        .list_indices()
        .await?
        .into_iter()
        .map(|i| i.columns.join(","))
        .collect();

    let want: &[(&str, Index)] = &[
        ("pos", Index::BTree(BTreeIndexBuilder::default())),
        ("chrom", Index::Bitmap(BitmapIndexBuilder::default())),
        ("sample_name", Index::Bitmap(BitmapIndexBuilder::default())),
    ];
    for (col, idx) in want {
        if existing.iter().any(|c| c == col) {
            continue;
        }
        // Using default builder (.train(true), .replace(true)) to build over
        // all current rows. Later ingests will re-train but that's cheap.
        table
            .create_index(&[(*col).to_string()], idx.clone())
            .execute()
            .await?;
    }
    Ok(())
}

#[derive(Default)]
struct SampleFormatFields {
    gt_raw: Option<String>,
    dp: Option<u32>,
    ad: Option<String>,
    gq: Option<u32>,
    pl: Option<String>,
    extra: Option<String>,
    /// Original FORMAT key order, colon-joined (e.g. "GT:GQ:DP:AD:VAF:PL").
    key_order: Option<String>,
}

fn append_variant_rows(
    record: &vcf::Record,
    header: &vcf::Header,
    sample_name: &str,
    buf: &mut VariantBuffer,
) -> Result<()> {
    let chrom = record.reference_sequence_name().to_string();
    let pos = match record.variant_start() {
        Some(p) => usize::from(p?) as u64,
        None => return Ok(()), // skip records without a position (symbolic etc.)
    };
    let ref_allele = record.reference_bases().to_string();
    let quality = record.quality_score().transpose().ok().flatten();

    let filter = {
        let filters = record.filters();
        let mut parts = Vec::new();
        for result in filters.iter(header) {
            parts.push(result?.to_string());
        }
        if parts.is_empty() {
            None
        } else {
            Some(parts.join(";"))
        }
    };

    let ids_str = {
        let ids_view = record.ids();
        let mut v = Vec::new();
        for id in ids_view.iter() {
            v.push(id.to_string());
        }
        if v.is_empty() {
            None
        } else {
            Some(v.join(";"))
        }
    };

    let fmt = extract_sample_format_fields(record, header)?;
    let info_af_per_alt = parse_info_af(record, header)?;
    let info_raw = serialize_info(record, header)?;

    let gt_indices: Option<Vec<Option<usize>>> = parse_gt_indices(fmt.gt_raw.as_deref());

    // Count real ALTs (skip "." / empty).
    let alt_list: Vec<String> = record
        .alternate_bases()
        .iter()
        .filter_map(|r| r.ok().map(|s| s.to_string()))
        .filter(|s| s != "." && !s.is_empty())
        .collect();
    let alt_count = alt_list.len() as u32;
    if alt_count == 0 {
        return Ok(());
    }

    for (alt_index, alt) in alt_list.into_iter().enumerate() {
        // Per-allele genotype: report whether the reference sample carries this alt.
        let gt_for_alt = gt_indices
            .as_ref()
            .map(|calls| format_gt_for_alt(calls, alt_index));

        let af = info_af_per_alt
            .as_ref()
            .and_then(|v| v.get(alt_index).copied().flatten());

        buf.sample_name.push(sample_name.to_string());
        buf.chrom.push(chrom.clone());
        buf.pos.push(pos);
        buf.ids.push(ids_str.clone());
        buf.ref_allele.push(ref_allele.clone());
        buf.alt_allele.push(alt);
        buf.alt_index.push(alt_index as u32);
        buf.alt_count.push(alt_count);
        buf.quality.push(quality);
        buf.filter.push(filter.clone());
        buf.genotype.push(gt_for_alt);
        buf.gt_raw.push(fmt.gt_raw.clone());
        buf.read_depth.push(fmt.dp);
        buf.format_ad.push(fmt.ad.clone());
        buf.format_gq.push(fmt.gq);
        buf.format_pl.push(fmt.pl.clone());
        buf.allele_freq.push(af);
        buf.format_extra.push(fmt.extra.clone());
        buf.format_key_order.push(fmt.key_order.clone());
        buf.info_raw.push(info_raw.clone());
    }

    Ok(())
}

/// Extract GT (as a raw display string), DP, AD, GQ, PL from the first
/// sample column. Values are returned as owned strings/u32 so nothing
/// borrows from the transient `samples()` view.
fn extract_sample_format_fields(
    record: &vcf::Record,
    header: &vcf::Header,
) -> Result<SampleFormatFields> {
    let samples = record.samples();
    let Some(sample) = samples.get_index(0) else {
        return Ok(SampleFormatFields::default());
    };

    let gt_raw = match sample.get(header, "GT") {
        Some(res) => res?.as_ref().map(value_to_string),
        None => None,
    };
    let dp = match sample.get(header, "DP") {
        Some(res) => match res? {
            Some(SampleValue::Integer(n)) if n >= 0 => Some(n as u32),
            _ => None,
        },
        None => None,
    };
    let ad = match sample.get(header, "AD") {
        Some(res) => res?.as_ref().map(value_to_string),
        None => None,
    };
    let gq = match sample.get(header, "GQ") {
        Some(res) => match res? {
            Some(SampleValue::Integer(n)) if n >= 0 => Some(n as u32),
            _ => None,
        },
        None => None,
    };
    let pl = match sample.get(header, "PL") {
        Some(res) => res?.as_ref().map(value_to_string),
        None => None,
    };

    // Collect everything else into a compact "KEY=VAL;KEY=VAL" blob so
    // `biolance export` can re-emit FORMAT fields we don't model directly
    // (VAF, MIN_DP, MED_DP, …). Known keys are skipped to avoid duplication.
    // Also capture the full key order (all keys, in record order) for faithful
    // FORMAT column reconstruction on export.
    const KNOWN: &[&str] = &["GT", "AD", "DP", "GQ", "PL"];
    let mut extra_parts: Vec<String> = Vec::new();
    let mut all_keys: Vec<String> = Vec::new();
    for result in sample.iter(header) {
        let (key, value) = result?;
        all_keys.push(key.to_string());
        if KNOWN.contains(&key) {
            continue;
        }
        match value {
            None => extra_parts.push(key.to_string()),
            Some(v) => extra_parts.push(format!("{key}={}", value_to_string(&v))),
        }
    }
    let extra = if extra_parts.is_empty() {
        None
    } else {
        Some(extra_parts.join(";"))
    };
    let key_order = if all_keys.is_empty() {
        None
    } else {
        Some(all_keys.join(":"))
    };

    Ok(SampleFormatFields {
        gt_raw,
        dp,
        ad,
        gq,
        pl,
        extra,
        key_order,
    })
}

fn value_to_string(v: &SampleValue<'_>) -> String {
    use noodles::vcf::variant::record::samples::series::value::{
        genotype::Phasing, Array as SampleArray,
    };
    match v {
        SampleValue::Integer(n) => n.to_string(),
        SampleValue::Float(n) => n.to_string(),
        SampleValue::Character(c) => c.to_string(),
        SampleValue::String(s) => s.to_string(),
        SampleValue::Genotype(g) => {
            let mut out = String::new();
            let mut first = true;
            for call in g.iter() {
                match call {
                    Ok((allele, phasing)) => {
                        if !first {
                            out.push(match phasing {
                                Phasing::Phased => '|',
                                Phasing::Unphased => '/',
                            });
                        }
                        first = false;
                        match allele {
                            Some(i) => out.push_str(&i.to_string()),
                            None => out.push('.'),
                        }
                    }
                    Err(_) => return String::from("."),
                }
            }
            out
        }
        SampleValue::Array(arr) => {
            // Flatten array-valued FORMAT fields (AD, PL, …) into "a,b,c".
            let parts: Vec<String> = match arr {
                SampleArray::Integer(a) => a
                    .iter()
                    .map(|r| match r {
                        Ok(Some(n)) => n.to_string(),
                        _ => ".".to_string(),
                    })
                    .collect(),
                SampleArray::Float(a) => a
                    .iter()
                    .map(|r| match r {
                        Ok(Some(n)) => n.to_string(),
                        _ => ".".to_string(),
                    })
                    .collect(),
                SampleArray::Character(a) => a
                    .iter()
                    .map(|r| match r {
                        Ok(Some(c)) => c.to_string(),
                        _ => ".".to_string(),
                    })
                    .collect(),
                SampleArray::String(a) => a
                    .iter()
                    .map(|r| match r {
                        Ok(Some(s)) => s.into_owned(),
                        _ => ".".to_string(),
                    })
                    .collect(),
            };
            parts.join(",")
        }
    }
}

/// Parse a genotype string like "0/1", "1|1", "./." into per-position allele
/// indices. Returns None if parsing fails or the input is None.
fn parse_gt_indices(gt: Option<&str>) -> Option<Vec<Option<usize>>> {
    let gt = gt?;
    let parts: Vec<Option<usize>> = gt
        .split(|c| c == '/' || c == '|')
        .map(|s| s.parse::<usize>().ok())
        .collect();
    Some(parts)
}

/// For a given ALT index (0-based into the ALT list), encode the genotype
/// as a per-allele biallelic call relative to reference. For multi-allelic
/// sites, other ALTs are treated as reference for this split row.
fn format_gt_for_alt(calls: &[Option<usize>], alt_index: usize) -> String {
    let target = alt_index + 1; // VCF allele indices: 0=REF, 1..=ALTs
    let encoded: Vec<String> = calls
        .iter()
        .map(|c| match c {
            Some(i) if *i == target => "1".to_string(),
            Some(_) => "0".to_string(),
            None => ".".to_string(),
        })
        .collect();
    encoded.join("/")
}

/// Try to read INFO/AF (one value per ALT). Returns Some(...) if present.
fn parse_info_af(record: &vcf::Record, header: &vcf::Header) -> Result<Option<Vec<Option<f32>>>> {
    use noodles::vcf::variant::record::info::field::{
        value::Array as InfoArray, Value as InfoValue,
    };

    let info = record.info();
    let Some(result) = info.get(header, "AF") else {
        return Ok(None);
    };
    let Some(value) = result? else {
        return Ok(None);
    };
    match value {
        InfoValue::Float(f) => Ok(Some(vec![Some(f)])),
        InfoValue::Array(InfoArray::Float(arr)) => {
            let mut out = Vec::new();
            for r in arr.iter() {
                out.push(r?);
            }
            Ok(Some(out))
        }
        _ => Ok(None),
    }
}

/// Serialize the full INFO column back to VCF text ("K=V;K=V" or just "K"
/// for flag fields). Stored once per site but duplicated across split
/// per-alt rows so export only has to read one.
fn serialize_info(record: &vcf::Record, header: &vcf::Header) -> Result<Option<String>> {
    use noodles::vcf::variant::record::info::field::{
        value::Array as InfoArray, Value as InfoValue,
    };

    let info = record.info();
    if info.is_empty() {
        return Ok(None);
    }

    let mut parts: Vec<String> = Vec::new();
    for result in info.iter(header) {
        let (key, value) = result?;
        match value {
            None => parts.push(key.to_string()),
            Some(v) => {
                let rendered = match v {
                    InfoValue::Flag => {
                        parts.push(key.to_string());
                        continue;
                    }
                    InfoValue::Integer(n) => n.to_string(),
                    InfoValue::Float(n) => n.to_string(),
                    InfoValue::Character(c) => c.to_string(),
                    InfoValue::String(s) => s.into_owned(),
                    InfoValue::Array(arr) => match arr {
                        InfoArray::Integer(a) => join_opt_iter(
                            a.iter()
                                .map(|r| r.map(|o: Option<i32>| o.map(|n| n.to_string()))),
                        ),
                        InfoArray::Float(a) => join_opt_iter(
                            a.iter()
                                .map(|r| r.map(|o: Option<f32>| o.map(|n| n.to_string()))),
                        ),
                        InfoArray::Character(a) => join_opt_iter(
                            a.iter()
                                .map(|r| r.map(|o: Option<char>| o.map(|c| c.to_string()))),
                        ),
                        InfoArray::String(a) => join_opt_iter(a.iter().map(|r| {
                            r.map(|o| o.map(|s: std::borrow::Cow<'_, str>| s.into_owned()))
                        })),
                    },
                };
                parts.push(format!("{key}={rendered}"));
            }
        }
    }
    if parts.is_empty() {
        Ok(None)
    } else {
        Ok(Some(parts.join(";")))
    }
}

fn join_opt_iter<I>(it: I) -> String
where
    I: Iterator<Item = std::io::Result<Option<String>>>,
{
    let parts: Vec<String> = it
        .map(|r| match r {
            Ok(Some(s)) => s,
            _ => ".".to_string(),
        })
        .collect();
    parts.join(",")
}

// ---------------------------------------------------------------------------
// ClinVar ingest
// ---------------------------------------------------------------------------

#[derive(Default)]
struct ClinVarBuffer {
    chrom: Vec<String>,
    pos: Vec<u64>,
    ref_allele: Vec<String>,
    alt_allele: Vec<String>,
    variation_id: Vec<Option<String>>,
    gene_symbol: Vec<Option<String>>,
    clinical_significance: Vec<Option<String>>,
    review_status: Vec<Option<String>>,
    disease_name: Vec<Option<String>>,
}

impl ClinVarBuffer {
    fn len(&self) -> usize {
        self.pos.len()
    }

    fn into_batch(self) -> Result<RecordBatch> {
        let schema = clinvar_schema();
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(self.chrom)),
                Arc::new(UInt64Array::from(self.pos)),
                Arc::new(StringArray::from(self.ref_allele)),
                Arc::new(StringArray::from(self.alt_allele)),
                Arc::new(StringArray::from(self.variation_id)),
                Arc::new(StringArray::from(self.gene_symbol)),
                Arc::new(StringArray::from(self.clinical_significance)),
                Arc::new(StringArray::from(self.review_status)),
                Arc::new(StringArray::from(self.disease_name)),
            ],
        )?;
        Ok(batch)
    }
}

async fn ingest_clinvar<R: std::io::BufRead>(
    conn: &Connection,
    reader: &mut vcf::io::Reader<R>,
    header: &vcf::Header,
) -> Result<()> {
    let table = open_or_create(conn, CLINVAR_TABLE, clinvar_schema()).await?;
    let progress = make_progress();
    let mut buffer = ClinVarBuffer::default();
    let mut total: u64 = 0;

    let (tx, mut rx) = mpsc::channel::<RecordBatch>(2);
    let writer_table = table.clone();
    let writer = tokio::spawn(async move {
        while let Some(batch) = rx.recv().await {
            writer_table.add(vec![batch]).execute().await?;
        }
        Ok::<_, anyhow::Error>(())
    });

    for result in reader.records() {
        let record = result?;
        append_clinvar_rows(&record, header, &mut buffer)?;

        if buffer.len() >= BATCH_SIZE {
            let n = buffer.len() as u64;
            let batch = std::mem::take(&mut buffer).into_batch()?;
            if tx.send(batch).await.is_err() {
                break;
            }
            total += n;
            progress.set_position(total);
        }
    }

    if buffer.len() > 0 {
        let n = buffer.len() as u64;
        let batch = buffer.into_batch()?;
        tx.send(batch).await.ok();
        total += n;
    }
    drop(tx);
    writer
        .await
        .map_err(|e| anyhow!("writer task panicked: {e}"))??;

    let done_msg = format!("ingested {total} ClinVar rows");
    progress.finish_with_message(done_msg.clone());
    eprintln!("{done_msg}");
    Ok(())
}

fn append_clinvar_rows(
    record: &vcf::Record,
    header: &vcf::Header,
    buf: &mut ClinVarBuffer,
) -> Result<()> {
    let chrom = record.reference_sequence_name().to_string();
    let Some(pos_result) = record.variant_start() else {
        return Ok(());
    };
    let pos = usize::from(pos_result?) as u64;
    let ref_allele = record.reference_bases().to_string();

    let mut ids = Vec::new();
    let ids_view = record.ids();
    for id in ids_view.iter() {
        ids.push(id.to_string());
    }
    let variation_id = if ids.is_empty() {
        None
    } else {
        Some(ids.join(";"))
    };

    let gene_symbol = info_scalar_string(record, header, "GENEINFO")?.map(|s| {
        // GENEINFO is formatted "SYMBOL:ID|SYMBOL:ID" → keep the first symbol.
        s.split('|')
            .next()
            .unwrap_or(&s)
            .split(':')
            .next()
            .unwrap_or(&s)
            .to_string()
    });

    let clinical_significance = info_scalar_string(record, header, "CLNSIG")?;
    let review_status = info_scalar_string(record, header, "CLNREVSTAT")?;
    let disease_name = info_scalar_string(record, header, "CLNDN")?;

    for alt_result in record.alternate_bases().iter() {
        let alt = alt_result?.to_string();
        if alt == "." || alt.is_empty() {
            continue;
        }
        buf.chrom.push(chrom.clone());
        buf.pos.push(pos);
        buf.ref_allele.push(ref_allele.clone());
        buf.alt_allele.push(alt);
        buf.variation_id.push(variation_id.clone());
        buf.gene_symbol.push(gene_symbol.clone());
        buf.clinical_significance
            .push(clinical_significance.clone());
        buf.review_status.push(review_status.clone());
        buf.disease_name.push(disease_name.clone());
    }

    Ok(())
}

fn info_scalar_string(
    record: &vcf::Record,
    header: &vcf::Header,
    key: &str,
) -> Result<Option<String>> {
    use noodles::vcf::variant::record::info::field::{
        value::Array as InfoArray, Value as InfoValue,
    };

    let info = record.info();
    let Some(result) = info.get(header, key) else {
        return Ok(None);
    };
    let Some(value) = result? else {
        return Ok(None);
    };
    Ok(Some(match value {
        InfoValue::String(s) => s.into_owned(),
        InfoValue::Integer(n) => n.to_string(),
        InfoValue::Float(f) => f.to_string(),
        InfoValue::Flag => "true".to_string(),
        InfoValue::Character(c) => c.to_string(),
        InfoValue::Array(arr) => match arr {
            InfoArray::String(a) => {
                let mut parts = Vec::new();
                for r in a.iter() {
                    if let Some(s) = r? {
                        parts.push(s.into_owned());
                    }
                }
                parts.join(",")
            }
            InfoArray::Integer(a) => {
                let mut parts = Vec::new();
                for r in a.iter() {
                    if let Some(n) = r? {
                        parts.push(n.to_string());
                    }
                }
                parts.join(",")
            }
            InfoArray::Float(a) => {
                let mut parts = Vec::new();
                for r in a.iter() {
                    if let Some(n) = r? {
                        parts.push(n.to_string());
                    }
                }
                parts.join(",")
            }
            InfoArray::Character(a) => {
                let mut parts = Vec::new();
                for r in a.iter() {
                    if let Some(c) = r? {
                        parts.push(c.to_string());
                    }
                }
                parts.join(",")
            }
        },
    }))
}

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

fn make_progress() -> ProgressBar {
    let pb = ProgressBar::new_spinner();
    pb.set_style(
        ProgressStyle::with_template("{spinner} {pos} rows {msg}")
            .unwrap_or_else(|_| ProgressStyle::default_spinner()),
    );
    pb.enable_steady_tick(std::time::Duration::from_millis(200));
    pb
}
