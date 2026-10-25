use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;
use std::sync::atomic::{AtomicU16, Ordering};
use std::sync::Arc;

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
    header::record::value::map::info::Number as InfoNumber,
    variant::record::{
        samples::series::Value as SampleValue, AlternateBases, Filters as FiltersTrait,
        Ids as IdsTrait,
    },
};
use tokio::sync::mpsc;

use genolance_core::schema::{
    clinvar_schema, samples_schema, variant_schema, CLINVAR_TABLE, SAMPLES_TABLE, VARIANTS_TABLE,
};
use genolance_core::store::Store;
use genolance_core::Error;

/// Number of rows to buffer before flushing a batch to the Lance table.
const BATCH_SIZE: usize = 10_000;

/// The 28-byte BGZF end-of-file marker block (an empty BGZF block).
/// Source: SAM/BAM specification (hts-specs, SAMv1.pdf), section 4.1.2
/// "End-of-file marker".
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// Options for [`run`]. Build with [`IngestOptions::new`] and the
/// `with_*` methods; new options can be added without breaking callers.
#[derive(Debug, Clone, Default)]
#[non_exhaustive]
pub struct IngestOptions {
    /// Name to store the sample under instead of the VCF header's name.
    /// Only valid when ingesting a single file.
    pub sample_override: Option<String>,
    /// Which sample column to ingest from a multi-sample VCF.
    pub sample_column: Option<String>,
    /// Replace a sample (or the ClinVar release) that is already stored.
    pub replace: bool,
    /// Ingest as a ClinVar release even if the header has no INFO/CLNSIG.
    pub clinvar: bool,
}

impl IngestOptions {
    /// Strict defaults: no renaming, one-sample VCFs only, no replacing.
    pub fn new() -> Self {
        Self::default()
    }
    pub fn with_sample_override(mut self, name: impl Into<String>) -> Self {
        self.sample_override = Some(name.into());
        self
    }
    pub fn with_sample_column(mut self, name: impl Into<String>) -> Self {
        self.sample_column = Some(name.into());
        self
    }
    pub fn with_replace(mut self, replace: bool) -> Self {
        self.replace = replace;
        self
    }
    pub fn with_clinvar(mut self, clinvar: bool) -> Self {
        self.clinvar = clinvar;
        self
    }
}

/// What an ingest run did. Counts of everything it did not store are here,
/// so a skipped record is never invisible.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
#[non_exhaustive]
pub struct IngestReport {
    /// Files fully ingested.
    pub files: usize,
    /// Rows written (one per sample and ALT allele, or per ClinVar ALT).
    pub rows: u64,
    /// Records skipped because POS is 0.
    pub skipped_no_position: u64,
    /// Records skipped because ALT is `.` (reference-only sites).
    pub skipped_no_alt: u64,
}

impl IngestReport {
    fn add(&mut self, other: &IngestReport) {
        self.files += other.files;
        self.rows += other.rows;
        self.skipped_no_position += other.skipped_no_position;
        self.skipped_no_alt += other.skipped_no_alt;
    }
}

/// Ingest one or more VCF files into the GenoLance store at `store_path`.
///
/// A file with a ClinVar header (INFO/CLNSIG), or any file when
/// `opts.clinvar` is set, is ingested into the `clinvar` annotation table.
/// Other files must have sample columns and are ingested into the
/// `variants` table, one row per (sample, chrom, pos, ref, alt).
///
/// Each file is all-or-nothing: if it fails part-way, the rows it wrote are
/// removed before the error is returned.
pub async fn run(store_path: &str, files: &[String], opts: &IngestOptions) -> Result<IngestReport> {
    if opts.sample_override.is_some() && files.len() > 1 {
        return Err(Error::InvalidArgument(format!(
            "--sample renames one sample, but {} files were given; ingest the files one at a time",
            files.len()
        ))
        .into());
    }
    let store = Store::open(store_path).await?;

    let mut report = IngestReport::default();
    for file in files {
        report.add(&ingest_one(&store, file, opts).await?);
    }
    Ok(report)
}

async fn ingest_one(store: &Store, path: &str, opts: &IngestOptions) -> Result<IngestReport> {
    let filename = Path::new(path)
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or(path);
    println!("[ingest] {filename}");

    let (mut reader, tail) = open_vcf_reader(path).with_context(|| format!("opening {path}"))?;
    let header = reader
        .read_header()
        .with_context(|| format!("reading VCF header of {path}"))?;

    let n_samples = header.sample_names().len();
    let is_clinvar = if opts.clinvar {
        if n_samples > 0 {
            return Err(Error::InvalidArgument(format!(
                "{path}: --clinvar was given but the VCF has {n_samples} sample column(s); a ClinVar release has none"
            ))
            .into());
        }
        true
    } else if n_samples == 0 {
        if header.infos().contains_key("CLNSIG") {
            true
        } else {
            return Err(Error::UnrecognizedInput {
                path: path.to_string(),
                reason: "VCF has no sample columns and no INFO/CLNSIG definition, so it is not recognised as ClinVar; pass --clinvar to ingest it as a ClinVar release".into(),
            }
            .into());
        }
    } else {
        false
    };

    if is_clinvar {
        let table = open_or_create(&store.variants, CLINVAR_TABLE, clinvar_schema()).await?;
        prepare_clinvar_table(&table, opts.replace).await?;
        let result = ingest_clinvar(&table, &mut reader, &header, path)
            .await
            .and_then(|r| tail.check_complete(path).map(|()| r));
        match result {
            Ok(report) => Ok(report),
            Err(e) => Err(rollback(delete_all(&table).await, e, "clinvar")),
        }
    } else {
        let col_idx = select_sample_column(path, &header, opts)?;
        let file_sample_name = header
            .sample_names()
            .get_index(col_idx)
            .cloned()
            .ok_or_else(|| anyhow!("{path}: sample column {col_idx} not in header"))?;
        let sample_name = opts.sample_override.clone().unwrap_or(file_sample_name);

        let table = open_or_create(&store.variants, VARIANTS_TABLE, variant_schema()).await?;
        prepare_sample(store, &table, &sample_name, opts.replace).await?;

        let result = ingest_variants(&table, &mut reader, &header, &sample_name, col_idx, path)
            .await
            .and_then(|r| tail.check_complete(path).map(|()| r));
        let result = match result {
            Ok(r) => register_sample(&store.conn, &sample_name, path, &header)
                .await
                .map(|()| r),
            Err(e) => Err(e),
        };
        let report = match result {
            Ok(r) => r,
            Err(e) => {
                let pred = format!("sample_name = '{}'", sql_escape(&sample_name));
                let cleanup = table
                    .delete(&pred)
                    .await
                    .map(|_| ())
                    .map_err(anyhow::Error::from);
                return Err(rollback(cleanup, e, &sample_name));
            }
        };
        ensure_variant_indices(&table).await?;
        Ok(report)
    }
}

/// Combine an ingest error with the outcome of removing its partial rows.
/// A failed cleanup is reported, never dropped.
fn rollback(cleanup: Result<()>, err: anyhow::Error, what: &str) -> anyhow::Error {
    match cleanup {
        Ok(()) => err.context(format!(
            "ingest of {what} failed; the rows it wrote were removed"
        )),
        Err(c) => err.context(format!(
            "ingest of {what} failed AND removing its partial rows also failed ({c:#}); the store may hold partial data for {what}"
        )),
    }
}

async fn delete_all(table: &Table) -> Result<()> {
    table.delete("pos >= 0").await?;
    Ok(())
}

fn sql_escape(s: &str) -> String {
    s.replace('\'', "''")
}

/// Choose the sample column to ingest. A multi-sample VCF needs an explicit
/// `--sample-column`; nothing is dropped silently.
fn select_sample_column(path: &str, header: &vcf::Header, opts: &IngestOptions) -> Result<usize> {
    let names = header.sample_names();
    match &opts.sample_column {
        Some(want) => names.get_index_of(want.as_str()).ok_or_else(|| {
            Error::UnknownSample {
                name: want.clone(),
                available: names.iter().cloned().collect(),
            }
            .into()
        }),
        None if names.len() == 1 => Ok(0),
        None => Err(Error::MultiSample {
            path: path.to_string(),
            samples: names.iter().cloned().collect(),
        }
        .into()),
    }
}

async fn table_exists(conn: &Connection, name: &str) -> Result<bool> {
    Ok(conn
        .table_names()
        .execute()
        .await?
        .iter()
        .any(|n| n == name))
}

/// Refuse to add a sample that is already stored (or that has rows left by
/// an interrupted ingest). With `replace`, remove the old rows first.
async fn prepare_sample(
    store: &Store,
    calls: &Table,
    sample_name: &str,
    replace: bool,
) -> Result<()> {
    let pred = format!("sample_name = '{}'", sql_escape(sample_name));
    let mut exists = calls.count_rows(Some(pred.clone())).await? > 0;
    let have_registry = table_exists(&store.conn, SAMPLES_TABLE).await?;
    let registry = if have_registry {
        Some(store.conn.open_table(SAMPLES_TABLE).execute().await?)
    } else {
        None
    };
    if let Some(reg) = &registry {
        exists |= reg.count_rows(Some(pred.clone())).await? > 0;
    }
    if !exists {
        return Ok(());
    }
    if !replace {
        return Err(Error::DuplicateSample {
            sample: sample_name.to_string(),
        }
        .into());
    }
    calls.delete(&pred).await?;
    if let Some(reg) = registry {
        reg.delete(&pred).await?;
    }
    Ok(())
}

async fn prepare_clinvar_table(table: &Table, replace: bool) -> Result<()> {
    let existing = table.count_rows(None).await?;
    if existing == 0 {
        return Ok(());
    }
    if !replace {
        return Err(Error::ClinVarReleasePresent { existing }.into());
    }
    delete_all(table).await
}

/// Remembers the last byte handed to the VCF parser, so a file that stops
/// mid-line is detected after parsing.
struct TailCheck(Arc<AtomicU16>);

const NO_BYTE: u16 = 256;

impl TailCheck {
    fn check_complete(&self, path: &str) -> Result<()> {
        match self.0.load(Ordering::SeqCst) {
            b if b == u16::from(b'\n') => Ok(()),
            _ => Err(Error::Truncated {
                path: path.to_string(),
                detail: "the last line has no terminating newline".into(),
            }
            .into()),
        }
    }
}

struct TrackingReader<R> {
    inner: R,
    last: Arc<AtomicU16>,
}

impl<R: Read> Read for TrackingReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = self.inner.read(buf)?;
        if n > 0 {
            self.last.store(u16::from(buf[n - 1]), Ordering::SeqCst);
        }
        Ok(n)
    }
}

impl<R: BufRead> BufRead for TrackingReader<R> {
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        self.inner.fill_buf()
    }
    fn consume(&mut self, amt: usize) {
        if amt > 0 {
            if let Ok(buf) = self.inner.fill_buf() {
                if amt <= buf.len() {
                    self.last.store(u16::from(buf[amt - 1]), Ordering::SeqCst);
                }
            }
        }
        self.inner.consume(amt);
    }
}

/// Open a VCF reader. For BGZF-compressed inputs (`.vcf.gz` / `.vcf.bgz`),
/// decompress in parallel via `bgzf::io::MultithreadedReader`; for plain
/// `.vcf` use a standard buffered file reader. The returned reader is
/// erased to `Box<dyn BufRead + Send>` so both paths share the same
/// downstream ingest code.
///
/// BGZF inputs must end with the BGZF EOF marker block; a file without it
/// was cut at a block boundary and is rejected. The returned [`TailCheck`]
/// reports, once parsing is done, whether the decompressed stream ended on
/// a newline.
///
/// As of noodles-bgzf 0.48, `MultithreadedReader` decompresses on rayon's
/// global thread pool rather than a caller-specified worker count; configure
/// pool size via `rayon::ThreadPoolBuilder` at the process level if needed.
fn open_vcf_reader(path: &str) -> Result<(vcf::io::Reader<Box<dyn BufRead + Send>>, TailCheck)> {
    let is_bgzf = path.ends_with(".gz") || path.ends_with(".bgz") || path.ends_with(".bcf");
    let mut file = File::open(path).with_context(|| format!("open {path}"))?;
    let last = Arc::new(AtomicU16::new(NO_BYTE));
    let inner: Box<dyn BufRead + Send> = if is_bgzf {
        let len = file.metadata()?.len();
        let mut tail = [0u8; 28];
        let has_eof = len >= 28 && {
            file.seek(SeekFrom::End(-28))?;
            file.read_exact(&mut tail)?;
            file.seek(SeekFrom::Start(0))?;
            tail == BGZF_EOF
        };
        if !has_eof {
            return Err(Error::Truncated {
                path: path.to_string(),
                detail: "the BGZF end-of-file marker block is missing (the writer did not finish)"
                    .into(),
            }
            .into());
        }
        Box::new(bgzf::io::MultithreadedReader::new(file))
    } else {
        Box::new(BufReader::new(file))
    };
    let tracked = TrackingReader {
        inner,
        last: Arc::clone(&last),
    };
    Ok((
        vcf::io::Reader::new(Box::new(tracked) as Box<dyn BufRead + Send>),
        TailCheck(last),
    ))
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
/// `genolance export` can reproduce it byte-for-byte modulo ordering.
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

/// Record the sample in the registry. Called only after its variants were
/// written and the file was checked complete, so a registry row means a
/// finished ingest.
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

/// Records the ingester could not turn into rows, counted by reason and
/// printed in the summary.
#[derive(Default)]
struct SkipCounts {
    /// POS 0 (no position; e.g. telomere markers).
    no_position: u64,
    /// ALT is `.`: reference-only sites such as hom-ref blocks.
    no_alt: u64,
}

impl SkipCounts {
    fn report(&self) {
        if self.no_position > 0 {
            eprintln!(
                "skipped {} record(s) with POS 0 (no position)",
                self.no_position
            );
        }
        if self.no_alt > 0 {
            eprintln!(
                "skipped {} record(s) with ALT '.' (reference-only sites are not stored)",
                self.no_alt
            );
        }
    }
}

async fn ingest_variants<R: std::io::BufRead>(
    table: &Table,
    reader: &mut vcf::io::Reader<R>,
    header: &vcf::Header,
    sample_name: &str,
    sample_col: usize,
    path: &str,
) -> Result<IngestReport> {
    let progress = make_progress();
    let mut buffer = VariantBuffer::default();
    let mut total: u64 = 0;
    let mut skipped = SkipCounts::default();
    let expected_samples = header.sample_names().len();

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

    let mut parse_err: Option<anyhow::Error> = None;
    for (n, result) in reader.records().enumerate() {
        let step = result
            .with_context(|| format!("{path}: reading record {}", n + 1))
            .and_then(|record| {
                append_variant_rows(
                    &record,
                    header,
                    sample_name,
                    sample_col,
                    expected_samples,
                    &mut buffer,
                    &mut skipped,
                )
                .with_context(|| format!("{path}: record {}", n + 1))
            });
        if let Err(e) = step {
            parse_err = Some(e);
            break;
        }

        if buffer.len() >= BATCH_SIZE {
            let rows = buffer.len() as u64;
            let batch = std::mem::take(&mut buffer).into_batch()?;
            if tx.send(batch).await.is_err() {
                break;
            }
            total += rows;
            progress.set_position(total);
        }
    }

    if parse_err.is_none() && buffer.len() > 0 {
        let rows = buffer.len() as u64;
        let batch = buffer.into_batch()?;
        if tx.send(batch).await.is_ok() {
            total += rows;
        }
    }
    drop(tx);
    let written = writer
        .await
        .map_err(|e| anyhow!("writer task panicked: {e}"))?;
    if let Some(e) = parse_err {
        // Report the parse error; a concurrent writer error would only be a
        // secondary symptom, but do not drop it.
        return Err(match written {
            Ok(()) => e,
            Err(w) => e.context(format!("(the table writer also failed: {w:#})")),
        });
    }
    written?;

    let done_msg = format!("ingested {total} rows for sample {sample_name}");
    progress.finish_with_message(done_msg.clone());
    eprintln!("{done_msg}");
    skipped.report();
    Ok(IngestReport {
        files: 1,
        rows: total,
        skipped_no_position: skipped.no_position,
        skipped_no_alt: skipped.no_alt,
    })
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

fn invalid_format(chrom: &str, pos: u64, key: &str, reason: impl ToString) -> anyhow::Error {
    Error::MalformedRecord {
        chrom: chrom.to_string(),
        pos,
        field: format!("FORMAT/{key}"),
        reason: reason.to_string(),
    }
    .into()
}

fn append_variant_rows(
    record: &vcf::Record,
    header: &vcf::Header,
    sample_name: &str,
    sample_col: usize,
    expected_samples: usize,
    buf: &mut VariantBuffer,
    skipped: &mut SkipCounts,
) -> Result<()> {
    let chrom = record.reference_sequence_name().to_string();
    let pos = match record.variant_start() {
        Some(p) => usize::from(p.context("invalid POS")?) as u64,
        None => {
            skipped.no_position += 1;
            return Ok(());
        }
    };

    let found_samples = record.samples().iter().count();
    if found_samples != expected_samples {
        return Err(Error::MalformedRecord {
            chrom,
            pos,
            field: "sample columns".into(),
            reason: format!(
                "record has {found_samples} sample column(s) but the header declares {expected_samples}"
            ),
        }
        .into());
    }

    let ref_allele = record.reference_bases().to_string();
    let quality = record
        .quality_score()
        .transpose()
        .with_context(|| format!("{chrom}:{pos}: invalid QUAL"))?;

    let filter = {
        let filters = record.filters();
        let mut parts = Vec::new();
        for result in filters.iter(header) {
            parts.push(
                result
                    .with_context(|| format!("{chrom}:{pos}: invalid FILTER"))?
                    .to_string(),
            );
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

    // ALT list. A lone "." means a reference-only site. Missing or empty
    // elements next to real alleles would shift the allele numbering, so
    // they are an error.
    let mut alt_list: Vec<String> = Vec::new();
    let mut missing_alts = 0usize;
    for r in record.alternate_bases().iter() {
        let s = r
            .with_context(|| format!("{chrom}:{pos}: invalid ALT"))?
            .to_string();
        if s == "." || s.is_empty() {
            missing_alts += 1;
        } else {
            alt_list.push(s);
        }
    }
    if missing_alts > 0 && !alt_list.is_empty() {
        return Err(Error::MalformedRecord {
            chrom,
            pos,
            field: "ALT".into(),
            reason: format!(
                "'{}' mixes missing or empty alleles with real alleles",
                record.alternate_bases().as_ref()
            ),
        }
        .into());
    }

    let fmt = extract_sample_format_fields(record, header, sample_col, &chrom, pos)?;
    let alt_count = alt_list.len();
    let info_raw = validate_and_raw_info(record, header, &chrom, pos)?;
    if alt_count == 0 {
        skipped.no_alt += 1;
        return Ok(());
    }
    let info_af_per_alt = parse_info_af(record, header, alt_count, &chrom, pos)?;
    let gt_indices = match fmt.gt_raw.as_deref() {
        Some(gt) => Some(parse_gt_indices(gt, alt_count, &chrom, pos)?),
        None => None,
    };
    let alt_count = alt_count as u32;

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

/// Escape a raw FORMAT value for storage in the `KEY=VAL;KEY=VAL`
/// `format_extra` column: `%`, `;` and `=` become `%25`, `%3B`, `%3D`.
pub(crate) fn encode_extra_value(raw: &str) -> String {
    let mut out = String::with_capacity(raw.len());
    for c in raw.chars() {
        match c {
            '%' => out.push_str("%25"),
            ';' => out.push_str("%3B"),
            '=' => out.push_str("%3D"),
            c => out.push(c),
        }
    }
    out
}

/// Inverse of [`encode_extra_value`].
pub(crate) fn decode_extra_value(stored: &str) -> String {
    stored
        .replace("%3B", ";")
        .replace("%3D", "=")
        .replace("%25", "%")
}

/// Extract GT, DP, AD, GQ, PL and the remaining FORMAT fields from the
/// chosen sample column. Values are kept as the raw text of the file so
/// that export writes back exactly what was read; the typed parse is only
/// used to reject malformed values.
fn extract_sample_format_fields(
    record: &vcf::Record,
    header: &vcf::Header,
    sample_col: usize,
    chrom: &str,
    pos: u64,
) -> Result<SampleFormatFields> {
    let samples = record.samples();
    let Some(sample) = samples.get_index(sample_col) else {
        return Ok(SampleFormatFields::default());
    };
    let raw_sample: &str = sample.as_ref();

    let mut out = SampleFormatFields::default();
    let mut extra_parts: Vec<String> = Vec::new();
    let mut all_keys: Vec<String> = Vec::new();
    let mut raw_values = raw_sample.split(':');
    for result in sample.iter(header) {
        let (key, value) = result.map_err(|e| invalid_format(chrom, pos, "*", e))?;
        let raw = raw_values.next().unwrap_or(".");
        all_keys.push(key.to_string());
        if let Some(v) = &value {
            check_sample_value(v).map_err(|e| invalid_format(chrom, pos, key, e))?;
        }
        let missing = raw == ".";
        match key {
            "GT" => {
                if !missing {
                    out.gt_raw = Some(raw.to_string());
                }
            }
            "AD" => {
                if !missing {
                    out.ad = Some(raw.to_string());
                }
            }
            "PL" => {
                if !missing {
                    out.pl = Some(raw.to_string());
                }
            }
            "DP" | "GQ" => {
                let slot = if key == "DP" {
                    &mut out.dp
                } else {
                    &mut out.gq
                };
                match &value {
                    None => {}
                    Some(SampleValue::Integer(n)) if *n >= 0 => *slot = Some(*n as u32),
                    // Not a non-negative integer (a Float-typed GQ, a
                    // negative sentinel): keep the text as written.
                    Some(_) => extra_parts.push(format!("{key}={}", encode_extra_value(raw))),
                }
            }
            _ => extra_parts.push(format!("{key}={}", encode_extra_value(raw))),
        }
    }
    out.extra = if extra_parts.is_empty() {
        None
    } else {
        Some(extra_parts.join(";"))
    };
    out.key_order = if all_keys.is_empty() {
        None
    } else {
        Some(all_keys.join(":"))
    };
    Ok(out)
}

/// Walk every element of an array or genotype so lazily-parsed errors
/// surface here instead of being turned into missing values.
fn check_sample_value(v: &SampleValue<'_>) -> std::io::Result<()> {
    use noodles::vcf::variant::record::samples::series::value::Array as SampleArray;
    match v {
        SampleValue::Genotype(g) => {
            for call in g.iter() {
                call?;
            }
        }
        SampleValue::Array(arr) => match arr {
            SampleArray::Integer(a) => {
                for r in a.iter() {
                    r?;
                }
            }
            SampleArray::Float(a) => {
                for r in a.iter() {
                    r?;
                }
            }
            SampleArray::Character(a) => {
                for r in a.iter() {
                    r?;
                }
            }
            SampleArray::String(a) => {
                for r in a.iter() {
                    r?;
                }
            }
        },
        _ => {}
    }
    Ok(())
}

/// Parse a genotype string like "0/1", "1|1", "./." into per-position allele
/// indices. Every allele must be `.` or an index within the ALT list.
fn parse_gt_indices(
    gt: &str,
    alt_count: usize,
    chrom: &str,
    pos: u64,
) -> Result<Vec<Option<usize>>> {
    let mut parts = Vec::new();
    for tok in gt.split(['/', '|']) {
        if tok == "." {
            parts.push(None);
            continue;
        }
        let n: usize = tok.parse().map_err(|_| Error::MalformedRecord {
            chrom: chrom.to_string(),
            pos,
            field: "FORMAT/GT".into(),
            reason: format!("malformed genotype '{gt}'"),
        })?;
        if n > alt_count {
            return Err(Error::MalformedRecord {
                chrom: chrom.to_string(),
                pos,
                field: "FORMAT/GT".into(),
                reason: format!(
                    "genotype '{gt}' refers to allele {n} but the record has {alt_count} ALT allele(s)"
                ),
            }
            .into());
        }
        parts.push(Some(n));
    }
    Ok(parts)
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

/// Read INFO/AF as one value per ALT allele, using the Number declared in
/// the header: `A` gives one value per ALT, `R` gives REF then one per ALT.
/// Any other layout, a non-Float type, or a wrong value count is an error.
fn parse_info_af(
    record: &vcf::Record,
    header: &vcf::Header,
    alt_count: usize,
    chrom: &str,
    pos: u64,
) -> Result<Option<Vec<Option<f32>>>> {
    use noodles::vcf::header::record::value::map::info::Type as InfoType;
    use noodles::vcf::variant::record::info::field::{
        value::Array as InfoArray, Value as InfoValue,
    };
    let af_err = |reason: String| -> anyhow::Error {
        Error::MalformedRecord {
            chrom: chrom.to_string(),
            pos,
            field: "INFO/AF".into(),
            reason,
        }
        .into()
    };

    let info = record.info();
    let Some(result) = info.get(header, "AF") else {
        return Ok(None);
    };
    let Some(value) = result.map_err(|e| af_err(e.to_string()))? else {
        return Ok(None);
    };
    let (number, ty) = match header.infos().get("AF") {
        Some(def) => (def.number(), def.ty()),
        None => (InfoNumber::A, InfoType::Float),
    };
    if ty != InfoType::Float {
        return Err(af_err(format!(
            "header declares AF as Type={ty:?}, not Float"
        )));
    }
    let values: Vec<Option<f32>> = match value {
        InfoValue::Float(f) => vec![Some(f)],
        InfoValue::Array(InfoArray::Float(arr)) => {
            let mut out = Vec::new();
            for r in arr.iter() {
                out.push(r.map_err(|e| af_err(e.to_string()))?);
            }
            out
        }
        other => return Err(af_err(format!("expected Float values, found {other:?}"))),
    };
    let per_alt = match number {
        InfoNumber::A => values,
        InfoNumber::Count(1) if alt_count == 1 => values,
        InfoNumber::R => values.get(1..).map(<[_]>::to_vec).unwrap_or_default(),
        other => {
            return Err(af_err(format!(
                "header declares AF with Number={other:?}; only A and R are understood"
            )))
        }
    };
    if per_alt.len() != alt_count {
        return Err(af_err(format!(
            "{} value(s) for {alt_count} ALT allele(s)",
            per_alt.len()
        )));
    }
    Ok(Some(per_alt))
}

/// Check that every INFO value parses (including each array element), and
/// return the INFO column exactly as written. `.` (no INFO) becomes None.
fn validate_and_raw_info(
    record: &vcf::Record,
    header: &vcf::Header,
    chrom: &str,
    pos: u64,
) -> Result<Option<String>> {
    use noodles::vcf::variant::record::info::field::{
        value::Array as InfoArray, Value as InfoValue,
    };
    let info = record.info();
    for result in info.iter(header) {
        let invalid = |key: &str, reason: String| -> anyhow::Error {
            Error::MalformedRecord {
                chrom: chrom.to_string(),
                pos,
                field: format!("INFO/{key}"),
                reason,
            }
            .into()
        };
        let (key, value) = result.map_err(|e| invalid("*", e.to_string()))?;
        if let Some(InfoValue::Array(arr)) = value {
            let res: std::io::Result<()> = (|| {
                match arr {
                    InfoArray::Integer(a) => {
                        for r in a.iter() {
                            r?;
                        }
                    }
                    InfoArray::Float(a) => {
                        for r in a.iter() {
                            r?;
                        }
                    }
                    InfoArray::Character(a) => {
                        for r in a.iter() {
                            r?;
                        }
                    }
                    InfoArray::String(a) => {
                        for r in a.iter() {
                            r?;
                        }
                    }
                }
                Ok(())
            })();
            res.map_err(|e| invalid(key, e.to_string()))?;
        }
    }
    let raw: &str = info.as_ref();
    Ok(if raw.is_empty() || raw == "." {
        None
    } else {
        Some(raw.to_string())
    })
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
    table: &Table,
    reader: &mut vcf::io::Reader<R>,
    header: &vcf::Header,
    path: &str,
) -> Result<IngestReport> {
    let progress = make_progress();
    let mut buffer = ClinVarBuffer::default();
    let mut total: u64 = 0;
    let mut skipped = SkipCounts::default();

    let (tx, mut rx) = mpsc::channel::<RecordBatch>(2);
    let writer_table = table.clone();
    let writer = tokio::spawn(async move {
        while let Some(batch) = rx.recv().await {
            writer_table.add(vec![batch]).execute().await?;
        }
        Ok::<_, anyhow::Error>(())
    });

    let mut parse_err: Option<anyhow::Error> = None;
    for (n, result) in reader.records().enumerate() {
        let step = result
            .with_context(|| format!("{path}: reading record {}", n + 1))
            .and_then(|record| {
                append_clinvar_rows(&record, header, &mut buffer, &mut skipped)
                    .with_context(|| format!("{path}: record {}", n + 1))
            });
        if let Err(e) = step {
            parse_err = Some(e);
            break;
        }

        if buffer.len() >= BATCH_SIZE {
            let rows = buffer.len() as u64;
            let batch = std::mem::take(&mut buffer).into_batch()?;
            if tx.send(batch).await.is_err() {
                break;
            }
            total += rows;
            progress.set_position(total);
        }
    }

    if parse_err.is_none() && buffer.len() > 0 {
        let rows = buffer.len() as u64;
        let batch = buffer.into_batch()?;
        if tx.send(batch).await.is_ok() {
            total += rows;
        }
    }
    drop(tx);
    let written = writer
        .await
        .map_err(|e| anyhow!("writer task panicked: {e}"))?;
    if let Some(e) = parse_err {
        return Err(match written {
            Ok(()) => e,
            Err(w) => e.context(format!("(the table writer also failed: {w:#})")),
        });
    }
    written?;

    let done_msg = format!("ingested {total} ClinVar rows");
    progress.finish_with_message(done_msg.clone());
    eprintln!("{done_msg}");
    skipped.report();
    Ok(IngestReport {
        files: 1,
        rows: total,
        skipped_no_position: skipped.no_position,
        skipped_no_alt: skipped.no_alt,
    })
}

fn append_clinvar_rows(
    record: &vcf::Record,
    header: &vcf::Header,
    buf: &mut ClinVarBuffer,
    skipped: &mut SkipCounts,
) -> Result<()> {
    let chrom = record.reference_sequence_name().to_string();
    let Some(pos_result) = record.variant_start() else {
        skipped.no_position += 1;
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

    let mut any_alt = false;
    for alt_result in record.alternate_bases().iter() {
        let alt = alt_result?.to_string();
        if alt == "." || alt.is_empty() {
            continue;
        }
        any_alt = true;
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
    if !any_alt {
        skipped.no_alt += 1;
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
