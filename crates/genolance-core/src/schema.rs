use arrow_schema::{DataType, Field, Schema};
use std::sync::Arc;

/// Arrow schema for a called variant record.
///
/// One row per ALT allele per sample. Multi-allelic sites are split
/// into one row per ALT so positional joins stay simple, but enough
/// raw fields are duplicated across split rows to allow a lossless
/// regrouped VCF export.
pub fn variant_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("sample_name", DataType::Utf8, false),
        Field::new("chrom", DataType::Utf8, false),
        Field::new("pos", DataType::UInt64, false), // 1-based
        Field::new("ids", DataType::Utf8, true),    // VCF ID column (`;`-joined)
        Field::new("ref_allele", DataType::Utf8, false),
        Field::new("alt_allele", DataType::Utf8, false),
        Field::new("alt_index", DataType::UInt32, false), // 0-based position in ALT list
        Field::new("alt_count", DataType::UInt32, false), // total ALT count at this site
        Field::new("quality", DataType::Float32, true),
        Field::new("filter", DataType::Utf8, true), // `;`-joined; None = "."
        // Per-alt split GT: biallelic encoding relative to this ALT ("0/1", "1/1", "./.").
        Field::new("genotype", DataType::Utf8, true),
        // Original, un-split GT ("1/2", "0|1", ...) - preserved for lossless export.
        Field::new("gt_raw", DataType::Utf8, true),
        Field::new("read_depth", DataType::UInt32, true), // FORMAT/DP
        Field::new("format_ad", DataType::Utf8, true),    // allele depths, e.g. "10,5"
        Field::new("format_gq", DataType::UInt32, true),  // genotype quality
        Field::new("format_pl", DataType::Utf8, true),    // phred likelihoods, e.g. "0,30,255"
        Field::new("allele_freq", DataType::Float32, true), // per-alt AF
        // The record's FORMAT column exactly as written ("GT:GQ:DP:AD"), and
        // the chosen sample's column exactly as written ("0/1:35:10:5,5").
        // `genolance export` writes them back unchanged. The typed columns
        // above (gt_raw, read_depth, ...) are derived from the same text.
        Field::new("format_keys", DataType::Utf8, true),
        Field::new("format_raw", DataType::Utf8, true),
        // The INFO column exactly as written. Duplicated across all per-alt
        // split rows of a record so any one can reconstruct it.
        Field::new("info_raw", DataType::Utf8, true),
        // 0-based position of the source record among all records of the
        // file it came from. Rows with the same (sample_name, source_ordinal)
        // are the split alleles of one record.
        Field::new("source_ordinal", DataType::UInt64, false),
        // Left-trimmed, right-trimmed form of (pos, ref, alt) used to join
        // against ClinVar and other annotation sources. The raw columns above
        // stay as written for export.
        Field::new("norm_pos", DataType::UInt64, false),
        Field::new("norm_ref", DataType::Utf8, false),
        Field::new("norm_alt", DataType::Utf8, false),
    ]))
}

/// Arrow schema for ClinVar annotation records.
pub fn clinvar_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("chrom", DataType::Utf8, false),
        Field::new("pos", DataType::UInt64, false),
        Field::new("ref_allele", DataType::Utf8, false),
        Field::new("alt_allele", DataType::Utf8, false),
        Field::new("variation_id", DataType::Utf8, true),
        // First gene symbol of GENEINFO.
        Field::new("gene_symbol", DataType::Utf8, true),
        Field::new("clinical_significance", DataType::Utf8, true),
        Field::new("review_status", DataType::Utf8, true),
        Field::new("disease_name", DataType::Utf8, true),
        // Every gene symbol of GENEINFO, in file order; null when the record
        // has no GENEINFO. `gene_symbol` is the first one.
        Field::new(
            "gene_symbols",
            DataType::List(Arc::new(Field::new("item", DataType::Utf8, true))),
            true,
        ),
        Field::new("norm_pos", DataType::UInt64, false),
        Field::new("norm_ref", DataType::Utf8, false),
        Field::new("norm_alt", DataType::Utf8, false),
    ]))
}

/// Arrow schema for the `samples` registry table. One row per ingested
/// sample keyed by `sample_name`. Holds the raw VCF header text so
/// `genolance export` can reconstruct a valid header on the way out.
pub fn samples_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("sample_name", DataType::Utf8, false),
        Field::new("source_path", DataType::Utf8, true),
        Field::new("vcf_header", DataType::Utf8, false),
        Field::new("ingested_at", DataType::Utf8, true), // ISO8601 UTC
        Field::new("reference", DataType::Utf8, true),   // from ##reference line
        // false while an ingest is running or after one was interrupted.
        Field::new("ingest_complete", DataType::Boolean, false),
        // Contig lengths from the header as "name:length,name:length" (only
        // contigs that declare a length). Used to check the reference against
        // the ClinVar release.
        Field::new("contigs", DataType::Utf8, true),
        // Reference assembly label ("GRCh37", "GRCh38", ...) inferred from the
        // header or declared with `--assembly`; null when unknown.
        Field::new("assembly", DataType::Utf8, true),
    ]))
}

/// Name of the variant calls table (in the variants-layer connection).
/// On disk: `<store>/variants/calls.lance/`
pub const VARIANTS_TABLE: &str = "calls";

/// Name of the ClinVar annotations table (in the variants-layer connection).
/// On disk: `<store>/variants/clinvar.lance/`
pub const CLINVAR_TABLE: &str = "clinvar";

/// Name of the sample registry table (in the root connection).
/// On disk: `<store>/samples.lance/`
pub const SAMPLES_TABLE: &str = "samples";

/// Name of the store metadata table (key/value, in the root connection).
/// On disk: `<store>/store_meta.lance/`
pub const STORE_META_TABLE: &str = "store_meta";

/// Arrow schema for the `store_meta` table.
pub fn store_meta_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        Field::new("key", DataType::Utf8, false),
        Field::new("value", DataType::Utf8, false),
    ]))
}

/// Metadata key: store format version.
pub const META_FORMAT_VERSION: &str = "format_version";
/// Metadata key: `complete` or `incomplete` for the stored ClinVar release.
pub const META_CLINVAR_STATE: &str = "clinvar_ingest";
/// Metadata key: ClinVar `##fileDate`.
pub const META_CLINVAR_RELEASE: &str = "clinvar_release";
/// Metadata key: ClinVar source file name.
pub const META_CLINVAR_SOURCE: &str = "clinvar_source";
/// Metadata key: ClinVar `##reference`.
pub const META_CLINVAR_REFERENCE: &str = "clinvar_reference";
/// Metadata key: ClinVar assembly label (absent when unknown).
pub const META_CLINVAR_ASSEMBLY: &str = "clinvar_assembly";
/// Metadata key: ClinVar contig lengths as "name:length,name:length".
pub const META_CLINVAR_CONTIGS: &str = "clinvar_contigs";

/// Store format version, `MAJOR.MINOR`, recorded in `store_meta`.
///
/// Compatibility rule:
/// - A different MAJOR is refused (newer: upgrade genolance; older: the store
///   cannot be upgraded in place and must be rebuilt).
/// - The same MAJOR with any MINOR can be read. A new MINOR only adds nullable
///   columns and metadata keys; a command that needs a column an older store
///   lacks stops with `Error::MissingColumn`, naming the table and column.
/// - Writing (ingest) needs a MINOR equal to the build's: a store with a newer
///   MINOR is read-only for this build.
/// - Changing or removing a column, or changing what a value means, bumps MAJOR.
pub const STORE_FORMAT_MAJOR: u32 = 1;
/// See [`STORE_FORMAT_MAJOR`].
pub const STORE_FORMAT_MINOR: u32 = 0;

/// The format version this build writes, as stored in `store_meta`.
pub fn store_format_version() -> String {
    format!("{STORE_FORMAT_MAJOR}.{STORE_FORMAT_MINOR}")
}
