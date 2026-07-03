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
        // Non-standard FORMAT fields we don't have dedicated columns for,
        // serialized as "KEY=VAL;KEY=VAL". Preserved so `genolance export`
        // can re-emit them after the standard GT:AD:DP:GQ:PL block.
        Field::new("format_extra", DataType::Utf8, true),
        // Original FORMAT key order as ":" joined string (e.g. "GT:GQ:DP:AD:VAF:PL").
        // Used by export to reproduce the exact FORMAT column order from the source VCF.
        Field::new("format_key_order", DataType::Utf8, true),
        // Full INFO column re-serialized as "K=V;K=V" text. Duplicated across
        // all per-alt split rows of a site so any one can reconstruct it.
        Field::new("info_raw", DataType::Utf8, true),
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
        Field::new("gene_symbol", DataType::Utf8, true),
        Field::new("clinical_significance", DataType::Utf8, true),
        Field::new("review_status", DataType::Utf8, true),
        Field::new("disease_name", DataType::Utf8, true),
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
