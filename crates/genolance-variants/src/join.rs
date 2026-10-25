use anyhow::{anyhow, Result};
use arrow_array::{Array, RecordBatch};
use futures::TryStreamExt;
use lancedb::query::{ExecutableQuery, QueryBase};

use genolance_core::schema::{CLINVAR_TABLE, VARIANTS_TABLE};
use genolance_core::store::Store;

use crate::common::{
    build_clinvar_map, carrier_state, ensure_clinvar_unique, fmt_dp, fmt_qual, normalize_chrom,
    opt_string, req_f32, req_str, req_u32, req_u64, significance_predicate, truncate,
    validate_filters, Carrier, Gate,
};
pub use crate::filters::Filters;

type Key = (String, u64, String, String);

/// Annotate variants in the store by joining against the `clinvar` table.
///
/// The join is performed in-memory on (chrom, pos, ref, alt). The
/// ClinVar side is filtered by significance (substring or exact-list)
/// and optionally by gene set (ACMG SF v3 or user-supplied). Calls whose
/// FILTER is not PASS, and calls with a missing QUAL/DP when a minimum is
/// set, are excluded and counted unless `Filters` says otherwise.
///
/// `annotation_vcf` is only used to load ClinVar when the store has none. If
/// the store already holds a release, passing a file is an error: the stored
/// release is what the join uses, and a different file would be silently
/// ignored. Load a new release with `genolance ingest --replace`.
pub async fn run(store_path: &str, annotation_vcf: Option<&str>, f: &Filters) -> Result<()> {
    validate_filters(f)?;
    let store = Store::open(store_path).await?;

    // Make sure the ClinVar table exists - ingest it on demand if missing.
    let tables = store.variants.table_names().execute().await?;
    let have_clinvar = tables.iter().any(|n| n == CLINVAR_TABLE);
    match (have_clinvar, annotation_vcf) {
        (true, Some(file)) => {
            return Err(anyhow!(
                "the store already holds a ClinVar release, so `join` would ignore {file}; \
                 drop --against to use the stored release, or load the file with \
                 `genolance ingest --replace {file}`"
            ));
        }
        (false, None) => {
            return Err(anyhow!(
                "the store has no '{CLINVAR_TABLE}' table; pass --against <clinvar.vcf.gz> or run `genolance ingest` on it first"
            ));
        }
        (false, Some(file)) => {
            println!("[join] '{CLINVAR_TABLE}' table not found; ingesting {file} ...");
            crate::ingest::run(
                store_path,
                std::slice::from_ref(&file.to_string()),
                &crate::ingest::IngestOptions::new().with_clinvar(true),
            )
            .await?;
        }
        (true, None) => {}
    }
    if !tables.iter().any(|n| n == VARIANTS_TABLE)
        && !store
            .conn
            .table_names()
            .execute()
            .await?
            .iter()
            .any(|n| n == VARIANTS_TABLE)
    {
        return Err(anyhow!(
            "store has no '{VARIANTS_TABLE}' table; ingest a sample VCF first"
        ));
    }

    let variants = store.variants.open_table(VARIANTS_TABLE).execute().await?;
    let clinvar = store.variants.open_table(CLINVAR_TABLE).execute().await?;
    ensure_clinvar_unique(&clinvar).await?;

    // ClinVar predicate: significance-exact OR substring; no filter matches everything.
    let cv_pred = build_clinvar_predicate(f)?;
    let cv_batches: Vec<RecordBatch> = clinvar
        .query()
        .only_if(cv_pred)
        .execute()
        .await?
        .try_collect()
        .await?;
    let cv_map = build_clinvar_map(&cv_batches, f)?;
    println!(
        "[join] {} ClinVar records in lookup (after filters)",
        cv_map.len()
    );
    if cv_map.is_empty() {
        return Ok(());
    }

    // Variants predicate: push min_qual / min_dp into SQL, keeping NULLs so
    // they can be counted (or kept) below.
    let mut pred_parts: Vec<String> = Vec::new();
    if let Some(q) = f.min_qual {
        pred_parts.push(format!("(quality IS NULL OR quality >= {q})"));
    }
    if let Some(d) = f.min_dp {
        pred_parts.push(format!("(read_depth IS NULL OR read_depth >= {d})"));
    }
    let mut vq = variants.query();
    if !pred_parts.is_empty() {
        vq = vq.only_if(pred_parts.join(" AND "));
    }
    let v_batches: Vec<RecordBatch> = vq.execute().await?.try_collect().await?;

    println!(
        "{:<16} {:<8} {:<12} {:<4} {:<4} {:<10} {:>5} {:>4} {:<8} {:<12} {:<30} disease",
        "sample",
        "chrom",
        "pos",
        "ref",
        "alt",
        "gt",
        "qual",
        "dp",
        "filter",
        "gene",
        "significance"
    );
    println!("{}", "-".repeat(150));
    let mut matches = 0usize;
    let mut gate = Gate::new(f);
    for b in &v_batches {
        let t = VARIANTS_TABLE;
        let sample = req_str(b, t, "sample_name")?;
        let chrom = req_str(b, t, "chrom")?;
        let pos = req_u64(b, t, "pos")?;
        let ref_a = req_str(b, t, "ref_allele")?;
        let alt_a = req_str(b, t, "alt_allele")?;
        let gt = req_str(b, t, "genotype")?;
        let qual = req_f32(b, t, "quality")?;
        let dp = req_u32(b, t, "read_depth")?;
        let filter = req_str(b, t, "filter")?;
        for i in 0..b.num_rows() {
            let key: Key = (
                normalize_chrom(chrom.value(i)),
                pos.value(i),
                ref_a.value(i).to_string(),
                alt_a.value(i).to_string(),
            );
            let Some(ann) = cv_map.get(&key) else {
                continue;
            };
            let gt_s = opt_string(gt, i);
            match carrier_state(gt_s.as_deref()) {
                Carrier::Yes => {}
                Carrier::No => continue,
                Carrier::NoCall => {
                    gate.no_call += 1;
                    continue;
                }
            }
            let q = (!qual.is_null(i)).then(|| qual.value(i));
            let d = (!dp.is_null(i)).then(|| dp.value(i));
            let fl = opt_string(filter, i);
            if !gate.admit(fl.as_deref(), q, d) {
                continue;
            }
            matches += 1;
            println!(
                "{:<16} {:<8} {:<12} {:<4} {:<4} {:<10} {:>5} {:>4} {:<8} {:<12} {:<30} {}",
                sample.value(i),
                chrom.value(i),
                pos.value(i),
                truncate(ref_a.value(i), 4),
                truncate(alt_a.value(i), 4),
                gt_s.as_deref().unwrap_or("."),
                fmt_qual(q),
                fmt_dp(d),
                truncate(fl.as_deref().unwrap_or("."), 8),
                ann.gene_symbol.as_deref().unwrap_or(""),
                truncate(ann.clinical_significance.as_deref().unwrap_or(""), 30),
                ann.disease_name.as_deref().unwrap_or(""),
            );
        }
    }
    println!("{}", "-".repeat(150));
    println!("{matches} annotated variants");
    gate.report("join");
    Ok(())
}

fn build_clinvar_predicate(f: &Filters) -> Result<String> {
    significance_predicate(f, None)
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod tests {
    use super::*;
    use arrow_array::{StringArray, UInt64Array};
    use arrow_schema::{DataType, Field, Schema};
    use std::sync::Arc;

    /// Minimal ClinVar-shaped batch covering only the columns `build_clinvar_map` reads.
    fn clinvar_batch(rows: &[(&str, u64, &str, &str, &str, &str, &str)]) -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("pos", DataType::UInt64, false),
            Field::new("ref_allele", DataType::Utf8, false),
            Field::new("alt_allele", DataType::Utf8, false),
            Field::new("gene_symbol", DataType::Utf8, true),
            Field::new("clinical_significance", DataType::Utf8, true),
            Field::new("disease_name", DataType::Utf8, true),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(
                    rows.iter().map(|r| r.0).collect::<Vec<_>>(),
                )),
                Arc::new(UInt64Array::from(
                    rows.iter().map(|r| r.1).collect::<Vec<_>>(),
                )),
                Arc::new(StringArray::from(
                    rows.iter().map(|r| r.2).collect::<Vec<_>>(),
                )),
                Arc::new(StringArray::from(
                    rows.iter().map(|r| r.3).collect::<Vec<_>>(),
                )),
                Arc::new(StringArray::from(
                    rows.iter().map(|r| r.4).collect::<Vec<_>>(),
                )),
                Arc::new(StringArray::from(
                    rows.iter().map(|r| r.5).collect::<Vec<_>>(),
                )),
                Arc::new(StringArray::from(
                    rows.iter().map(|r| r.6).collect::<Vec<_>>(),
                )),
            ],
        )
        .unwrap()
    }

    #[test]
    fn build_clinvar_predicate_defaults_to_match_all_without_filters() {
        assert_eq!(build_clinvar_predicate(&Filters::default()).unwrap(), "1=1");
    }

    #[test]
    fn build_clinvar_predicate_uses_substring_filter() {
        let mut f = Filters::default();
        f.significance_substring = Some("Pathogenic".to_string());
        assert_eq!(
            build_clinvar_predicate(&f).unwrap(),
            "lower(clinical_significance) LIKE '%pathogenic%'"
        );
    }

    #[test]
    fn build_clinvar_predicate_prefers_exact_list_over_substring() {
        let mut f = Filters::default();
        f.significance_substring = Some("ignored".to_string());
        f.significance_exact = Some(vec![
            "Pathogenic".to_string(),
            "Likely_pathogenic".to_string(),
        ]);
        assert_eq!(
            build_clinvar_predicate(&f).unwrap(),
            "(lower(clinical_significance) = 'pathogenic' OR lower(clinical_significance) = 'likely_pathogenic')"
        );
    }

    #[test]
    fn build_clinvar_predicate_escapes_quotes() {
        let mut f = Filters::default();
        f.significance_substring = Some("O'Brien".to_string());
        assert_eq!(
            build_clinvar_predicate(&f).unwrap(),
            "lower(clinical_significance) LIKE '%o''brien%'"
        );
    }

    #[test]
    fn build_clinvar_map_normalizes_chrom_and_keys_by_variant() {
        let batch = clinvar_batch(&[
            ("chr1", 100, "A", "T", "BRCA1", "Pathogenic", "Cancer"),
            ("2", 200, "C", "G", "", "Benign", ""),
        ]);
        let map = build_clinvar_map(&[batch], &Filters::default()).unwrap();
        assert_eq!(map.len(), 2);

        let ann = map
            .get(&("1".to_string(), 100, "A".to_string(), "T".to_string()))
            .unwrap();
        assert_eq!(ann.gene_symbol.as_deref(), Some("BRCA1"));
        assert_eq!(ann.clinical_significance.as_deref(), Some("Pathogenic"));

        // Empty gene_symbol is normalized to None, not Some("").
        let ann2 = map
            .get(&("2".to_string(), 200, "C".to_string(), "G".to_string()))
            .unwrap();
        assert_eq!(ann2.gene_symbol, None);
    }

    #[test]
    fn build_clinvar_map_acmg_only_filters_out_non_acmg_genes() {
        let batch = clinvar_batch(&[
            ("1", 100, "A", "T", "BRCA1", "Pathogenic", "Cancer"), // in ACMG_SF_V3
            ("1", 200, "C", "G", "NOTAREALGENE", "Pathogenic", "Other"),
        ]);
        let mut f = Filters::default();
        f.acmg_only = true;
        let map = build_clinvar_map(&[batch], &f).unwrap();
        assert_eq!(map.len(), 1);
        assert!(map.contains_key(&("1".to_string(), 100, "A".to_string(), "T".to_string())));
    }

    #[test]
    fn build_clinvar_map_gene_filter_is_case_insensitive() {
        let batch = clinvar_batch(&[
            ("1", 100, "A", "T", "brca1", "Pathogenic", "Cancer"),
            ("1", 200, "C", "G", "TP53", "Pathogenic", "Other"),
        ]);
        let mut f = Filters::default();
        f.gene_filter = Some(vec!["BRCA1".to_string()]);
        let map = build_clinvar_map(&[batch], &f).unwrap();
        assert_eq!(map.len(), 1);
        assert!(map.contains_key(&("1".to_string(), 100, "A".to_string(), "T".to_string())));
    }

    #[test]
    fn empty_exact_list_is_an_error_not_match_all() {
        let mut f = Filters::default();
        f.significance_exact = Some(vec![]);
        assert!(build_clinvar_predicate(&f).is_err());
    }

    #[test]
    fn conflicting_duplicate_clinvar_rows_are_an_error() {
        let batch = clinvar_batch(&[
            ("1", 100, "A", "T", "BRCA1", "Pathogenic", "Cancer"),
            ("chr1", 100, "A", "T", "BRCA1", "Benign", "Cancer"),
        ]);
        let err = build_clinvar_map(&[batch], &Filters::default()).unwrap_err();
        assert!(format!("{err}").contains("more than one release"), "{err}");
    }
}
