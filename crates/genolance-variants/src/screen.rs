//! `genolance screen` - combined carrier-screen across N samples with a
//! ClinVar pathogenicity filter. Surfaces sites where every listed
//! sample carries at least one ALT allele AND the site is annotated
//! with a matching ClinVar clinical-significance string.
//!
//! This is the combined workflow a couple or trio typically cares
//! about (compound carrier screening for recessive disease); doing it
//! as a first-class command avoids the previous two-step pipeline of
//! `compare --mode carrier-screen` + `join --significance pathogenic`.

use std::collections::{HashMap, HashSet};

use anyhow::{anyhow, Result};
use arrow_array::{Array, RecordBatch};
use futures::TryStreamExt;
use lancedb::query::{ExecutableQuery, QueryBase};

use genolance_core::schema::{CLINVAR_TABLE, VARIANTS_TABLE};
use genolance_core::store::Store;

use crate::common::{
    build_clinvar_map, carrier_state, ensure_clinvar_unique, fmt_dp, fmt_qual, normalize_chrom,
    opt_string, req_f32, req_str, req_u32, req_u64, significance_predicate, sql_escape, truncate,
    validate_filters, Annotation, Carrier, Gate, Key,
};
pub use crate::filters::Filters;

/// Default substring needle when the user doesn't pass one.
const DEFAULT_SIGNIFICANCE: &str = "pathogenic";

pub async fn run(store_path: &str, samples: &[String], f: &Filters) -> Result<()> {
    if samples.len() < 2 {
        return Err(anyhow!(
            "`screen` needs at least two samples (compound carrier screen). Got {}",
            samples.len()
        ));
    }
    validate_filters(f)?;

    let store = Store::open(store_path).await?;
    let tables = store.variants.table_names().execute().await?;
    if !tables.iter().any(|n| n == VARIANTS_TABLE) {
        return Err(anyhow!(
            "store has no '{VARIANTS_TABLE}' table; ingest sample VCFs first"
        ));
    }
    if !tables.iter().any(|n| n == CLINVAR_TABLE) {
        return Err(anyhow!(
            "store has no '{CLINVAR_TABLE}' table; ingest clinvar.vcf.gz first"
        ));
    }

    // 1. Pull filtered ClinVar into an in-memory map keyed on (chrom, pos, ref, alt).
    let clinvar = store.variants.open_table(CLINVAR_TABLE).execute().await?;
    ensure_clinvar_unique(&clinvar).await?;
    let cv_pred = build_clinvar_predicate(f)?;
    let cv_batches: Vec<RecordBatch> = clinvar
        .query()
        .only_if(cv_pred)
        .execute()
        .await?
        .try_collect()
        .await?;

    let cv_map = build_clinvar_map(&cv_batches, f)?;
    eprintln!(
        "[screen] {} ClinVar records in lookup (after filters)",
        cv_map.len()
    );
    if cv_map.is_empty() {
        return Ok(());
    }

    // 2. Pull variants for the listed samples. Thresholds are pushed into SQL
    // with NULLs kept, so missing QUAL/DP can be counted below.
    let sample_pred = samples
        .iter()
        .map(|s| format!("sample_name = '{}'", sql_escape(s)))
        .collect::<Vec<_>>()
        .join(" OR ");
    let mut pred = format!("({sample_pred})");
    if let Some(q) = f.min_qual {
        pred.push_str(&format!(" AND (quality IS NULL OR quality >= {q})"));
    }
    if let Some(d) = f.min_dp {
        pred.push_str(&format!(" AND (read_depth IS NULL OR read_depth >= {d})"));
    }
    let variants = store.variants.open_table(VARIANTS_TABLE).execute().await?;
    let v_batches: Vec<RecordBatch> = variants
        .query()
        .only_if(pred)
        .execute()
        .await?
        .try_collect()
        .await?;

    // 3. Index carriers per site: Key -> { sample -> (gt, qual, dp) }.
    type PerSample = HashMap<String, (String, Option<f32>, Option<u32>)>;
    let mut carriers: HashMap<Key, PerSample> = HashMap::new();
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
            if !cv_map.contains_key(&key) {
                continue;
            }
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
            carriers.entry(key).or_default().insert(
                sample.value(i).to_string(),
                (gt_s.unwrap_or_default(), q, d),
            );
        }
    }

    // 4. Emit sites where every requested sample is a carrier.
    let want: HashSet<&str> = samples.iter().map(String::as_str).collect();
    let mut hits: Vec<(&Key, &Annotation, &PerSample)> = Vec::new();
    for (key, ann) in cv_map.iter() {
        let Some(c) = carriers.get(key) else { continue };
        if want.iter().all(|w| c.contains_key(*w)) {
            hits.push((key, ann, c));
        }
    }
    hits.sort_by(|a, b| (&a.0 .0, a.0 .1).cmp(&(&b.0 .0, b.0 .1)));

    println!(
        "{:<8} {:<12} {:<4} {:<4} {:<12} {:<30} disease",
        "chrom", "pos", "ref", "alt", "gene", "significance"
    );
    println!("{}", "-".repeat(120));
    for (key, ann, per) in &hits {
        println!(
            "{:<8} {:<12} {:<4} {:<4} {:<12} {:<30} {}",
            key.0,
            key.1,
            truncate(&key.2, 4),
            truncate(&key.3, 4),
            ann.gene_symbol.as_deref().unwrap_or(""),
            truncate(ann.clinical_significance.as_deref().unwrap_or(""), 30),
            ann.disease_name.as_deref().unwrap_or(""),
        );
        let mut rows: Vec<_> = per.iter().collect();
        rows.sort_by(|a, b| a.0.cmp(b.0));
        for (sname, (gt, q, d)) in rows {
            println!(
                "    -> {:<12} GT={:<5} QUAL={:>5} DP={:>4}",
                sname,
                gt,
                fmt_qual(*q),
                fmt_dp(*d)
            );
        }
    }
    println!("{}", "-".repeat(120));
    println!(
        "{} sites where all of [{}] carry a matching ClinVar variant",
        hits.len(),
        samples.join(", ")
    );
    gate.report("screen");
    Ok(())
}

fn build_clinvar_predicate(f: &Filters) -> Result<String> {
    significance_predicate(f, Some(DEFAULT_SIGNIFICANCE))
}

#[cfg(test)]
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
    fn build_clinvar_predicate_defaults_to_pathogenic_substring() {
        // Unlike join's predicate, screen defaults to "pathogenic" rather
        // than matching everything when no filter is given.
        assert_eq!(
            build_clinvar_predicate(&Filters::default()).unwrap(),
            "lower(clinical_significance) LIKE '%pathogenic%'"
        );
    }

    #[test]
    fn build_clinvar_predicate_prefers_exact_list_over_substring() {
        let f = Filters {
            significance_exact: Some(vec!["Pathogenic".to_string()]),
            ..Default::default()
        };
        assert_eq!(
            build_clinvar_predicate(&f).unwrap(),
            "(lower(clinical_significance) = 'pathogenic')"
        );
    }

    #[test]
    fn build_clinvar_map_normalizes_chrom_and_keys_by_variant() {
        let batch = clinvar_batch(&[("chr1", 100, "A", "T", "BRCA1", "Pathogenic", "Cancer")]);
        let map = build_clinvar_map(&[batch], &Filters::default()).unwrap();
        assert_eq!(map.len(), 1);
        assert!(map.contains_key(&("1".to_string(), 100, "A".to_string(), "T".to_string())));
    }

    #[test]
    fn build_clinvar_map_acmg_only_filters_out_non_acmg_genes() {
        let batch = clinvar_batch(&[
            ("1", 100, "A", "T", "BRCA1", "Pathogenic", "Cancer"), // in ACMG_SF_V3
            ("1", 200, "C", "G", "NOTAREALGENE", "Pathogenic", "Other"),
        ]);
        let f = Filters {
            acmg_only: true,
            ..Default::default()
        };
        let map = build_clinvar_map(&[batch], &f).unwrap();
        assert_eq!(map.len(), 1);
        assert!(map.contains_key(&("1".to_string(), 100, "A".to_string(), "T".to_string())));
    }
}
