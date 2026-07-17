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
use arrow_array::{Array, Float32Array, RecordBatch, StringArray, UInt32Array, UInt64Array};
use futures::TryStreamExt;
use lancedb::query::{ExecutableQuery, QueryBase};

use genolance_core::gene_lists::{gene_in_set, ACMG_SF_V3};
use genolance_core::schema::{CLINVAR_TABLE, VARIANTS_TABLE};
use genolance_core::store::Store;

type Key = (String, u64, String, String); // chrom-normalized, pos, ref, alt

/// Default substring needle when the user doesn't pass one.
const DEFAULT_SIGNIFICANCE: &str = "pathogenic";

#[derive(Default, Clone)]
pub struct Filters {
    pub significance_substring: Option<String>,
    pub significance_exact: Option<Vec<String>>,
    pub min_qual: Option<f32>,
    pub min_dp: Option<u32>,
    pub acmg_only: bool,
    pub gene_filter: Option<Vec<String>>,
}

pub async fn run(store_path: &str, samples: &[String], f: &Filters) -> Result<()> {
    if samples.len() < 2 {
        return Err(anyhow!(
            "`screen` needs at least two samples (compound carrier screen). Got {}",
            samples.len()
        ));
    }

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
    let cv_pred = build_clinvar_predicate(f);
    let cv_batches: Vec<RecordBatch> = clinvar
        .query()
        .only_if(cv_pred)
        .execute()
        .await?
        .try_collect()
        .await?;

    let cv_map = build_clinvar_map(&cv_batches, f);
    eprintln!(
        "[screen] {} ClinVar records in lookup (after filters)",
        cv_map.len()
    );
    if cv_map.is_empty() {
        return Ok(());
    }

    // 2. Pull variants for the listed samples with quality gates applied in SQL.
    let sample_pred = samples
        .iter()
        .map(|s| format!("sample_name = '{}'", sql_escape(s)))
        .collect::<Vec<_>>()
        .join(" OR ");
    let mut pred = format!("({sample_pred})");
    if let Some(q) = f.min_qual {
        pred.push_str(&format!(" AND quality >= {}", q));
    }
    if let Some(d) = f.min_dp {
        pred.push_str(&format!(" AND read_depth >= {}", d));
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
    type PerSample = HashMap<String, (String, f32, u32)>;
    let mut carriers: HashMap<Key, PerSample> = HashMap::new();
    for b in &v_batches {
        let sample = col_str(b, "sample_name");
        let chrom = col_str(b, "chrom");
        let pos = col_u64(b, "pos");
        let ref_a = col_str(b, "ref_allele");
        let alt_a = col_str(b, "alt_allele");
        let gt = col_str(b, "genotype");
        let qual = col_f32(b, "quality");
        let dp = col_u32(b, "read_depth");
        let (Some(sample), Some(chrom), Some(pos), Some(ref_a), Some(alt_a)) =
            (sample, chrom, pos, ref_a, alt_a)
        else {
            continue;
        };
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
            let gt_s = gt.map(|g| g.value(i)).unwrap_or("");
            if !gt_s.chars().any(|c| c == '1') {
                continue;
            }
            carriers.entry(key).or_default().insert(
                sample.value(i).to_string(),
                (
                    gt_s.to_string(),
                    qual.map(|a| a.value(i)).unwrap_or(0.0),
                    dp.map(|a| a.value(i)).unwrap_or(0),
                ),
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
        let mut rows: Vec<(&String, &(String, f32, u32))> = per.iter().collect();
        rows.sort_by(|a, b| a.0.cmp(b.0));
        for (sname, (gt, q, d)) in rows {
            println!(
                "    ↳ {:<12} GT={:<5} QUAL={:>5.1} DP={:>4}",
                sname, gt, q, d
            );
        }
    }
    println!("{}", "-".repeat(120));
    println!(
        "{} sites where all of [{}] carry a matching ClinVar variant",
        hits.len(),
        samples.join(", ")
    );
    Ok(())
}

struct Annotation {
    gene_symbol: Option<String>,
    clinical_significance: Option<String>,
    disease_name: Option<String>,
}

fn build_clinvar_predicate(f: &Filters) -> String {
    if let Some(list) = &f.significance_exact {
        let ors: Vec<String> = list
            .iter()
            .map(|s| {
                format!(
                    "lower(clinical_significance) = '{}'",
                    sql_escape(&s.to_lowercase())
                )
            })
            .collect();
        if ors.is_empty() {
            "1=1".to_string()
        } else {
            format!("({})", ors.join(" OR "))
        }
    } else {
        let sub = f
            .significance_substring
            .as_deref()
            .unwrap_or(DEFAULT_SIGNIFICANCE);
        format!(
            "lower(clinical_significance) LIKE '%{}%'",
            sql_escape(&sub.to_lowercase())
        )
    }
}

fn build_clinvar_map(batches: &[RecordBatch], f: &Filters) -> HashMap<Key, Annotation> {
    let mut map = HashMap::new();
    for b in batches {
        let chrom = col_str(b, "chrom");
        let pos = col_u64(b, "pos");
        let ref_a = col_str(b, "ref_allele");
        let alt_a = col_str(b, "alt_allele");
        let gene = col_str(b, "gene_symbol");
        let sig = col_str(b, "clinical_significance");
        let dis = col_str(b, "disease_name");
        let (Some(chrom), Some(pos), Some(ref_a), Some(alt_a)) = (chrom, pos, ref_a, alt_a) else {
            continue;
        };
        for i in 0..b.num_rows() {
            let gene_s = gene
                .and_then(|a| (!a.is_null(i)).then(|| a.value(i).to_string()))
                .unwrap_or_default();
            if f.acmg_only && !gene_in_set(&gene_s, ACMG_SF_V3) {
                continue;
            }
            if let Some(list) = &f.gene_filter {
                if !list.iter().any(|g| g.eq_ignore_ascii_case(&gene_s)) {
                    continue;
                }
            }
            let key: Key = (
                normalize_chrom(chrom.value(i)),
                pos.value(i),
                ref_a.value(i).to_string(),
                alt_a.value(i).to_string(),
            );
            map.insert(
                key,
                Annotation {
                    gene_symbol: if gene_s.is_empty() {
                        None
                    } else {
                        Some(gene_s)
                    },
                    clinical_significance: sig
                        .and_then(|a| (!a.is_null(i)).then(|| a.value(i).to_string())),
                    disease_name: dis.and_then(|a| (!a.is_null(i)).then(|| a.value(i).to_string())),
                },
            );
        }
    }
    map
}

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

fn normalize_chrom(c: &str) -> String {
    c.strip_prefix("chr").unwrap_or(c).to_string()
}
fn sql_escape(s: &str) -> String {
    s.replace('\'', "''")
}
fn truncate(s: &str, n: usize) -> String {
    if s.chars().count() > n {
        let cut: String = s.chars().take(n.saturating_sub(1)).collect();
        format!("{cut}...")
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
    fn normalize_chrom_strips_prefix() {
        assert_eq!(normalize_chrom("chr1"), "1");
        assert_eq!(normalize_chrom("1"), "1");
    }

    #[test]
    fn sql_escape_doubles_single_quotes() {
        assert_eq!(sql_escape("O'Brien"), "O''Brien");
    }

    #[test]
    fn truncate_long_string_gets_ellipsis() {
        assert_eq!(truncate("abcdefgh", 4), "abc...");
    }

    #[test]
    fn build_clinvar_predicate_defaults_to_pathogenic_substring() {
        // Unlike join's predicate, screen defaults to "pathogenic" rather
        // than matching everything when no filter is given.
        assert_eq!(
            build_clinvar_predicate(&Filters::default()),
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
            build_clinvar_predicate(&f),
            "(lower(clinical_significance) = 'pathogenic')"
        );
    }

    #[test]
    fn build_clinvar_map_normalizes_chrom_and_keys_by_variant() {
        let batch = clinvar_batch(&[("chr1", 100, "A", "T", "BRCA1", "Pathogenic", "Cancer")]);
        let map = build_clinvar_map(&[batch], &Filters::default());
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
        let map = build_clinvar_map(&[batch], &f);
        assert_eq!(map.len(), 1);
        assert!(map.contains_key(&("1".to_string(), 100, "A".to_string(), "T".to_string())));
    }
}
