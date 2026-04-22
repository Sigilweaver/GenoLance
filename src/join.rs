use std::collections::HashMap;

use anyhow::{anyhow, Result};
use arrow_array::{Array, Float32Array, RecordBatch, StringArray, UInt32Array, UInt64Array};
use futures::TryStreamExt;
use lancedb::query::{ExecutableQuery, QueryBase};

use crate::gene_lists::{gene_in_set, ACMG_SF_V3};
use crate::schema::{CLINVAR_TABLE, VARIANTS_TABLE};
use crate::store::Store;

type Key = (String, u64, String, String);

#[derive(Default, Clone)]
pub struct Filters {
    pub significance_substring: Option<String>,
    pub significance_exact: Option<Vec<String>>,
    pub min_qual: Option<f32>,
    pub min_dp: Option<u32>,
    pub acmg_only: bool,
    pub gene_filter: Option<Vec<String>>,
}

/// Annotate variants in the store by joining against the `clinvar` table.
///
/// The join is performed in-memory on (chrom, pos, ref, alt). The
/// ClinVar side is filtered by significance (substring or exact-list)
/// and optionally by gene set (ACMG SF v3 or user-supplied). The
/// variants side is filtered by `min_qual` / `min_dp` in SQL so we
/// stream fewer rows out of storage.
pub async fn run(
    store_path: &str,
    annotation_vcf: &str,
    f: &Filters,
) -> Result<()> {
    let store = Store::open(store_path).await?;

    // Make sure the ClinVar table exists — ingest it on demand if missing.
    let tables = store.conn.table_names().execute().await?;
    if !tables.iter().any(|n| n == CLINVAR_TABLE) {
        println!("[join] '{CLINVAR_TABLE}' table not found; ingesting {annotation_vcf} …");
        crate::ingest::run(
            store_path,
            std::slice::from_ref(&annotation_vcf.to_string()),
            None,
        )
        .await?;
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

    let variants = store.conn.open_table(VARIANTS_TABLE).execute().await?;
    let clinvar = store.conn.open_table(CLINVAR_TABLE).execute().await?;

    // ClinVar predicate: significance-exact OR substring, default substring "pathogenic"
    // only applies when no filter is given at all.
    let cv_pred = build_clinvar_predicate(f);
    let cv_batches: Vec<RecordBatch> = clinvar
        .query()
        .only_if(cv_pred)
        .execute()
        .await?
        .try_collect()
        .await?;
    let cv_map = build_clinvar_map(&cv_batches, f);
    println!("[join] {} ClinVar records in lookup (after filters)", cv_map.len());
    if cv_map.is_empty() {
        return Ok(());
    }

    // Variants predicate: push min_qual / min_dp into SQL.
    let mut pred_parts: Vec<String> = Vec::new();
    if let Some(q) = f.min_qual {
        pred_parts.push(format!("quality >= {}", q));
    }
    if let Some(d) = f.min_dp {
        pred_parts.push(format!("read_depth >= {}", d));
    }
    let mut vq = variants.query();
    if !pred_parts.is_empty() {
        vq = vq.only_if(pred_parts.join(" AND "));
    }
    let v_batches: Vec<RecordBatch> = vq.execute().await?.try_collect().await?;

    println!(
        "{:<16} {:<8} {:<12} {:<4} {:<4} {:<10} {:>5} {:>4} {:<12} {:<30} {}",
        "sample", "chrom", "pos", "ref", "alt", "gt", "qual", "dp", "gene", "significance", "disease"
    );
    println!("{}", "-".repeat(140));
    let mut matches = 0usize;
    for b in &v_batches {
        let sample = col_str(b, "sample_name");
        let chrom = col_str(b, "chrom");
        let pos = col_u64(b, "pos");
        let ref_a = col_str(b, "ref_allele");
        let alt_a = col_str(b, "alt_allele");
        let gt = col_str(b, "genotype");
        let qual = col_f32(b, "quality");
        let dp = col_u32(b, "read_depth");
        for i in 0..b.num_rows() {
            let (Some(s), Some(c), Some(p), Some(r), Some(a)) = (sample, chrom, pos, ref_a, alt_a)
            else {
                continue;
            };
            let key: Key = (
                normalize_chrom(c.value(i)),
                p.value(i),
                r.value(i).to_string(),
                a.value(i).to_string(),
            );
            if let Some(ann) = cv_map.get(&key) {
                let gt_s = gt.map(|g| g.value(i)).unwrap_or("");
                // Skip purely-reference genotypes so we only surface carriers.
                if gt_s.chars().any(|c| c == '1') || gt_s.is_empty() {
                    matches += 1;
                    println!(
                        "{:<16} {:<8} {:<12} {:<4} {:<4} {:<10} {:>5.1} {:>4} {:<12} {:<30} {}",
                        s.value(i),
                        c.value(i),
                        p.value(i),
                        truncate(r.value(i), 4),
                        truncate(a.value(i), 4),
                        gt_s,
                        qual.map(|a| a.value(i)).unwrap_or(0.0),
                        dp.map(|a| a.value(i)).unwrap_or(0),
                        ann.gene_symbol.as_deref().unwrap_or(""),
                        truncate(ann.clinical_significance.as_deref().unwrap_or(""), 30),
                        ann.disease_name.as_deref().unwrap_or(""),
                    );
                }
            }
        }
    }
    println!("{}", "-".repeat(140));
    println!("{matches} annotated variants");
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
    } else if let Some(sub) = &f.significance_substring {
        format!(
            "lower(clinical_significance) LIKE '%{}%'",
            sql_escape(&sub.to_lowercase())
        )
    } else {
        // No filter: match everything.
        "1=1".to_string()
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
                    gene_symbol: if gene_s.is_empty() { None } else { Some(gene_s) },
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

/// Normalize chromosome to "chr"-stripped form (so "chr1" and "1" compare equal).
fn normalize_chrom(c: &str) -> String {
    c.strip_prefix("chr").unwrap_or(c).to_string()
}

fn sql_escape(s: &str) -> String {
    s.replace('\'', "''")
}

fn truncate(s: &str, n: usize) -> String {
    if s.chars().count() > n {
        let cut: String = s.chars().take(n.saturating_sub(1)).collect();
        format!("{cut}…")
    } else {
        s.to_string()
    }
}
