//! `biolance compound-het` — flag genes where a single sample carries
//! two or more pathogenic/likely-pathogenic heterozygous variants, or
//! at least one homozygous P/LP variant.
//!
//! Short-read WGS can't phase variants that are more than a read-pair
//! apart, so we can't distinguish *in cis* from *in trans*. This
//! command surfaces genes that are worth resolving via family phasing
//! or Sanger follow-up; it does NOT diagnose compound heterozygosity.

use std::collections::HashMap;

use anyhow::{anyhow, Result};
use arrow_array::{Array, Float32Array, RecordBatch, StringArray, UInt32Array, UInt64Array};
use futures::TryStreamExt;
use lancedb::query::{ExecutableQuery, QueryBase};

use biolance_core::gene_lists::{gene_in_set, ACMG_SF_V3};
use biolance_core::schema::{CLINVAR_TABLE, VARIANTS_TABLE};
use biolance_core::store::Store;

type Key = (String, u64, String, String);

#[derive(Clone)]
struct Hit {
    chrom: String,
    pos: u64,
    ref_allele: String,
    alt_allele: String,
    significance: String,
    disease: String,
    gt: String,
    qual: f32,
    dp: u32,
}

#[derive(Default, Clone)]
pub struct Filters {
    pub significance_exact: Option<Vec<String>>,
    pub significance_substring: Option<String>,
    pub min_qual: Option<f32>,
    pub min_dp: Option<u32>,
    pub acmg_only: bool,
    pub gene_filter: Option<Vec<String>>,
}

/// Default substring needle when no filter is supplied.
const DEFAULT_SIGNIFICANCE: &str = "pathogenic";

pub async fn run(store_path: &str, sample: &str, f: &Filters) -> Result<()> {
    let store = Store::open(store_path).await?;
    let tables = store.variants.table_names().execute().await?;
    if !tables.iter().any(|n| n == VARIANTS_TABLE) {
        return Err(anyhow!("store has no '{VARIANTS_TABLE}' table"));
    }
    if !tables.iter().any(|n| n == CLINVAR_TABLE) {
        return Err(anyhow!("store has no '{CLINVAR_TABLE}' table"));
    }

    // 1. Build ClinVar lookup with the requested filter.
    let clinvar = store.variants.open_table(CLINVAR_TABLE).execute().await?;
    let cv_pred = build_clinvar_predicate(f);
    let mut cv_q = clinvar.query();
    if let Some(p) = cv_pred.as_deref() {
        cv_q = cv_q.only_if(p.to_string());
    }
    let cv_batches: Vec<RecordBatch> = cv_q.execute().await?.try_collect().await?;
    let cv_map = build_clinvar_map(&cv_batches, f);
    eprintln!(
        "[compound-het] {} ClinVar records in lookup (after filters)",
        cv_map.len()
    );
    if cv_map.is_empty() {
        return Ok(());
    }

    // 2. Pull the sample's variants.
    let variants = store.variants.open_table(VARIANTS_TABLE).execute().await?;
    let mut pred = format!("sample_name = '{}'", sql_escape(sample));
    if let Some(q) = f.min_qual {
        pred.push_str(&format!(" AND quality >= {}", q));
    }
    if let Some(d) = f.min_dp {
        pred.push_str(&format!(" AND read_depth >= {}", d));
    }
    let v_batches: Vec<RecordBatch> = variants
        .query()
        .only_if(pred)
        .execute()
        .await?
        .try_collect()
        .await?;

    // 3. Per-gene collection of carrier hits.
    let mut by_gene: HashMap<String, Vec<Hit>> = HashMap::new();
    for b in &v_batches {
        let chrom = col_str(b, "chrom");
        let pos = col_u64(b, "pos");
        let ref_a = col_str(b, "ref_allele");
        let alt_a = col_str(b, "alt_allele");
        let gt = col_str(b, "genotype");
        let qual = col_f32(b, "quality");
        let dp = col_u32(b, "read_depth");
        let (Some(chrom), Some(pos), Some(ref_a), Some(alt_a)) = (chrom, pos, ref_a, alt_a) else {
            continue;
        };
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
            let gt_s = gt.map(|g| g.value(i)).unwrap_or("").to_string();
            if !gt_s.chars().any(|c| c == '1') {
                continue;
            }
            let gene = ann.gene.clone();
            if gene.is_empty() {
                continue;
            }
            by_gene.entry(gene.clone()).or_default().push(Hit {
                chrom: chrom.value(i).to_string(),
                pos: pos.value(i),
                ref_allele: ref_a.value(i).to_string(),
                alt_allele: alt_a.value(i).to_string(),
                significance: ann.significance.clone(),
                disease: ann.disease.clone(),
                gt: gt_s,
                qual: qual.map(|a| a.value(i)).unwrap_or(0.0),
                dp: dp.map(|a| a.value(i)).unwrap_or(0),
            });
        }
    }

    // 4. Report: genes with >=2 hets, or >=1 hom-alt.
    let mut genes: Vec<(String, Vec<Hit>)> = by_gene
        .into_iter()
        .filter(|(_, hits)| {
            let hets = hits.iter().filter(|h| is_het(&h.gt)).count();
            let homs = hits.iter().filter(|h| is_hom_alt(&h.gt)).count();
            hets >= 2 || homs >= 1
        })
        .collect();
    genes.sort_by(|a, b| a.0.cmp(&b.0));

    if genes.is_empty() {
        println!("[compound-het] no genes flagged for {sample}");
        return Ok(());
    }

    println!("Compound-het / homozygous P/LP screen — sample {sample}");
    println!("{}", "-".repeat(100));
    for (gene, hits) in &genes {
        let hets = hits.iter().filter(|h| is_het(&h.gt)).count();
        let homs = hits.iter().filter(|h| is_hom_alt(&h.gt)).count();
        let disposition = if homs > 0 {
            format!("HOMOZYGOUS ({} hom{})", homs, plural(homs))
        } else {
            format!(
                "POSSIBLE COMPOUND HET ({} het P/LP — phase unknown)",
                hets
            )
        };
        println!("\n{} — {}", gene, disposition);
        for h in hits {
            println!(
                "  {}:{} {}>{}  GT={:<5} QUAL={:>5.1} DP={:>3}  [{}]  {}",
                h.chrom,
                h.pos,
                truncate(&h.ref_allele, 8),
                truncate(&h.alt_allele, 8),
                h.gt,
                h.qual,
                h.dp,
                truncate(&h.significance, 30),
                truncate(&h.disease, 60),
            );
        }
    }
    println!("\n{}", "-".repeat(100));
    println!("{} gene(s) flagged", genes.len());
    println!(
        "note: heterozygous pairs could be *in cis* (carrier) or *in trans* (affected). \
         Resolve with family phasing or Sanger."
    );
    Ok(())
}

struct Ann {
    gene: String,
    significance: String,
    disease: String,
}

fn build_clinvar_predicate(f: &Filters) -> Option<String> {
    let mut parts: Vec<String> = Vec::new();
    if let Some(list) = &f.significance_exact {
        // Exact case-insensitive match on any of the provided strings.
        let ors: Vec<String> = list
            .iter()
            .map(|s| format!("lower(clinical_significance) = '{}'", sql_escape(&s.to_lowercase())))
            .collect();
        parts.push(format!("({})", ors.join(" OR ")));
    } else {
        let sub = f.significance_substring.as_deref().unwrap_or(DEFAULT_SIGNIFICANCE);
        parts.push(format!(
            "lower(clinical_significance) LIKE '%{}%'",
            sql_escape(&sub.to_lowercase())
        ));
    }
    Some(parts.join(" AND "))
}

fn build_clinvar_map(batches: &[RecordBatch], f: &Filters) -> HashMap<Key, Ann> {
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
            // Gene filters are applied client-side (gene_symbol is a free-form
            // GENEINFO string in ClinVar, not easily SQL-filterable).
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
                Ann {
                    gene: gene_s,
                    significance: sig
                        .and_then(|a| (!a.is_null(i)).then(|| a.value(i).to_string()))
                        .unwrap_or_default(),
                    disease: dis
                        .and_then(|a| (!a.is_null(i)).then(|| a.value(i).to_string()))
                        .unwrap_or_default(),
                },
            );
        }
    }
    map
}

// ---- shared helpers (kept local to avoid crate-wide API churn) --------

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
        format!("{cut}…")
    } else {
        s.to_string()
    }
}
fn is_het(gt: &str) -> bool {
    // "0/1", "1/0", "0|1", "1|0", "1/2", etc. — exactly one copy of '1' *or*
    // any heterozygous call (counts mismatch between the two alleles).
    let alleles: Vec<&str> = gt.split(|c| c == '/' || c == '|').collect();
    if alleles.len() != 2 {
        return false;
    }
    let a = alleles[0];
    let b = alleles[1];
    a != b && (a == "1" || b == "1" || a == "0" || b == "0")
}
fn is_hom_alt(gt: &str) -> bool {
    matches!(gt, "1/1" | "1|1")
}
fn plural(n: usize) -> &'static str {
    if n == 1 {
        ""
    } else {
        "s"
    }
}
