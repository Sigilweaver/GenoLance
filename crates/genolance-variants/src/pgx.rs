//! `genolance pgx` - pharmacogenomic screening.
//!
//! This is intentionally narrow. Real star-allele diplotype calling
//! (CYP2D6 especially) needs phased haplotypes plus structural variant
//! calls; we don't attempt that. Instead, for a curated list of
//! pharmacogenomically-relevant genes we intersect the sample's
//! variants with ClinVar entries whose `clinical_significance` or
//! `disease_name` indicates a drug-response annotation, then report
//! each carried variant with the drug/effect text.
//!
//! Think of it as `genolance join --significance "drug response"`
//! scoped to PGx-relevant genes. Useful for surfacing review-worthy
//! variants for a clinician; *not* a diplotype call and *not*
//! a prescribing recommendation.

use std::collections::HashMap;

use anyhow::{anyhow, Result};
use arrow_array::{Array, RecordBatch, StringArray, UInt64Array};
use futures::TryStreamExt;
use lancedb::query::{ExecutableQuery, QueryBase};

use genolance_core::schema::{CLINVAR_TABLE, VARIANTS_TABLE};
use genolance_core::store::Store;

/// Genes with the strongest CPIC Tier 1 / Tier 2 guidance. Extend via
/// `--genes` on the CLI if you want to include others.
const DEFAULT_PGX_GENES: &[&str] = &[
    "CYP2C19", "CYP2C9", "CYP2D6", "CYP3A5", "CYP4F2", "DPYD", "G6PD", "IFNL3", "NUDT15",
    "SLCO1B1", "TPMT", "UGT1A1", "VKORC1",
];

/// Substrings (case-insensitive) we match against ClinVar
/// `clinical_significance` and `disease_name` to identify
/// pharmacogenomic annotations. ClinVar's CLNSIG uses `drug_response`
/// with an underscore; the disease-name column uses human prose
/// ("response to X", "X response -").
const PGX_NEEDLES: &[&str] = &["drug_response", "drug response", "response to "];

type Key = (String, u64, String, String);

pub async fn run(store_path: &str, sample: &str, extra_genes: &[String]) -> Result<()> {
    let store = Store::open(store_path).await?;
    let tables = store.variants.table_names().execute().await?;
    if !tables.iter().any(|n| n == VARIANTS_TABLE) {
        return Err(anyhow!("store has no '{VARIANTS_TABLE}' table"));
    }
    if !tables.iter().any(|n| n == CLINVAR_TABLE) {
        return Err(anyhow!(
            "store has no '{CLINVAR_TABLE}' table; ingest clinvar.vcf.gz first"
        ));
    }

    let mut genes: Vec<String> = DEFAULT_PGX_GENES.iter().map(|s| s.to_string()).collect();
    for g in extra_genes {
        if !genes.iter().any(|x| x.eq_ignore_ascii_case(g)) {
            genes.push(g.clone());
        }
    }

    // 1. Pull ClinVar entries for PGx genes, filtered to drug-response rows.
    let clinvar = store.variants.open_table(CLINVAR_TABLE).execute().await?;
    let gene_in = genes
        .iter()
        .map(|g| format!("'{}'", sql_escape(g)))
        .collect::<Vec<_>>()
        .join(", ");
    let needle_pred = PGX_NEEDLES
        .iter()
        .flat_map(|n| {
            let n = sql_escape(&n.to_lowercase());
            [
                format!("lower(clinical_significance) LIKE '%{n}%'"),
                format!("lower(disease_name) LIKE '%{n}%'"),
            ]
        })
        .collect::<Vec<_>>()
        .join(" OR ");
    let pred = format!("gene_symbol IN ({gene_in}) AND ({needle_pred})");
    let cv_batches: Vec<RecordBatch> = clinvar
        .query()
        .only_if(pred)
        .execute()
        .await?
        .try_collect()
        .await?;
    let cv_map = build_clinvar_map(&cv_batches);
    eprintln!(
        "[pgx] {} PGx ClinVar variants across {} genes",
        cv_map.len(),
        genes.len()
    );
    if cv_map.is_empty() {
        println!("(no PGx ClinVar variants in store - is clinvar.vcf.gz current?)");
        return Ok(());
    }

    // 2. Pull the sample's variants and intersect by (chrom, pos, ref, alt).
    let variants = store.variants.open_table(VARIANTS_TABLE).execute().await?;
    let v_batches: Vec<RecordBatch> = variants
        .query()
        .only_if(format!("sample_name = '{}'", sql_escape(sample)))
        .execute()
        .await?
        .try_collect()
        .await?;

    let mut hits: Vec<Hit> = Vec::new();
    for b in &v_batches {
        let chrom = col_str(b, "chrom");
        let pos = col_u64(b, "pos");
        let ref_a = col_str(b, "ref_allele");
        let alt_a = col_str(b, "alt_allele");
        let gt = col_str(b, "genotype");
        let (Some(chrom), Some(pos), Some(ref_a), Some(alt_a)) = (chrom, pos, ref_a, alt_a) else {
            continue;
        };
        for i in 0..b.num_rows() {
            let gt_str = gt
                .and_then(|g| (!g.is_null(i)).then(|| g.value(i).to_string()))
                .unwrap_or_default();
            if !gt_str.chars().any(|c| c == '1') {
                continue;
            }
            let key: Key = (
                normalize_chrom(chrom.value(i)),
                pos.value(i),
                ref_a.value(i).to_string(),
                alt_a.value(i).to_string(),
            );
            if let Some(ann) = cv_map.get(&key) {
                hits.push(Hit {
                    gene: ann.gene_symbol.clone().unwrap_or_default(),
                    chrom: key.0.clone(),
                    pos: key.1,
                    ref_a: key.2.clone(),
                    alt_a: key.3.clone(),
                    genotype: gt_str,
                    significance: ann.clinical_significance.clone().unwrap_or_default(),
                    phenotype: ann.disease_name.clone().unwrap_or_default(),
                });
            }
        }
    }

    hits.sort_by(|a, b| {
        a.gene
            .cmp(&b.gene)
            .then_with(|| a.chrom.cmp(&b.chrom))
            .then_with(|| a.pos.cmp(&b.pos))
    });

    println!(
        "{:<10} {:<8} {:<12} {:<4} {:<4} {:<6} {:<24} drug / phenotype",
        "gene", "chrom", "pos", "ref", "alt", "gt", "significance"
    );
    println!("{}", "-".repeat(120));
    for h in &hits {
        println!(
            "{:<10} {:<8} {:<12} {:<4} {:<4} {:<6} {:<24} {}",
            trunc(&h.gene, 10),
            h.chrom,
            h.pos,
            trunc(&h.ref_a, 4),
            trunc(&h.alt_a, 4),
            h.genotype,
            trunc(&h.significance, 24),
            h.phenotype,
        );
    }
    println!("{}", "-".repeat(120));
    println!("{} PGx-relevant variants carried by {}", hits.len(), sample);
    println!(
        "NOTE: screening only - not a diplotype call, not a prescribing \
         recommendation. Review hits with a clinician against current CPIC \
         guidance (https://cpicpgx.org/)."
    );
    Ok(())
}

struct Hit {
    gene: String,
    chrom: String,
    pos: u64,
    ref_a: String,
    alt_a: String,
    genotype: String,
    significance: String,
    phenotype: String,
}

struct Annotation {
    gene_symbol: Option<String>,
    clinical_significance: Option<String>,
    disease_name: Option<String>,
}

fn build_clinvar_map(batches: &[RecordBatch]) -> HashMap<Key, Annotation> {
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
            let key: Key = (
                normalize_chrom(chrom.value(i)),
                pos.value(i),
                ref_a.value(i).to_string(),
                alt_a.value(i).to_string(),
            );
            map.insert(
                key,
                Annotation {
                    gene_symbol: gene.and_then(|a| (!a.is_null(i)).then(|| a.value(i).to_string())),
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
fn normalize_chrom(c: &str) -> String {
    c.strip_prefix("chr").unwrap_or(c).to_string()
}
fn sql_escape(s: &str) -> String {
    s.replace('\'', "''")
}
fn trunc(s: &str, n: usize) -> String {
    if s.chars().count() > n {
        let mut out: String = s.chars().take(n.saturating_sub(1)).collect();
        out.push('...');
        out
    } else {
        s.to_string()
    }
}
