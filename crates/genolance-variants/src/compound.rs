//! `genolance compound-het` - flag genes where a single sample carries
//! two or more pathogenic/likely-pathogenic heterozygous variants, or
//! at least one homozygous P/LP variant.
//!
//! Short-read WGS can't phase variants that are more than a read-pair
//! apart, so we can't distinguish *in cis* from *in trans*. This
//! command surfaces genes that are worth resolving via family phasing
//! or Sanger follow-up; it does NOT diagnose compound heterozygosity.

use std::collections::HashMap;

use anyhow::{anyhow, Result};
use arrow_array::{Array, RecordBatch};
use futures::TryStreamExt;
use lancedb::query::{ExecutableQuery, QueryBase};

use genolance_core::schema::{CLINVAR_TABLE, VARIANTS_TABLE};
use genolance_core::store::Store;

use crate::common::{
    build_clinvar_map, carrier_state, ensure_clinvar_unique, fmt_dp, fmt_qual, normalize_chrom,
    opt_string, req_f32, req_str, req_u32, req_u64, significance_predicate, sql_escape, truncate,
    validate_filters, Carrier, Gate, Key,
};
pub use crate::filters::Filters;

#[derive(Clone)]
struct Hit {
    chrom: String,
    pos: u64,
    ref_allele: String,
    alt_allele: String,
    significance: String,
    disease: String,
    gt: String,
    qual: Option<f32>,
    dp: Option<u32>,
}

/// Default substring needle when no filter is supplied.
const DEFAULT_SIGNIFICANCE: &str = "pathogenic";

pub async fn run(store_path: &str, sample: &str, f: &Filters) -> Result<()> {
    validate_filters(f)?;
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
        "[compound-het] {} ClinVar records in lookup (after filters)",
        cv_map.len()
    );
    if cv_map.is_empty() {
        return Ok(());
    }

    // 2. Pull the sample's variants (NULL QUAL/DP kept so they can be counted).
    let variants = store.variants.open_table(VARIANTS_TABLE).execute().await?;
    let mut pred = format!("sample_name = '{}'", sql_escape(sample));
    if let Some(q) = f.min_qual {
        pred.push_str(&format!(" AND (quality IS NULL OR quality >= {q})"));
    }
    if let Some(d) = f.min_dp {
        pred.push_str(&format!(" AND (read_depth IS NULL OR read_depth >= {d})"));
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
    let mut gate = Gate::new(f);
    let mut no_gene = 0usize;
    for b in &v_batches {
        let t = VARIANTS_TABLE;
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
            let Some(gene) = ann.gene_symbol.clone() else {
                no_gene += 1;
                continue;
            };
            by_gene.entry(gene).or_default().push(Hit {
                chrom: chrom.value(i).to_string(),
                pos: pos.value(i),
                ref_allele: ref_a.value(i).to_string(),
                alt_allele: alt_a.value(i).to_string(),
                significance: ann.clinical_significance.clone().unwrap_or_default(),
                disease: ann.disease_name.clone().unwrap_or_default(),
                gt: gt_s.unwrap_or_default(),
                qual: q,
                dp: d,
            });
        }
    }

    // 4. Report: genes with >=2 hets, or >=1 hom-alt / hemizygous-alt.
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
        gate.report("compound-het");
        report_no_gene(no_gene);
        return Ok(());
    }

    println!("Compound-het / homozygous P/LP screen - sample {sample}");
    println!("{}", "-".repeat(100));
    for (gene, hits) in &genes {
        let hets = hits.iter().filter(|h| is_het(&h.gt)).count();
        let homs = hits.iter().filter(|h| is_hom_alt(&h.gt)).count();
        let disposition = if homs > 0 {
            format!("HOMOZYGOUS/HEMIZYGOUS ({} call{})", homs, plural(homs))
        } else {
            format!("POSSIBLE COMPOUND HET ({} het P/LP - phase unknown)", hets)
        };
        println!("\n{} - {}", gene, disposition);
        for h in hits {
            println!(
                "  {}:{} {}>{}  GT={:<5} QUAL={:>5} DP={:>3}  [{}]  {}",
                h.chrom,
                h.pos,
                truncate(&h.ref_allele, 8),
                truncate(&h.alt_allele, 8),
                h.gt,
                fmt_qual(h.qual),
                fmt_dp(h.dp),
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
    gate.report("compound-het");
    report_no_gene(no_gene);
    Ok(())
}

fn report_no_gene(n: usize) {
    if n > 0 {
        eprintln!("[compound-het] {n} matching call(s) skipped because the ClinVar record has no gene symbol");
    }
}

fn build_clinvar_predicate(f: &Filters) -> Result<String> {
    significance_predicate(f, Some(DEFAULT_SIGNIFICANCE))
}

fn is_het(gt: &str) -> bool {
    // "0/1", "1/0", "0|1", "1|0", "1/2", etc. - any diploid call whose two
    // alleles differ.
    let alleles: Vec<&str> = gt.split(['/', '|']).collect();
    if alleles.len() != 2 {
        return false;
    }
    alleles[0] != alleles[1] && !alleles.contains(&".")
}
/// All called alleles are the ALT: `1/1`, `1|1`, or a haploid `1`
/// (hemizygous, e.g. chrX in males).
fn is_hom_alt(gt: &str) -> bool {
    let alleles: Vec<&str> = gt.split(['/', '|']).collect();
    !alleles.is_empty() && alleles.iter().all(|a| *a == "1")
}
fn plural(n: usize) -> &'static str {
    if n == 1 {
        ""
    } else {
        "s"
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hemizygous_alt_is_flagged_like_hom_alt() {
        assert!(is_hom_alt("1"));
        assert!(is_hom_alt("1/1"));
        assert!(is_hom_alt("1|1"));
        assert!(!is_hom_alt("0/1"));
        assert!(!is_hom_alt("0"));
        assert!(!is_hom_alt("./."));
    }

    #[test]
    fn het_needs_two_called_differing_alleles() {
        assert!(is_het("0/1"));
        assert!(is_het("1|0"));
        assert!(!is_het("1/1"));
        assert!(!is_het("1"));
        assert!(!is_het("./1"));
    }
}
