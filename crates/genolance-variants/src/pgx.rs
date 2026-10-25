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

use arrow_array::{Array, RecordBatch};
use futures::TryStreamExt;
use genolance_core::{Error, Result};
use lancedb::query::{ExecutableQuery, QueryBase, Select};

use genolance_core::schema::{CLINVAR_TABLE, VARIANTS_TABLE};
use genolance_core::store::Store;

use crate::common::{
    build_clinvar_map, carrier_state, check_reference, ensure_clinvar_unique, gene_match_predicate,
    normalize_chrom, opt_string, req_str, req_u64, sql_escape, Carrier, Gate, Key,
};
use crate::filters::Filters;

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

/// Screen `sample` against ClinVar drug-response annotations. `extra_genes`
/// are added to the default list (case-insensitive); a gene with no ClinVar
/// rows at all is an error. Calls whose FILTER is not PASS are excluded and
/// counted unless `f.include_filtered` is set.
pub async fn run(
    store_path: &str,
    sample: &str,
    extra_genes: &[String],
    f: &Filters,
) -> Result<()> {
    let store = Store::open(store_path).await?;
    let tables = store.variants.table_names().execute().await?;
    if !tables.iter().any(|n| n == VARIANTS_TABLE) {
        return Err(Error::missing_table(VARIANTS_TABLE, "ingest it first"));
    }
    if !tables.iter().any(|n| n == CLINVAR_TABLE) {
        return Err(Error::missing_table(
            CLINVAR_TABLE,
            "ingest clinvar.vcf.gz first",
        ));
    }

    let mut genes: Vec<String> = DEFAULT_PGX_GENES.iter().map(|s| s.to_string()).collect();
    let mut requested: Vec<String> = Vec::new();
    for g in extra_genes {
        let g = g.trim().to_uppercase();
        if g.is_empty() {
            return Err(Error::InvalidArgument(
                "--genes contains an empty gene symbol".into(),
            ));
        }
        if !genes.iter().any(|x| x.eq_ignore_ascii_case(&g)) {
            genes.push(g.clone());
        }
        requested.push(g);
    }

    // 1. Pull ClinVar entries for PGx genes, filtered to drug-response rows.
    let clinvar = store.variants.open_table(CLINVAR_TABLE).execute().await?;
    ensure_clinvar_unique(&clinvar).await?;
    check_reference(
        &store,
        Some(std::slice::from_ref(&sample.to_string())),
        f,
        "pgx",
    )
    .await?;
    if !requested.is_empty() {
        let found: Vec<RecordBatch> = clinvar
            .query()
            .only_if(gene_match_predicate(&requested)?)
            .select(Select::columns(&["gene_symbols"]))
            .execute()
            .await?
            .try_collect()
            .await?;
        let mut present: std::collections::HashSet<String> = std::collections::HashSet::new();
        for b in &found {
            let col = crate::common::req_list(b, CLINVAR_TABLE, "gene_symbols")?;
            for i in 0..col.len() {
                for g in crate::common::list_strings(col, i) {
                    present.insert(g.to_uppercase());
                }
            }
        }
        let missing: Vec<&String> = requested.iter().filter(|g| !present.contains(*g)).collect();
        if !missing.is_empty() {
            return Err(Error::unknown_gene(
                missing.iter().map(|g| g.to_string()).collect(),
            ));
        }
    }
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
    let pred = format!("{} AND ({needle_pred})", gene_match_predicate(&genes)?);
    let cv_batches: Vec<RecordBatch> = clinvar
        .query()
        .only_if(pred)
        .execute()
        .await?
        .try_collect()
        .await?;
    let cv_map = build_clinvar_map(&cv_batches, &Filters::default())?;
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
    let mut gate = Gate::new(f);
    for b in &v_batches {
        let t = VARIANTS_TABLE;
        let chrom = req_str(b, t, "chrom")?;
        let npos = req_u64(b, t, "norm_pos")?;
        let nref = req_str(b, t, "norm_ref")?;
        let nalt = req_str(b, t, "norm_alt")?;
        let gt = req_str(b, t, "genotype")?;
        let filter = req_str(b, t, "filter")?;
        for i in 0..b.num_rows() {
            let gt_str = opt_string(gt, i);
            let key: Key = (
                normalize_chrom(chrom.value(i)),
                npos.value(i),
                nref.value(i).to_string(),
                nalt.value(i).to_string(),
            );
            let Some(ann) = cv_map.get(&key) else {
                continue;
            };
            match carrier_state(gt_str.as_deref()) {
                Carrier::Yes => {}
                Carrier::No => continue,
                Carrier::NoCall => {
                    gate.no_call += 1;
                    continue;
                }
            }
            let fl = opt_string(filter, i);
            if !gate.admit(fl.as_deref(), None, None) {
                continue;
            }
            hits.push(Hit {
                gene: ann.gene_symbol.clone().unwrap_or_default(),
                chrom: key.0.clone(),
                pos: key.1,
                ref_a: key.2.clone(),
                alt_a: key.3.clone(),
                genotype: gt_str.unwrap_or_default(),
                significance: ann.clinical_significance.clone().unwrap_or_default(),
                phenotype: ann.disease_name.clone().unwrap_or_default(),
            });
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
    gate.report("pgx");
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

fn trunc(s: &str, n: usize) -> String {
    if s.chars().count() > n {
        let mut out: String = s.chars().take(n.saturating_sub(1)).collect();
        out.push_str("...");
        out
    } else {
        s.to_string()
    }
}
