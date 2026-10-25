use std::collections::HashSet;

use arrow_array::{Array, RecordBatch};
use futures::TryStreamExt;
use genolance_core::{Error, Result};
use lancedb::query::{ExecutableQuery, QueryBase, Select};

use genolance_core::schema::{CLINVAR_TABLE, VARIANTS_TABLE};
use genolance_core::store::Store;

use crate::common::{
    chrom_predicate, fmt_dp, fmt_qual, gene_match_predicate, normalize_chrom, opt_string, req_f32,
    req_str, req_u32, req_u64,
};

/// Output format of `genolance query`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum OutputFormat {
    /// Human-readable table (first 100 rows).
    #[default]
    Table,
    /// One JSON object per line (NDJSON), one per row.
    Json,
}

impl std::str::FromStr for OutputFormat {
    type Err = Error;
    fn from_str(s: &str) -> std::result::Result<Self, Error> {
        match s {
            "table" => Ok(Self::Table),
            "json" => Ok(Self::Json),
            other => Err(Error::InvalidArgument(format!(
                "unknown output format '{other}'; expected table or json"
            ))),
        }
    }
}

/// Query variants from a GenoLance store.
///
/// If `gene` is provided the store must also contain a `clinvar` table -
/// ClinVar positions for that gene are looked up and used as the region
/// filter.
pub async fn run(
    store_path: &str,
    gene: Option<&str>,
    chrom: Option<&str>,
    start: Option<u64>,
    end: Option<u64>,
    output_format: OutputFormat,
) -> Result<()> {
    let store = Store::open(store_path).await?;

    let table_names = store.variants.table_names().execute().await?;
    if !table_names.iter().any(|n| n == VARIANTS_TABLE) {
        return Err(invalid!(
            "store {store_path} has no '{VARIANTS_TABLE}' table - did you run `genolance ingest` first?"
        ));
    }
    let variants = store.variants.open_table(VARIANTS_TABLE).execute().await?;

    let mut predicates: Vec<String> = Vec::new();

    // Region filters from CLI args.
    if let Some(c) = chrom {
        predicates.push(chrom_predicate(c));
    }
    if let Some(s) = start {
        predicates.push(format!("pos >= {s}"));
    }
    if let Some(e) = end {
        predicates.push(format!("pos <= {e}"));
    }

    // Gene -> positions via ClinVar (normalized positions, so a left-trimmed
    // call matches ClinVar's representation).
    let mut members: Option<HashSet<(String, u64)>> = None;
    if let Some(g) = gene {
        if !table_names.iter().any(|n| n == CLINVAR_TABLE) {
            return Err(Error::missing_table(
                CLINVAR_TABLE,
                "--gene needs ClinVar; ingest clinvar.vcf.gz first",
            ));
        }
        let cv = store.variants.open_table(CLINVAR_TABLE).execute().await?;
        let stream = cv
            .query()
            .only_if(gene_match_predicate(&[g.to_string()])?)
            .select(Select::columns(&["chrom", "norm_pos"]))
            .execute()
            .await?;
        let batches: Vec<RecordBatch> = stream.try_collect().await?;
        let positions = collect_chrom_pos(&batches)?;
        if positions.is_empty() {
            return Err(Error::unknown_gene(vec![g.to_string()]));
        }
        predicates.push(region_predicate(&positions));
        members = Some(
            positions
                .iter()
                .map(|(c, p)| (normalize_chrom(c), *p))
                .collect(),
        );
    }

    let mut q = variants.query();
    if !predicates.is_empty() {
        q = q.only_if(predicates.join(" AND "));
    }
    let stream = q.execute().await?;
    let batches: Vec<RecordBatch> = stream.try_collect().await?;

    match output_format {
        OutputFormat::Json => print_json(&batches, members.as_ref())?,
        OutputFormat::Table => print_table(&batches, members.as_ref())?,
    }
    Ok(())
}

fn collect_chrom_pos(batches: &[RecordBatch]) -> Result<Vec<(String, u64)>> {
    let mut out = Vec::new();
    for b in batches {
        let chrom = req_str(b, CLINVAR_TABLE, "chrom")?;
        let pos = req_u64(b, CLINVAR_TABLE, "norm_pos")?;
        for i in 0..b.num_rows() {
            out.push((chrom.value(i).to_string(), pos.value(i)));
        }
    }
    out.sort();
    out.dedup();
    Ok(out)
}

/// A coarse SQL predicate for the gene: per chromosome, the span from the
/// first to the last ClinVar position. The exact positions are checked while
/// printing, so a very large gene does not become a huge `IN` list.
fn region_predicate(positions: &[(String, u64)]) -> String {
    use std::collections::BTreeMap;
    let mut by_chrom: BTreeMap<String, (u64, u64)> = BTreeMap::new();
    for (c, p) in positions {
        let e = by_chrom.entry(normalize_chrom(c)).or_insert((*p, *p));
        e.0 = e.0.min(*p);
        e.1 = e.1.max(*p);
    }
    let clauses: Vec<String> = by_chrom
        .into_iter()
        .map(|(c, (lo, hi))| {
            format!(
                "({} AND norm_pos BETWEEN {lo} AND {hi})",
                chrom_predicate(&c)
            )
        })
        .collect();
    format!("({})", clauses.join(" OR "))
}

// ---------------------------------------------------------------------------
// Output formatting
// ---------------------------------------------------------------------------

fn print_table(batches: &[RecordBatch], members: Option<&HashSet<(String, u64)>>) -> Result<()> {
    let mut total = 0usize;
    println!(
        "{:<16} {:<8} {:<12} {:<6} {:<6} {:<8} {:<10} {:<6}",
        "sample", "chrom", "pos", "ref", "alt", "qual", "gt", "dp"
    );
    println!("{}", "-".repeat(80));
    for b in batches {
        let t = VARIANTS_TABLE;
        let sample = req_str(b, t, "sample_name")?;
        let chrom = req_str(b, t, "chrom")?;
        let pos = req_u64(b, t, "pos")?;
        let ref_a = req_str(b, t, "ref_allele")?;
        let alt_a = req_str(b, t, "alt_allele")?;
        let qual = req_f32(b, t, "quality")?;
        let gt = req_str(b, t, "genotype")?;
        let dp = req_u32(b, t, "read_depth")?;
        let npos = req_u64(b, t, "norm_pos")?;

        for i in 0..b.num_rows() {
            if !is_member(members, chrom.value(i), npos.value(i)) {
                continue;
            }
            total += 1;
            if total > 100 {
                continue;
            }
            println!(
                "{:<16} {:<8} {:<12} {:<6} {:<6} {:<8} {:<10} {:<6}",
                sample.value(i),
                chrom.value(i),
                pos.value(i),
                truncate(ref_a.value(i), 6),
                truncate(alt_a.value(i), 6),
                fmt_qual((!qual.is_null(i)).then(|| qual.value(i))),
                opt_string(gt, i).unwrap_or_else(|| ".".into()),
                fmt_dp((!dp.is_null(i)).then(|| dp.value(i))),
            );
        }
    }
    println!("{}", "-".repeat(80));
    if total > 100 {
        println!("showing first 100 of {total} rows");
    } else {
        println!("{total} rows");
    }
    Ok(())
}

/// One JSON object per row, properly escaped; a missing genotype is `null`.
fn print_json(batches: &[RecordBatch], members: Option<&HashSet<(String, u64)>>) -> Result<()> {
    for b in batches {
        let t = VARIANTS_TABLE;
        let sample = req_str(b, t, "sample_name")?;
        let chrom = req_str(b, t, "chrom")?;
        let pos = req_u64(b, t, "pos")?;
        let ref_a = req_str(b, t, "ref_allele")?;
        let alt_a = req_str(b, t, "alt_allele")?;
        let gt = req_str(b, t, "genotype")?;
        let npos = req_u64(b, t, "norm_pos")?;
        for i in 0..b.num_rows() {
            if !is_member(members, chrom.value(i), npos.value(i)) {
                continue;
            }
            let row = serde_json::json!({
                "sample": sample.value(i),
                "chrom": chrom.value(i),
                "pos": pos.value(i),
                "ref": ref_a.value(i),
                "alt": alt_a.value(i),
                "gt": opt_string(gt, i),
            });
            println!("{row}");
        }
    }
    Ok(())
}

fn is_member(members: Option<&HashSet<(String, u64)>>, chrom: &str, npos: u64) -> bool {
    members.is_none_or(|m| m.contains(&(normalize_chrom(chrom), npos)))
}

fn truncate(s: &str, n: usize) -> String {
    if s.chars().count() > n {
        let cut: String = s.chars().take(n.saturating_sub(1)).collect();
        format!("{cut}...")
    } else {
        s.to_string()
    }
}
