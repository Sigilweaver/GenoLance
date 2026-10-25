use arrow_array::{Array, RecordBatch, StringArray, UInt64Array};
use futures::TryStreamExt;
use genolance_core::{Error, Result};
use lancedb::query::{ExecutableQuery, QueryBase, Select};

use genolance_core::schema::{CLINVAR_TABLE, VARIANTS_TABLE};
use genolance_core::store::Store;

use crate::common::{
    chrom_forms, chrom_predicate, fmt_dp, fmt_qual, opt_string, req_f32, req_str, req_u32, req_u64,
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

    // Gene → positions via ClinVar.
    if let Some(g) = gene {
        if !table_names.iter().any(|n| n == CLINVAR_TABLE) {
            return Err(invalid!(
                "--gene requires a '{CLINVAR_TABLE}' table; ingest clinvar.vcf.gz first"
            ));
        }
        let cv = store.variants.open_table(CLINVAR_TABLE).execute().await?;
        let stream = cv
            .query()
            .only_if(format!(
                "lower(gene_symbol) = '{}'",
                sql_escape(&g.to_lowercase())
            ))
            .select(Select::columns(&["chrom", "pos"]))
            .execute()
            .await?;
        let batches: Vec<RecordBatch> = stream.try_collect().await?;
        let positions = collect_chrom_pos(&batches);
        if positions.is_empty() {
            return Err(invalid!(
                "no ClinVar rows for gene '{g}' (symbols match case-insensitively); check the spelling"
            ));
        }
        predicates.push(region_predicate(&positions));
    }

    let mut q = variants.query();
    if !predicates.is_empty() {
        q = q.only_if(predicates.join(" AND "));
    }
    let stream = q.execute().await?;
    let batches: Vec<RecordBatch> = stream.try_collect().await?;

    match output_format {
        OutputFormat::Json => print_json(&batches)?,
        OutputFormat::Table => print_table(&batches)?,
    }
    Ok(())
}

fn sql_escape(s: &str) -> String {
    s.replace('\'', "''")
}

fn collect_chrom_pos(batches: &[RecordBatch]) -> Vec<(String, u64)> {
    let mut out = Vec::new();
    for b in batches {
        let chrom = b
            .column_by_name("chrom")
            .and_then(|c| c.as_any().downcast_ref::<StringArray>());
        let pos = b
            .column_by_name("pos")
            .and_then(|c| c.as_any().downcast_ref::<UInt64Array>());
        if let (Some(chrom), Some(pos)) = (chrom, pos) {
            for i in 0..b.num_rows() {
                out.push((chrom.value(i).to_string(), pos.value(i)));
            }
        }
    }
    out.sort();
    out.dedup();
    out
}

/// Build a SQL predicate for `(chrom, pos)` pairs.
///
/// To keep the predicate compact we group positions by chromosome. ClinVar
/// uses unprefixed chromosome names ("1") while most sample VCFs use
/// "chr1"; we match either form.
fn region_predicate(positions: &[(String, u64)]) -> String {
    use std::collections::BTreeMap;
    let mut by_chrom: BTreeMap<String, Vec<u64>> = BTreeMap::new();
    for (c, p) in positions {
        by_chrom.entry(c.clone()).or_default().push(*p);
    }
    let clauses: Vec<String> = by_chrom
        .into_iter()
        .map(|(c, ps)| {
            let in_list: Vec<String> = ps.iter().map(|p| p.to_string()).collect();
            let names: Vec<String> = chrom_forms(&c)
                .iter()
                .map(|f| format!("chrom = '{}'", sql_escape(f)))
                .collect();
            format!(
                "(({}) AND pos IN ({}))",
                names.join(" OR "),
                in_list.join(",")
            )
        })
        .collect();
    format!("({})", clauses.join(" OR "))
}

// ---------------------------------------------------------------------------
// Output formatting
// ---------------------------------------------------------------------------

fn print_table(batches: &[RecordBatch]) -> Result<()> {
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

        for i in 0..b.num_rows() {
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
fn print_json(batches: &[RecordBatch]) -> Result<()> {
    for b in batches {
        let t = VARIANTS_TABLE;
        let sample = req_str(b, t, "sample_name")?;
        let chrom = req_str(b, t, "chrom")?;
        let pos = req_u64(b, t, "pos")?;
        let ref_a = req_str(b, t, "ref_allele")?;
        let alt_a = req_str(b, t, "alt_allele")?;
        let gt = req_str(b, t, "genotype")?;
        for i in 0..b.num_rows() {
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

fn truncate(s: &str, n: usize) -> String {
    if s.chars().count() > n {
        let cut: String = s.chars().take(n.saturating_sub(1)).collect();
        format!("{cut}...")
    } else {
        s.to_string()
    }
}
