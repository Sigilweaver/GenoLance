use anyhow::{anyhow, Result};
use arrow_array::{Array, Float32Array, RecordBatch, StringArray, UInt32Array, UInt64Array};
use futures::TryStreamExt;
use lancedb::query::{ExecutableQuery, QueryBase, Select};

use genolance_core::schema::{CLINVAR_TABLE, VARIANTS_TABLE};
use genolance_core::store::Store;

/// Query variants from a GenoLance store.
///
/// If `gene` is provided the store must also contain a `clinvar` table —
/// ClinVar positions for that gene are looked up and used as the region
/// filter.
pub async fn run(
    store_path: &str,
    gene: Option<&str>,
    chrom: Option<&str>,
    start: Option<u64>,
    end: Option<u64>,
    output_format: &str,
) -> Result<()> {
    let store = Store::open(store_path).await?;

    let table_names = store.variants.table_names().execute().await?;
    if !table_names.iter().any(|n| n == VARIANTS_TABLE) {
        return Err(anyhow!(
            "store {store_path} has no '{VARIANTS_TABLE}' table — did you run `genolance ingest` first?"
        ));
    }
    let variants = store.variants.open_table(VARIANTS_TABLE).execute().await?;

    let mut predicates: Vec<String> = Vec::new();

    // Region filters from CLI args.
    if let Some(c) = chrom {
        predicates.push(format!("chrom = '{}'", sql_escape(c)));
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
            return Err(anyhow!(
                "--gene requires a '{CLINVAR_TABLE}' table; ingest clinvar.vcf.gz first"
            ));
        }
        let cv = store.variants.open_table(CLINVAR_TABLE).execute().await?;
        let stream = cv
            .query()
            .only_if(format!("gene_symbol = '{}'", sql_escape(g)))
            .select(Select::columns(&["chrom", "pos"]))
            .execute()
            .await?;
        let batches: Vec<RecordBatch> = stream.try_collect().await?;
        let positions = collect_chrom_pos(&batches);
        if positions.is_empty() {
            println!("(no ClinVar positions matched gene {g})");
            return Ok(());
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
        "json" => print_json(&batches)?,
        "arrow" => print_arrow(&batches),
        _ => print_table(&batches),
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
            let bare = c.trim_start_matches("chr");
            let prefixed = if c.starts_with("chr") {
                c.clone()
            } else {
                format!("chr{c}")
            };
            format!(
                "((chrom = '{}' OR chrom = '{}') AND pos IN ({}))",
                sql_escape(bare),
                sql_escape(&prefixed),
                in_list.join(",")
            )
        })
        .collect();
    format!("({})", clauses.join(" OR "))
}

// ---------------------------------------------------------------------------
// Output formatting
// ---------------------------------------------------------------------------

fn print_table(batches: &[RecordBatch]) {
    let mut total = 0usize;
    println!(
        "{:<16} {:<8} {:<12} {:<6} {:<6} {:<8} {:<10} {:<6}",
        "sample", "chrom", "pos", "ref", "alt", "qual", "gt", "dp"
    );
    println!("{}", "-".repeat(80));
    for b in batches {
        let sample = col_str(b, "sample_name");
        let chrom = col_str(b, "chrom");
        let pos = col_u64(b, "pos");
        let ref_a = col_str(b, "ref_allele");
        let alt_a = col_str(b, "alt_allele");
        let qual = col_f32(b, "quality");
        let gt = col_str(b, "genotype");
        let dp = col_u32(b, "read_depth");

        for i in 0..b.num_rows() {
            total += 1;
            if total > 100 {
                continue;
            }
            println!(
                "{:<16} {:<8} {:<12} {:<6} {:<6} {:<8} {:<10} {:<6}",
                sample.as_ref().map(|a| a.value(i)).unwrap_or(""),
                chrom.as_ref().map(|a| a.value(i)).unwrap_or(""),
                pos.as_ref()
                    .map(|a| a.value(i).to_string())
                    .unwrap_or_default(),
                truncate(ref_a.as_ref().map(|a| a.value(i)).unwrap_or(""), 6),
                truncate(alt_a.as_ref().map(|a| a.value(i)).unwrap_or(""), 6),
                fmt_opt_f32(qual.as_ref(), i),
                gt.as_ref().map(|a| a.value(i)).unwrap_or(""),
                fmt_opt_u32(dp.as_ref(), i),
            );
        }
    }
    println!("{}", "-".repeat(80));
    if total > 100 {
        println!("showing first 100 of {total} rows");
    } else {
        println!("{total} rows");
    }
}

fn print_json(batches: &[RecordBatch]) -> Result<()> {
    for b in batches {
        let sample = col_str(b, "sample_name");
        let chrom = col_str(b, "chrom");
        let pos = col_u64(b, "pos");
        let ref_a = col_str(b, "ref_allele");
        let alt_a = col_str(b, "alt_allele");
        let gt = col_str(b, "genotype");
        for i in 0..b.num_rows() {
            println!(
                "{{\"sample\":\"{}\",\"chrom\":\"{}\",\"pos\":{},\"ref\":\"{}\",\"alt\":\"{}\",\"gt\":\"{}\"}}",
                sample.as_ref().map(|a| a.value(i)).unwrap_or(""),
                chrom.as_ref().map(|a| a.value(i)).unwrap_or(""),
                pos.as_ref().map(|a| a.value(i)).unwrap_or(0),
                ref_a.as_ref().map(|a| a.value(i)).unwrap_or(""),
                alt_a.as_ref().map(|a| a.value(i)).unwrap_or(""),
                gt.as_ref().map(|a| a.value(i)).unwrap_or(""),
            );
        }
    }
    Ok(())
}

fn print_arrow(batches: &[RecordBatch]) {
    let rows: usize = batches.iter().map(|b| b.num_rows()).sum();
    let cols = batches.first().map(|b| b.num_columns()).unwrap_or(0);
    println!(
        "Arrow RecordBatches: {} ({} rows x {} cols)",
        batches.len(),
        rows,
        cols
    );
    if let Some(first) = batches.first() {
        println!("schema: {}", first.schema());
    }
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

fn fmt_opt_f32(a: Option<&&Float32Array>, i: usize) -> String {
    match a {
        Some(arr) if !arr.is_null(i) => format!("{:.1}", arr.value(i)),
        _ => String::from("."),
    }
}
fn fmt_opt_u32(a: Option<&&UInt32Array>, i: usize) -> String {
    match a {
        Some(arr) if !arr.is_null(i) => arr.value(i).to_string(),
        _ => String::from("."),
    }
}

fn truncate(s: &str, n: usize) -> String {
    if s.len() > n {
        format!("{}…", &s[..n.saturating_sub(1)])
    } else {
        s.to_string()
    }
}
