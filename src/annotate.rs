//! `biolance annotate` — join a sample's variants against any annotation
//! VCF on-the-fly and emit configurable INFO fields as output columns.
//!
//! Unlike `join`, which reads the pre-ingested `clinvar` table, this
//! command streams the annotation VCF directly so it works with any
//! VCF — gnomAD, COSMIC, dbSNP, a local call set, anything. No ingest
//! or schema change needed.
//!
//! For large annotations (e.g. whole-genome gnomAD) pair with
//! `--chrom/--start/--end` to cap the region scanned. If the
//! annotation VCF has a tabix/CSI index noodles will still do a linear
//! scan here — position-indexed random-access is a future optimization.

use std::collections::HashMap;

use anyhow::{anyhow, Context, Result};
use arrow_array::{Array, RecordBatch, StringArray, UInt64Array};
use futures::TryStreamExt;
use lancedb::query::{ExecutableQuery, QueryBase};
use noodles::vcf::{
    self as vcf,
    variant::record::{info::field::Value as InfoValue, AlternateBases},
};

use crate::schema::VARIANTS_TABLE;
use crate::store::Store;

type Key = (String, u64, String, String);

pub async fn run(
    store_path: &str,
    sample_name: &str,
    annotation_vcf: &str,
    info_fields: &[String],
    chrom: Option<&str>,
    start: Option<u64>,
    end: Option<u64>,
) -> Result<()> {
    if info_fields.is_empty() {
        return Err(anyhow!(
            "--info is required (comma-separated INFO keys to extract)"
        ));
    }

    let store = Store::open(store_path).await?;
    let tables = store.conn.table_names().execute().await?;
    if !tables.iter().any(|n| n == VARIANTS_TABLE) {
        return Err(anyhow!("store has no '{VARIANTS_TABLE}' table"));
    }

    // 1. Pull sample variants for the region.
    let variants = store.conn.open_table(VARIANTS_TABLE).execute().await?;
    let mut preds = vec![format!("sample_name = '{}'", sql_escape(sample_name))];
    if let Some(c) = chrom {
        // Accept either "chr17" or "17" — match both forms in the store.
        let bare = c.strip_prefix("chr").unwrap_or(c);
        preds.push(format!(
            "(chrom = '{}' OR chrom = 'chr{}')",
            sql_escape(bare),
            sql_escape(bare)
        ));
    }
    if let Some(v) = start {
        preds.push(format!("pos >= {v}"));
    }
    if let Some(v) = end {
        preds.push(format!("pos <= {v}"));
    }
    let batches: Vec<RecordBatch> = variants
        .query()
        .only_if(preds.join(" AND "))
        .execute()
        .await?
        .try_collect()
        .await?;

    let carriers = build_carrier_map(&batches);
    eprintln!("[annotate] {} carrier sites in query region", carriers.len());
    if carriers.is_empty() {
        return Ok(());
    }

    // 2. Stream the annotation VCF, emit matched rows with requested INFO fields.
    let mut reader = vcf::io::reader::Builder::default()
        .build_from_path(annotation_vcf)
        .with_context(|| format!("opening {annotation_vcf}"))?;
    let header = reader
        .read_header()
        .with_context(|| format!("reading header of {annotation_vcf}"))?;

    // Print column header.
    let mut hdr = vec![
        "chrom".to_string(),
        "pos".to_string(),
        "ref".to_string(),
        "alt".to_string(),
        "gt".to_string(),
    ];
    for f in info_fields {
        hdr.push(f.clone());
    }
    println!("{}", hdr.join("\t"));

    let mut matched = 0usize;
    for result in reader.records() {
        let record = result?;
        let chrom_s = record.reference_sequence_name().to_string();
        let Some(pos_r) = record.variant_start() else {
            continue;
        };
        let pos_v = usize::from(pos_r?) as u64;
        // Quick region short-circuit (linear scan; annotation VCF is not indexed here).
        if let Some(c) = chrom {
            if normalize_chrom(&chrom_s) != normalize_chrom(c) {
                continue;
            }
        }
        if let Some(s) = start {
            if pos_v < s {
                continue;
            }
        }
        if let Some(e) = end {
            if pos_v > e {
                // If sorted by pos we could break, but the file may span multiple
                // chroms — a full skip only when chrom also matches the filter.
                if chrom.is_some() {
                    break;
                } else {
                    continue;
                }
            }
        }

        let ref_a = record.reference_bases().to_string();
        let alts: Vec<String> = record
            .alternate_bases()
            .iter()
            .filter_map(|a| a.ok().map(|s| s.to_string()))
            .collect();

        for (alt_idx, alt) in alts.iter().enumerate() {
            let key: Key = (
                normalize_chrom(&chrom_s),
                pos_v,
                ref_a.clone(),
                alt.clone(),
            );
            let Some(gt) = carriers.get(&key) else {
                continue;
            };
            let mut row = vec![key.0.clone(), key.1.to_string(), key.2.clone(), alt.clone(), gt.clone()];
            for field in info_fields {
                row.push(extract_info(&record, &header, field, alt_idx));
            }
            println!("{}", row.join("\t"));
            matched += 1;
        }
    }
    eprintln!("[annotate] {matched} annotated variants emitted");
    Ok(())
}

fn build_carrier_map(batches: &[RecordBatch]) -> HashMap<Key, String> {
    let mut map = HashMap::new();
    for b in batches {
        let chrom = col_str(b, "chrom");
        let pos = col_u64(b, "pos");
        let ref_a = col_str(b, "ref_allele");
        let alt_a = col_str(b, "alt_allele");
        let gt = col_str(b, "genotype");
        let (Some(chrom), Some(pos), Some(ref_a), Some(alt_a)) = (chrom, pos, ref_a, alt_a) else {
            continue;
        };
        for i in 0..b.num_rows() {
            let gt_s = gt
                .and_then(|g| (!g.is_null(i)).then(|| g.value(i).to_string()))
                .unwrap_or_default();
            // Skip sites where the sample isn't a carrier.
            if !gt_s.chars().any(|c| c == '1') {
                continue;
            }
            let key: Key = (
                normalize_chrom(chrom.value(i)),
                pos.value(i),
                ref_a.value(i).to_string(),
                alt_a.value(i).to_string(),
            );
            map.insert(key, gt_s);
        }
    }
    map
}

/// Render an INFO field value for the given ALT index. If the field is
/// `Number=A` or `Number=R` we pick the per-alt scalar; otherwise we
/// emit all values comma-joined. This heuristic works for the common
/// population-genetics fields (AF, AC, AN) without having to look up
/// the Number definition in the header.
fn extract_info(
    record: &vcf::Record,
    header: &vcf::Header,
    field: &str,
    alt_idx: usize,
) -> String {
    use vcf::variant::record::info::field::value::Array;
    let info = record.info();
    let Some(Ok(Some(value))) = info.get(header, field) else {
        return ".".to_string();
    };
    match value {
        InfoValue::String(s) => s.to_string(),
        InfoValue::Integer(n) => n.to_string(),
        InfoValue::Float(f) => f.to_string(),
        InfoValue::Flag => "1".to_string(),
        InfoValue::Character(c) => c.to_string(),
        InfoValue::Array(arr) => {
            let vals: Vec<String> = match arr {
                Array::Integer(it) => it
                    .iter()
                    .filter_map(|r| r.ok())
                    .map(|o| o.map(|v| v.to_string()).unwrap_or_else(|| ".".into()))
                    .collect(),
                Array::Float(it) => it
                    .iter()
                    .filter_map(|r| r.ok())
                    .map(|o| o.map(|v| v.to_string()).unwrap_or_else(|| ".".into()))
                    .collect(),
                Array::Character(it) => it
                    .iter()
                    .filter_map(|r| r.ok())
                    .map(|o| o.map(|v| v.to_string()).unwrap_or_else(|| ".".into()))
                    .collect(),
                Array::String(it) => it
                    .iter()
                    .filter_map(|r| r.ok())
                    .map(|o| o.map(|v| v.to_string()).unwrap_or_else(|| ".".into()))
                    .collect(),
            };
            if vals.is_empty() {
                ".".to_string()
            } else if vals.len() == 1 {
                vals.into_iter().next().unwrap()
            } else if alt_idx < vals.len() {
                // Ambiguous A vs R; per-alt (A-style) is the common case.
                vals[alt_idx].clone()
            } else {
                vals.join(",")
            }
        }
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
fn normalize_chrom(c: &str) -> String {
    c.strip_prefix("chr").unwrap_or(c).to_string()
}
fn sql_escape(s: &str) -> String {
    s.replace('\'', "''")
}
