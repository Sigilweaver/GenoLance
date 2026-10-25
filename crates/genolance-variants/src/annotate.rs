//! `genolance annotate` - join a sample's variants against any annotation
//! VCF on-the-fly and emit configurable INFO fields as output columns.
//!
//! Unlike `join`, which reads the pre-ingested `clinvar` table, this
//! command streams the annotation VCF directly so it works with any
//! VCF - gnomAD, COSMIC, dbSNP, a local call set, anything. No ingest
//! or schema change needed.
//!
//! For large annotations (e.g. whole-genome gnomAD) pair with
//! `--chrom/--start/--end` to cap the region scanned. If the
//! annotation VCF has a tabix/CSI index noodles will still do a linear
//! scan here - position-indexed random-access is a future optimization.

use std::collections::HashMap;

use arrow_array::RecordBatch;
use futures::TryStreamExt;
use genolance_core::{Error, Result, ResultExt};
use lancedb::query::{ExecutableQuery, QueryBase};
use noodles::vcf::{
    self as vcf,
    variant::record::{info::field::Value as InfoValue, AlternateBases},
};

use genolance_core::schema::VARIANTS_TABLE;
use genolance_core::store::Store;

use crate::common::{
    carrier_state, chrom_predicate, normalize_alleles, normalize_chrom, opt_string, req_str,
    req_u64, require_samples, sql_escape, Carrier,
};

type Key = (String, u64, String, String);

const VARIANTS: &str = VARIANTS_TABLE;

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
        return Err(invalid!(
            "--info is required (comma-separated INFO keys to extract)"
        ));
    }

    let store = Store::open(store_path).await?;
    let tables = store.variants.table_names().execute().await?;
    if !tables.iter().any(|n| n == VARIANTS_TABLE) {
        return Err(Error::missing_table(VARIANTS_TABLE, "ingest it first"));
    }

    require_samples(&store, std::slice::from_ref(&sample_name.to_string())).await?;

    // Open the annotation file first: a truncated or wrongly compressed file,
    // or an --info key it does not declare, must fail before any output.
    let (mut reader, tail) = crate::ingest::open_vcf_reader(annotation_vcf)
        .with_context(|| format!("opening {annotation_vcf}"))?;
    let header = reader
        .read_header()
        .with_context(|| format!("reading header of {annotation_vcf}"))?;
    let undeclared: Vec<&String> = info_fields
        .iter()
        .filter(|f| !header.infos().contains_key(f.as_str()))
        .collect();
    if !undeclared.is_empty() {
        return Err(Error::unrecognized_input(
            annotation_vcf,
            format!(
                "does not declare INFO key(s) {}; declared keys: {}",
                undeclared
                    .iter()
                    .map(|s| s.as_str())
                    .collect::<Vec<_>>()
                    .join(", "),
                header
                    .infos()
                    .keys()
                    .map(String::as_str)
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
        ));
    }

    // The annotation file must be on the sample's assembly.
    {
        let ann = crate::ingest::header_assembly(&header, None)?;
        let reg = crate::common::load_registry(&store).await?;
        let sample_assembly = reg
            .iter()
            .find(|r| r.name == sample_name)
            .and_then(|r| r.assembly);
        match (sample_assembly, ann.assembly) {
            (Some(a), Some(b)) if a != b => {
                return Err(Error::reference_mismatch(
                    sample_name,
                    format!("the sample is {a} but {annotation_vcf} is {b}; use an annotation file for {a}"),
                ));
            }
            (Some(_), Some(_)) => {}
            (sa, aa) => eprintln!(
                "[annotate] WARNING: cannot verify the assembly (sample: {}, {annotation_vcf}: {}); positions are assumed to match",
                sa.map(|a| a.to_string()).unwrap_or_else(|| "unknown".into()),
                aa.map(|a| a.to_string()).unwrap_or_else(|| "unknown".into()),
            ),
        }
    }

    // 1. Pull sample variants for the region.
    let variants = store.variants.open_table(VARIANTS_TABLE).execute().await?;
    let mut preds = vec![format!("sample_name = '{}'", sql_escape(sample_name))];
    if let Some(c) = chrom {
        // Accept "chr17" or "17" (and chrM, M, MT) - match every form in the store.
        preds.push(chrom_predicate(c));
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

    let (carriers, no_calls) = build_carrier_map(&batches)?;
    eprintln!(
        "[annotate] {} carrier sites in query region",
        carriers.len()
    );
    if no_calls > 0 {
        eprintln!(
            "[annotate] {no_calls} no-call row(s) in the region were not treated as carriers"
        );
    }
    if carriers.is_empty() {
        return Ok(());
    }

    // 2. Stream the annotation VCF, emit matched rows with requested INFO fields.

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
    for (n, result) in reader.records().enumerate() {
        let record =
            result.with_context(|| format!("{annotation_vcf}: reading record {}", n + 1))?;
        let chrom_s = record.reference_sequence_name().to_string();
        let Some(pos_r) = record.variant_start() else {
            continue;
        };
        let pos_v = usize::from(
            pos_r.with_context(|| format!("{annotation_vcf}: record {}: invalid POS", n + 1))?,
        ) as u64;
        // Region filter. The annotation VCF is scanned linearly and in full:
        // stopping at the first record past --end would silently skip
        // records of an unsorted or concatenated file.
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
                continue;
            }
        }

        let ref_a = record.reference_bases().to_string();
        let mut alts: Vec<String> = Vec::new();
        for a in record.alternate_bases().iter() {
            alts.push(
                a.with_context(|| format!("{annotation_vcf}: {chrom_s}:{pos_v}: invalid ALT"))?
                    .to_string(),
            );
        }

        for (alt_idx, alt) in alts.iter().enumerate() {
            let (npos, nref, nalt) = normalize_alleles(pos_v, &ref_a, alt);
            let key: Key = (normalize_chrom(&chrom_s), npos, nref, nalt);
            let Some(gt) = carriers.get(&key) else {
                continue;
            };
            let mut row = vec![
                normalize_chrom(&chrom_s),
                pos_v.to_string(),
                ref_a.clone(),
                alt.clone(),
                gt.clone(),
            ];
            for field in info_fields {
                row.push(
                    extract_info(&record, &header, field, alt_idx, alts.len())
                        .with_context(|| format!("{annotation_vcf}: {chrom_s}:{pos_v}"))?,
                );
            }
            println!("{}", row.join("\t"));
            matched += 1;
        }
    }
    tail.check_complete(annotation_vcf)?;
    eprintln!("[annotate] {matched} annotated variants emitted");
    Ok(())
}

fn build_carrier_map(batches: &[RecordBatch]) -> Result<(HashMap<Key, String>, usize)> {
    let mut map = HashMap::new();
    let mut no_calls = 0usize;
    for b in batches {
        let chrom = req_str(b, VARIANTS, "chrom")?;
        let npos = req_u64(b, VARIANTS, "norm_pos")?;
        let nref = req_str(b, VARIANTS, "norm_ref")?;
        let nalt = req_str(b, VARIANTS, "norm_alt")?;
        let gt = req_str(b, VARIANTS, "genotype")?;
        for i in 0..b.num_rows() {
            let gt_s = opt_string(gt, i);
            match carrier_state(gt_s.as_deref()) {
                Carrier::Yes => {}
                Carrier::No => continue,
                Carrier::NoCall => {
                    no_calls += 1;
                    continue;
                }
            }
            let key: Key = (
                normalize_chrom(chrom.value(i)),
                npos.value(i),
                nref.value(i).to_string(),
                nalt.value(i).to_string(),
            );
            map.insert(key, gt_s.unwrap_or_default());
        }
    }
    Ok((map, no_calls))
}

/// Render an INFO field for the given ALT, using the Number declared in the
/// header: `A` picks the value for this ALT, `R` picks the value after the
/// REF value, anything else (`.`, a fixed count, `G`) emits every value.
/// A per-allele field whose value count does not match the ALT list, or an
/// element that does not parse, is an error.
fn extract_info(
    record: &vcf::Record,
    header: &vcf::Header,
    field: &str,
    alt_idx: usize,
    alt_count: usize,
) -> Result<String> {
    use vcf::header::record::value::map::info::Number;
    use vcf::variant::record::info::field::value::Array;
    let chrom = record.reference_sequence_name().to_string();
    let pos = record
        .variant_start()
        .and_then(|p| p.ok())
        .map(|p| usize::from(p) as u64)
        .unwrap_or(0);
    let bad = |reason: String| -> Error {
        Error::malformed_record(chrom.clone(), pos, format!("INFO/{field}"), reason)
    };
    let info = record.info();
    let value = match info.get(header, field) {
        None => return Ok(".".to_string()),
        Some(r) => r.map_err(|e| bad(e.to_string()))?,
    };
    let Some(value) = value else {
        return Ok(".".to_string());
    };
    let vals: Vec<String> = match value {
        InfoValue::String(s) => vec![s.to_string()],
        InfoValue::Integer(n) => vec![n.to_string()],
        InfoValue::Float(f) => vec![f.to_string()],
        InfoValue::Flag => vec!["1".to_string()],
        InfoValue::Character(c) => vec![c.to_string()],
        InfoValue::Array(arr) => {
            fn collect<T: ToString>(
                it: impl Iterator<Item = std::io::Result<Option<T>>>,
            ) -> std::io::Result<Vec<String>> {
                it.map(|r| r.map(|o| o.map(|v| v.to_string()).unwrap_or_else(|| ".".into())))
                    .collect()
            }
            let r = match arr {
                Array::Integer(it) => collect(it.iter()),
                Array::Float(it) => collect(it.iter()),
                Array::Character(it) => collect(it.iter()),
                Array::String(it) => collect(it.iter()),
            };
            r.map_err(|e| bad(e.to_string()))?
        }
    };
    let number = header.infos().get(field).map(|d| d.number());
    match number {
        Some(Number::AlternateBases) => {
            if vals.len() != alt_count {
                return Err(bad(format!(
                    "Number=A but {} value(s) for {alt_count} ALT allele(s)",
                    vals.len()
                )));
            }
            Ok(vals[alt_idx].clone())
        }
        Some(Number::ReferenceAlternateBases) => {
            if vals.len() != alt_count + 1 {
                return Err(bad(format!(
                    "Number=R but {} value(s) for REF plus {alt_count} ALT allele(s)",
                    vals.len()
                )));
            }
            Ok(vals[alt_idx + 1].clone())
        }
        _ => Ok(vals.join(",")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use arrow_array::{StringArray, UInt64Array};
    use arrow_schema::{DataType, Field, Schema};
    use std::sync::Arc;

    /// Minimal variants-shaped batch covering only the columns `build_carrier_map` reads.
    fn variants_batch(rows: &[(&str, u64, &str, &str, &str)]) -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("pos", DataType::UInt64, false),
            Field::new("ref_allele", DataType::Utf8, false),
            Field::new("alt_allele", DataType::Utf8, false),
            Field::new("genotype", DataType::Utf8, true),
            Field::new("norm_pos", DataType::UInt64, false),
            Field::new("norm_ref", DataType::Utf8, false),
            Field::new("norm_alt", DataType::Utf8, false),
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
                Arc::new(UInt64Array::from(
                    rows.iter().map(|r| r.1).collect::<Vec<_>>(),
                )),
                Arc::new(StringArray::from(
                    rows.iter().map(|r| r.2).collect::<Vec<_>>(),
                )),
                Arc::new(StringArray::from(
                    rows.iter().map(|r| r.3).collect::<Vec<_>>(),
                )),
            ],
        )
        .unwrap()
    }

    #[test]
    fn build_carrier_map_normalizes_chrom_and_keys_by_variant() {
        let batch = variants_batch(&[("chr17", 100, "A", "T", "0/1")]);
        let (map, _) = build_carrier_map(&[batch]).unwrap();
        assert_eq!(map.len(), 1);
        assert_eq!(
            map.get(&("17".to_string(), 100, "A".to_string(), "T".to_string()))
                .map(String::as_str),
            Some("0/1")
        );
    }

    #[test]
    fn build_carrier_map_skips_homozygous_reference_calls() {
        let batch = variants_batch(&[
            ("1", 100, "A", "T", "0/0"),
            ("1", 200, "C", "G", "0/1"),
            ("1", 300, "G", "A", "./."),
        ]);
        let (map, _) = build_carrier_map(&[batch]).unwrap();
        assert_eq!(map.len(), 1);
        assert!(map.contains_key(&("1".to_string(), 200, "C".to_string(), "G".to_string())));
    }
}
