use std::collections::{HashMap, HashSet};

use anyhow::{anyhow, Result};
use arrow_array::{Array, RecordBatch, StringArray, UInt64Array};
use futures::TryStreamExt;
use lancedb::query::{ExecutableQuery, QueryBase};

use biolance_core::schema::VARIANTS_TABLE;
use biolance_core::store::Store;

type Key = (String, u64, String, String);

/// Compare variants across two or more samples.
pub async fn run(store_path: &str, samples: &[String], mode: &str) -> Result<()> {
    if samples.len() < 2 {
        return Err(anyhow!("compare needs at least 2 samples"));
    }

    let store = Store::open(store_path).await?;
    let tables = store.variants.table_names().execute().await?;
    if !tables.iter().any(|n| n == VARIANTS_TABLE) {
        return Err(anyhow!(
            "store has no '{VARIANTS_TABLE}' table; ingest samples first"
        ));
    }
    let table = store.variants.open_table(VARIANTS_TABLE).execute().await?;

    // Pull only the relevant samples using an IN-list.
    let in_list = samples
        .iter()
        .map(|s| format!("'{}'", s.replace('\'', "''")))
        .collect::<Vec<_>>()
        .join(",");
    let batches: Vec<RecordBatch> = table
        .query()
        .only_if(format!("sample_name IN ({in_list})"))
        .execute()
        .await?
        .try_collect()
        .await?;

    // Index: variant key -> map of sample_name -> genotype
    let mut index: HashMap<Key, HashMap<String, String>> = HashMap::new();
    for b in &batches {
        let sample = col_str(b, "sample_name");
        let chrom = col_str(b, "chrom");
        let pos = col_u64(b, "pos");
        let ref_a = col_str(b, "ref_allele");
        let alt_a = col_str(b, "alt_allele");
        let gt = col_str(b, "genotype");
        let (Some(sample), Some(chrom), Some(pos), Some(ref_a), Some(alt_a)) =
            (sample, chrom, pos, ref_a, alt_a)
        else {
            continue;
        };
        for i in 0..b.num_rows() {
            let s = sample.value(i).to_string();
            let g = gt
                .and_then(|g| (!g.is_null(i)).then(|| g.value(i).to_string()))
                .unwrap_or_else(|| "./.".to_string());
            // Skip homozygous-reference calls — they're not really variants.
            if g == "0/0" || g == "0|0" {
                continue;
            }
            let key: Key = (
                chrom.value(i).to_string(),
                pos.value(i),
                ref_a.value(i).to_string(),
                alt_a.value(i).to_string(),
            );
            index.entry(key).or_default().insert(s, g);
        }
    }

    let all_samples: HashSet<&str> = samples.iter().map(String::as_str).collect();

    match mode {
        "concordance" => print_concordance(&index, samples),
        "carrier-screen" => print_carrier_screen(&index, &all_samples, samples),
        "private" => print_private(&index, samples),
        other => return Err(anyhow!("unknown compare mode: {other}")),
    }
    Ok(())
}

fn print_concordance(index: &HashMap<Key, HashMap<String, String>>, samples: &[String]) {
    let mut shared_all = 0usize;
    let mut private_counts: HashMap<&str, usize> = HashMap::new();
    let mut shared_pair: HashMap<(&str, &str), usize> = HashMap::new();
    let total = index.len();

    for (_k, per_sample) in index {
        if per_sample.len() == samples.len() {
            shared_all += 1;
        } else if per_sample.len() == 1 {
            if let Some(s) = per_sample.keys().next() {
                *private_counts
                    .entry(s.as_str())
                    .or_insert(0) += 1;
            }
        }
        for i in 0..samples.len() {
            for j in (i + 1)..samples.len() {
                if per_sample.contains_key(&samples[i]) && per_sample.contains_key(&samples[j]) {
                    *shared_pair
                        .entry((samples[i].as_str(), samples[j].as_str()))
                        .or_insert(0) += 1;
                }
            }
        }
    }

    println!("Concordance across {} samples", samples.len());
    println!("  total union of variant sites : {total}");
    println!("  shared by all                : {shared_all}");
    for s in samples {
        let n = private_counts.get(s.as_str()).copied().unwrap_or(0);
        println!("  private to {:<16}     : {n}", s);
    }
    if samples.len() > 2 {
        println!("  pairwise shared:");
        for ((a, b), n) in shared_pair {
            println!("    {a} ∩ {b} = {n}");
        }
    } else if samples.len() == 2 {
        let key = (samples[0].as_str(), samples[1].as_str());
        let n = shared_pair.get(&key).copied().unwrap_or(0);
        println!(
            "  {} ∩ {} = {} (Jaccard ≈ {:.4})",
            samples[0],
            samples[1],
            n,
            jaccard(n, total)
        );
    }
}

fn jaccard(intersection: usize, union: usize) -> f64 {
    if union == 0 {
        0.0
    } else {
        intersection as f64 / union as f64
    }
}

fn print_carrier_screen(
    index: &HashMap<Key, HashMap<String, String>>,
    all: &HashSet<&str>,
    samples: &[String],
) {
    println!(
        "Variants carried by ALL {} samples (het or hom-alt)",
        samples.len()
    );
    println!(
        "{:<8} {:<12} {:<6} {:<6} {}",
        "chrom", "pos", "ref", "alt", "genotypes"
    );
    println!("{}", "-".repeat(60));
    let mut n = 0usize;
    for ((c, p, r, a), per) in index {
        if all.iter().all(|s| per.contains_key(*s)) {
            n += 1;
            if n <= 200 {
                let gts: Vec<String> = samples
                    .iter()
                    .map(|s| format!("{}={}", s, per.get(s).cloned().unwrap_or_default()))
                    .collect();
                println!(
                    "{:<8} {:<12} {:<6} {:<6} {}",
                    c,
                    p,
                    truncate(r, 6),
                    truncate(a, 6),
                    gts.join(" ")
                );
            }
        }
    }
    if n > 200 {
        println!("… (showing first 200 of {n})");
    } else {
        println!("{n} sites");
    }
}

fn print_private(index: &HashMap<Key, HashMap<String, String>>, samples: &[String]) {
    let mut by_sample: HashMap<&str, Vec<&Key>> = HashMap::new();
    for (k, per) in index {
        if per.len() == 1 {
            if let Some(s) = per.keys().next() {
                by_sample.entry(s.as_str()).or_default().push(k);
            }
        }
    }
    for s in samples {
        let v = by_sample.get(s.as_str());
        let n = v.map(|v| v.len()).unwrap_or(0);
        println!("private to {s}: {n}");
        if let Some(v) = v {
            for k in v.iter().take(20) {
                println!("  {}\t{}\t{}\t{}", k.0, k.1, truncate(&k.2, 6), truncate(&k.3, 6));
            }
            if v.len() > 20 {
                println!("  … (+{} more)", v.len() - 20);
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

fn truncate(s: &str, n: usize) -> String {
    if s.len() > n {
        format!("{}…", &s[..n.saturating_sub(1)])
    } else {
        s.to_string()
    }
}
