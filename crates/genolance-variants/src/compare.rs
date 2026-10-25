use std::collections::{HashMap, HashSet};

use arrow_array::RecordBatch;
use futures::TryStreamExt;
use genolance_core::{Error, Result};
use lancedb::query::{ExecutableQuery, QueryBase};

use genolance_core::schema::VARIANTS_TABLE;
use genolance_core::store::Store;

use crate::common::{carrier_state, opt_string, req_str, req_u64, require_samples, Carrier};

type Key = (String, u64, String, String);

/// Compare variants across two or more samples.
pub async fn run(store_path: &str, samples: &[String], mode: &str) -> Result<()> {
    if samples.len() < 2 {
        return Err(invalid!("compare needs at least 2 samples"));
    }
    if !matches!(mode, "concordance" | "carrier-screen" | "private") {
        return Err(invalid!(
            "unknown compare mode: {mode} (expected concordance, carrier-screen or private)"
        ));
    }

    let store = Store::open(store_path).await?;
    require_samples(&store, samples).await?;
    let tables = store.variants.table_names().execute().await?;
    if !tables.iter().any(|n| n == VARIANTS_TABLE) {
        return Err(Error::missing_table(VARIANTS_TABLE, "ingest samples first"));
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

    // Index: variant key -> map of sample_name -> genotype. Only carriers
    // are indexed; no-calls are unknown, not variants, and are counted.
    let mut index: HashMap<Key, HashMap<String, String>> = HashMap::new();
    let mut no_calls: HashMap<String, usize> = HashMap::new();
    for b in &batches {
        let t = VARIANTS_TABLE;
        let sample = req_str(b, t, "sample_name")?;
        let chrom = req_str(b, t, "chrom")?;
        let pos = req_u64(b, t, "pos")?;
        let ref_a = req_str(b, t, "ref_allele")?;
        let alt_a = req_str(b, t, "alt_allele")?;
        let gt = req_str(b, t, "genotype")?;
        for i in 0..b.num_rows() {
            let s = sample.value(i).to_string();
            let g = opt_string(gt, i);
            match carrier_state(g.as_deref()) {
                Carrier::Yes => {}
                Carrier::No => continue,
                Carrier::NoCall => {
                    *no_calls.entry(s).or_insert(0) += 1;
                    continue;
                }
            }
            let key: Key = (
                chrom.value(i).to_string(),
                pos.value(i),
                ref_a.value(i).to_string(),
                alt_a.value(i).to_string(),
            );
            index
                .entry(key)
                .or_default()
                .insert(s, g.unwrap_or_default());
        }
    }

    let all_samples: HashSet<&str> = samples.iter().map(String::as_str).collect();

    match mode {
        "concordance" => print_concordance(&index, samples),
        "carrier-screen" => print_carrier_screen(&index, &all_samples, samples),
        "private" => print_private(&index, samples),
        other => return Err(invalid!("unknown compare mode: {other}")),
    }
    for s in samples {
        if let Some(n) = no_calls.get(s) {
            eprintln!("[compare] {n} no-call genotype row(s) for {s} were not counted as variants");
        }
    }
    Ok(())
}

fn print_concordance(index: &HashMap<Key, HashMap<String, String>>, samples: &[String]) {
    let mut shared_all = 0usize;
    let mut private_counts: HashMap<&str, usize> = HashMap::new();
    let mut shared_pair: HashMap<(&str, &str), usize> = HashMap::new();
    let total = index.len();

    for per_sample in index.values() {
        if per_sample.len() == samples.len() {
            shared_all += 1;
        } else if per_sample.len() == 1 {
            if let Some(s) = per_sample.keys().next() {
                *private_counts.entry(s.as_str()).or_insert(0) += 1;
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
        "{:<8} {:<12} {:<6} {:<6} genotypes",
        "chrom", "pos", "ref", "alt"
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
        println!("... (showing first 200 of {n})");
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
                println!(
                    "  {}\t{}\t{}\t{}",
                    k.0,
                    k.1,
                    truncate(&k.2, 6),
                    truncate(&k.3, 6)
                );
            }
            if v.len() > 20 {
                println!("  ... (+{} more)", v.len() - 20);
            }
        }
    }
}

fn truncate(s: &str, n: usize) -> String {
    if s.len() > n {
        format!("{}...", &s[..n.saturating_sub(1)])
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn jaccard_of_empty_union_is_zero() {
        assert_eq!(jaccard(0, 0), 0.0);
    }

    #[test]
    fn jaccard_of_full_overlap_is_one() {
        assert_eq!(jaccard(10, 10), 1.0);
    }

    #[test]
    fn jaccard_is_intersection_over_union() {
        assert!((jaccard(1, 4) - 0.25).abs() < 1e-12);
    }

    #[test]
    fn truncate_short_string_is_unchanged() {
        assert_eq!(truncate("abc", 10), "abc");
    }

    #[test]
    fn truncate_long_string_gets_ellipsis() {
        assert_eq!(truncate("abcdefgh", 4), "abc...");
    }

    #[test]
    fn truncate_at_exact_length_is_unchanged() {
        assert_eq!(truncate("abcd", 4), "abcd");
    }
}
