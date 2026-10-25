//! Helpers shared by the clinical commands (`join`, `screen`,
//! `compound-het`, `pgx`, `compare`, `annotate`, `query`).

use std::collections::{HashMap, HashSet};

use arrow_array::{Array, Float32Array, RecordBatch, StringArray, UInt32Array, UInt64Array};
use futures::TryStreamExt;
use genolance_core::gene_lists::{gene_in_set, ACMG_SF_V3};
use genolance_core::schema::{
    CLINVAR_TABLE, META_CLINVAR_CONTIGS, META_CLINVAR_RELEASE, META_CLINVAR_SOURCE, SAMPLES_TABLE,
};
use genolance_core::store::Store;
use genolance_core::{Error, Result};
use lancedb::query::{ExecutableQuery, QueryBase, Select};

use crate::filters::Filters;

pub(crate) fn sql_escape(s: &str) -> String {
    s.replace('\'', "''")
}

pub(crate) fn truncate(s: &str, n: usize) -> String {
    if s.chars().count() > n {
        let cut: String = s.chars().take(n.saturating_sub(1)).collect();
        format!("{cut}...")
    } else {
        s.to_string()
    }
}

/// Canonical chromosome key: `chr` prefix stripped, mitochondrial `M`
/// folded into `MT`, so `chrM`, `M` and `MT` compare equal.
pub(crate) fn normalize_chrom(c: &str) -> String {
    let bare = c.strip_prefix("chr").unwrap_or(c);
    if bare == "M" {
        "MT".to_string()
    } else {
        bare.to_string()
    }
}

/// Every spelling of a chromosome that a store might use.
pub(crate) fn chrom_forms(c: &str) -> Vec<String> {
    let norm = normalize_chrom(c);
    let bases: Vec<&str> = if norm == "MT" {
        vec!["MT", "M"]
    } else {
        vec![norm.as_str()]
    };
    let mut out = Vec::new();
    for b in bases {
        out.push(b.to_string());
        out.push(format!("chr{b}"));
    }
    out
}

/// SQL predicate matching a chromosome in any of its spellings.
pub(crate) fn chrom_predicate(c: &str) -> String {
    let ors: Vec<String> = chrom_forms(c)
        .iter()
        .map(|f| format!("chrom = '{}'", sql_escape(f)))
        .collect();
    format!("({})", ors.join(" OR "))
}

/// Whether a sample carries the ALT allele of a split row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Carrier {
    /// At least one called allele is the ALT (`1` in the per-ALT genotype).
    Yes,
    /// Called, and no allele is the ALT.
    No,
    /// Missing, empty or all-`.` genotype: unknown, never counted as a carrier.
    NoCall,
}

pub(crate) fn carrier_state(gt: Option<&str>) -> Carrier {
    let Some(gt) = gt else {
        return Carrier::NoCall;
    };
    let mut called = false;
    for tok in gt.split(['/', '|']) {
        match tok {
            "1" => return Carrier::Yes,
            "" | "." => {}
            _ => called = true,
        }
    }
    if called {
        Carrier::No
    } else {
        Carrier::NoCall
    }
}

/// Fail with `Error::UnknownSample` unless every name is in the `samples`
/// registry.
pub(crate) async fn require_samples(store: &Store, names: &[String]) -> Result<()> {
    let root_tables = store.conn.table_names().execute().await?;
    let mut available: Vec<String> = Vec::new();
    if root_tables.iter().any(|n| n == SAMPLES_TABLE) {
        let samples = store.conn.open_table(SAMPLES_TABLE).execute().await?;
        let batches: Vec<RecordBatch> = samples.query().execute().await?.try_collect().await?;
        for b in &batches {
            let col = req_str(b, SAMPLES_TABLE, "sample_name")?;
            for i in 0..col.len() {
                available.push(col.value(i).to_string());
            }
        }
    }
    for n in names {
        if !available.iter().any(|a| a == n) {
            return Err(Error::unknown_sample(n.clone(), available));
        }
    }
    Ok(())
}

/// Trim a REF/ALT pair the way `bcftools norm` does: drop trailing bases the
/// two share, then leading bases (moving POS), always keeping one base.
/// Symbolic (`<DEL>`), breakend and `*` alleles are returned unchanged.
pub(crate) fn normalize_alleles(pos: u64, r: &str, a: &str) -> (u64, String, String) {
    if a.starts_with('<') || a.contains('[') || a.contains(']') || a == "*" || a == "." {
        return (pos, r.to_string(), a.to_string());
    }
    let (mut rb, mut ab) = (r.as_bytes(), a.as_bytes());
    while rb.len() > 1 && ab.len() > 1 && rb.last() == ab.last() {
        rb = &rb[..rb.len() - 1];
        ab = &ab[..ab.len() - 1];
    }
    let mut p = pos;
    while rb.len() > 1 && ab.len() > 1 && rb[0] == ab[0] {
        rb = &rb[1..];
        ab = &ab[1..];
        p += 1;
    }
    (
        p,
        String::from_utf8_lossy(rb).into_owned(),
        String::from_utf8_lossy(ab).into_owned(),
    )
}

/// SQL predicate matching ClinVar rows whose `gene_symbols` contain any of
/// the given symbols as a whole symbol.
pub(crate) fn gene_match_predicate(genes: &[String]) -> String {
    let ors: Vec<String> = genes
        .iter()
        .map(|g| format!("gene_symbols LIKE '%|{}|%'", sql_escape(g)))
        .collect();
    format!("({})", ors.join(" OR "))
}

pub(crate) fn fmt_qual(q: Option<f32>) -> String {
    q.map(|v| format!("{v:.1}")).unwrap_or_else(|| ".".into())
}

pub(crate) fn fmt_dp(d: Option<u32>) -> String {
    d.map(|v| v.to_string()).unwrap_or_else(|| ".".into())
}

pub(crate) fn missing_column(table: &str, column: &str) -> Error {
    Error::missing_column(table, column)
}

/// Required columns: absent or retyped columns are errors, not empty output.
pub(crate) fn req_str<'a>(b: &'a RecordBatch, table: &str, name: &str) -> Result<&'a StringArray> {
    b.column_by_name(name)
        .and_then(|c| c.as_any().downcast_ref::<StringArray>())
        .ok_or_else(|| missing_column(table, name))
}
pub(crate) fn req_u64<'a>(b: &'a RecordBatch, table: &str, name: &str) -> Result<&'a UInt64Array> {
    b.column_by_name(name)
        .and_then(|c| c.as_any().downcast_ref::<UInt64Array>())
        .ok_or_else(|| missing_column(table, name))
}
pub(crate) fn req_u32<'a>(b: &'a RecordBatch, table: &str, name: &str) -> Result<&'a UInt32Array> {
    b.column_by_name(name)
        .and_then(|c| c.as_any().downcast_ref::<UInt32Array>())
        .ok_or_else(|| missing_column(table, name))
}
pub(crate) fn req_f32<'a>(b: &'a RecordBatch, table: &str, name: &str) -> Result<&'a Float32Array> {
    b.column_by_name(name)
        .and_then(|c| c.as_any().downcast_ref::<Float32Array>())
        .ok_or_else(|| missing_column(table, name))
}

pub(crate) fn opt_string(a: &StringArray, i: usize) -> Option<String> {
    use arrow_array::Array;
    (!a.is_null(i)).then(|| a.value(i).to_string())
}

/// Fail if the ClinVar table holds two rows for one (chrom, pos, ref, alt),
/// which happens when two releases were ingested into it. Looking one up
/// would then return whichever row the scan produced last.
pub(crate) async fn ensure_clinvar_unique(table: &lancedb::Table) -> Result<()> {
    let batches: Vec<RecordBatch> = table
        .query()
        .select(Select::columns(&[
            "chrom", "norm_pos", "norm_ref", "norm_alt",
        ]))
        .execute()
        .await?
        .try_collect()
        .await?;
    let mut seen: HashSet<(String, u64, String, String)> = HashSet::new();
    for b in &batches {
        let t = genolance_core::schema::CLINVAR_TABLE;
        let chrom = req_str(b, t, "chrom")?;
        let npos = req_u64(b, t, "norm_pos")?;
        let nref = req_str(b, t, "norm_ref")?;
        let nalt = req_str(b, t, "norm_alt")?;
        for i in 0..b.num_rows() {
            let key = (
                normalize_chrom(chrom.value(i)),
                npos.value(i),
                nref.value(i).to_string(),
                nalt.value(i).to_string(),
            );
            if !seen.insert(key.clone()) {
                return Err(Error::clinvar_conflict(
                    key.0,
                    key.1,
                    key.2,
                    key.3,
                    "the clinvar table has more than one row for this variant (it holds more than one release); re-ingest ClinVar with --replace",
                ));
            }
        }
    }
    Ok(())
}

/// ClinVar significance predicate: exact list, else substring, else
/// `default_substring`, else match everything. An empty exact list is an error.
pub(crate) fn significance_predicate(
    f: &Filters,
    default_substring: Option<&str>,
) -> Result<String> {
    if let Some(list) = &f.significance_exact {
        if list.is_empty() {
            return Err(Error::InvalidArgument(
                "--significance-exact has no values; refusing to match everything".into(),
            ));
        }
        let ors: Vec<String> = list
            .iter()
            .map(|s| {
                format!(
                    "lower(clinical_significance) = '{}'",
                    sql_escape(&s.to_lowercase())
                )
            })
            .collect();
        return Ok(format!("({})", ors.join(" OR ")));
    }
    match f.significance_substring.as_deref().or(default_substring) {
        Some(sub) => Ok(format!(
            "lower(clinical_significance) LIKE '%{}%'",
            sql_escape(&sub.to_lowercase())
        )),
        None => Ok("1=1".to_string()),
    }
}

pub(crate) type Key = (String, u64, String, String);

#[derive(Debug)]
pub(crate) struct Annotation {
    /// Every gene symbol of the record, `|`-joined (None if it has none).
    pub gene_symbol: Option<String>,
    pub genes: Vec<String>,
    pub clinical_significance: Option<String>,
    pub disease_name: Option<String>,
}

pub(crate) fn build_clinvar_map(
    batches: &[RecordBatch],
    f: &Filters,
) -> Result<HashMap<Key, Annotation>> {
    let mut map: HashMap<Key, Annotation> = HashMap::new();
    for b in batches {
        let t = CLINVAR_TABLE;
        let chrom = req_str(b, t, "chrom")?;
        let npos = req_u64(b, t, "norm_pos")?;
        let nref = req_str(b, t, "norm_ref")?;
        let nalt = req_str(b, t, "norm_alt")?;
        let gene = req_str(b, t, "gene_symbols")?;
        let sig = req_str(b, t, "clinical_significance")?;
        let dis = req_str(b, t, "disease_name")?;
        for i in 0..b.num_rows() {
            let genes: Vec<String> = opt_string(gene, i)
                .unwrap_or_default()
                .split('|')
                .filter(|g| !g.is_empty())
                .map(str::to_string)
                .collect();
            if f.acmg_only && !genes.iter().any(|g| gene_in_set(g, ACMG_SF_V3)) {
                continue;
            }
            if let Some(list) = &f.gene_filter {
                if !genes
                    .iter()
                    .any(|g| list.iter().any(|w| w.eq_ignore_ascii_case(g)))
                {
                    continue;
                }
            }
            let key: Key = (
                normalize_chrom(chrom.value(i)),
                npos.value(i),
                nref.value(i).to_string(),
                nalt.value(i).to_string(),
            );
            let ann = Annotation {
                gene_symbol: if genes.is_empty() {
                    None
                } else {
                    Some(genes.join("|"))
                },
                genes,
                clinical_significance: opt_string(sig, i),
                disease_name: opt_string(dis, i),
            };
            if let Some(prev) = map.get(&key) {
                if prev.clinical_significance != ann.clinical_significance {
                    let detail = format!(
                        "ClinVar has two rows with different significance ({:?} and {:?}); the store holds more than one release. Re-ingest with --replace",
                        prev.clinical_significance, ann.clinical_significance
                    );
                    return Err(Error::clinvar_conflict(key.0, key.1, key.2, key.3, detail));
                }
                continue;
            }
            map.insert(key, ann);
        }
    }
    Ok(map)
}

/// Parse "name:length,name:length" (names may contain ':'; the length is
/// after the last one) into normalized-chromosome -> length.
fn parse_contigs(s: &str) -> HashMap<String, u64> {
    s.split(',')
        .filter_map(|e| {
            let (name, len) = e.rsplit_once(':')?;
            Some((normalize_chrom(name), len.parse::<u64>().ok()?))
        })
        .collect()
}

/// Print the stored ClinVar release and check that every sample in scope
/// was made against the same assembly as ClinVar, by comparing the lengths
/// of the contigs both headers declare. A mismatch, or no contig to compare,
/// is an error unless `f.skip_reference_check` is set. `samples` of `None`
/// means every registered sample. Also fails for a sample that is not
/// registered.
pub(crate) async fn check_reference(
    store: &Store,
    samples: Option<&[String]>,
    f: &Filters,
    cmd: &str,
) -> Result<()> {
    let release = store
        .meta(META_CLINVAR_RELEASE)
        .await?
        .unwrap_or_else(|| "unknown".into());
    let source = store.meta(META_CLINVAR_SOURCE).await?.unwrap_or_default();
    eprintln!("[{cmd}] ClinVar release: {release} ({source})");

    let cv_contigs = parse_contigs(&store.meta(META_CLINVAR_CONTIGS).await?.unwrap_or_default());
    let reg = store.conn.open_table(SAMPLES_TABLE).execute().await?;
    let batches: Vec<RecordBatch> = reg.query().execute().await?.try_collect().await?;
    let mut registered: Vec<(String, HashMap<String, u64>)> = Vec::new();
    for b in &batches {
        let names = req_str(b, SAMPLES_TABLE, "sample_name")?;
        let contigs = req_str(b, SAMPLES_TABLE, "contigs")?;
        for i in 0..b.num_rows() {
            registered.push((
                names.value(i).to_string(),
                parse_contigs(&opt_string(contigs, i).unwrap_or_default()),
            ));
        }
    }
    if let Some(wanted) = samples {
        for w in wanted {
            if !registered.iter().any(|(n, _)| n == w) {
                return Err(Error::unknown_sample(
                    w.clone(),
                    registered.iter().map(|(n, _)| n.clone()).collect(),
                ));
            }
        }
    }
    if f.skip_reference_check {
        return Ok(());
    }
    for (name, contigs) in &registered {
        if samples.is_some_and(|w| !w.contains(name)) {
            continue;
        }
        let mut shared = 0usize;
        for (chrom, len) in contigs {
            if let Some(cv_len) = cv_contigs.get(chrom) {
                shared += 1;
                if cv_len != len {
                    return Err(Error::reference_mismatch(name.clone(), format!(
                            "contig {chrom} is {len} bases in the sample's header but {cv_len} in the ClinVar header; they are different assemblies (ClinVar release {release})"
                        )));
                }
            }
        }
        if shared == 0 {
            return Err(Error::reference_mismatch(name.clone(), "no contig with a declared length is shared with the ClinVar header, so the assembly cannot be verified; pass --skip-reference-check to proceed anyway"));
        }
    }
    Ok(())
}

/// Check the list-valued filters: an empty list would silently match
/// everything (exact significance) or nothing (genes).
pub(crate) fn validate_filters(f: &Filters) -> Result<()> {
    if let Some(list) = &f.significance_exact {
        if list.is_empty() {
            return Err(Error::InvalidArgument(
                "--significance-exact has no values after splitting; refusing to treat it as match-all".into(),
            ));
        }
    }
    if f.significance_substring
        .as_deref()
        .is_some_and(|s| s.trim().is_empty())
    {
        return Err(Error::InvalidArgument(
            "--significance is empty; an empty substring would match every record".into(),
        ));
    }
    if let Some(list) = &f.gene_filter {
        if list.is_empty() {
            return Err(Error::InvalidArgument(
                "--genes has no values after splitting; refusing to treat it as match-none or match-all".into(),
            ));
        }
    }
    Ok(())
}

/// Row-level admission for hits: FILTER, and QUAL/DP thresholds. Rows it
/// turns away for a reason other than a threshold being missed are counted
/// and reported, never dropped silently.
pub(crate) struct Gate<'a> {
    f: &'a Filters,
    pub excluded_filtered: usize,
    pub excluded_missing: usize,
    pub no_call: usize,
}

impl<'a> Gate<'a> {
    pub(crate) fn new(f: &'a Filters) -> Self {
        Self {
            f,
            excluded_filtered: 0,
            excluded_missing: 0,
            no_call: 0,
        }
    }

    /// `filter` is the stored FILTER text (None = "."). Returns whether the
    /// call may be reported.
    pub(crate) fn admit(
        &mut self,
        filter: Option<&str>,
        qual: Option<f32>,
        dp: Option<u32>,
    ) -> bool {
        if !self.f.include_filtered {
            if let Some(v) = filter {
                if v != "PASS" {
                    self.excluded_filtered += 1;
                    return false;
                }
            }
        }
        if let Some(min) = self.f.min_qual {
            match qual {
                Some(q) if q >= min => {}
                Some(_) => return false,
                None if self.f.keep_missing => {}
                None => {
                    self.excluded_missing += 1;
                    return false;
                }
            }
        }
        if let Some(min) = self.f.min_dp {
            match dp {
                Some(d) if d >= min => {}
                Some(_) => return false,
                None if self.f.keep_missing => {}
                None => {
                    self.excluded_missing += 1;
                    return false;
                }
            }
        }
        true
    }

    pub(crate) fn report(&self, cmd: &str) {
        if self.excluded_filtered > 0 {
            eprintln!(
                "[{cmd}] {} matching call(s) excluded because FILTER is not PASS (use --include-filtered to include them)",
                self.excluded_filtered
            );
        }
        if self.excluded_missing > 0 {
            eprintln!(
                "[{cmd}] {} matching call(s) excluded because QUAL or DP is missing and a minimum was set (use --keep-missing to include them)",
                self.excluded_missing
            );
        }
        if self.no_call > 0 {
            eprintln!(
                "[{cmd}] {} matching call(s) skipped because the genotype is a no-call",
                self.no_call
            );
        }
    }
}

#[cfg(test)]
#[allow(clippy::field_reassign_with_default)]
mod tests {
    use super::*;

    #[test]
    fn carrier_states() {
        assert_eq!(carrier_state(Some("0/1")), Carrier::Yes);
        assert_eq!(carrier_state(Some("1")), Carrier::Yes);
        assert_eq!(carrier_state(Some("1|1")), Carrier::Yes);
        assert_eq!(carrier_state(Some("0/0")), Carrier::No);
        assert_eq!(carrier_state(Some("0")), Carrier::No);
        assert_eq!(carrier_state(Some("0/.")), Carrier::No);
        assert_eq!(carrier_state(Some("./.")), Carrier::NoCall);
        assert_eq!(carrier_state(Some("")), Carrier::NoCall);
        assert_eq!(carrier_state(None), Carrier::NoCall);
    }

    #[test]
    fn allele_trimming() {
        assert_eq!(
            normalize_alleles(100, "GCA", "GA"),
            (100, "GC".into(), "G".into())
        );
        assert_eq!(
            normalize_alleles(100, "GCA", "G"),
            (100, "GCA".into(), "G".into())
        );
        assert_eq!(
            normalize_alleles(100, "AT", "AC"),
            (101, "T".into(), "C".into())
        );
        assert_eq!(
            normalize_alleles(100, "A", "G"),
            (100, "A".into(), "G".into())
        );
        assert_eq!(
            normalize_alleles(100, "A", "<DEL>"),
            (100, "A".into(), "<DEL>".into())
        );
    }

    #[test]
    fn chromosome_aliases() {
        assert_eq!(normalize_chrom("chrM"), "MT");
        assert_eq!(normalize_chrom("MT"), "MT");
        assert_eq!(normalize_chrom("chr17"), "17");
        assert_eq!(chrom_forms("chrM"), vec!["MT", "chrMT", "M", "chrM"]);
        assert_eq!(chrom_forms("17"), vec!["17", "chr17"]);
    }

    #[test]
    fn gate_counts_exclusions() {
        let mut f = Filters::default();
        f.min_qual = Some(30.0);
        let mut g = Gate::new(&f);
        assert!(g.admit(Some("PASS"), Some(40.0), None));
        assert!(g.admit(None, Some(40.0), None));
        assert!(!g.admit(Some("LowQual"), Some(40.0), None));
        assert!(!g.admit(Some("PASS"), None, None));
        assert!(!g.admit(Some("PASS"), Some(10.0), None));
        assert_eq!((g.excluded_filtered, g.excluded_missing), (1, 1));
        f.include_filtered = true;
        f.keep_missing = true;
        let mut g = Gate::new(&f);
        assert!(g.admit(Some("LowQual"), None, None));
    }
}

#[cfg(test)]
pub(crate) mod testutil {
    use std::sync::Arc;

    use arrow_array::{RecordBatch, StringArray, UInt64Array};
    use genolance_core::schema::clinvar_schema;

    type Row<'a> = (&'a str, u64, &'a str, &'a str, &'a str, &'a str, &'a str);

    /// ClinVar-shaped batch from `(chrom, pos, ref, alt, gene, significance,
    /// disease)` rows; `norm_*` equal the raw alleles, `gene` becomes
    /// `gene_symbols` (`|GENE|`, or null when empty).
    pub(crate) fn clinvar_batch(rows: &[Row]) -> RecordBatch {
        let col = |f: fn(&Row) -> String| {
            Arc::new(StringArray::from(rows.iter().map(f).collect::<Vec<_>>()))
        };
        let genes: Vec<Option<String>> = rows
            .iter()
            .map(|r| (!r.4.is_empty()).then(|| format!("|{}|", r.4)))
            .collect();
        let first: Vec<Option<String>> = rows
            .iter()
            .map(|r| (!r.4.is_empty()).then(|| r.4.to_string()))
            .collect();
        RecordBatch::try_new(
            clinvar_schema(),
            vec![
                col(|r| r.0.to_string()),
                Arc::new(UInt64Array::from(
                    rows.iter().map(|r| r.1).collect::<Vec<_>>(),
                )),
                col(|r| r.2.to_string()),
                col(|r| r.3.to_string()),
                Arc::new(StringArray::from(vec![None::<String>; rows.len()])),
                Arc::new(StringArray::from(first)),
                col(|r| r.5.to_string()),
                Arc::new(StringArray::from(vec![None::<String>; rows.len()])),
                col(|r| r.6.to_string()),
                Arc::new(StringArray::from(genes)),
                Arc::new(UInt64Array::from(
                    rows.iter().map(|r| r.1).collect::<Vec<_>>(),
                )),
                col(|r| r.2.to_string()),
                col(|r| r.3.to_string()),
            ],
        )
        .unwrap()
    }
}
