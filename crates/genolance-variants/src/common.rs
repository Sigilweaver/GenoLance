//! Helpers shared by the clinical commands (`join`, `screen`,
//! `compound-het`, `pgx`, `compare`, `annotate`, `query`).

use std::collections::{HashMap, HashSet};

use arrow_array::{Array, Float32Array, RecordBatch, StringArray, UInt32Array, UInt64Array};
use futures::TryStreamExt;
use genolance_core::gene_lists::{gene_in_set, ACMG_SF_V3};
use genolance_core::schema::{
    CLINVAR_TABLE, META_CLINVAR_ASSEMBLY, META_CLINVAR_CONTIGS, META_CLINVAR_RELEASE,
    META_CLINVAR_SOURCE, SAMPLES_TABLE,
};
use genolance_core::store::Store;
use genolance_core::{Error, Result};
use lancedb::query::{ExecutableQuery, QueryBase, Select};

use crate::assembly::Assembly;
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
/// Bases are upper-cased first (VCF bases are case-insensitive, so a
/// soft-masked `acgt` must match ClinVar's `ACGT`).
///
/// The result is parsimonious (minimal), not left-aligned: without the
/// reference sequence a caller's right-aligned deletion in a repeat gets a
/// different key from ClinVar's left-aligned one. Normalise inputs with
/// `bcftools norm -f <reference.fa>` first.
pub(crate) fn normalize_alleles(pos: u64, r: &str, a: &str) -> (u64, String, String) {
    let (r, a) = (r.to_ascii_uppercase(), a.to_ascii_uppercase());
    if is_unnormalizable(&a) {
        return (pos, r, a);
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

/// Gene symbols as ClinVar writes them: letters, digits and `- . @ + '`.
/// Anything else cannot be a real symbol (and would be a LIKE wildcard).
fn plausible_symbol(g: &str) -> bool {
    !g.is_empty()
        && g.chars()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, '-' | '.' | '@' | '+' | '\''))
}

/// SQL predicate matching ClinVar rows whose `gene_symbols` list contains any
/// of the given symbols, ignoring case. Fails for a symbol that cannot occur in
/// ClinVar, naming it, instead of building a pattern from it.
pub(crate) fn gene_match_predicate(genes: &[String]) -> Result<String> {
    let mut ors = Vec::new();
    for g in genes {
        if !plausible_symbol(g) {
            return Err(Error::unknown_gene(vec![g.clone()]));
        }
        ors.push(format!(
            "lower(concat('|', array_to_string(gene_symbols, '|'), '|')) LIKE '%|{}|%'",
            sql_escape(&g.to_lowercase())
        ));
    }
    Ok(format!("({})", ors.join(" OR ")))
}

/// The symbols of a `List<Utf8>` cell.
pub(crate) fn list_strings(a: &arrow_array::ListArray, i: usize) -> Vec<String> {
    if a.is_null(i) {
        return Vec::new();
    }
    let values = a.value(i);
    match values.as_any().downcast_ref::<StringArray>() {
        Some(sa) => (0..sa.len()).map(|j| sa.value(j).to_string()).collect(),
        None => Vec::new(),
    }
}

pub(crate) fn req_list<'a>(
    b: &'a RecordBatch,
    table: &str,
    name: &str,
) -> Result<&'a arrow_array::ListArray> {
    b.column_by_name(name)
        .and_then(|c| c.as_any().downcast_ref::<arrow_array::ListArray>())
        .ok_or_else(|| missing_column(table, name))
}

/// Symbolic (`<DEL>`), breakend and `*` alleles cannot be trimmed and are
/// joined as written; ingest counts them.
pub(crate) fn is_unnormalizable(a: &str) -> bool {
    a.starts_with('<') || a.contains('[') || a.contains(']') || a == "*" || a == "."
}

/// The error for a BCF input: BCF is not supported by any command.
pub(crate) fn bcf_unsupported(path: &str) -> Error {
    Error::unrecognized_input(
        path,
        "BCF input is not supported; convert it to VCF first (for example `bcftools view -Oz`)",
    )
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
    /// Every gene symbol of the record, `,`-joined for display (None if it has none).
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
        let gene = req_list(b, t, "gene_symbols")?;
        let sig = req_str(b, t, "clinical_significance")?;
        let dis = req_str(b, t, "disease_name")?;
        for i in 0..b.num_rows() {
            let genes: Vec<String> = list_strings(gene, i);
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
                    Some(genes.join(","))
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

/// A registered sample as far as the reference checks need it.
pub(crate) struct Registered {
    pub name: String,
    pub contigs: HashMap<String, u64>,
    pub assembly: Option<Assembly>,
}

pub(crate) async fn load_registry(store: &Store) -> Result<Vec<Registered>> {
    let reg = store.conn.open_table(SAMPLES_TABLE).execute().await?;
    let batches: Vec<RecordBatch> = reg.query().execute().await?.try_collect().await?;
    let mut out = Vec::new();
    for b in &batches {
        let names = req_str(b, SAMPLES_TABLE, "sample_name")?;
        let contigs = req_str(b, SAMPLES_TABLE, "contigs")?;
        let assembly = req_str(b, SAMPLES_TABLE, "assembly")?;
        for i in 0..b.num_rows() {
            out.push(Registered {
                name: names.value(i).to_string(),
                contigs: parse_contigs(&opt_string(contigs, i).unwrap_or_default()),
                assembly: opt_string(assembly, i).and_then(|l| Assembly::from_label(&l)),
            });
        }
    }
    Ok(out)
}

fn label(a: Option<Assembly>) -> String {
    a.map(|a| a.to_string()).unwrap_or_else(|| "unknown".into())
}

/// Print the stored ClinVar release and check that every sample in scope was
/// made against the same assembly as ClinVar.
///
/// The assembly of each side was inferred at ingest from its contig lengths
/// and `##reference` line (see the `assembly` module) or declared with
/// `--assembly`. Equal labels pass, different labels fail. If either is
/// unknown, the lengths of contigs both headers declare are compared; if there
/// are none, the check fails and says how to declare the assembly. With
/// `f.skip_reference_check` the same check runs and its outcome is printed as
/// a warning, but never stops the command. `samples` of `None` means every
/// registered sample. A sample that is not registered is an error.
pub(crate) async fn check_reference(
    store: &Store,
    samples: Option<&[String]>,
    f: &Filters,
    cmd: &str,
) -> Result<()> {
    let release = store
        .meta(META_CLINVAR_RELEASE)
        .await?
        .unwrap_or_else(|| "unknown (the file has no ##fileDate)".into());
    let source = store.meta(META_CLINVAR_SOURCE).await?.unwrap_or_default();
    let cv_assembly = store
        .meta(META_CLINVAR_ASSEMBLY)
        .await?
        .and_then(|l| Assembly::from_label(&l));
    eprintln!(
        "[{cmd}] ClinVar release: {release} ({source}), assembly {}",
        label(cv_assembly)
    );

    let cv_contigs = parse_contigs(&store.meta(META_CLINVAR_CONTIGS).await?.unwrap_or_default());
    let registered = load_registry(store).await?;
    if let Some(wanted) = samples {
        for w in wanted {
            if !registered.iter().any(|r| &r.name == w) {
                return Err(Error::unknown_sample(
                    w.clone(),
                    registered.iter().map(|r| r.name.clone()).collect(),
                ));
            }
        }
    }
    let mut failure: Option<(String, String)> = None;
    for r in &registered {
        if samples.is_some_and(|w| !w.contains(&r.name)) {
            continue;
        }
        if let Err(detail) = compare_with_clinvar(r, cv_assembly, &cv_contigs, &release) {
            failure = Some((r.name.clone(), detail));
            break;
        }
    }
    match (failure, f.skip_reference_check) {
        (None, false) => Ok(()),
        (None, true) => {
            eprintln!(
                "[{cmd}] WARNING: --skip-reference-check is set; the reference check was not enforced (it would have passed)"
            );
            Ok(())
        }
        (Some((sample, detail)), true) => {
            eprintln!(
                "[{cmd}] WARNING: --skip-reference-check is set; the reference check was not enforced. For sample '{sample}' it would have failed: {detail}"
            );
            Ok(())
        }
        (Some((sample, detail)), false) => Err(Error::reference_mismatch(sample, detail)),
    }
}

fn compare_with_clinvar(
    r: &Registered,
    cv_assembly: Option<Assembly>,
    cv_contigs: &HashMap<String, u64>,
    release: &str,
) -> std::result::Result<(), String> {
    if let (Some(sa), Some(ca)) = (r.assembly, cv_assembly) {
        return if sa == ca {
            Ok(())
        } else {
            Err(format!(
                "the sample is {sa} but the ClinVar release {release} is {ca}; use a ClinVar release for {sa}"
            ))
        };
    }
    let mut shared = 0usize;
    for (chrom, len) in &r.contigs {
        if let Some(cv_len) = cv_contigs.get(chrom) {
            shared += 1;
            if cv_len != len {
                return Err(format!(
                    "contig {chrom} is {len} bases in the sample's header but {cv_len} in the ClinVar header; they are different assemblies"
                ));
            }
        }
    }
    if shared > 0 {
        return Ok(());
    }
    Err(format!(
        "cannot tell whether the sample and ClinVar use the same assembly (sample: {}, ClinVar: {}). Declare the assembly of the one that is unknown by re-ingesting it with `--assembly GRCh37` or `--assembly GRCh38` (and `--replace`), or pass --skip-reference-check",
        label(r.assembly),
        label(cv_assembly)
    ))
}

/// Check that the listed samples were made against one assembly. Known and
/// different is an error; unknown is a printed warning (the commands using
/// this never read ClinVar, so there is nothing to declare against).
pub(crate) async fn check_samples_agree(
    store: &Store,
    samples: &[String],
    cmd: &str,
) -> Result<()> {
    let registered = load_registry(store).await?;
    let picked: Vec<&Registered> = registered
        .iter()
        .filter(|r| samples.contains(&r.name))
        .collect();
    let mut first: Option<(&str, Assembly)> = None;
    for r in &picked {
        match (r.assembly, first) {
            (Some(a), None) => first = Some((&r.name, a)),
            (Some(a), Some((n0, a0))) if a != a0 => {
                return Err(Error::reference_mismatch(
                    r.name.clone(),
                    format!("sample '{}' is {a} but sample '{n0}' is {a0}; their positions are not comparable", r.name),
                ));
            }
            (None, _) => eprintln!(
                "[{cmd}] WARNING: the assembly of sample '{}' is unknown; its positions are assumed to match the other samples'",
                r.name
            ),
            _ => {}
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
        let genes = {
            use arrow_array::builder::{ListBuilder, StringBuilder};
            let mut b = ListBuilder::new(StringBuilder::new());
            for r in rows {
                if r.4.is_empty() {
                    b.append(false);
                } else {
                    b.values().append_value(r.4);
                    b.append(true);
                }
            }
            b.finish()
        };
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
                Arc::new(genes),
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
