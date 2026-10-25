//! Filters shared by `join`, `screen` and `compound-het`.

/// Filters on ClinVar annotations and on the sample's calls.
///
/// Construct with `Filters::default()` and assign the fields you need.
#[derive(Default, Clone, Debug)]
#[non_exhaustive]
pub struct Filters {
    /// Case-insensitive substring of `clinical_significance`.
    pub significance_substring: Option<String>,
    /// Exact (case-insensitive) `clinical_significance` values. An empty
    /// list is an error.
    pub significance_exact: Option<Vec<String>>,
    /// Drop calls with QUAL below this value.
    pub min_qual: Option<f32>,
    /// Drop calls with FORMAT/DP below this value.
    pub min_dp: Option<u32>,
    /// Restrict to the ACMG secondary-findings gene list.
    pub acmg_only: bool,
    /// Restrict to these genes. An empty list is an error.
    pub gene_filter: Option<Vec<String>>,
    /// Include calls whose FILTER is not PASS. By default they are excluded
    /// and counted.
    pub include_filtered: bool,
    /// Keep calls whose QUAL/DP is missing when `min_qual`/`min_dp` is set.
    /// By default they are excluded and counted.
    pub keep_missing: bool,
}
