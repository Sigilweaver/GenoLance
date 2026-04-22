//! Curated gene lists used for clinical filtering.
//!
//! These are short enough to keep inline rather than ship as data files.
//! Add new lists here as they become relevant.

/// ACMG Secondary Findings v3.2 (Miller et al., 2023).
///
/// These are the genes the ACMG recommends laboratories report on an
/// opt-out basis when performing clinical exome/genome sequencing,
/// regardless of the original test indication, because they harbour
/// variants with well-established clinical actionability.
///
/// Current as of v3.2 (published 2023). This list is updated every
/// ~2 years — verify against the current recommendation before using
/// in a clinical context.
pub const ACMG_SF_V3: &[&str] = &[
    // Cancer (hereditary tumor syndromes)
    "APC", "BMPR1A", "BRCA1", "BRCA2", "BRIP1", "CDH1", "CDK4", "CDKN2A", "DICER1", "EPCAM",
    "HOXB13", "MAX", "MEN1", "MLH1", "MSH2", "MSH6", "MUTYH", "NF2", "PALB2", "PMS2", "PTEN",
    "RAD51C", "RAD51D", "RB1", "RET", "SDHAF2", "SDHB", "SDHC", "SDHD", "SMAD4", "STK11", "TMEM127",
    "TP53", "TSC1", "TSC2", "VHL", "WT1",
    // Cardiovascular (cardiomyopathy, arrhythmia, aortopathy, FH, etc.)
    "ACTA2", "ACTC1", "APOB", "CASQ2", "COL3A1", "DES", "DSC2", "DSG2", "DSP", "FBN1", "FLNC",
    "GLA", "HNF1A", "KCNH2", "KCNQ1", "LDLR", "LMNA", "MYBPC3", "MYH11", "MYH7", "MYL2", "MYL3",
    "PCSK9", "PKP2", "PRKAG2", "RYR2", "SCN5A", "SMAD3", "TGFBR1", "TGFBR2", "TMEM43", "TNNC1",
    "TNNI3", "TNNT2", "TPM1", "TRDN", "TTN", "TTR",
    // Inborn errors / misc (malignant hyperthermia, Ehlers-Danlos, Wilson)
    "ATP7B", "BAG3", "CALM1", "CALM2", "CALM3", "GAA", "OTC", "PALLD", "RPE65", "RYR1", "TTR",
];

/// Split a comma-separated significance list on commas and semicolons, trim
/// whitespace, lowercase, de-duplicate.
pub fn parse_significance_list(raw: &str) -> Vec<String> {
    raw.split([',', ';'])
        .map(|s| s.trim().to_lowercase())
        .filter(|s| !s.is_empty())
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect()
}

/// Lowercase, whitespace-trim helper for a comma-separated gene list.
pub fn parse_gene_list(raw: &str) -> Vec<String> {
    raw.split([',', ';', '\n'])
        .map(|s| s.trim().to_uppercase())
        .filter(|s| !s.is_empty())
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect()
}

/// Is the gene in a case-insensitive gene set?
pub fn gene_in_set(gene: &str, set: &[&str]) -> bool {
    set.iter().any(|g| g.eq_ignore_ascii_case(gene))
}
