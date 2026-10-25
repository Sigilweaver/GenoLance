//! Which reference assembly a VCF header describes, so a sample made against
//! GRCh37 is not joined with a GRCh38 ClinVar release.
//!
//! Real ClinVar VCFs carry no `##contig` lines, only `##reference=GRCh38`
//! (or `GRCh37`), so the assembly is taken from two sources:
//!
//! 1. the lengths of the primary-assembly contigs the header declares, matched
//!    against the tables below;
//! 2. the name in the `##reference` line (`GRCh37`, `hg19`, `GRCh38`, `hg38`,
//!    `T2T-CHM13`).
//!
//! When the two disagree, or neither is present, the assembly is unknown.

use std::fmt;

/// A reference assembly this tool can tell apart.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Assembly {
    Grch37,
    Grch38,
    /// Named in `##reference` only: no length table is built in.
    T2tChm13,
}

impl fmt::Display for Assembly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl Assembly {
    pub(crate) fn as_str(&self) -> &'static str {
        match self {
            Assembly::Grch37 => "GRCh37",
            Assembly::Grch38 => "GRCh38",
            Assembly::T2tChm13 => "T2T-CHM13",
        }
    }

    /// Parse a label as stored or given on the command line.
    pub(crate) fn from_label(s: &str) -> Option<Assembly> {
        match s.trim().to_ascii_lowercase().as_str() {
            "grch37" | "hg19" | "b37" => Some(Assembly::Grch37),
            "grch38" | "hg38" => Some(Assembly::Grch38),
            "t2t-chm13" | "chm13" | "t2t" => Some(Assembly::T2tChm13),
            _ => None,
        }
    }
}

/// Primary-assembly chromosome lengths (chr1-22, X, Y), GRCh38.p14.
/// Source: NCBI Assembly GCF_000001405.40 (GRCh38.p14) assembly report,
/// column "Sequence-Length" of the assembled-molecule rows
/// (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_assembly_report.txt).
/// The mitochondrial contig (16569 in both assemblies) does not tell them apart
/// and is left out.
const GRCH38: &[(&str, u64)] = &[
    ("1", 248_956_422),
    ("2", 242_193_529),
    ("3", 198_295_559),
    ("4", 190_214_555),
    ("5", 181_538_259),
    ("6", 170_805_979),
    ("7", 159_345_973),
    ("8", 145_138_636),
    ("9", 138_394_717),
    ("10", 133_797_422),
    ("11", 135_086_622),
    ("12", 133_275_309),
    ("13", 114_364_328),
    ("14", 107_043_718),
    ("15", 101_991_189),
    ("16", 90_338_345),
    ("17", 83_257_441),
    ("18", 80_373_285),
    ("19", 58_617_616),
    ("20", 64_444_167),
    ("21", 46_709_983),
    ("22", 50_818_468),
    ("X", 156_040_895),
    ("Y", 57_227_415),
];

/// Primary-assembly chromosome lengths (chr1-22, X, Y), GRCh37.p13.
/// Source: NCBI Assembly GCF_000001405.25 (GRCh37.p13) assembly report,
/// column "Sequence-Length" of the assembled-molecule rows
/// (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_assembly_report.txt).
const GRCH37: &[(&str, u64)] = &[
    ("1", 249_250_621),
    ("2", 243_199_373),
    ("3", 198_022_430),
    ("4", 191_154_276),
    ("5", 180_915_260),
    ("6", 171_115_067),
    ("7", 159_138_663),
    ("8", 146_364_022),
    ("9", 141_213_431),
    ("10", 135_534_747),
    ("11", 135_006_516),
    ("12", 133_851_895),
    ("13", 115_169_878),
    ("14", 107_349_540),
    ("15", 102_531_392),
    ("16", 90_354_753),
    ("17", 81_195_210),
    ("18", 78_077_248),
    ("19", 59_128_983),
    ("20", 63_025_520),
    ("21", 48_129_895),
    ("22", 51_304_566),
    ("X", 155_270_560),
    ("Y", 59_373_566),
];

fn count_matches(table: &[(&str, u64)], contigs: &[(String, u64)]) -> usize {
    contigs
        .iter()
        .filter(|(name, len)| table.iter().any(|(n, l)| n == name && l == len))
        .count()
}

/// The assembly named by a `##reference` value, if any.
fn from_reference_text(text: &str) -> Option<Assembly> {
    let t = text.to_ascii_lowercase();
    if t.contains("grch38") || t.contains("hg38") {
        Some(Assembly::Grch38)
    } else if t.contains("grch37") || t.contains("hg19") || t.contains("b37") {
        Some(Assembly::Grch37)
    } else if t.contains("chm13") || t.contains("t2t") {
        Some(Assembly::T2tChm13)
    } else {
        None
    }
}

/// What a header says about its assembly.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Inferred {
    pub assembly: Option<Assembly>,
    /// Why the assembly is what it is, or why it is unknown.
    pub why: String,
}

/// Infer the assembly from the declared contigs (normalized name, length) and
/// the `##reference` text.
pub(crate) fn infer(contigs: &[(String, u64)], reference: Option<&str>) -> Inferred {
    let m38 = count_matches(GRCH38, contigs);
    let m37 = count_matches(GRCH37, contigs);
    let by_length = match (m38, m37) {
        (0, 0) => None,
        (n, 0) => Some((Assembly::Grch38, n)),
        (0, n) => Some((Assembly::Grch37, n)),
        _ => {
            return Inferred {
                assembly: None,
                why: "contig lengths match both GRCh37 and GRCh38".into(),
            }
        }
    };
    let by_name = reference.and_then(from_reference_text);
    match (by_length, by_name) {
        (Some((a, n)), Some(b)) if a != b => Inferred {
            assembly: None,
            why: format!(
                "contig lengths ({n} contigs) say {a} but ##reference names {b}"
            ),
        },
        (Some((a, n)), _) => Inferred {
            assembly: Some(a),
            why: format!("{n} contig length(s) match {a}"),
        },
        (None, Some(b)) => Inferred {
            assembly: Some(b),
            why: format!("##reference names {b}"),
        },
        (None, None) => Inferred {
            assembly: None,
            why: match reference {
                Some(r) => format!(
                    "no contig lengths match GRCh37 or GRCh38 and ##reference ({r}) does not name an assembly"
                ),
                None => "no ##contig lengths match GRCh37 or GRCh38 and there is no ##reference line"
                    .into(),
            },
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn contigs(t: &[(&str, u64)]) -> Vec<(String, u64)> {
        t.iter().map(|(n, l)| (n.to_string(), *l)).collect()
    }

    #[test]
    fn lengths_pick_the_assembly() {
        let c38 = contigs(&[("1", 248_956_422), ("MT", 16_569)]);
        assert_eq!(infer(&c38, None).assembly, Some(Assembly::Grch38));
        let c37 = contigs(&[("1", 249_250_621), ("17", 81_195_210)]);
        assert_eq!(infer(&c37, None).assembly, Some(Assembly::Grch37));
    }

    #[test]
    fn real_clinvar_header_shape_uses_the_reference_line() {
        // Real NCBI ClinVar VCFs have `##reference=GRCh38` and no ##contig lines.
        let i = infer(&[], Some("GRCh38"));
        assert_eq!(i.assembly, Some(Assembly::Grch38));
        assert_eq!(infer(&[], Some("GRCh37")).assembly, Some(Assembly::Grch37));
        assert_eq!(
            infer(&[], Some("file:///ref/hg19.fa")).assembly,
            Some(Assembly::Grch37)
        );
    }

    #[test]
    fn disagreement_and_absence_are_unknown() {
        let c37 = contigs(&[("1", 249_250_621)]);
        assert_eq!(infer(&c37, Some("GRCh38")).assembly, None);
        assert_eq!(infer(&[], None).assembly, None);
        assert_eq!(
            infer(&contigs(&[("1", 1000)]), Some("custom.fa")).assembly,
            None
        );
    }

    #[test]
    fn labels_parse() {
        assert_eq!(Assembly::from_label("hg38"), Some(Assembly::Grch38));
        assert_eq!(Assembly::from_label("GRCh37"), Some(Assembly::Grch37));
        assert_eq!(Assembly::from_label("mm10"), None);
    }
}
