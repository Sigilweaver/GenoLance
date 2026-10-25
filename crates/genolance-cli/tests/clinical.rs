//! End-to-end tests of the clinical commands through the `genolance` binary.

use std::path::PathBuf;
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU32, Ordering};

static COUNTER: AtomicU32 = AtomicU32::new(0);

struct Tmp(PathBuf);
impl Tmp {
    fn new() -> Self {
        let n = COUNTER.fetch_add(1, Ordering::SeqCst);
        let p = std::env::temp_dir().join(format!("genolance-cli-test-{}-{n}", std::process::id()));
        std::fs::create_dir_all(&p).unwrap();
        Tmp(p)
    }
    fn path(&self, n: &str) -> String {
        self.0.join(n).to_string_lossy().into_owned()
    }
    fn store(&self) -> String {
        self.path("store")
    }
}
impl Drop for Tmp {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn gl(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_genolance"))
        .args(args)
        .output()
        .unwrap()
}
fn ok(args: &[&str]) -> (String, String) {
    let o = gl(args);
    let (so, se) = (
        String::from_utf8_lossy(&o.stdout).into_owned(),
        String::from_utf8_lossy(&o.stderr).into_owned(),
    );
    assert!(o.status.success(), "{args:?} failed: {se}");
    (so, se)
}
fn fails(args: &[&str]) -> String {
    let o = gl(args);
    assert!(!o.status.success(), "{args:?} should have failed");
    String::from_utf8_lossy(&o.stderr).into_owned()
}

const SAMPLE_HEADER: &str = "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=100000>\n##contig=<ID=chrX,length=100000>\n##contig=<ID=chrM,length=20000>\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"af\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"dp\">\n";

fn sample_vcf(name: &str, lines: &[&str]) -> String {
    let mut s = String::from(SAMPLE_HEADER);
    s.push_str(&format!(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{name}\n"
    ));
    for l in lines {
        s.push_str(l);
        s.push('\n');
    }
    s
}

fn clinvar_vcf(lines: &[&str]) -> String {
    let mut s = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=100000>\n##contig=<ID=X,length=100000>\n##contig=<ID=MT,length=20000>\n\
##INFO=<ID=CLNSIG,Number=.,Type=String,Description=\"sig\">\n\
##INFO=<ID=GENEINFO,Number=1,Type=String,Description=\"gene\">\n\
##INFO=<ID=CLNDN,Number=.,Type=String,Description=\"dis\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    for l in lines {
        s.push_str(l);
        s.push('\n');
    }
    s
}

fn setup(t: &Tmp, sample_lines: &[&str], clinvar_lines: &[&str]) {
    std::fs::write(t.path("s.vcf"), sample_vcf("S1", sample_lines)).unwrap();
    std::fs::write(t.path("cv.vcf"), clinvar_vcf(clinvar_lines)).unwrap();
    ok(&["ingest", "-s", &t.store(), &t.path("s.vcf")]);
    ok(&["ingest", "-s", &t.store(), &t.path("cv.vcf")]);
}

const CV_1: &str = "1\t100\t1\tA\tG\t.\t.\tCLNSIG=Pathogenic;GENEINFO=GENE1:1;CLNDN=Disease_X";

// S1-6, S1-7, S1-15
#[test]
fn join_skips_no_calls_excludes_filtered_and_prints_missing_as_dot() {
    let t = Tmp::new();
    setup(
        &t,
        &[
            "chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT:DP\t0/1:12",
            "chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT\t./.",
            "chr1\t300\t.\tG\tA\t50\tLowQual\t.\tGT\t0/1",
        ],
        &[
            CV_1,
            "1\t200\t2\tC\tT\t.\t.\tCLNSIG=Pathogenic;GENEINFO=GENE2:2",
            "1\t300\t3\tG\tA\t.\t.\tCLNSIG=Pathogenic;GENEINFO=GENE3:3",
        ],
    );
    let (out, err) = ok(&["join", "-s", &t.store()]);
    assert!(out.contains("GENE1"), "carrier listed:\n{out}");
    assert!(!out.contains("GENE2"), "no-call must not be listed:\n{out}");
    assert!(
        !out.contains("GENE3"),
        "filtered call must not be listed:\n{out}"
    );
    assert!(out.contains("1 annotated variants"), "{out}");
    let line = out.lines().find(|l| l.contains("GENE1")).unwrap();
    assert!(
        !line.contains("0.0"),
        "missing QUAL must print as '.': {line}"
    );
    assert!(err.contains("FILTER is not PASS"), "{err}");
    assert!(err.contains("no-call"), "{err}");

    let (out, _) = ok(&["join", "-s", &t.store(), "--include-filtered"]);
    assert!(out.contains("GENE3") && out.contains("LowQual"), "{out}");
}

// S1-16
#[test]
fn empty_significance_exact_is_an_error() {
    let t = Tmp::new();
    setup(&t, &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"], &[CV_1]);
    for val in ["", ","] {
        let err = fails(&["join", "-s", &t.store(), "--significance-exact", val]);
        assert!(err.contains("significance-exact"), "{err}");
    }
}

// S2-12
#[test]
fn min_qual_counts_rows_with_missing_qual() {
    let t = Tmp::new();
    setup(&t, &["chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1"], &[CV_1]);
    let (out, err) = ok(&["join", "-s", &t.store(), "--min-qual", "30"]);
    assert!(!out.contains("GENE1"));
    assert!(err.contains("QUAL or DP is missing"), "{err}");
    let (out, _) = ok(&[
        "join",
        "-s",
        &t.store(),
        "--min-qual",
        "30",
        "--keep-missing",
    ]);
    assert!(out.contains("GENE1"), "{out}");
}

// S2-4
#[test]
fn chr_m_matches_clinvar_mt() {
    let t = Tmp::new();
    setup(
        &t,
        &["chrM\t3243\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"],
        &["MT\t3243\t9\tA\tG\t.\t.\tCLNSIG=Pathogenic;GENEINFO=MT-TL1:1"],
    );
    let (out, _) = ok(&["join", "-s", &t.store()]);
    assert!(out.contains("MT-TL1"), "{out}");
}

// S2-10
#[test]
fn hemizygous_alt_on_chrx_is_flagged_by_compound_het() {
    let t = Tmp::new();
    setup(
        &t,
        &["chrX\t500\t.\tC\tT\t50\tPASS\t.\tGT\t1"],
        &["X\t500\t7\tC\tT\t.\t.\tCLNSIG=Pathogenic;GENEINFO=OTC:1;CLNDN=OTC_deficiency"],
    );
    let (out, _) = ok(&["compound-het", "-s", &t.store(), "--sample", "S1"]);
    assert!(out.contains("OTC") && out.contains("HEMIZYGOUS"), "{out}");
}

// S2-10: ClinVar hits with no gene are counted
#[test]
fn compound_het_reports_hits_without_gene() {
    let t = Tmp::new();
    setup(
        &t,
        &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t1/1"],
        &["1\t100\t1\tA\tG\t.\t.\tCLNSIG=Pathogenic"],
    );
    let (_, err) = ok(&["compound-het", "-s", &t.store(), "--sample", "S1"]);
    assert!(err.contains("no gene symbol"), "{err}");
}

// S2-5
#[test]
fn query_chrom_matches_either_spelling() {
    let t = Tmp::new();
    std::fs::write(
        t.path("s.vcf"),
        sample_vcf("S1", &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"]),
    )
    .unwrap();
    ok(&["ingest", "-s", &t.store(), &t.path("s.vcf")]);
    for c in ["chr1", "1"] {
        let (out, _) = ok(&["query", "-s", &t.store(), "--chrom", c]);
        assert!(out.contains("1 rows"), "{c}: {out}");
    }
}

// S3-1
#[test]
fn unknown_query_output_is_a_usage_error() {
    let t = Tmp::new();
    std::fs::write(
        t.path("s.vcf"),
        sample_vcf("S1", &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"]),
    )
    .unwrap();
    ok(&["ingest", "-s", &t.store(), &t.path("s.vcf")]);
    let o = gl(&["query", "-s", &t.store(), "--output", "jsno"]);
    assert_eq!(o.status.code(), Some(2));
}

// S2-9
#[test]
fn pgx_genes_are_split_uppercased_and_checked() {
    let t = Tmp::new();
    setup(&t, &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"], &[CV_1]);
    let err = fails(&[
        "pgx",
        "-s",
        &t.store(),
        "--sample",
        "S1",
        "--genes",
        "gene1,nosuchgene",
    ]);
    assert!(err.contains("NOSUCHGENE"), "{err}");
    assert!(!err.contains("GENE1,"), "{err}");
    let (_, err) = ok(&[
        "pgx",
        "-s",
        &t.store(),
        "--sample",
        "S1",
        "--genes",
        "gene1",
    ]);
    assert!(err.contains("PGx ClinVar variants"), "{err}");
}

// S1-6 for compare
#[test]
fn compare_does_not_count_no_calls_as_variants() {
    let t = Tmp::new();
    std::fs::write(
        t.path("a.vcf"),
        sample_vcf(
            "A",
            &[
                "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1",
                "chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT\t./.",
            ],
        ),
    )
    .unwrap();
    std::fs::write(
        t.path("b.vcf"),
        sample_vcf(
            "B",
            &[
                "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1",
                "chr1\t300\t.\tC\tT\t50\tPASS\t.\tGT\t0",
            ],
        ),
    )
    .unwrap();
    ok(&["ingest", "-s", &t.store(), &t.path("a.vcf")]);
    ok(&["ingest", "-s", &t.store(), &t.path("b.vcf")]);
    let (out, err) = ok(&["compare", "-s", &t.store(), "A", "B"]);
    assert!(out.contains("total union of variant sites : 1"), "{out}");
    assert!(err.contains("no-call"), "{err}");
}

// S1-13
#[test]
fn annotate_uses_declared_number_for_per_alt_values() {
    let t = Tmp::new();
    std::fs::write(
        t.path("s.vcf"),
        sample_vcf("S1", &["chr1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t0/1"]),
    )
    .unwrap();
    ok(&["ingest", "-s", &t.store(), &t.path("s.vcf")]);
    let ann = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n\
##INFO=<ID=RA,Number=R,Type=Integer,Description=\"r\">\n\
##INFO=<ID=CSQ,Number=.,Type=String,Description=\"csq\">\n\
##INFO=<ID=BAD,Number=A,Type=Integer,Description=\"a\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t100\t.\tA\tG,T\t.\t.\tRA=1,2,3;CSQ=x,y;BAD=7,8\n";
    std::fs::write(t.path("ann.vcf"), ann).unwrap();
    let (out, _) = ok(&[
        "annotate",
        "-s",
        &t.store(),
        "--sample",
        "S1",
        "--vcf",
        &t.path("ann.vcf"),
        "--info",
        "RA,CSQ,BAD",
    ]);
    let line = out.lines().find(|l| l.starts_with("1\t100")).unwrap();
    let cols: Vec<&str> = line.split('\t').collect();
    assert_eq!(&cols[3..], &["T", "0/1", "3", "x,y", "8"], "{line}");
}

// S1: join never silently ignores the file it was given.
#[test]
fn join_against_is_not_ignored_when_a_release_is_stored() {
    let t = Tmp::new();
    setup(&t, &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"], &[CV_1]);
    std::fs::write(
        t.path("new.vcf"),
        clinvar_vcf(&["1\t100\t1\tA\tG\t.\t.\tCLNSIG=Benign;GENEINFO=GENE1:1"]),
    )
    .unwrap();
    let err = fails(&["join", "-s", &t.store(), "--against", &t.path("new.vcf")]);
    assert!(err.contains("ingest --replace"), "{err}");
    // Without --against and without a table there is nothing to join against.
    let t2 = Tmp::new();
    std::fs::write(
        t2.path("s.vcf"),
        sample_vcf("S1", &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"]),
    )
    .unwrap();
    ok(&["ingest", "-s", &t2.store(), &t2.path("s.vcf")]);
    let err = fails(&["join", "-s", &t2.store()]);
    assert!(err.contains("--against"), "{err}");
    // With a file and no table it ingests the file and joins.
    std::fs::write(t2.path("cv.vcf"), clinvar_vcf(&[CV_1])).unwrap();
    let (out, _) = ok(&["join", "-s", &t2.store(), "--against", &t2.path("cv.vcf")]);
    assert!(out.contains("GENE1"), "{out}");
}

fn annotation_vcf(lines: &[&str], final_newline: bool) -> String {
    let mut s = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"af\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    s.push_str(&lines.join("\n"));
    if final_newline {
        s.push('\n');
    }
    s
}

// S3: annotate checks the same things as the other commands.
#[test]
fn annotate_rejects_truncated_files_unknown_samples_and_undeclared_keys() {
    let t = Tmp::new();
    std::fs::write(
        t.path("s.vcf"),
        sample_vcf("S1", &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"]),
    )
    .unwrap();
    ok(&["ingest", "-s", &t.store(), &t.path("s.vcf")]);
    let line = "chr1\t100\t.\tA\tG\t.\t.\tAF=0.25";
    std::fs::write(t.path("full.vcf"), annotation_vcf(&[line], true)).unwrap();
    std::fs::write(t.path("cut.vcf"), annotation_vcf(&[line], false)).unwrap();
    let base = |vcf: String, sample: &str, info: &str| {
        let args = [
            "annotate",
            "-s",
            &t.store(),
            "--sample",
            sample,
            "--vcf",
            &vcf,
            "--info",
            info,
        ];
        gl(&args)
    };
    let o = base(t.path("full.vcf"), "S1", "AF");
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    let o = base(t.path("cut.vcf"), "S1", "AF");
    assert!(!o.status.success());
    assert!(String::from_utf8_lossy(&o.stderr).contains("truncated"));
    let o = base(t.path("full.vcf"), "Typo", "AF");
    assert!(!o.status.success());
    assert!(String::from_utf8_lossy(&o.stderr).contains("not found"));
    let o = base(t.path("full.vcf"), "S1", "AFX");
    assert!(!o.status.success());
    assert!(String::from_utf8_lossy(&o.stderr).contains("does not declare"));
}

// S4
#[test]
fn compare_rejects_unknown_samples_and_modes() {
    let t = Tmp::new();
    for (n, name) in [("a", "A"), ("b", "B")] {
        std::fs::write(
            t.path(&format!("{n}.vcf")),
            sample_vcf(name, &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"]),
        )
        .unwrap();
        ok(&["ingest", "-s", &t.store(), &t.path(&format!("{n}.vcf"))]);
    }
    let err = fails(&["compare", "-s", &t.store(), "A", "Typo"]);
    assert!(err.contains("not found"), "{err}");
    let err = fails(&["compare", "-s", &t.store(), "A", "B", "--mode", "nope"]);
    assert!(err.contains("unknown compare mode"), "{err}");
}

// S6, S7
#[test]
fn query_gene_is_case_insensitive_json_is_valid_and_arrow_is_gone() {
    let t = Tmp::new();
    setup(
        &t,
        &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"],
        &["1\t100\t1\tA\tG\t.\t.\tCLNSIG=Pathogenic;GENEINFO=C1orf112:1"],
    );
    let (out, _) = ok(&["query", "-s", &t.store(), "--gene", "c1ORF112"]);
    assert!(out.contains("1 rows"), "{out}");
    let err = fails(&["query", "-s", &t.store(), "--gene", "NOPE1"]);
    assert!(err.contains("no ClinVar rows"), "{err}");
    let o = gl(&["query", "-s", &t.store(), "--output", "arrow"]);
    assert_eq!(o.status.code(), Some(2));

    // JSON escaping: a sample name with a quote and a backslash.
    let t2 = Tmp::new();
    std::fs::write(
        t2.path("s.vcf"),
        sample_vcf("S1", &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t./."]),
    )
    .unwrap();
    ok(&[
        "ingest",
        "-s",
        &t2.store(),
        "--sample",
        "s\"Q\\x",
        &t2.path("s.vcf"),
    ]);
    let (out, _) = ok(&["query", "-s", &t2.store(), "--output", "json"]);
    let line = out.lines().next().unwrap();
    assert!(line.contains(r#""sample":"s\"Q\\x""#), "{line}");
    assert!(line.contains(r#""gt":"./.""#), "{line}");
}

// S14: an empty --significance would match everything
#[test]
fn empty_significance_substring_is_an_error() {
    let t = Tmp::new();
    setup(&t, &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"], &[CV_1]);
    let err = fails(&["join", "-s", &t.store(), "--significance", ""]);
    assert!(err.contains("--significance is empty"), "{err}");
}

// ---- reference assembly, ClinVar release, genes, alleles, store version ----

/// The header of the real NCBI ClinVar GRCh38 VCF (release 2026-10-19), cut to
/// the lines these tests need: it has `##reference=GRCh38` and NO `##contig`
/// lines. Source: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
const REAL_CLINVAR_HEADER: &str = "##fileformat=VCFv4.1\n\
##fileDate=2026-10-19\n\
##source=ClinVar\n\
##reference=GRCh38\n\
##ID=<Description=\"ClinVar Variation ID\">\n\
##INFO=<ID=CLNDN,Number=.,Type=String,Description=\"ClinVar's preferred disease name for the concept specified by disease identifiers in CLNDISDB\">\n\
##INFO=<ID=CLNSIG,Number=.,Type=String,Description=\"Aggregate germline classification for this single variant; multiple values are separated by a vertical bar\">\n\
##INFO=<ID=GENEINFO,Number=1,Type=String,Description=\"Gene(s) for the variant reported as gene symbol:gene id. The gene symbol and id are delimited by a colon (:) and each pair is delimited by a vertical bar (|)\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

fn real_clinvar(lines: &[&str]) -> String {
    let mut s = REAL_CLINVAR_HEADER.to_string();
    for l in lines {
        s.push_str(l);
        s.push('\n');
    }
    s
}

fn sample_with_contig(name: &str, contig: &str, lines: &[&str]) -> String {
    let mut s = format!(
        "##fileformat=VCFv4.2\n{contig}\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{name}\n"
    );
    for l in lines {
        s.push_str(l);
        s.push('\n');
    }
    s
}

const S38: &str = "##contig=<ID=chr1,length=248956422>\n";
const S37: &str = "##contig=<ID=chr1,length=249250621>\n";
const HIT: &str = "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1";
const CV_HIT: &str = "1\t100\t1\tA\tG\t.\t.\tCLNSIG=Pathogenic;GENEINFO=GENE1:1;CLNDN=Disease_X";

// B1: the real ClinVar header shape (no ##contig lines) works.
#[test]
fn real_clinvar_header_passes_for_a_grch38_sample() {
    let t = Tmp::new();
    std::fs::write(t.path("s.vcf"), sample_with_contig("S1", S38, &[HIT])).unwrap();
    std::fs::write(t.path("cv.vcf"), real_clinvar(&[CV_HIT])).unwrap();
    ok(&["ingest", "-s", &t.store(), &t.path("s.vcf")]);
    ok(&["ingest", "-s", &t.store(), &t.path("cv.vcf")]);
    let (out, err) = ok(&["join", "-s", &t.store()]);
    assert!(out.contains("GENE1"), "{out}");
    assert!(err.contains("ClinVar release: 2026-10-19"), "{err}");
    assert!(err.contains("assembly GRCh38"), "{err}");
}

// B1 / S1-14: a GRCh37 sample against GRCh38 ClinVar is refused, and the
// skip flag prints what the check would have said.
#[test]
fn grch37_sample_against_grch38_clinvar_is_refused() {
    let t = Tmp::new();
    std::fs::write(t.path("s.vcf"), sample_with_contig("S1", S37, &[HIT])).unwrap();
    std::fs::write(t.path("cv.vcf"), real_clinvar(&[CV_HIT])).unwrap();
    ok(&["ingest", "-s", &t.store(), &t.path("s.vcf")]);
    ok(&["ingest", "-s", &t.store(), &t.path("cv.vcf")]);
    let store = t.store();
    for cmd in [
        vec!["join", "-s", &store],
        vec!["compound-het", "-s", &store, "--sample", "S1"],
        vec!["pgx", "-s", &store, "--sample", "S1"],
    ] {
        let o = gl(&cmd);
        assert_eq!(o.status.code(), Some(3), "{cmd:?}");
        let err = String::from_utf8_lossy(&o.stderr);
        assert!(err.contains("GRCh37") && err.contains("GRCh38"), "{err}");
    }
    let (out, err) = ok(&["join", "-s", &t.store(), "--skip-reference-check"]);
    assert!(out.contains("GENE1"), "{out}");
    assert!(err.contains("WARNING: --skip-reference-check"), "{err}");
    assert!(err.contains("would have failed"), "{err}");
}

// B1: when the assembly cannot be told, the error names the fix, and
// declaring it at ingest resolves it.
#[test]
fn unknown_assembly_asks_for_a_declaration() {
    let t = Tmp::new();
    std::fs::write(t.path("s.vcf"), sample_with_contig("S1", "", &[HIT])).unwrap();
    std::fs::write(t.path("cv.vcf"), real_clinvar(&[CV_HIT])).unwrap();
    ok(&["ingest", "-s", &t.store(), &t.path("s.vcf")]);
    ok(&["ingest", "-s", &t.store(), &t.path("cv.vcf")]);
    let o = gl(&["join", "-s", &t.store()]);
    assert_eq!(o.status.code(), Some(3));
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("--assembly"), "{err}");
    ok(&[
        "ingest",
        "-s",
        &t.store(),
        "--replace",
        "--assembly",
        "GRCh38",
        &t.path("s.vcf"),
    ]);
    let (out, _) = ok(&["join", "-s", &t.store()]);
    assert!(out.contains("GENE1"), "{out}");
    let o = gl(&[
        "ingest",
        "-s",
        &t.path("x"),
        "--assembly",
        "mm10",
        &t.path("s.vcf"),
    ]);
    assert_eq!(o.status.code(), Some(2));
}

// S3-2: a ClinVar file without ##fileDate is shown as such, not as "unknown".
#[test]
fn release_without_a_file_date_is_stated_plainly() {
    let t = Tmp::new();
    std::fs::write(t.path("s.vcf"), sample_with_contig("S1", S38, &[HIT])).unwrap();
    let cv = real_clinvar(&[CV_HIT]).replace("##fileDate=2026-10-19\n", "");
    std::fs::write(t.path("cv.vcf"), cv).unwrap();
    ok(&["ingest", "-s", &t.store(), &t.path("s.vcf")]);
    ok(&["ingest", "-s", &t.store(), &t.path("cv.vcf")]);
    let (_, err) = ok(&["join", "-s", &t.store()]);
    assert!(err.contains("no ##fileDate"), "{err}");
}

// S1: `--against` with the file that is already loaded is fine.
#[test]
fn join_against_the_loaded_file_is_accepted() {
    let t = Tmp::new();
    std::fs::write(t.path("s.vcf"), sample_with_contig("S1", S38, &[HIT])).unwrap();
    std::fs::write(t.path("cv.vcf"), real_clinvar(&[CV_HIT])).unwrap();
    ok(&["ingest", "-s", &t.store(), &t.path("s.vcf")]);
    ok(&["ingest", "-s", &t.store(), &t.path("cv.vcf")]);
    let (out, _) = ok(&["join", "-s", &t.store(), "--against", &t.path("cv.vcf")]);
    assert!(out.contains("GENE1"), "{out}");
}

// S2-6 / S13: every gene of a record is searchable, ignoring case.
#[test]
fn every_gene_of_a_clinvar_record_is_searchable() {
    let t = Tmp::new();
    std::fs::write(t.path("s.vcf"), sample_with_contig("S1", S38, &[HIT])).unwrap();
    std::fs::write(
        t.path("cv.vcf"),
        real_clinvar(&["1\t100\t1\tA\tG\t.\t.\tCLNSIG=Pathogenic;GENEINFO=GENEX:1|C1orf112:2"]),
    )
    .unwrap();
    ok(&["ingest", "-s", &t.store(), &t.path("s.vcf")]);
    ok(&["ingest", "-s", &t.store(), &t.path("cv.vcf")]);
    for g in ["C1orf112", "c1ORF112", "genex"] {
        let (out, _) = ok(&["query", "-s", &t.store(), "--gene", g]);
        assert!(out.contains("1 rows"), "{g}: {out}");
    }
    let err = fails(&["query", "-s", &t.store(), "--gene", "NOPE_1"]);
    assert!(err.contains("NOPE_1"), "{err}");
    let (out, _) = ok(&["join", "-s", &t.store(), "--genes", "c1orf112"]);
    assert!(out.contains("1 annotated variants"), "{out}");
    assert!(out.contains("GENEX,C1orf112"), "{out}");
}

// S2-7 / S10: a padded multi-allelic split and lower-case bases still match.
#[test]
fn trimmed_and_lowercase_alleles_join_clinvar() {
    let t = Tmp::new();
    std::fs::write(
        t.path("s.vcf"),
        sample_with_contig(
            "S1",
            S38,
            &["chr1\t100\t.\tgca\tg,ga\t50\tPASS\t.\tGT\t0/2"],
        ),
    )
    .unwrap();
    std::fs::write(
        t.path("cv.vcf"),
        real_clinvar(&["1\t100\t1\tGC\tG\t.\t.\tCLNSIG=Pathogenic;GENEINFO=GENE1:1"]),
    )
    .unwrap();
    ok(&["ingest", "-s", &t.store(), &t.path("s.vcf")]);
    ok(&["ingest", "-s", &t.store(), &t.path("cv.vcf")]);
    let (out, _) = ok(&["join", "-s", &t.store()]);
    assert!(out.contains("1 annotated variants"), "{out}");
}

// Exit code 2 for a missing store (a wrong path is a request error).
#[test]
fn missing_store_exits_with_code_2_and_creates_nothing() {
    let t = Tmp::new();
    let o = gl(&["query", "-s", &t.path("nope")]);
    assert_eq!(o.status.code(), Some(2));
    assert!(!std::path::Path::new(&t.path("nope")).exists());
}

// BCF is not supported anywhere: say so instead of failing on a header parse.
#[test]
fn bcf_input_is_rejected_with_a_clear_message() {
    let t = Tmp::new();
    std::fs::write(t.path("x.bcf"), b"BCF\x02\x02").unwrap();
    let o = gl(&["ingest", "-s", &t.store(), &t.path("x.bcf")]);
    assert_eq!(o.status.code(), Some(3));
    let err = String::from_utf8_lossy(&o.stderr);
    assert!(err.contains("BCF input is not supported"), "{err}");
}

// Repeating --sample merges; --merge no longer exists.
#[test]
fn repeated_sample_merges_and_merge_flag_is_gone() {
    let t = Tmp::new();
    for (n, name) in [("a", "A"), ("b", "B")] {
        std::fs::write(
            t.path(&format!("{n}.vcf")),
            sample_vcf(name, &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"]),
        )
        .unwrap();
        ok(&["ingest", "-s", &t.store(), &t.path(&format!("{n}.vcf"))]);
    }
    let (out, _) = ok(&["export", "-s", &t.store(), "--sample", "A", "--sample", "B"]);
    let chrom = out.lines().find(|l| l.starts_with("#CHROM")).unwrap();
    assert!(chrom.ends_with("FORMAT\tA\tB"), "{chrom}");
    let o = gl(&["export", "-s", &t.store(), "--sample", "A", "--merge"]);
    assert_eq!(o.status.code(), Some(2));
}

// Exit codes: 2 for a request that cannot be honoured, 3 for rejected data.
#[test]
fn exit_codes_follow_the_error_kind() {
    let t = Tmp::new();
    setup(&t, &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"], &[CV_1]);
    let o = gl(&["join", "-s", &t.store(), "--significance-exact", ""]);
    assert_eq!(o.status.code(), Some(2), "empty filter list");
    let o = gl(&["ingest", "-s", &t.store(), &t.path("s.vcf")]);
    assert_eq!(o.status.code(), Some(2), "duplicate sample");
    let o = gl(&["compare", "-s", &t.store(), "S1", "Typo"]);
    assert_eq!(o.status.code(), Some(2), "unknown sample");
    let o = gl(&["query", "-s", &t.store(), "--gene", "NOPE1"]);
    assert_eq!(o.status.code(), Some(2), "unknown gene");
    let o = gl(&["export", "-s", &t.store(), "--sample", "S1", "--chrom", "9"]);
    assert_eq!(o.status.code(), Some(2), "empty region");
    // Rejected data: a VCF cut mid-line.
    let mut cut = sample_vcf("S2", &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"]);
    cut.pop();
    std::fs::write(t.path("cut.vcf"), cut).unwrap();
    let o = gl(&["ingest", "-s", &t.store(), &t.path("cut.vcf")]);
    assert_eq!(o.status.code(), Some(3), "truncated input");
}
