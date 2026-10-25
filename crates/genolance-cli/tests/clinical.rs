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
    let (out, err) = ok(&["join", "-s", &t.store(), "--against", &t.path("cv.vcf")]);
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

    let (out, _) = ok(&[
        "join",
        "-s",
        &t.store(),
        "--against",
        &t.path("cv.vcf"),
        "--include-filtered",
    ]);
    assert!(out.contains("GENE3") && out.contains("LowQual"), "{out}");
}

// S1-16
#[test]
fn empty_significance_exact_is_an_error() {
    let t = Tmp::new();
    setup(&t, &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"], &[CV_1]);
    for val in ["", ","] {
        let err = fails(&[
            "join",
            "-s",
            &t.store(),
            "--against",
            &t.path("cv.vcf"),
            "--significance-exact",
            val,
        ]);
        assert!(err.contains("significance-exact"), "{err}");
    }
}

// S2-12
#[test]
fn min_qual_counts_rows_with_missing_qual() {
    let t = Tmp::new();
    setup(&t, &["chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1"], &[CV_1]);
    let (out, err) = ok(&[
        "join",
        "-s",
        &t.store(),
        "--against",
        &t.path("cv.vcf"),
        "--min-qual",
        "30",
    ]);
    assert!(!out.contains("GENE1"));
    assert!(err.contains("QUAL or DP is missing"), "{err}");
    let (out, _) = ok(&[
        "join",
        "-s",
        &t.store(),
        "--against",
        &t.path("cv.vcf"),
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
    let (out, _) = ok(&["join", "-s", &t.store(), "--against", &t.path("cv.vcf")]);
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
