mod common;

use common::*;
use genolance_core::Error;
use genolance_variants::ingest::{self, IngestOptions};

const GOOD: &str = "chr1\t100\t.\tA\tG\t50\tPASS\tAF=0.5\tGT:AD:DP:GQ\t0/1:5,5:10:40";

#[tokio::test]
async fn good_file_ingests_plain_and_bgzf() {
    let t = Tmp::new();
    let text = vcf("S1", &[GOOD]);
    write_plain(&t.path("a.vcf"), &text);
    write_bgzf(&t.path("b.vcf.gz"), text.replace("S1", "S2").as_bytes());
    ingest_file(&t.store(), &t.path("a.vcf")).await.unwrap();
    ingest_file(&t.store(), &t.path("b.vcf.gz")).await.unwrap();
    assert_eq!(calls_count(&t.store(), None).await, 2);
}

// S1-3
#[tokio::test]
async fn bgzf_without_eof_block_is_rejected() {
    let t = Tmp::new();
    write_bgzf(&t.path("a.vcf.gz"), vcf("S1", &[GOOD]).as_bytes());
    let mut bytes = std::fs::read(t.path("a.vcf.gz")).unwrap();
    bytes.truncate(bytes.len() - 28);
    std::fs::write(t.path("a.vcf.gz"), bytes).unwrap();
    let e = ingest_file(&t.store(), &t.path("a.vcf.gz"))
        .await
        .unwrap_err();
    assert!(matches!(ingest_error(&e), Error::Truncated { .. }));
    assert_eq!(calls_count(&t.store(), None).await, 0);
}

// S1-3
#[tokio::test]
async fn partial_last_line_is_rejected_plain_and_bgzf() {
    let t = Tmp::new();
    // Cut inside the sample column: "0/1:3" would otherwise be stored as GQ.
    let mut text = vcf("S1", &[GOOD]);
    text.pop(); // drop final newline
    text.truncate(text.len() - 10);
    write_plain(&t.path("a.vcf"), &text);
    write_bgzf(&t.path("a.vcf.gz"), text.as_bytes());
    for f in ["a.vcf", "a.vcf.gz"] {
        let e = ingest_file(&t.store(), &t.path(f)).await.unwrap_err();
        assert!(
            matches!(ingest_error(&e), Error::Truncated { .. }),
            "{f}: {e:#}"
        );
    }
    assert_eq!(calls_count(&t.store(), None).await, 0);
    assert_eq!(samples_count(&t.store()).await, 0);
}

// S1-3: a record cut before the sample column still ends in a newline.
#[tokio::test]
async fn record_with_missing_sample_column_is_rejected() {
    let t = Tmp::new();
    let text = vcf("S1", &["chr1\t100\t.\tA\tG\t50\tPASS\tAF=0.5"]);
    write_plain(&t.path("a.vcf"), &text);
    let e = ingest_file(&t.store(), &t.path("a.vcf")).await.unwrap_err();
    let msg = format!("{e:#}");
    assert!(
        matches!(
            ingest_error(&e),
            Error::MalformedRecord { field, .. } if field == "sample columns"
        ),
        "{msg}"
    );
}

// S1-4
#[tokio::test]
async fn same_sample_twice_is_rejected_unless_replace() {
    let t = Tmp::new();
    write_plain(&t.path("a.vcf"), &vcf("S1", &[GOOD]));
    ingest_file(&t.store(), &t.path("a.vcf")).await.unwrap();
    let e = ingest_file(&t.store(), &t.path("a.vcf")).await.unwrap_err();
    assert!(matches!(ingest_error(&e), Error::DuplicateSample { .. }));
    assert_eq!(calls_count(&t.store(), None).await, 1);

    let opts = IngestOptions::new().with_replace(true);
    ingest::run(&t.store(), &[t.path("a.vcf")], &opts)
        .await
        .unwrap();
    assert_eq!(calls_count(&t.store(), None).await, 1);
    assert_eq!(samples_count(&t.store()).await, 1);
}

// S1-4
#[tokio::test]
async fn sample_override_with_two_files_is_rejected() {
    let t = Tmp::new();
    write_plain(&t.path("a.vcf"), &vcf("S1", &[GOOD]));
    write_plain(&t.path("b.vcf"), &vcf("S2", &[GOOD]));
    let opts = IngestOptions::new().with_sample_override("X");
    let e = ingest::run(&t.store(), &[t.path("a.vcf"), t.path("b.vcf")], &opts)
        .await
        .unwrap_err();
    assert!(matches!(ingest_error(&e), Error::InvalidArgument(_)));
    assert_eq!(calls_count(&t.store(), None).await, 0);
}

// S1-4: two files whose header names collide.
#[tokio::test]
async fn colliding_header_names_are_rejected() {
    let t = Tmp::new();
    write_plain(&t.path("a.vcf"), &vcf("default", &[GOOD]));
    write_plain(&t.path("b.vcf"), &vcf("default", &[GOOD]));
    let e = ingest::run(
        &t.store(),
        &[t.path("a.vcf"), t.path("b.vcf")],
        &IngestOptions::new(),
    )
    .await
    .unwrap_err();
    assert!(matches!(ingest_error(&e), Error::DuplicateSample { .. }));
}

// S2-13: a failure part-way leaves nothing behind and a retry works.
#[tokio::test]
async fn failed_ingest_is_rolled_back() {
    let t = Tmp::new();
    let bad = vcf("S1", &[GOOD, "chr1\t200\t.\tA\tG,T\t50\tPASS\t.\tGT\t0/3"]);
    write_plain(&t.path("a.vcf"), &bad);
    let e = ingest_file(&t.store(), &t.path("a.vcf")).await.unwrap_err();
    assert!(matches!(
        ingest_error(&e),
        Error::MalformedRecord { field, .. } if field == "FORMAT/GT"
    ));
    assert_eq!(calls_count(&t.store(), None).await, 0);
    assert_eq!(samples_count(&t.store()).await, 0);

    write_plain(&t.path("a.vcf"), &vcf("S1", &[GOOD]));
    ingest_file(&t.store(), &t.path("a.vcf")).await.unwrap();
    assert_eq!(calls_count(&t.store(), None).await, 1);
}

// S2-1
#[tokio::test]
async fn multi_sample_vcf_needs_explicit_column() {
    let t = Tmp::new();
    let mut text = String::from(HEADER_1);
    text.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tkid\tmom\n");
    text.push_str("chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t1/1\n");
    write_plain(&t.path("trio.vcf"), &text);
    let e = ingest_file(&t.store(), &t.path("trio.vcf"))
        .await
        .unwrap_err();
    assert!(matches!(ingest_error(&e), Error::MultiSample { .. }));

    let opts = IngestOptions::new().with_sample_column("mom");
    ingest::run(&t.store(), &[t.path("trio.vcf")], &opts)
        .await
        .unwrap();
    assert_eq!(
        calls_count(&t.store(), Some("sample_name = 'mom'")).await,
        1
    );
    assert_eq!(
        calls_count(&t.store(), Some("sample_name = 'mom' AND gt_raw = '1/1'")).await,
        1
    );

    let opts = IngestOptions::new().with_sample_column("dad");
    let e = ingest::run(&t.store(), &[t.path("trio.vcf")], &opts)
        .await
        .unwrap_err();
    assert!(matches!(ingest_error(&e), Error::UnknownSample { .. }));
}

// S1-10, S1-11, S2-15: malformed values are errors, not missing values.
#[tokio::test]
async fn malformed_values_are_errors() {
    let cases: &[(&str, &str)] = &[
        ("GT=0/1x", "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1x"),
        ("GT=0-1", "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0\\1"),
        (
            "AD abc",
            "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:AD\t0/1:12,abc",
        ),
        ("INFO AC", "chr1\t100\t.\tA\tG\t50\tPASS\tAC=1,x\tGT\t0/1"),
        ("QUAL", "chr1\t100\t.\tA\tG\thigh\tPASS\t.\tGT\t0/1"),
        ("empty ALT", "chr1\t100\t.\tA\tG,,T\t50\tPASS\t.\tGT\t0/2"),
        (
            "AF count",
            "chr1\t100\t.\tA\tG,T\t50\tPASS\tAF=0.5\tGT\t0/1",
        ),
    ];
    for (name, line) in cases {
        let t = Tmp::new();
        write_plain(&t.path("a.vcf"), &vcf("S1", &[line]));
        let r = ingest_file(&t.store(), &t.path("a.vcf")).await;
        assert!(r.is_err(), "{name}: expected an error");
        assert_eq!(calls_count(&t.store(), None).await, 0, "{name}");
    }
}

// S1-12: AF Number=R skips the REF value.
#[tokio::test]
async fn af_number_r_is_offset_by_one() {
    let t = Tmp::new();
    let text = vcf(
        "S1",
        &["chr1\t100\t.\tA\tG,T\t50\tPASS\tAF=0.5,0.25,0.25\tGT\t0/1"],
    )
    .replace("ID=AF,Number=A", "ID=AF,Number=R");
    write_plain(&t.path("a.vcf"), &text);
    ingest_file(&t.store(), &t.path("a.vcf")).await.unwrap();
    assert_eq!(
        calls_count(&t.store(), Some("alt_allele = 'G' AND allele_freq = 0.25")).await,
        1
    );
}

// S2-2: Float-typed GQ and DP are kept as written.
#[tokio::test]
async fn float_typed_gq_is_preserved() {
    let t = Tmp::new();
    let text = vcf("S1", &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:GQ\t0/1:33.5"])
        .replace("ID=GQ,Number=1,Type=Integer", "ID=GQ,Number=1,Type=Float");
    write_plain(&t.path("a.vcf"), &text);
    ingest_file(&t.store(), &t.path("a.vcf")).await.unwrap();
    let lines = export_lines(&t.store(), "S1", &t.0.join("out.vcf")).await;
    assert!(lines[0].ends_with("GT:GQ\t0/1:33.5"), "{}", lines[0]);
}

// S1-8, S1-9, S1-18, S1-19: INFO and FORMAT text round-trips exactly.
#[tokio::test]
async fn info_and_format_text_roundtrip_exactly() {
    let t = Tmp::new();
    let line = "chr1\t100\t.\tA\tG\t50\tPASS\tNOTE=a%3Bb;MQRankSum=.;AF=0.50\tGT:FT:AD\t0/1:LowGQ;LowDP:1,2";
    write_plain(&t.path("a.vcf"), &vcf("S1", &[line]));
    ingest_file(&t.store(), &t.path("a.vcf")).await.unwrap();
    let lines = export_lines(&t.store(), "S1", &t.0.join("out.vcf")).await;
    assert_eq!(lines, vec![line.to_string()]);
}

// S2-3: records without ALT are skipped, the rest are stored.
#[tokio::test]
async fn reference_only_records_are_skipped_and_rest_kept() {
    let t = Tmp::new();
    let text = vcf("S1", &["chr1\t50\t.\tA\t.\t50\tPASS\t.\tGT\t0/0", GOOD]);
    write_plain(&t.path("a.vcf"), &text);
    ingest_file(&t.store(), &t.path("a.vcf")).await.unwrap();
    assert_eq!(calls_count(&t.store(), None).await, 1);
}

fn clinvar_text(sig: &str) -> String {
    format!(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n\
##INFO=<ID=CLNSIG,Number=.,Type=String,Description=\"sig\">\n\
##INFO=<ID=GENEINFO,Number=1,Type=String,Description=\"gene\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t100\t1\tA\tG\t.\t.\tCLNSIG={sig};GENEINFO=GENE1:1\n"
    )
}

// S1-5, S2-11
#[tokio::test]
async fn clinvar_reingest_requires_replace() {
    let t = Tmp::new();
    write_plain(&t.path("old.vcf"), &clinvar_text("Pathogenic"));
    write_plain(&t.path("new.vcf"), &clinvar_text("Benign"));
    ingest_file(&t.store(), &t.path("old.vcf")).await.unwrap();
    let e = ingest_file(&t.store(), &t.path("new.vcf"))
        .await
        .unwrap_err();
    assert!(matches!(
        ingest_error(&e),
        Error::ClinVarReleasePresent { existing: 1, .. }
    ));
    assert_eq!(clinvar_count(&t.store()).await, 1);

    let opts = IngestOptions::new().with_replace(true);
    ingest::run(&t.store(), &[t.path("new.vcf")], &opts)
        .await
        .unwrap();
    assert_eq!(clinvar_count(&t.store()).await, 1);
}

// S2-11: a sites-only VCF that is not ClinVar is not put in the clinvar table.
#[tokio::test]
async fn sites_only_non_clinvar_vcf_is_rejected() {
    let t = Tmp::new();
    let text = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"af\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t100\t.\tA\tG\t.\t.\tAF=0.1\n";
    write_plain(&t.path("gnomad.vcf"), text);
    let e = ingest_file(&t.store(), &t.path("gnomad.vcf"))
        .await
        .unwrap_err();
    assert!(matches!(ingest_error(&e), Error::UnrecognizedInput { .. }));
    assert_eq!(clinvar_count(&t.store()).await, 0);

    let opts = IngestOptions::new().with_clinvar(true);
    ingest::run(&t.store(), &[t.path("gnomad.vcf")], &opts)
        .await
        .unwrap();
    assert_eq!(clinvar_count(&t.store()).await, 1);
}

// S2-11: a sample VCF with "clinvar" in its name is still a sample.
#[tokio::test]
async fn sample_vcf_named_clinvar_is_a_sample() {
    let t = Tmp::new();
    write_plain(
        &t.path("patient1.clinvar_annotated.vcf"),
        &vcf("S1", &[GOOD]),
    );
    ingest_file(&t.store(), &t.path("patient1.clinvar_annotated.vcf"))
        .await
        .unwrap();
    assert_eq!(calls_count(&t.store(), None).await, 1);
    assert_eq!(clinvar_count(&t.store()).await, 0);
}

// S1-3 for ClinVar: truncated release.
#[tokio::test]
async fn truncated_clinvar_is_rejected_and_rolled_back() {
    let t = Tmp::new();
    let mut text = clinvar_text("Pathogenic");
    text.pop();
    write_bgzf(&t.path("clinvar.vcf.gz"), text.as_bytes());
    let e = ingest_file(&t.store(), &t.path("clinvar.vcf.gz"))
        .await
        .unwrap_err();
    assert!(matches!(ingest_error(&e), Error::Truncated { .. }));
    assert_eq!(clinvar_count(&t.store()).await, 0);
}

// S12: plain gzip and non-gzip files named .gz get their own message.
#[tokio::test]
async fn plain_gzip_and_fake_gz_are_not_reported_as_truncated() {
    let t = Tmp::new();
    // An empty gzip member: magic, deflate, no FEXTRA, stored empty block.
    let plain_gzip: [u8; 20] = [
        0x1f, 0x8b, 0x08, 0x00, 0, 0, 0, 0, 0, 0xff, 0x03, 0x00, 0, 0, 0, 0, 0, 0, 0, 0,
    ];
    std::fs::write(t.path("plain.vcf.gz"), plain_gzip).unwrap();
    let e = ingest_file(&t.store(), &t.path("plain.vcf.gz"))
        .await
        .unwrap_err();
    let msg = format!("{e:#}");
    assert!(msg.contains("plain gzip"), "{msg}");
    assert!(msg.contains("bgzip"), "{msg}");
    assert!(!msg.contains("truncated"), "{msg}");

    write_plain(&t.path("fake.vcf.gz"), &vcf("S1", &[GOOD]));
    let e = ingest_file(&t.store(), &t.path("fake.vcf.gz"))
        .await
        .unwrap_err();
    assert!(format!("{e:#}").contains("gzip magic"), "{e:#}");
}
