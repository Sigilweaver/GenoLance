mod common;

use common::*;
use genolance_core::Error;
use genolance_variants::export;

async fn ingest_text(t: &Tmp, name: &str, text: &str) {
    let path = t.path(name);
    write_plain(&path, text);
    ingest_file(&t.store(), &path).await.unwrap();
}

fn cols(line: &str) -> Vec<&str> {
    line.split('\t').collect()
}

// S1-1: records split by `bcftools norm -m-` share (chrom, pos, ref).
#[tokio::test]
async fn split_biallelic_records_are_not_merged() {
    let t = Tmp::new();
    let l1 = "chr1\t100\t.\tA\tG\t50\tPASS\tAC=1\tGT:DP\t0/1:10";
    let l2 = "chr1\t100\t.\tA\tT\t30\tLowQual\tAC=2\tGT:DP\t1/1:8";
    ingest_text(&t, "a.vcf", &vcf("S1", &[l1, l2])).await;
    let mut lines = export_lines(&t.store(), "S1", &t.0.join("out.vcf")).await;
    lines.sort();
    let mut want = vec![l1.to_string(), l2.to_string()];
    want.sort();
    assert_eq!(lines, want);
}

#[tokio::test]
async fn multiallelic_record_roundtrips_as_one_line() {
    let t = Tmp::new();
    let l = "chr1\t200\t.\tA\tG,T\t50\tPASS\tAC=1,1\tGT:AD:PL\t1/2:0,5,6:9,8,7,6,5,4";
    ingest_text(&t, "a.vcf", &vcf("S1", &[l])).await;
    let lines = export_lines(&t.store(), "S1", &t.0.join("out.vcf")).await;
    assert_eq!(lines, vec![l.to_string()]);
}

// S1-1: two multi-allelic records with the same shape at one position are
// told apart by their source ordinal.
#[tokio::test]
async fn same_shape_multiallelic_records_stay_separate() {
    let t = Tmp::new();
    let l1 = "chr1\t200\t.\tA\tG,T\t50\tPASS\t.\tGT\t0/1";
    let l2 = "chr1\t200\t.\tA\tG,T\t60\tPASS\t.\tGT\t0/2";
    ingest_text(&t, "a.vcf", &vcf("S1", &[l1, l2])).await;
    let lines = export_lines(&t.store(), "S1", &t.0.join("out.vcf")).await;
    assert_eq!(lines, vec![l1.to_string(), l2.to_string()]);
}

// S1-17: the #CHROM line names the exported sample.
#[tokio::test]
async fn chrom_line_uses_stored_sample_name() {
    let t = Tmp::new();
    let path = t.path("a.vcf");
    write_plain(
        &path,
        &vcf("HG002", &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"]),
    );
    let opts = genolance_variants::ingest::IngestOptions::new().with_sample_override("sampleA");
    genolance_variants::ingest::run(&t.store(), &[path], &opts)
        .await
        .unwrap();
    let out = t.0.join("out.vcf");
    export::run(&t.store(), "sampleA", None, None, None, out.to_str())
        .await
        .unwrap();
    let text = std::fs::read_to_string(&out).unwrap();
    let chrom_line = text.lines().find(|l| l.starts_with("#CHROM")).unwrap();
    assert!(chrom_line.ends_with("FORMAT\tsampleA"), "{chrom_line}");
}

// S2-8
#[tokio::test]
async fn unknown_sample_is_an_error() {
    let t = Tmp::new();
    ingest_text(
        &t,
        "a.vcf",
        &vcf("SampleB", &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"]),
    )
    .await;
    let out = t.0.join("out.vcf");
    let e = export::run(&t.store(), "sampleB", None, None, None, out.to_str())
        .await
        .unwrap_err();
    assert!(matches!(ingest_error(&e), Error::UnknownSample { .. }));
    assert!(!out.exists());

    let names = vec!["SampleB".to_string(), "typo".to_string()];
    let e = export::run_merge(&t.store(), &names, None, None, None, out.to_str())
        .await
        .unwrap_err();
    assert!(matches!(ingest_error(&e), Error::UnknownSample { .. }));
}

// S1-2
#[tokio::test]
async fn merge_remaps_ad_pl_and_merges_filter_info() {
    let t = Tmp::new();
    ingest_text(
        &t,
        "a.vcf",
        &vcf(
            "A",
            &["chr1\t100\t.\tA\tT\t50\tPASS\tAC=1\tGT:AD:DP\t0/1:5,6:11"],
        ),
    )
    .await;
    ingest_text(
        &t,
        "b.vcf",
        &vcf(
            "B",
            &["chr1\t100\t.\tA\tG,T\t40\tLowQual\tAC=1,1\tGT:AD:PL\t1/2:3,10,12:0,10,20,30,40,50"],
        ),
    )
    .await;
    let out = t.0.join("m.vcf");
    export::run_merge(
        &t.store(),
        &["A".to_string(), "B".to_string()],
        None,
        None,
        None,
        out.to_str(),
    )
    .await
    .unwrap();
    let text = std::fs::read_to_string(&out).unwrap();
    let line = text.lines().find(|l| !l.starts_with('#')).unwrap();
    let c = cols(line);
    assert_eq!(c[4], "T,G");
    assert_eq!(c[5], ".", "QUAL differs between samples");
    assert_eq!(c[6], "LowQual", "non-PASS filter wins");
    assert_eq!(
        c[7], ".",
        "INFO differs, so it is not copied from one sample"
    );
    assert_eq!(c[9], "0/1:5,6,.:11:.:.");
    assert_eq!(c[10], "2/1:3,12,10:.:.:0,30,50,10,40,20");
}

// S1-1 for merge: split records of one sample cannot be merged unambiguously.
#[tokio::test]
async fn merge_rejects_split_records_of_one_sample() {
    let t = Tmp::new();
    ingest_text(
        &t,
        "a.vcf",
        &vcf(
            "A",
            &[
                "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1",
                "chr1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t0/1",
            ],
        ),
    )
    .await;
    ingest_text(
        &t,
        "b.vcf",
        &vcf("B", &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"]),
    )
    .await;
    let out = t.0.join("m.vcf");
    let e = export::run_merge(
        &t.store(),
        &["A".to_string(), "B".to_string()],
        None,
        None,
        None,
        out.to_str(),
    )
    .await
    .unwrap_err();
    assert!(matches!(ingest_error(&e), Error::MalformedRecord { .. }));
    assert!(!out.exists());
}

// S5: --chrom matches either spelling, and an empty region is an error.
#[tokio::test]
async fn export_chrom_matches_either_spelling_and_empty_region_errors() {
    let t = Tmp::new();
    let l = "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1";
    ingest_text(&t, "a.vcf", &vcf("S1", &[l])).await;
    for c in ["chr1", "1"] {
        let out = t.0.join(format!("out-{c}.vcf"));
        export::run(&t.store(), "S1", Some(c), None, None, out.to_str())
            .await
            .unwrap();
        let text = std::fs::read_to_string(&out).unwrap();
        assert!(text.lines().any(|x| x == l), "{c}");
    }
    let out = t.0.join("none.vcf");
    let e = export::run(&t.store(), "S1", Some("chr9"), None, None, out.to_str())
        .await
        .unwrap_err();
    assert!(format!("{e:#}").contains("region filter"), "{e:#}");
    assert!(!out.exists());
}
