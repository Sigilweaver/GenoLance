mod common;

use std::sync::Arc;

use arrow_array::{RecordBatch, StringArray, UInt64Array};
use common::*;
use genolance_core::schema::clinvar_schema;
use genolance_variants::filters::Filters;

// S1-5: a store that already holds two ClinVar releases (written by an
// older version) must not answer from whichever row the scan returns last.
#[tokio::test]
async fn store_with_two_clinvar_releases_is_an_error() {
    let t = Tmp::new();
    write_plain(
        &t.path("s.vcf"),
        &vcf("S1", &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"]),
    );
    ingest_file(&t.store(), &t.path("s.vcf")).await.unwrap();
    let cv = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n\
##INFO=<ID=CLNSIG,Number=.,Type=String,Description=\"s\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t100\t1\tA\tG\t.\t.\tCLNSIG=Pathogenic\n";
    write_plain(&t.path("cv.vcf"), cv);
    ingest_file(&t.store(), &t.path("cv.vcf")).await.unwrap();

    // Simulate the older behavior: append a second release with a new class.
    let conn = lancedb::connect(&format!("{}/variants", t.store()))
        .execute()
        .await
        .unwrap();
    let table = conn.open_table("clinvar").execute().await.unwrap();
    let opt = |v: Option<&str>| Arc::new(StringArray::from(vec![v])) as _;
    let batch = RecordBatch::try_new(
        clinvar_schema(),
        vec![
            Arc::new(StringArray::from(vec!["chr1"])),
            Arc::new(UInt64Array::from(vec![100u64])),
            Arc::new(StringArray::from(vec!["A"])),
            Arc::new(StringArray::from(vec!["G"])),
            opt(Some("1")),
            opt(None),
            opt(Some("Benign")),
            opt(None),
            opt(None),
        ],
    )
    .unwrap();
    table.add(vec![batch]).execute().await.unwrap();

    let e = genolance_variants::join::run(&t.store(), None, &Filters::default())
        .await
        .unwrap_err();
    assert!(format!("{e:#}").contains("more than one release"), "{e:#}");
}
