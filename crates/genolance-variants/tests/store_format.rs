mod common;

use std::path::Path;

use common::*;
use genolance_core::schema::{samples_schema, META_FORMAT_VERSION};
use genolance_core::store::Store;
use genolance_core::Error;
use genolance_variants::ingest::{self, IngestOptions};

async fn open_err(path: &str) -> Error {
    match Store::open(path).await {
        Err(e) => e,
        Ok(_) => panic!("opening {path} should have failed"),
    }
}

#[tokio::test]
async fn opening_a_missing_path_fails_and_creates_nothing() {
    let t = Tmp::new();
    let p = t.path("nope");
    let e = open_err(&p).await;
    assert!(matches!(e.root(), Error::StoreNotFound { .. }), "{e}");
    assert_eq!(e.exit_code(), 4);
    assert!(!Path::new(&p).exists());
}

#[tokio::test]
async fn a_store_without_a_format_version_must_be_rebuilt() {
    let t = Tmp::new();
    let p = t.store();
    std::fs::create_dir_all(format!("{p}/variants")).unwrap();
    let conn = lancedb::connect(&p).execute().await.unwrap();
    conn.create_empty_table("samples", samples_schema())
        .execute()
        .await
        .unwrap();
    let e = open_err(&p).await;
    assert!(matches!(e.root(), Error::StoreTooOld { .. }), "{e}");
    assert!(format!("{e}").contains("Rebuild it"), "{e}");
    assert_eq!(e.exit_code(), 4);
    // Ingest refuses it too, rather than writing new-format rows next to old ones.
    write_plain(
        &t.path("a.vcf"),
        &vcf("S1", &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"]),
    );
    let e = ingest_file(&p, &t.path("a.vcf")).await.unwrap_err();
    assert!(matches!(e.root(), Error::StoreTooOld { .. }), "{e}");
}

#[tokio::test]
async fn a_newer_store_is_rejected() {
    let t = Tmp::new();
    write_plain(
        &t.path("a.vcf"),
        &vcf("S1", &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"]),
    );
    ingest_file(&t.store(), &t.path("a.vcf")).await.unwrap();
    let store = Store::open_or_create(&t.store()).await.unwrap();
    store.set_meta(META_FORMAT_VERSION, "99").await.unwrap();
    drop(store);
    let e = open_err(&t.store()).await;
    assert!(matches!(e.root(), Error::StoreTooNew { .. }), "{e}");
    assert_eq!(e.exit_code(), 4);
}

// S2-13: an interrupted ingest leaves a flag that every command refuses.
#[tokio::test]
async fn an_incomplete_sample_blocks_opening_until_replaced() {
    let t = Tmp::new();
    write_plain(
        &t.path("a.vcf"),
        &vcf("S1", &["chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1"]),
    );
    ingest_file(&t.store(), &t.path("a.vcf")).await.unwrap();
    assert!(Store::open(&t.store()).await.is_ok());

    let conn = lancedb::connect(&t.store()).execute().await.unwrap();
    let samples = conn.open_table("samples").execute().await.unwrap();
    samples
        .update()
        .only_if("sample_name = 'S1'")
        .column("ingest_complete", "false")
        .execute()
        .await
        .unwrap();
    let e = open_err(&t.store()).await;
    assert!(matches!(e.root(), Error::IncompleteIngest { .. }), "{e}");

    let opts = IngestOptions::new().with_replace(true);
    ingest::run(&t.store(), &[t.path("a.vcf")], &opts)
        .await
        .unwrap();
    assert!(Store::open(&t.store()).await.is_ok());
}

// S1-8/9/18: the FORMAT and sample columns come back as written, including
// trailing fields the record leaves out.
#[tokio::test]
async fn format_and_sample_columns_roundtrip_as_written() {
    let t = Tmp::new();
    let l1 = "chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT:DP:GQ\t0/1:5";
    let l2 = "chr1\t200\t.\tC\tT\t50\tPASS\t.\tGT:FT\t0/1:LowGQ;LowDP";
    write_plain(&t.path("a.vcf"), &vcf("S1", &[l1, l2]));
    ingest_file(&t.store(), &t.path("a.vcf")).await.unwrap();
    let lines = export_lines(&t.store(), "S1", &t.0.join("out.vcf")).await;
    assert_eq!(lines, vec![l1.to_string(), l2.to_string()]);
}
