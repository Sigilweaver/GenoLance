//! Shared helpers for integration tests: temp directories, synthetic VCF
//! text, BGZF writing and store inspection.
#![allow(dead_code)]

use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU32, Ordering};

use genolance_core::Error;
use genolance_variants::ingest::{self, IngestOptions};

static COUNTER: AtomicU32 = AtomicU32::new(0);

pub struct Tmp(pub PathBuf);

impl Tmp {
    pub fn new() -> Self {
        let n = COUNTER.fetch_add(1, Ordering::SeqCst);
        let p = std::env::temp_dir().join(format!("genolance-test-{}-{n}", std::process::id()));
        std::fs::create_dir_all(&p).unwrap();
        Tmp(p)
    }
    pub fn path(&self, name: &str) -> String {
        self.0.join(name).to_string_lossy().into_owned()
    }
    pub fn store(&self) -> String {
        self.path("store")
    }
}

impl Drop for Tmp {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

pub const HEADER_1: &str = "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=100000>\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"af\">\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"ac\">\n\
##INFO=<ID=NOTE,Number=1,Type=String,Description=\"note\">\n\
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"mq\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"ad\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"dp\">\n\
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"gq\">\n\
##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"pl\">\n\
##FORMAT=<ID=FT,Number=1,Type=String,Description=\"ft\">\n";

/// A VCF with one sample column named `name` and the given data lines
/// (each without a trailing newline).
pub fn vcf(name: &str, lines: &[&str]) -> String {
    let mut s = String::from(HEADER_1);
    s.push_str(&format!(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{name}\n"
    ));
    for l in lines {
        s.push_str(l);
        s.push('\n');
    }
    s
}

pub fn write_plain(path: &str, text: &str) {
    std::fs::write(path, text).unwrap();
}

/// Write `bytes` as a BGZF file that ends with the EOF marker block.
pub fn write_bgzf(path: &str, bytes: &[u8]) {
    let f = std::fs::File::create(path).unwrap();
    let mut w = noodles::bgzf::io::Writer::new(f);
    w.write_all(bytes).unwrap();
    w.finish().unwrap();
}

pub async fn ingest_file(store: &str, file: &str) -> genolance_core::Result<()> {
    ingest::run(store, &[file.to_string()], &IngestOptions::new())
        .await
        .map(|_| ())
}

pub fn ingest_error(e: &Error) -> &Error {
    e.root()
}

async fn open(store: &str, sub: &str, table: &str) -> Option<lancedb::Table> {
    let conn = lancedb::connect(&format!("{store}{sub}"))
        .execute()
        .await
        .ok()?;
    conn.open_table(table).execute().await.ok()
}

/// Number of rows in `variants/calls` matching `filter` (0 if the table is missing).
pub async fn calls_count(store: &str, filter: Option<&str>) -> usize {
    match open(store, "/variants", "calls").await {
        Some(t) => t.count_rows(filter.map(str::to_string)).await.unwrap(),
        None => 0,
    }
}

pub async fn clinvar_count(store: &str) -> usize {
    match open(store, "/variants", "clinvar").await {
        Some(t) => t.count_rows(None).await.unwrap(),
        None => 0,
    }
}

pub async fn samples_count(store: &str) -> usize {
    match open(store, "", "samples").await {
        Some(t) => t.count_rows(None).await.unwrap(),
        None => 0,
    }
}

/// Export a sample to a file and return the data lines (header dropped).
pub async fn export_lines(store: &str, sample: &str, out: &Path) -> Vec<String> {
    genolance_variants::export::run(store, sample, None, None, None, out.to_str())
        .await
        .unwrap();
    std::fs::read_to_string(out)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_string)
        .collect()
}
