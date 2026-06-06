//! Build scalar indices on a GenoLance store so region / sample / gene
//! lookups stop doing full table scans.
//!
//! Lance supports two scalar index flavors we care about:
//!
//! * **BTree** — sorted index. Best on high-cardinality columns with
//!   range / equality queries. Used here for `pos` and `variation_id`.
//! * **Bitmap** — one bitset per distinct value. Best on low-cardinality
//!   columns (dozens to low hundreds of distinct values). Used here for
//!   `chrom`, `sample_name`, `gene_symbol`, `clinical_significance`.
//!
//! All GenoLance query paths issue SQL-like predicates against these
//! columns (`chrom = 'chr17' AND pos >= … AND pos <= …`,
//! `sample_name = 'Nathan'`, `gene_symbol IN (…)`), so indexing them
//! gives every existing subcommand a speedup with no caller changes.

use anyhow::{Context, Result};
use lancedb::index::{
    scalar::{BTreeIndexBuilder, BitmapIndexBuilder},
    Index,
};

use genolance_core::schema::{CLINVAR_TABLE, SAMPLES_TABLE, VARIANTS_TABLE};
use genolance_core::store::Store;

/// Build all recommended indices. Safe to run repeatedly: Lance's
/// `create_index` replaces any existing index on the same column(s).
pub async fn run(store_path: &str) -> Result<()> {
    let store = Store::open(store_path).await?;
    let var_tables = store.variants.table_names().execute().await?;
    let root_tables = store.conn.table_names().execute().await?;

    if var_tables.iter().any(|n| n == VARIANTS_TABLE) {
        let t = store
            .variants
            .open_table(VARIANTS_TABLE)
            .execute()
            .await
            .with_context(|| format!("opening {VARIANTS_TABLE}"))?;
        eprintln!("[index] {VARIANTS_TABLE}: building BTree(pos)");
        t.create_index(&["pos"], Index::BTree(BTreeIndexBuilder::default()))
            .execute()
            .await?;
        eprintln!("[index] {VARIANTS_TABLE}: building Bitmap(chrom)");
        t.create_index(&["chrom"], Index::Bitmap(BitmapIndexBuilder::default()))
            .execute()
            .await?;
        eprintln!("[index] {VARIANTS_TABLE}: building Bitmap(sample_name)");
        t.create_index(
            &["sample_name"],
            Index::Bitmap(BitmapIndexBuilder::default()),
        )
        .execute()
        .await?;
    } else {
        eprintln!("[index] no {VARIANTS_TABLE} table; skipping");
    }

    if var_tables.iter().any(|n| n == CLINVAR_TABLE) {
        let t = store
            .variants
            .open_table(CLINVAR_TABLE)
            .execute()
            .await
            .with_context(|| format!("opening {CLINVAR_TABLE}"))?;
        eprintln!("[index] {CLINVAR_TABLE}: building BTree(pos)");
        t.create_index(&["pos"], Index::BTree(BTreeIndexBuilder::default()))
            .execute()
            .await?;
        eprintln!("[index] {CLINVAR_TABLE}: building Bitmap(chrom)");
        t.create_index(&["chrom"], Index::Bitmap(BitmapIndexBuilder::default()))
            .execute()
            .await?;
        eprintln!("[index] {CLINVAR_TABLE}: building Bitmap(gene_symbol)");
        t.create_index(
            &["gene_symbol"],
            Index::Bitmap(BitmapIndexBuilder::default()),
        )
        .execute()
        .await?;
        eprintln!("[index] {CLINVAR_TABLE}: building Bitmap(clinical_significance)");
        t.create_index(
            &["clinical_significance"],
            Index::Bitmap(BitmapIndexBuilder::default()),
        )
        .execute()
        .await?;
    } else {
        eprintln!("[index] no {CLINVAR_TABLE} table; skipping");
    }

    if root_tables.iter().any(|n| n == SAMPLES_TABLE) {
        let t = store.conn.open_table(SAMPLES_TABLE).execute().await?;
        eprintln!("[index] {SAMPLES_TABLE}: building Bitmap(sample_name)");
        t.create_index(
            &["sample_name"],
            Index::Bitmap(BitmapIndexBuilder::default()),
        )
        .execute()
        .await?;
    }

    // Enumerate what ended up on disk so the user sees the final state.
    for name in var_tables.iter().chain(root_tables.iter()) {
        let t = if var_tables.contains(name) {
            store.variants.open_table(name).execute().await?
        } else {
            store.conn.open_table(name).execute().await?
        };
        let indices = t.list_indices().await?;
        if indices.is_empty() {
            eprintln!("[index] {name}: no indices");
        } else {
            for idx in indices {
                eprintln!(
                    "[index] {name}: {} on {:?} (type = {:?})",
                    idx.name, idx.columns, idx.index_type
                );
            }
        }
    }

    Ok(())
}
