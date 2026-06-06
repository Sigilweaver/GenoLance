use anyhow::Result;
use lancedb::Connection;

/// A handle to an open GenoLance store.
///
/// The on-disk layout uses one LanceDB connection per omics layer:
///
/// ```text
/// <store>/
///   samples.lance/        ← root connection: sample registry (shared spine)
///   variants/
///     calls.lance/        ← variants connection: called variant rows
///     clinvar.lance/      ← variants connection: ClinVar annotation rows
///   rna/                  ← future: genolance-rna layer
///   methyl/               ← future: genolance-methyl layer
/// ```
pub struct Store {
    /// Root connection — holds the `samples` registry table.
    pub conn: Connection,
    /// Variants-layer connection — holds `calls` and `clinvar` tables.
    pub variants: Connection,
}

impl Store {
    /// Open (or create) a GenoLance store at `root`.
    ///
    /// Creates `<root>/variants/` if it doesn't exist.
    pub async fn open(root: &str) -> Result<Self> {
        let variants_path = format!("{}/variants", root);
        // Ensure the subdirectory exists before LanceDB tries to connect.
        std::fs::create_dir_all(&variants_path)?;
        let conn = lancedb::connect(root).execute().await?;
        let variants = lancedb::connect(&variants_path).execute().await?;
        Ok(Self { conn, variants })
    }

    /// List the names of all tables across all layers.
    pub async fn table_names(&self) -> Result<Vec<String>> {
        let mut names: Vec<String> = self
            .conn
            .table_names()
            .execute()
            .await?
            .into_iter()
            .collect();
        for name in self.variants.table_names().execute().await? {
            names.push(format!("variants/{}", name));
        }
        Ok(names)
    }
}
