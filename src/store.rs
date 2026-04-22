use anyhow::Result;
use lancedb::Connection;

/// A handle to an open BioLance store (a LanceDB database directory).
pub struct Store {
    pub conn: Connection,
}

impl Store {
    /// Open an existing store or create a new one at `path`.
    pub async fn open(path: &str) -> Result<Self> {
        let conn = lancedb::connect(path).execute().await?;
        Ok(Self { conn })
    }

    /// List the names of all tables in this store.
    pub async fn table_names(&self) -> Result<Vec<String>> {
        Ok(self.conn.table_names().execute().await?)
    }
}
