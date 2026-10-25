use std::path::Path;
use std::sync::Arc;

use arrow_array::{Array, RecordBatch, StringArray};
use futures::TryStreamExt;
use lancedb::query::{ExecutableQuery, QueryBase};
use lancedb::Connection;

use crate::error::{Error, Result};
use crate::schema::{
    store_meta_schema, META_CLINVAR_STATE, META_FORMAT_VERSION, SAMPLES_TABLE,
    STORE_FORMAT_VERSION, STORE_META_TABLE,
};

/// A handle to an open GenoLance store.
///
/// The on-disk layout uses one LanceDB connection per omics layer:
///
/// ```text
/// <store>/
///   store_meta.lance/     <- root connection: format version, ClinVar release
///   samples.lance/        <- root connection: sample registry (shared spine)
///   variants/
///     calls.lance/        <- variants connection: called variant rows
///     clinvar.lance/      <- variants connection: ClinVar annotation rows
///   rna/                  <- future: genolance-rna layer
///   methyl/               <- future: genolance-methyl layer
/// ```
///
/// The store records its format version in `store_meta`. A store with
/// another version is rejected; it is not upgraded in place (see
/// [`Error::StoreTooOld`]).
pub struct Store {
    /// Root connection - holds the `samples` registry table.
    pub conn: Connection,
    /// Variants-layer connection - holds `calls` and `clinvar` tables.
    pub variants: Connection,
}

impl Store {
    /// Open an existing GenoLance store for reading.
    ///
    /// Fails if `root` is not a store, if its format version is not the one
    /// this build writes, or if a sample (or the ClinVar release) has an
    /// incomplete ingest. Nothing is created.
    pub async fn open(root: &str) -> Result<Self> {
        if !Path::new(root).is_dir() {
            return Err(Error::StoreNotFound {
                path: root.to_string(),
            });
        }
        let store = Self::connect(root).await?;
        match store.read_version(root).await? {
            Version::Current => {}
            Version::Empty => {
                return Err(Error::StoreNotFound {
                    path: root.to_string(),
                })
            }
        }
        store.check_complete().await?;
        Ok(store)
    }

    /// Open the store at `root` for ingest, creating it if it does not exist.
    /// An existing store must have the current format version. Incomplete
    /// samples are allowed here so `ingest --replace` can repair them.
    pub async fn open_or_create(root: &str) -> Result<Self> {
        // Ensure the subdirectory exists before LanceDB tries to connect.
        std::fs::create_dir_all(format!("{root}/variants"))?;
        let store = Self::connect(root).await?;
        if let Version::Empty = store.read_version(root).await? {
            store
                .set_meta(META_FORMAT_VERSION, &STORE_FORMAT_VERSION.to_string())
                .await?;
        }
        Ok(store)
    }

    async fn connect(root: &str) -> Result<Self> {
        let variants_path = format!("{root}/variants");
        if !Path::new(&variants_path).is_dir() {
            return Err(Error::StoreNotFound {
                path: root.to_string(),
            });
        }
        let conn = lancedb::connect(root).execute().await?;
        let variants = lancedb::connect(&variants_path).execute().await?;
        Ok(Self { conn, variants })
    }

    async fn read_version(&self, root: &str) -> Result<Version> {
        let root_tables = self.conn.table_names().execute().await?;
        if root_tables.iter().any(|n| n == STORE_META_TABLE) {
            let found = self.meta(META_FORMAT_VERSION).await?;
            let Some(found) = found else {
                return Err(Error::StoreTooOld {
                    path: root.to_string(),
                    found: "missing".into(),
                    supported: STORE_FORMAT_VERSION,
                });
            };
            return match found.parse::<u32>() {
                Ok(v) if v == STORE_FORMAT_VERSION => Ok(Version::Current),
                Ok(v) if v > STORE_FORMAT_VERSION => Err(Error::StoreTooNew {
                    found,
                    supported: STORE_FORMAT_VERSION,
                }),
                _ => Err(Error::StoreTooOld {
                    path: root.to_string(),
                    found,
                    supported: STORE_FORMAT_VERSION,
                }),
            };
        }
        let variant_tables = self.variants.table_names().execute().await?;
        if root_tables.is_empty() && variant_tables.is_empty() {
            Ok(Version::Empty)
        } else {
            // Tables but no store_meta: written before format versions.
            Err(Error::StoreTooOld {
                path: root.to_string(),
                found: "unversioned".into(),
                supported: STORE_FORMAT_VERSION,
            })
        }
    }

    /// Fail if a sample or the ClinVar release has an unfinished ingest.
    async fn check_complete(&self) -> Result<()> {
        let root_tables = self.conn.table_names().execute().await?;
        if root_tables.iter().any(|n| n == SAMPLES_TABLE) {
            let samples = self.conn.open_table(SAMPLES_TABLE).execute().await?;
            let batches: Vec<RecordBatch> = samples
                .query()
                .only_if("ingest_complete = false")
                .execute()
                .await?
                .try_collect()
                .await?;
            for b in &batches {
                let names = b
                    .column_by_name("sample_name")
                    .and_then(|c| c.as_any().downcast_ref::<StringArray>());
                let Some(names) = names else {
                    return Err(Error::missing_column(SAMPLES_TABLE, "sample_name"));
                };
                if names.len() > 0 {
                    return Err(Error::IncompleteIngest {
                        what: format!("sample '{}'", names.value(0)),
                    });
                }
            }
            // The flag column must exist with the right type, or the filter
            // above would have failed; check the type explicitly too.
            let schema = samples.schema().await?;
            let ok = schema
                .field_with_name("ingest_complete")
                .map(|f| f.data_type() == &arrow_schema::DataType::Boolean)
                .unwrap_or(false);
            if !ok {
                return Err(Error::missing_column(SAMPLES_TABLE, "ingest_complete"));
            }
        }
        if self.meta(META_CLINVAR_STATE).await?.as_deref() == Some("incomplete") {
            return Err(Error::IncompleteIngest {
                what: "the ClinVar release".into(),
            });
        }
        Ok(())
    }

    /// Read a metadata value.
    pub async fn meta(&self, key: &str) -> Result<Option<String>> {
        let root_tables = self.conn.table_names().execute().await?;
        if !root_tables.iter().any(|n| n == STORE_META_TABLE) {
            return Ok(None);
        }
        let t = self.conn.open_table(STORE_META_TABLE).execute().await?;
        let batches: Vec<RecordBatch> = t
            .query()
            .only_if(format!("key = '{}'", key.replace('\'', "''")))
            .execute()
            .await?
            .try_collect()
            .await?;
        for b in &batches {
            let vals = b
                .column_by_name("value")
                .and_then(|c| c.as_any().downcast_ref::<StringArray>())
                .ok_or_else(|| Error::missing_column(STORE_META_TABLE, "value"))?;
            if vals.len() > 0 {
                return Ok(Some(vals.value(0).to_string()));
            }
        }
        Ok(None)
    }

    /// Set a metadata value, replacing any previous value for the key.
    pub async fn set_meta(&self, key: &str, value: &str) -> Result<()> {
        let root_tables = self.conn.table_names().execute().await?;
        let t = if root_tables.iter().any(|n| n == STORE_META_TABLE) {
            self.conn.open_table(STORE_META_TABLE).execute().await?
        } else {
            self.conn
                .create_empty_table(STORE_META_TABLE, store_meta_schema())
                .execute()
                .await?
        };
        t.delete(&format!("key = '{}'", key.replace('\'', "''")))
            .await?;
        let batch = RecordBatch::try_new(
            store_meta_schema(),
            vec![
                Arc::new(StringArray::from(vec![key.to_string()])),
                Arc::new(StringArray::from(vec![value.to_string()])),
            ],
        )?;
        t.add(vec![batch]).execute().await?;
        Ok(())
    }

    /// Remove a metadata key (no error if it is absent).
    pub async fn clear_meta(&self, key: &str) -> Result<()> {
        let root_tables = self.conn.table_names().execute().await?;
        if root_tables.iter().any(|n| n == STORE_META_TABLE) {
            let t = self.conn.open_table(STORE_META_TABLE).execute().await?;
            t.delete(&format!("key = '{}'", key.replace('\'', "''")))
                .await?;
        }
        Ok(())
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

enum Version {
    Current,
    Empty,
}
