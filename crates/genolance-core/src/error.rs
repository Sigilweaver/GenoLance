//! Error type shared by the GenoLance crates.
//!
//! Library functions return [`Result`]. Match on the variant (use
//! [`Error::root`] to look through added context) instead of on message
//! text. New kinds are added over time.

use std::fmt::Display;

/// Result alias used by the GenoLance libraries.
pub type Result<T> = std::result::Result<T, Error>;

/// Why an input or a request was rejected, or what went wrong underneath.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// The file ends early: a missing BGZF end-of-file block, or a last
    /// line without its newline.
    #[error("{path}: truncated input: {detail}")]
    Truncated { path: String, detail: String },

    /// A record (or one field of it) is not what the decoder expects.
    #[error("{chrom}:{pos}: malformed {field}: {reason}")]
    MalformedRecord {
        chrom: String,
        pos: u64,
        field: String,
        reason: String,
    },

    /// A sample with this name is already stored.
    #[error("sample '{sample}' is already in the store; pass --replace to overwrite it")]
    DuplicateSample { sample: String },

    /// A multi-sample VCF without an explicit choice of sample column.
    #[error("{path}: VCF has {} samples ({}); pass --sample-column <name> to choose the one to ingest", .samples.len(), .samples.join(", "))]
    MultiSample { path: String, samples: Vec<String> },

    /// A requested sample is not present.
    #[error("sample '{name}' not found; available: {}", .available.join(", "))]
    UnknownSample {
        name: String,
        available: Vec<String>,
    },

    /// The file is not a format the command can use as given.
    #[error("{path}: {reason}")]
    UnrecognizedInput { path: String, reason: String },

    /// The store already holds a ClinVar release.
    #[error("the clinvar table already holds {existing} rows; ingesting a second release would leave both in the store. Pass --replace to replace the stored release")]
    ClinVarReleasePresent { existing: usize },

    /// A sample's ingest did not finish (interrupted or failed).
    #[error(
        "sample '{sample}' has an incomplete ingest; re-run `genolance ingest --replace` for it"
    )]
    IncompleteSample { sample: String },

    /// The sample and ClinVar were not made against the same reference, or
    /// that could not be checked.
    #[error("reference check failed for sample '{sample}': {detail}")]
    ReferenceMismatch { sample: String, detail: String },

    /// The path is not a GenoLance store.
    #[error("{path}: not a GenoLance store (run `genolance ingest` to create one)")]
    StoreNotFound { path: String },

    /// The store was written by a newer GenoLance.
    #[error("store format version {found} is newer than this build supports ({supported}); upgrade genolance")]
    StoreTooNew { found: String, supported: u32 },

    /// The store predates format versions or uses an older one.
    #[error("store {path} has format version {found}, this build needs {supported}; it cannot be upgraded in place. Rebuild it: ingest your VCFs and the ClinVar release into a new store directory with `genolance ingest -s <new-dir> ...`")]
    StoreTooOld {
        path: String,
        found: String,
        supported: u32,
    },

    /// A request that cannot be honoured as written.
    #[error("{0}")]
    InvalidArgument(String),

    /// Added context around another error.
    #[error("{context}: {source}")]
    Context {
        context: String,
        #[source]
        source: Box<Error>,
    },

    /// A file system error.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// An error from the storage engine or the VCF reader.
    #[error("{0}")]
    Store(String),
}

impl Error {
    /// The innermost error, looking through [`Error::Context`] layers.
    pub fn root(&self) -> &Error {
        match self {
            Error::Context { source, .. } => source.root(),
            other => other,
        }
    }

    /// Wrap this error with a description of what was being done.
    pub fn context(self, context: impl Display) -> Error {
        Error::Context {
            context: context.to_string(),
            source: Box::new(self),
        }
    }

    /// Process exit code for the command-line tool.
    ///
    /// 1 internal, I/O or storage-engine error; 2 request cannot be honoured
    /// as written; 3 input data or store contents rejected; 4 store format
    /// version or layout. (clap itself exits 2 for usage errors.)
    pub fn exit_code(&self) -> i32 {
        match self.root() {
            Error::InvalidArgument(_) => 2,
            Error::Truncated { .. }
            | Error::MalformedRecord { .. }
            | Error::DuplicateSample { .. }
            | Error::MultiSample { .. }
            | Error::UnknownSample { .. }
            | Error::UnrecognizedInput { .. }
            | Error::ClinVarReleasePresent { .. }
            | Error::IncompleteSample { .. }
            | Error::ReferenceMismatch { .. } => 3,
            Error::StoreNotFound { .. } | Error::StoreTooNew { .. } | Error::StoreTooOld { .. } => {
                4
            }
            _ => 1,
        }
    }
}

impl From<lancedb::Error> for Error {
    fn from(e: lancedb::Error) -> Self {
        Error::Store(e.to_string())
    }
}

impl From<arrow_schema::ArrowError> for Error {
    fn from(e: arrow_schema::ArrowError) -> Self {
        Error::Store(e.to_string())
    }
}

impl From<std::string::FromUtf8Error> for Error {
    fn from(e: std::string::FromUtf8Error) -> Self {
        Error::Store(e.to_string())
    }
}

/// Add context to the error of a `Result`, like `anyhow::Context`.
pub trait ResultExt<T> {
    fn context<C: Display>(self, context: C) -> Result<T>;
    fn with_context<C: Display, F: FnOnce() -> C>(self, f: F) -> Result<T>;
}

impl<T, E: Into<Error>> ResultExt<T> for std::result::Result<T, E> {
    fn context<C: Display>(self, context: C) -> Result<T> {
        self.map_err(|e| e.into().context(context))
    }
    fn with_context<C: Display, F: FnOnce() -> C>(self, f: F) -> Result<T> {
        self.map_err(|e| e.into().context(f()))
    }
}
