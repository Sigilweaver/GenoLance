//! Error type shared by the GenoLance crates.
//!
//! Library functions return [`Result`]. Match on the variant (use
//! [`Error::root`] to look through added context) instead of on message
//! text. The enum and every struct variant are `#[non_exhaustive]`, so kinds
//! and fields can be added without breaking callers.
//!
//! [`Error::exit_code`] is the single place that maps a kind to a process
//! exit code: 1 internal, I/O or storage engine; 2 a request that cannot be
//! honoured as written (bad argument, unknown sample or gene, missing table,
//! empty region, duplicate or conflicting target); 3 input data or stored
//! data rejected as malformed, truncated or inconsistent; 4 the store's
//! format version or layout is not supported.

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
    #[non_exhaustive]
    Truncated { path: String, detail: String },

    /// A record (or one field of it) is not what the decoder expects.
    #[error("{chrom}:{pos}: malformed {field}: {reason}")]
    #[non_exhaustive]
    MalformedRecord {
        chrom: String,
        pos: u64,
        field: String,
        reason: String,
    },

    /// The file is not a format the command can use as given.
    #[error("{path}: {reason}")]
    #[non_exhaustive]
    UnrecognizedInput { path: String, reason: String },

    /// A sample with this name is already stored.
    #[error("sample '{sample}' is already in the store; pass --replace to overwrite it")]
    #[non_exhaustive]
    DuplicateSample { sample: String },

    /// A multi-sample VCF without an explicit choice of sample column.
    #[error("{path}: VCF has {} samples ({}); pass --sample-column <name> to choose the one to ingest", .samples.len(), .samples.join(", "))]
    #[non_exhaustive]
    MultiSample { path: String, samples: Vec<String> },

    /// A requested sample is not present.
    #[error("sample '{name}' not found; available: {}", .available.join(", "))]
    #[non_exhaustive]
    UnknownSample {
        name: String,
        available: Vec<String>,
    },

    /// A requested gene has no ClinVar rows.
    #[error("no ClinVar rows for gene(s): {} (symbols match case-insensitively); check the spelling", .symbols.join(", "))]
    #[non_exhaustive]
    UnknownGene { symbols: Vec<String> },

    /// The store already holds a ClinVar release.
    #[error("the clinvar table already holds {existing} rows; ingesting a second release would leave both in the store. Pass --replace to replace the stored release")]
    #[non_exhaustive]
    ClinVarReleasePresent { existing: usize },

    /// The ClinVar table has two rows for one variant.
    #[error("{chrom}:{pos} {reference}>{alt}: {detail}")]
    #[non_exhaustive]
    ClinVarConflict {
        chrom: String,
        pos: u64,
        reference: String,
        alt: String,
        detail: String,
    },

    /// A command needs a table the store does not have yet.
    #[error("the store has no {what}; {hint}")]
    #[non_exhaustive]
    MissingData { what: String, hint: String },

    /// A region filter matched nothing.
    #[error("no variants of {what} match the region filter ({region}); check the contig name and positions")]
    #[non_exhaustive]
    RegionEmpty { what: String, region: String },

    /// An ingest did not finish (interrupted or failed).
    #[error("{what} has an incomplete ingest (it was interrupted or failed); re-run `genolance ingest --replace` for it")]
    #[non_exhaustive]
    IncompleteIngest { what: String },

    /// The sample and ClinVar were not made against the same reference, or
    /// that could not be checked.
    #[error("reference check failed for sample '{sample}': {detail}")]
    #[non_exhaustive]
    ReferenceMismatch { sample: String, detail: String },

    /// The path is not a GenoLance store.
    #[error("{path}: not a GenoLance store (run `genolance ingest` to create one)")]
    #[non_exhaustive]
    StoreNotFound { path: String },

    /// The store was written by a newer GenoLance.
    #[error("store format version {found} is newer than this build supports ({supported}); upgrade genolance (a store with a newer minor version can be read but not written to)")]
    #[non_exhaustive]
    StoreTooNew { found: String, supported: String },

    /// The store predates format versions or uses an older one.
    #[error("store {path} has format version {found}, this build needs {supported}; it cannot be upgraded in place. Rebuild it: ingest your VCFs and the ClinVar release into a new store directory with `genolance ingest -s <new-dir> ...`")]
    #[non_exhaustive]
    StoreTooOld {
        path: String,
        found: String,
        supported: String,
    },

    /// A table lacks a column this version reads (or has it with another
    /// type): the store layout does not match this build.
    #[error("table '{table}' has no column '{column}' of the expected type; the store layout does not match this version of genolance")]
    #[non_exhaustive]
    MissingColumn { table: String, column: String },

    /// A request that cannot be honoured as written (bad argument value).
    #[error("{0}")]
    InvalidArgument(String),

    /// Added context around another error.
    #[error("{context}: {source}")]
    #[non_exhaustive]
    Context {
        context: String,
        #[source]
        source: Box<Error>,
    },

    /// A file system error.
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    /// An error from the storage engine (Lance, Arrow) or a background task;
    /// the original error is the `source()`.
    #[error("storage error: {0}")]
    Storage(#[source] Box<dyn std::error::Error + Send + Sync>),
}

impl Error {
    /// A storage-engine error with a message and no underlying error type.
    pub fn storage(message: impl Into<String>) -> Error {
        Error::Storage(message.into().into())
    }

    /// A required table is missing; `hint` says what to run.
    pub fn missing_table(table: &str, hint: &str) -> Error {
        Error::MissingData {
            what: format!("'{table}' table"),
            hint: hint.to_string(),
        }
    }

    // Constructors: the struct variants are `#[non_exhaustive]`, so other
    // crates build them through these.

    pub fn truncated(path: impl Into<String>, detail: impl Into<String>) -> Error {
        Error::Truncated {
            path: path.into(),
            detail: detail.into(),
        }
    }

    pub fn malformed_record(
        chrom: impl Into<String>,
        pos: u64,
        field: impl Into<String>,
        reason: impl Into<String>,
    ) -> Error {
        Error::MalformedRecord {
            chrom: chrom.into(),
            pos,
            field: field.into(),
            reason: reason.into(),
        }
    }

    pub fn unrecognized_input(path: impl Into<String>, reason: impl Into<String>) -> Error {
        Error::UnrecognizedInput {
            path: path.into(),
            reason: reason.into(),
        }
    }

    pub fn duplicate_sample(sample: impl Into<String>) -> Error {
        Error::DuplicateSample {
            sample: sample.into(),
        }
    }

    pub fn multi_sample(path: impl Into<String>, samples: Vec<String>) -> Error {
        Error::MultiSample {
            path: path.into(),
            samples,
        }
    }

    pub fn unknown_sample(name: impl Into<String>, available: Vec<String>) -> Error {
        Error::UnknownSample {
            name: name.into(),
            available,
        }
    }

    pub fn unknown_gene(symbols: Vec<String>) -> Error {
        Error::UnknownGene { symbols }
    }

    pub fn clinvar_release_present(existing: usize) -> Error {
        Error::ClinVarReleasePresent { existing }
    }

    pub fn clinvar_conflict(
        chrom: impl Into<String>,
        pos: u64,
        reference: impl Into<String>,
        alt: impl Into<String>,
        detail: impl Into<String>,
    ) -> Error {
        Error::ClinVarConflict {
            chrom: chrom.into(),
            pos,
            reference: reference.into(),
            alt: alt.into(),
            detail: detail.into(),
        }
    }

    pub fn region_empty(what: impl Into<String>, region: impl Into<String>) -> Error {
        Error::RegionEmpty {
            what: what.into(),
            region: region.into(),
        }
    }

    pub fn incomplete_ingest(what: impl Into<String>) -> Error {
        Error::IncompleteIngest { what: what.into() }
    }

    pub fn reference_mismatch(sample: impl Into<String>, detail: impl Into<String>) -> Error {
        Error::ReferenceMismatch {
            sample: sample.into(),
            detail: detail.into(),
        }
    }

    pub fn store_not_found(path: impl Into<String>) -> Error {
        Error::StoreNotFound { path: path.into() }
    }

    pub fn store_too_new(found: impl Into<String>, supported: impl Into<String>) -> Error {
        Error::StoreTooNew {
            found: found.into(),
            supported: supported.into(),
        }
    }

    pub fn store_too_old(
        path: impl Into<String>,
        found: impl Into<String>,
        supported: impl Into<String>,
    ) -> Error {
        Error::StoreTooOld {
            path: path.into(),
            found: found.into(),
            supported: supported.into(),
        }
    }

    pub fn missing_column(table: impl Into<String>, column: impl Into<String>) -> Error {
        Error::MissingColumn {
            table: table.into(),
            column: column.into(),
        }
    }

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

    /// Process exit code for the command-line tool (see the module docs).
    pub fn exit_code(&self) -> i32 {
        match self.root() {
            Error::InvalidArgument(_)
            | Error::UnknownSample { .. }
            | Error::UnknownGene { .. }
            | Error::DuplicateSample { .. }
            | Error::MultiSample { .. }
            | Error::ClinVarReleasePresent { .. }
            | Error::MissingData { .. }
            | Error::RegionEmpty { .. }
            | Error::StoreNotFound { .. } => 2,
            Error::Truncated { .. }
            | Error::MalformedRecord { .. }
            | Error::UnrecognizedInput { .. }
            | Error::ClinVarConflict { .. }
            | Error::IncompleteIngest { .. }
            | Error::ReferenceMismatch { .. } => 3,
            Error::StoreTooNew { .. } | Error::StoreTooOld { .. } | Error::MissingColumn { .. } => {
                4
            }
            _ => 1,
        }
    }
}

impl From<lancedb::Error> for Error {
    fn from(e: lancedb::Error) -> Self {
        Error::Storage(Box::new(e))
    }
}

impl From<arrow_schema::ArrowError> for Error {
    fn from(e: arrow_schema::ArrowError) -> Self {
        Error::Storage(Box::new(e))
    }
}

impl From<std::string::FromUtf8Error> for Error {
    fn from(e: std::string::FromUtf8Error) -> Self {
        Error::Storage(Box::new(e))
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exit_codes_follow_the_documented_classes() {
        // 2: the request cannot be honoured as written.
        assert_eq!(Error::InvalidArgument("x".into()).exit_code(), 2);
        assert_eq!(Error::unknown_sample("s", vec![]).exit_code(), 2);
        assert_eq!(Error::unknown_gene(vec!["G".into()]).exit_code(), 2);
        assert_eq!(Error::duplicate_sample("s").exit_code(), 2);
        assert_eq!(Error::missing_table("calls", "ingest").exit_code(), 2);
        assert_eq!(Error::region_empty("sample 's'", "chrom 9").exit_code(), 2);
        assert_eq!(Error::store_not_found("p").exit_code(), 2);
        // 3: input or stored data rejected.
        assert_eq!(Error::truncated("f", "x").exit_code(), 3);
        assert_eq!(Error::malformed_record("1", 5, "GT", "x").exit_code(), 3);
        assert_eq!(Error::unrecognized_input("f", "x").exit_code(), 3);
        assert_eq!(
            Error::clinvar_conflict("1", 5, "A", "G", "x").exit_code(),
            3
        );
        // 4: store format or layout.
        assert_eq!(Error::store_too_new("9.0", "1.0").exit_code(), 4);
        assert_eq!(
            Error::store_too_old("p", "unversioned", "1.0").exit_code(),
            4
        );
        assert_eq!(Error::missing_column("calls", "pos").exit_code(), 4);
        // 1: everything else, including context-wrapped storage errors.
        assert_eq!(Error::storage("boom").exit_code(), 1);
        assert_eq!(
            Error::storage("boom").context("while reading").exit_code(),
            1
        );
    }

    #[test]
    fn context_keeps_the_kind_reachable() {
        let e = Error::truncated("f.vcf.gz", "no EOF block").context("ingesting f.vcf.gz");
        assert!(matches!(e.root(), Error::Truncated { .. }));
        assert_eq!(e.exit_code(), 3);
        assert!(e.to_string().contains("ingesting f.vcf.gz"));
    }
}
