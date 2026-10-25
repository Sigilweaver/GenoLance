//! Error type shared by the GenoLance crates.
//!
//! Library code returns these as the source of an `anyhow::Error`, so callers
//! can `downcast_ref::<genolance_core::Error>()` and match on the kind instead
//! of on message text. New kinds are added over time.

/// Why an input or a request was rejected.
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

    /// The store already holds data that this request would duplicate.
    #[error("the clinvar table already holds {existing} rows; ingesting a second release would leave both in the store. Pass --replace to replace the stored release")]
    ClinVarReleasePresent { existing: usize },

    /// A request that cannot be honoured as written.
    #[error("{0}")]
    InvalidArgument(String),
}
