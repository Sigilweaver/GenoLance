/// Build an `Error::InvalidArgument` from a format string.
macro_rules! invalid {
    ($($t:tt)*) => {
        genolance_core::Error::InvalidArgument(format!($($t)*))
    };
}

pub mod annotate;
mod assembly;
mod common;
pub mod compare;
pub mod compound;
pub mod export;
pub mod filters;
pub mod index;
pub mod ingest;
pub mod join;
pub mod pgx;
pub mod query;
pub mod screen;
