# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Behavior changes

- `genolance ingest` now returns an error instead of storing a wrong or
  partial result when:
  - a `.gz`/`.bgz` file lacks the BGZF end-of-file block, or the last line
    has no newline (truncated file); a plain-gzip file (not BGZF) and a file
    that is not gzip at all get their own messages;
  - a record has fewer or more sample columns than the header, or lists a
    FORMAT key twice;
  - a genotype is malformed or names an allele beyond the ALT list, a FORMAT
    or INFO array element does not parse, QUAL does not parse, the ALT list
    mixes `.`/empty alleles with real ones, or INFO/AF does not match its
    declared Number and the ALT count;
  - the sample is already stored (use `--replace`), `--sample` is given with
    more than one file, or the VCF has more than one sample column (use
    `--sample-column`);
  - the ClinVar table already holds a release (use `--replace`), or a
    sites-only VCF has no INFO/CLNSIG header (use `--clinvar`).
- ClinVar is recognised by its INFO/CLNSIG header, not by "clinvar" in the
  filename.
- A failed ingest removes the rows it wrote. A sample is registered only
  after its file was read completely.
- INFO and FORMAT values are stored as written in the file: percent-encoded
  text, `;` inside FORMAT strings, float precision, and `K=.` round-trip
  through `export` unchanged. Float-typed or negative DP/GQ are kept in
  `format_extra` instead of being dropped. `format_extra` values now escape
  `%`, `;` and `=` as `%25`, `%3B`, `%3D`.
- Records with POS 0 or ALT `.` are still skipped, and the ingest summary now
  prints how many.

### API breaks

- `ingest::run` takes `&IngestOptions` instead of `Option<&str>` and returns
  an `IngestReport` (skipped-record counts). Added `IngestOptions`,
  `IngestReport` and `genolance_core::Error` (all `#[non_exhaustive]`);
  rejected input surfaces as `genolance_core::Error` inside the
  `anyhow::Error`.

### Security

- Updated the transitive TLS dependencies to rustls 0.23.45 and
  rustls-webpki 0.103.15 to fix RUSTSEC-2026-0285.
- Updated `quinn-proto` to 0.11.19, `rustls` to 0.23.46, and `rustls-webpki`
  to 0.103.16 in the lockfile, clearing RUSTSEC-2026-0338 through -0342,
  RUSTSEC-2026-0346 through -0349, and RUSTSEC-2026-0382. Also updated
  `tokio` (RUSTSEC-2026-0386), `memmap2` (RUSTSEC-2026-0186),
  `event-listener` (RUSTSEC-2026-0221), and `chacha20` (yanked). No code or
  API changes.

### Testing

- Added the first unit-test suite for `genolance-variants`: 28
  synthetic-fixture tests covering the pure/transform-heavy helpers in
  `join`, `compare`, `screen`, and `annotate` (chromosome normalization,
  SQL predicate/escaping, ClinVar-significance filtering including
  `acmg_only` and `gene_filter`, carrier-genotype filtering, and the
  `jaccard`/`truncate` helpers). CI previously ran `cargo test` over an
  empty test set. (#1, contributed by @Nabejo)

## [0.2.0] - 2026-05-22

First publication-ready cut. Moves GenoLance onto the same release and
metadata conventions as the rest of the Sigilweaver suite.

### Changed

- Relicensed from MIT to **Apache-2.0** (no prior release was ever
  published under MIT; this brings GenoLance in line with the rest of
  the Sigilweaver suite, including the explicit patent grant the
  Apache-2.0 license adds for downstream users).
- Moved crate metadata to `[workspace.package]` so `version`,
  `edition`, `rust-version`, `license`, `repository`, `homepage`,
  `documentation`, `readme`, `keywords`, and `categories` are
  declared once and inherited by all three crates.
- Workspace MSRV pinned at `1.87` to match the rest of the suite.
- Added `unsafe_code = "forbid"` workspace lint.

### Added

- `LICENSE` file (Apache-2.0).
- `CONTRIBUTING.md`.
- `CHANGELOG.md` (this file).
- README badges (CI, license, MSRV, docs).
- GitHub Actions CI workflow: `cargo fmt`, `cargo clippy`, `cargo test`.
- `homepage = "https://sigilweaver.app/genolance/"` and
  `documentation = "https://sigilweaver.app/genolance/docs/"` for
  crates.io and docs.rs discovery.

## [0.1.0] - 2025-XX-XX

- Initial layered store layout: `genolance-core`, `genolance-variants`,
  `genolance-cli`.
- CLI subcommands: `ingest`, `query`, `join`, `compare`, `screen`,
  `compound-het`, `export`.
- Validated against DeepVariant WGS data (GRCh37, ~5.8M sites/sample).
- Verified byte-exact VCF roundtrip via `format_key_order`.
