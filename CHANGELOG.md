# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

- _No unreleased changes yet._

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
