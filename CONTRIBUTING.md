# Contributing to GenoLance

Thanks for your interest in GenoLance. This is a small, single-maintainer
project that ships [Apache-2.0](LICENSE) Rust tooling for fast,
columnar multi-sample variant analysis on top of [Lance].

Crates in this repo: `genolance-core`, `genolance-variants`, `genolance-cli`.

## Before you open a PR

- Open an issue first if the change is non-trivial (new API surface,
  schema change, new ingest/export semantics, dependency bump beyond a
  patch). For small fixes - typos, docs, minor bug fixes, additional
  tests - go straight to a PR.
- Run `cargo fmt --all` and `cargo clippy --all-targets -- -D warnings` locally.
  CI will run them too.
- Run `cargo test --workspace`.
- Update [CHANGELOG.md](CHANGELOG.md) under `## [Unreleased]` with a
  short bullet describing the user-visible change.
- Keep commits small and prefer [Conventional Commits](https://www.conventionalcommits.org/)
  (`feat:`, `fix:`, `docs:`, `refactor:`, `test:`, `chore:`).
- Code is ASCII only and `#![forbid(unsafe_code)]`.

## Working with VCFs and ClinVar

GenoLance is validated against real DeepVariant WGS data (GRCh37,
~5.8M sites per sample) and several ClinVar builds. When changing
ingest, query, or join semantics, please add a fixture-backed test
that exercises the affected code path. Real-data smoke tests live
under `fixtures/` and `tests/`.

ClinVar VCFs are auto-detected by filename. Chromosome naming
mismatch (`chr1` vs `1`) is normalized automatically, so GRCh37
sample VCFs join correctly against GRCh37 ClinVar builds - if you
touch this, please test both directions.

## Roundtrip guarantee

`genolance export` is verified byte-exact (modulo float QUAL
trailing-zero formatting) on DeepVariant VCFs. If you change anything
that touches the export path or per-row `format_key_order` storage,
the roundtrip test in `tests/` must continue to pass.

## Security

Please report security vulnerabilities privately via GitHub Security
Advisories - see [SECURITY.md](SECURITY.md). Do not open public issues
for vulnerabilities.

## DCO

By submitting a contribution you certify that you have the right
to submit the work under the project license (Apache-2.0) and
agree to the
[Developer Certificate of Origin](https://developercertificate.org/).

## License

By submitting a PR you agree that your contribution is licensed under
the Apache License 2.0, the same terms as the rest of the project.

[Lance]: https://lancedb.github.io/lance/
