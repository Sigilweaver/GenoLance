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

- `genolance export` returns an error instead of a wrong or empty VCF when:
  the sample is not in the store (single and `--merge`), the store lacks an
  expected column, rows at one position cannot be assigned to source records
  (several multi-allelic records with the same ALT count), or `--merge` finds
  several records for one sample at a site (split-biallelic input).
- `export --chrom` (single and merge) matches `chr17` and `17`, and a region
  filter that matches no rows is an error instead of an empty VCF. Merge
  export treats `chr1` and `1` as one contig and writes the first sample's
  spelling.
- `export` keeps records that share a position separate (for example VCFs
  split with `bcftools norm -m-`) instead of merging them into one multi-ALT
  line that kept only the first record's genotype and fields.
- The exported `#CHROM` line names the exported sample (the stored name, which
  `--sample` on ingest may have changed).
- `export --merge` re-indexes AD (Number=R) and PL (Number=G) to the merged
  allele order, writes `.` for QUAL when samples disagree and for INFO unless
  all samples agree on it and on the ALT list, and uses the union of the
  samples' non-PASS FILTER values.
- `join`, `screen` and `compound-het` (and `pgx` for FILTER) now:
  - list only calls that carry the ALT allele; no-call genotypes (`./.`,
    missing GT) are no longer treated as carriers and are counted;
  - exclude calls whose FILTER is not PASS and print how many were excluded
    (`--include-filtered` keeps them); the hit lines of `join` show FILTER;
  - print `.` instead of `0` for a missing QUAL or DP, and count calls dropped
    because QUAL/DP is missing under `--min-qual`/`--min-dp`
    (`--keep-missing` keeps them);
  - fail on `--significance-exact` or `--genes` that is empty after splitting
    (it used to match everything);
  - fail when the clinvar table holds two rows for one variant (two releases
    ingested by an older version); re-ingest with `--replace`;
  - fail when a table column is missing or retyped instead of printing zero
    results.
- `chrM`, `M` and `MT` match each other in `join`, `screen`, `compound-het`,
  `pgx`, `annotate` and `query`; `query --chrom` and `annotate --chrom` match
  both `chr17` and `17`.
- `compound-het` flags a hemizygous `GT=1` call (chrX in males) like a
  homozygous one and counts ClinVar hits that have no gene symbol.
- `compare` counts only calls that carry the ALT; `./.`, `0/.` and haploid `0`
  are no longer counted as variants, and no-call rows are reported.
- `pgx --genes` takes a comma-separated list, upper-cases symbols, and errors
  for a gene with no ClinVar rows.
- `annotate` reads the declared Number of each INFO field (`A` and `R` pick the
  value for the ALT; others print all values), errors on a wrong value count or
  an unparsable element, and scans the whole annotation file instead of
  stopping at the first record past `--end`.
- `query --output` accepts only `table` or `json`; anything else is a usage
  error. JSON is one correctly escaped object per line with `null` for a
  missing genotype (it used to break on quotes and print `""` for null).
  `--output arrow` is removed: it printed a summary, not Arrow.
- `query --gene` matches case-insensitively and fails when ClinVar has no rows
  for the symbol (it used to print "no positions matched" and exit 0).
- `join` never ignores the file it is given: with a ClinVar release already in
  the store, `--against FILE` is an error that points at `ingest --replace`;
  without `--against` and without a release it is an error too.
- `annotate` checks its annotation file for truncation and compression like
  `ingest`, fails for an unregistered sample, and fails for an `--info` key the
  annotation header does not declare. `compare` fails for an unregistered
  sample or an unknown `--mode` before reading data. An empty `--significance`
  is an error.

- The command-line tool exits with 1 for internal, I/O or storage-engine
  errors, 2 for usage errors and requests that cannot be honoured as written
  (unknown sample or gene, missing table, empty region, duplicate sample,
  missing store), 3 when input data or stored data is rejected as malformed,
  truncated or inconsistent, and 4 when the store's format version or layout
  is not supported.

### API breaks

- `ingest::run` takes `&IngestOptions` and returns an `IngestReport`
  (skipped-record counts). Added `IngestOptions` and `IngestReport` (both
  `#[non_exhaustive]`).
- Library functions return `genolance_core::Result<T>` with
  `genolance_core::Error`. The enum and each of its struct variants are
  `#[non_exhaustive]`; other crates build the variants through constructor
  functions such as `Error::truncated(..)`. `Error::root`, `Error::context`,
  `Error::exit_code` and the `ResultExt` trait are the helpers. The kinds are
  typed (`MissingColumn`, `ClinVarConflict`, `MissingData`, `RegionEmpty`,
  `UnknownSample`, `UnknownGene`, ...); `InvalidArgument` is only for bad
  argument values.
- `join::Filters`, `screen::Filters` and `compound::Filters` are one
  `filters::Filters` (`#[non_exhaustive]`, new `include_filtered` and
  `keep_missing`). `pgx::run` takes a `&Filters`. `query::run` takes an
  `OutputFormat` instead of a string; `join::run` takes the ClinVar file as an
  `Option`.

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
