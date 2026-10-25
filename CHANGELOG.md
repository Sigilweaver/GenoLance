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
- INFO and FORMAT values are stored as written in the file (the `info_raw`,
  `format_keys` and `format_raw` columns): percent-encoded text, `;` inside
  FORMAT strings, float precision, trailing fields a record leaves out, and
  `K=.` round-trip through `export` unchanged. Float-typed or negative DP/GQ
  stay in the raw sample column and are counted in the ingest summary; the
  typed `read_depth`/`format_gq` columns stay empty for them.
- Records with POS 0 or ALT `.` are still skipped, and the ingest summary now
  prints how many.

- `genolance export` returns an error instead of a wrong or empty VCF when:
  the sample is not in the store (single and merged export), the store lacks an
  expected column, rows at one position cannot be assigned to source records
  (several multi-allelic records with the same ALT count), or a merged export
  finds several records for one sample at a site (split-biallelic input).
- `export --chrom` (single and merge) matches `chr17` and `17`, and a region
  filter that matches no rows is an error instead of an empty VCF. Merge
  export treats `chr1` and `1` as one contig and writes the first sample's
  spelling.
- `export` keeps records that share a position separate (for example VCFs
  split with `bcftools norm -m-`) instead of merging them into one multi-ALT
  line that kept only the first record's genotype and fields.
- The exported `#CHROM` line names the exported sample (the stored name, which
  `--sample` on ingest may have changed).
- Merged export (repeated `--sample`) re-indexes AD (Number=R) and PL (Number=G) to the merged
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

- The store records a format version, `MAJOR.MINOR`. A store with another
  major is refused with exit code 4. A store from before format versions (all
  earlier releases) cannot be upgraded in place: rebuild it by ingesting the
  VCFs and the ClinVar release into a new store directory. The same major is
  readable at any minor and a newer minor is read-only; a new minor only adds
  nullable columns and metadata keys, and a command that needs a column an
  older store lacks stops with an error naming the table and column. Commands
  that only read a store no longer create the directory when the path is
  wrong; they fail with exit code 2.
- A sample (or the ClinVar release) whose ingest was interrupted is registered
  as incomplete, and every command refuses the store until it is re-ingested
  with `--replace`.
- `genolance export` writes the FORMAT and sample columns exactly as they were
  ingested, and tells records that share a position apart by their source
  record, so several multi-allelic records at one position export correctly.
- Alleles are joined to ClinVar by an upper-cased, trimmed form of (pos, ref,
  alt), so a multi-allelic split that keeps the shared REF (`GCA>GA`) matches
  ClinVar's `GC>G`, and lower-case (soft-masked) bases match. The key is
  parsimonious, not left-aligned: normalise inputs with
  `bcftools norm -f <reference.fa>` so a deletion in a repeat is written the
  way ClinVar writes it. The ingest summary counts the alleles whose key
  differs and the symbolic, breakend and `*` alleles, which cannot be trimmed.
- Every gene of a ClinVar `GENEINFO` record is stored (a list column) and
  searched, ignoring case (`query --gene`, `--genes`, `--acmg`, `pgx`,
  `compound-het`), not only the first. `query --gene` joins on the normalized
  position and fails for a symbol ClinVar does not have.
- The reference assembly of each sample and of the ClinVar release is
  recorded at ingest: from the declared `##contig` lengths (matched against
  built-in GRCh37 and GRCh38 primary-contig lengths) and the `##reference`
  line, or declared with `ingest --assembly GRCh37|GRCh38|T2T-CHM13` (real
  ClinVar files have `##reference=GRCh38` and no `##contig` lines).
  `join`, `screen`, `compound-het` and `pgx` print the stored ClinVar release
  and assembly and refuse a sample on another assembly (exit code 3). If an
  assembly cannot be determined on either side they stop and say to declare it
  with `--assembly`. `--skip-reference-check` lets the command run and prints
  a warning, with what the check would have said, every time. `compare` and
  merged `export` refuse samples on different assemblies; `annotate` checks
  the annotation file the same way and warns when it cannot tell. A sample
  name that is not registered is an error in these commands.
- Float-typed or negative DP/GQ no longer need special handling: the sample
  column is stored as written, and the ingest summary counts the values that
  stay out of the typed columns.

- BCF input is not supported and is now rejected with a clear message (exit
  code 3); the README and command help no longer claim it, and the unused
  `noodles` `bcf` and `csi` features are dropped.
- `export --merge` is removed; repeat `--sample` to merge samples. A single
  `--sample` always exports that sample alone.

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
- `Store::open` no longer creates the store and checks its version; use
  `Store::open_or_create` (ingest) to create one. Added `Store::meta`,
  `set_meta` and `clear_meta`, `IngestReport::normalized_alleles` and
  `unparsed_dp_gq`, and `Filters::skip_reference_check`.
- Store schema: `variants` drops `format_extra` and `format_key_order` for
  `format_keys` and `format_raw`, and gains `source_ordinal` and `norm_*`;
  `clinvar` gains `gene_symbols` and `norm_*`; `samples` gains
  `ingest_complete` and `contigs`; new `store_meta` table.
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
