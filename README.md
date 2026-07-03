# GenoLance

[![CI](https://github.com/Sigilweaver/GenoLance/actions/workflows/ci.yml/badge.svg)](https://github.com/Sigilweaver/GenoLance/actions/workflows/ci.yml)
[![License: Apache-2.0](https://img.shields.io/badge/license-Apache--2.0-blue.svg)](LICENSE)
[![Rust MSRV](https://img.shields.io/badge/rust-1.87%2B-orange.svg)](https://www.rust-lang.org)
[![Docs](https://img.shields.io/badge/docs-sigilweaver.app-blue.svg)](https://sigilweaver.app/genolance/docs/)

A fast, columnar multi-sample variant store powered by [Lance].

GenoLance ingests sample VCFs (and annotation VCFs like ClinVar) into a
LanceDB directory where every field is an Arrow column. Once ingested,
variants can be queried, joined against ClinVar, and compared across
samples entirely out of a single process - no Spark, no cloud, no
GenomicsDB.

## Status

Validated against real DeepVariant WGS data (GRCh37, ~5.8M sites per
sample). It currently supports:

- `genolance ingest` - parse a VCF/BCF and write one Arrow row per
  (sample, chrom, pos, ref, alt). Multi-allelic sites are split.
  ClinVar VCFs are auto-detected by filename and go into a separate
  `clinvar` table.
- `genolance query` - region queries (`--chrom --start --end`) plus
  gene lookups that resolve gene → positions via ClinVar before
  filtering the variants table.
- `genolance join` - annotate variants against ClinVar with an
  optional `--significance pathogenic`-style substring filter, an
  exact-match alternative (`--significance-exact`), quality gates
  (`--min-qual`, `--min-dp`), and gene-set restriction (`--acmg`,
  `--genes`). Chromosome naming mismatch (`chr1` vs `1`) is normalized
  automatically, so GRCh37 sample VCFs join correctly against the
  GRCh37 ClinVar build.
- `genolance compare` - `concordance`, `carrier-screen`, and
  `private` modes across two or more samples.
- `genolance screen` - combined carrier-screen + ClinVar pathogenicity
  filter: sites where *all* listed samples carry at least one ALT AND
  the site matches a ClinVar significance substring. The compound
  carrier-screening question in a single command. Supports exact
  significance matching (`--significance-exact`), quality gates
  (`--min-qual`, `--min-dp`), and gene-set restriction (`--acmg`,
  `--genes`).
- `genolance compound-het` - single-sample screen that flags genes with
  two or more heterozygous P/LP variants (possible compound het, phase
  unknown) or any homozygous P/LP variant. Same filter flags as `screen`.
- `genolance export` - reconstruct a VCF for a given sample (with
  optional region filter). **Verified byte-exact roundtrip** (modulo
  float QUAL trailing-zero formatting) on DeepVariant VCFs:
  FORMAT key order, all field values, INFO, multi-allelic grouping all
  preserved. FORMAT key order is stored per-row via `format_key_order`
  so the original source order (`GT:GQ:DP:AD:VAF:PL`) is reproduced
  faithfully rather than reordered to a canonical form.

## Building

This project depends on `lance-encoding`, which needs `protoc`. The easiest way to get it is through the bundled [pixi] environment:

```sh
cd GenoLance
pixi install                    # installs bcftools, samtools, libprotobuf
cargo build --release           # picks up .cargo/config.toml → $PROTOC
```

## Example

```sh
STORE=/tmp/genolance-demo
rm -rf $STORE

# Ingest two samples and a ClinVar build into one store.
# Match the ClinVar assembly to your sample VCFs (GRCh37 vs GRCh38).
genolance ingest -s $STORE --sample sampleA path/to/sampleA.vcf.gz
genolance ingest -s $STORE --sample sampleB path/to/sampleB.vcf.gz
genolance ingest -s $STORE path/to/clinvar.vcf.gz

# Region query
genolance query -s $STORE --chrom chr17 --start 43044000 --end 43125000

# Gene query (uses ClinVar position index)
genolance query -s $STORE --gene BRCA1

# ClinVar join, only likely-/pathogenic calls the sample actually carries
genolance join -s $STORE --significance pathogenic

# Stricter: exact significance match + quality gates + ACMG SF v3 gene list
genolance join -s $STORE \
  --significance-exact "Pathogenic,Likely_pathogenic,Pathogenic/Likely_pathogenic" \
  --min-qual 20 --min-dp 10 --acmg

# Compare samples
genolance compare -s $STORE sampleA sampleB --mode concordance
genolance compare -s $STORE sampleA sampleB --mode carrier-screen

# Compound screen: sites where ALL listed samples carry a pathogenic ClinVar variant
genolance screen -s $STORE sampleA sampleB
genolance screen -s $STORE sampleA sampleB --significance "likely_pathogenic"
# Clinical-grade variant: exact P/LP match, QUAL>=20, DP>=10, ACMG genes only
genolance screen -s $STORE sampleA sampleB \
  --significance-exact "Pathogenic,Likely_pathogenic,Pathogenic/Likely_pathogenic" \
  --min-qual 20 --min-dp 10 --acmg

# Compound-het / hom-alt screen for one sample: flags genes with 2+ het P/LP
# variants (phase unknown) or any hom-alt P/LP variant. Useful for recessive disease.
genolance compound-het -s $STORE --sample sampleA \
  --significance-exact "Pathogenic,Likely_pathogenic,Pathogenic/Likely_pathogenic" \
  --min-qual 20 --min-dp 10

# Export back to VCF - verified perfect roundtrip on DeepVariant VCFs
genolance export -s $STORE --sample sampleA -o sampleA.roundtrip.vcf
genolance export -s $STORE --sample sampleA --chrom chr17 --start 43044000 --end 43125000
```

## Schema

The `variants` table:

| column        | arrow type | notes                                                  |
| ------------- | ---------- | ------------------------------------------------------ |
| `sample_name` | Utf8       | per-sample label                                       |
| `chrom`       | Utf8       | as written in the VCF (e.g. `chr1` or `1`)             |
| `pos`         | UInt64     | 1-based                                                |
| `ids`         | Utf8       | VCF ID column (`;`-joined), nullable                   |
| `ref_allele`  | Utf8       |                                                        |
| `alt_allele`  | Utf8       | one row per ALT                                        |
| `alt_index`   | UInt32     | 0-based position in the original ALT list              |
| `alt_count`   | UInt32     | total ALTs at this site                                |
| `quality`     | Float32    | QUAL                                                   |
| `filter`      | Utf8       | `;`-joined FILTER field                                |
| `genotype`    | Utf8       | per-ALT encoded GT (`0/1`, `1/1`, `./.`, ...)            |
| `gt_raw`      | Utf8       | original un-split GT (`1/2`, `0\|1`, ...)                |
| `read_depth`  | UInt32     | FORMAT/DP                                              |
| `format_ad`   | Utf8       | FORMAT/AD, comma-joined ("10,5")                       |
| `format_gq`   | UInt32     | FORMAT/GQ                                              |
| `format_pl`   | Utf8       | FORMAT/PL, comma-joined                                |
| `allele_freq` | Float32    | INFO/AF for this ALT (if present)                      |
| `format_extra`| Utf8       | other FORMAT fields (VAF, MIN_DP, ...) as `K=V;K=V`      |
| `format_key_order`| Utf8   | original FORMAT key order (`GT:GQ:DP:AD:VAF:PL`), used by export |
| `info_raw`    | Utf8       | full INFO column re-serialized; duplicated across ALTs |

The `clinvar` table adds `variation_id`, `gene_symbol`,
`clinical_significance`, `review_status`, and `disease_name` columns.

The `samples` table (registry, one row per ingested sample) stores
`sample_name`, `source_path`, `vcf_header` (raw text), `ingested_at`,
and `reference`. It backs `genolance export`'s header reconstruction.

## Known limitations

- Ingest is single-threaded per file; parallelism is on the CLI
  `files` list. For 100+ samples this will be the first thing to
  optimize.
- `--gene` uses exact ClinVar `GENEINFO` symbol matching.
- Per-ALT split genotype (`genotype` column) is lossy for >2-allele
  sites; use `gt_raw` for the original call.
- Chromosome naming (`chr1` vs `1`) is normalized in join/query/screen,
  so a store with `chr`-prefixed variants joins correctly against a
  ClinVar build that uses bare chromosome names (e.g. GRCh37). No
  rewriting on ingest.
- **`genolance export` is verified-lossless on DeepVariant VCFs**:
  all 10,000-site stress tests produce zero diffs (after float
  QUAL normalization). FORMAT key order is stored per-row and reproduced
  exactly. The only intentional semantic difference: QUAL `52.10` is
  emitted as `52.1` (equivalent float, `bcftools`-transparent).
- Scalar indices are **not** built automatically - run
  `genolance index -s $STORE` once after ingest to create BTree
  indices on `pos` and Bitmap indices on `chrom`, `sample_name`,
  `gene_symbol`, and `clinical_significance`. Safe to re-run after
  further ingests to refresh. Gene queries go from ~5s to sub-second
  on the fixture data once indexed.
- Stores created before the `format_key_order` schema change (April 2026)
  need to be re-ingested; column additions are not auto-migrated.

[Lance]: https://lancedb.github.io/lance/
[pixi]: https://pixi.sh
