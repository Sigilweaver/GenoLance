use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(
    name = "genolance",
    about = "A fast, columnar multi-sample variant store powered by Lance",
    version
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Ingest one or more VCF files into a GenoLance store
    Ingest {
        /// Path to the GenoLance store (created if it doesn't exist)
        #[arg(short, long)]
        store: String,

        /// VCF files to ingest (plain or BGZF-compressed)
        #[arg(required = true)]
        files: Vec<String>,

        /// Override the sample name (defaults to the sample name in the VCF header).
        /// Only valid with a single file.
        #[arg(long)]
        sample: Option<String>,

        /// Sample column to ingest from a multi-sample VCF
        #[arg(long = "sample-column")]
        sample_column: Option<String>,

        /// Replace a sample (or the ClinVar release) that is already stored
        #[arg(long)]
        replace: bool,

        /// Ingest the file as a ClinVar release even without an INFO/CLNSIG header
        #[arg(long)]
        clinvar: bool,

        /// Reference assembly of the file (GRCh37, GRCh38, T2T-CHM13), for
        /// files whose header does not make it recognisable (no ##contig
        /// lengths and no ##reference naming it)
        #[arg(long)]
        assembly: Option<String>,
    },

    /// Query variants from a GenoLance store
    Query {
        /// Path to the GenoLance store
        #[arg(short, long)]
        store: String,

        /// Filter by gene symbol
        #[arg(long)]
        gene: Option<String>,

        /// Filter by chromosome (e.g. chr1, chrX)
        #[arg(long)]
        chrom: Option<String>,

        /// Filter by start position (1-based, inclusive)
        #[arg(long)]
        start: Option<u64>,

        /// Filter by end position (1-based, inclusive)
        #[arg(long)]
        end: Option<u64>,

        /// Output format: table | json
        #[arg(long, default_value = "table", value_parser = clap::value_parser!(genolance_variants::query::OutputFormat))]
        output: genolance_variants::query::OutputFormat,
    },

    /// Annotate variants by joining against a reference VCF (e.g. ClinVar)
    Join {
        /// Path to the GenoLance store
        #[arg(short, long)]
        store: String,

        /// ClinVar VCF to load when the store has no ClinVar table yet. Not
        /// allowed when the store already holds a release (use
        /// `ingest --replace` to load a new one).
        #[arg(long)]
        against: Option<String>,

        /// Filter by clinical_significance as a case-insensitive substring
        /// (e.g. "pathogenic" also matches "Likely_pathogenic" and "Conflicting...").
        #[arg(long, conflicts_with = "significance_exact")]
        significance: Option<String>,

        /// Filter by clinical_significance as an EXACT, case-insensitive match.
        /// Comma-separated; e.g. "Pathogenic,Likely_pathogenic,Pathogenic/Likely_pathogenic".
        #[arg(long = "significance-exact")]
        significance_exact: Option<String>,

        /// Drop calls with QUAL below this threshold
        #[arg(long = "min-qual")]
        min_qual: Option<f32>,

        /// Drop calls with FORMAT/DP below this threshold
        #[arg(long = "min-dp")]
        min_dp: Option<u32>,

        /// Restrict to the ACMG SF v3 secondary-findings gene list (~80 genes)
        #[arg(long)]
        acmg: bool,

        /// Restrict to a comma-separated gene list
        #[arg(long = "genes")]
        gene_filter: Option<String>,

        /// Include calls whose FILTER is not PASS (excluded and counted by default)
        #[arg(long = "include-filtered")]
        include_filtered: bool,

        /// Keep calls with missing QUAL/DP when --min-qual/--min-dp is set
        /// (excluded and counted by default)
        #[arg(long = "keep-missing")]
        keep_missing: bool,

        /// Skip the check that the samples and the ClinVar release use the
        /// same assembly (compares contig lengths in the headers)
        #[arg(long = "skip-reference-check")]
        skip_reference_check: bool,
    },

    /// Combined carrier-screen + ClinVar pathogenicity filter across N samples.
    ///
    /// Surfaces sites where every listed sample carries at least one ALT
    /// AND the site matches a ClinVar clinical_significance substring
    /// (default "pathogenic", case-insensitive).
    Screen {
        /// Path to the GenoLance store
        #[arg(short, long)]
        store: String,

        /// Sample names to require as carriers (>=2)
        #[arg(required = true)]
        samples: Vec<String>,

        /// ClinVar significance substring to match (default: pathogenic)
        #[arg(long, conflicts_with = "significance_exact")]
        significance: Option<String>,

        /// Exact clinical_significance match (comma-separated). Use this
        /// to exclude "Conflicting_classifications_of_pathogenicity" etc.
        #[arg(long = "significance-exact")]
        significance_exact: Option<String>,

        /// Drop calls with QUAL below this threshold (default: off)
        #[arg(long = "min-qual")]
        min_qual: Option<f32>,

        /// Drop calls with FORMAT/DP below this threshold
        #[arg(long = "min-dp")]
        min_dp: Option<u32>,

        /// Restrict to the ACMG SF v3 secondary-findings gene list
        #[arg(long)]
        acmg: bool,

        /// Restrict to a comma-separated gene list
        #[arg(long = "genes")]
        gene_filter: Option<String>,

        /// Include calls whose FILTER is not PASS (excluded and counted by default)
        #[arg(long = "include-filtered")]
        include_filtered: bool,

        /// Keep calls with missing QUAL/DP when --min-qual/--min-dp is set
        /// (excluded and counted by default)
        #[arg(long = "keep-missing")]
        keep_missing: bool,

        /// Skip the check that the samples and the ClinVar release use the
        /// same assembly (compares contig lengths in the headers)
        #[arg(long = "skip-reference-check")]
        skip_reference_check: bool,
    },

    /// Flag genes where a single sample carries two or more P/LP heterozygous
    /// variants (possible compound het - phase unknown with short-read WGS)
    /// or at least one homozygous P/LP variant.
    CompoundHet {
        /// Path to the GenoLance store
        #[arg(short, long)]
        store: String,

        /// Sample name
        #[arg(long)]
        sample: String,

        /// Significance substring (default: pathogenic)
        #[arg(long, conflicts_with = "significance_exact")]
        significance: Option<String>,

        /// Exact clinical_significance match (comma-separated)
        #[arg(long = "significance-exact")]
        significance_exact: Option<String>,

        /// Drop calls with QUAL below this threshold
        #[arg(long = "min-qual")]
        min_qual: Option<f32>,

        /// Drop calls with FORMAT/DP below this threshold
        #[arg(long = "min-dp")]
        min_dp: Option<u32>,

        /// Restrict to the ACMG SF v3 secondary-findings gene list
        #[arg(long)]
        acmg: bool,

        /// Restrict to a comma-separated gene list
        #[arg(long = "genes")]
        gene_filter: Option<String>,

        /// Include calls whose FILTER is not PASS (excluded and counted by default)
        #[arg(long = "include-filtered")]
        include_filtered: bool,

        /// Keep calls with missing QUAL/DP when --min-qual/--min-dp is set
        /// (excluded and counted by default)
        #[arg(long = "keep-missing")]
        keep_missing: bool,

        /// Skip the check that the samples and the ClinVar release use the
        /// same assembly (compares contig lengths in the headers)
        #[arg(long = "skip-reference-check")]
        skip_reference_check: bool,
    },

    /// Annotate a sample's variants with arbitrary INFO fields from any VCF.
    /// Generic alternative to `join` - match on (chrom, pos, ref, alt) against
    /// any annotation VCF (gnomAD, COSMIC, dbSNP, ...) and emit selected INFO
    /// fields as output columns.
    Annotate {
        /// Path to the GenoLance store
        #[arg(short, long)]
        store: String,

        /// Sample name whose variants drive the join
        #[arg(long)]
        sample: String,

        /// Annotation VCF to match against
        #[arg(long)]
        vcf: String,

        /// Comma-separated INFO field keys to extract (e.g. AF,AF_popmax)
        #[arg(long, value_delimiter = ',')]
        info: Vec<String>,

        /// Restrict to a chromosome
        #[arg(long)]
        chrom: Option<String>,
        /// Restrict to positions >= start (1-based)
        #[arg(long)]
        start: Option<u64>,
        /// Restrict to positions <= end (1-based)
        #[arg(long)]
        end: Option<u64>,
    },

    /// Pharmacogenomic screening: intersect a sample's variants with
    /// ClinVar drug-response annotations in known PGx genes. Screening
    /// only - not a diplotype call.
    Pgx {
        /// Path to the GenoLance store
        #[arg(short, long)]
        store: String,

        /// Sample name to screen
        #[arg(long)]
        sample: String,

        /// Additional gene symbols beyond the default PGx list (comma-separated)
        #[arg(long, value_delimiter = ',')]
        genes: Vec<String>,

        /// Include calls whose FILTER is not PASS (excluded and counted by default)
        #[arg(long = "include-filtered")]
        include_filtered: bool,
        /// Skip the check that the samples and the ClinVar release use the
        /// same assembly (compares contig lengths in the headers)
        #[arg(long = "skip-reference-check")]
        skip_reference_check: bool,
    },

    /// Export variants for a sample back to VCF (reconstructs original header).
    /// Repeat `--sample` to emit a merged multi-sample VCF.
    Export {
        /// Path to the GenoLance store
        #[arg(short, long)]
        store: String,

        /// Sample name to export. Repeat for a merged multi-sample VCF.
        #[arg(long)]
        sample: Vec<String>,

        /// Optional chromosome filter (e.g. chr1)
        #[arg(long)]
        chrom: Option<String>,

        /// Optional start position (1-based, inclusive)
        #[arg(long)]
        start: Option<u64>,

        /// Optional end position (1-based, inclusive)
        #[arg(long)]
        end: Option<u64>,

        /// Output path (default: stdout)
        #[arg(short, long)]
        output: Option<String>,
    },

    /// Build scalar indices on (chrom, pos) and other hot columns so
    /// region / sample / gene lookups stop full-scanning. Safe to re-run;
    /// rebuilds any existing indices on the same columns.
    Index {
        /// Path to the GenoLance store
        #[arg(short, long)]
        store: String,
    },

    /// Compare variants across two or more samples
    Compare {
        /// Path to the GenoLance store
        #[arg(short, long)]
        store: String,

        /// Sample names to compare
        #[arg(required = true)]
        samples: Vec<String>,

        /// Comparison mode: concordance | carrier-screen | private
        #[arg(long, default_value = "concordance")]
        mode: String,
    },
}

/// Exit codes: 0 success; 1 internal, I/O or storage error; 2 usage error or a
/// request that cannot be honoured as written; 3 input data or store contents
/// rejected; 4 store format version or layout not supported.
#[tokio::main]
async fn main() {
    let cli = Cli::parse();
    if let Err(e) = run(cli).await {
        eprintln!("Error: {e}");
        std::process::exit(e.exit_code());
    }
}

// `Filters` is non-exhaustive, so it is built field by field from its default.
#[allow(clippy::field_reassign_with_default)]
async fn run(cli: Cli) -> genolance_core::Result<()> {
    match cli.command {
        Commands::Ingest {
            store,
            files,
            sample,
            sample_column,
            replace,
            clinvar,
            assembly,
        } => {
            let mut opts = genolance_variants::ingest::IngestOptions::new()
                .with_replace(replace)
                .with_clinvar(clinvar);
            if let Some(s) = sample {
                opts = opts.with_sample_override(s);
            }
            if let Some(c) = sample_column {
                opts = opts.with_sample_column(c);
            }
            if let Some(a) = assembly {
                opts = opts.with_assembly(a);
            }
            genolance_variants::ingest::run(&store, &files, &opts).await?;
        }
        Commands::Query {
            store,
            gene,
            chrom,
            start,
            end,
            output,
        } => {
            genolance_variants::query::run(
                &store,
                gene.as_deref(),
                chrom.as_deref(),
                start,
                end,
                output,
            )
            .await?;
        }
        Commands::Join {
            store,
            against,
            significance,
            significance_exact,
            min_qual,
            min_dp,
            acmg,
            gene_filter,
            include_filtered,
            keep_missing,
            skip_reference_check,
        } => {
            let mut f = genolance_variants::filters::Filters::default();
            f.significance_substring = significance;
            f.significance_exact = significance_exact
                .as_deref()
                .map(genolance_core::gene_lists::parse_significance_list);
            f.min_qual = min_qual;
            f.min_dp = min_dp;
            f.acmg_only = acmg;
            f.gene_filter = gene_filter
                .as_deref()
                .map(genolance_core::gene_lists::parse_gene_list);
            f.include_filtered = include_filtered;
            f.keep_missing = keep_missing;
            f.skip_reference_check = skip_reference_check;
            genolance_variants::join::run(&store, against.as_deref(), &f).await?;
        }
        Commands::Screen {
            store,
            samples,
            significance,
            significance_exact,
            min_qual,
            min_dp,
            acmg,
            gene_filter,
            include_filtered,
            keep_missing,
            skip_reference_check,
        } => {
            let mut f = genolance_variants::filters::Filters::default();
            f.significance_substring = significance;
            f.significance_exact = significance_exact
                .as_deref()
                .map(genolance_core::gene_lists::parse_significance_list);
            f.min_qual = min_qual;
            f.min_dp = min_dp;
            f.acmg_only = acmg;
            f.gene_filter = gene_filter
                .as_deref()
                .map(genolance_core::gene_lists::parse_gene_list);
            f.include_filtered = include_filtered;
            f.keep_missing = keep_missing;
            f.skip_reference_check = skip_reference_check;
            genolance_variants::screen::run(&store, &samples, &f).await?;
        }
        Commands::CompoundHet {
            store,
            sample,
            significance,
            significance_exact,
            min_qual,
            min_dp,
            acmg,
            gene_filter,
            include_filtered,
            keep_missing,
            skip_reference_check,
        } => {
            let mut f = genolance_variants::filters::Filters::default();
            f.significance_substring = significance;
            f.significance_exact = significance_exact
                .as_deref()
                .map(genolance_core::gene_lists::parse_significance_list);
            f.min_qual = min_qual;
            f.min_dp = min_dp;
            f.acmg_only = acmg;
            f.gene_filter = gene_filter
                .as_deref()
                .map(genolance_core::gene_lists::parse_gene_list);
            f.include_filtered = include_filtered;
            f.keep_missing = keep_missing;
            f.skip_reference_check = skip_reference_check;
            genolance_variants::compound::run(&store, &sample, &f).await?;
        }
        Commands::Export {
            store,
            sample,
            chrom,
            start,
            end,
            output,
        } => {
            if sample.len() > 1 {
                genolance_variants::export::run_merge(
                    &store,
                    &sample,
                    chrom.as_deref(),
                    start,
                    end,
                    output.as_deref(),
                )
                .await?;
            } else {
                let s = sample.into_iter().next().ok_or_else(|| {
                    genolance_core::Error::InvalidArgument("--sample is required".into())
                })?;
                genolance_variants::export::run(
                    &store,
                    &s,
                    chrom.as_deref(),
                    start,
                    end,
                    output.as_deref(),
                )
                .await?;
            }
        }
        Commands::Annotate {
            store,
            sample,
            vcf,
            info,
            chrom,
            start,
            end,
        } => {
            genolance_variants::annotate::run(
                &store,
                &sample,
                &vcf,
                &info,
                chrom.as_deref(),
                start,
                end,
            )
            .await?;
        }
        Commands::Pgx {
            store,
            sample,
            genes,
            include_filtered,
            skip_reference_check,
        } => {
            let mut f = genolance_variants::filters::Filters::default();
            f.include_filtered = include_filtered;
            f.skip_reference_check = skip_reference_check;
            genolance_variants::pgx::run(&store, &sample, &genes, &f).await?;
        }
        Commands::Index { store } => {
            genolance_variants::index::run(&store).await?;
        }
        Commands::Compare {
            store,
            samples,
            mode,
        } => {
            genolance_variants::compare::run(&store, &samples, &mode).await?;
        }
    }

    Ok(())
}
