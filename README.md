# GLOWpipeline: Multi-Region Orchestration for GLOW Variant-Set Scans

## Overview

**GLOWpipeline** is the orchestration-and-summary layer of the GLOW software
family. It scales the single-region GLOW variant-set test
([GLOWr](https://github.com/ZWuLab/GLOWr)) up to a genome-wide scan: it splits a
genome into chunks for a compute-cluster array, runs the GLOWr single-region
test over each chunk's regions under an output policy, and aggregates the
per-chunk results into a genome-wide table with significance thresholding and
overlapping-window locus merging.

GLOWpipeline is the STAARpipeline analog in the
GLOWpipeline -> GLOWr -> GFisher stack. It is deliberately cohort- and
cluster-agnostic (resolved objects in, tables out), so a workflow script binds
it to a specific dataset and scheduler.

## Key Features

- **Genome chunking for cluster arrays**: `define_genome_chunks()` splits a
  chromosome set into contiguous sliding-window chunks, one chunk per array task,
  with a no-duplicate / no-gap boundary rule.
- **Single-chunk scan engine**: `glow_scan_chunk()` runs the GLOWr single-region
  test (`GLOWr::glow_region_test()`) over a chunk's windows under a
  `glow_output_policy()`, returning a flat per-window results table.
- **Scan-unit runner**: `run_scan_unit()` builds a unit's regions and filters
  (`build_scan_regions()`), runs the per-window scan, and persists the results.
- **Run preparation**: `prepare_scan_run()` resolves the shared setup of a
  multi-region run (null model, weight inputs, output directories) once, so each
  chunk task reuses it.
- **Genome-wide aggregation**: `aggregate_scan_results()` concatenates per-chunk
  window results into a genome-wide table, flags significant windows, and merges
  overlapping significant windows into loci; `aggregate_coding_results()` does
  the gene-centric coding-scan analog.
- **Diagnostics and summaries**: per-test genomic-control lambda and QQ plots
  (`summarize_qq_and_lambdas()`), completeness checking
  (`check_scan_completeness()`), run-log parsing (`parse_chr_log()`,
  `parse_chunk_log()`), and error drill-downs (`write_error_drilldowns()`).
- **Built on GLOWr (and GFisher)**: the per-region statistics come from GLOWr,
  which in turn uses GFisher for the generalized Fisher-type combination
  p-values.

## Installation

### From GitHub (recommended)

```r
# install.packages("remotes")
remotes::install_github("ZWuLab/GLOWpipeline")
```

GLOWpipeline depends on **GLOWr** (`ZWuLab/GLOWr`), which itself depends on
**GFisher** (`ZWuLab/GFisher`). The GitHub install resolves this transitive
chain automatically via the `Remotes:` field in `DESCRIPTION`; no separate step
is needed. For a from-source install you may prefer to install the dependencies
first, GFisher before GLOWr:

```r
remotes::install_github("ZWuLab/GFisher")
remotes::install_github("ZWuLab/GLOWr")
remotes::install_github("ZWuLab/GLOWpipeline")
```

### System Requirements

- R >= 3.5.0
- The Bioconductor package **SeqArray** (a hard dependency, for GDS/aGDS I/O)

Install the Bioconductor dependency with:

```r
# install.packages("BiocManager")
BiocManager::install("SeqArray")
```

### Optional Dependencies

Several capabilities rely on Suggested packages installed on demand:

- **Coding-scan and aGDS paths**: `STAAR`, `gdsfmt`, `SNPRelate` (Bioconductor).
- **Fast on-disk results IO**: `fst`.

## Quick Start

A genome-wide GLOW scan has three phases: chunk the genome, scan each chunk
(typically one cluster-array task per chunk), then aggregate. The snippet below
shows the call shapes for the first and last phases; the per-chunk scan
(`glow_scan_chunk()` / `run_scan_unit()`) needs a GDS file plus a fitted null
model and weight inputs, which a workflow script supplies (see the vignette).

```r
library(GLOWpipeline)

# 1. Chunk the genome into sliding-window array tasks. Each row is one chunk:
#    a chromosome and a contiguous run of window start positions.
chrom_spans <- data.frame(
  chr   = c("21", "22"),
  start = c(5010000L, 10510000L),
  end   = c(48120000L, 50810000L)
)
chunks <- define_genome_chunks(
  chrom_spans,
  window_size = 2000, step_size = 1000,
  chunk_strategy = "fixed_mb", chunk_size_mb = 5
)
head(chunks)

# 2. Per chunk (one cluster-array task), run the scan over its windows.
#    glow_scan_chunk() / run_scan_unit() take a GDS path plus the resolved
#    null model + weight inputs from prepare_scan_run(); see the vignette.

# 3. Aggregate the per-chunk result files into a genome-wide table, flag
#    significant windows, and merge overlapping windows into loci.
# agg <- aggregate_scan_results(chunk_files, primary_test = "GLOW_Omni",
#                               alpha = 0.05)
# agg$scan_results_all   # genome-wide flat table
# agg$scan_sig_windows   # significant windows
```

For the full workflow (resolving the shared setup with `prepare_scan_run()`,
running the scan, and aggregating), see the package vignette.

## Documentation

- **Vignette**: `vignette("glowpipeline-overview", package = "GLOWpipeline")` —
  the orchestration overview (chunk -> scan -> aggregate -> summarize).
- **Function help**: `?define_genome_chunks`, `?glow_scan_chunk`,
  `?run_scan_unit`, `?aggregate_scan_results`, `?prepare_scan_run`, etc.

## Testing

The package includes a `testthat` suite:

```r
devtools::test()
```

## Package Structure

```
GLOWpipeline/
├── R/                  # orchestration + summary functions
├── man/                # function help (Rd)
├── tests/testthat/     # unit tests
└── vignettes/          # the orchestration overview vignette
```

## License

GPL-3

## Citation

If you use GLOWpipeline, please cite the GLOW methodology paper:

> Zhang, H., Liu, M., Landers, J. E., and Wu, Z. Integrated Weighted Association
> Test with Application to Genetic Association Studies. *Annals of Applied
> Statistics* (in revision).

See `inst/CITATION` (or run `citation("GLOWpipeline")`).

## Related software / acknowledgments

GLOWpipeline builds on and interoperates with:

- **FAVOR** — functional annotations underpinning the aGDS annotation paths
  consumed by the scan. Zhou, H., Arapoglou, T., Li, X., et al. (2023). FAVOR:
  functional annotation of variants online resource and annotator for variation
  across the human genome. *Nucleic Acids Research*, 51(D1), D1300-D1311.
  doi:10.1093/nar/gkac966
- **STAAR / STAARpipeline** — the STAARpipeline analog this package mirrors; the
  coding aggregator merges native STAAR omnibi and uses STAARpipeline
  `genes_info`. Li, X., Li, Z., Zhou, H., et al. (2020). Dynamic incorporation of
  multiple in silico functional annotations empowers rare variant association
  analysis of large whole-genome sequencing studies at scale. *Nature Genetics*,
  52, 969-983. doi:10.1038/s41588-020-0676-4. Li, Z., Li, X., Zhou, H., et al.
  (2022). A framework for detecting noncoding rare-variant associations of
  large-scale whole-genome sequencing studies. *Nature Methods*, 19, 1599-1611.
  doi:10.1038/s41592-022-01640-x

## Contact

Maintainer: **Zheyang Wu** (zheyangwu@wpi.edu). Please file bugs and questions at
<https://github.com/ZWuLab/GLOWpipeline/issues>.

## Authors

- **Ming Liu** (mliu5@wpi.edu) — Author
- **Hong Zhang** (consistencyzhang@gmail.com) — Author
- **Zheyang Wu** (zheyangwu@wpi.edu) — Author, Maintainer

> Developed with AI assistance; see [AI-USE.md](AI-USE.md).
