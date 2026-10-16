# GLOWpipeline 0.2.0

Requires GLOWr 0.3.0.

- **The fourth region type, `noncoding`.** `build_scan_regions("noncoding", chr, config,
  gds_path = ...)` returns one regions table per gene-centric non-coding mask of
  `config$noncoding_categories` (default the eight masks), built once from the chromosome's
  annotations by `GLOWr::define_regions_noncoding()`, with the matching
  `GLOWr::noncoding_filter()` map. The new trailing `gds_path` argument is required for this
  type and ignored by the others. `run_scan_unit(region_type = "noncoding")` runs each
  mask's table as it runs a coding category (the label check, `max_regions` and the
  scheduler apply per table) and tags the rows with `category`. It writes the per-chromosome
  table gene-first, ordered by the gene anchor of the new `noncoding_gene_anchor()`, so that a
  gene's masks sit together although their candidate spans differ. Beside the table it writes
  an errors CSV (zero rows when no cell errored) and a per-mask summary CSV (the cells
  attempted, ok, skipped and errored) on every completed unit. `aggregate_coding_results()` serves the type through
  `result_prefix`, `type_label`, `gene_anchor` and `family_alpha` (the across-family cutoff computed as
  `family_alpha / n_all` over the tested cells, in place of the fixed `gw_alpha`). It
  returns a `thresholds` table (every cutoff with its denominator and the published
  constants of `reference_cutoffs` beside it) and a `completeness_table` with one row per
  requested chromosome (its status, the cells, the tallies of the summary sidecar).
  The coding defaults are unchanged.
- **Size-aware scheduling of a unit's regions at `n_cores > 1`.** `run_scan_unit()` no longer
  splits a chromosome's regions into blocks of 25 consecutive regions, which left most workers idle
  while one of them held several neighboring heavy genes. It now reads every region's GDS record
  count from the unit's variant index (no genotype read), takes the squared count as the region's
  cost, runs every region whose cost reaches the unit's total cost divided by 8 times `n_cores`
  in a block of its own, packs the remaining regions in decreasing cost into blocks of at most
  25 regions whose summed cost stays within that threshold, and hands the blocks to the forked
  workers heaviest first.
  The merged output is restored to the regions table's order by label, so the written table, the
  sidecars and the error table are the single-core ones. The gene, window and coding types are all
  scheduled this way, and for the coding type each category's pass over the gene table is. The log
  carries one scheduler line per regions table (regions, blocks, how many run alone, the heaviest
  region, the budget) and one line per finished block (its region count, heaviest record count,
  tally and wall time). On the 1000 Genomes chromosome 22 the table is identical to the 0.1.2 one
  at four and at twenty cores. At twenty workers on a node allocated to the job alone, where the
  old schedule's heaviest block set the wall time, the chromosome took 55 and 56 percent of the
  0.1.2 time in two repeats; at four workers, where the old schedule's heaviest block was already
  below the per-core share of the work, the two schedules were within about 1 percent of each
  other.
- **A repeated region label stops the run.** `run_scan_unit()` stops, at every `n_cores`, when a
  unit's regions table carries a repeated label, naming the duplicates. The sidecars are keyed by
  label and the error table carries only labels, so a duplicate already lost information silently.
  The bundled gene table has no repeated symbol within a chromosome. A user-supplied gene list
  through `config$resolved_gene_list` is what the check is for.

# GLOWpipeline 0.1.2

Requires GLOWr 0.2.3.

- **Large-gene segmentation.** For `region_type = "gene"`, `build_scan_regions()` counts the
  GDS records in each gene's span from the unit's variant index and tiles every gene above
  `config$gene_max_records` (10,000 in the GLOWanalyses defaults) into non-overlapping
  segments of `config$gene_segment_bp` (default 10,000 bp) from the gene start. Each segment is
  its own row, labeled `<gene>_seg<k>`, tested and reported like a gene, with no gene-level
  combination. This is how the GLOW methodology paper handled genes above 1,000 variants. The
  count is taken before any filter, so the split depends only on the file. Without a threshold
  or an index nothing is segmented. `run_scan_unit()` logs how many genes were tiled.
- **The variant index, built once per unit.** `run_scan_unit()` builds GLOWr's variant index
  for its GDS and passes it to `build_scan_regions()` and `glow_scan_chunk()`, which also
  accepts it as `variant_index` (and builds one itself otherwise), so no region re-reads the
  chromosome's variant table. Results are unchanged.
- **`n_cores`** in `run_scan_unit()` (argument, else `config$n_cores`, default 1): the unit's
  regions run in blocks of 25 on forked workers, each with its own GDS handle, and the block
  outputs are merged in order, so the written table equals the single-core one. For one task
  per chromosome on a cloud platform. Not available on Windows.
- `prepare_scan_run()` counts the variants whose QC value is missing (`NA` in the
  `qc_label` node, which is how `SeqArray::seqVCF2GDS()` stores a VCF `FILTER` of
  `"."`) before anything is fitted. GLOWr 0.2.1 excludes such variants instead of
  stopping the region, so the run now warns with the per-chromosome counts, and
  stops when a chromosome would contribute no variant at all.

# GLOWpipeline 0.1.1

- `parse_chr_log()` now reads the per-chromosome summary line written by
  `run_scan_unit()` (`"Done gene unit chrN. n=..., n_ok=..., n_skip_empty=...,
  n_staar_skip_single=..., n_error=..., <t>s"`), in addition to the legacy
  `"Done. genes_attempted=..."` line it already understood. The 2026-06 scan
  unification replaced the legacy line, after which the parser filled every
  count with `NA` on current logs, and a downstream summation over those counts
  failed. Counts map as `n` to `attempted` (and to `gene_list_n` when no
  gene-list line is present), `n_ok + n_staar_skip_single` to `results`,
  `n_skip_empty` to `skipped`, `n_error` to `errors`, and the trailing seconds
  to `runtime_sec`. `results` is the number of result rows the unit wrote: a
  gene whose STAAR sub-step was skipped still has its GLOW result and its row,
  so `attempted == results + skipped + errors` holds as it did for the legacy
  line. Legacy logs parse exactly as before.

# GLOWpipeline 0.1.0

Initial public release.

- Genome chunking for compute-cluster arrays: `define_genome_chunks()` splits a
  chromosome set into contiguous sliding-window chunks (one chunk per array
  task) under a no-duplicate / no-gap boundary rule.
- Single-chunk scan engine: `glow_scan_chunk()` runs the GLOWr single-region
  test over a chunk's windows under a `glow_output_policy()`, returning a flat
  per-window results table.
- Scan-unit runner and Stage-1 setup: `run_scan_unit()` builds a unit's regions
  and filters (`build_scan_regions()`) and persists results; `prepare_scan_run()`
  resolves the shared spine of a multi-region run once for reuse across tasks.
- Genome-wide aggregation: `aggregate_scan_results()` concatenates per-chunk
  window results, flags significant windows, and merges overlapping windows into
  loci; `aggregate_coding_results()` is the gene-centric coding-scan analog.
- Diagnostics and summaries: per-test genomic-control lambda and QQ plots
  (`summarize_qq_and_lambdas()`), scan completeness checking
  (`check_scan_completeness()`), run-log parsing (`parse_chr_log()`,
  `parse_chunk_log()`), and per-error variant drill-downs
  (`write_error_drilldowns()`).
- Built on the GLOWr methods layer (and, transitively, the GFisher backend) for
  the per-region variant-set statistics.
