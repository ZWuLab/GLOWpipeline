# GLOWpipeline 0.3.0

Requires R 4.0.0 or later (it was 3.5.0), GLOWr 0.3.5 or later, and the CRAN package `yaml`,
which is new among the imports. The
scan functions and their results are unchanged. This release adds the functions that the
workflow templates of GLOWanalyses 0.4.0 use to read one YAML configuration, to place every
output folder, and to keep a record of each unit of work. A user of GLOWanalyses meets them
through its `glow` command and does not call them.

- **One YAML configuration, checked against a registry of settings.** `read_glow_config()`
  reads a configuration file and the files its `extends` lines name, down to a base file, and
  returns every setting converted to its registered type. A setting that the registry does not
  know is refused with the nearest name, a value of the wrong type with the type that is
  expected, and a value outside a fixed list with the list. Every scalar is read as the text
  the user wrote and converted by the type of its setting, so YAML's own guesses (a `no` read
  as false, a version number read as a decimal) play no part. `write_base_config()` writes the
  base file from the registry, each setting with its default and its description, and
  `write_resolved_config()` writes what a run used.
- **The output folders of a study.** `glow_layout()` gives, from the three names
  `output_root`, `dataset` and `analysis`, the folder of the data of a dataset and of each step
  of an analysis. `glow_absolute_paths()` makes the paths of a configuration absolute.
- **Completion records.** A unit of work (a part of a step, or one chromosome or chunk of it)
  writes a record when it has finished: the settings, the input files and the package versions
  it depended on, and the files it wrote. `glow_unit_begin()` locks a unit and decides whether
  it is complete, has to run, or would replace results made with other dependencies, which it
  refuses without `force`. `glow_unit_complete()` writes the record, `glow_unit_state()` and
  `glow_unit_record()` read it, and `glow_dependency_differences()` names what differs. A
  record can name the process that produced the unit (`producer`, for example the job of a
  scheduler), so that a caller that started a job can tell the record that this job wrote.
  `glow_read_parts()` reads the table that says what each part of a workflow depends on, and
  `glow_setting_text()`, `glow_file_fingerprint()` and `glow_run_versions()` give the three
  kinds of dependency as text. `glow_seed_unit()` supplies shipped files in place of a run, with
  a record of the settings they were built with.
- **A lock that survives a run that died.** The lock of a unit records its owner (host, process
  and SLURM job). It is taken over only when the owner is known to be gone, and two runs cannot
  hold it together, also when one of them acts late. `glow_lock()` gives the same lock under a
  name, which the SLURM mode of the `glow` command holds while it submits.
- **The per-family QQ page draws at most one point per grid cell.**
  `summarize_qq_and_lambdas()` drew one vector point per row for every test of its per-family
  page. R's PDF device holds the drawing of a page in memory when it closes the page, about
  0.24 GB per million points, so the page of a whole-genome window scan (2.67 million windows,
  16 tests, 43 million points) took about 10 GB beside the table. The aggregation task of
  such a run on Terra peaked at 14.8 GB on a 16 GB machine (reported by Alberto Brusati,
  Landers lab, UMass Chan Medical School). With `thin = TRUE` (the default of the new
  argument) the points of each test are mapped to the grid that `GLOWr::plot_qq()` uses, half
  the marker's diameter to a cell, and one point per cell is drawn, the highest. The omnibus
  page takes the same argument. The lambdas, the expected quantiles and the axes come from
  the complete input, so the numbers and the RDS file are unchanged, and at the marker's
  scale the picture is the full set's. On a synthetic table of that size each of the two
  per-family pages went from 19.6 MB to 0.13 MB. `thin = FALSE` restores the previous drawing.

# GLOWpipeline 0.2.3

Requires GLOWr 0.3.5.

- **The variant-preprocessing settings have one source.** `glow_scan_chunk()` defaults
  `ld_threshold`, `mac_threshold` and `collapse_method` to `GLOWr::glow_preprocessing_defaults()`
  (0.9 on the absolute correlation, 10, "mean"). `prepare_scan_run()` and `run_scan_unit()`
  resolve a configuration's fields the same way. A missing field takes that default. An invalid
  value (an LD threshold outside (0, 1], a negative or fractional MAC threshold, an unknown
  collapse method) stops the run with a message naming the field before any output is written. Up to 0.2.2 a configuration without the fields passed `NULL` down to the engine,
  which failed inside `GLOWr::filter_variants_ld()` after the null model was fitted. The
  Stage-1 snapshot carries the resolved values. Results are unchanged for a configuration that
  set the fields. The GLOWanalyses base configuration moves from 0.95 and 11 to these values in
  GLOWanalyses 0.3.4 (its NEWS).
- The GLOWr floor rises to 0.3.5, the first version with the defaults accessor.

# GLOWpipeline 0.2.2

Requires GLOWr 0.3.4.

- **`debug_errors` threaded to the engine.** `glow_scan_chunk()` and `run_scan_unit()` take
  `debug_errors` (default FALSE) and pass it to `GLOWr::glow_region_test()`, whose error
  rows then carry the failing call and the stack of function names in their `message`.
  Results unchanged. The GLOWr floor rises to 0.3.4, the first version that accepts the
  argument.
- The optimal-weight p-values of the set-based tests change with GLOWr 0.3.4 (its NEWS):
  the `_N` and `_sparse` schemes of `GLOW_Fisher` and the omnibus columns. The equal-weight
  tests, Burden and the single-variant scan are unchanged.

# GLOWpipeline 0.2.1

Requires GLOWr 0.3.0.

- **The GLOW methods paper is published.** `citation("GLOWpipeline")` and the README cite
  Zhang, Liu, Landers and Wu (2026), The Annals of Applied Statistics 20(3), 2074-2097,
  doi:10.1214/26-AOAS2223, in place of "in revision". Documentation only.
- **The package check passes without GLOWr's optional promoter-mask packages.** One test
  of `build_scan_regions("noncoding")` built the eight default masks and so needed
  `GenomicFeatures` and `TxDb.Hsapiens.UCSC.hg38.knownGene`, which a plain `R CMD check`
  does not install; the public check of 0.2.0 failed on that test. The test now asserts
  the eight masks when the two packages are present and the named error when they are
  absent, and the two packages are listed as suggested, as GLOWr lists them, so that a
  check environment installs them and the promoter masks are exercised there too. No code
  change.

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
