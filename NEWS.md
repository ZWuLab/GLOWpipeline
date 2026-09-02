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
