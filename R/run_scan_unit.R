# ==============================================================================
# Unified Per-Unit Scan Body for GLOW SNV-Set Scans (run_scan_unit)
# ==============================================================================
#
# Generalized into a package function so the gene, sliding-window and coding
# scans share ONE per-unit scan body. Given resolved inputs
# (gds_path + the Stage-1 bundle + region_type + the unit selector + config +
# output paths + an output_policy + an optional staar_context), it:
#   1. build_scan_regions(region_type, chr|chunk, config) -> regions + filter(s)
#   2. loop glow_scan_chunk ONCE (gene/window, single filter) or ONCE PER
#      category (coding, category-outer; tag each result table with `category`
#      then rbind, gene-first column order AND a gene-first (start, category) ROW
#      sort so the on-disk per-chr order is gene-first -- design Q6)
#   3. write the per-unit flat results CSV + policy-gated sidecars + a log line.
# With n_cores > 1 each regions table is scheduled over forked workers by size
# (heavy regions alone, cost-budgeted blocks, heaviest first) and the outputs
# are restored to the table's order by label; labels must be unique.
# region_summary is turned ON for region_type %in% c("gene") (covers gene and
# custom-gene-list scans, which dispatch via region_type = "gene") and OFF
# otherwise, reproducing the gene scan's 7 summary columns while leaving the
# sliding/coding schemas unchanged.
#
# This function does NOT parse commandArgs (the thin per-type entry scripts do)
# and does NOT load the Stage-1 snapshot itself (the entry scripts resolve it).
# It is purely additive and changes no existing behavior.
#
# EXPORTED FUNCTIONS:
#   - run_scan_unit()   Run + write one scan unit (chr or chunk) under a policy
#
# INTERNAL HELPER FUNCTIONS:
#   - .rsu_get()        Read one config field from an env/list.
#   - .rsu_null_or()    a if non-NULL else b.
#   - .rsu_write_table() Write a results table as CSV or fst.
#   - .rsu_check_unique_labels() Stop when a regions table repeats a label.
#   - .rsu_region_records() GDS records per region (the cost proxy's input).
#   - .rsu_form_blocks() The size-aware block former (heavy regions alone,
#                        cost-budgeted blocks, dispatch by block cost).
#   - .rsu_merge_chunk_outputs() Merge per-block engine outputs (n_cores > 1)
#                        and restore the regions table's order by label.
#   - .rsu_stack_summaries() Stack per-category summaries into one row.


#################### EXPORTED MAIN FUNCTIONS ####################

#' Run and Persist One GLOW Scan Unit
#'
#' Runs one scan unit: builds the unit's regions and filter(s) with
#' \code{\link{build_scan_regions}}, runs \code{\link{glow_scan_chunk}} (once for
#' \code{gene}/\code{window}; once per category for \code{coding}), and writes
#' the per-unit flat results table plus the policy-gated sidecars, emitting log
#' lines through \code{logger}. A scan "unit" is one chromosome (gene/coding) or
#' one genome chunk (window). The thin per-type entry scripts resolve the
#' prepared inputs and call this with \code{region_type} set.
#'
#' @details
#' \strong{Compute order (coding).} Categories are looped \emph{outer} — one
#' \code{glow_scan_chunk} call per category over the gene table, with the
#' category's filter as the only difference — so each call is exactly a
#' gene-style call and coding reuses the scan engine unchanged. Each call's result
#' table is tagged with a \code{category} column and the tables are row-bound;
#' the combined table leads with \code{label, chr, category} (gene-first column
#' order) and is sorted \code{(start, category)} so the on-disk ROW order is
#' gene-first (each gene's categories together at its genomic position),
#' matching the gene-first aggregation output.
#'
#' \strong{Region summary.} \code{region_summary} is enabled for
#' \code{region_type \%in\% c("gene")} (gene and custom-gene-list scans both
#' dispatch as \code{"gene"}) and disabled for \code{window}/\code{coding}, so the
#' gene flat table carries the 7 region-summary columns and the other schemas are
#' unchanged.
#'
#' \strong{Writes.} The flat results table is written to \code{results_path} (CSV
#' by default, or fst when \code{config$output_format == "fst"}). Policy-gated
#' sidecars (\code{evidence}, \code{staar_detail}, \code{result_objects}) are
#' written under \code{output_dir} into their own on-demand folders, only when
#' both the policy flag is on and the scan engine
#' produced the artifact (so no orphan/empty folders).
#'
#' @param gds_path Character. Path to the unit's GDS/aGDS file.
#' @param null_model,pi_models,pi_features,reference_medians The resolved
#'   prepared inputs (from \code{\link{prepare_scan_run}}), passed to the scan engine.
#'   \code{reference_medians} is this unit's chromosome's annotation-median cache.
#' @param b_func,b_model The B source (exactly one non-NULL), passed to the
#'   engine.
#' @param region_type Character: \code{"gene"}, \code{"window"},
#'   \code{"coding"} or \code{"noncoding"} (the gene-centric non-coding masks,
#'   see the section below).
#' @param config An environment or named list of resolved config values. Used for
#'   the region/filter build (see \code{\link{build_scan_regions}}) and the engine
#'   knobs \code{ld_threshold}, \code{mac_threshold}, \code{collapse_method},
#'   \code{output_format}.
#' @param results_path Character. Path the flat results table is written to.
#' @param output_dir Character. Directory under which sidecar folders
#'   (\code{evidence/}, \code{staar_detail/}, \code{result_objects/}) are created.
#' @param output_policy A \code{\link{glow_output_policy}}.
#' @param chr Integer or character chromosome (required for
#'   \code{gene}/\code{coding}; for \code{window} it is taken from \code{chunk}).
#'   Default NULL.
#' @param chunk A single-row genome-chunk data.frame (required for
#'   \code{window}). Default NULL.
#' @param unit_id Character or integer used in the sidecar filenames and log
#'   lines to identify the unit (e.g. the chr number or chunk id). Default
#'   derived from \code{chr} / \code{chunk$chunk_id}.
#' @param staar_context A STAAR context or NULL (STAAR off). Default NULL.
#' @param use_spa,z_scale Threaded to the engine (see
#'   \code{\link{glow_scan_chunk}}). Defaults NULL / 1.
#' @param max_regions Integer or NA. If not NA, cap the regions to the first
#'   \code{max_regions} (smoke-test; the \code{--max-genes}/\code{--max-windows}
#'   parallel). Default NA.
#' @param progress_every Integer threaded to the engine. Default 0L.
#' @param logger A \code{function(msg)} for log lines, or NULL (silent). Default
#'   NULL.
#' @param n_cores Integer or NULL (default). Forked workers for this unit. NULL
#'   reads \code{config$n_cores} and falls back to 1. With more than one, the
#'   regions are grouped into blocks by size (see the section "Scheduling over
#'   the workers"), each block runs \code{\link{glow_scan_chunk}} in a forked
#'   worker with its own GDS handle (\code{parallel::mclapply}, dynamic
#'   dispatch), and the block outputs are merged and restored to the regions
#'   table's order by label, so the written table equals the single-core one.
#'   The progress callback is off inside the workers; the parent logs one
#'   scheduler line per regions table and one line per finished block. Forking
#'   is unavailable on Windows (a value above 1 warns and runs serially).
#'   Memory scales with the workers, since each holds its own genotype and LD
#'   matrices.
#' @param verbose Integer threaded to the engine. Default 0.
#'
#' @section Gene-centric non-coding masks:
#' For \code{region_type = "noncoding"} the unit is a chromosome and
#' \code{\link{build_scan_regions}} returns one regions table per mask of
#' \code{config$noncoding_categories} (default the eight masks of
#' \code{GLOWr::noncoding_masks()}), built once from the chromosome's
#' annotations by \code{GLOWr::define_regions_noncoding}. Each mask's table is
#' run like a coding category (the label check, \code{max_regions} and the
#' scheduler apply per table; \code{region_summary} is off), the rows are
#' tagged with \code{category}, and the combined table leads with
#' \code{label, chr, category} and is sorted gene-first by the gene anchor of
#' \code{\link{noncoding_gene_anchor}} (a gene's masks together, because their
#' candidate spans differ), then by mask. A row's \code{start}, \code{end} and
#' \code{midpoint} describe the candidate span, which is wide for an enhancer
#' far from its gene. Two sidecars are written beside the flat table on every
#' completed unit, for the aggregation's completeness table: the error table as
#' \code{<results_path stem>_errors.csv} (zero rows when no cell errored, so a
#' rerun never leaves an earlier run's errors beside a clean table) and the
#' per-mask tallies (cells attempted, tested, skipped for an empty set or a
#' single rare variant, errored) as \code{<results_path stem>_summary.csv}.
#' The log carries one line per mask with the genes, the selected records and
#' the links.
#'
#' @section Large genes:
#' For \code{region_type = "gene"} the unit's variant index is passed to
#' \code{\link{build_scan_regions}}, which tiles every gene with more GDS
#' records than \code{config$gene_max_records} into segments of
#' \code{config$gene_segment_bp} (see that function). The log states how many
#' genes were segmented and into how many segments.
#'
#' @section Scheduling over the workers:
#' With \code{n_cores > 1} the regions of a table are not run in genomic order.
#' Heavy genes sit next to each other on the genome, and a block of consecutive
#' regions that holds several of them keeps one worker busy long after the
#' others have finished. Instead, every region's GDS record count \eqn{n_i} is
#' read from the unit's variant index (two binary searches per region, no
#' genotype read; the gene segmentation's \code{n_records} is reused where
#' present) and its cost is taken as \eqn{c_i = n_i^2}, a heuristic that lets a
#' few heavy regions dominate the sums. With \eqn{C = \sum_i c_i} and
#' \eqn{\tau = C / (8\,\mathrm{n\_cores})}, a region with \eqn{c_i \ge \tau}
#' forms a block of its own; the remaining regions, in decreasing cost, are
#' packed into blocks that close at 25 regions or when the next region would
#' push the summed cost above \eqn{\tau}; and the blocks are dispatched in
#' decreasing summed cost, so the heaviest work starts first and the light
#' blocks fill the tail. When \eqn{C = 0} the regions keep their order in blocks
#' of 25. Each region is computed independently of its block, so the merged
#' rows are the serial ones; the merge puts the flat rows, the error rows and
#' the sidecar entries back in the order of the labels in the regions table.
#' This needs unique labels, which the function enforces at every
#' \code{n_cores}: a repeated label stops the run with the duplicates named,
#' because the sidecars and the error table are keyed by label and a duplicate
#' already loses information silently in the serial path. For the coding type
#' each category's pass over the gene table is scheduled this way, one category
#' after another. The log's scheduler line gives the number of blocks, how many
#' are a heavy region alone, the heaviest region and \eqn{\tau}; each block line
#' gives the block's region count, its heaviest record count, its tally and
#' its wall time, so the reader sees where the unit's time went.
#'
#' @return Invisibly, a named list with \code{results} (the written flat table),
#'   \code{summary} (a 1-row tally; the engine summary for gene/window, or a
#'   per-category-stacked summary for coding), \code{errors} (a data.frame), and
#'   \code{results_path}.
#'
#' @seealso \code{\link{build_scan_regions}}, \code{\link{glow_scan_chunk}},
#'   \code{\link{glow_output_policy}}, \code{\link{prepare_scan_run}}
#'
#' @examples
#' \dontrun{
#' bundle <- prepare_scan_run(env, "runs/myrun/outputs", chroms = 22L)
#' run_scan_unit(
#'   gds_path = "data/chr22.gds",
#'   null_model = bundle$null_model, pi_models = bundle$pi_models,
#'   pi_features = bundle$pi_features, reference_medians = meds22,
#'   b_func = bundle$b_func, b_model = bundle$b_model,
#'   region_type = "gene", chr = 22, config = env,
#'   results_path = "runs/myrun/outputs/results/glow_chr22.csv",
#'   output_dir = "runs/myrun/outputs",
#'   output_policy = glow_output_policy(), logger = message)
#' }
#' @export
run_scan_unit <- function(gds_path,
                          null_model, pi_models, pi_features, reference_medians,
                          b_func = NULL, b_model = NULL,
                          region_type, config,
                          results_path, output_dir, output_policy,
                          chr = NULL, chunk = NULL, unit_id = NULL,
                          staar_context = NULL,
                          use_spa = NULL, z_scale = 1,
                          max_regions = NA_integer_,
                          progress_every = 0L,
                          logger = NULL, n_cores = NULL, verbose = 0) {

  region_type <- match.arg(region_type, c("gene", "window", "coding", "noncoding"))
  log_msg <- if (is.null(logger)) function(...) invisible(NULL) else logger

  # The comparator's annotation columns must be among the extracted ones
  # (pi_features); otherwise every region would error at the comparator step
  # and the unit would end with no comparator value at all.
  if (!is.null(staar_context) && !is.null(staar_context$anno_cols)) {
    miss <- setdiff(staar_context$anno_cols, pi_features)
    if (length(miss) > 0L) {
      stop("staar_context$anno_cols names annotation column(s) that are not among ",
           "pi_features, the columns extracted for every region: ",
           paste(miss, collapse = ", "),
           ". Pass anno_cols as a subset of pi_features.", call. = FALSE)
    }
  }

  # ---- Cores: the argument, else the config field, else 1 ----
  if (is.null(n_cores)) n_cores <- .rsu_null_or(.rsu_get(config, "n_cores"), 1L)
  if (!is.numeric(n_cores) || length(n_cores) != 1L || is.na(n_cores) ||
      n_cores < 1) {
    stop("`n_cores` must be a single integer >= 1.", call. = FALSE)
  }
  n_cores <- as.integer(n_cores)
  if (n_cores > 1L && .Platform$OS.type == "windows") {
    warning("n_cores > 1 needs forked workers, which Windows does not support; ",
            "running with one core.")
    n_cores <- 1L
  }

  # ---- Resolve the unit selector (chr | chunk) ----
  if (region_type == "window") {
    if (is.null(chunk) || !is.data.frame(chunk) || nrow(chunk) != 1L) {
      stop("`region_type = \"window\"` requires a single-row `chunk`.",
           call. = FALSE)
    }
    chr <- as.character(chunk$chr)
    if (is.null(unit_id))
      unit_id <- if (!is.null(chunk$chunk_id)) chunk$chunk_id else chr
  } else {
    if (is.null(chr)) {
      stop("`region_type = \"", region_type, "\"` requires `chr`.",
           call. = FALSE)
    }
    if (is.null(unit_id)) unit_id <- chr
  }

  # ---- region_summary: on for gene scans (dispatched as "gene"), off else ----
  region_summary <- region_type %in% c("gene")

  # ---- Engine knobs (from config) ----
  ld_threshold    <- .rsu_get(config, "ld_threshold")
  mac_threshold   <- .rsu_get(config, "mac_threshold")
  collapse_method <- .rsu_get(config, "collapse_method")
  output_format   <- .rsu_null_or(.rsu_get(config, "output_format"), "csv")

  # ---- The variant index: built once per unit ----
  # Positions, ids and file positions of the unit's GDS, read once. The
  # segmentation counts records with it and every region's lookup uses it.
  variant_index <- GLOWr::build_variant_index(gds_path)

  # ---- Build regions + filter(s) for this unit ----
  built <- build_scan_regions(region_type, chr, config, chunk = chunk,
                              variant_index = variant_index, gds_path = gds_path)
  if (is.data.frame(built$regions)) {
    n_seg_genes <- attr(built$regions, "n_genes_segmented")
    if (!is.null(n_seg_genes)) {
      log_msg(sprintf(paste0("gene segmentation: %d gene(s) above gene_max_records = %s ",
                             "records tiled into %d segment(s) of %s bp"),
                      n_seg_genes, format(.rsu_get(config, "gene_max_records")),
                      attr(built$regions, "n_segments"),
                      format(.rsu_null_or(.rsu_get(config, "gene_segment_bp"), 10000L))))
    }
  } else {                                     # noncoding: one table per mask
    for (k in names(built$regions)) {
      r <- built$regions[[k]]
      log_msg(sprintf(paste0("noncoding mask %-14s: %d gene(s) with linked records, ",
                             "%d record(s) selected, %d (gene, record) link(s)"),
                      k, nrow(r), .rsu_null_or(attr(r, "n_selected_records"), NA_integer_),
                      .rsu_null_or(attr(r, "n_linked_records"), NA_integer_)))
    }
  }
  # Labels are the key of the sidecars, the error table and the parallel merge,
  # so a repeated label stops the run at every n_cores (plan 12-size-aware, D6).
  # max_regions caps a regions table: the unit's table, or each mask's table.
  prepare_regions <- function(regs, what) {
    if (!"label" %in% names(regs)) {           # the engine's default label
      regs$label <- paste0("chr", regs$chr, ":", regs$start, "-", regs$end)
    }
    .rsu_check_unique_labels(regs)
    if (!is.na(max_regions) && nrow(regs) > max_regions) {
      regs <- regs[seq_len(max_regions), , drop = FALSE]
      log_msg(sprintf("max_regions %d: limiting %s to %d regions",
                      max_regions, what, nrow(regs)))
    }
    regs
  }
  if (is.data.frame(built$regions)) {
    regions <- prepare_regions(built$regions, "the unit")
  } else {
    built$regions <- lapply(stats::setNames(names(built$regions), names(built$regions)),
                            function(k) prepare_regions(built$regions[[k]], paste("mask", k)))
    regions <- NULL
  }

  # One engine call over a regions table; the progress callback only when the
  # call runs in this process (a forked worker must not write to the parent's
  # log connection).
  engine_call <- function(regs, filter_spec, in_parent) {
    glow_scan_chunk(
      gds_path, regs, null_model, pi_models, pi_features, reference_medians,
      b_func = b_func, b_model = b_model, filter_spec = filter_spec,
      ld_threshold = ld_threshold, mac_threshold = mac_threshold,
      collapse_method = collapse_method,
      use_spa = use_spa, z_scale = z_scale,
      region_summary = region_summary,
      output_policy = output_policy, staar_context = staar_context,
      progress_every = if (in_parent) progress_every else 0L,
      progress_callback = if (in_parent && !is.null(logger))
        function(m) log_msg(m) else NULL,
      variant_index = variant_index,
      verbose = verbose)
  }

  # The schedule of one regions table over the workers (the roxygen section
  # "Scheduling over the workers"): NULL means the serial path. Computed once
  # per table and logged once, so the coding categories share one schedule.
  schedule_of <- function(regs) {
    if (n_cores == 1L || nrow(regs) <= 1L) return(NULL)
    n_rec  <- .rsu_region_records(regs, variant_index)
    blocks <- .rsu_form_blocks(n_rec^2, n_cores)
    heaviest <- which.max(n_rec)
    log_msg(sprintf(paste0("scheduler (n_cores = %d): %d regions in %d blocks, %d ",
                           "of them a heavy region alone; heaviest %s with %d ",
                           "records; cost = records^2, block budget tau = %.4g"),
                    n_cores, nrow(regs), length(blocks),
                    attr(blocks, "n_singletons"), as.character(regs$label[heaviest]),
                    as.integer(n_rec[heaviest]), attr(blocks, "tau")))
    list(blocks = blocks, n_records = n_rec)
  }

  run_one <- function(regs, filter_spec, sched = NULL) {
    if (is.null(sched)) {
      return(engine_call(regs, filter_spec, in_parent = TRUE))
    }
    # One forked worker per block, the blocks handed out heaviest first as
    # workers free up (mc.preschedule = FALSE). Each region is computed
    # independently of its block, so the merge below only has to put the rows
    # back in the table's order.
    blocks <- sched$blocks
    t_par <- proc.time()[["elapsed"]]
    outs <- parallel::mclapply(
      blocks,
      function(idx) engine_call(regs[idx, , drop = FALSE], filter_spec,
                                in_parent = FALSE),
      mc.cores = n_cores, mc.preschedule = FALSE)
    failed <- vapply(outs, function(o) inherits(o, "try-error"), logical(1))
    if (any(failed)) {
      k <- which(failed)[1L]
      stop(sprintf("run_scan_unit: worker for block %d (%d region(s): %s) failed: %s",
                   k, length(blocks[[k]]),
                   paste(utils::head(as.character(regs$label[blocks[[k]]]), 5L),
                         collapse = ", "),
                   as.character(outs[[k]])), call. = FALSE)
    }
    for (k in seq_along(outs)) {
      sm <- outs[[k]]$summary
      log_msg(sprintf(paste0("block %d/%d: %d region(s), heaviest %d records, ",
                             "n_ok=%d n_skip_empty=%d n_error=%d %.1fs"),
                      k, length(outs), length(blocks[[k]]),
                      as.integer(max(sched$n_records[blocks[[k]]])),
                      sm$n_ok, sm$n_skip_empty, sm$n_error, sm$runtime_secs))
    }
    .rsu_merge_chunk_outputs(outs, wall_secs = proc.time()[["elapsed"]] - t_par,
                             region_summary = region_summary,
                             labels = as.character(regs$label))
  }

  # ---- Run: ONCE (gene/window) or ONCE PER CATEGORY (coding, noncoding) ----
  evidence_acc <- list(); staar_detail_acc <- list(); result_objects_acc <- list()

  if (region_type %in% c("coding", "noncoding")) {
    category_labels <- names(built$filters)
    acc <- list(); err_rows <- list(); summaries <- list()
    if (region_type == "coding") {
      log_msg(sprintf("Starting coding chr%s: %d genes, %d categories",
                      chr, nrow(regions), length(category_labels)))
      sched_shared <- schedule_of(regions)     # the categories share the gene table
    } else {
      log_msg(sprintf("Starting noncoding chr%s: %d masks (%s), %d (gene, mask) cells",
                      chr, length(category_labels), paste(category_labels, collapse = ", "),
                      sum(vapply(built$regions, nrow, integer(1)))))
    }
    for (category in category_labels) {
      regs  <- if (region_type == "coding") regions else built$regions[[category]]
      sched <- if (region_type == "coding") sched_shared else schedule_of(regs)
      out <- run_one(regs, built$filters[[category]], sched)
      if (nrow(out$results) > 0L) out$results$category <- category
      acc[[category]] <- out$results
      summaries[[category]] <- out$summary
      if (isTRUE(output_policy$write_evidence)       && length(out$evidence)       > 0L)
        evidence_acc[[category]]       <- out$evidence
      if (isTRUE(output_policy$write_staar_detail)   && length(out$staar_detail)   > 0L)
        staar_detail_acc[[category]]   <- out$staar_detail
      if (isTRUE(output_policy$write_result_objects) && length(out$result_objects) > 0L)
        result_objects_acc[[category]] <- out$result_objects
      s <- out$summary
      log_msg(sprintf(paste0("category %-20s: %d cells (n_ok=%d n_skip_empty=%d ",
                             "n_staar_skip_single=%d n_error=%d %.1fs)"),
                      category, nrow(out$results), s$n_ok, s$n_skip_empty,
                      s$n_staar_skip_single, s$n_error, s$runtime_secs))
      if (nrow(out$errors) > 0L) {
        out$errors$category <- category
        err_rows[[category]] <- out$errors
      }
    }
    res <- do.call(rbind, acc)
    if (is.null(res) || nrow(res) == 0L) {
      res <- data.frame(category = character(0L))
      log_msg("WARNING: no (gene, category) cell tested on this unit (all empty).")
    } else {
      front <- c("label", "chr", "category")           # gene-first column order
      res <- res[, c(front, setdiff(names(res), front)), drop = FALSE]
      # Gene-first ROW order on disk (design Q6): sort (gene_start, category) so a
      # gene's categories sit together at the gene's genomic position, matching
      # the gene-first aggregation output. The compute loop is still category-outer
      # (the rbind above is category-grouped); this re-sort only changes ROW ORDER
      # (sub-phase 2c -- the one intended change vs 2b's category-outer on-disk
      # order; verified content-identical after a (label, category) key sort).
      if (region_type == "coding" && "start" %in% names(res)) {
        res <- res[order(res$start, res$category), , drop = FALSE]
      } else if (region_type == "noncoding") {
        # A gene's masks have different candidate spans, so the gene-first
        # order uses the gene anchor (the bundled gene table's start for a
        # protein-coding gene, the first candidate position for an ncRNA gene).
        anchor <- noncoding_gene_anchor(res)
        a <- anchor$anchor[match(paste(res$chr, res$label), paste(anchor$chr, anchor$gene))]
        res <- res[order(a, res$label, res$category), , drop = FALSE]
      }
      rownames(res) <- NULL
    }
    errors <- if (length(err_rows) > 0L) do.call(rbind, err_rows) else
      data.frame(category = character(0L), label = character(0L),
                 message = character(0L), stringsAsFactors = FALSE)
    errors <- errors[, c("category", "label", "message"), drop = FALSE]  # stable schema
    rownames(errors) <- NULL
    # One stacked summary row (sum the per-category tallies, max runtime), and
    # the per-category tallies themselves (the noncoding summary sidecar).
    summary <- .rsu_stack_summaries(summaries)
    summary_by_category <- do.call(rbind, lapply(names(summaries), function(k) {
      s <- summaries[[k]]
      data.frame(category = k, n_attempted = as.integer(s$n_windows),
                 n_ok = as.integer(s$n_ok), n_skip_empty = as.integer(s$n_skip_empty),
                 n_staar_skip_single = as.integer(s$n_staar_skip_single),
                 n_error = as.integer(s$n_error), runtime_secs = s$runtime_secs,
                 stringsAsFactors = FALSE)
    }))
  } else {
    out <- run_one(regions, built$filters, schedule_of(regions))
    res <- out$results
    rownames(res) <- NULL
    errors  <- out$errors
    summary <- out$summary
    if (isTRUE(output_policy$write_evidence)       && length(out$evidence)       > 0L)
      evidence_acc       <- out$evidence
    if (isTRUE(output_policy$write_staar_detail)   && length(out$staar_detail)   > 0L)
      staar_detail_acc   <- out$staar_detail
    if (isTRUE(output_policy$write_result_objects) && length(out$result_objects) > 0L)
      result_objects_acc <- out$result_objects
    log_msg(sprintf(paste0("Done %s unit %s. n=%d, n_ok=%d, n_skip_empty=%d, ",
                           "n_staar_skip_single=%d, n_error=%d, %.1fs"),
                    region_type, as.character(unit_id), summary$n_windows,
                    summary$n_ok, summary$n_skip_empty,
                    summary$n_staar_skip_single, summary$n_error,
                    summary$runtime_secs))
  }

  # ---- Write the flat results table ----
  dir.create(dirname(results_path), recursive = TRUE, showWarnings = FALSE)
  .rsu_write_table(res, results_path, output_format)
  log_msg(sprintf("Wrote %s (%d rows)", results_path, nrow(res)))
  # The noncoding sidecars for the aggregation's completeness table: the error
  # rows and the per-mask tallies (cells attempted, tested, skipped, errored).
  # Both are written on every completed unit, the errors table with zero rows
  # when no cell errored, so a rerun never leaves an earlier run's errors file
  # beside a clean table.
  if (region_type == "noncoding") {
    stem  <- sub("\\.(csv|fst)$", "", basename(results_path))
    p_err <- file.path(dirname(results_path), paste0(stem, "_errors.csv"))
    utils::write.csv(errors, p_err, row.names = FALSE)
    log_msg(sprintf("Wrote %s (%d error row(s))", p_err, nrow(errors)))
    p_sum <- file.path(dirname(results_path), paste0(stem, "_summary.csv"))
    utils::write.csv(summary_by_category, p_sum, row.names = FALSE)
    log_msg(sprintf("Wrote %s (%d mask row(s))", p_sum, nrow(summary_by_category)))
  }

  # ---- Policy-gated sidecars (own on-demand folders; no orphans) ----
  if (isTRUE(output_policy$write_evidence) && length(evidence_acc) > 0L) {
    d <- file.path(output_dir, "evidence")
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    p <- file.path(d, sprintf("evidence_%s.rds", as.character(unit_id)))
    saveRDS(evidence_acc, p)
    log_msg(sprintf("Wrote %s", p))
  }
  if (isTRUE(output_policy$write_staar_detail) && length(staar_detail_acc) > 0L) {
    d <- file.path(output_dir, "staar_detail")
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    p <- file.path(d, sprintf("staar_detail_%s.rds", as.character(unit_id)))
    saveRDS(staar_detail_acc, p)
    log_msg(sprintf("Wrote %s", p))
  }
  if (isTRUE(output_policy$write_result_objects) &&
      length(result_objects_acc) > 0L) {
    d <- file.path(output_dir, "result_objects")
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    p <- file.path(d, sprintf("result_objects_%s.rds", as.character(unit_id)))
    saveRDS(result_objects_acc, p)
    log_msg(sprintf("Wrote %s", p))
  }

  if (nrow(errors) > 0L) {
    log_msg(sprintf("%d error(s) on this unit:", nrow(errors)))
    for (k in seq_len(nrow(errors))) {
      lab <- if ("category" %in% names(errors))
        paste0(errors$category[k], "/", errors$label[k]) else errors$label[k]
      log_msg(sprintf("  %s: %s", lab, errors$message[k]))
    }
  }

  invisible(list(results = res, summary = summary, errors = errors,
                 results_path = results_path))
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Read One Config Field from an Environment or List
#' @keywords internal
#' @noRd
.rsu_get <- function(config, name) {
  if (is.environment(config)) {
    get0(name, envir = config, inherits = FALSE)
  } else if (is.list(config)) {
    config[[name]]
  } else {
    stop("`config` must be an environment or a named list.")
  }
}

#' Return a if Non-NULL, Else b
#' @keywords internal
#' @noRd
.rsu_null_or <- function(a, b) if (is.null(a)) b else a

#' Write a Results Table as CSV (default) or fst
#' @keywords internal
#' @noRd
.rsu_write_table <- function(df, path, output_format = "csv") {
  if (identical(output_format, "fst")) {
    if (!requireNamespace("fst", quietly = TRUE)) {
      stop("output_format = 'fst' but the fst package is not installed.")
    }
    fst::write_fst(df, path)
  } else {
    utils::write.csv(df, path, row.names = FALSE)
  }
}

#' Stop When a Regions Table Repeats a Label
#'
#' The sidecars of \code{glow_scan_chunk()} are keyed by label and its error
#' table carries only labels, so a repeated label already loses information
#' silently in the serial path, and the parallel merge restores the table's
#' order by label. The bundled gene table has no repeated symbol within a
#' chromosome, window labels encode their coordinates and segment labels carry
#' their index; a user-supplied gene list is what this check is for.
#'
#' @param regions A regions table with a \code{label} column.
#' @return Invisibly TRUE; stops otherwise, naming the duplicates.
#' @keywords internal
#' @noRd
.rsu_check_unique_labels <- function(regions) {
  lab <- as.character(regions$label)
  dup <- unique(lab[duplicated(lab)])
  if (length(dup) > 0L) {
    stop(sprintf(paste0("run_scan_unit: the unit's regions carry %d repeated ",
                        "label(s): %s%s. Labels must be unique, because the ",
                        "sidecars and the error table are keyed by label."),
                 length(dup), paste(utils::head(dup, 10L), collapse = ", "),
                 if (length(dup) > 10L) ", ..." else ""), call. = FALSE)
  }
  invisible(TRUE)
}

#' GDS Records of Every Region (the Scheduler's Cost Input)
#'
#' Reuses the table's \code{n_records} column where the gene segmentation
#' attached it (or the gene-centric non-coding builder, whose candidate count
#' is the right size for a set whose descriptive span is wide), and otherwise
#' counts the records in each span from the variant index with
#' \code{GLOWr::count_index_records()}: two binary searches per region when
#' the chromosome's positions are sorted, no genotype read. A missing count
#' is taken as 0.
#'
#' @param regions A regions table (\code{chr}, \code{start}, \code{end}).
#' @param variant_index A \code{glow_variant_index} of the unit's GDS.
#' @return A numeric vector of counts, one per row.
#' @keywords internal
#' @noRd
.rsu_region_records <- function(regions, variant_index) {
  if (nrow(regions) == 0L) return(numeric(0))
  n <- if ("n_records" %in% names(regions)) as.numeric(regions$n_records) else
    as.numeric(GLOWr::count_index_records(variant_index, regions$chr,
                                          regions$start, regions$end))
  n[is.na(n)] <- 0
  n
}

#' Form the Size-Aware Blocks of a Regions Table
#'
#' The block former of the scheduler (plan 12-size-aware-scheduling, section
#' 4.3). With \eqn{C = \sum_i c_i > 0} and \eqn{\tau = C / (8 m)} for
#' \code{m = n_cores}: the regions are sorted by decreasing cost (ties by row
#' index), every region with \eqn{c_i \ge \tau} becomes a block of its own, the
#' rest are packed in sorted order into blocks that close at \code{block_size}
#' regions or when the next region would push the summed cost above
#' \eqn{\tau}, and the blocks are ordered by decreasing summed cost (ties by
#' formation order). With \eqn{C = 0} the rows keep their order in blocks of
#' \code{block_size}. Every row lands in exactly one block.
#'
#' @param cost Numeric vector of region costs (non-negative; NA counts as 0).
#' @param n_cores Integer, the workers.
#' @param block_size Integer, the cap on a packed block (25).
#' @return A list of integer row-index vectors, heaviest block first, with the
#'   attributes \code{tau}, \code{n_singletons} and \code{block_cost}.
#' @keywords internal
#' @noRd
.rsu_form_blocks <- function(cost, n_cores, block_size = 25L) {
  n <- length(cost)
  cost <- as.numeric(cost)
  cost[is.na(cost)] <- 0
  if (n == 0L) {
    return(structure(list(), tau = 0, n_singletons = 0L, block_cost = numeric(0)))
  }
  total <- sum(cost)
  if (total <= 0) {
    blocks <- unname(split(seq_len(n), ceiling(seq_len(n) / block_size)))
    return(structure(blocks, tau = 0, n_singletons = 0L,
                     block_cost = numeric(length(blocks))))
  }
  tau <- total / (8 * n_cores)
  ord <- order(-cost, seq_len(n))                 # decreasing cost, ties by row
  heavy <- cost[ord] >= tau
  singles <- as.list(ord[heavy])                  # one heavy region per block
  packed <- list()
  cur <- integer(0); cur_cost <- 0
  for (i in ord[!heavy]) {                        # light regions, sorted order
    if (length(cur) > 0L &&
        (length(cur) >= block_size || cur_cost + cost[i] > tau)) {
      packed[[length(packed) + 1L]] <- cur
      cur <- integer(0); cur_cost <- 0
    }
    cur <- c(cur, i); cur_cost <- cur_cost + cost[i]
  }
  if (length(cur) > 0L) packed[[length(packed) + 1L]] <- cur
  blocks <- c(singles, packed)
  block_cost <- vapply(blocks, function(b) sum(cost[b]), numeric(1))
  o <- order(-block_cost, seq_along(blocks))      # heaviest block first
  structure(blocks[o], tau = tau, n_singletons = length(singles),
            block_cost = block_cost[o])
}

#' Merge the Per-Block Outputs of glow_scan_chunk() Into One
#'
#' Row-binds the flat tables, concatenates the evidence, STAAR-detail and
#' result-object lists, row-binds the error tables, and sums the summary
#' tally, with the wall time of the whole parallel call as \code{runtime_secs}.
#' With \code{labels} (the regions table's labels, in its order) the flat rows,
#' the error rows and the sidecar entries are put back in that order, so the
#' merged output equals the serial one whatever the block membership; without
#' it, block order is preserved.
#'
#' @param outs A list of \code{glow_scan_chunk()} outputs, in block order.
#' @param wall_secs Numeric wall time of the parallel call.
#' @param region_summary Logical, for the empty-table schema.
#' @param labels Character, the unique labels of the regions table in its
#'   order, or NULL.
#' @return One \code{glow_scan_chunk()}-shaped list.
#' @keywords internal
#' @noRd
.rsu_merge_chunk_outputs <- function(outs, wall_secs, region_summary = FALSE,
                                     labels = NULL) {
  tabs <- lapply(outs, `[[`, "results")
  tabs <- tabs[vapply(tabs, nrow, integer(1)) > 0L]
  results <- if (length(tabs) == 0L) .empty_scan_results(region_summary) else
    do.call(rbind, tabs)
  rownames(results) <- NULL
  # The order of the regions table, by label (unique by the D6 check).
  by_label <- function(x, key) {
    if (is.null(labels) || length(key) == 0L) return(x)
    x[order(match(as.character(key), labels)), , drop = FALSE]
  }
  cat_lists <- function(field) {
    parts <- lapply(outs, `[[`, field)
    if (all(vapply(parts, is.null, logical(1)))) return(NULL)
    out <- do.call(c, parts[!vapply(parts, is.null, logical(1))])
    if (!is.null(labels) && length(out) > 0L && !is.null(names(out)))
      out <- out[order(match(names(out), labels))]
    out
  }
  results <- by_label(results, results$label)
  rownames(results) <- NULL
  errors <- do.call(rbind, lapply(outs, `[[`, "errors"))
  errors <- by_label(errors, errors$label)
  rownames(errors) <- NULL
  sums <- do.call(rbind, lapply(outs, `[[`, "summary"))
  summary <- data.frame(
    n_windows           = sum(sums$n_windows),
    n_ok                = sum(sums$n_ok),
    n_skip_empty        = sum(sums$n_skip_empty),
    n_staar_skip_single = sum(sums$n_staar_skip_single),
    n_error             = sum(sums$n_error),
    runtime_secs        = round(wall_secs, 3),
    stringsAsFactors = FALSE)
  list(results = results,
       evidence = cat_lists("evidence"),
       staar_detail = cat_lists("staar_detail"),
       result_objects = cat_lists("result_objects"),
       errors = errors, summary = summary)
}

#' Stack Per-Category Engine Summaries Into One Row
#'
#' Sums the count fields across categories and takes the max runtime, so the
#' coding path returns a single comparable 1-row summary.
#' @keywords internal
#' @noRd
.rsu_stack_summaries <- function(summaries) {
  if (length(summaries) == 0L) {
    return(data.frame(n_windows = 0L, n_ok = 0L, n_skip_empty = 0L,
                      n_staar_skip_single = 0L, n_error = 0L, runtime_secs = 0,
                      stringsAsFactors = FALSE))
  }
  tbl <- do.call(rbind, summaries)
  data.frame(
    n_windows           = sum(tbl$n_windows),
    n_ok                = sum(tbl$n_ok),
    n_skip_empty        = sum(tbl$n_skip_empty),
    n_staar_skip_single = sum(tbl$n_staar_skip_single),
    n_error             = sum(tbl$n_error),
    runtime_secs        = max(tbl$runtime_secs),
    stringsAsFactors = FALSE
  )
}
