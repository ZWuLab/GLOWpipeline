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
#   - .rsu_merge_chunk_outputs() Merge per-block engine outputs (n_cores > 1).
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
#' @param region_type Character: \code{"gene"}, \code{"window"}, or
#'   \code{"coding"}.
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
#'   regions are split into blocks of 25 consecutive regions, each block runs
#'   \code{\link{glow_scan_chunk}} in a forked worker with its own GDS handle
#'   (\code{parallel::mclapply}, dynamic scheduling so a block with a large
#'   gene does not hold the others), and the block outputs are merged in block
#'   order, so the written table equals the single-core one. The progress
#'   callback is off inside the workers; the parent logs one line per finished
#'   block. Forking is unavailable on Windows (a value above 1 warns and runs
#'   serially). Memory scales with the workers, since each holds its own
#'   genotype and LD matrices.
#' @param verbose Integer threaded to the engine. Default 0.
#'
#' @section Large genes:
#' For \code{region_type = "gene"} the unit's variant index is passed to
#' \code{\link{build_scan_regions}}, which tiles every gene with more GDS
#' records than \code{config$gene_max_records} into segments of
#' \code{config$gene_segment_bp} (see that function). The log states how many
#' genes were segmented and into how many segments.
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

  region_type <- match.arg(region_type, c("gene", "window", "coding"))
  log_msg <- if (is.null(logger)) function(...) invisible(NULL) else logger

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
  built   <- build_scan_regions(region_type, chr, config, chunk = chunk,
                                variant_index = variant_index)
  regions <- built$regions
  n_seg_genes <- attr(regions, "n_genes_segmented")
  if (!is.null(n_seg_genes)) {
    log_msg(sprintf(paste0("gene segmentation: %d gene(s) above gene_max_records = %s ",
                           "records tiled into %d segment(s) of %s bp"),
                    n_seg_genes, format(.rsu_get(config, "gene_max_records")),
                    attr(regions, "n_segments"),
                    format(.rsu_null_or(.rsu_get(config, "gene_segment_bp"), 10000L))))
  }
  if (!is.na(max_regions)) {
    regions <- regions[seq_len(min(max_regions, nrow(regions))), , drop = FALSE]
    log_msg(sprintf("max_regions %d: limiting to %d regions",
                    max_regions, nrow(regions)))
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

  run_one <- function(filter_spec) {
    if (n_cores == 1L || nrow(regions) <= 1L) {
      return(engine_call(regions, filter_spec, in_parent = TRUE))
    }
    # Blocks of 25 consecutive regions, one forked worker per block, dynamic
    # scheduling (mc.preschedule = FALSE) so a block holding a large gene does
    # not hold the rest. The outputs come back in block order and are merged
    # in that order, so the table is the serial one.
    block_size <- 25L
    blocks <- split(seq_len(nrow(regions)),
                    ceiling(seq_len(nrow(regions)) / block_size))
    log_msg(sprintf("n_cores = %d: %d regions in %d blocks of up to %d",
                    n_cores, nrow(regions), length(blocks), block_size))
    t_par <- proc.time()[["elapsed"]]
    outs <- parallel::mclapply(
      blocks,
      function(idx) engine_call(regions[idx, , drop = FALSE], filter_spec,
                                in_parent = FALSE),
      mc.cores = n_cores, mc.preschedule = FALSE)
    failed <- vapply(outs, function(o) inherits(o, "try-error"), logical(1))
    if (any(failed)) {
      k <- which(failed)[1L]
      stop(sprintf("run_scan_unit: worker for block %d (regions %d-%d) failed: %s",
                   k, min(blocks[[k]]), max(blocks[[k]]),
                   as.character(outs[[k]])), call. = FALSE)
    }
    for (k in seq_along(outs)) {
      sm <- outs[[k]]$summary
      log_msg(sprintf("block %d/%d: regions %d-%d, n_ok=%d n_skip_empty=%d n_error=%d %.1fs",
                      k, length(outs), min(blocks[[k]]), max(blocks[[k]]),
                      sm$n_ok, sm$n_skip_empty, sm$n_error, sm$runtime_secs))
    }
    .rsu_merge_chunk_outputs(outs, wall_secs = proc.time()[["elapsed"]] - t_par,
                             region_summary = region_summary)
  }

  # ---- Run: ONCE (gene/window) or ONCE PER CATEGORY (coding) ----
  evidence_acc <- list(); staar_detail_acc <- list(); result_objects_acc <- list()

  if (region_type == "coding") {
    category_labels <- names(built$filters)
    acc <- list(); err_rows <- list(); summaries <- list()
    log_msg(sprintf("Starting coding chr%s: %d genes, %d categories",
                    chr, nrow(regions), length(category_labels)))
    for (category in category_labels) {
      out <- run_one(built$filters[[category]])
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
      if ("start" %in% names(res)) {
        res <- res[order(res$start, res$category), , drop = FALSE]
      }
      rownames(res) <- NULL
    }
    errors <- if (length(err_rows) > 0L) do.call(rbind, err_rows) else
      data.frame(category = character(0L), label = character(0L),
                 message = character(0L))
    # One stacked summary row (sum the per-category tallies, max runtime).
    summary <- .rsu_stack_summaries(summaries)
  } else {
    out <- run_one(built$filters)
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

#' Merge the Per-Block Outputs of glow_scan_chunk() Into One
#'
#' Row-binds the flat tables, concatenates the evidence, STAAR-detail and
#' result-object lists, row-binds the error tables, and sums the summary
#' tally, with the wall time of the whole parallel call as \code{runtime_secs}.
#' Block order is preserved, so the merged table equals the serial one.
#'
#' @param outs A list of \code{glow_scan_chunk()} outputs, in block order.
#' @param wall_secs Numeric wall time of the parallel call.
#' @param region_summary Logical, for the empty-table schema.
#' @return One \code{glow_scan_chunk()}-shaped list.
#' @keywords internal
#' @noRd
.rsu_merge_chunk_outputs <- function(outs, wall_secs, region_summary = FALSE) {
  tabs <- lapply(outs, `[[`, "results")
  tabs <- tabs[vapply(tabs, nrow, integer(1)) > 0L]
  results <- if (length(tabs) == 0L) .empty_scan_results(region_summary) else
    do.call(rbind, tabs)
  rownames(results) <- NULL
  cat_lists <- function(field) {
    parts <- lapply(outs, `[[`, field)
    if (all(vapply(parts, is.null, logical(1)))) return(NULL)
    do.call(c, parts[!vapply(parts, is.null, logical(1))])
  }
  errors <- do.call(rbind, lapply(outs, `[[`, "errors"))
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
