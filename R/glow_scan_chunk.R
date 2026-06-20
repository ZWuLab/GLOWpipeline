# ==============================================================================
# Chunk Scan Engine — Run the GLOW Single-Region Test over One Chunk's Windows
# ==============================================================================
#
# The counterpart of STAARpipeline's `Sliding_Window(type = "multiple")`. Opens
# the chunk's GDS once, loops GLOWr::glow_region_test() over a regions table
# under an output policy, and assembles one always-on flat results table plus the
# policy-gated sidecars. Cohort- and cluster-agnostic: resolved objects + a
# regions table in, tables out. The workflow scripts bind it to a cohort, a
# config, and a compute-task id.
#
# Engine-added flat-table columns (on top of the single-region test's canonical row):
#   - midpoint  = as.integer((start + end) / 2)   [Manhattan x-axis, design §10.2]
#   - is_single = (n_after_collapse == 1L)
#   - the 8 STAAR family columns (all-NA when STAAR off or a window staar_skip_single)
# so the flat-table schema is uniform across the chunk whether or not STAAR runs.
#
# EXPORTED FUNCTIONS:
#   - glow_scan_chunk()           Scan one chunk's windows under a policy
#
# INTERNAL HELPER FUNCTIONS:
#   - .scan_flat_row()            Primitive row + engine columns -> one flat row
#   - .empty_scan_results()       0-row flat table with the full schema
#   - .staar_block()              1-row STAAR 8-col block (values or all-NA)


# The 16 GLOW omnibus p-value columns (glow_test output) and the 8 STAAR family
# columns, in their canonical order. These pin the flat-table schema for the
# empty-chunk case (no row to read columns from); a real-GDS test asserts the
# live single-region test emits exactly these names, guarding against drift.
.GLOW_PVALUE_COLS <- c(
  "GLOW_SKAT_BE_N", "GLOW_SKAT_APR_N", "GLOW_SKAT_BE_sparse",
  "GLOW_SKAT_APR_sparse", "GLOW_SKAT_equ",
  "GLOW_Burden_BE", "GLOW_Burden_APR", "GLOW_Burden_equ",
  "GLOW_Fisher_BE_N", "GLOW_Fisher_APR_N", "GLOW_Fisher_BE_sparse",
  "GLOW_Fisher_APR_sparse", "GLOW_Fisher_equ",
  "GLOW_BSF_Omni", "GLOW_SNV_CCT", "GLOW_Omni")

.STAAR_FAMILY_COLS <- c(
  "STAAR_O", "ACAT_O",
  "STAAR_S_1_25", "STAAR_S_1_1",
  "STAAR_B_1_25", "STAAR_B_1_1",
  "STAAR_A_1_25", "STAAR_A_1_1")

.SCAN_ID_COLS    <- c("label", "chr", "start", "end")
.SCAN_COUNT_COLS <- c("n_variants", "n_after_annotation", "n_after_filter",
                      "n_after_ld", "n_after_collapse", "cMAC")

# The 7 region-summary columns the GLOWr single-region test emits when region_summary =
# TRUE (gene scans). Placed after is_single, before the GLOW_* columns. Empty
# when the scan ran with region_summary = FALSE (sliding/coding default).
.SCAN_SUMMARY_COLS <- c("n_collapsed_units", "min_MAF", "max_MAF",
                        "min_MAC", "max_MAC", "min_p_standard", "min_p_SPA")


#################### EXPORTED FUNCTIONS ####################

#' Scan One Genomic Chunk's Sliding Windows with GLOW
#'
#' Opens the chunk's GDS once and runs the single-region GLOW test
#' ([GLOWr::glow_region_test()]) over every window in a regions table,
#' accumulating the always-on flat results table and the artifacts requested by
#' \code{output_policy}. This is the orchestration core of the sliding-window
#' scan; the workflow scripts wrap it per chunk.
#'
#' @details
#' \strong{STAAR on/off} is decided solely by \code{staar_context}: non-NULL
#' turns STAAR on for every window (passed straight to the single-region test); NULL
#' leaves the 8 STAAR columns all-NA. \code{output_policy$write_staar_detail}
#' therefore requires a non-NULL \code{staar_context}, which is validated here
#' (the single cross-check that [glow_output_policy()] deliberately omits).
#'
#' \strong{Per-window status handling} (from the single-region test's \code{$status}):
#' \describe{
#'   \item{\code{"ok"} / \code{"staar_skip_single"}}{A flat row is appended; STAAR
#'     columns carry the single-region test's p-values (\code{"ok"} with STAAR on) or NA
#'     (\code{staar_skip_single}, or STAAR off). Evidence / staar_detail /
#'     result_objects are stashed if the policy enables them and the single-region test
#'     produced them.}
#'   \item{\code{"skip_empty"}}{Tallied; no row.}
#'   \item{\code{"error"}}{Tallied; the message recorded in \code{$errors}; no
#'     row. One window's error never aborts the chunk.}
#' }
#'
#' \strong{Empty chunk} (every window \code{skip_empty} — common at fine scale on
#' sparse data): \code{$results} is a 0-row data.frame carrying the full schema,
#' so a downstream file-presence check stays uniform.
#'
#' @param gds_path Path to the chunk's GDS/aGDS file. Opened read-only once and
#'   reused across windows (closed on exit).
#' @param regions A regions table (the \code{GLOWr::define_regions_window()}
#'   shape: \code{chr}, \code{start}, \code{end}, and ideally \code{label}). All
#'   rows must be on a single chromosome (chunks are per-chromosome, and
#'   \code{reference_medians} is that chromosome's cache).
#' @param null_model,pi_models,pi_features,reference_medians,b_func,b_model,filter_spec,ld_threshold,mac_threshold,collapse_method
#'   Passed through to [GLOWr::glow_region_test()] for every window (exactly one
#'   of \code{b_func}/\code{b_model}). \code{reference_medians} is the chunk
#'   chromosome's cached annotation medians.
#' @param use_spa Logical or NULL (default NULL). Threaded to
#'   [GLOWr::glow_region_test()]: selects GLOW's single-variant score statistics
#'   (NULL = auto: SPA for binary, standard for continuous; TRUE = force SPA,
#'   binary only; FALSE = force standard). GLOW-only — STAAR is decided solely by
#'   \code{staar_context} and always runs non-SPA. Validated up front (a
#'   \code{TRUE} + continuous mismatch fails once, before the GDS is opened).
#' @param z_scale Numeric scalar (default 1). Threaded to
#'   [GLOWr::glow_region_test()]: the M2 single-variant confounding correction
#'   (\code{1/sqrt(inflation_factor)}) applied to each unit's Z before the GLOW
#'   combination. \code{z_scale = 1} is a no-op.
#' @param region_summary Logical (default FALSE). Pass-through to
#'   [GLOWr::glow_region_test()]: when TRUE the flat table gains 7 region-summary
#'   columns (\code{n_collapsed_units}, \code{min_MAF}, \code{max_MAF},
#'   \code{min_MAC}, \code{max_MAC}, \code{min_p_standard}, \code{min_p_SPA})
#'   placed after \code{is_single}. These are the gene scan's per-gene
#'   diagnostic columns; the gene scans turn this on. Default FALSE leaves
#'   the sliding/coding flat schema unchanged.
#' @param output_policy A \code{glow_output_policy} from [glow_output_policy()].
#' @param staar_context A STAAR context (see [GLOWr::glow_region_test()]) or NULL
#'   (STAAR off).
#' @param progress_every Integer. If positive, emit a progress line every N
#'   windows (after window N, 2N, ...). Default 0L (off). Useful for long
#'   chunks where a silent log otherwise looks indistinguishable from a hang.
#' @param progress_callback Optional function(msg) taking the formatted
#'   progress string. Default NULL routes to [message()]. Pass a logger (e.g.
#'   the workflow's `log_msg`) to write progress into the run's per-chunk log.
#' @param verbose Integer (default 0). Threaded to the single-region test.
#'
#' @return A list with:
#' \describe{
#'   \item{results}{A \code{data.frame}, one row per non-skipped window, with the
#'     flat-table schema: \code{label, chr, start, end, midpoint};
#'     \code{n_variants, n_after_annotation, n_after_filter, n_after_ld,
#'     n_after_collapse, cMAC} (where \code{n_after_annotation} = after
#'     QC+type+annotation, before MAF/MAC; in the coding analysis = the
#'     per-category raw SNV count); \code{is_single}; (when
#'     \code{region_summary = TRUE}) the 7 region-summary columns; the 16
#'     \code{GLOW_*} columns; the 8 \code{STAAR_*} columns. 0-row (full schema)
#'     for an empty chunk.}
#'   \item{evidence}{A named list of per-window evidence packets (keyed by window
#'     label) for kept windows, or NULL when \code{write_evidence} is off.}
#'   \item{staar_detail}{A named list of full \code{STAAR()} outputs (keyed by
#'     label) for windows where STAAR ran, or NULL when \code{write_staar_detail}
#'     is off.}
#'   \item{result_objects}{A named list of per-window result rows (keyed by
#'     label) when \code{write_result_objects} is on, else NULL. (The single-region test
#'     exposes the flat row, not the raw \code{glow_test_result} S3 object.)}
#'   \item{errors}{A \code{data.frame(label, message)} of windows that errored
#'     (0-row when none).}
#'   \item{summary}{A 1-row \code{data.frame}: \code{n_windows, n_ok,
#'     n_skip_empty, n_staar_skip_single, n_error, runtime_secs}.}
#' }
#'
#' @seealso [glow_output_policy()], [define_genome_chunks()],
#'   [GLOWr::glow_region_test()]
#'
#' @export
glow_scan_chunk <- function(gds_path, regions,
                            null_model, pi_models, pi_features, reference_medians,
                            b_func = NULL, b_model = NULL,
                            filter_spec, ld_threshold = 0.9, mac_threshold = 10L,
                            collapse_method = "mean",
                            use_spa = NULL,
                            z_scale = 1,
                            region_summary = FALSE,
                            output_policy,
                            staar_context = NULL,
                            progress_every = 0L,
                            progress_callback = NULL,
                            verbose = 0) {

  t0 <- proc.time()[["elapsed"]]

  # ---- 0. Validate (fail fast, before opening the GDS) ----
  if (!inherits(output_policy, "glow_output_policy")) {
    stop("`output_policy` must be a glow_output_policy (see glow_output_policy()).",
         call. = FALSE)
  }
  if (isTRUE(output_policy$write_staar_detail) && is.null(staar_context)) {
    stop("output_policy$write_staar_detail = TRUE requires a non-NULL ",
         "`staar_context` (STAAR must be on to write its detail). Pass a STAAR ",
         "context, or set write_staar_detail = FALSE.", call. = FALSE)
  }
  if (sum(!is.null(b_func), !is.null(b_model)) != 1L) {
    stop("Exactly one of `b_func`, `b_model` must be provided.", call. = FALSE)
  }
  if (!is.numeric(z_scale) || length(z_scale) != 1L || !is.finite(z_scale) ||
      z_scale <= 0) {
    stop("`z_scale` must be a single finite positive number.", call. = FALSE)
  }
  if (!is.null(use_spa) &&
      (!is.logical(use_spa) || length(use_spa) != 1L || is.na(use_spa))) {
    stop("`use_spa` must be NULL, TRUE, or FALSE.", call. = FALSE)
  }
  if (isTRUE(use_spa) && null_model$trait != "binary") {
    stop("`use_spa = TRUE` requires a binary trait (null_model$trait = '",
         null_model$trait, "').", call. = FALSE)
  }
  if (!is.logical(region_summary) || length(region_summary) != 1L ||
      is.na(region_summary)) {
    stop("`region_summary` must be a single non-NA logical (TRUE/FALSE).",
         call. = FALSE)
  }
  if (!is.data.frame(regions) ||
      !all(c("chr", "start", "end") %in% names(regions))) {
    stop("`regions` must be a data.frame with at least chr, start, end columns ",
         "(the define_regions_window() shape).", call. = FALSE)
  }
  if (is.null(reference_medians)) {
    stop("`reference_medians` is NULL (the chunk chromosome's annotation-median ",
         "cache is missing; see the Stage-1 reference-median cache).",
         call. = FALSE)
  }
  if (!"label" %in% names(regions)) {
    regions$label <- paste0("chr", regions$chr, ":", regions$start, "-",
                            regions$end)
  }
  if (nrow(regions) > 0L && length(unique(as.character(regions$chr))) > 1L) {
    stop("`regions` spans multiple chromosomes (",
         paste(unique(as.character(regions$chr)), collapse = ", "),
         "); a chunk is single-chromosome and `reference_medians` is one ",
         "chromosome's cache.", call. = FALSE)
  }
  progress_every <- as.integer(progress_every)
  if (is.na(progress_every) || progress_every < 0L) {
    stop("`progress_every` must be a non-negative integer.", call. = FALSE)
  }
  if (!is.null(progress_callback) && !is.function(progress_callback)) {
    stop("`progress_callback` must be a function(msg) or NULL.", call. = FALSE)
  }
  emit_progress <- if (progress_every > 0L) {
    fn <- if (is.null(progress_callback)) function(msg) message(msg)
          else progress_callback
    function(i, n) fn(sprintf("progress: %d/%d windows", i, n))
  } else NULL

  # ---- 1. Open the GDS once for the whole chunk ----
  gds <- seqOpen(gds_path, readonly = TRUE)
  on.exit(seqClose(gds), add = TRUE)

  # ---- 2. Loop windows ----
  n_windows <- nrow(regions)
  row_list      <- vector("list", n_windows)   # flat rows (NULL for skips/errors)
  evidence      <- if (isTRUE(output_policy$write_evidence)) list() else NULL
  staar_detail  <- if (isTRUE(output_policy$write_staar_detail)) list() else NULL
  result_objs   <- if (isTRUE(output_policy$write_result_objects)) list() else NULL
  err_labels    <- character(0)
  err_messages  <- character(0)
  n_ok <- n_skip_empty <- n_staar_skip_single <- n_error <- 0L

  for (i in seq_len(n_windows)) {
    region <- as.list(regions[i, , drop = FALSE])
    label  <- as.character(region$label)

    r <- glow_region_test(
      gds, region, filter_spec, null_model,
      pi_models = pi_models, pi_features = pi_features,
      reference_medians = reference_medians,
      b_func = b_func, b_model = b_model,
      ld_threshold = ld_threshold, mac_threshold = mac_threshold,
      collapse_method = collapse_method,
      use_spa = use_spa,
      z_scale = z_scale,
      staar = staar_context,
      return_evidence = isTRUE(output_policy$write_evidence),
      region_summary = region_summary,
      verbose = verbose
    )

    if (r$status %in% c("ok", "staar_skip_single")) {
      if (identical(r$status, "ok")) n_ok <- n_ok + 1L else
        n_staar_skip_single <- n_staar_skip_single + 1L
      row_list[[i]] <- .scan_flat_row(r$result, r$staar, region_summary)
      if (!is.null(evidence) && !is.null(r$evidence))
        evidence[[label]] <- r$evidence
      if (!is.null(staar_detail) && !is.null(r$staar_detail))
        staar_detail[[label]] <- r$staar_detail
      if (!is.null(result_objs))
        result_objs[[label]] <- r$result
    } else if (identical(r$status, "skip_empty")) {
      n_skip_empty <- n_skip_empty + 1L
    } else {                                   # "error"
      n_error <- n_error + 1L
      err_labels   <- c(err_labels, label)
      err_messages <- c(err_messages,
                        if (is.null(r$message)) NA_character_ else r$message)
    }

    if (!is.null(emit_progress) && (i %% progress_every == 0L)) {
      emit_progress(i, n_windows)
    }
  }

  # ---- 3. Assemble outputs ----
  kept <- row_list[!vapply(row_list, is.null, logical(1))]
  results <- if (length(kept) == 0L) .empty_scan_results(region_summary) else
    do.call(rbind, kept)
  rownames(results) <- NULL

  summary <- data.frame(
    n_windows     = n_windows,
    n_ok          = n_ok,
    n_skip_empty  = n_skip_empty,
    n_staar_skip_single = n_staar_skip_single,
    n_error       = n_error,
    runtime_secs  = round(proc.time()[["elapsed"]] - t0, 3),
    stringsAsFactors = FALSE
  )
  errors <- data.frame(label = err_labels, message = err_messages,
                       stringsAsFactors = FALSE)

  list(results = results, evidence = evidence, staar_detail = staar_detail,
       result_objects = result_objs, errors = errors, summary = summary)
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Build one flat-table row from a single-region test result + STAAR p-values
#'
#' Reorders the single-region test's canonical row into the flat-table schema, inserting
#' the engine-added \code{midpoint} (after \code{end}) and \code{is_single}
#' (after \code{cMAC}) columns and appending the 8-column STAAR block. When
#' \code{region_summary} is TRUE the single-region test's row also carries the 7
#' region-summary columns; they are placed after \code{is_single} and before the
#' \code{GLOW_*} p-value columns (and excluded from \code{glow_cols} so they are
#' not appended twice).
#'
#' @param result A 1-row single-region test \code{$result} data.frame.
#' @param staar The single-region test's \code{$staar} (named list of 8) or NULL.
#' @param region_summary Logical. Whether the single-region test emitted the 7
#'   region-summary columns (so they are carried in their canonical slot).
#' @return A 1-row data.frame in the flat-table schema.
#' @keywords internal
#' @noRd
.scan_flat_row <- function(result, staar, region_summary = FALSE) {
  glow_cols <- setdiff(names(result),
                       c(.SCAN_ID_COLS, .SCAN_COUNT_COLS, .SCAN_SUMMARY_COLS))

  front <- data.frame(
    label    = as.character(result$label),
    chr      = as.character(result$chr),
    start    = as.integer(result$start),
    end      = as.integer(result$end),
    midpoint = as.integer((as.numeric(result$start) + as.numeric(result$end)) / 2),
    stringsAsFactors = FALSE
  )
  counts <- result[, .SCAN_COUNT_COLS, drop = FALSE]
  flag   <- data.frame(is_single = isTRUE(result$n_after_collapse == 1L))
  glow   <- result[, glow_cols, drop = FALSE]

  # NULL summary_block must be dropped from the cbind (cbind() treats a NULL
  # argument as a 0-row operand and errors on the row mismatch), so assemble the
  # parts list and filter NULL before binding.
  parts <- list(front, counts, flag,
                if (isTRUE(region_summary))
                  result[, .SCAN_SUMMARY_COLS, drop = FALSE] else NULL,
                glow, .staar_block(staar))
  do.call(cbind, parts[!vapply(parts, is.null, logical(1))])
}


#' Build the 1-row STAAR 8-column block
#'
#' @param staar A named list of 8 STAAR family p-values, or NULL (STAAR off).
#' @return A 1-row data.frame with the 8 \code{.STAAR_FAMILY_COLS} (all-NA when
#'   \code{staar} is NULL).
#' @keywords internal
#' @noRd
.staar_block <- function(staar) {
  if (is.null(staar)) {
    vals <- as.list(rep(NA_real_, length(.STAAR_FAMILY_COLS)))
  } else {
    vals <- lapply(.STAAR_FAMILY_COLS, function(nm) {
      v <- staar[[nm]]
      if (is.null(v)) NA_real_ else as.numeric(v)
    })
  }
  out <- as.data.frame(stats::setNames(vals, .STAAR_FAMILY_COLS),
                       stringsAsFactors = FALSE)
  out
}


#' A 0-row flat results table carrying the full schema
#'
#' @param region_summary Logical. When TRUE, include the 7 region-summary columns
#'   (after \code{is_single}) so the empty-chunk schema matches the populated one.
#' @return A 0-row data.frame with every flat-table column, correctly typed.
#' @keywords internal
#' @noRd
.empty_scan_results <- function(region_summary = FALSE) {
  base <- data.frame(
    label              = character(0),
    chr                = character(0),
    start              = integer(0),
    end                = integer(0),
    midpoint           = integer(0),
    n_variants         = integer(0),
    n_after_annotation = integer(0),
    n_after_filter     = integer(0),
    n_after_ld         = integer(0),
    n_after_collapse   = integer(0),
    cMAC               = numeric(0),
    is_single          = logical(0),
    stringsAsFactors = FALSE
  )
  if (isTRUE(region_summary)) {
    summary_block <- data.frame(
      n_collapsed_units = integer(0),
      min_MAF           = numeric(0),
      max_MAF           = numeric(0),
      min_MAC           = numeric(0),
      max_MAC           = numeric(0),
      min_p_standard    = numeric(0),
      min_p_SPA         = numeric(0),
      stringsAsFactors = FALSE
    )
    base <- cbind(base, summary_block)
  }
  pcols <- c(.GLOW_PVALUE_COLS, .STAAR_FAMILY_COLS)
  pblock <- as.data.frame(
    stats::setNames(rep(list(numeric(0)), length(pcols)), pcols),
    stringsAsFactors = FALSE
  )
  cbind(base, pblock)
}
