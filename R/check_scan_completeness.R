# ==============================================================================
# Genome-Scan Completeness Gate
# ==============================================================================
#
# A genome-wide scan submits one compute task per genome chunk. The genome-wide
# summary must NOT silently aggregate a PARTIAL scan: the Bonferroni threshold is
# alpha / n_windows over ALL tested windows, so a missing chunk shrinks the
# denominator and inflates the count of "significant" windows.
# check_scan_completeness() verifies every expected chunk produced its result
# file, reusing parse_chunk_log() for the per-chunk Done/error signal.
#
# EXPORTED FUNCTIONS:
#   - check_scan_completeness()   Expected-vs-present chunk completeness gate
#
# INTERNAL HELPER FUNCTIONS:
#   - collapse_int_ranges()       c(3,17,18,19) -> "3, 17-19" (not exported)


#################### EXPORTED MAIN FUNCTIONS ####################

#' Check Whether a Genome Scan Is Complete
#'
#' Verifies that every expected genome chunk produced its result file, so the
#' genome-wide summary does not silently aggregate a partial scan (which would
#' shrink the Bonferroni denominator and inflate the count of significant
#' windows).
#'
#' @details
#' Expected chunk ids come from \code{<output_dir>/shared/01-genome_chunks.rds};
#' present chunks are those with a
#' \code{results/scan_chunk_<id>.<ext>} file. For every expected chunk, the
#' per-chunk log (if present) is parsed via \code{\link{parse_chunk_log}} to
#' detect whether the run reached "Done" and how many per-window errors it had.
#'
#' @param output_dir Character. The run's \code{outputs/} directory.
#' @param output_format Character. \code{"csv"} (default) or \code{"fst"};
#'   selects the result-file extension.
#'
#' @return A list with: \code{complete} (TRUE iff every expected chunk has a
#'   result file), \code{n_expected}, \code{n_present}, \code{missing_ids},
#'   \code{not_done_ids} (result exists but log never reached "Done";
#'   advisory), \code{n_window_errors} (total \code{n_error} across finished
#'   chunks), and \code{status} (per-chunk data.frame: \code{chunk_id},
#'   \code{has_result}, \code{reached_done}, \code{n_error}).
#'
#' @seealso \code{\link{parse_chunk_log}}
#' @export
check_scan_completeness <- function(output_dir, output_format = "csv") {
  shared_dir  <- file.path(output_dir, "shared")
  results_dir <- file.path(output_dir, "results")
  logs_dir    <- file.path(output_dir, "logs")

  ct_path <- file.path(shared_dir, "01-genome_chunks.rds")
  if (!file.exists(ct_path))
    stop("Chunk table not found: ", ct_path,
         "\n  Run Stage 1 (01-prepare-shared.R) for this run first.")
  chunk_table  <- readRDS(ct_path)
  expected_ids <- sort(as.integer(chunk_table$chunk_id))
  n_expected   <- length(expected_ids)

  ext       <- if (identical(output_format, "fst")) "fst" else "csv"
  res_files <- list.files(results_dir,
                          pattern = sprintf("^scan_chunk_\\d+\\.%s$", ext))
  present_ids <- if (length(res_files) == 0L) integer(0L) else
    sort(as.integer(sub(sprintf("^scan_chunk_0*(\\d+)\\.%s$", ext), "\\1",
                        res_files)))
  missing_ids <- setdiff(expected_ids, present_ids)

  has_result   <- expected_ids %in% present_ids
  reached_done <- rep(NA, n_expected)
  n_error      <- rep(NA_integer_, n_expected)
  for (i in seq_along(expected_ids)) {
    lp <- file.path(logs_dir, sprintf("log_chunk_%04d.txt", expected_ids[i]))
    if (!file.exists(lp)) next
    pl <- tryCatch(parse_chunk_log(lp), error = function(e) NULL)
    if (is.null(pl)) next
    reached_done[i] <- !is.na(pl$summary$n_windows)  # NA n_windows => no "Done"
    n_error[i]      <- pl$summary$n_error
  }
  not_done_ids    <- expected_ids[has_result & !is.na(reached_done) & !reached_done]
  n_window_errors <- sum(n_error, na.rm = TRUE)

  list(complete        = length(missing_ids) == 0L,
       n_expected      = n_expected,
       n_present       = length(present_ids),
       missing_ids     = missing_ids,
       not_done_ids    = not_done_ids,
       n_window_errors = n_window_errors,
       status          = data.frame(chunk_id     = expected_ids,
                                    has_result   = has_result,
                                    reached_done = reached_done,
                                    n_error      = n_error,
                                    stringsAsFactors = FALSE))
}


#################### INTERNAL HELPER FUNCTIONS ####################

# Collapse a sorted integer vector into compact range strings for readable
# reporting: c(3, 17, 18, 19, 200) -> "3, 17-19, 200". "" for empty input.
# Internal (not exported).
#' @noRd
collapse_int_ranges <- function(ids) {
  if (length(ids) == 0L) return("")
  ids    <- sort(unique(as.integer(ids)))
  brk    <- diff(ids) != 1L
  starts <- ids[c(TRUE, brk)]
  ends   <- ids[c(brk, TRUE)]
  parts  <- ifelse(starts == ends, as.character(starts),
                   sprintf("%d-%d", starts, ends))
  paste(parts, collapse = ", ")
}
