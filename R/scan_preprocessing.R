# ==============================================================================
# File: scan_preprocessing.R
# Purpose: Resolve and validate the three variant-preprocessing settings a scan reads from
#   its configuration (ld_threshold, mac_threshold, collapse_method). A missing field takes
#   GLOWr's default (GLOWr::glow_preprocessing_defaults(), the one source every layer
#   shares); an invalid value is rejected with a message naming the field, before any
#   output is written. Used by prepare_scan_run() (Stage 1) and run_scan_unit() (Stage 2).
#
# INTERNAL HELPERS:
#   - .resolve_preprocessing()   NULL -> default; validate; typed list

#################### INTERNAL HELPER FUNCTIONS ####################

#' Resolve and Validate the Variant-Preprocessing Settings
#'
#' A NULL (absent) field takes the value of \code{GLOWr::glow_preprocessing_defaults()}.
#' Checks: \code{ld_threshold} one finite number in (0, 1] (the prune is on the absolute
#' correlation); \code{mac_threshold} one non-negative whole number (0 disables collapsing);
#' \code{collapse_method} one of "mean", "max", "sum".
#'
#' @param ld_threshold,mac_threshold,collapse_method The configuration's values, or NULL.
#' @param where Character naming the source, for the error message.
#' @return A list with the three fields, typed (numeric, integer, character).
#' @keywords internal
#' @noRd
.resolve_preprocessing <- function(ld_threshold = NULL, mac_threshold = NULL,
                                   collapse_method = NULL, where = "the configuration") {
  d <- GLOWr::glow_preprocessing_defaults()
  if (is.null(ld_threshold))    ld_threshold    <- d$ld_threshold
  if (is.null(mac_threshold))   mac_threshold   <- d$mac_threshold
  if (is.null(collapse_method)) collapse_method <- d$collapse_method
  if (!is.numeric(ld_threshold) || length(ld_threshold) != 1L || !is.finite(ld_threshold) ||
      ld_threshold <= 0 || ld_threshold > 1) {
    stop("`ld_threshold` in ", where, " must be one finite number in (0, 1], the absolute ",
         "correlation above which one variant of a pair is pruned (got ",
         deparse(ld_threshold), ").", call. = FALSE)
  }
  if (!is.numeric(mac_threshold) || length(mac_threshold) != 1L || !is.finite(mac_threshold) ||
      mac_threshold < 0 || mac_threshold != round(mac_threshold) ||
      mac_threshold > .Machine$integer.max) {
    stop("`mac_threshold` in ", where, " must be one non-negative whole number, the cohort ",
         "MAC below which variants are collapsed (0 disables collapsing; got ",
         deparse(mac_threshold), ").", call. = FALSE)
  }
  if (!is.character(collapse_method) || length(collapse_method) != 1L ||
      !collapse_method %in% c("mean", "max", "sum")) {
    stop("`collapse_method` in ", where, " must be \"mean\", \"max\" or \"sum\", the ",
         "aggregation of the B and PI weights of a collapsed group (got ",
         deparse(collapse_method), ").", call. = FALSE)
  }
  list(ld_threshold = as.numeric(ld_threshold), mac_threshold = as.integer(mac_threshold),
       collapse_method = collapse_method)
}
