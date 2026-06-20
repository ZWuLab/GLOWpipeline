# ==============================================================================
# Output-Policy Constructor for the GLOW Window Scan
# ==============================================================================
#
# A small, validated value object that records WHICH per-window artifacts a scan
# persists — passed as ONE argument to glow_scan_chunk() instead of four loose
# flags, and printed (S3) into the run log as a one-line retention record. It
# does NOT inspect window count or output size. What keeps a fine-scale run (up
# to ~1.4M windows) from writing hundreds of GB of sidecars is that every sidecar
# flag DEFAULTS OFF and glow_scan_chunk() gates *accumulation* on it — with the
# default policy the engine never even builds the per-window evidence/detail/
# §7.2, §9.2). This constructor's OWN checks are only: each flag is a scalar
# non-NA logical, and write_flat_table must be TRUE (the one load-bearing output
# cannot be disabled). The size safety lives in the defaults + the engine, not
# in a check here.
#
# The policy governs PERSISTENCE ONLY. Whether STAAR runs is signalled
# separately by passing a `staar_context` to glow_scan_chunk() (the single
# source of truth for "STAAR on"); the cross-check that `write_staar_detail`
# requires such a context lives in glow_scan_chunk(), not here, because the
# context is engine input, not a policy field.
#
# EXPORTED FUNCTIONS:
#   - glow_output_policy()        Construct/validate the policy object
#   - print.glow_output_policy()  Readable one-screen summary (S3)
#
# INTERNAL HELPER FUNCTIONS:
#   - .check_scalar_logical()     Validate a length-1, non-NA logical flag


#################### EXPORTED FUNCTIONS ####################

#' Construct a GLOW Window-Scan Output Policy
#'
#' Builds and validates the object that tells [glow_scan_chunk()] which
#' per-window artifacts to persist. The always-on flat results table is
#' mandatory (it is the one load-bearing output); every other artifact is
#' opt-in and defaults off, because at fine window scale the per-window sidecars
#' would dominate storage.
#'
#' @details
#' This object records persistence intent only. It deliberately does \emph{not}
#' carry a STAAR switch: whether STAAR runs is decided by passing a
#' \code{staar_context} to [glow_scan_chunk()] (a non-NULL context turns STAAR
#' on for both \code{glow_scan_chunk()} and the single-region GLOW test). Consequently
#' \code{write_staar_detail = TRUE} is only \emph{meaningful} alongside a STAAR
#' context, and that cross-check is enforced by [glow_scan_chunk()] — not here —
#' so there is a single source of truth for "STAAR on". Constructing a policy
#' with \code{write_staar_detail = TRUE} and no context is therefore allowed and
#' does not error at construction time.
#'
#' The size safety at fine window scale comes from the sidecar flags
#' \emph{defaulting off} plus [glow_scan_chunk()] gating in-memory accumulation on
#' them (with the default policy \code{glow_scan_chunk()} never builds the per-window sidecar
#' lists) — not from any size or window-count check here. This constructor only
#' validates each flag is a scalar non-NA logical and enforces
#' \code{write_flat_table = TRUE}.
#'
#' @param write_flat_table Logical (default \code{TRUE}). The compact one-row-
#'   per-window results table. Must be \code{TRUE}; \code{FALSE} is rejected.
#' @param write_evidence Logical (default \code{FALSE}). Per-window evidence
#'   packets (post-filter variant table + per-unit B/PI/Z/weights). ~1.6 KB per
#'   window, so tens of GB at fine WGS scale — keep off there.
#' @param write_staar_detail Logical (default \code{FALSE}). The full per-window
#'   \code{STAAR::STAAR} output grid. Requires a \code{staar_context} at the
#'   engine (validated in [glow_scan_chunk()]).
#' @param write_result_objects Logical (default \code{FALSE}). The per-window
#'   \code{glow_test_result} S3 objects (RDS). Redundant with the flat table for
#'   most uses.
#'
#' @return An S3 object of class \code{"glow_output_policy"}: a named list of the
#'   four logical flags.
#'
#' @seealso [glow_scan_chunk()]
#'
#' @examples
#' # Lean default: only the flat table.
#' glow_output_policy()
#'
#' # Detailed policy: flat table + evidence packets (window count modest).
#' glow_output_policy(write_evidence = TRUE)
#'
#' @export
glow_output_policy <- function(write_flat_table     = TRUE,
                               write_evidence       = FALSE,
                               write_staar_detail   = FALSE,
                               write_result_objects = FALSE) {

  .check_scalar_logical(write_flat_table,     "write_flat_table")
  .check_scalar_logical(write_evidence,       "write_evidence")
  .check_scalar_logical(write_staar_detail,   "write_staar_detail")
  .check_scalar_logical(write_result_objects, "write_result_objects")

  if (!isTRUE(write_flat_table)) {
    stop("`write_flat_table` must be TRUE: the flat results table is the one ",
         "load-bearing output of a window scan and cannot be disabled. ",
         "(All other artifacts are opt-in.)", call. = FALSE)
  }

  structure(
    list(write_flat_table     = TRUE,
         write_evidence       = write_evidence,
         write_staar_detail   = write_staar_detail,
         write_result_objects = write_result_objects),
    class = "glow_output_policy"
  )
}


#' @rdname glow_output_policy
#' @param x A \code{glow_output_policy} object.
#' @param ... Ignored.
#' @export
print.glow_output_policy <- function(x, ...) {
  on_off <- function(flag) if (isTRUE(flag)) "on" else "off"
  cat("<glow_output_policy>\n")
  cat("  flat_table     :", on_off(x$write_flat_table),
      "(always; the load-bearing output)\n")
  cat("  evidence       :", on_off(x$write_evidence), "\n")
  cat("  staar_detail   :", on_off(x$write_staar_detail),
      "(needs a staar_context at the engine)\n")
  cat("  result_objects :", on_off(x$write_result_objects), "\n")
  invisible(x)
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Validate a length-1, non-NA logical flag
#'
#' @param x The value to check.
#' @param name The argument name (for the error message).
#' @keywords internal
#' @noRd
.check_scalar_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", name, "` must be a single non-NA logical (TRUE/FALSE), not ",
         if (length(x) != 1L) paste0("length ", length(x), " ") else "",
         "`", paste(format(x), collapse = ", "), "`.", call. = FALSE)
  }
  invisible(TRUE)
}
