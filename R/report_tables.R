# ==============================================================================
# Report-Table Builders
# ==============================================================================
#
# Shared computation behind the post-analysis report tables: resolve the GLOW
# genomic-control calibration fields read by the coding diagnostics / compare
# reports, and build the per-category lambda_GC markdown table (common-support
# view) those reports render. Factored here so the reports cannot drift in HOW
# they read the snapshot / present the table.
#
# EXPORTED FUNCTIONS:
#   - gc_calibration_fields()   Resolve GC-calibration report fields from config
#   - build_lambda_table_md()   Per-category lambda_GC markdown table


#################### EXPORTED MAIN FUNCTIONS ####################

#' Resolve GLOW Genomic-Control Calibration Report Fields
#'
#' Resolves the GLOW-only genomic-control (GC) calibration fields that the
#' coding diagnostics and compare reports both read for their run-identity
#' "GC calibration (GLOW-only)" row. The data-prep is identical in both reports
#' (the rendered cell text differs, so each caller still formats its own cell);
#' this factors only the shared computation so the two reports cannot drift in
#' how they read the snapshot.
#'
#' @param config The resolved run configuration
#'   (\code{\link{load_config_snapshot}}).
#' @param null_or A NULL-coalescing \code{function(a, b) -> if (is.null(a)) b
#'   else a}, passed by the caller. Kept as an argument (rather than relying on
#'   a sourced global) so this helper is self-contained.
#'
#' @return A list with \code{inflation} (\code{config$inflation_factor} or
#'   \code{1.0} when absent; \code{1.0} means calibration OFF), \code{method}
#'   (\code{config$calibration$method[1]} or \code{"manual"} when absent), and
#'   \code{z_scale} (\code{1 / sqrt(inflation)}).
#'
#' @export
gc_calibration_fields <- function(config, null_or) {
  infl <- null_or(config$inflation_factor, 1.0)
  cal_meth <- if (!is.null(config$calibration$method))
    as.character(config$calibration$method)[1] else "manual"
  list(inflation = infl, method = cal_meth, z_scale = 1 / sqrt(infl))
}


#' Build the Per-Category Lambda_GC Markdown Table
#'
#' Builds the markdown per-category \eqn{\lambda_{\mathrm{GC}}} table (the
#' common-support view) that the coding diagnostics and compare reports render
#' from \code{lambda_by_category.csv}. The shared core (filter to common
#' support, fixed method-column order, category rows sorted with ALL appended
#' last, the negative-control tag) is factored here and the per-report
#' differences are parameterised by \code{show_n} / \code{append_all_row} /
#' \code{extra_note}.
#'
#' @details
#' Returns the table block from the header row through the trailing note(s): the
#' caller still emits the section header, leading blank line, the
#' \code{lam}-is-NULL else-branch, and the trailing blank line, so the
#' surrounding structure of each report is unchanged.
#'
#' @param lam data.frame read from \code{lambda_by_category.csv}. Required
#'   columns: \code{category}, \code{method}, \code{support}, \code{n},
#'   \code{lambda_gc}. Must be non-NULL.
#' @param methods Character vector of method columns, in display order. Only
#'   those also present in the common-support rows of \code{lam} are shown
#'   (intersect, preserving this order).
#' @param neg_control Category name tagged with the negative-control marker.
#'   Default \code{"synonymous"}.
#' @param support Value of \code{lam$support} selecting the common-support rows.
#'   Default \code{"common_ge2"}.
#' @param show_n Logical. If \code{TRUE}, cells read \code{"X.XXX (n=NN)"} and
#'   missing lambdas read \code{"NA (n=NN)"}; if \code{FALSE} (default), cells
#'   read \code{"X.XXX"} / \code{"NA"}.
#' @param append_all_row Logical. If \code{TRUE} (default), an \code{"ALL"} row
#'   is appended after the sorted non-ALL categories.
#' @param extra_note Optional character vector of trailing note line(s) appended
#'   after a blank line below the table. \code{NULL} (default) appends no note.
#'
#' @return Character vector of markdown lines (header, separator, one line per
#'   category, then a blank line and the \code{extra_note} lines when
#'   \code{extra_note} is non-NULL).
#'
#' @export
build_lambda_table_md <- function(lam, methods,
                                   neg_control = "synonymous",
                                   support = "common_ge2",
                                   show_n = FALSE,
                                   append_all_row = TRUE,
                                   extra_note = NULL) {
  lc <- lam[lam$support == support, , drop = FALSE]
  meths <- intersect(methods, unique(lc$method))

  out <- c(
    sprintf("| category | %s |", paste(meths, collapse = " | ")),
    sprintf("|---|%s", paste(rep("---|", length(meths)), collapse = ""))
  )

  cats <- sort(setdiff(unique(lc$category), "ALL"))
  if (append_all_row) cats <- c(cats, "ALL")
  for (cat0 in cats) {
    cells <- vapply(meths, function(m) {
      v <- lc$lambda_gc[lc$category == cat0 & lc$method == m]
      if (show_n) {
        n <- lc$n[lc$category == cat0 & lc$method == m]
        if (!length(v) || is.na(v[1])) sprintf("NA (n=%s)", if (length(n)) n[1] else 0)
        else sprintf("%.3f (n=%d)", v[1], n[1])
      } else {
        if (!length(v) || is.na(v[1])) "NA" else sprintf("%.3f", v[1])
      }
    }, character(1))
    tag <- if (identical(cat0, neg_control)) " **[neg ctrl]**" else ""
    out <- c(out, sprintf("| %s%s | %s |", cat0, tag, paste(cells, collapse = " | ")))
  }

  if (!is.null(extra_note)) out <- c(out, "", extra_note)
  out
}
