# ==============================================================================
# QQ Plots + Per-Test Lambda Summary
# ==============================================================================
#
# Computes per-test genomic-control lambda values for a vector of test columns,
# writes a single-series QQ plot for the headline omnibus, optionally writes a
# multi-panel per-family QQ plot, and (optionally) saves the lambdas to an RDS
# for downstream consumption. All written artifacts are tagged with `label` so
# the caller can run the function twice without filename collisions. This is the
# one side-effecting summary writer in the summary layer.
#
# EXPORTED FUNCTIONS:
#   - summarize_qq_and_lambdas()   Per-test lambda + QQ PDFs + lambda RDS


#################### EXPORTED MAIN FUNCTIONS ####################

#' Summarize QQ Plots and Per-Test Lambda Values
#'
#' Computes per-test \eqn{\lambda_{\mathrm{GC}}} values for a vector of test
#' columns, writes a single-series QQ plot for the headline omnibus, optionally
#' writes a multi-panel per-family QQ plot, and (optionally) saves the lambdas
#' to an RDS for downstream consumption. All written artifacts are tagged with
#' \code{label} so the caller can run the function twice (e.g. on the full set
#' and on a subset) without filename collisions.
#'
#' @details
#' This helper does no row filtering; callers pass an already-filtered
#' \code{df}. Per-test lambda is \code{\link[GLOWr]{compute_lambda_gc}} of each
#' column (\code{NA} flows through, e.g. for a column with fewer than 50 non-NA
#' p-values). The omnibus QQ uses \code{\link[GLOWr]{plot_qq}}.
#'
#' @param df data.frame of results, one row per analysis unit; must contain
#'   every column in \code{test_columns} and \code{omni_test}.
#' @param label Short string used as the filename / title suffix.
#' @param output_dir Directory to write the PDFs and RDS into (created if
#'   missing).
#' @param test_columns Character vector of column names whose p-values get a
#'   lambda entry. May be \code{NULL} only if \code{families} is given (then
#'   derived as \code{unlist(families)}).
#' @param omni_test Character. Name of the headline omnibus column; must be in
#'   \code{test_columns}.
#' @param families Optional named list mapping family label -> character vector
#'   of test columns to draw on that panel. If \code{NULL}, the per-family PDF
#'   is skipped. If supplied alongside \code{test_columns}, the two sets must
#'   match.
#' @param subtitle Optional string appended to the omni title and drawn as outer
#'   mtext on the family grid.
#' @param rds_basename File-stem for the lambdas RDS (default
#'   \code{".family_lambdas"}; final name \code{"<basename>_<label>.rds"}).
#' @param pdf_basename_omni File-stem for the qq_omni PDF. Default
#'   \code{"qq_omni"}.
#' @param pdf_basename_family File-stem for the qq_by_family PDF. Default
#'   \code{"qq_by_family"}.
#' @param family_palette Optional character vector of colors for the family-grid
#'   points; \code{NULL} uses \code{hcl.colors(max_family_size, "Dark 3")}.
#' @param palette_pch Optional integer vector of plotting symbols cycled by
#'   within-family position. Default \code{c(16, 17, 15, 18, 8, 4)}.
#' @param family_title Optional \code{function(fam) -> character} for the panel
#'   title; \code{NULL} (default) uses \code{"<fam> family"}.
#' @param write_rds Logical. If \code{TRUE} (default), save the per-test lambdas
#'   to \code{"<rds_basename>_<label>.rds"}.
#' @param verbose Logical. If \code{TRUE} (default), cats a timestamped line per
#'   written file.
#'
#' @return Invisibly, a list with \code{lambda_omni}, \code{lambdas} (named per
#'   \code{test_columns}), \code{n} (\code{nrow(df)}), and \code{paths}
#'   (\code{qq_omni}; \code{qq_by_family}, which is \code{NA} when no
#'   \code{families} are supplied; and \code{family_lambdas_rds}, which is
#'   \code{NA} when \code{write_rds} is \code{FALSE}).
#'
#' @seealso \code{\link[GLOWr]{plot_qq}}, \code{\link[GLOWr]{compute_lambda_gc}}
#'
#' @importFrom GLOWr plot_qq compute_lambda_gc
#' @importFrom grDevices pdf dev.off hcl.colors
#' @importFrom graphics par plot abline points legend mtext
#' @importFrom stats setNames
#' @export
summarize_qq_and_lambdas <- function(
    df, label, output_dir,
    test_columns = NULL,
    omni_test,
    families = NULL,
    subtitle = NULL,
    rds_basename = ".family_lambdas",
    pdf_basename_omni = "qq_omni",
    pdf_basename_family = "qq_by_family",
    family_palette = NULL,
    palette_pch = c(16, 17, 15, 18, 8, 4),
    family_title = NULL,
    write_rds = TRUE,
    verbose = TRUE) {

  if (!is.data.frame(df))
    stop("`df` must be a data.frame")
  if (!is.character(label) || length(label) != 1L || !nzchar(label))
    stop("`label` must be a single non-empty string")
  if (!is.character(omni_test) || length(omni_test) != 1L)
    stop("`omni_test` must be a single string")
  if (!is.null(family_title) && !is.function(family_title))
    stop("`family_title` must be a function(fam) -> character, or NULL")

  if (is.null(test_columns)) {
    if (is.null(families))
      stop("Must supply `test_columns` or `families` (or both)")
    test_columns <- unlist(families, use.names = FALSE)
  } else if (!is.null(families)) {
    fam_tests <- unlist(families, use.names = FALSE)
    if (!setequal(fam_tests, test_columns)) {
      missing_in_tc  <- setdiff(fam_tests,    test_columns)
      missing_in_fam <- setdiff(test_columns, fam_tests)
      stop("`families` and `test_columns` must contain the same tests.\n",
           if (length(missing_in_tc))
             paste0("  In families but not test_columns: ",
                    paste(missing_in_tc, collapse = ", "), "\n"),
           if (length(missing_in_fam))
             paste0("  In test_columns but not families: ",
                    paste(missing_in_fam, collapse = ", ")))
    }
  }

  missing_cols <- setdiff(c(test_columns, omni_test), names(df))
  if (length(missing_cols))
    stop("Columns not found in `df`: ", paste(missing_cols, collapse = ", "))
  if (!omni_test %in% test_columns)
    stop("`omni_test` ('", omni_test, "') must be one of `test_columns`")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  n <- nrow(df)
  tag <- function(stem, ext)
    file.path(output_dir, sprintf("%s_%s%s", stem, label, ext))

  # Per-test lambda_GC. NA flows through (e.g., a column with fewer than
  # 50 non-NA p-values; compute_lambda_gc returns NA in that case).
  lambdas <- setNames(
    lapply(test_columns, function(t) compute_lambda_gc(df[[t]])),
    test_columns
  )
  lambda_omni <- lambdas[[omni_test]]

  # ---- 1. Omnibus QQ ----
  omni_path <- tag(pdf_basename_omni, ".pdf")
  if (verbose) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                            basename(omni_path)))

  title_omni <- sprintf("QQ %s  (lambda_GC = %.2f%s)",
                         omni_test, lambda_omni,
                         if (!is.null(subtitle)) sprintf("; %s", subtitle)
                         else "")
  pdf(omni_path, width = 6, height = 6)
  plot_qq(df, pvalue_col = omni_test, title = title_omni)
  dev.off()

  # ---- 2. Family grid (optional) ----
  fam_path <- NA_character_
  if (!is.null(families)) {
    fam_path <- tag(pdf_basename_family, ".pdf")
    if (verbose) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"),
                              basename(fam_path)))

    .logp_max <- function(v) {
      v <- v[is.finite(v) & v > 0 & v <= 1]
      if (length(v) == 0L) 0 else max(-log10(v))
    }
    y_max <- max(vapply(test_columns,
                        function(t) .logp_max(df[[t]]), numeric(1L)))
    x_max <- -log10(0.5 / n)

    n_fam  <- length(families)
    n_cols <- ceiling(sqrt(n_fam))
    n_rows <- ceiling(n_fam / n_cols)

    max_fam_size <- max(lengths(families))
    fam_palette  <- family_palette %||% hcl.colors(max_fam_size,
                                                    palette = "Dark 3")

    pdf(fam_path, width = 5 * n_cols, height = 5 * n_rows)
    op <- par(mfrow = c(n_rows, n_cols),
              mar  = c(4.0, 4.0, 3.5, 0.5),
              oma  = c(0, 0, 1.5, 0))
    for (fam in names(families)) {
      tests <- families[[fam]]
      # Cycle colour and pch by within-family position; recycle pch if a
      # family has more tests than supplied glyphs.
      pch_cycle    <- palette_pch[((seq_along(tests) - 1L) %%
                                     length(palette_pch)) + 1L]
      test_palette <- setNames(fam_palette[seq_along(tests)], tests)
      test_pch     <- setNames(pch_cycle, tests)

      plot(NA, NA,
           xlim = c(0, x_max), ylim = c(0, y_max),
           xlab = "Expected -log10(p)",
           ylab = "Observed -log10(p)",
           main = if (is.null(family_title)) sprintf("%s family", fam)
                  else family_title(fam))
      abline(0, 1, lty = 2, col = "grey60")

      for (t in tests) {
        v <- df[[t]]
        v <- v[is.finite(v) & v > 0 & v <= 1]
        if (length(v) == 0L) next
        obs_q  <- -log10(sort(v))
        this_e <- -log10((seq_len(length(obs_q)) - 0.5) / length(obs_q))
        points(this_e, obs_q,
               pch = test_pch[t], cex = 0.5, col = test_palette[t])
      }
      legend("topleft",
             legend = sprintf("%s  (lambda = %.2f)",
                              tests, unlist(lambdas[tests])),
             col = unlist(test_palette),
             pch = unlist(test_pch),
             bty = "n", cex = 0.80)
    }
    if (!is.null(subtitle))
      mtext(subtitle, side = 3, outer = TRUE, line = 0.2, cex = 0.95)
    par(op)
    dev.off()
  }

  # ---- 3. Save RDS (optional; 06-run-diagnostics consumes it where on) ----
  rds_path <- NA_character_
  if (write_rds) {
    rds_path <- tag(rds_basename, ".rds")
    saveRDS(lambdas, rds_path)
  }

  invisible(list(
    lambda_omni = lambda_omni,
    lambdas     = lambdas,
    n           = n,
    paths       = list(qq_omni            = omni_path,
                       qq_by_family       = fam_path,
                       family_lambdas_rds = rds_path)
  ))
}
