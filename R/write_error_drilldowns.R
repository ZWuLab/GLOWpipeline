# ==============================================================================
# Per-Error CSV Drill-down Writer
# ==============================================================================
#
# Pure, resolver-agnostic per-error CSV drill-down writer. For each error row it
# asks a caller-supplied resolver where the analyzed region(s) are, re-opens the
# chr GDS, runs GLOWr::extract_variant_set(), and writes per-variant cohort
# genotype summaries; on any non-"ok" resolver status it writes a note-only CSV
# that explains itself. Shared by the gene-based and sliding-window diagnostics
# stages.
#
# EXPORTED FUNCTIONS:
#   - write_error_drilldowns()   Resolver-agnostic per-error CSV drill-down writer


#################### EXPORTED MAIN FUNCTIONS ####################

#' Write Per-Error Variant Drill-down CSVs
#'
#' Pure, resolver-agnostic per-error CSV drill-down writer. For each row of
#' \code{err_rows} it asks \code{resolve_regions(err_row)} where the analyzed
#' region(s) are; on \code{status = "ok"} it re-opens the chr GDS, runs
#' \code{\link[GLOWr]{extract_variant_set}}, and writes
#' \code{<err_dir>/<label>-variants.csv} with per-variant cohort genotype
#' summaries (\code{n_g0}, \code{n_g1}, \code{n_g2}, \code{n_missing},
#' variance, \code{is_constant}); on any non-"ok" status it writes a note-only
#' CSV carrying a \code{status} code plus a plain reason (so an empty CSV
#' explains itself).
#'
#' @details
#' Two scripts call this helper with different resolvers: the gene-based
#' diagnostics resolves by matching the error's label/chr against the aggregated
#' results; the window-based diagnostics parses the analyzed region directly out
#' of the window label \code{"chrN:start-end"}. Multiple regions for the same
#' error are emitted as \code{<label>-<chr>_<start>_<end>-variants.csv} (no
#' silent overwrite).
#'
#' @param err_rows data.frame with at least \code{chr}, \code{label}. Optional
#'   \code{source}, \code{message} (carried into note CSVs when present).
#' @param resolve_regions \code{function(err_row) -> list}. Required key
#'   \code{$status}; \code{"ok"} proceeds to GDS re-extraction, any other value
#'   writes a note CSV with that status code. On \code{"ok"}, also required:
#'   \code{$regions}, a list of lists \code{{chr, start, end, label}}. On non-ok
#'   statuses, optional \code{$note} (note text) and \code{$location} (copied
#'   into the returned \code{location} column).
#' @param null_model The run's fitted null-model list (must include \code{$sample_id}).
#' @param filter A \code{glow_variant_filter} (from \code{variant_filter()}).
#' @param gds_dir Character. Directory holding the per-chr GDS files.
#' @param gds_pattern Character. Per-chr GDS filename pattern with a literal
#'   \code{\{chr\}} placeholder; the path is
#'   \code{file.path(gds_dir, gsub("\{chr\}", err_row$chr, gds_pattern))}.
#' @param pi_features Character vector of annotation feature names passed to
#'   \code{extract_variant_set()}.
#' @param err_dir Output directory for the per-error CSVs (created if missing).
#' @param file_basename_fn Optional \code{function(err_row) -> CSV basename}.
#'   Default \code{sprintf("\%s-variants.csv", err_row$label)}.
#' @param link_path_prefix Character. Prefix for the markdown links written into
#'   \code{err_rows$error_diagnostics}. Default \code{"error_diagnostics"}.
#' @param verbose Logical. If \code{TRUE} (default), cats one timestamped line
#'   per error processed.
#'
#' @return \code{err_rows} with two appended columns: \code{location}
#'   ("chr:start-end" of the analyzed region(s), semicolon-separated when
#'   multiple) and \code{error_diagnostics} (markdown link(s) to the written
#'   CSV(s)).
#'
#' @importFrom GLOWr extract_variant_set
#' @importFrom SeqArray seqOpen seqClose
#' @importFrom stats var na.omit
#' @export
write_error_drilldowns <- function(err_rows, resolve_regions,
                                    null_model, filter,
                                    gds_dir, gds_pattern, pi_features,
                                    err_dir,
                                    file_basename_fn = NULL,
                                    link_path_prefix = "error_diagnostics",
                                    verbose = TRUE) {

  required <- c("chr", "label")
  miss <- setdiff(required, names(err_rows))
  if (length(miss) > 0L)
    stop("err_rows is missing required column(s): ",
         paste(miss, collapse = ", "))
  if (!is.function(resolve_regions))
    stop("`resolve_regions` must be a function(err_row) -> list")

  if (is.null(file_basename_fn))
    file_basename_fn <- function(r) sprintf("%s-variants.csv", r$label)

  dir.create(err_dir, recursive = TRUE, showWarnings = FALSE)

  # Note-only CSV: status + plain reason so an empty file explains itself.
  write_note <- function(path, status, note) {
    utils::write.csv(data.frame(status = status, note = note,
                                stringsAsFactors = FALSE),
                     path, row.names = FALSE)
  }

  if (is.null(err_rows$source))  err_rows$source  <- NA_character_
  if (is.null(err_rows$message)) err_rows$message <- NA_character_
  err_rows$location          <- NA_character_
  err_rows$error_diagnostics <- NA_character_

  log_one <- function(fmt, ...) if (isTRUE(verbose))
    cat(sprintf("[%s]   %s\n",
                format(Sys.time(), "%H:%M:%S"),
                sprintf(fmt, ...)))

  default_note <- function(status, r)
    sprintf("Resolver returned status %s for %s (chr %s). source=%s; log message: %s",
            status, r$label, r$chr,
            if (is.na(r$source))  "NA" else r$source,
            if (is.na(r$message)) "NA" else r$message)

  for (k in seq_len(nrow(err_rows))) {
    r <- err_rows[k, ]
    base    <- file_basename_fn(r)
    out_csv <- file.path(err_dir, base)
    err_rows$error_diagnostics[k] <-
      sprintf("[`%s`](%s/%s)", base, link_path_prefix, base)

    res <- resolve_regions(r)
    if (!is.list(res) || is.null(res$status))
      stop("`resolve_regions` must return a list with $status (and $regions when status='ok')")

    if (!identical(res$status, "ok")) {
      write_note(out_csv, res$status,
                 if (!is.null(res$note)) res$note
                 else default_note(res$status, r))
      if (!is.null(res$location))
        err_rows$location[k] <- res$location
      log_one("%s: %s", r$label, res$status)
      next
    }

    regions <- res$regions
    if (!is.list(regions) || length(regions) == 0L)
      stop("`resolve_regions` returned status='ok' but no $regions for ",
           r$label)

    err_rows$location[k] <- paste(
      vapply(regions, function(rg) sprintf("%s:%d-%d",
                                            rg$chr, rg$start, rg$end),
             character(1L)), collapse = "; ")

    gds_path <- file.path(gds_dir, gsub("\\{chr\\}", r$chr, gds_pattern))
    if (!file.exists(gds_path)) {
      warning("GDS not found for chr", r$chr, ": ", gds_path)
      write_note(out_csv, "GDS_NOT_FOUND",
                 sprintf("GDS file missing: %s", gds_path))
      next
    }

    multi <- length(regions) > 1L
    if (multi) {
      files <- vapply(regions, function(rg)
        sprintf("%s-%s_%d_%d-variants.csv",
                r$label, rg$chr, rg$start, rg$end),
        character(1L))
      err_rows$error_diagnostics[k] <- paste(
        sprintf("[`%s`](%s/%s)",
                files, link_path_prefix, files),
        collapse = "<br>")
    }

    for (rg in regions) {
      this_csv <- if (multi)
        file.path(err_dir, sprintf("%s-%s_%d_%d-variants.csv",
                                    r$label, rg$chr, rg$start, rg$end))
      else out_csv

      gds  <- seqOpen(gds_path, readonly = TRUE)
      vset <- tryCatch(
        extract_variant_set(gds, rg, filter,
                            sample_id        = null_model$sample_id,
                            annotation_names = pi_features,
                            verbose          = 0),
        error = function(e) list(.extract_error = conditionMessage(e)))
      seqClose(gds)

      if (is.null(vset)) {
        write_note(this_csv, "UNEXPECTED_EMPTY",
                   sprintf("Region %s:%d-%d re-extracted to NULL (no qualifying variants now). Likely a region-lookup or code-drift bug -- investigate; do NOT read as an empty unit. source=%s; log message: %s",
                           rg$chr, rg$start, rg$end,
                           if (is.na(r$source))  "NA" else r$source,
                           if (is.na(r$message)) "NA" else r$message))
        warning("UNEXPECTED_EMPTY: ", r$label, " (", rg$chr, ":",
                rg$start, "-", rg$end, ") re-extracted to NULL")
        log_one("%s: UNEXPECTED NULL (%s:%d-%d)",
                r$label, rg$chr, rg$start, rg$end)
        next
      }
      if (!is.null(vset$.extract_error)) {
        write_note(this_csv, "EXTRACT_ERROR",
                   sprintf("extract_variant_set() errored: %s. source=%s; log message: %s",
                           vset$.extract_error,
                           if (is.na(r$source))  "NA" else r$source,
                           if (is.na(r$message)) "NA" else r$message))
        log_one("%s: extract_variant_set errored", r$label)
        next
      }

      G <- vset$G
      n_g0 <- colSums(G == 0L, na.rm = TRUE)
      n_g1 <- colSums(G == 1L, na.rm = TRUE)
      n_g2 <- colSums(G == 2L, na.rm = TRUE)
      n_miss_geno <- colSums(is.na(G))
      var_in_coh  <- apply(G, 2L, function(v) var(v, na.rm = TRUE))
      is_const <- apply(G, 2L, function(v) {
        u <- unique(stats::na.omit(v)); length(u) <= 1L
      })
      const_val <- vapply(seq_len(ncol(G)), function(j) {
        u <- unique(stats::na.omit(G[, j]))
        if (length(u) == 1L) u else NA_real_
      }, numeric(1L))

      diag_df <- cbind(
        vset$variant_info,
        data.frame(
          n_g0_homo_alt         = as.integer(n_g0),
          n_g1_het              = as.integer(n_g1),
          n_g2_homo_ref         = as.integer(n_g2),
          n_missing_genotypes   = as.integer(n_miss_geno),
          variance_in_cohort    = as.numeric(var_in_coh),
          is_constant_in_cohort = as.logical(is_const),
          constant_value        = const_val,
          stringsAsFactors      = FALSE))
      utils::write.csv(diag_df, this_csv, row.names = FALSE)

      log_one("%s (%s:%d-%d): %d variants, %d constant in cohort",
              r$label, rg$chr, rg$start, rg$end,
              nrow(diag_df), sum(is_const))
    }
  }

  err_rows
}
