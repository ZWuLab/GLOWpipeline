# ==============================================================================
# Genome-Wide Aggregation of Per-Chunk Window-Scan Results
# ==============================================================================
#
# Concatenates the per-chunk flat tables in genomic order, tags genome-wide
# significant windows at alpha / n_windows, and merges overlapping/abutting
# significant windows into loci by interval-union. This is the counterpart of
# STAARpipeline's Sliding_Window_Results_Summary() (genome-wide row-bind +
# significant-results path).
#
# This function only reads a run's per-chunk results and returns tables; the
# run's summary scripts write them to disk. The locus-merge runs on the
# assembled genome-wide table, so loci that straddle chunk boundaries merge
# correctly.
#
# EXPORTED FUNCTIONS:
#   - aggregate_scan_results()    Genome rbind + threshold + locus-merge
#
# INTERNAL HELPER FUNCTIONS:
#   - .read_chunk_table()         Read one chunk file (csv | fst), chr as char
#   - .merge_windows_to_loci()    Interval-union of significant windows
#   - .pvalue_cols()              The test p-value columns (non-id/count/flag)
#   - .empty_loci()               0-row loci table with the correct schema


#################### EXPORTED FUNCTIONS ####################

#' Aggregate Per-Chunk Window-Scan Results Genome-Wide
#'
#' Reads the per-chunk flat tables a window scan produced, concatenates them in
#' genomic order, applies the genome-wide significance threshold
#' \code{alpha / n_windows}, and merges contiguous/overlapping significant
#' windows into loci.
#'
#' @details
#' \strong{Multiple-testing denominator.} \code{n_windows} is the number of
#' tested windows after the optional \code{cmac_cutoff} filter; the threshold is
#' \code{alpha / n_windows}, matching the STAARpipeline convention. Overlapping
#' windows are correlated, so this Bonferroni denominator is conservative — a
#' deliberate convention, not a defect.
#'
#' \strong{Locus-merging.} Significant windows are sorted by \code{(chr, start)};
#' a window joins the running locus when it is on the same chromosome and its
#' \code{start} is \code{<= running_end + merge_gap}, else it opens a new locus.
#' Each locus reports its position union, the number of member windows, and the
#' \emph{peak} (minimum \code{primary_test} p-value) member's label and full
#' p-value row (prefixed \code{peak_}). Running on the assembled genome-wide
#' table makes chunk-boundary loci merge correctly.
#'
#' @param chunk_files Character vector of per-chunk flat-table file paths.
#' @param primary_test Column name driving significance and locus peaks
#'   (default \code{"GLOW_Omni"}).
#' @param alpha Family-wise error rate (default \code{0.05}).
#' @param merge_gap Integer bp (default \code{0L}); windows within this gap of a
#'   locus merge into it. \code{0} merges only directly overlapping/abutting
#'   windows.
#' @param cmac_cutoff Optional numeric; windows with \code{cMAC <= cmac_cutoff}
#'   are dropped before thresholding (STAARpipeline-style). \code{NULL} (default)
#'   keeps all windows.
#' @param output_format Format of the \emph{input} chunk files: \code{"csv"}
#'   (default) or \code{"fst"} (requires the \pkg{fst} package).
#'
#' @return A list with:
#' \describe{
#'   \item{scan_results_all}{The genome-wide flat table (tested windows, ordered
#'     by chr then start).}
#'   \item{scan_sig_windows}{The raw significant windows
#'     (\code{primary_test < alpha / n_windows}); kept so window-to-window
#'     comparisons stay possible.}
#'   \item{scan_loci}{Merged significant loci (see Details for columns).}
#'   \item{summary}{A 1-row \code{data.frame}: \code{n_windows}, \code{threshold},
#'     \code{n_sig_windows}, \code{n_loci}, \code{primary_test}, \code{alpha},
#'     \code{n_dropped_cmac}.}
#' }
#'
#' @seealso [glow_scan_chunk()], [define_genome_chunks()]
#'
#' @examples
#' \dontrun{
#' files <- list.files("outputs/results", "^scan_chunk_.*\\.csv$",
#'                     full.names = TRUE)
#' agg <- aggregate_scan_results(files, primary_test = "GLOW_Omni", alpha = 0.05)
#' agg$summary
#' head(agg$scan_loci)
#' }
#'
#' @export
aggregate_scan_results <- function(chunk_files,
                                   primary_test = "GLOW_Omni",
                                   alpha = 0.05,
                                   merge_gap = 0L,
                                   cmac_cutoff = NULL,
                                   output_format = "csv") {

  if (length(chunk_files) == 0L) {
    stop("`chunk_files` is empty: no per-chunk tables to aggregate.",
         call. = FALSE)
  }
  if (!output_format %in% c("csv", "fst")) {
    stop("`output_format` must be \"csv\" or \"fst\" (got \"", output_format,
         "\").", call. = FALSE)
  }
  if (identical(output_format, "fst") &&
      !requireNamespace("fst", quietly = TRUE)) {
    stop("output_format = \"fst\" requires the fst package. Install fst or use ",
         "\"csv\".", call. = FALSE)
  }

  # ---- 1. Read + rbind ----
  tables <- lapply(chunk_files, .read_chunk_table, output_format = output_format)
  all <- do.call(rbind, tables)
  rownames(all) <- NULL

  if (!primary_test %in% names(all)) {
    stop("`primary_test` column \"", primary_test, "\" is not in the chunk ",
         "tables. Available test columns: ",
         paste(.pvalue_cols(all), collapse = ", "), ".", call. = FALSE)
  }

  # ---- 2. Order genomically; optional cMAC filter ----
  all$chr <- as.character(all$chr)
  all <- all[order(.chr_rank(all$chr), all$start), , drop = FALSE]
  rownames(all) <- NULL

  n_dropped_cmac <- 0L
  if (!is.null(cmac_cutoff)) {
    if (!"cMAC" %in% names(all)) {
      stop("`cmac_cutoff` was supplied but the tables have no `cMAC` column.",
           call. = FALSE)
    }
    keep <- !is.na(all$cMAC) & all$cMAC > cmac_cutoff
    n_dropped_cmac <- sum(!keep)
    all <- all[keep, , drop = FALSE]
    rownames(all) <- NULL
  }

  # ---- 3-5. Threshold ----
  n_windows <- nrow(all)
  threshold <- if (n_windows > 0L) alpha / n_windows else NA_real_
  p <- all[[primary_test]]
  sig_idx <- !is.na(p) & p < threshold
  scan_sig_windows <- all[sig_idx, , drop = FALSE]
  rownames(scan_sig_windows) <- NULL

  # ---- 6. Merge significant windows into loci ----
  scan_loci <- .merge_windows_to_loci(scan_sig_windows, merge_gap, primary_test)

  summary <- data.frame(
    n_windows      = n_windows,
    threshold      = threshold,
    n_sig_windows  = nrow(scan_sig_windows),
    n_loci         = nrow(scan_loci),
    primary_test   = primary_test,
    alpha          = alpha,
    n_dropped_cmac = n_dropped_cmac,
    stringsAsFactors = FALSE
  )

  list(scan_results_all = all,
       scan_sig_windows = scan_sig_windows,
       scan_loci        = scan_loci,
       summary          = summary)
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Read one per-chunk flat table
#'
#' @param path File path.
#' @param output_format \code{"csv"} or \code{"fst"}.
#' @return A data.frame with \code{chr} coerced to character (so genomic
#'   ordering and locus-merge are not confused by numeric chromosome reads).
#' @keywords internal
#' @noRd
.read_chunk_table <- function(path, output_format) {
  if (identical(output_format, "fst")) {
    df <- fst::read_fst(path)
  } else {
    df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                          colClasses = c(chr = "character"))
  }
  df$chr <- as.character(df$chr)
  df
}


#' The test p-value columns of a flat scan table
#'
#' Every column that is not a window-identity, variant-count, or flag column —
#' i.e. the GLOW and STAAR p-value columns.
#'
#' @param df A flat scan table.
#' @return Character vector of p-value column names.
#' @keywords internal
#' @noRd
.pvalue_cols <- function(df) {
  non_pval <- c("label", "chr", "start", "end", "midpoint",
                "n_variants", "n_after_filter", "n_after_ld",
                "n_after_collapse", "cMAC", "is_single")
  setdiff(names(df), non_pval)
}


#' Merge significant windows into loci (interval-union)
#'
#' @param sig Significant windows (full flat schema), or 0-row.
#' @param merge_gap Integer bp join gap.
#' @param primary_test The test column whose minimum marks the locus peak.
#' @return A loci data.frame (see [aggregate_scan_results()]).
#' @keywords internal
#' @noRd
.merge_windows_to_loci <- function(sig, merge_gap, primary_test) {
  pcols <- .pvalue_cols(sig)
  if (nrow(sig) == 0L) return(.empty_loci(pcols, primary_test))

  sig <- sig[order(.chr_rank(sig$chr), sig$start), , drop = FALSE]

  make_locus <- function(members) {
    peak_i   <- which.min(members[[primary_test]])
    peak_row <- members[peak_i, pcols, drop = FALSE]
    names(peak_row) <- paste0("peak_", names(peak_row))
    cbind(
      data.frame(
        chr         = as.character(members$chr[1]),
        locus_start = as.integer(min(members$start)),
        locus_end   = as.integer(max(members$end)),
        n_windows   = nrow(members),
        peak_label  = as.character(members$label[peak_i]),
        stringsAsFactors = FALSE
      ),
      peak_row
    )
  }

  loci   <- list()
  open   <- 1L                      # row index that opened the running locus
  cur    <- c(open)                 # member row indices of the running locus
  cur_chr <- sig$chr[1]
  cur_end <- sig$end[1]

  for (k in seq_len(nrow(sig))[-1]) {
    if (identical(sig$chr[k], cur_chr) &&
        sig$start[k] <= cur_end + merge_gap) {
      cur     <- c(cur, k)
      cur_end <- max(cur_end, sig$end[k])
    } else {
      loci[[length(loci) + 1L]] <- make_locus(sig[cur, , drop = FALSE])
      cur     <- k
      cur_chr <- sig$chr[k]
      cur_end <- sig$end[k]
    }
  }
  loci[[length(loci) + 1L]] <- make_locus(sig[cur, , drop = FALSE])

  out <- do.call(rbind, loci)
  out <- cbind(locus_id = seq_len(nrow(out)), out)
  rownames(out) <- NULL
  out
}


#' A 0-row loci table with the correct schema
#'
#' @param pcols The p-value column names (to build the \code{peak_} block).
#' @param primary_test Unused beyond documentation symmetry.
#' @return A 0-row loci data.frame.
#' @keywords internal
#' @noRd
.empty_loci <- function(pcols, primary_test) {
  base <- data.frame(
    locus_id    = integer(0),
    chr         = character(0),
    locus_start = integer(0),
    locus_end   = integer(0),
    n_windows   = integer(0),
    peak_label  = character(0),
    stringsAsFactors = FALSE
  )
  if (length(pcols) == 0L) return(base)
  peak_block <- as.data.frame(
    stats::setNames(rep(list(numeric(0)), length(pcols)),
                    paste0("peak_", pcols)),
    stringsAsFactors = FALSE
  )
  cbind(base, peak_block)
}
