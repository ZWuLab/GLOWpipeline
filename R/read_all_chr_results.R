# ==============================================================================
# Concatenate Per-Chromosome Result CSVs
# ==============================================================================
#
# Reads the 22 per-chr gene-level result CSVs and concatenates them in genomic
# order (chr 1..22, then ascending start), with a schema-consistency check.
#
# EXPORTED FUNCTIONS:
#   - read_all_chr_results()   Concat 22 per-chr result CSVs, genomic order


#################### EXPORTED MAIN FUNCTIONS ####################

#' Concatenate the Per-Chromosome Result CSVs in Genomic Order
#'
#' Concatenates the 22 per-chr gene-level result CSVs into a single data.frame
#' in genomic order (chr 1..22, then ascending start). Errors loudly if any chr
#' is missing or if column schemas differ across chrs.
#'
#' @details
#' The per-chr filenames are built as \code{sprintf(pattern, chroms)} under
#' \code{results_dir}. The default \code{chroms = 1:22} + \code{pattern}
#' reproduces the gene-module \code{glow_chr\{N\}.csv} genome-wide behaviour
#' exactly; pass \code{chroms} to aggregate only a subset (e.g. a single-chr
#' smoke run). After reading, every file's column schema is compared against the
#' first and the rows are returned ordered by \code{(as.integer(chr), start)}.
#'
#' @param results_dir Character. Directory holding the per-chr CSVs. Required.
#' @param chroms Integer/character vector of chromosomes to read. Default
#'   \code{1:22}. Each requested chromosome's file must exist (errors otherwise).
#' @param pattern Character. \code{sprintf} pattern (one \code{\%d}) for the
#'   per-chr filename. Default \code{"glow_chr\%d.csv"}.
#'
#' @return A data.frame of the requested per-chr rows concatenated and sorted
#'   into genomic order.
#'
#' @export
read_all_chr_results <- function(results_dir, chroms = 1:22,
                                 pattern = "glow_chr%d.csv") {
  files <- file.path(results_dir, sprintf(pattern, chroms))
  missing <- files[!file.exists(files)]
  if (length(missing) > 0L) {
    stop("Missing per-chr result CSVs:\n  ",
         paste(missing, collapse = "\n  "))
  }

  parts <- lapply(files, utils::read.csv, stringsAsFactors = FALSE)

  # Schema consistency check (every file vs the first requested one)
  cols1 <- colnames(parts[[1L]])
  for (i in seq_along(parts)[-1L]) {
    if (!identical(colnames(parts[[i]]), cols1)) {
      stop("Per-chr CSV column schemas differ. ",
           "the first file has ", length(cols1), " columns; file ", i,
           " has ", ncol(parts[[i]]), " columns.")
    }
  }

  all <- do.call(rbind, parts)
  # Genomic order
  all[order(as.integer(all$chr), all$start), , drop = FALSE]
}
