# ==============================================================================
# Genome Chunking for the SLURM Window-Scan Array
# ==============================================================================
#
# Splits a chromosome set into fixed-Mb chunks, one per compute array task, so
# the window scan saturates many CPUs and bounds per-task wall-time. Returns a
# lookup table the scan indexes by array-task id — arithmetic-free at run time.
#
# BOUNDARY RULE ("assign by window start"; refined here so the
# load-bearing invariant holds by construction):
#   A sliding window is OWNED by the chunk that contains its START. Because
#   sliding windows overlap (window_size > step_size), a window owned by a chunk
#   may END past chunk_end — that is fine; ownership is by start only, so no
#   window is shared between chunks and none is dropped.
#
#   For the union of all chunks' windows to EQUAL a single full-span slide (the
#   "no duplicates, no gaps" guarantee), every chunk must generate windows
#   on the SAME global step grid anchored at the chromosome's span start. We
#   therefore make each chunk own a contiguous run of `windows_per_chunk =
#   floor(chunk_size_bp / step_size)` grid positions, so chunk_start is always a
#   grid-aligned window start. The scan then generates a chunk's
#   windows with
#       define_regions_window(chr, chunk_start, span_end, window_size, step_size)
#   and keeps the rows with start <= chunk_end. (Generating to `span_end`, not
#   `chunk_end`, is what lets a window start inside the chunk and spill past it,
#   exactly as the full-span slide would — generating only to `chunk_end` would
#   silently drop every window whose start lies within `window_size` of a chunk
#   edge, leaving ~window_size coverage gaps at every boundary.) `span_end` is
#   carried in the table so the scan needs only this one lookup.
#
#   For all realistic configs (step_size divides chunk_size_bp, e.g. 5 Mb / 10
#   kb), chunk_start / chunk_end equal the naive fixed-Mb boundaries exactly;
#   only when chunk_size_bp is not a whole number of steps are chunks rounded
#   down to floor(chunk_size_bp/step)*step bp to keep the grid aligned.
#
# EXPORTED FUNCTIONS:
#   - define_genome_chunks()      Build the per-chunk lookup table
#
# INTERNAL HELPER FUNCTIONS:
#   - .chr_rank()                 Genomic-order rank for a chromosome label
#   - .n_windows_in_span()        Full-span window count (define_regions.R:127)


#################### EXPORTED FUNCTIONS ####################

#' Split a Chromosome Set into Window-Scan Chunks
#'
#' Builds the chunk lookup table for the genomic-chunk SLURM array that drives a
#' sliding-window GLOW scan. Each chunk is one array task; the table maps
#' \code{SLURM_ARRAY_TASK_ID} to a chromosome and a contiguous run of window
#' start positions.
#'
#' @details
#' \strong{Boundary rule.} A window is owned by the chunk containing its
#' \emph{start}. Chunks own contiguous runs of
#' \code{windows_per_chunk = floor(chunk_size_mb * 1e6 / step_size)} positions on
#' the per-chromosome step grid (anchored at the span start), so \code{chunk_start}
#' is always a grid-aligned window start and the union of all chunks' windows
#' equals one full-span slide — no duplicates, no gaps. The scan builds
#' a chunk's windows with
#' \code{define_regions_window(chr, chunk_start, span_end, window_size, step_size)}
#' and keeps rows with \code{start <= chunk_end}; generating to \code{span_end}
#' (not \code{chunk_end}) lets a window start in the chunk and spill past its end,
#' matching the full-span slide. \code{span_end} is returned for that one lookup.
#'
#' For configs where \code{step_size} divides \code{chunk_size_mb * 1e6} (the
#' usual case, e.g. 5 Mb / 10 kb), \code{chunk_start}/\code{chunk_end} equal the
#' naive fixed-Mb boundaries; otherwise chunks are rounded down to a whole number
#' of steps to keep the grid aligned.
#'
#' @param chrom_spans A \code{data.frame} with columns \code{chr}, \code{start},
#'   \code{end} — one row per chromosome, giving its inclusive position extent
#'   (typically the per-chr GDS min/max position, computed during preparation). Rows
#'   are processed in genomic order regardless of input order.
#' @param window_size Integer window width in bp (e.g. \code{100000L}).
#' @param step_size Integer step between consecutive window starts (e.g.
#'   \code{10000L}). Governs the window count.
#' @param chunk_strategy Character. Only \code{"fixed_mb"} is supported in v1;
#'   any other value is an error.
#' @param chunk_size_mb Numeric chunk width in Mb (default \code{5}, the
#'   STAARpipeline-proven value). Must be large enough to hold at least one step
#'   (\code{chunk_size_mb * 1e6 >= step_size}).
#' @param gds_paths Reserved for a future variant-count-balancing strategy;
#'   unused in v1 (\code{"fixed_mb"}).
#'
#' @return A \code{data.frame}, one row per chunk in genomic order, with columns:
#' \describe{
#'   \item{chunk_id}{Integer, 1-based global id; the scan selects its
#'     chunk by \code{SLURM_ARRAY_TASK_ID == chunk_id}.}
#'   \item{chr}{Chromosome label.}
#'   \item{chunk_start}{Integer, inclusive lower bound of owned window starts
#'     (grid-aligned). Passed as \code{define_regions_window}'s \code{start}.}
#'   \item{chunk_end}{Integer, inclusive upper bound of owned window starts; the
#'     scan keeps windows with \code{start <= chunk_end}.}
#'   \item{span_end}{Integer, the chromosome's max position; passed as
#'     \code{define_regions_window}'s \code{end} so boundary-spanning windows are
#'     generated and window ends are capped at the chromosome extent.}
#'   \item{est_n_windows}{Integer, the number of window starts owned by the chunk
#'     (exact for fixed geometry; for SLURM load inspection).}
#' }
#'
#' @seealso \code{GLOWr::define_regions_window}, [glow_scan_chunk()]
#'
#' @examples
#' spans <- data.frame(chr = c("21", "22"),
#'                     start = c(1L, 1L),
#'                     end = c(12000000L, 8000000L))
#' chunks <- define_genome_chunks(spans, window_size = 100000L,
#'                                step_size = 10000L, chunk_size_mb = 5)
#' chunks
#'
#' @export
define_genome_chunks <- function(chrom_spans,
                                  window_size, step_size,
                                  chunk_strategy = "fixed_mb",
                                  chunk_size_mb = 5,
                                  gds_paths = NULL) {

  # ---- Validate ----
  if (!identical(chunk_strategy, "fixed_mb")) {
    stop("Only chunk_strategy = \"fixed_mb\" is supported in v1 (got \"",
         chunk_strategy, "\"). Variant-count balancing is a future option.",
         call. = FALSE)
  }
  req <- c("chr", "start", "end")
  if (!is.data.frame(chrom_spans) || !all(req %in% names(chrom_spans))) {
    stop("`chrom_spans` must be a data.frame with columns: ",
         paste(req, collapse = ", "), ".", call. = FALSE)
  }
  if (nrow(chrom_spans) == 0L) {
    stop("`chrom_spans` has no rows.", call. = FALSE)
  }
  window_size <- as.integer(window_size)
  step_size   <- as.integer(step_size)
  if (is.na(window_size) || window_size <= 0L)
    stop("`window_size` must be a positive integer.", call. = FALSE)
  if (is.na(step_size) || step_size <= 0L)
    stop("`step_size` must be a positive integer.", call. = FALSE)
  if (!is.numeric(chunk_size_mb) || length(chunk_size_mb) != 1L ||
      is.na(chunk_size_mb) || chunk_size_mb <= 0)
    stop("`chunk_size_mb` must be a positive number.", call. = FALSE)

  chunk_size_bp <- chunk_size_mb * 1e6
  windows_per_chunk <- floor(chunk_size_bp / step_size)
  if (windows_per_chunk < 1) {
    stop("chunk_size_mb * 1e6 (", format(chunk_size_bp, scientific = FALSE),
         " bp) is smaller than step_size (", step_size,
         " bp): a chunk must hold at least one window start. ",
         "Increase chunk_size_mb.", call. = FALSE)
  }
  windows_per_chunk <- as.integer(windows_per_chunk)

  # ---- Order chromosomes genomically ----
  ord <- order(.chr_rank(chrom_spans$chr), chrom_spans$start)
  spans <- chrom_spans[ord, , drop = FALSE]

  # ---- Build chunks per chromosome ----
  per_chr <- lapply(seq_len(nrow(spans)), function(i) {
    chr        <- as.character(spans$chr[i])
    span_start <- as.integer(spans$start[i])
    span_end   <- as.integer(spans$end[i])
    if (is.na(span_start) || is.na(span_end) || span_end < span_start) {
      stop("chrom_spans row ", i, " (chr ", chr,
           ") has invalid extent [", span_start, ", ", span_end, "].",
           call. = FALSE)
    }

    N <- .n_windows_in_span(span_start, span_end, window_size, step_size)
    n_chunks <- as.integer(ceiling(N / windows_per_chunk))

    j0 <- (seq_len(n_chunks) - 1L) * windows_per_chunk      # 0-based first idx
    j1 <- pmin(seq_len(n_chunks) * windows_per_chunk, N) - 1L # 0-based last idx

    chunk_start <- as.integer(span_start + j0 * step_size)
    # Owned-start upper bound: one below the next chunk's first start; the last
    # chunk extends to the chromosome's span end.
    chunk_end <- as.integer(span_start + (j0 + windows_per_chunk) * step_size - 1L)
    chunk_end[n_chunks] <- span_end

    data.frame(
      chr           = rep(chr, n_chunks),
      chunk_start   = chunk_start,
      chunk_end     = chunk_end,
      span_end      = rep(span_end, n_chunks),
      est_n_windows = as.integer(j1 - j0 + 1L),
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, per_chr)
  out <- cbind(chunk_id = seq_len(nrow(out)), out)
  rownames(out) <- NULL
  out
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Genomic-order rank for a chromosome label
#'
#' Numeric chromosomes 1..22 sort first (ascending), then X, Y, M/MT, then any
#' other label alphabetically — so a global \code{chunk_id} runs in genomic
#' order regardless of how \code{chrom_spans} rows were supplied. Accepts labels
#' with or without a \code{"chr"} prefix.
#'
#' @param chr Character (or coercible) vector of chromosome labels.
#' @return Numeric rank vector (lower = earlier in the genome).
#' @keywords internal
#' @noRd
.chr_rank <- function(chr) {
  raw <- sub("^chr", "", as.character(chr), ignore.case = TRUE)
  num <- suppressWarnings(as.numeric(raw))
  special <- c(X = 100, Y = 101, M = 102, MT = 102)
  rank <- num                                   # numeric chrs keep their value
  is_special <- is.na(num) & toupper(raw) %in% names(special)
  rank[is_special] <- special[toupper(raw[is_special])]
  # Anything else (unknown contigs) sorts after, alphabetically and stably.
  is_other <- is.na(rank)
  if (any(is_other)) {
    rank[is_other] <- 1000 + rank(raw[is_other], ties.method = "min")
  }
  rank
}


#' Full-span window count for a chromosome
#'
#' Mirrors \code{GLOWr::define_regions_window}'s row-count formula
#' so \code{est_n_windows} matches what the scan's
#' \code{define_regions_window} call will produce.
#'
#' @param span_start,span_end Integer inclusive chromosome extent.
#' @param window_size,step_size Integer window geometry.
#' @return Integer window count (>= 1).
#' @keywords internal
#' @noRd
.n_windows_in_span <- function(span_start, span_end, window_size, step_size) {
  max(1L, as.integer(floor((span_end - span_start + 1L - window_size) /
                             step_size)) + 1L)
}
