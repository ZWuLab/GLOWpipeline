# ==============================================================================
# Tests for define_genome_chunks() — the SLURM-array chunk lookup table
# ==============================================================================
#
# the UNION of every chunk's driver-generated window starts equals a single
# full-span slide — no duplicates, no gaps. That guards the boundary rule
# against silently dropping or double-counting windows at chunk edges.


# Reproduce the Stage-2 driver's per-chunk window generation: generate from
# chunk_start to span_end, keep starts <= chunk_end (design §4.3.1).
.driver_chunk_starts <- function(chunks, chr, window_size, step_size) {
  rows <- chunks[chunks$chr == chr, , drop = FALSE]
  starts <- integer(0)
  for (i in seq_len(nrow(rows))) {
    w <- GLOWr::define_regions_window(chr, rows$chunk_start[i], rows$span_end[i],
                                      window_size, step_size)
    w <- w[w$start <= rows$chunk_end[i], , drop = FALSE]
    starts <- c(starts, as.integer(w$start))
  }
  starts
}

.full_span_starts <- function(chr, span_start, span_end, window_size, step_size) {
  as.integer(
    GLOWr::define_regions_window(chr, span_start, span_end,
                                 window_size, step_size)$start)
}


test_that("single 12 Mb chr at 5 Mb chunks -> 3 chunks (5, 5, 2 Mb), capped", {
  spans <- data.frame(chr = "1", start = 1L, end = 12000000L)
  ch <- define_genome_chunks(spans, window_size = 100000L, step_size = 10000L,
                             chunk_size_mb = 5)

  expect_equal(nrow(ch), 3L)
  expect_equal(ch$chunk_id, 1:3)
  expect_equal(ch$chunk_start, c(1L, 5000001L, 10000001L))
  expect_equal(ch$chunk_end,   c(5000000L, 10000000L, 12000000L))
  expect_true(all(ch$span_end == 12000000L))          # ends carried for capping
  expect_equal(ch$est_n_windows, c(500L, 500L, 191L))
})


test_that("multiple chrs -> global, genomic-ordered chunk_id (regardless of input order)", {
  # Supplied out of genomic order (chr 2 before chr 1).
  spans <- data.frame(chr = c("2", "1"),
                      start = c(1L, 1L),
                      end = c(6000000L, 12000000L))
  ch <- define_genome_chunks(spans, window_size = 100000L, step_size = 10000L,
                             chunk_size_mb = 5)

  # chr 1 (3 chunks) sorts before chr 2 (2 chunks); chunk_id runs 1..5.
  expect_equal(ch$chunk_id, 1:5)
  expect_equal(ch$chr, c("1", "1", "1", "2", "2"))
  expect_equal(sum(ch$chr == "1"), 3L)
  expect_equal(sum(ch$chr == "2"), 2L)
})


test_that("BOUNDARY: union of chunk window starts == full-span slide (no dup, no gap)", {
  spans <- data.frame(chr = "1", start = 1L, end = 12000000L)
  win <- 100000L; step <- 10000L
  ch <- define_genome_chunks(spans, window_size = win, step_size = step,
                             chunk_size_mb = 5)

  got  <- .driver_chunk_starts(ch, "1", win, step)
  full <- .full_span_starts("1", 1L, 12000000L, win, step)

  expect_equal(anyDuplicated(got), 0L)                # no duplicates
  expect_equal(sort(got), sort(full))                 # no gaps, exact set
  expect_equal(length(got), length(full))
})


test_that("BOUNDARY holds when chunk size is NOT a whole number of steps", {
  # chunk_bp = 50000 is not a multiple of step 7000 -> wpc = floor(50000/7000) = 7;
  # chunks are rounded to 7 steps (49000 bp) to keep the grid aligned.
  spans <- data.frame(chr = "3", start = 1L, end = 200000L)
  win <- 10000L; step <- 7000L
  ch <- define_genome_chunks(spans, window_size = win, step_size = step,
                             chunk_size_mb = 0.05)

  # chunk_start must stay on the global step grid (chunk_start - span_start
  # divisible by step).
  expect_true(all((ch$chunk_start - 1L) %% step == 0L))

  got  <- .driver_chunk_starts(ch, "3", win, step)
  full <- .full_span_starts("3", 1L, 200000L, win, step)
  expect_equal(anyDuplicated(got), 0L)
  expect_equal(sort(got), sort(full))
})


test_that("est_n_windows matches the driver's define_regions_window row count per chunk", {
  spans <- data.frame(chr = c("1", "2"),
                      start = c(1L, 1L),
                      end = c(12000000L, 6000000L))
  win <- 100000L; step <- 10000L
  ch <- define_genome_chunks(spans, window_size = win, step_size = step,
                             chunk_size_mb = 5)

  for (i in seq_len(nrow(ch))) {
    w <- GLOWr::define_regions_window(ch$chr[i], ch$chunk_start[i],
                                      ch$span_end[i], win, step)
    w <- w[w$start <= ch$chunk_end[i], , drop = FALSE]
    expect_equal(ch$est_n_windows[i], nrow(w),
                 info = paste("chunk", i))
  }
  # Sum of owned windows equals the genome-wide full-span window count.
  full1 <- length(.full_span_starts("1", 1L, 12000000L, win, step))
  full2 <- length(.full_span_starts("2", 1L, 6000000L, win, step))
  expect_equal(sum(ch$est_n_windows), full1 + full2)
})


test_that("a tiny span yields a single one-window chunk", {
  spans <- data.frame(chr = "22", start = 1000L, end = 1500L)
  ch <- define_genome_chunks(spans, window_size = 2000L, step_size = 1000L,
                             chunk_size_mb = 5)
  expect_equal(nrow(ch), 1L)
  expect_equal(ch$est_n_windows, 1L)             # define_regions_window's max(1,)
  expect_equal(ch$chunk_start, 1000L)
  expect_equal(ch$chunk_end, 1500L)
})


test_that("non-fixed_mb strategy and a too-small chunk are errors", {
  spans <- data.frame(chr = "1", start = 1L, end = 1000000L)
  expect_error(
    define_genome_chunks(spans, 100000L, 10000L, chunk_strategy = "variant_balanced"),
    "fixed_mb")
  # chunk_size_bp (5000) < step_size (10000) -> a chunk cannot hold a window start
  expect_error(
    define_genome_chunks(spans, 100000L, 10000L, chunk_size_mb = 0.005),
    "smaller than step_size")
})


test_that("invalid chrom_spans is rejected", {
  expect_error(define_genome_chunks(data.frame(a = 1), 1000L, 500L),
               "columns")
  expect_error(
    define_genome_chunks(data.frame(chr = character(0), start = integer(0),
                                    end = integer(0)), 1000L, 500L),
    "no rows")
})
