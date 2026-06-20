# ==============================================================================
# Tests for aggregate_scan_results() — genome rbind + threshold + locus merge
# ==============================================================================
#
# Synthetic per-chunk flat tables written to tempdir(); no GDS needed. The
# load-bearing case is the chunk-boundary merge: overlapping significant windows
# split across two chunk FILES must merge into one locus (design §10.3).


# A representative flat-table subset: identity + count/flag + three p-value
# columns (GLOW_Omni is the default primary_test).
.mk_chunk <- function(chr, start, end, glow_omni, cMAC = 50) {
  n <- length(start)
  data.frame(
    label              = paste0("chr", chr, ":", start, "-", end),
    chr                = as.character(chr),
    start              = as.integer(start),
    end                = as.integer(end),
    midpoint           = as.integer((start + end) / 2),
    n_after_annotation = rep(8L, n),
    n_after_collapse   = rep(5L, n),
    cMAC               = as.numeric(cMAC),
    is_single          = rep(FALSE, n),
    GLOW_Omni          = as.numeric(glow_omni),
    GLOW_SKAT_equ      = as.numeric(glow_omni),
    STAAR_O            = as.numeric(glow_omni),
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )
}

.write_chunk <- function(df, id, dir = tempfile("aggchunks_")) {
  if (!dir.exists(dir)) dir.create(dir)
  path <- file.path(dir, sprintf("scan_chunk_%04d.csv", id))
  utils::write.csv(df, path, row.names = FALSE)
  path
}


test_that("concatenation puts chunks in genomic order; n_windows is total", {
  d <- tempfile("agg_order_"); dir.create(d)
  f1 <- .write_chunk(.mk_chunk("2", c(100, 200), c(150, 250), 0.5), 1, d)
  f2 <- .write_chunk(.mk_chunk("1", c(1000, 2000), c(1050, 2050), 0.5), 2, d)
  f3 <- .write_chunk(.mk_chunk("1", 500, 550, 0.5), 3, d)

  agg <- aggregate_scan_results(c(f1, f2, f3))
  expect_equal(nrow(agg$scan_results_all), 5L)
  expect_equal(agg$summary$n_windows, 5L)
  # chr 1 (sorted by start) precedes chr 2.
  expect_equal(agg$scan_results_all$chr, c("1", "1", "1", "2", "2"))
  expect_equal(agg$scan_results_all$start[1:3], c(500L, 1000L, 2000L))
})


test_that("thresholding selects exactly the windows below alpha / n_windows", {
  d <- tempfile("agg_thr_"); dir.create(d)
  # 4 windows -> threshold = 0.05 / 4 = 0.0125. One window below it.
  f <- .write_chunk(
    .mk_chunk("1", c(100, 200, 300, 400), c(150, 250, 350, 450),
              glow_omni = c(1e-5, 0.5, 0.5, 0.5)), 1, d)

  agg <- aggregate_scan_results(f, primary_test = "GLOW_Omni", alpha = 0.05)
  expect_equal(agg$summary$n_windows, 4L)
  expect_equal(agg$summary$threshold, 0.05 / 4)
  expect_equal(nrow(agg$scan_sig_windows), 1L)
  expect_equal(agg$scan_sig_windows$start, 100L)
})


test_that("overlapping/abutting significant windows merge into one locus", {
  d <- tempfile("agg_overlap_"); dir.create(d)
  # 3 windows, all significant (each p << 0.05/3); peak is the middle one.
  f <- .write_chunk(
    .mk_chunk("1", c(100, 200, 350), c(300, 400, 500),
              glow_omni = c(1e-4, 1e-8, 1e-4)), 1, d)

  agg <- aggregate_scan_results(f, merge_gap = 0L)
  expect_equal(nrow(agg$scan_sig_windows), 3L)
  expect_equal(nrow(agg$scan_loci), 1L)

  locus <- agg$scan_loci
  expect_equal(locus$chr, "1")
  expect_equal(locus$locus_start, 100L)
  expect_equal(locus$locus_end, 500L)
  expect_equal(locus$n_windows, 3L)
  expect_equal(locus$peak_label, "chr1:200-400")        # min-p member
  expect_equal(locus$peak_GLOW_Omni, 1e-8)
  # peak block carries every p-value column, prefixed peak_.
  expect_true(all(c("peak_GLOW_Omni", "peak_GLOW_SKAT_equ", "peak_STAAR_O") %in%
                  names(locus)))
})


test_that("significant windows separated by more than merge_gap stay separate", {
  d <- tempfile("agg_gap_"); dir.create(d)
  f <- .write_chunk(
    .mk_chunk("1", c(100, 100000), c(200, 100100), glow_omni = c(1e-6, 1e-6)),
    1, d)
  agg <- aggregate_scan_results(f, merge_gap = 0L)
  expect_equal(nrow(agg$scan_loci), 2L)
  expect_equal(agg$scan_loci$locus_start, c(100L, 100000L))
})


test_that("overlapping significant windows split across two chunk files merge", {
  d <- tempfile("agg_boundary_"); dir.create(d)
  # Same chr, overlapping, but in different chunk files -> only an assembled-
  # table merge can join them (design §10.3 chunk-boundary guarantee).
  fa <- .write_chunk(.mk_chunk("1", 1000, 2000, 1e-4), 1, d)
  fb <- .write_chunk(.mk_chunk("1", 1500, 2500, 1e-4), 2, d)

  agg <- aggregate_scan_results(c(fa, fb))
  expect_equal(agg$summary$n_windows, 2L)
  expect_equal(nrow(agg$scan_loci), 1L)
  expect_equal(agg$scan_loci$locus_start, 1000L)
  expect_equal(agg$scan_loci$locus_end, 2500L)
  expect_equal(agg$scan_loci$n_windows, 2L)
})


test_that("no significant windows: 0-row sig + loci tables with correct columns", {
  d <- tempfile("agg_none_"); dir.create(d)
  f <- .write_chunk(.mk_chunk("1", c(100, 200), c(150, 250), 0.9), 1, d)

  agg <- aggregate_scan_results(f)
  expect_equal(nrow(agg$scan_sig_windows), 0L)
  expect_equal(nrow(agg$scan_loci), 0L)
  expect_equal(agg$summary$n_loci, 0L)
  expect_true(all(c("locus_id", "chr", "locus_start", "locus_end", "n_windows",
                    "peak_label", "peak_GLOW_Omni") %in% names(agg$scan_loci)))
})


test_that("cMAC cutoff drops low-cMAC windows before thresholding", {
  d <- tempfile("agg_cmac_"); dir.create(d)
  f <- .write_chunk(
    .mk_chunk("1", c(100, 200, 300), c(150, 250, 350),
              glow_omni = c(0.5, 0.5, 0.5), cMAC = c(5, 50, 100)), 1, d)

  agg <- aggregate_scan_results(f, cmac_cutoff = 10)
  expect_equal(agg$summary$n_dropped_cmac, 1L)           # the cMAC == 5 window
  expect_equal(agg$summary$n_windows, 2L)                # denominator after drop
  expect_false(5 %in% agg$scan_results_all$cMAC)
})


test_that("empty chunk_files and unknown primary_test are errors", {
  d <- tempfile("agg_err_"); dir.create(d)
  f <- .write_chunk(.mk_chunk("1", 100, 150, 0.5), 1, d)
  expect_error(aggregate_scan_results(character(0)), "empty")
  expect_error(aggregate_scan_results(f, primary_test = "NOPE"),
               "not in the chunk")
})


test_that("fst round-trips chunk tables (skipped when fst is absent)", {
  skip_if_not_installed("fst")
  d <- tempfile("agg_fst_"); dir.create(d)
  df <- .mk_chunk("1", c(100, 200), c(150, 250), c(1e-6, 0.5))
  p1 <- file.path(d, "scan_chunk_0001.fst")
  fst::write_fst(df, p1)

  agg <- aggregate_scan_results(p1, output_format = "fst")
  expect_equal(nrow(agg$scan_results_all), 2L)
  expect_equal(nrow(agg$scan_sig_windows), 1L)
})
