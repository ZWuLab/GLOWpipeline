# ==============================================================================
# Tests for read_all_chr_results.R
# ==============================================================================
#
# Tests for read_all_chr_results() (concat 22 per-chr CSVs in genomic order).


# Write 22 per-chr CSVs (two rows each, deliberately out of start order) under a
# temp dir using a given filename pattern; return the dir.
write_chr_csvs <- function(pattern = "glow_chr%d.csv") {
  d <- tempfile("results-")
  dir.create(d)
  for (chr in 1:22) {
    df <- data.frame(
      label = c(sprintf("g%d_b", chr), sprintf("g%d_a", chr)),
      chr   = c(chr, chr),
      start = c(2000L, 1000L),          # descending so the sort is testable
      GLOW_Omni = c(0.5, 0.2),
      stringsAsFactors = FALSE)
    write.csv(df, file.path(d, sprintf(pattern, chr)), row.names = FALSE)
  }
  d
}

test_that("read_all_chr_results concatenates all 22 chrs", {
  d <- write_chr_csvs()
  res <- read_all_chr_results(d)
  expect_equal(nrow(res), 44L)
  expect_setequal(unique(res$chr), 1:22)
})

test_that("read_all_chr_results returns genomic order (chr, then start)", {
  d <- write_chr_csvs()
  res <- read_all_chr_results(d)
  # chr 1 rows come first, ascending start within chr
  expect_equal(res$chr[1:2], c(1L, 1L))
  expect_equal(res$start[1:2], c(1000L, 2000L))
  # last rows are chr 22
  expect_equal(tail(res$chr, 1), 22L)
})

test_that("read_all_chr_results default pattern reproduces glow_chr%d.csv", {
  expect_equal(formals(read_all_chr_results)$pattern, "glow_chr%d.csv")
  d <- write_chr_csvs("glow_chr%d.csv")
  # default-call output is identical to explicit-pattern call
  res_default <- read_all_chr_results(d)
  res_explicit <- read_all_chr_results(d, pattern = "glow_chr%d.csv")
  expect_identical(res_default, res_explicit)
})

test_that("read_all_chr_results honors a custom pattern", {
  d <- write_chr_csvs("scan_chr%d.csv")
  res <- read_all_chr_results(d, pattern = "scan_chr%d.csv")
  expect_equal(nrow(res), 44L)
})

test_that("read_all_chr_results errors on a missing chr CSV", {
  d <- write_chr_csvs()
  file.remove(file.path(d, "glow_chr7.csv"))
  expect_error(read_all_chr_results(d), "Missing per-chr result CSVs")
})

test_that("read_all_chr_results errors on inconsistent schemas", {
  d <- write_chr_csvs()
  # rewrite chr5 with an extra column
  bad <- data.frame(label = "x", chr = 5L, start = 1L, GLOW_Omni = 0.1,
                    extra = 1, stringsAsFactors = FALSE)
  write.csv(bad, file.path(d, "glow_chr5.csv"), row.names = FALSE)
  expect_error(read_all_chr_results(d), "column schemas differ")
})

test_that("read_all_chr_results default chroms is 1:22", {
  expect_equal(eval(formals(read_all_chr_results)$chroms), 1:22)
})

test_that("read_all_chr_results honors a chroms subset", {
  d <- write_chr_csvs()
  res <- read_all_chr_results(d, chroms = c(1L, 22L))
  expect_equal(nrow(res), 4L)
  expect_setequal(unique(res$chr), c(1L, 22L))
  # genomic order preserved within the subset
  expect_equal(res$chr, c(1L, 1L, 22L, 22L))
})

test_that("read_all_chr_results(chroms=) ignores files outside the subset", {
  d <- write_chr_csvs()
  file.remove(file.path(d, "glow_chr7.csv"))   # not requested -> no error
  res <- read_all_chr_results(d, chroms = c(1L, 22L))
  expect_equal(nrow(res), 4L)
})

test_that("read_all_chr_results(chroms=) errors when a REQUESTED chr is missing", {
  d <- write_chr_csvs()
  file.remove(file.path(d, "glow_chr22.csv"))
  expect_error(read_all_chr_results(d, chroms = c(1L, 22L)),
               "Missing per-chr result CSVs")
})
