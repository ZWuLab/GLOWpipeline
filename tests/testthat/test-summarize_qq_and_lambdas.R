# ==============================================================================
# Tests for summarize_qq_and_lambdas.R
# ==============================================================================
#
# Tests for summarize_qq_and_lambdas() on a tiny results frame: assert the
# computed lambda values match GLOWr::compute_lambda_gc and that the PDF/RDS
# artifacts are written (PDF bytes are NOT compared -- they carry timestamps).


# Tiny results frame: 300 rows, two GLOW-style test columns. Inflate one column
# so its lambda is clearly above 1, keep the omni near-null.
make_results <- function() {
  set.seed(20260606)
  n <- 300L
  z_null <- rnorm(n)
  z_infl <- rnorm(n) * 1.3
  data.frame(
    label     = paste0("g", seq_len(n)),
    GLOW_Omni = stats::pnorm(abs(z_null), lower.tail = FALSE) * 2,
    GLOW_SKAT = stats::pnorm(abs(z_infl), lower.tail = FALSE) * 2,
    stringsAsFactors = FALSE)
}

test_that("summarize_qq_and_lambdas computes per-test lambdas matching GLOWr", {
  df <- make_results()
  out_dir <- tempfile("qq-")
  res <- summarize_qq_and_lambdas(
    df, label = "all", output_dir = out_dir,
    test_columns = c("GLOW_Omni", "GLOW_SKAT"),
    omni_test = "GLOW_Omni", verbose = FALSE)

  expect_equal(res$n, nrow(df))
  expect_equal(res$lambda_omni,
               GLOWr::compute_lambda_gc(df$GLOW_Omni))
  expect_equal(res$lambdas$GLOW_SKAT,
               GLOWr::compute_lambda_gc(df$GLOW_SKAT))
  # the inflated column has the larger lambda
  expect_gt(res$lambdas$GLOW_SKAT, res$lambdas$GLOW_Omni)
})

test_that("summarize_qq_and_lambdas writes the omni PDF and lambda RDS", {
  df <- make_results()
  out_dir <- tempfile("qq-")
  res <- summarize_qq_and_lambdas(
    df, label = "all", output_dir = out_dir,
    test_columns = c("GLOW_Omni", "GLOW_SKAT"),
    omni_test = "GLOW_Omni", verbose = FALSE)

  expect_true(file.exists(res$paths$qq_omni))
  expect_match(basename(res$paths$qq_omni), "^qq_omni_all\\.pdf$")
  expect_true(file.exists(res$paths$family_lambdas_rds))
  # the RDS holds the same per-test lambdas
  saved <- readRDS(res$paths$family_lambdas_rds)
  expect_equal(saved, res$lambdas)
  # no family grid was requested
  expect_true(is.na(res$paths$qq_by_family))
})

test_that("summarize_qq_and_lambdas writes a family grid when families given", {
  df <- make_results()
  out_dir <- tempfile("qq-")
  res <- summarize_qq_and_lambdas(
    df, label = "fam", output_dir = out_dir,
    families = list(SKAT = "GLOW_SKAT", Omni = "GLOW_Omni"),
    omni_test = "GLOW_Omni", verbose = FALSE)
  expect_true(file.exists(res$paths$qq_by_family))
  expect_match(basename(res$paths$qq_by_family), "^qq_by_family_fam\\.pdf$")
})

test_that("summarize_qq_and_lambdas skips the RDS when write_rds = FALSE", {
  df <- make_results()
  out_dir <- tempfile("qq-")
  res <- summarize_qq_and_lambdas(
    df, label = "all", output_dir = out_dir,
    test_columns = c("GLOW_Omni", "GLOW_SKAT"),
    omni_test = "GLOW_Omni", write_rds = FALSE, verbose = FALSE)
  expect_true(is.na(res$paths$family_lambdas_rds))
})

test_that("summarize_qq_and_lambdas validates omni_test membership", {
  df <- make_results()
  expect_error(
    summarize_qq_and_lambdas(
      df, label = "all", output_dir = tempfile("qq-"),
      test_columns = "GLOW_Omni", omni_test = "GLOW_SKAT",
      verbose = FALSE),
    "must be one of")
})

# ------------------------------------------------------------------------------
# Display thinning of the QQ pages
# ------------------------------------------------------------------------------

# A larger frame: 200,000 rows, three tests in two families, a few small p-values in each.
make_large_results <- function(n = 200000L) {
  set.seed(20261117)
  df <- data.frame(label = paste0("w", seq_len(n)),
                   GLOW_Omni = stats::runif(n), GLOW_SKAT_a = stats::runif(n),
                   GLOW_SKAT_b = stats::runif(n), stringsAsFactors = FALSE)
  df$GLOW_Omni[1:10]   <- 10^-seq(6, 15, length.out = 10)
  df$GLOW_SKAT_a[1:10] <- 10^-seq(5, 12, length.out = 10)
  df
}

test_that("thinning leaves the lambdas unchanged and bounds the size of both QQ files", {
  df <- make_large_results()
  fams <- list(SKAT = c("GLOW_SKAT_a", "GLOW_SKAT_b"), Omni = "GLOW_Omni")
  d1 <- tempfile("qq-thin-"); d2 <- tempfile("qq-full-")
  on.exit(unlink(c(d1, d2), recursive = TRUE), add = TRUE)
  thinned <- summarize_qq_and_lambdas(df, label = "all", output_dir = d1, families = fams,
                                      omni_test = "GLOW_Omni", verbose = FALSE)
  full    <- summarize_qq_and_lambdas(df, label = "all", output_dir = d2, families = fams,
                                      omni_test = "GLOW_Omni", verbose = FALSE, thin = FALSE)
  # the numbers are of the complete input either way
  expect_identical(thinned$lambdas, full$lambdas)
  expect_identical(readRDS(thinned$paths$family_lambdas_rds), readRDS(full$paths$family_lambdas_rds))
  # the per-family page: 600,000 points when nothing is thinned, a few thousand otherwise
  size <- function(res, which) file.size(res$paths[[which]])
  expect_lt(size(thinned, "qq_by_family"), size(full, "qq_by_family") / 10)
  expect_lt(size(thinned, "qq_by_family"), 2e5)
  # the omnibus page follows the same argument
  expect_lt(size(thinned, "qq_omni"), size(full, "qq_omni") / 3)
})

test_that("the points of a QQ curve that stand apart are all drawn when thinning", {
  df <- make_large_results()
  v <- df$GLOW_Omni
  obs <- -log10(sort(v)); expd <- -log10((seq_along(obs) - 0.5) / length(obs))
  grDevices::pdf(NULL, width = 10, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, 2))
  plot(NA, NA, xlim = c(0, max(expd)), ylim = c(0, max(obs)))
  keep <- .thin_points(expd, obs, .thin_cells(NULL, cex = 0.5))
  # the ten planted p-values are the ten largest observed values, each in a cell of its own
  expect_true(all(keep[1:10]))
  # the bulk of the curve is thinned: far fewer points than rows, and more than the tail alone
  expect_lt(sum(keep), 2000L)
  expect_gt(sum(keep), 100L)
})

test_that("the thinning helpers give what GLOWr's give on the same points", {
  ns <- asNamespace("GLOWr")
  skip_if_not(exists(".thin_points", envir = ns, inherits = FALSE) &&
                exists(".thin_cells", envir = ns, inherits = FALSE),
              "this GLOWr has no thinning helpers to compare with")
  grDevices::pdf(NULL, width = 6, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)
  set.seed(11)
  x <- stats::runif(5e4, 0, 6); y <- stats::runif(5e4, 0, 12)
  x[1:4] <- c(-1, 7, NA, 3); y[1:4] <- c(1, 1, 1, 13)        # outside the frame, or missing
  plot(c(0, 6), c(0, 12), type = "n")
  expect_identical(.thin_points(x, y, 57), get(".thin_points", envir = ns)(x, y, 57))
  plot(c(1, 1000), c(0.001, 12), type = "n", log = "xy")       # logarithmic axes
  xl <- 10^stats::runif(2e4, 0, 3)
  expect_identical(.thin_points(xl, y[seq_along(xl)] + 0.001, 30),
                   get(".thin_points", envir = ns)(xl, y[seq_along(xl)] + 0.001, 30))
  for (cex in c(0.3, 0.5, 1))
    expect_identical(.thin_cells(NULL, cex), get(".thin_cells", envir = ns)(NULL, cex))
  expect_identical(.thin_cells(100, 0.5), 100)
  expect_error(.thin_cells(-1, 0.5))
})

test_that("summarize_qq_and_lambdas refuses a `thin` that is not TRUE or FALSE", {
  df <- make_results()
  for (bad in list(NA, "yes", c(TRUE, FALSE), 1))
    expect_error(
      summarize_qq_and_lambdas(df, label = "all", output_dir = tempfile("qq-"),
                               test_columns = c("GLOW_Omni", "GLOW_SKAT"),
                               omni_test = "GLOW_Omni", verbose = FALSE, thin = bad),
      "`thin` must be TRUE or FALSE")
})
