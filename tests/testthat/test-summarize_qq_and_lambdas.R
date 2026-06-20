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
