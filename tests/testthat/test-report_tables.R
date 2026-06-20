# ==============================================================================
# Tests for report_tables.R
# ==============================================================================
#
# Tests for gc_calibration_fields() and build_lambda_table_md().


`%||%` <- function(a, b) if (is.null(a)) b else a   # caller-supplied null_or


# -------------------- gc_calibration_fields() --------------------

test_that("gc_calibration_fields resolves inflation/method/z_scale", {
  cfg <- list(inflation_factor = 1.44,
              calibration = list(method = "ldsc_intercept"))
  res <- gc_calibration_fields(cfg, `%||%`)
  expect_equal(res$inflation, 1.44)
  expect_equal(res$method, "ldsc_intercept")
  expect_equal(res$z_scale, 1 / sqrt(1.44))
})

test_that("gc_calibration_fields defaults to 1.0 / manual when absent", {
  res <- gc_calibration_fields(list(), `%||%`)
  expect_equal(res$inflation, 1.0)
  expect_equal(res$method, "manual")
  expect_equal(res$z_scale, 1.0)
})


# -------------------- build_lambda_table_md() --------------------

make_lambda_df <- function() {
  # two methods, three categories (incl ALL), common-support rows + an off-support
  # row that must be filtered out.
  rbind(
    data.frame(category = "missense",   method = "GLOW_Omni",
               support = "common_ge2", n = 100L, lambda_gc = 1.05),
    data.frame(category = "missense",   method = "STAAR_O_native",
               support = "common_ge2", n = 100L, lambda_gc = 1.20),
    data.frame(category = "synonymous", method = "GLOW_Omni",
               support = "common_ge2", n = 80L,  lambda_gc = 0.98),
    data.frame(category = "synonymous", method = "STAAR_O_native",
               support = "common_ge2", n = 80L,  lambda_gc = NA_real_),
    data.frame(category = "ALL",        method = "GLOW_Omni",
               support = "common_ge2", n = 180L, lambda_gc = 1.02),
    data.frame(category = "ALL",        method = "STAAR_O_native",
               support = "common_ge2", n = 180L, lambda_gc = 1.10),
    data.frame(category = "missense",   method = "GLOW_Omni",
               support = "rare_only",  n = 50L,  lambda_gc = 9.99),
    stringsAsFactors = FALSE)
}

test_that("build_lambda_table_md emits header + category rows in order", {
  lam <- make_lambda_df()
  md <- build_lambda_table_md(lam, methods = c("GLOW_Omni", "STAAR_O_native"))
  expect_equal(md[1], "| category | GLOW_Omni | STAAR_O_native |")
  expect_match(md[2], "^\\|---\\|")
  # categories sorted (missense, synonymous), then ALL appended last
  cat_lines <- md[-(1:2)]
  expect_match(cat_lines[1], "^\\| missense ")
  expect_match(cat_lines[2], "^\\| synonymous")
  expect_match(cat_lines[3], "^\\| ALL ")
})

test_that("build_lambda_table_md tags the negative control and renders NA", {
  lam <- make_lambda_df()
  md <- build_lambda_table_md(lam, methods = c("GLOW_Omni", "STAAR_O_native"))
  syn_line <- grep("synonymous", md, value = TRUE)
  expect_match(syn_line, "\\*\\*\\[neg ctrl\\]\\*\\*")
  expect_match(syn_line, "0.980")
  expect_match(syn_line, "\\| NA \\|")   # the NA STAAR cell
})

test_that("build_lambda_table_md filters to the requested support level", {
  lam <- make_lambda_df()
  md <- build_lambda_table_md(lam, methods = "GLOW_Omni", support = "rare_only")
  # rare_only has only the missense row; ALL still appended (no rare ALL -> NA)
  expect_match(grep("missense", md, value = TRUE), "9.990")
})

test_that("build_lambda_table_md show_n adds (n=NN) to cells", {
  lam <- make_lambda_df()
  md <- build_lambda_table_md(lam, methods = "GLOW_Omni", show_n = TRUE)
  mis_line <- grep("missense", md, value = TRUE)
  expect_match(mis_line, "1.050 \\(n=100\\)")
})

test_that("build_lambda_table_md appends extra_note after a blank line", {
  lam <- make_lambda_df()
  md <- build_lambda_table_md(lam, methods = "GLOW_Omni",
                              extra_note = "note line")
  expect_equal(md[length(md)], "note line")
  expect_equal(md[length(md) - 1L], "")
})

test_that("build_lambda_table_md can suppress the ALL row", {
  lam <- make_lambda_df()
  md <- build_lambda_table_md(lam, methods = "GLOW_Omni",
                              append_all_row = FALSE)
  expect_false(any(grepl("^\\| ALL ", md)))
})
