# ==============================================================================
# Tests for glow_output_policy() — the artifact-retention value object
# ==============================================================================


test_that("defaults: flat table on, all sidecars off, correct class", {
  p <- glow_output_policy()
  expect_s3_class(p, "glow_output_policy")
  expect_true(p$write_flat_table)
  expect_false(p$write_evidence)
  expect_false(p$write_staar_detail)
  expect_false(p$write_result_objects)
})


test_that("write_flat_table = FALSE is rejected", {
  expect_error(glow_output_policy(write_flat_table = FALSE),
               "must be TRUE")
})


test_that("non-logical or NA flags are rejected", {
  expect_error(glow_output_policy(write_evidence = "yes"),
               "single non-NA logical")
  expect_error(glow_output_policy(write_evidence = NA),
               "single non-NA logical")
  expect_error(glow_output_policy(write_result_objects = c(TRUE, FALSE)),
               "single non-NA logical")
  expect_error(glow_output_policy(write_flat_table = 1L),
               "single non-NA logical")
})


test_that("write_staar_detail = TRUE alone constructs cleanly (engine validates the context)", {
  p <- expect_no_error(glow_output_policy(write_staar_detail = TRUE))
  expect_true(p$write_staar_detail)
  expect_s3_class(p, "glow_output_policy")
})


test_that("print method runs and is invisible", {
  p <- glow_output_policy(write_evidence = TRUE)
  expect_output(print(p), "glow_output_policy")
  expect_output(print(p), "evidence")
  expect_invisible(print(p))
})
