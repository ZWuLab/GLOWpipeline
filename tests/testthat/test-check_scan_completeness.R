# ==============================================================================
# Tests for check_scan_completeness.R
# ==============================================================================
#
# Tests for check_scan_completeness() (expected-vs-present chunk gate) and the
# internal collapse_int_ranges() reporting util. Builds a temp run dir with a
# Stage-1 chunk table, some result files, and some chunk logs.


# Build a temp run dir with expected chunk ids; create result files for
# `present`, and chunk logs (with optional Done) per `logs` spec.
build_run_dir <- function(expected_ids, present_ids, logs = list(),
                          ext = "csv") {
  d <- tempfile("run-")
  dir.create(file.path(d, "shared"), recursive = TRUE)
  dir.create(file.path(d, "results"))
  dir.create(file.path(d, "logs"))
  saveRDS(data.frame(chunk_id = expected_ids),
          file.path(d, "shared", "01-genome_chunks.rds"))
  for (id in present_ids) {
    writeLines("dummy",
               file.path(d, "results", sprintf("scan_chunk_%d.%s", id, ext)))
  }
  for (id in names(logs)) {
    writeLines(logs[[id]],
               file.path(d, "logs",
                         sprintf("log_chunk_%04d.txt", as.integer(id))))
  }
  d
}

done_log <- function(chunk_id, n_error = 0L) c(
  sprintf("[ts] Starting chunk %d: chr1 [1000, 5000], est_n_windows=4", chunk_id),
  sprintf(paste0("[ts] Done chunk %d: n_windows=4 n_ok=%d n_skip_empty=0 ",
                 "n_staar_skip_single=0 n_error=%d engine_runtime=1s wall=1s"),
          chunk_id, 4L - n_error, n_error))

started_only_log <- function(chunk_id) c(
  sprintf("[ts] Starting chunk %d: chr1 [1000, 5000], est_n_windows=4", chunk_id))


test_that("check_scan_completeness reports complete when all present", {
  d <- build_run_dir(1:3, 1:3,
                     logs = list("1" = done_log(1), "2" = done_log(2),
                                 "3" = done_log(3)))
  res <- check_scan_completeness(d)
  expect_true(res$complete)
  expect_equal(res$n_expected, 3L)
  expect_equal(res$n_present, 3L)
  expect_equal(length(res$missing_ids), 0L)
})

test_that("check_scan_completeness flags missing chunks", {
  d <- build_run_dir(1:3, c(1L, 3L),
                     logs = list("1" = done_log(1), "3" = done_log(3)))
  res <- check_scan_completeness(d)
  expect_false(res$complete)
  expect_equal(res$missing_ids, 2L)
  expect_equal(res$n_present, 2L)
})

test_that("check_scan_completeness flags not-done chunks (result but no Done)", {
  # chunk 2 has a result file but its log never reached Done
  d <- build_run_dir(1:2, 1:2,
                     logs = list("1" = done_log(1),
                                 "2" = started_only_log(2)))
  res <- check_scan_completeness(d)
  expect_true(res$complete)             # all results present
  expect_equal(res$not_done_ids, 2L)    # but chunk 2 never finished
})

test_that("check_scan_completeness tallies window errors", {
  d <- build_run_dir(1:2, 1:2,
                     logs = list("1" = done_log(1, n_error = 1L),
                                 "2" = done_log(2, n_error = 2L)))
  res <- check_scan_completeness(d)
  expect_equal(res$n_window_errors, 3L)
})

test_that("check_scan_completeness errors without a chunk table", {
  d <- tempfile("run-empty-")
  dir.create(file.path(d, "results"), recursive = TRUE)
  expect_error(check_scan_completeness(d), "Chunk table not found")
})

test_that("check_scan_completeness status data.frame has the right schema", {
  d <- build_run_dir(1:2, 1L, logs = list("1" = done_log(1)))
  res <- check_scan_completeness(d)
  expect_equal(colnames(res$status),
               c("chunk_id", "has_result", "reached_done", "n_error"))
  expect_equal(res$status$has_result, c(TRUE, FALSE))
})


# -------------------- internal collapse_int_ranges() --------------------

test_that("collapse_int_ranges collapses consecutive runs", {
  f <- GLOWpipeline:::collapse_int_ranges
  expect_equal(f(c(3, 17, 18, 19, 200)), "3, 17-19, 200")
  expect_equal(f(integer(0)), "")
  expect_equal(f(c(5L)), "5")
  expect_equal(f(c(1L, 2L, 3L)), "1-3")
})
