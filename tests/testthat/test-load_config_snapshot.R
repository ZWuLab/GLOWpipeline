# ==============================================================================
# Tests for load_config_snapshot.R
# ==============================================================================
#
# Tests for load_config_snapshot() (reads <output_dir>/shared/01-config_used.rds).


test_that("load_config_snapshot reads the Stage-1 snapshot", {
  out_dir <- tempfile("cfg-")
  dir.create(file.path(out_dir, "shared"), recursive = TRUE)
  cfg <- list(trait = "binary", inflation_factor = 1.05)
  saveRDS(cfg, file.path(out_dir, "shared", "01-config_used.rds"))

  res <- load_config_snapshot(out_dir)
  expect_equal(res, cfg)
})

test_that("load_config_snapshot errors when the snapshot is missing", {
  out_dir <- tempfile("cfg-missing-")
  dir.create(out_dir)
  expect_error(load_config_snapshot(out_dir),
               "config snapshot not found")
})

test_that("load_config_snapshot requires output_dir (no default)", {
  expect_true("output_dir" %in% names(formals(load_config_snapshot)))
  expect_identical(formals(load_config_snapshot)$output_dir,
                   quote(expr = ))  # empty default -> required arg
})
