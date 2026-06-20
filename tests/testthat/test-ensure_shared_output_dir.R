# ==============================================================================
# Tests for ensure_shared_output_dir() — portable no-op cases only (CI-safe)
# ==============================================================================
#
# ensure_shared_output_dir() is the sanctioned born-shared hook (shared-storage
# impl §8/§8a). It is a deliberate, portable NO-OP in every "can't / shouldn't"
# case, falling back to a plain local dir created by the caller. These tests
# exercise ONLY those no-op cases — they never touch /project and never append
# to the tracked manifest:
#   - a tempfile() path (outside any git repo)        -> FALSE, no symlink
#   - a path whose PARENT does not exist              -> FALSE, no symlink
#
# The born-shared SUCCESS path (shelling to scaffold_shared_output.sh -> /project
# + appending tools/shared-storage/shared-paths.tsv) is intentionally NOT tested
# here: that would pollute the tracked manifest. It is exercised end-to-end by
# prepare_scan_run's own coverage / the opt-in smoke in the adoption task.


# A path is "an existing symlink" iff Sys.readlink() returns a non-NA, non-empty
# string. Sys.readlink() returns NA for a path that does not exist at all, and
# nzchar(NA) is TRUE, so a bare nzchar() check would mis-fire on a missing path.
.is_symlink <- function(p) {
  rl <- Sys.readlink(p)
  !is.na(rl) && nzchar(rl)
}

test_that("ensure_shared_output_dir is a no-op outside a git repo (tempfile path)", {
  # A tempfile() under /tmp is not inside the package repo, so the hook must
  # decline (no scaffold helper reachable) and create nothing — leaving the
  # caller free to make a plain local dir.
  out <- tempfile("eso-norepo-")           # parent (tempdir()) exists; leaf does not
  withr::defer(unlink(out, recursive = TRUE))

  made <- ensure_shared_output_dir(out, share_root = NULL, verbose = 0)

  expect_false(made)                        # returned FALSE
  expect_false(file.exists(out))            # created nothing
  expect_false(dir.exists(out))             # created no directory
  expect_false(.is_symlink(out))            # created no symlink
})

test_that("ensure_shared_output_dir is a no-op with an explicit share_root too", {
  # The explicit-share_root arm must ALSO decline outside a repo: the share_root
  # is only consulted when a repo + scaffold helper are present.
  out <- tempfile("eso-norepo-share-")
  withr::defer(unlink(out, recursive = TRUE))

  made <- ensure_shared_output_dir(out, share_root = file.path(tempdir(), "ss_root"),
                                   verbose = 0)

  expect_false(made)
  expect_false(file.exists(out))
  expect_false(dir.exists(out))
  expect_false(.is_symlink(out))
})

test_that("ensure_shared_output_dir is a no-op when the parent does not exist", {
  # No parent dir -> nowhere to host a symlink -> the hook declines immediately,
  # before any repo / helper probing.
  base <- tempfile("eso-noparent-")         # does NOT exist
  out  <- file.path(base, "missing", "outputs")
  withr::defer(unlink(base, recursive = TRUE))

  expect_false(dir.exists(dirname(out)))    # precondition: parent absent
  made <- ensure_shared_output_dir(out, share_root = NULL, verbose = 0)

  expect_false(made)
  expect_false(file.exists(out))
  expect_false(dir.exists(out))
  expect_false(.is_symlink(out))
})
