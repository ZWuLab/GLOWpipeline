# ==============================================================================
# Tests for parse_run_log.R
# ==============================================================================
#
# Tests for parse_chr_log() (per-chr gene scan) and parse_chunk_log()
# (per-chunk sliding-window scan). Writes tiny synthetic logs into a temp dir
# and asserts the parsed $summary fields and $errors rows.


# -------------------- parse_chr_log() --------------------

write_chr_log <- function(lines, chr = 5L) {
  d <- tempfile("chrlog-")
  dir.create(d)
  p <- file.path(d, sprintf("log_chr%d.txt", chr))
  writeLines(lines, p)
  p
}

test_that("parse_chr_log extracts the summary fields", {
  p <- write_chr_log(c(
    "[ts] n_samples = 1234",
    "[ts] trait = binary",
    "[ts] Gene list (chr5): 678 coding genes",
    "[ts] Filter: rare_maf<0.01",
    "[ts] Done. genes_attempted=600 results=590 skipped=8 errors=2 runtime=42.5s"
  ))
  res <- parse_chr_log(p)
  s <- res$summary
  expect_equal(s$chr, 5L)
  expect_equal(s$n_samples, 1234L)
  expect_equal(s$trait, "binary")
  expect_equal(s$gene_list_n, 678L)
  expect_match(s$filter, "Filter:")
  expect_equal(s$attempted, 600L)
  expect_equal(s$results, 590L)
  expect_equal(s$skipped, 8L)
  expect_equal(s$errors, 2L)
  expect_equal(s$runtime_sec, 42.5)
})

test_that("parse_chr_log accepts the bare '(\\d+) genes' wording", {
  p <- write_chr_log(c("[ts] Gene list (chr5): 42 genes"))
  expect_equal(parse_chr_log(p)$summary$gene_list_n, 42L)
})

test_that("parse_chr_log parses the run_scan_unit 'Done gene unit' vocabulary", {
  p <- write_chr_log(c(
    "[ts] n_samples = 48198",
    "[ts] trait = binary",
    "[ts] Filter: rare_maf_cutoff=0.5, variant_type=SNV, min_variants=1",
    "[ts] Done gene unit chr5. n=427, n_ok=306, n_skip_empty=67, n_staar_skip_single=54, n_error=0, 689.2s"
  ))
  s <- parse_chr_log(p)$summary
  expect_equal(s$gene_list_n, 427L)   # unit count stands in for the absent gene-list line
  expect_equal(s$attempted, 427L)
  expect_equal(s$results, 360L)       # n_ok + n_staar_skip_single: both statuses write a row
  expect_equal(s$skipped, 67L)
  expect_equal(s$attempted, s$results + s$skipped + s$errors)   # the legacy identity holds
  expect_equal(s$errors, 0L)
  expect_equal(s$runtime_sec, 689.2)
})

test_that("parse_chr_log NA-fills when there is no Done line", {
  p <- write_chr_log(c("[ts] n_samples = 10"))
  s <- parse_chr_log(p)$summary
  expect_true(is.na(s$attempted))
  expect_true(is.na(s$runtime_sec))
})

test_that("parse_chr_log parses error rows and the source prefix", {
  p <- write_chr_log(c(
    "[ts] Done. genes_attempted=3 results=1 skipped=0 errors=2 runtime=1s",
    "[ts]   STAAR ERROR in GENE_A: singular matrix",
    "[ts]   ERROR in GENE_B: something failed"
  ))
  e <- parse_chr_log(p)$errors
  expect_equal(nrow(e), 2L)
  expect_equal(e$label, c("GENE_A", "GENE_B"))
  expect_equal(e$source, c("STAAR", "GLOW"))   # absent prefix -> GLOW
  expect_equal(e$message, c("singular matrix", "something failed"))
  expect_true(all(e$chr == 5L))
})

test_that("parse_chr_log returns 0-row errors df when no errors", {
  p <- write_chr_log(c(
    "[ts] Done. genes_attempted=3 results=3 skipped=0 errors=0 runtime=1s"))
  e <- parse_chr_log(p)$errors
  expect_equal(nrow(e), 0L)
  expect_equal(colnames(e), c("chr", "label", "source", "message"))
})

test_that("parse_chr_log errors when the log file is absent", {
  expect_error(parse_chr_log(tempfile("nope_log_chr1.txt")),
               "Log file not found")
})


# -------------------- parse_chunk_log() --------------------

write_chunk_log <- function(lines, chunk_id = 3L) {
  d <- tempfile("chunklog-")
  dir.create(d)
  p <- file.path(d, sprintf("log_chunk_%04d.txt", chunk_id))
  writeLines(lines, p)
  p
}

test_that("parse_chunk_log extracts the summary fields (with run name)", {
  p <- write_chunk_log(c(
    "[ts] Starting chunk 3 [run alpha]: chr7 [1000, 5000], est_n_windows=12",
    "[ts] Done chunk 3: n_windows=12 n_ok=10 n_skip_empty=1 n_staar_skip_single=0 n_error=1 engine_runtime=30.0s wall=33.0s"
  ))
  s <- parse_chunk_log(p)$summary
  expect_equal(s$chunk_id, 3L)
  expect_equal(s$run_name, "alpha")
  expect_equal(s$chr, "7")
  expect_equal(s$chunk_start, 1000L)
  expect_equal(s$chunk_end, 5000L)
  expect_equal(s$est_n_windows, 12L)
  expect_equal(s$n_windows, 12L)
  expect_equal(s$n_ok, 10L)
  expect_equal(s$n_skip_empty, 1L)
  expect_equal(s$n_staar_skip_single, 0L)
  expect_equal(s$n_error, 1L)
  expect_equal(s$engine_runtime_secs, 30.0)
  expect_equal(s$wall_secs, 33.0)
})

test_that("parse_chunk_log accepts the older header without a run name", {
  p <- write_chunk_log(c(
    "[ts] Starting chunk 3: chr7 [1000, 5000], est_n_windows=4",
    "[ts] Done chunk 3: n_windows=4 n_ok=4 n_skip_empty=0 n_staar_skip_single=0 n_error=0 engine_runtime=1.0s wall=1.0s"
  ))
  s <- parse_chunk_log(p)$summary
  expect_true(is.na(s$run_name))
  expect_equal(s$chr, "7")
})

test_that("parse_chunk_log parses per-window error rows", {
  p <- write_chunk_log(c(
    "[ts] Starting chunk 3: chr7 [1000, 9000], est_n_windows=3",
    "[ts] Done chunk 3: n_windows=3 n_ok=1 n_skip_empty=0 n_staar_skip_single=0 n_error=2 engine_runtime=1s wall=1s",
    "[ts] 2 window error(s):",
    "[ts]   chr7:1000-2000: boom one",
    "[ts]   chr7:3000-4000: boom two"
  ))
  e <- parse_chunk_log(p)$errors
  expect_equal(nrow(e), 2L)
  expect_equal(e$label, c("chr7:1000-2000", "chr7:3000-4000"))
  expect_equal(e$chr, c("7", "7"))
  expect_equal(e$message, c("boom one", "boom two"))
  expect_true(all(e$chunk_id == 3L))
})

test_that("parse_chunk_log NA-fills an incomplete chunk (no Done line)", {
  p <- write_chunk_log(c(
    "[ts] Starting chunk 3: chr7 [1000, 5000], est_n_windows=4"))
  s <- parse_chunk_log(p)$summary
  expect_true(is.na(s$n_windows))
  expect_true(is.na(s$wall_secs))
})
