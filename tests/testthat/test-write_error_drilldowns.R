# ==============================================================================
# Tests for write_error_drilldowns.R
# ==============================================================================
#
# Tests for write_error_drilldowns() using the chr22 synthetic aGDS fixture
# (helper-create-test-agds.R). Covers: the note-only path (non-ok resolver
# status) and the GDS re-extraction path (status="ok" -> per-variant CSV).


test_that("write_error_drilldowns writes a note-only CSV on a non-ok status", {
  err_rows <- data.frame(chr = "22", label = "GENE_X",
                         source = "GLOW", message = "boom",
                         stringsAsFactors = FALSE)
  resolve <- function(r) list(status = "NO_RESULT_ROW")

  err_dir <- tempfile("errdir-")
  out <- write_error_drilldowns(
    err_rows, resolve_regions = resolve,
    null_model = list(sample_id = "x"), filter = NULL,
    gds_dir = tempdir(), gds_pattern = "absent_chr{chr}.gds",
    pi_features = character(0), err_dir = err_dir, verbose = FALSE)

  csv <- file.path(err_dir, "GENE_X-variants.csv")
  expect_true(file.exists(csv))
  note <- read.csv(csv, stringsAsFactors = FALSE)
  expect_equal(note$status, "NO_RESULT_ROW")
  expect_match(note$note, "NO_RESULT_ROW")
  # the markdown link is recorded on the returned rows
  expect_match(out$error_diagnostics[1], "GENE_X-variants.csv")
})

test_that("write_error_drilldowns re-extracts variants on status='ok'", {
  skip_if_not_installed("SeqArray")
  skip_if_not_installed("SNPRelate")
  skip_if_not_installed("gdsfmt")

  # Build the chr22 fixture at the {chr}-substituted path.
  gds_dir <- tempfile("gds-")
  dir.create(gds_dir)
  gds_path <- file.path(gds_dir, "agds_chr22.gds")
  create_test_agds(gds_path)
  withr::defer(unlink(gds_path))

  feats <- c("cadd_phred", "linsight", "fathmm_xf")
  spec  <- GLOWr::variant_filter(rare_maf_cutoff = 0.5, variant_type = "SNV",
                                 min_mac = 1L, min_variants = 1L)
  gds0  <- SeqArray::seqOpen(gds_path, readonly = TRUE)
  sids  <- SeqArray::seqGetData(gds0, "sample.id")
  SeqArray::seqClose(gds0)
  null_model <- list(sample_id = sids)

  # GENE_A spans pos 1000-5000 on chr22.
  err_rows <- data.frame(chr = "22", label = "GENE_A",
                         stringsAsFactors = FALSE)
  resolve <- function(r) list(
    status = "ok",
    regions = list(list(chr = "22", start = 1000L, end = 5000L,
                        label = "GENE_A")))

  err_dir <- tempfile("errdir-ok-")
  out <- write_error_drilldowns(
    err_rows, resolve_regions = resolve,
    null_model = null_model, filter = spec,
    gds_dir = gds_dir, gds_pattern = "agds_chr{chr}.gds",
    pi_features = feats, err_dir = err_dir, verbose = FALSE)

  csv <- file.path(err_dir, "GENE_A-variants.csv")
  expect_true(file.exists(csv))
  d <- read.csv(csv, stringsAsFactors = FALSE)
  expect_gt(nrow(d), 0L)
  # the per-variant cohort genotype summary columns are present
  expect_true(all(c("n_g0_homo_alt", "n_g1_het", "n_g2_homo_ref",
                    "n_missing_genotypes", "variance_in_cohort",
                    "is_constant_in_cohort", "constant_value") %in%
                  colnames(d)))
  # location recorded as chr:start-end
  expect_equal(out$location[1], "22:1000-5000")
})

test_that("write_error_drilldowns errors on missing required columns", {
  expect_error(
    write_error_drilldowns(
      data.frame(chr = "22"), resolve_regions = function(r) list(status = "x"),
      null_model = list(sample_id = "x"), filter = NULL,
      gds_dir = tempdir(), gds_pattern = "p{chr}.gds",
      pi_features = character(0), err_dir = tempfile()),
    "missing required column")
})

test_that("write_error_drilldowns writes GDS_NOT_FOUND when the GDS is absent", {
  err_rows <- data.frame(chr = "22", label = "GENE_Y",
                         stringsAsFactors = FALSE)
  resolve <- function(r) list(
    status = "ok",
    regions = list(list(chr = "22", start = 1L, end = 2L, label = "GENE_Y")))

  err_dir <- tempfile("errdir-nogds-")
  expect_warning(
    out <- write_error_drilldowns(
      err_rows, resolve_regions = resolve,
      null_model = list(sample_id = "x"), filter = NULL,
      gds_dir = tempdir(), gds_pattern = "definitely_absent_chr{chr}.gds",
      pi_features = character(0), err_dir = err_dir, verbose = FALSE),
    "GDS not found")
  note <- read.csv(file.path(err_dir, "GENE_Y-variants.csv"),
                   stringsAsFactors = FALSE)
  expect_equal(note$status, "GDS_NOT_FOUND")
})
