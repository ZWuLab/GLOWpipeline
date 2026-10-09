# ==============================================================================
# Tests for prepare_scan_run() — the shared Stage-1 spine
# ==============================================================================
#
# Exercises the spine on the chr22 synthetic aGDS fixture (helper-create-test-
# agds.R): output-dir creation, GDS-presence check, pheno load + sample
# alignment + null-model fit, PI-model load + PI-feature resolve, B-source XOR
# dispatch, input checksums, and the base config_used snapshot fields.


# ---- Shared fixtures (built once for the file) ----

# Synthetic aGDS named so a {chr} pattern resolves to it for chr22.
.psr_dir <- tempfile("psr-gds-")
dir.create(.psr_dir)
.psr_gds <- file.path(.psr_dir, "chr22_test.gds")
create_test_agds(.psr_gds)
withr::defer(unlink(.psr_dir, recursive = TRUE), teardown_env())

# Sample order + a prepared-pheno object matching the fixture's 50 samples.
.psr_g0   <- SeqArray::seqOpen(.psr_gds, readonly = TRUE)
.psr_sids <- SeqArray::seqGetData(.psr_g0, "sample.id")
SeqArray::seqClose(.psr_g0)
.psr_n    <- length(.psr_sids)

set.seed(20260606)
.psr_prep <- list(
  sample_id = sample(.psr_sids),                 # deliberately shuffled vs GDS
  Y         = rbinom(.psr_n, 1, 0.3),
  X         = cbind(intercept = 1, cov1 = rnorm(.psr_n)),
  trait     = "binary"
)
.psr_pheno <- tempfile(fileext = ".rds")
saveRDS(.psr_prep, .psr_pheno)
withr::defer(unlink(.psr_pheno), teardown_env())

# A 2-model PI ensemble written as model_*.rds (what load_PI_models reads).
.psr_pi_dir <- tempfile("psr-pi-")
dir.create(.psr_pi_dir)
.psr_make_pi <- function(seed) {
  set.seed(seed)
  feats <- c("cadd_phred", "linsight", "fathmm_xf")
  xx <- matrix(runif(300 * length(feats)), 300, length(feats),
               dimnames = list(NULL, feats))
  yy <- rbinom(300, 1, plogis(xx %*% rep(0.4, length(feats))))
  glm(yy ~ xx, family = "binomial")
}
saveRDS(.psr_make_pi(1), file.path(.psr_pi_dir, "model_1.rds"))
saveRDS(.psr_make_pi(2), file.path(.psr_pi_dir, "model_2.rds"))
withr::defer(unlink(.psr_pi_dir, recursive = TRUE), teardown_env())

# A trained B model envelope (the structure load_B_model validates) for the
# b_model_path arm. The validator needs a glow_B_model with a predictable lm
# under $models$beta_method on the (X, fX, logX, logfX) design, plus metadata
# and format_version. Built directly so the test is self-contained.
.psr_b_path <- tempfile(fileext = ".rds")
.psr_b_built <- tryCatch({
  set.seed(7)
  X    <- runif(500, 0.001, 0.5)
  fX   <- X * (1 - X)
  logX <- log(X); logfX <- log(fX)
  b    <- sqrt(pmax(-0.14307 * log(fX) - 0.19522, 1e-6))
  lm_obj <- lm(b ~ X + fX + logX + logfX,
               data = data.frame(b, X, fX, logX, logfX))
  b_model <- structure(
    list(models = list(beta_method = lm_obj, pvalue_method = NULL)),
    class = "glow_B_model")
  envelope <- list(
    model = b_model,
    metadata = list(glowr_version = as.character(utils::packageVersion("GLOWr"))),
    format_version = "1.0")
  saveRDS(envelope, .psr_b_path)
  TRUE
}, error = function(e) FALSE)

# A common config skeleton (env) used by most tests; b source filled per-test.
.psr_base_config <- function(b = c("func", "model")) {
  b <- match.arg(b)
  env <- new.env()
  env$gds_dir         <- .psr_dir
  env$gds_pattern     <- "chr{chr}_test.gds"
  env$pheno_path      <- .psr_pheno
  env$pi_model_dir    <- .psr_pi_dir
  env$filter_spec     <- list(rare_maf_cutoff = 0.5, variant_type = "SNV",
                              min_mac = 1L, min_variants = 1L)
  env$ld_threshold    <- 0.9
  env$mac_threshold   <- 10L
  env$collapse_method <- "mean"
  env$pi_features     <- c("cadd_phred", "linsight", "fathmm_xf")
  if (b == "func") {
    env$b_func <- function(maf) rep(0.3, length(maf))
  } else {
    env$b_model_path <- .psr_b_path
  }
  env
}


# ---- Tests ----

test_that("prepare_scan_run runs the spine end-to-end (b_func arm)", {
  out  <- tempfile("psr-out-"); withr::defer(unlink(out, recursive = TRUE))
  cfg  <- .psr_base_config("func")
  bnd  <- prepare_scan_run(cfg, output_dir = out, chroms = 22L,
                           config_path = NULL, verbose = 0)

  # Bundle has the documented components.
  expect_named(bnd, c("null_model", "pi_models", "pi_features", "b_func",
                      "b_model", "filter_spec", "sample_ids", "checksums",
                      "base_snapshot", "gds_files"))

  # Null model fitted on the aligned cohort, in GDS sample order.
  expect_equal(bnd$null_model$n, .psr_n)
  expect_equal(bnd$null_model$trait, "binary")
  expect_equal(bnd$sample_ids, .psr_sids)             # GDS order, not pheno order
  expect_equal(bnd$null_model$sample_id, .psr_sids)

  # PI resolve: explicit features echoed back.
  expect_equal(bnd$pi_features, c("cadd_phred", "linsight", "fathmm_xf"))
  expect_equal(bnd$pi_models$n_models, 2L)

  # B source: func arm sets b_func, leaves b_model NULL.
  expect_true(is.function(bnd$b_func))
  expect_null(bnd$b_model)

  # Checksums present for each input class.
  expect_named(bnd$checksums,
               c("pheno_path", "b_model_path", "pi_model_files", "gds_files"))
  expect_true(is.na(bnd$checksums$b_model_path))      # no path in func arm
  expect_match(unname(bnd$checksums$pheno_path), "^[0-9a-f]{32}$")
})

test_that("prepare_scan_run writes the spine-owned shared artifacts", {
  out <- tempfile("psr-out-"); withr::defer(unlink(out, recursive = TRUE))
  bnd <- prepare_scan_run(.psr_base_config("func"), output_dir = out,
                          chroms = 22L)
  expect_true(file.exists(file.path(out, "shared", "01-null_model.rds")))
  expect_true(file.exists(file.path(out, "shared", "01-input_checksums.rds")))
  # The written null model equals the returned one.
  expect_equal(readRDS(file.path(out, "shared", "01-null_model.rds")),
               bnd$null_model)
  expect_equal(readRDS(file.path(out, "shared", "01-input_checksums.rds")),
               bnd$checksums)
})

test_that("prepare_scan_run base_snapshot carries the shared config fields", {
  out <- tempfile("psr-out-"); withr::defer(unlink(out, recursive = TRUE))
  bnd <- prepare_scan_run(.psr_base_config("func"), output_dir = out,
                          chroms = 22L, config_path = NULL)
  bs  <- bnd$base_snapshot

  for (f in c("config_path", "config_source", "gds_dir", "gds_pattern",
              "pheno_path", "pi_model_dir", "b_model_path", "b_func",
              "filter_spec", "ld_threshold", "mac_threshold", "collapse_method",
              "pi_features", "output_dir", "trait", "glow_test_used_spa",
              "git_sha", "session_info", "prepared_at")) {
    expect_true(f %in% names(bs), info = paste("missing base field:", f))
  }

  expect_equal(bs$gds_dir, .psr_dir)
  expect_equal(bs$output_dir, out)
  expect_equal(bs$trait, "binary")
  expect_true(bs$glow_test_used_spa)           # binary trait
  expect_null(bs$b_model_path)                 # func arm -> no path field
  expect_true(is.function(bs$b_func))
  expect_equal(bs$pi_features, c("cadd_phred", "linsight", "fathmm_xf"))
})

test_that("config_path NULL -> NULL snapshot path/source; a real path is read", {
  out <- tempfile("psr-out-"); withr::defer(unlink(out, recursive = TRUE))
  b1  <- prepare_scan_run(.psr_base_config("func"), output_dir = out,
                          chroms = 22L, config_path = NULL)
  expect_null(b1$base_snapshot$config_path)
  expect_null(b1$base_snapshot$config_source)

  cfile <- tempfile(fileext = ".R")
  writeLines(c("# a config", "x <- 1"), cfile)
  withr::defer(unlink(cfile))
  out2 <- tempfile("psr-out-"); withr::defer(unlink(out2, recursive = TRUE))
  b2 <- prepare_scan_run(.psr_base_config("func"), output_dir = out2,
                         chroms = 22L, config_path = cfile)
  expect_equal(b2$base_snapshot$config_path, cfile)
  expect_equal(b2$base_snapshot$config_source, readLines(cfile))
})

test_that("PI features default to GLOWr defaults when config pi_features is NULL", {
  out <- tempfile("psr-out-"); withr::defer(unlink(out, recursive = TRUE))
  cfg <- .psr_base_config("func")
  cfg$pi_features <- NULL
  bnd <- prepare_scan_run(cfg, output_dir = out, chroms = 22L)
  expect_equal(bnd$pi_features, GLOWr:::.default_PI_features())
  expect_true(length(bnd$pi_features) > 3L)
})

test_that("B-source XOR is enforced (both set / neither set both error)", {
  out <- tempfile("psr-out-"); withr::defer(unlink(out, recursive = TRUE))
  both <- .psr_base_config("func"); both$b_model_path <- .psr_b_path
  expect_error(prepare_scan_run(both, output_dir = out, chroms = 22L),
               "exactly one of")
  neither <- .psr_base_config("func"); neither$b_func <- NULL
  expect_error(prepare_scan_run(neither, output_dir = out, chroms = 22L),
               "exactly one of")
})

test_that("a non-finite b_func probe is rejected", {
  out <- tempfile("psr-out-"); withr::defer(unlink(out, recursive = TRUE))
  cfg <- .psr_base_config("func")
  cfg$b_func <- function(maf) rep(NA_real_, length(maf))
  expect_error(prepare_scan_run(cfg, output_dir = out, chroms = 22L),
               "finite scalar numeric")
})

test_that("missing per-chr GDS files error informatively", {
  out <- tempfile("psr-out-"); withr::defer(unlink(out, recursive = TRUE))
  cfg <- .psr_base_config("func")
  # chr1 GDS does not exist under the fixture dir.
  expect_error(prepare_scan_run(cfg, output_dir = out, chroms = c(1L, 22L)),
               "Missing per-chr GDS files")
})

test_that("make_dirs creates the requested subdirectories", {
  out <- tempfile("psr-out-"); withr::defer(unlink(out, recursive = TRUE))
  prepare_scan_run(.psr_base_config("func"), output_dir = out, chroms = 22L,
                   make_dirs = c("results", "logs"))
  expect_true(dir.exists(file.path(out, "shared")))
  expect_true(dir.exists(file.path(out, "results")))
  expect_true(dir.exists(file.path(out, "logs")))
})

test_that("the b_model_path arm loads the model and checksums the path", {
  skip_if_not(.psr_b_built, "B model fixture could not be built in this env")
  out <- tempfile("psr-out-"); withr::defer(unlink(out, recursive = TRUE))
  bnd <- prepare_scan_run(.psr_base_config("model"), output_dir = out,
                          chroms = 22L)
  expect_null(bnd$b_func)
  expect_false(is.null(bnd$b_model))
  expect_equal(bnd$base_snapshot$b_model_path, .psr_b_path)
  expect_match(unname(bnd$checksums$b_model_path), "^[0-9a-f]{32}$")
})



# ---- .psr_check_qc_node(): missing QC values (a VCF FILTER of ".") ----

# A three-variant VCF converted with SeqArray; returns the GDS path.
.psr_vcf_gds <- function(filters) {
  vcf <- tempfile(fileext = ".vcf")
  writeLines(c("##fileformat=VCFv4.2", "##contig=<ID=22>",
               '##FILTER=<ID=LowQual,Description="Low quality">',
               '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
               sprintf("22\t%d\trs%d\tA\tG\t.\t%s\t.\tGT\t0/1\t0/0",
                       c(100L, 200L, 300L), 1:3, filters)), vcf)
  out <- tempfile(fileext = ".gds")
  SeqArray::seqVCF2GDS(vcf, out, verbose = FALSE)
  unlink(vcf)
  out
}

test_that(".psr_check_qc_node() is silent when no QC value is missing", {
  expect_silent(res <- GLOWpipeline:::.psr_check_qc_node(c("22" = .psr_gds),
                                                          qc_label = "annotation/filter"))
  expect_equal(res$n_missing_qc, 0L)
})

test_that(".psr_check_qc_node() warns with the counts when some QC values are missing", {
  g <- .psr_vcf_gds(c(".", "PASS", "LowQual"))
  on.exit(unlink(g), add = TRUE)
  expect_warning(res <- GLOWpipeline:::.psr_check_qc_node(c("22" = g), qc_label = "annotation/filter"),
                 "1 variant\\(s\\) have a missing value")
  expect_equal(res$n_variants, 3L)
  expect_equal(res$n_missing_qc, 1L)
})

test_that(".psr_check_qc_node() stops when every QC value of a chromosome is missing", {
  g <- .psr_vcf_gds(c(".", ".", "."))
  on.exit(unlink(g), add = TRUE)
  expect_error(GLOWpipeline:::.psr_check_qc_node(c("22" = g), qc_label = "annotation/filter"),
               "Every variant on chr 22")
})
