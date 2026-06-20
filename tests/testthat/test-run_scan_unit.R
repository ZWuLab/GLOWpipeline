# ==============================================================================
# Tests for run_scan_unit() — the unified Stage-2 driver body
# ==============================================================================
#
# Runs the synthetic aGDS fixture (helper-create-test-agds.R) through the
# unified driver for each region type and asserts the written per-unit table.
# Kept fast: the gene/coding gene table is the two multi-variant fixture genes
# (GENE_A, GENE_B) supplied via config$resolved_gene_list (gene) or by mocking
# GLOWr::define_regions_gene() to those two rows (coding).


# ---- Shared fixtures (built once for the file) ----

.agds <- tempfile(fileext = ".gds")
create_test_agds(.agds)
withr::defer(unlink(.agds), teardown_env())

.feats <- c("cadd_phred", "linsight", "fathmm_xf")
.spec_args <- list(rare_maf_cutoff = 0.5, variant_type = "SNV",
                   min_mac = 1L, min_variants = 1L)
.spec <- do.call(GLOWr::variant_filter, .spec_args)

.gds0 <- SeqArray::seqOpen(.agds, readonly = TRUE)
.sids <- SeqArray::seqGetData(.gds0, "sample.id")
.n    <- length(.sids)
.meds <- GLOWr::compute_annotation_medians(.gds0, .feats, .spec,
                                           sample_id = .sids, verbose = 0)
SeqArray::seqClose(.gds0)

set.seed(20260526)
.X  <- cbind(1, rnorm(.n))
.Y  <- rbinom(.n, 1, 0.3)
.nm <- GLOWr::fit_null_model(.X, .Y, trait = "binary", sample_id = .sids)

.make_glm_pi <- function(seed) {
  set.seed(seed)
  xx <- matrix(runif(300 * length(.feats)), 300, length(.feats))
  yy <- rbinom(300, 1, plogis(xx %*% rep(0.4, length(.feats))))
  glm(y ~ x, data = list(x = xx, y = yy), family = "binomial")
}
.pi_models <- list(models = list(.make_glm_pi(11), .make_glm_pi(12)),
                   model_type = "GLM", n_models = 2L)
.b_func <- function(maf) rep(0.3, length(maf))

# The two multi-variant fixture genes.
.fixt_genes <- data.frame(
  region_id = c("GENE_A", "GENE_B"), chr = "22",
  start = c(1000L, 10000L), end = c(5000L, 15000L),
  label = c("GENE_A", "GENE_B"), stringsAsFactors = FALSE)

# A base config carrying the engine knobs every type needs.
.base_config <- function(extra = list()) {
  env <- new.env()
  env$filter_spec     <- .spec_args
  env$ld_threshold    <- 0.95
  env$mac_threshold   <- 11L
  env$collapse_method <- "mean"
  for (nm in names(extra)) assign(nm, extra[[nm]], envir = env)
  env
}

# Convenience: bind the fixed engine inputs.
.run <- function(region_type, config, results_path, output_dir, ...) {
  run_scan_unit(.agds, .nm, .pi_models, .feats, .meds,
                b_func = .b_func, b_model = NULL,
                region_type = region_type, config = config,
                results_path = results_path, output_dir = output_dir,
                output_policy = glow_output_policy(), logger = NULL, ...)
}

.SUMMARY_COLS <- c("n_collapsed_units", "min_MAF", "max_MAF",
                   "min_MAC", "max_MAC", "min_p_standard", "min_p_SPA")


# ---- gene: flat table WITH the 7 summary columns ----

test_that("gene-type run writes a flat table WITH the 7 region-summary columns", {
  od <- withr::local_tempdir()
  env <- .base_config(list(resolved_gene_list = .fixt_genes))
  rp <- file.path(od, "results", "glow_chr22.csv")

  ret <- .run("gene", env, rp, od, chr = 22)
  expect_true(file.exists(rp))

  tbl <- read.csv(rp, stringsAsFactors = FALSE)
  expect_equal(nrow(tbl), 2L)                          # GENE_A, GENE_B
  expect_setequal(tbl$label, c("GENE_A", "GENE_B"))
  expect_true(all(.SUMMARY_COLS %in% names(tbl)))      # the 7 summary columns
  # Engine flat schema: ...cMAC, is_single, [7 summary], GLOW_*. The summary
  # block sits after is_single and immediately before the GLOW_* columns.
  nms <- names(tbl)
  expect_equal(which(nms == "n_collapsed_units"), which(nms == "is_single") + 1L)
  expect_equal(which(nms == GLOWpipeline:::.GLOW_PVALUE_COLS[1L]),
               which(nms == "min_p_SPA") + 1L)
  expect_true("midpoint" %in% nms)
  expect_equal(ret$results_path, rp)
})


# ---- window: per-unit window table, NO summary columns ----

test_that("window-type run writes the expected per-unit window table", {
  od <- withr::local_tempdir()
  env <- .base_config(list(window_size = 2000L, step_size = 1000L))
  # span_end past the fixture variants so multiple windows fall on real data.
  chunk <- data.frame(chr = "22", chunk_id = 1L,
                      chunk_start = 1000L, chunk_end = 15000L, span_end = 16000L)
  rp <- file.path(od, "results", "scan_chunk_0001.csv")

  ret <- .run("window", env, rp, od, chunk = chunk)
  expect_true(file.exists(rp))

  tbl <- read.csv(rp, stringsAsFactors = FALSE)
  expect_true(nrow(tbl) > 0L)
  expect_true("midpoint" %in% names(tbl))
  expect_false(any(.SUMMARY_COLS %in% names(tbl)))     # window: no summary cols
  # The flat schema's GLOW columns are present.
  expect_true(all(GLOWpipeline:::.GLOW_PVALUE_COLS %in% names(tbl)))
})


# ---- coding: tag `category`, stack categories, gene-first column order ----

test_that("coding-type run tags `category` and stacks categories", {
  testthat::local_mocked_bindings(
    define_regions_gene = function(genes_info = NULL, chr = NULL, extend = 0L)
      .fixt_genes,
    .package = "GLOWr")

  od <- withr::local_tempdir()
  # Two categories that hit the fixture's GENE_A/GENE_B variants. GENE_A has
  # nonsynonymous + synonymous SNVs; GENE_B has nonsynonymous SNVs.
  env <- .base_config(list(categories = c("missense", "synonymous")))
  rp <- file.path(od, "results-glow", "glow_coding_chr22.csv")

  ret <- .run("coding", env, rp, od, chr = 22)
  expect_true(file.exists(rp))

  tbl <- read.csv(rp, stringsAsFactors = FALSE)
  expect_true("category" %in% names(tbl))
  # Gene-first column order: label, chr, category lead the table.
  expect_equal(names(tbl)[1:3], c("label", "chr", "category"))
  # Both categories produced at least one (gene, category) cell -> stacked.
  expect_true(nrow(tbl) >= 2L)
  expect_setequal(unique(tbl$category), c("missense", "synonymous"))
  # No region-summary columns on the coding schema.
  expect_false(any(.SUMMARY_COLS %in% names(tbl)))
})
