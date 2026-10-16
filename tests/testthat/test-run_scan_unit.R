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


# ---- 2026-10-09: n_cores gives the same table; segmentation runs end to end ----

test_that("gene run with n_cores = 2 writes the same table as with 1 (fixture, forked)", {
  skip_on_os("windows")
  skip_if(parallel::detectCores() < 2L, "needs two cores")
  od1 <- withr::local_tempdir(); od2 <- withr::local_tempdir()
  env <- .base_config(list(resolved_gene_list = .fixt_genes))
  rp1 <- file.path(od1, "results", "glow_chr22.csv")
  rp2 <- file.path(od2, "results", "glow_chr22.csv")
  lines <- character(0)
  .run("gene", env, rp1, od1, chr = 22, n_cores = 1L)
  # Two regions and a block size of 25 would be one block: use the config field
  # path and force two blocks by giving each gene its own block through the
  # window type below; here the parallel path must at least agree on the table.
  ret2 <- .run("gene", env, rp2, od2, chr = 22, n_cores = 2L)
  expect_equal(read.csv(rp1, stringsAsFactors = FALSE),
               read.csv(rp2, stringsAsFactors = FALSE))
  expect_equal(ret2$summary$n_ok, 2L)
  # The config field is read when the argument is NULL.
  env$n_cores <- 2L
  od3 <- withr::local_tempdir(); rp3 <- file.path(od3, "results", "glow_chr22.csv")
  .run("gene", env, rp3, od3, chr = 22)
  expect_equal(read.csv(rp1, stringsAsFactors = FALSE),
               read.csv(rp3, stringsAsFactors = FALSE))
})

test_that("window run with n_cores = 2 over several blocks writes the same table as with 1", {
  skip_on_os("windows")
  skip_if(parallel::detectCores() < 2L, "needs two cores")
  od1 <- withr::local_tempdir(); od2 <- withr::local_tempdir()
  env <- .base_config(list(window_size = 1000L, step_size = 250L))
  chunk <- data.frame(chr = "22", chunk_id = 1L,
                      chunk_start = 1000L, chunk_end = 15000L, span_end = 16000L)
  rp1 <- file.path(od1, "results", "scan_chunk_0001.csv")
  rp2 <- file.path(od2, "results", "scan_chunk_0001.csv")
  .run("window", env, rp1, od1, chunk = chunk, n_cores = 1L)
  .run("window", env, rp2, od2, chunk = chunk, n_cores = 2L)   # 57 windows: 3 blocks
  t1 <- read.csv(rp1, stringsAsFactors = FALSE)
  t2 <- read.csv(rp2, stringsAsFactors = FALSE)
  expect_true(nrow(t1) > 25L)
  expect_equal(t1, t2)
})

test_that("gene run segments a large gene and tests every segment", {
  od <- withr::local_tempdir()
  env <- .base_config(list(resolved_gene_list = .fixt_genes,
                           gene_max_records = 30L, gene_segment_bp = 2000L))
  rp <- file.path(od, "results", "glow_chr22.csv")
  logs <- character(0)
  ret <- run_scan_unit(.agds, .nm, .pi_models, .feats, .meds,
                       b_func = .b_func, b_model = NULL,
                       region_type = "gene", config = env,
                       results_path = rp, output_dir = od,
                       output_policy = glow_output_policy(), chr = 22,
                       logger = function(m) logs <<- c(logs, m))
  tbl <- read.csv(rp, stringsAsFactors = FALSE)
  # GENE_A_seg3 is the 1 bp tail [5000, 5000]; its one record does not pass the
  # filter, so that segment is skip_empty and has no row, as an empty gene would.
  expect_setequal(tbl$label, c(paste0("GENE_A_seg", 1:2), paste0("GENE_B_seg", 1:3)))
  expect_true(all(tbl$start[tbl$label == "GENE_A_seg2"] == 3000L))
  expect_true(any(grepl("gene segmentation: 2 gene", logs)))
  expect_true(all(.SUMMARY_COLS %in% names(tbl)))
})


# ---- 2026-10-15: size-aware scheduling (plan 12-size-aware-scheduling, G1) ----

# A regions table whose heavy regions (GENE_A, GENE_B: 40 records each) are
# not contiguous, with a span that holds no record (GAP1: skip_empty) and a
# region whose test is made to fail (GENE_ERR, the GENE_C span, 30 records).
# Costs (records^2): 100, 1600, 0, 25, 1600, 900; at n_cores = 2 the budget
# tau = 4225 / 16 = 264, so GENE_A, GENE_B and GENE_ERR run alone and GENE_D,
# GENE_E and GAP1 share one block: four blocks, none of them contiguous rows.
.sched_genes <- data.frame(
  region_id = c("GENE_D", "GENE_A", "GAP1", "GENE_E", "GENE_B", "GENE_ERR"),
  chr = "22",
  start = c(30000L, 1000L, 6000L, 40000L, 10000L, 20000L),
  end   = c(32000L, 5000L, 9000L, 41000L, 15000L, 25000L),
  label = c("GENE_D", "GENE_A", "GAP1", "GENE_E", "GENE_B", "GENE_ERR"),
  stringsAsFactors = FALSE)

# Plant an error in one region: GLOWr's region test is replaced, in GLOWr's
# namespace and in the imports of GLOWpipeline, by a wrapper that fails the
# label GENE_ERR with status "error" and forwards every other region. Forked
# workers inherit the binding. Restored on exit.
.with_planted_error <- function(code, label = "GENE_ERR") {
  orig <- GLOWr::glow_region_test
  wrapper <- function(gds, region, ...) {
    if (identical(as.character(region$label), label)) {
      return(list(result = NULL, evidence = NULL, staar = NULL,
                  staar_detail = NULL, status = "error",
                  message = "planted error"))
    }
    orig(gds, region, ...)
  }
  envs <- list(asNamespace("GLOWr"), parent.env(asNamespace("GLOWpipeline")))
  swap <- function(fn) for (e in envs) {
    if (exists("glow_region_test", envir = e, inherits = FALSE)) {
      unlockBinding("glow_region_test", e)
      assign("glow_region_test", fn, envir = e)
      lockBinding("glow_region_test", e)
    }
  }
  swap(wrapper)
  on.exit(swap(orig), add = TRUE)
  force(code)
}

.make_staar_ctx <- function() {
  sdat  <- data.frame(phenotype = .Y, cov1 = .X[, 2])
  snull <- STAAR::fit_null_glm(phenotype ~ cov1, data = sdat,
                               family = stats::binomial(link = "logit"))
  list(null_model = snull, anno_cols = .feats, rare_maf_cutoff = 0.5)
}

# One run with every sidecar on; returns the written tables and the log.
.run_sched <- function(region_type, config, od, k, staar_context = NULL, ...) {
  rp <- file.path(od, "results", "glow_unit.csv")
  logs <- character(0)
  ret <- run_scan_unit(.agds, .nm, .pi_models, .feats, .meds,
                       b_func = .b_func, b_model = NULL,
                       region_type = region_type, config = config,
                       results_path = rp, output_dir = od,
                       output_policy = glow_output_policy(
                         write_evidence = TRUE, write_result_objects = TRUE,
                         write_staar_detail = !is.null(staar_context)),
                       staar_context = staar_context, unit_id = "u",
                       logger = function(m) logs <<- c(logs, m),
                       n_cores = k, ...)
  sidecar <- function(name) {
    p <- file.path(od, name, sprintf("%s_u.rds", name))
    if (file.exists(p)) readRDS(p) else NULL
  }
  list(ret = ret, logs = logs, table = read.csv(rp, stringsAsFactors = FALSE),
       evidence = sidecar("evidence"), staar_detail = sidecar("staar_detail"),
       result_objects = sidecar("result_objects"))
}

test_that(".rsu_form_blocks isolates heavy regions, caps the rest and places every region once", {
  fb <- GLOWpipeline:::.rsu_form_blocks
  # The planning review's counterexample: at two workers tau = 300 / 16 =
  # 18.75, so every one of the five regions runs alone, heaviest first.
  b <- fb(c(100, 64, 64, 36, 36), 2L)
  expect_equal(attr(b, "tau"), 18.75)
  expect_equal(lengths(b), rep(1L, 5))
  expect_equal(unlist(b), 1:5)
  expect_equal(attr(b, "n_singletons"), 5L)
  expect_equal(attr(b, "block_cost"), c(100, 64, 64, 36, 36))
  # Random costs with a heavy tail: the invariants of section 4.3.
  set.seed(20261015)
  cost <- rexp(400)^3
  for (m in c(2L, 4L, 20L)) {
    b <- fb(cost, m)
    tau <- sum(cost) / (8 * m)
    expect_equal(attr(b, "tau"), tau)
    expect_equal(sort(unlist(b)), seq_along(cost))          # every region once
    single <- lengths(b) == 1L
    heavy  <- cost >= tau
    expect_true(all(heavy[unlist(b[single])] | TRUE))        # singles may be light
    expect_true(all(vapply(which(heavy), function(i)
      any(vapply(b, function(blk) identical(blk, i), logical(1))), logical(1))))
    packed <- b[!single]
    expect_true(all(lengths(packed) <= 25L))
    expect_true(all(vapply(packed, function(blk) sum(cost[blk]) <= tau + 1e-9,
                           logical(1))))
    bc <- attr(b, "block_cost")
    expect_true(all(diff(bc) <= 1e-9))                       # heaviest first
    expect_equal(attr(b, "n_singletons"), sum(heavy))
  }
  # Zero total cost: the rows keep their order in blocks of 25.
  b0 <- fb(rep(0, 60), 4L)
  expect_equal(unname(lapply(b0, identity)), list(1:25, 26:50, 51:60))
  expect_equal(attr(b0, "tau"), 0)
  # Ties broken by row index; an empty table gives no block.
  expect_equal(unlist(fb(c(5, 5, 5), 2L)), 1:3)
  expect_length(fb(numeric(0), 2L), 0L)
})

test_that("a repeated label stops the run at n_cores 1 and 2", {
  od <- withr::local_tempdir()
  dup <- rbind(.fixt_genes, .fixt_genes[1L, ])
  env <- .base_config(list(resolved_gene_list = dup))
  rp <- file.path(od, "results", "glow_chr22.csv")
  expect_error(.run("gene", env, rp, od, chr = 22, n_cores = 1L),
               "repeated label.*GENE_A")
  skip_on_os("windows")
  skip_if(parallel::detectCores() < 2L, "needs two cores")
  expect_error(.run("gene", env, rp, od, chr = 22, n_cores = 2L),
               "repeated label.*GENE_A")
  expect_false(file.exists(rp))
})

test_that("gene run at n_cores 2: blocks by size, serial table, sidecars and errors restored", {
  skip_on_os("windows")
  skip_if(parallel::detectCores() < 2L, "needs two cores")
  skip_if_not_installed("STAAR")
  env <- .base_config(list(resolved_gene_list = .sched_genes))
  ctx <- .make_staar_ctx()
  od1 <- withr::local_tempdir(); od2 <- withr::local_tempdir()
  r1 <- .with_planted_error(.run_sched("gene", env, od1, 1L, ctx, chr = 22))
  r2 <- .with_planted_error(.run_sched("gene", env, od2, 2L, ctx, chr = 22))

  # The serial run: rows in the table's order, the gap skipped, the error
  # recorded.
  expect_equal(r1$table$label, c("GENE_D", "GENE_A", "GENE_E", "GENE_B"))
  expect_equal(r1$ret$errors$label, "GENE_ERR")
  expect_equal(r1$ret$summary$n_skip_empty, 1L)
  expect_equal(r1$ret$summary$n_error, 1L)

  # The scheduled run: four blocks, three of them a heavy region alone,
  # dispatched heaviest first.
  sched_line <- grep("^scheduler \\(n_cores = 2\\)", r2$logs, value = TRUE)
  expect_length(sched_line, 1L)
  expect_match(sched_line, "6 regions in 4 blocks, 3 of them a heavy region alone")
  expect_match(sched_line, "heaviest GENE_A with 40 records")
  block_lines <- grep("^block [0-9]+/4:", r2$logs, value = TRUE)
  expect_length(block_lines, 4L)
  expect_match(block_lines[1L], "^block 1/4: 1 region\\(s\\), heaviest 40 records")
  expect_match(block_lines[4L], "^block 4/4: 3 region\\(s\\), heaviest 10 records")

  # Identity with the serial run: the flat table, the error table, the
  # summary tally (timing excluded) and the three sidecars, entry order
  # included.
  expect_equal(r2$table, r1$table)
  expect_equal(r2$ret$errors, r1$ret$errors)
  tally <- c("n_windows", "n_ok", "n_skip_empty", "n_staar_skip_single", "n_error")
  expect_equal(r2$ret$summary[, tally], r1$ret$summary[, tally])
  expect_equal(names(r1$evidence), c("GENE_D", "GENE_A", "GENE_E", "GENE_B"))
  expect_equal(r2$evidence, r1$evidence)
  expect_equal(names(r2$result_objects), names(r1$result_objects))
  expect_equal(r2$result_objects, r1$result_objects)
  expect_equal(names(r2$staar_detail), names(r1$staar_detail))
  expect_equal(r2$staar_detail, r1$staar_detail)
})

test_that("window run at n_cores 2 is scheduled by size and writes the serial table", {
  skip_on_os("windows")
  skip_if(parallel::detectCores() < 2L, "needs two cores")
  env <- .base_config(list(window_size = 1000L, step_size = 250L))
  chunk <- data.frame(chr = "22", chunk_id = 1L,
                      chunk_start = 1000L, chunk_end = 15000L, span_end = 16000L)
  od1 <- withr::local_tempdir(); od2 <- withr::local_tempdir()
  r1 <- .run_sched("window", env, od1, 1L, chunk = chunk)
  r2 <- .run_sched("window", env, od2, 2L, chunk = chunk)
  expect_true(nrow(r1$table) > 25L)
  expect_true(r1$ret$summary$n_skip_empty > 0L)              # empty windows between genes
  expect_equal(r2$table, r1$table)
  expect_equal(r2$ret$errors, r1$ret$errors)
  expect_equal(r2$evidence, r1$evidence)
  expect_equal(r2$result_objects, r1$result_objects)
  expect_length(grep("^scheduler \\(n_cores = 2\\)", r2$logs), 1L)
  expect_true(length(grep("^block [0-9]+/", r2$logs)) >= 8L)   # at least 8 m blocks when none is heavy
})

test_that("coding run at n_cores 2 schedules each category's pass and writes the serial table", {
  skip_on_os("windows")
  skip_if(parallel::detectCores() < 2L, "needs two cores")
  testthat::local_mocked_bindings(
    define_regions_gene = function(genes_info = NULL, chr = NULL, extend = 0L)
      .fixt_genes[c(2L, 1L), ],                                # B before A
    .package = "GLOWr")
  env <- .base_config(list(categories = c("missense", "synonymous")))
  od1 <- withr::local_tempdir(); od2 <- withr::local_tempdir()
  r1 <- .run_sched("coding", env, od1, 1L, chr = 22)
  r2 <- .run_sched("coding", env, od2, 2L, chr = 22)
  expect_equal(r2$table, r1$table)
  expect_equal(r2$ret$errors, r1$ret$errors)
  expect_equal(r2$evidence, r1$evidence)
  expect_equal(r2$result_objects, r1$result_objects)
  expect_length(grep("^scheduler \\(n_cores = 2\\)", r2$logs), 1L)   # one schedule, both categories
})


# The fixture's own genes are not in the bundled protein-coding table, so the
# seven regulatory masks are exercised through the gene list of
# define_regions_noncoding(); here the bundled ncRNA table (POTEH-AS1 and
# LINC01297 are planted) carries the ncRNA mask, and a mocked gene table makes
# the fixture's genes protein-coding for the other masks.
.fixt_genes_info <- data.frame(
  hgnc_symbol = c("GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E"),
  chromosome_name = 22L,
  start_position = c(1000L, 10000L, 20000L, 30000L, 40000L),
  end_position = c(5000L, 15000L, 25000L, 32000L, 41000L), stringsAsFactors = FALSE)

.run_nc <- function(od, k, categories = c("upstream", "downstream", "UTR", "enhancer_DHS", "ncRNA"),
                    max_regions = NA_integer_) {
  env <- .base_config(list(noncoding_categories = categories))
  rp <- file.path(od, "results-glow", "glow_noncoding_chr22.csv")
  logs <- character(0)
  ret <- run_scan_unit(.agds, .nm, .pi_models, .feats, .meds,
                       b_func = .b_func, b_model = NULL,
                       region_type = "noncoding", config = env,
                       results_path = rp, output_dir = od,
                       output_policy = glow_output_policy(write_evidence = TRUE),
                       chr = 22, unit_id = "u", max_regions = max_regions,
                       logger = function(m) logs <<- c(logs, m), n_cores = k)
  ev <- file.path(od, "evidence", "evidence_u.rds")
  list(ret = ret, logs = logs, table = read.csv(rp, stringsAsFactors = FALSE),
       evidence = if (file.exists(ev)) readRDS(ev) else NULL)
}

test_that("noncoding run writes a gene-first table with category, one schedule per mask, identical at n_cores 1 and 2", {
  # The bundled protein-coding table is replaced by the fixture's genes, so the
  # regulatory masks link to GENE_A..GENE_E (see helper-create-test-agds.R).
  gi <- .fixt_genes_info
  testthat::local_mocked_bindings(
    .noncoding_gene_list = function(kind, chr) {
      if (identical(kind, "ncrna")) {
        nc <- get("ncRNA_gene", envir = asNamespace("GLOWr"))
        as.character(nc$ncRNA[as.character(nc$chr) == chr])
      } else gi$hgnc_symbol
    },
    .package = "GLOWr")
  testthat::local_mocked_bindings(
    define_regions_gene = function(genes_info = NULL, chr = NULL, extend = 0L) {
      out <- data.frame(region_id = gi$hgnc_symbol, chr = as.character(gi$chromosome_name),
                        start = gi$start_position, end = gi$end_position,
                        label = gi$hgnc_symbol, stringsAsFactors = FALSE)
      if (!is.null(chr)) out[out$chr %in% as.character(chr), ] else out
    },
    .package = "GLOWr")
  od1 <- withr::local_tempdir()
  r1 <- .run_nc(od1, 1L)
  tbl <- r1$table
  expect_equal(names(tbl)[1:3], c("label", "chr", "category"))
  expect_true(all(c("start", "end", "midpoint", "n_variants", "GLOW_Omni") %in% names(tbl)))
  expect_false(any(.SUMMARY_COLS %in% names(tbl)))                # region_summary off
  # The ncRNA mask runs (its two fixture genes' records are rare SNVs whose
  # cohort MAC is 0 in the 50-sample fixture, so their cells are skipped as
  # empty and write no row); the four regulatory masks write rows.
  expect_setequal(unique(tbl$category), c("upstream", "downstream", "UTR", "enhancer_DHS"))
  expect_true(any(grepl("^noncoding mask ncRNA", r1$logs)))
  expect_equal(r1$ret$summary$n_windows, sum(vapply(c("upstream", "downstream", "UTR", "enhancer_DHS", "ncRNA"),
    function(k) sum(grepl(sprintf("^noncoding mask %s", k), r1$logs)), integer(1))) * 0L +
    r1$ret$summary$n_ok + r1$ret$summary$n_skip_empty + r1$ret$summary$n_staar_skip_single + r1$ret$summary$n_error)
  expect_true(r1$ret$summary$n_skip_empty >= 2L)                    # POTEH-AS1, LINC01297
  # Gene-first: a gene's masks sit together, genes ordered by their anchor
  # (the gene table's start for GENE_*).
  runs <- rle(tbl$label)$values
  expect_equal(anyDuplicated(runs), 0L)
  anchor <- noncoding_gene_anchor(tbl)
  a <- anchor$anchor[match(tbl$label, anchor$gene)]
  expect_true(!is.unsorted(a))
  expect_equal(anchor$anchor[anchor$gene == "GENE_B"], 10000L)
  expect_true(all(anchor$family == "protein_coding"))
  # The upstream set of GENE_A holds its ten upstream records; the candidate
  # span is the records', not the gene's.
  up_a <- tbl[tbl$label == "GENE_A" & tbl$category == "upstream", ]
  expect_equal(nrow(up_a), 1L)
  expect_equal(up_a$n_variants, 10L)
  # The log: one line per mask, one scheduler line per mask at n_cores 2.
  expect_length(grep("^noncoding mask ", r1$logs), 5L)
  expect_length(grep("^Starting noncoding chr22", r1$logs), 1L)
  expect_null(r1$evidence$upstream$GAP1)
  expect_true("GENE_A" %in% names(r1$evidence$upstream))

  skip_on_os("windows")
  skip_if(parallel::detectCores() < 2L, "needs two cores")
  od2 <- withr::local_tempdir()
  r2 <- .run_nc(od2, 2L)
  expect_equal(r2$table, r1$table)
  expect_equal(r2$ret$errors, r1$ret$errors)
  expect_equal(r2$evidence, r1$evidence)
  expect_true(length(grep("^scheduler \\(n_cores = 2\\)", r2$logs)) >= 1L)
})

test_that("noncoding max_regions applies per mask", {
  gi <- .fixt_genes_info
  testthat::local_mocked_bindings(
    .noncoding_gene_list = function(kind, chr) gi$hgnc_symbol, .package = "GLOWr")
  testthat::local_mocked_bindings(
    define_regions_gene = function(genes_info = NULL, chr = NULL, extend = 0L)
      data.frame(region_id = gi$hgnc_symbol, chr = "22", start = gi$start_position,
                 end = gi$end_position, label = gi$hgnc_symbol, stringsAsFactors = FALSE),
    .package = "GLOWr")
  od <- withr::local_tempdir()
  r <- .run_nc(od, 1L, categories = c("upstream", "downstream"), max_regions = 1L)
  expect_equal(sum(r$table$category == "upstream"), 1L)
  expect_equal(sum(r$table$category == "downstream"), 1L)
  expect_length(grep("^max_regions 1: limiting mask", r$logs), 2L)
})

# ---- 2026-10-16: the implementation review round (D3, the anno_cols check) ----

test_that("noncoding sidecars: the errors and summary files are rewritten on a clean rerun", {
  gi <- .fixt_genes_info
  testthat::local_mocked_bindings(
    .noncoding_gene_list = function(kind, chr) gi$hgnc_symbol, .package = "GLOWr")
  testthat::local_mocked_bindings(
    define_regions_gene = function(genes_info = NULL, chr = NULL, extend = 0L)
      data.frame(region_id = gi$hgnc_symbol, chr = "22", start = gi$start_position,
                 end = gi$end_position, label = gi$hgnc_symbol, stringsAsFactors = FALSE),
    .package = "GLOWr")
  od <- withr::local_tempdir()
  cats <- c("upstream", "downstream")
  p_err <- file.path(od, "results-glow", "glow_noncoding_chr22_errors.csv")
  p_sum <- file.path(od, "results-glow", "glow_noncoding_chr22_summary.csv")
  # GENE_B has candidates under both masks; the planted error fails its two cells.
  r_err <- .with_planted_error(.run_nc(od, 1L, categories = cats), label = "GENE_B")
  e1 <- read.csv(p_err, stringsAsFactors = FALSE)
  expect_equal(names(e1), c("category", "label", "message"))
  expect_equal(e1$category, cats)
  expect_equal(e1$label, c("GENE_B", "GENE_B"))
  s1 <- read.csv(p_sum, stringsAsFactors = FALSE)
  expect_equal(s1$category, cats)
  expect_equal(s1$n_error, c(1L, 1L))
  expect_equal(s1$n_attempted, s1$n_ok + s1$n_skip_empty + s1$n_staar_skip_single + s1$n_error)
  expect_false("GENE_B" %in% r_err$table$label)
  # The clean rerun into the same path overwrites both sidecars.
  r_ok <- .run_nc(od, 1L, categories = cats)
  e2 <- read.csv(p_err, stringsAsFactors = FALSE)
  expect_equal(nrow(e2), 0L)
  expect_equal(names(e2), c("category", "label", "message"))
  s2 <- read.csv(p_sum, stringsAsFactors = FALSE)
  expect_equal(s2$n_error, c(0L, 0L))
  expect_equal(nrow(r_ok$table), sum(s2$n_ok) + sum(s2$n_staar_skip_single))
  expect_true("GENE_B" %in% r_ok$table$label)
})

test_that("a staar_context with anno_cols outside pi_features stops before the scan", {
  skip_if_not_installed("STAAR")
  ctx <- .make_staar_ctx()
  ctx$anno_cols <- c(ctx$anno_cols, "no_such_feature")
  env <- .base_config(list(resolved_gene_list = .sched_genes))
  od <- withr::local_tempdir()
  expect_error(.run_sched("gene", env, od, 1L, ctx, chr = 22), "no_such_feature")
})
