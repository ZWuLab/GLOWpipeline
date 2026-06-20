# ==============================================================================
# Tests for glow_scan_chunk() — the per-chunk scan engine
# ==============================================================================
#
# Uses the synthetic aGDS (helper-create-test-agds.R, copied from GLOWr) and the
# same prediction-compatible fixtures as GLOWr's test-glow_region_test.R, so the
# real primitive runs end-to-end through the engine. Regions are chosen from the
# fixture's gene layout for reproducible per-window statuses:
#   - GENE_A [1000,5000] / GENE_B [10000,15000] / GENE_C [20000,25000] : "ok"
#   - [1000,1000]        : single post-filter variant -> single / staar_skip_single
#   - [50000,60000] etc. : beyond all genes -> "skip_empty"


# ---- Shared fixtures (built once for the file) ----

.agds <- tempfile(fileext = ".gds")
create_test_agds(.agds)
withr::defer(unlink(.agds), teardown_env())

.feats <- c("cadd_phred", "linsight", "fathmm_xf")
.spec  <- GLOWr::variant_filter(rare_maf_cutoff = 0.5, variant_type = "SNV",
                                min_mac = 1L, min_variants = 1L)

# Transient open to build the null model, sample order, and median cache; the
# engine re-opens the path itself.
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

.make_staar_ctx <- function() {
  sdat  <- data.frame(phenotype = .Y, cov1 = .X[, 2])
  snull <- STAAR::fit_null_glm(phenotype ~ cov1, data = sdat,
                               family = stats::binomial(link = "logit"))
  list(null_model = snull, anno_cols = .feats, rare_maf_cutoff = 0.5)
}

.regions <- function(starts, ends) {
  data.frame(region_id = paste0("w", seq_along(starts)),
             chr = "22", start = as.integer(starts), end = as.integer(ends),
             label = paste0("chr22:", starts, "-", ends),
             stringsAsFactors = FALSE)
}
.reg_ok3   <- .regions(c(1000, 10000, 20000), c(5000, 15000, 25000)) # A,B,C ok
.reg_mix   <- .regions(c(1000, 50000), c(5000, 60000))               # ok + empty
.reg_single<- .regions(1000, 1000)                                   # one variant
.reg_empty <- .regions(c(50000, 70000), c(60000, 80000))             # both empty

# Convenience wrapper binding the fixed fixtures.
.scan <- function(regions, output_policy = glow_output_policy(),
                  staar_context = NULL, null_model = .nm,
                  b_func = .b_func, b_model = NULL, z_scale = 1, use_spa = NULL) {
  glow_scan_chunk(.agds, regions, null_model, .pi_models, .feats, .meds,
                  b_func = b_func, b_model = b_model, filter_spec = .spec,
                  ld_threshold = 0.95, mac_threshold = 11L,
                  collapse_method = "mean", use_spa = use_spa, z_scale = z_scale,
                  output_policy = output_policy,
                  staar_context = staar_context, verbose = 0)
}

.SCHEMA <- c("label", "chr", "start", "end", "midpoint",
             "n_variants", "n_after_annotation", "n_after_filter",
             "n_after_ld", "n_after_collapse", "cMAC", "is_single",
             GLOWpipeline:::.GLOW_PVALUE_COLS, GLOWpipeline:::.STAAR_FAMILY_COLS)


# ---- Step-0 validation ----

test_that("write_staar_detail without a staar_context errors before opening the GDS", {
  expect_error(
    .scan(.reg_mix, output_policy = glow_output_policy(write_staar_detail = TRUE),
          staar_context = NULL),
    "write_staar_detail")
})

test_that("supplying both or neither B source is an error", {
  expect_error(.scan(.reg_mix, b_func = .b_func, b_model = 1), "Exactly one")
  expect_error(.scan(.reg_mix, b_func = NULL,    b_model = NULL), "Exactly one")
})

test_that("a non-policy object is rejected", {
  expect_error(.scan(.reg_mix, output_policy = list(write_flat_table = TRUE)),
               "glow_output_policy")
})

test_that("multi-chromosome regions are rejected (single-chr median cache)", {
  bad <- rbind(.reg_ok3, .regions(1000, 5000))
  bad$chr[nrow(bad)] <- "21"
  expect_error(.scan(bad), "multiple chromosomes")
})


# ---- Accumulation ----

test_that("ok + skip_empty: one row, correct summary tallies, STAAR cols all NA", {
  out <- .scan(.reg_mix)                                  # STAAR off

  expect_equal(out$summary$n_windows, 2L)
  expect_equal(out$summary$n_ok, 1L)
  expect_equal(out$summary$n_skip_empty, 1L)
  expect_equal(out$summary$n_staar_skip_single, 0L)
  expect_equal(out$summary$n_error, 0L)

  expect_equal(nrow(out$results), 1L)
  expect_false(out$results$is_single)                    # GENE_A is multi-variant
  expect_true(all(is.na(out$results[, GLOWpipeline:::.STAAR_FAMILY_COLS])))
  expect_null(out$evidence)                              # lean policy
  expect_equal(nrow(out$errors), 0L)
})

test_that("all windows error (malformed null model): 0-row results, errors recorded", {
  nm_bad <- .nm
  nm_bad$n <- 999L                                       # dimension mismatch downstream
  out <- .scan(.reg_ok3, null_model = nm_bad)

  expect_equal(out$summary$n_error, 3L)
  expect_equal(out$summary$n_ok, 0L)
  expect_equal(nrow(out$results), 0L)                    # 0-row, full schema
  expect_equal(names(out$results), .SCHEMA)
  expect_equal(nrow(out$errors), 3L)
  expect_true(all(nzchar(out$errors$message)))
})

test_that("empty chunk: 0-row results carrying the full schema", {
  out <- .scan(.reg_empty)
  expect_equal(out$summary$n_skip_empty, 2L)
  expect_equal(out$summary$n_ok, 0L)
  expect_equal(nrow(out$results), 0L)
  expect_equal(names(out$results), .SCHEMA)
})


# ---- Output schema (the always-on flat table) ----

test_that("flat-table schema, midpoint, and GLOW/STAAR column drift guard", {
  out <- .scan(.reg_ok3)
  expect_equal(nrow(out$results), 3L)
  expect_equal(names(out$results), .SCHEMA)

  # midpoint is the integer window centre.
  expect_equal(out$results$midpoint,
               as.integer((out$results$start + out$results$end) / 2))

  # The live primitive emits exactly the GLOW columns the engine hard-codes for
  # the empty-chunk schema (drift guard).
  glow_live <- grep("^GLOW_", names(out$results), value = TRUE)
  expect_setequal(glow_live, GLOWpipeline:::.GLOW_PVALUE_COLS)

  # GLOW p-values are valid probabilities.
  pv <- unlist(out$results[, GLOWpipeline:::.GLOW_PVALUE_COLS])
  expect_true(all(pv >= 0 & pv <= 1))
})


# ---- Evidence gating ----

test_that("evidence is gated by the policy and keyed by window label", {
  expect_null(.scan(.reg_ok3)$evidence)                  # off by default

  out <- .scan(.reg_ok3, output_policy = glow_output_policy(write_evidence = TRUE))
  expect_type(out$evidence, "list")
  expect_length(out$evidence, 3L)
  expect_setequal(names(out$evidence), .reg_ok3$label)
  expect_true(all(c("region", "post_filter", "test_units") %in%
                  names(out$evidence[[1]])))
})


# ---- STAAR-on paths (skipped if STAAR is unavailable) ----

test_that("single-variant window with STAAR: staar_skip_single, NA STAAR block, is_single", {
  skip_if_not_installed("STAAR")
  out <- .scan(.reg_single, staar_context = .make_staar_ctx())

  expect_equal(out$summary$n_staar_skip_single, 1L)
  expect_equal(nrow(out$results), 1L)
  expect_true(out$results$is_single)
  expect_true(all(is.na(out$results[, GLOWpipeline:::.STAAR_FAMILY_COLS])))
})

# ---- Progress hook ----

test_that("progress_every fires the callback at the requested cadence", {
  buf <- character(0)
  cb  <- function(msg) buf <<- c(buf, msg)

  # 3 windows, progress_every = 1 -> 3 lines.
  .scan(.reg_ok3, output_policy = glow_output_policy())  # warm-up only
  buf <- character(0)
  glow_scan_chunk(.agds, .reg_ok3, .nm, .pi_models, .feats, .meds,
                  b_func = .b_func, filter_spec = .spec,
                  ld_threshold = 0.95, mac_threshold = 11L,
                  collapse_method = "mean",
                  output_policy = glow_output_policy(),
                  progress_every = 1L, progress_callback = cb, verbose = 0)
  expect_length(buf, 3L)
  expect_match(buf, "^progress: \\d+/3 windows$", all = TRUE)

  # progress_every = 2 over 3 windows -> 1 line (after window 2).
  buf <- character(0)
  glow_scan_chunk(.agds, .reg_ok3, .nm, .pi_models, .feats, .meds,
                  b_func = .b_func, filter_spec = .spec,
                  ld_threshold = 0.95, mac_threshold = 11L,
                  collapse_method = "mean",
                  output_policy = glow_output_policy(),
                  progress_every = 2L, progress_callback = cb, verbose = 0)
  expect_length(buf, 1L)
  expect_match(buf, "^progress: 2/3 windows$")

  # Default (off) emits nothing.
  buf <- character(0)
  glow_scan_chunk(.agds, .reg_ok3, .nm, .pi_models, .feats, .meds,
                  b_func = .b_func, filter_spec = .spec,
                  ld_threshold = 0.95, mac_threshold = 11L,
                  collapse_method = "mean",
                  output_policy = glow_output_policy(),
                  progress_callback = cb, verbose = 0)
  expect_length(buf, 0L)
})


test_that("STAAR on: family p-values populate the block and staar_detail is gated", {
  skip_if_not_installed("STAAR")
  ctx <- .make_staar_ctx()

  out_lean <- .scan(.reg_ok3, staar_context = ctx)        # STAAR runs, no detail
  expect_false(any(is.na(out_lean$results$STAAR_O)))      # 3 multi-variant -> non-NA
  expect_null(out_lean$staar_detail)

  out_det <- .scan(.reg_ok3,
                   output_policy = glow_output_policy(write_staar_detail = TRUE),
                   staar_context = ctx)
  expect_type(out_det$staar_detail, "list")
  expect_length(out_det$staar_detail, 3L)
  expect_setequal(names(out_det$staar_detail), .reg_ok3$label)
})


# ---- M2 z_scale passthrough ----

test_that("z_scale is threaded to the primitive (default 1 is a no-op)", {
  base <- .scan(.reg_ok3)                    # default z_scale = 1
  z1   <- .scan(.reg_ok3, z_scale = 1)
  zc   <- .scan(.reg_ok3, z_scale = 0.5)     # halve every unit's Z

  gcols <- GLOWpipeline:::.GLOW_PVALUE_COLS
  # default == explicit 1 (no-op)
  expect_equal(z1$results[, gcols], base$results[, gcols])
  # z_scale reaches glow_region_test(): shrinking Z changes the GLOW p-values,
  # and (as a uniform shrink) makes the omnibus no more significant.
  expect_false(isTRUE(all.equal(zc$results[, gcols], base$results[, gcols])))
  expect_true(all(zc$results$GLOW_Omni >= base$results$GLOW_Omni - 1e-12))
})

test_that("invalid z_scale is rejected before opening the GDS", {
  expect_error(.scan(.reg_ok3, z_scale = 0),       "z_scale")
  expect_error(.scan(.reg_ok3, z_scale = -2),      "z_scale")
  expect_error(.scan(.reg_ok3, z_scale = c(1, 2)), "z_scale")
})


# ---- use_spa passthrough ----

test_that("use_spa is threaded to the primitive (default NULL = SPA for binary)", {
  base <- .scan(.reg_ok3)                     # default NULL -> auto SPA (binary)
  spa  <- .scan(.reg_ok3, use_spa = TRUE)     # explicit SPA
  std  <- .scan(.reg_ok3, use_spa = FALSE)    # standard score stats

  gcols <- GLOWpipeline:::.GLOW_PVALUE_COLS
  # binary auto-default IS SPA -> NULL == TRUE (no-op default)
  expect_equal(spa$results[, gcols], base$results[, gcols])
  # use_spa reaches glow_region_test(): the non-SPA dispatch changes GLOW p-values
  expect_false(isTRUE(all.equal(std$results[, gcols], base$results[, gcols])))
})

test_that("use_spa = TRUE on a continuous trait errors before the loop", {
  set.seed(101)
  nm_cont <- GLOWr::fit_null_model(.X, rnorm(.n), trait = "continuous",
                                   sample_id = .sids)
  expect_error(.scan(.reg_ok3, null_model = nm_cont, use_spa = TRUE), "binary")
})

test_that("invalid use_spa is rejected before opening the GDS", {
  expect_error(.scan(.reg_ok3, use_spa = "yes"),          "use_spa")
  expect_error(.scan(.reg_ok3, use_spa = c(TRUE, FALSE)), "use_spa")
  expect_error(.scan(.reg_ok3, use_spa = NA),             "use_spa")
})
