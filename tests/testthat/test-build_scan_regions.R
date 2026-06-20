# ==============================================================================
# Tests for build_scan_regions() — region dispatch for the unified scan
# ==============================================================================
#
# build_scan_regions() does not open a GDS; it only constructs the per-unit
# regions table and the variant filter(s) from config, so these tests run on
# plain config lists (no aGDS fixture needed). The coding path calls
# GLOWr::define_regions_gene() internally, mocked here to the fixture genes so
# the returned regions are deterministic and small.


.fs <- list(rare_maf_cutoff = 0.5, variant_type = "SNV",
            min_mac = 1L, min_variants = 2L)


# ---- gene ----

test_that("gene returns a gene table + one filter (bundled genes_info)", {
  env <- new.env(); env$filter_spec <- .fs
  out <- build_scan_regions("gene", chr = 22, config = env)

  expect_true(all(c("regions", "filters") %in% names(out)))
  expect_s3_class(out$regions, "data.frame")
  expect_true(all(c("chr", "start", "end", "label") %in% names(out$regions)))
  expect_true(nrow(out$regions) > 0L)
  expect_true(all(as.character(out$regions$chr) == "22"))
  # One filter object (not a list of filters).
  expect_s3_class(out$filters, "glow_variant_filter")
  expect_equal(out$filters$rare_maf_cutoff, .fs$rare_maf_cutoff)
})

test_that("gene uses config$resolved_gene_list subset to chr (custom gene-list plug-in)", {
  gl <- data.frame(
    region_id = c("G22a", "G22b", "G1x"),
    chr = c("22", "22", "1"),
    start = c(1000L, 5000L, 2000L),
    end = c(2000L, 6000L, 3000L),
    label = c("G22a", "G22b", "G1x"),
    stringsAsFactors = FALSE)
  env <- new.env(); env$filter_spec <- .fs; env$resolved_gene_list <- gl
  out <- build_scan_regions("gene", chr = 22, config = env)

  expect_equal(nrow(out$regions), 2L)               # only the two chr22 rows
  expect_setequal(out$regions$label, c("G22a", "G22b"))
  expect_s3_class(out$filters, "glow_variant_filter")
})


# ---- window ----

test_that("window applies the sliding boundary trim (start <= chunk_end)", {
  env <- new.env()
  env$filter_spec <- .fs; env$window_size <- 2000L; env$step_size <- 1000L
  # Generate to span_end (16000) but keep only windows owned by this chunk
  # (start in [1000, 5000]).
  chunk <- data.frame(chr = "22", chunk_id = 1L,
                      chunk_start = 1000L, chunk_end = 5000L, span_end = 16000L)
  out <- build_scan_regions("window", chr = 22, config = env, chunk = chunk)

  expect_true(nrow(out$regions) > 0L)
  expect_true(all(out$regions$start <= 5000L))      # the boundary rule
  expect_true(all(out$regions$start >= 1000L))
  expect_equal(max(out$regions$start), 5000L)       # last owned window start
  expect_s3_class(out$filters, "glow_variant_filter")
})

test_that("window requires a single-row chunk with the chunk-table columns", {
  env <- new.env()
  env$filter_spec <- .fs; env$window_size <- 2000L; env$step_size <- 1000L
  expect_error(build_scan_regions("window", chr = 22, config = env),
               "requires a single-row")
  bad <- data.frame(chr = "22", chunk_start = 1L)   # missing chunk_end/span_end
  expect_error(build_scan_regions("window", chr = 22, config = env, chunk = bad),
               "missing column")
})


# ---- coding ----

test_that("coding returns the gene table + a NAMED category -> filter map", {
  fixt <- data.frame(region_id = c("GENE_A", "GENE_B"), chr = "22",
                     start = c(1000L, 10000L), end = c(5000L, 15000L),
                     label = c("GENE_A", "GENE_B"), stringsAsFactors = FALSE)
  testthat::local_mocked_bindings(
    define_regions_gene = function(genes_info = NULL, chr = NULL, extend = 0L) fixt,
    .package = "GLOWr")

  env <- new.env()
  env$filter_spec <- .fs
  env$categories <- c("synonymous", "missense")     # builtin names
  out <- build_scan_regions("coding", chr = 22, config = env)

  expect_equal(out$regions$label, c("GENE_A", "GENE_B"))
  # A named map, one filter per category, with the right names.
  expect_type(out$filters, "list")
  expect_named(out$filters, c("synonymous", "missense"))
  expect_true(all(vapply(out$filters,
                         function(f) inherits(f, "glow_variant_filter"),
                         logical(1))))
  # The run's MAF axis flows into every category filter.
  expect_equal(out$filters$synonymous$rare_maf_cutoff, .fs$rare_maf_cutoff)
})

test_that("coding resolves custom DNF clauses and builtin names in one named list", {
  fixt <- data.frame(region_id = "GENE_A", chr = "22", start = 1000L,
                     end = 5000L, label = "GENE_A", stringsAsFactors = FALSE)
  testthat::local_mocked_bindings(
    define_regions_gene = function(genes_info = NULL, chr = NULL, extend = 0L) fixt,
    .package = "GLOWr")

  env <- new.env()
  env$filter_spec <- .fs
  env$categories <- list(
    synonymous = "synonymous",                       # builtin name -> coding_filter
    nonsyn = list(list(genecode_comprehensive_exonic_category = "nonsynonymous SNV"))
  )                                                  # DNF clauses -> variant_filter
  out <- build_scan_regions("coding", chr = 22, config = env)

  expect_named(out$filters, c("synonymous", "nonsyn"))
  # The custom set carries the DNF annotation clauses; the builtin one does not
  # (its mask comes from coding_filter's category, with a category_name attr).
  expect_false(is.null(out$filters$nonsyn$annotation_clauses))
  expect_equal(attr(out$filters$synonymous, "category_name"), "synonymous")
})

test_that("coding requires named categories when given as a list", {
  fixt <- data.frame(region_id = "GENE_A", chr = "22", start = 1L, end = 2L,
                     label = "GENE_A", stringsAsFactors = FALSE)
  testthat::local_mocked_bindings(
    define_regions_gene = function(genes_info = NULL, chr = NULL, extend = 0L) fixt,
    .package = "GLOWr")
  env <- new.env(); env$filter_spec <- .fs
  env$categories <- list("synonymous", "missense")   # unnamed list
  expect_error(build_scan_regions("coding", chr = 22, config = env),
               "must be NAMED")
})


# ---- general validation ----

test_that("filter_spec is required and region_type is validated", {
  env <- new.env()
  expect_error(build_scan_regions("gene", chr = 22, config = env),
               "filter_spec")
  env$filter_spec <- .fs
  expect_error(build_scan_regions("nonsense", chr = 22, config = env))
})
