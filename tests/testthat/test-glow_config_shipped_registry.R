# ==============================================================================
# File: test-glow_config_shipped_registry.R
# Purpose: Tests of the loader on the settings registry that GLOWanalyses ships
#   (packages/GLOWanalyses/base-config/settings-registry.yml), against the functions
#   of GLOWr that consume the values. They run where the registry is next to the
#   package, in the GLOWpipeline monorepo, and are skipped elsewhere, because the
#   registry is not part of this package.
#
#   The cases are the ones a type table cannot catch on its own: a value that loads
#   without an error and then selects other records or variants than the R
#   configuration it replaces

shipped <- function() {
  dir <- testthat::test_path("..", "..", "..", "GLOWanalyses", "base-config")
  list(registry = file.path(dir, "settings-registry.yml"), base = file.path(dir, "base_config.yml"))
}
skip_without_registry <- function() {
  testthat::skip_if(in_release_tree(), "a release tree is tested without GLOWanalyses")   # helper-glow-stage.R
  p <- shipped()
  testthat::skip_if_not(file.exists(p$registry) && file.exists(p$base), "the GLOWanalyses registry is not next to the package")
}
# Read a configuration given as text against the shipped registry.
shipped_read <- function(text, dir) {
  p <- shipped()
  f <- file.path(dir, "user.yml")
  writeLines(text, f)
  read_glow_config(f, base = p$base, registry = p$registry, root = dir)
}

test_that("the shipped base is current, loads, and round-trips", {
  skip_without_registry()
  p <- shipped()
  expect_identical(write_base_config(p$registry, NULL), readLines(p$base, encoding = "UTF-8"))   # the base is not stale
  d <- withr::local_tempdir()
  cfg <- read_glow_config(NULL, base = p$base, registry = p$registry, root = d)
  expect_identical(cfg$dataset, "example")
  expect_identical(cfg$shared$chroms, 1:22)
  back <- read_glow_config(write_resolved_config(cfg, d), base = p$base, registry = p$registry, root = d)
  expect_identical(unclass(back)[names(cfg)], unclass(cfg)[names(cfg)], ignore_attr = TRUE)
})

test_that("text that names a column or a file prefix stays text", {
  skip_without_registry()
  cfg <- shipped_read(c(
    "training:",
    "  b_qc_filters: {duplicate_priority: \"007\"}",
    "  pi_qc_filters: {duplicate_priority: \"007\", remove_na_cols: [rsID, \"010\"]}",
    "  b_filter_include: {N: [1e3], study: [\"007\"]}",
    "data_prep:",
    "  plink_to_gds_opts: {output_prefix: \"007\", split_by_chr: true, chromosomes: [21, 22]}"), withr::local_tempdir())
  expect_identical(cfg$training$b_qc_filters, list(duplicate_priority = "007"))     # a list of your own, whole
  expect_identical(cfg$training$pi_qc_filters$duplicate_priority, "007")
  expect_identical(cfg$training$pi_qc_filters$remove_na_cols, c("rsID", "010"))
  # filter targets are kept as written; the stage converts them for a numeric column
  expect_identical(cfg$training$b_filter_include, list(N = "1e3", study = "007"))
  expect_identical(cfg$data_prep$plink_to_gds_opts,
                   list(split_by_chr = TRUE, chromosomes = c("21", "22"), output_prefix = "007"))
})

test_that("an annotation test loads as the predicate GLOWr builds from the R configuration", {
  skip_without_registry()
  cfg <- shipped_read(c(
    "snv_set:",
    "  filter_spec:",
    "    annotation_clauses:",
    "      - marker: {op: in, value: \"007\"}",
    "        cadd_phred: {op: lt, value: .inf}",
    "        score: {op: ge, value: 20}",
    "        genehancer: {op: nonempty}",
    "        category: [exonic, splicing]",
    "  categories:",
    "    lof: plof",
    "    rare_missense: [{category: [missense], metasvm: {op: in, value: [D]}}]"), withr::local_tempdir())
  cl <- cfg$snv_set$filter_spec$annotation_clauses[[1L]]
  pred <- function(p) GLOWr::annotation_predicate(p$op, p$value)
  expect_identical(pred(cl$marker), GLOWr::annotation_predicate("in", "007"))           # not the number 7
  expect_identical(pred(cl$cadd_phred), GLOWr::annotation_predicate("lt", Inf))
  expect_identical(pred(cl$score), GLOWr::annotation_predicate("ge", 20))
  expect_identical(pred(cl$genehancer), GLOWr::annotation_predicate("nonempty"))
  expect_identical(cl$category, c("exonic", "splicing"))
  expect_identical(cfg$snv_set$categories$lof, "plof")
  expect_identical(pred(cfg$snv_set$categories$rare_missense[[1L]]$metasvm), GLOWr::annotation_predicate("in", "D"))
  # a test without what its operator needs is refused by the loader
  d <- withr::local_tempdir()
  for (bad in c("{op: gt, value: high}", "{op: gt}", "{op: nonempty, value: 1}", "{value: 5}", "{op: between, value: 1}"))
    expect_error(shipped_read(c("snv_set:", "  filter_spec:", "    annotation_clauses:", paste0("      - x: ", bad)), d),
                 "annotation_clauses", info = bad)
})

test_that("the nulls that the stages accept load", {
  skip_without_registry()
  cfg <- shipped_read(c(
    "data_prep:",
    "  covariates: ~",
    "  pcs: {id_col: ~}",
    "training:",
    "  b_trait_type: ~",
    "  b_training_trait: ~",
    "  b_target_maf: {from: 0.1, to: 0.3, n: ~}"), withr::local_tempdir())
  expect_null(cfg$data_prep$covariates)
  expect_null(cfg$data_prep$pcs$id_col)
  expect_null(cfg$training$b_trait_type)
  expect_null(cfg$training$b_training_trait)
  expect_identical(cfg$training$b_target_maf, list(from = 0.1, to = 0.3, n = NULL))
  expect_null(shipped_read(c("data_prep:", "  pcs: ~"), withr::local_tempdir())$data_prep$pcs)
})

test_that("every form of the shipped settings that take several reads back unchanged", {
  skip_without_registry()
  p <- shipped()
  d <- withr::local_tempdir()
  forms <- list(
    c("snv_set:", "  categories: [plof, missense]"),
    c("snv_set:", "  categories: {lof: plof, custom: [{category: [missense]}]}"),
    c("training:", "  b_target_maf: [0.01, 0.05]"),
    c("training:", "  b_target_maf: {from: 0.001, to: 0.5, n: 50}"),
    c("data_prep:", "  outcome: {col: y, map: {\"1\": 0, \"2\": 1}}"),
    c("data_prep:", "  outcome: {col: y, map: {a: 0.5, b: 1}}"),
    c("data_prep:", "  outcome: {node: sample.annotation/status, map: \"function(v) as.numeric(v == 2)\"}"),
    c("data_prep:", "  covariates: {sex: {col: SEX, map: {M: 0, F: 1}}, age: AGE}"),
    c("data_prep:", "  plink_to_gds_opts: {regions: {chr: \"21\", start: 1000000, end: 5000000}, parallel: 4}"),
    c("data_prep:", "  plink_to_gds_opts: {regions: [{chr: \"21\", start: 1, end: 2}, {chr: X, start: 3, end: 4}], parallel: true}"),
    c("snv_set:", "  filter_spec:", "    annotation_clauses:", "      - {a: {op: not_in, value: [x, y]}, b: {op: le, value: 0.5}, c: {op: empty}}"))
  strip <- function(x) if (is.function(x)) attr(x, "glow_source") else if (is.list(x)) lapply(x, strip) else x
  for (txt in forms) {
    cfg <- shipped_read(txt, d)
    back <- read_glow_config(write_resolved_config(cfg, d), base = p$base, registry = p$registry, root = d)
    expect_identical(strip(unclass(back)[names(cfg)]), strip(unclass(cfg)[names(cfg)]), info = paste(txt, collapse = " "))
  }
})
