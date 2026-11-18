# ==============================================================================
# File: test-glow_config_read.R
# Purpose: Tests of read_glow_config() (R/glow_config_read.R): the chain of files,
#   the merge, the name checks and the study names

# The message of a failed read, or "" when the configuration loads.
gc_error <- function(s, text, name = "user.yml") {
  tryCatch({ gc_read(s, text, name); "" }, error = function(e) conditionMessage(e))
}

test_that("the base alone resolves to the registry defaults", {
  s <- gc_setup(withr::local_tempdir())
  cfg <- read_glow_config(NULL, base = s$base, registry = s$registry, root = s$dir)
  expect_s3_class(cfg, "glow_config")
  expect_identical(cfg$dataset, "example")
  expect_identical(cfg$shared$chroms, 1:22)
  expect_true("pi_features" %in% names(cfg$shared))        # a null setting is present
  expect_null(cfg$shared$pi_features)
  expect_identical(cfg$snv_set$filter_spec$rare_maf_cutoff, 0.5)
  expect_identical(nrow(attr(cfg, "glow_chain")), 1L)
  expect_identical(attr(cfg, "glow_root"), normalizePath(s$dir))
  # what a later step needs from the registry travels with the configuration
  reg <- attr(cfg, "glow_registry")
  expect_true(reg$entries[["shared.gds_dir"]]$derived)
  expect_identical(reg$entries[["snv_set.ld_threshold"]]$spec$type, "number")
  # reading the base as the user file is the base alone
  again <- read_glow_config(s$base, base = s$base, registry = s$registry, root = s$dir)
  expect_identical(nrow(attr(again, "glow_chain")), 1L)
})

test_that("a user file overrides only what it names", {
  s <- gc_setup(withr::local_tempdir())
  cfg <- gc_read(s, c("dataset: cohort0203", "shared:", "  chroms: [21, 22]", "snv_set:", "  ld_threshold: 0.95"))
  expect_identical(cfg$dataset, "cohort0203")
  expect_identical(cfg$analysis, "baseline")
  expect_identical(cfg$shared$chroms, c(21L, 22L))
  expect_identical(cfg$snv_set$ld_threshold, 0.95)
  expect_identical(cfg$snv_set$mac_threshold, 10L)
  expect_identical(nrow(attr(cfg, "glow_chain")), 2L)
})

test_that("a chain of three files resolves, the nearest file winning", {
  s <- gc_setup(withr::local_tempdir())
  dir.create(file.path(s$dir, "sub"))
  gc_file(s, "study.yml", c("dataset: cohort0203", "shared:", "  pi_features: [cadd, linsight, fathmm]",
                            "snv_set:", "  mac_threshold: 11"))
  # `extends` is relative to the file that holds it, and may stand anywhere in the file
  child <- gc_file(s, "sub/change_pi.yml", c("analysis: change_pi", "shared:", "  pi_features: [cadd]",
                                              "extends: ../study.yml"))
  cfg <- read_glow_config(child, base = s$base, registry = s$registry, root = s$dir)
  expect_identical(cfg$dataset, "cohort0203")
  expect_identical(cfg$analysis, "change_pi")
  expect_identical(cfg$shared$pi_features, "cadd")          # a list is replaced whole
  expect_identical(cfg$snv_set$mac_threshold, 11L)
  expect_identical(basename(attr(cfg, "glow_chain")$file), c("base_config.yml", "study.yml", "change_pi.yml"))
  # an absolute path works too
  abs <- gc_file(s, "abs.yml", c(paste0("extends: ", normalizePath(file.path(s$dir, "study.yml"))), "analysis: a2"))
  expect_identical(read_glow_config(abs, base = s$base, registry = s$registry, root = s$dir)$dataset, "cohort0203")
})

test_that("a file may name the base by its bare file name, and by nothing else that is missing", {
  s <- gc_setup(withr::local_tempdir())
  dir.create(file.path(s$dir, "elsewhere"))
  f <- gc_file(s, "elsewhere/u.yml", c("extends: base_config.yml", "dataset: d1"))
  expect_identical(read_glow_config(f, base = s$base, registry = s$registry, root = s$dir)$dataset, "d1")
  # a path that does not exist is an error even when it ends in the base's file name
  for (ext in c("/definitely/absent/base_config.yml", "missing/base_config.yml")) {
    g <- gc_file(s, "elsewhere/v.yml", c(paste("extends:", ext), "dataset: d1"))
    expect_error(read_glow_config(g, base = s$base, registry = s$registry, root = s$dir), "not found.*v.yml", info = ext)
  }
})

test_that("a missing file, a malformed `extends` and a cycle stop at once, with the chain", {
  s <- gc_setup(withr::local_tempdir())
  rd <- function(f) read_glow_config(f, base = s$base, registry = s$registry, root = s$dir)
  expect_error(rd(gc_file(s, "u.yml", "extends: nowhere.yml")), "not found.*u.yml")
  expect_error(rd(gc_file(s, "l.yml", "extends: [a.yml, b.yml]")), "has to name one configuration file")
  expect_error(rd(gc_file(s, "e.yml", "extends: \"\"")), "has to name one configuration file")
  expect_identical(rd(gc_file(s, "n.yml", c("extends: ~", "dataset: d1")))$dataset, "d1")   # a null: builds on the base
  gc_file(s, "a.yml", "extends: b.yml")
  expect_error(rd(gc_file(s, "b.yml", "extends: a.yml")), "cycle")
  expect_error(rd(gc_file(s, "self.yml", "extends: self.yml")), "cycle")
  expect_error(rd(file.path(s$dir, "absent.yml")), "not found")
  # a base that extends another file
  writeLines(c("extends: a.yml", readLines(s$base)), s$base)
  expect_error(rd(NULL), "base configuration cannot extend")
})

test_that("a cycle through a symbolic link is found", {
  s <- gc_setup(withr::local_tempdir())
  a <- gc_file(s, "a.yml", "extends: link.yml")
  skip_if_not(isTRUE(suppressWarnings(file.symlink(a, file.path(s$dir, "link.yml")))), "no symbolic links here")
  expect_error(read_glow_config(a, base = s$base, registry = s$registry, root = s$dir), "cycle")
})

test_that("a map with declared fields merges field by field, in a section and at the top level", {
  s <- gc_setup(withr::local_tempdir())
  cfg <- gc_read(s, c("limits: {b: 5}", "snv_set:", "  filter_spec: {rare_maf_cutoff: 0.01}",
                      "data_prep:", "  covariates: {pc1: PC1}"))
  expect_identical(cfg$snv_set$filter_spec$rare_maf_cutoff, 0.01)
  expect_identical(cfg$snv_set$filter_spec$variant_type, "SNV")       # inherited from the base
  expect_identical(cfg$snv_set$filter_spec$min_mac, 1L)
  expect_identical(cfg$limits, list(a = 1L, b = 5L))                  # a top-level map merges too
  expect_identical(cfg$data_prep$covariates, list(pc1 = "PC1"))       # a free-key map replaces: sex and age are gone
  # two files that set different fields of the same map
  gc_file(s, "p.yml", c("snv_set:", "  filter_spec: {min_mac: 3}"))
  child <- gc_file(s, "c.yml", c("extends: p.yml", "snv_set:", "  filter_spec: {variant_type: Indel}"))
  fs <- read_glow_config(child, base = s$base, registry = s$registry, root = s$dir)$snv_set$filter_spec
  expect_identical(fs[c("rare_maf_cutoff", "variant_type", "min_mac")], list(rare_maf_cutoff = 0.5, variant_type = "Indel", min_mac = 3L))
  # a value that takes one of several forms is replaced whole
  gc_file(s, "g.yml", c("training:", "  b_target_maf: {from: 0.1, to: 0.2, n: 10}"))
  expect_match(tryCatch({ read_glow_config(gc_file(s, "h.yml", c("extends: g.yml", "training:", "  b_target_maf: {n: 20}")),
                                           base = s$base, registry = s$registry, root = s$dir); "" },
                        error = function(e) conditionMessage(e)), "needs the field `from`")
})

test_that("a map marked replace_whole in the registry is replaced whole", {
  s <- gc_setup(withr::local_tempdir())
  expect_identical(read_glow_config(NULL, base = s$base, registry = s$registry, root = s$dir)$data_prep$source,
                   list(col = "outcome"))
  cfg <- gc_read(s, c("data_prep:", "  source: {node: pheno/status}"))
  expect_identical(cfg$data_prep$source, list(node = "pheno/status"))     # `col` is not inherited
})

test_that("an explicit null replaces a value and survives two levels of the chain", {
  s <- gc_setup(withr::local_tempdir())
  gc_file(s, "p.yml", c("shared:", "  pi_features: [cadd]", "data_prep:", "  outcome: {map: ~}"))
  gc_file(s, "c.yml", c("extends: p.yml", "shared:", "  pi_features: ~"))
  g <- gc_file(s, "g.yml", c("extends: c.yml", "dataset: d1"))
  cfg <- read_glow_config(g, base = s$base, registry = s$registry, root = s$dir)
  expect_true("pi_features" %in% names(cfg$shared))
  expect_null(cfg$shared$pi_features)
  expect_true("map" %in% names(cfg$data_prep$outcome))                # a null field of a map
  expect_null(cfg$data_prep$outcome$map)
  expect_identical(cfg$data_prep$outcome$col, "outcome")
})

test_that("values with parts resolve to what the stages need", {
  s <- gc_setup(withr::local_tempdir())
  cfg <- gc_read(s, c(
    "data_prep:",
    "  outcome: {col: status, map: {\"1\": 0, \"2\": 1}}",
    "  covariates:",
    "    sex: {col: SEX, map: {M: 0, F: 1}}",
    "    age: AGE",
    "    pc1: PC1",
    "snv_set:",
    "  categories: {lof: plof, mis: missense}",
    "  b_func: \"function(maf) sqrt(maf * (1 - maf))\"",
    "  filter_spec:",
    "    annotation_clauses:",
    "      - {cadd_phred: {op: gt, value: 20}, category: [exonic, splicing]}",
    "      - {genehancer: {op: nonempty}}",
    "training:",
    "  b_target_maf: {n: 50, to: 0.5, from: 0.001}"))
  expect_identical(cfg$data_prep$outcome$map, c("1" = 0L, "2" = 1L))
  expect_identical(names(cfg$data_prep$covariates), c("sex", "age", "pc1"))     # the order of the design matrix
  expect_identical(cfg$data_prep$covariates$sex, list(col = "SEX", map = c(M = 0, F = 1)))
  expect_identical(cfg$snv_set$categories, list(lof = "plof", mis = "missense"))
  expect_equal(cfg$snv_set$b_func(0.5), 0.5)
  cl <- cfg$snv_set$filter_spec$annotation_clauses
  expect_identical(cl[[1L]], list(cadd_phred = list(op = "gt", value = 20), category = c("exonic", "splicing")))
  expect_identical(cl[[2L]], list(genehancer = list(op = "nonempty")))
  expect_identical(cfg$training$b_target_maf, list(from = 0.001, to = 0.5, n = 50L))   # the declared order of the fields
  cfg2 <- gc_read(s, c("training:", "  b_target_maf: [0.01, 0.05]"), name = "u2.yml")
  expect_identical(cfg2$training$b_target_maf, c(0.01, 0.05))
})

test_that("an unknown name is an error with the nearest registered name", {
  s <- gc_setup(withr::local_tempdir())
  expect_error(gc_read(s, c("snv_set:", "  ld_treshold: 0.8")), "`snv_set.ld_treshold` is not a registered setting; did you mean `ld_threshold`\\?")
  expect_error(gc_read(s, "snvset: {}"), "`snvset` is not a setting or a section; did you mean `snv_set`\\?")
  # a setting written in the wrong section, or at the top level, names where it belongs
  expect_error(gc_read(s, c("shared:", "  staar_enabled: true")), "it is a setting of the section `snv_set`")
  expect_error(gc_read(s, "chroms: 1-22"), "it is a setting of the section `shared`")
  expect_error(gc_read(s, c("snv_set:", "  filter_spec: {min_mc: 2}")), "not a field of `snv_set.filter_spec`; did you mean `min_mac`\\?")
  expect_error(gc_read(s, c("snv_set:", "  zzzzzz: 1")), "`snv_set.zzzzzz` is not a registered setting\\.")   # nothing near: no suggestion
  expect_error(gc_read(s, c("snv_set:", "  alpha: 0.05", "  alpha: 0.01")), "Duplicate map key")
})

test_that("a section is a mapping, or nothing", {
  s <- gc_setup(withr::local_tempdir())
  base <- read_glow_config(NULL, base = s$base, registry = s$registry, root = s$dir)
  same <- function(cfg) expect_identical(unclass(cfg)[names(base)], unclass(base)[names(base)], ignore_attr = TRUE)
  # an empty file, a section with nothing after its colon and an empty mapping set nothing
  for (txt in list("", "# only a comment", "snv_set:", "snv_set: ~", "snv_set: {}")) same(gc_read(s, txt))
  expect_error(gc_read(s, "snv_set: [1, 2]"), "the section `snv_set` has to be a mapping of settings")
  expect_error(gc_read(s, "snv_set: []"), "the section `snv_set` has to be a mapping of settings")
  expect_error(gc_read(s, "snv_set: text"), "the section `snv_set` has to be a mapping of settings")
})

test_that("every problem is reported with its file, all in one message", {
  s <- gc_setup(withr::local_tempdir())
  err <- gc_error(s, c("dataset: ../x", "shared:", "  chroms: 1:22", "snv_set:", "  staar_enabled: yes",
                       "  ld_treshold: 0.8", "  alpha: abc"))
  expect_match(err, "has 5 problem")
  lines <- strsplit(err, "\n")[[1L]][-1L]
  expect_length(lines, 5L)
  expect_true(all(grepl("user.yml: ", lines, fixed = TRUE)))        # each problem carries its origin
  for (p in c("ld_treshold", "write it as a-b", "expects true or false", "expects a number", "folder name"))
    expect_true(any(grepl(p, lines, fixed = TRUE)), info = p)
})

test_that("each file of a chain is held to the registry on its own", {
  s <- gc_setup(withr::local_tempdir())
  gc_file(s, "parent.yml", c("dataset: bad/name", "snv_set:", "  staar_enabled: maybe", "  filter_spec: {min_mac: many}"))
  child <- gc_file(s, "child.yml", c("extends: parent.yml", "dataset: fine", "snv_set:", "  staar_enabled: true",
                                     "  filter_spec: {min_mac: 2}"))
  err <- tryCatch({ read_glow_config(child, base = s$base, registry = s$registry, root = s$dir); "" },
                  error = function(e) conditionMessage(e))
  # the child replaces both values, and the parent's are still reported, in the parent
  lines <- strsplit(err, "\n")[[1L]][-1L]
  expect_length(lines, 2L)
  expect_true(all(grepl("parent.yml: ", lines, fixed = TRUE)))
  expect_true(any(grepl("staar_enabled` expects true or false", lines, fixed = TRUE)))
  expect_true(any(grepl("filter_spec.min_mac` expects a number", lines, fixed = TRUE)))
  # a study name that is in force names the file it came from
  gc_file(s, "p2.yml", "dataset: bad/name")
  err2 <- tryCatch({ read_glow_config(gc_file(s, "c2.yml", c("extends: p2.yml", "analysis: a1")),
                                      base = s$base, registry = s$registry, root = s$dir); "" },
                   error = function(e) conditionMessage(e))
  expect_match(err2, "p2.yml: `dataset` is used as a folder name")
})

test_that("the study names are folder names", {
  s <- gc_setup(withr::local_tempdir())
  expect_error(gc_read(s, "dataset: .."), "folder name")
  expect_error(gc_read(s, "dataset: a/b"), "folder name")
  expect_error(gc_read(s, "analysis: data"), "cannot be `data`")
  expect_identical(gc_read(s, "dataset: cohort02.03_v2-qc")$dataset, "cohort02.03_v2-qc")
  expect_identical(gc_read(s, "dataset: 007", name = "n.yml")$dataset, "007")       # text, as written
})

test_that("a base older than the registry is reported", {
  s <- gc_setup(withr::local_tempdir())
  writeLines(c(gc_test_registry_text(),
               "- {name: new_one, section: snv_set, type: integer, default: 1, description: New.}"), s$registry)
  expect_error(read_glow_config(NULL, base = s$base, registry = s$registry, root = s$dir),
               "`snv_set.new_one` is set by no file")
})
