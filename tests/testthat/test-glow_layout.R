# ==============================================================================
# File: test-glow_layout.R
# Purpose: Unit tests of the derived run layout (R/glow_layout.R): glow_layout()
#   and glow_absolute_paths(), on a small registry written by the tests.

layout_registry_text <- function() c(
  "- {name: dataset, section: top, type: string, default: example, description: The dataset.}",
  "- {name: analysis, section: top, type: string, default: baseline, description: The analysis.}",
  "- {name: output_root, section: top, type: path, default: runs, description: The output root.}",
  "- {name: pheno_csv, section: data, type: path, default: data-example/pheno.csv, description: The phenotype table.}",
  "- {name: favor_match_method, section: shared, type: string, default: flexible, description: The matching method.}",
  "- {name: pheno_path, section: shared, type: path, null_ok: true, derived: true, default: ~, description: The bundle.}",
  "- name: pcs",
  "  section: data_prep",
  "  type: map",
  "  null_ok: true",
  "  fields: {path: {type: path, null_ok: true}, id_col: {type: string}}",
  "  default: {path: ~, id_col: sample.id}",
  "  description: The PC table.",
  "- name: tables",
  "  section: data_prep",
  "  type: map",
  "  values: {type: path}",
  "  default: {}",
  "  description: Extra tables by name.",
  "- name: target",
  "  section: data_prep",
  "  type: one_of",
  "  one_of:",
  "    - {type: number_list}",
  "    - {type: map, required: [file], fields: {file: {type: path}, n: {type: integer}}}",
  "  default: [0.1, 0.2]",
  "  description: A grid or a file.",
  "- {name: region_type, section: snv_set, type: string_list, default: [gene], description: The region types.}")

test_that("glow_layout derives every folder from the three names", {
  d <- withr::local_tempdir()
  s <- gc_setup(d, layout_registry_text())
  cfg <- gc_read(s, c("dataset: cohort0203", "analysis: change_pi", "output_root: /project/lab/runs",
                      "snv_set:", "  region_type: [gene, window]"))
  lay <- glow_layout(cfg)
  base <- "/project/lab/runs/cohort0203"
  expect_identical(lay$output_root, "/project/lab/runs")
  expect_identical(lay$dataset, "cohort0203")
  expect_identical(lay$analysis, "change_pi")
  expect_identical(lay$dataset_dir, base)
  expect_identical(lay$data, list(
    dir = file.path(base, "data"),
    gds = file.path(base, "data/gds"),
    gds_favor = file.path(base, "data/gds_favor"),
    gds_favor_gds = file.path(base, "data/gds_favor/flexible/gds"),
    gds_favor_csv = file.path(base, "data/gds_favor/flexible/csv"),
    pcs = file.path(base, "data/pcs")))
  expect_identical(lay$analysis_dir, file.path(base, "change_pi"))
  expect_identical(lay$steps$data_prep, file.path(base, "change_pi/data_prep"))
  expect_identical(lay$steps$training, file.path(base, "change_pi/training"))
  expect_identical(lay$steps$training_b, file.path(base, "change_pi/training/b"))
  expect_identical(lay$steps$training_pi, file.path(base, "change_pi/training/pi"))
  expect_identical(lay$steps$single_variant, file.path(base, "change_pi/single_variant"))
  expect_identical(lay$steps$snv_set_dir, file.path(base, "change_pi/snv_set"))
  expect_identical(lay$steps$snv_set, c(gene = file.path(base, "change_pi/snv_set/gene"),
                                        window = file.path(base, "change_pi/snv_set/window")))
})

test_that("a relative output_root is read from the root of the configuration", {
  d <- withr::local_tempdir()
  s <- gc_setup(d, layout_registry_text())
  cfg <- gc_read(s, "dataset: mine")
  lay <- glow_layout(cfg)
  expect_identical(lay$output_root, file.path(normalizePath(d), "runs"))
  expect_identical(lay$steps$snv_set, c(gene = file.path(normalizePath(d), "runs/mine/baseline/snv_set/gene")))
  # the function builds paths and creates nothing
  expect_false(dir.exists(file.path(d, "runs")))
})

test_that("the match level and the region types can be given by the caller", {
  d <- withr::local_tempdir()
  s <- gc_setup(d, layout_registry_text())
  cfg <- gc_read(s, "output_root: /r")
  lay <- glow_layout(cfg, match = "exact", region_types = c("coding", "noncoding"))
  expect_identical(lay$data$gds_favor_gds, "/r/example/data/gds_favor/exact/gds")
  expect_identical(names(lay$steps$snv_set), c("coding", "noncoding"))
  expect_identical(unname(glow_layout(cfg, region_types = character(0))$steps$snv_set), character(0))
})

test_that("glow_layout refuses what cannot name a folder", {
  d <- withr::local_tempdir()
  s <- gc_setup(d, layout_registry_text())
  cfg <- gc_read(s, "output_root: /r")
  expect_error(glow_layout(list(dataset = "a")), "read by read_glow_config")
  expect_error(glow_layout(cfg, match = "a/b"), "cannot name a folder")
  expect_error(glow_layout(cfg, match = ".."), "cannot name a folder")
  expect_error(glow_layout(cfg, match = NULL), "`match` has to be one non-empty text value")
  expect_error(glow_layout(cfg, region_types = c("gene", "gene")), "distinct, non-empty names")
  expect_error(glow_layout(cfg, region_types = "x/y"), "cannot name a folder")
  broken <- cfg; broken$dataset <- NULL
  expect_error(glow_layout(broken), "no usable `dataset`")
})

test_that("glow_absolute_paths makes every path absolute, also inside a value with parts", {
  d <- withr::local_tempdir()
  s <- gc_setup(d, layout_registry_text())
  root <- normalizePath(d)
  cfg <- gc_read(s, c("data:", "  pheno_csv: tables/pheno.csv",
                      "data_prep:", "  pcs: {path: pcs/pcs.rds}",
                      "  tables: {a: x/a.csv, b: /abs/b.csv}",
                      "  target: {file: grid.csv, n: 5}"))
  abs <- glow_absolute_paths(cfg)
  expect_identical(abs$output_root, file.path(root, "runs"))
  expect_identical(abs$data$pheno_csv, file.path(root, "tables/pheno.csv"))
  expect_identical(abs$data_prep$pcs$path, file.path(root, "pcs/pcs.rds"))
  expect_identical(abs$data_prep$pcs$id_col, "sample.id")
  expect_identical(abs$data_prep$tables, list(a = file.path(root, "x/a.csv"), b = "/abs/b.csv"))
  expect_identical(abs$data_prep$target$file, file.path(root, "grid.csv"))
  expect_identical(abs$data_prep$target$n, 5L)
  # a null stays null, and nothing but paths changes
  expect_null(abs$shared$pheno_path)
  expect_true("pheno_path" %in% names(abs$shared))
  expect_identical(abs$shared$favor_match_method, "flexible")
  expect_identical(class(abs), class(cfg))
  expect_identical(attr(abs, "glow_root"), attr(cfg, "glow_root"))
  # the other form of the one_of holds no path and is unchanged
  plain <- glow_absolute_paths(gc_read(s, "dataset: other"))
  expect_identical(plain$data_prep$target, c(0.1, 0.2))
  # an absolute configuration is unchanged by a second pass, and it can be written
  expect_identical(glow_absolute_paths(abs), abs)
  out <- withr::local_tempdir()
  write_resolved_config(abs, out)
  back <- read_glow_config(file.path(out, "config-resolved.yml"), base = s$base, registry = s$registry, root = out)
  expect_identical(back$data_prep$pcs$path, file.path(root, "pcs/pcs.rds"))
  expect_error(glow_absolute_paths(list()), "read by read_glow_config")
})

test_that("a path with a drive letter is refused where it is not an absolute path", {
  skip_on_os("windows")
  # (review 2026-11-16) On Linux and macOS `C:/runs` is a relative path, so calling it absolute
  # would leave the output place to the working directory.
  expect_error(GLOWpipeline:::.gl_absolute("C:/runs", "/root"), "starts with a drive letter, which is a Windows path")
  expect_error(GLOWpipeline:::.gl_absolute("d:\\runs", "/root"), "starts with a drive letter")
  expect_identical(GLOWpipeline:::.gl_absolute("/data/runs", "/root"), "/data/runs")
  expect_identical(GLOWpipeline:::.gl_absolute("runs", "/root"), "/root/runs")
  expect_identical(GLOWpipeline:::.gl_absolute("C/runs", "/root"), "/root/C/runs")      # a folder named C
  expect_identical(GLOWpipeline:::.gl_absolute("~/runs", "/root"), path.expand("~/runs"))
})

test_that("a setting that names a file inside a folder of the layout cannot lead out of it", {
  # (review 2026-11-16) the registry type `file_name`
  problem <- GLOWpipeline:::.gc_file_name_problem
  for (ok in c("pheno_covar.rds", "chr{chr}_hg38_favor.gds", "pi_models", "a b.csv", "..hidden"))
    expect_identical(problem(ok), "", info = ok)
  for (bad in c("../x.rds", "sub/x.rds", "/abs/x.rds", "..", ".", "", "a\\b", "~x", "../../other/data_prep/pheno_covar.rds"))
    expect_match(problem(bad), "expects the name of one file or folder", info = bad)
  spec <- list(type = "file_name")
  expect_identical(GLOWpipeline:::.gc_valid("pheno_covar.rds", spec), "")
  expect_match(GLOWpipeline:::.gc_valid("../x.rds", spec), "is not the name of one file or folder")
  expect_match(GLOWpipeline:::.gc_valid(c("a", "b"), spec), "is not a single text value")
  # through the shipped registry: the names of the files that the data preparation and the training write
  skip_without_helper()
  d <- withr::local_tempdir()
  read <- function(lines) {
    f <- file.path(d, "c.yml"); writeLines(lines, f)
    read_glow_config(f, base = file.path(ga_dir(), "base-config", "base_config.yml"),
                     registry = file.path(ga_dir(), "base-config", "settings-registry.yml"), root = ga_dir())
  }
  expect_identical(read(c("data_prep:", "  pheno_out_name: bundle.rds"))$data_prep$pheno_out_name, "bundle.rds")
  expect_error(read(c("data_prep:", "  pheno_out_name: ../../other/data_prep/pheno_covar.rds")),
               "`data_prep.pheno_out_name` expects the name of one file or folder")
  for (line in c("data_prep:\n  gds_pattern: \"../chr{chr}.gds\"", "data_prep:\n  favor_gds_pattern: \"x/chr{chr}.gds\"",
                 "data_prep:\n  favor_csv_pattern: \"/tmp/chr{chr}.csv\"", "training:\n  pi_models_subdir: ../../elsewhere",
                 "shared:\n  analysis_gds_pattern: \"../chr{chr}.gds\""))
    expect_error(read(strsplit(line, "\n")[[1]]), "expects the name of one file or folder", info = line)
  expect_null(read("dataset: mine")$shared$analysis_gds_pattern)            # null stays allowed where it was
})
