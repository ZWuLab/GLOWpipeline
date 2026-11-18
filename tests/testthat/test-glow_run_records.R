# ==============================================================================
# File: test-glow_run_records.R
# Purpose: Unit tests of the run records (R/glow_run_records.R): the table of
#   parts, the text and fingerprint of a dependency, the completion record, the
#   overwrite guard, the lock and the seed record.

deps_of <- function(settings = c("snv_set.ld_threshold" = "0.90000000000000002"),
                    inputs = character(0), files = character(0),
                    versions = c(GFisher = "0.3.2", GLOWr = "0.3.6", GLOWpipeline = "0.3.0", GLOWanalyses = "0.4.0")) {
  list(settings = settings, inputs = inputs, files = files, versions = versions)
}
# Whether the lock of a unit is held now, and the folders of its generations.
lock_base <- function(d, stem) file.path(d, "records", paste0(stem, ".lock"))
lock_held <- function(d, stem) GLOWpipeline:::.rr_lock_state(lock_base(d, stem))$held
lock_folders <- function(d, stem) grep(paste0("^", stem, "\\.lock"), list.files(file.path(d, "records")), value = TRUE)
# One finished unit with one output file.
finish_unit <- function(folder, part = "scan", unit = "chr21", deps = deps_of(), file = "results/chr21.csv", force = FALSE) {
  h <- glow_unit_begin(folder, part, unit, deps, force = force)
  if (h$status != "run") return(h)
  out <- file.path(folder, file)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  writeLines(format(Sys.time(), "%H:%M:%OS6"), out)
  glow_unit_complete(h, out)
  h
}

test_that("the table of parts is read and checked", {
  d <- withr::local_tempdir()
  f <- file.path(d, "parts.yml")
  good <- c(
    "- part: prepare",
    "  step: snv_set",
    "  script: 03-snv-set/prepare.R",
    "  folder: steps.snv_set",
    "  level: region_type",
    "  units: none",
    "  settings: [snv_set.ld_threshold, shared.chroms]",
    "  files: []",
    "  description: The preparation.",
    "- part: scan",
    "  step: snv_set",
    "  script: 03-snv-set/run-gene.R",
    "  folder: steps.snv_set",
    "  level: region_type",
    "  units: region",
    "  forwarded: [snv_set.ld_threshold]",
    "  region_settings:",
    "    gene: [snv_set.max_genes]",
    "    window: [snv_set.max_windows, snv_set.window_size]",
    "  inputs:",
    "    - {part: prepare, setting: snv_set.calibration.report}",
    "    - {part: prepare, unit: \"chr{shared.chroms}\", when: [{setting: snv_set.staar_enabled, is: \"true\"}, {region_type: [coding]}]}")
  writeLines(good, f)
  parts <- glow_read_parts(f)
  expect_identical(names(parts), c("snv_set.prepare", "snv_set.scan"))      # a part is identified with its step
  prepare <- parts[["snv_set.prepare"]]; scan <- parts[["snv_set.scan"]]
  expect_identical(prepare$id, "snv_set.prepare"); expect_identical(prepare$part, "prepare")
  expect_identical(prepare$settings, c("snv_set.ld_threshold", "shared.chroms"))
  expect_identical(prepare$files, list())
  expect_identical(prepare$units, "none")
  expect_identical(prepare$inputs, list())
  expect_identical(prepare$region_settings, list())
  expect_identical(scan$forwarded, "snv_set.ld_threshold")
  expect_identical(scan$units, "region")
  expect_identical(scan$region_settings, list(gene = "snv_set.max_genes", window = c("snv_set.max_windows", "snv_set.window_size")))
  # an input of the same step is named by its part and resolved to step.part
  expect_identical(scan$inputs[[1]], list(part = "snv_set.prepare", setting = "snv_set.calibration.report", unit = NULL,
                                          units_from = NULL, pattern = NULL, partial = NULL, when = list()))
  expect_identical(scan$inputs[[2]]$unit, "chr{shared.chroms}")
  # conditions have one form, in `when` and in `run_when`
  expect_identical(scan$inputs[[2]]$when, list(list(setting = "snv_set.staar_enabled", is = "true"),
                                               list(region_type = "coding")))
  expect_identical(scan$level, "region_type")
  # two steps may have a part of the same name, and an input may name a part of another step
  writeLines(c(good, "- {part: scan, step: single_variant, script: s.R, folder: steps.single_variant, level: analysis, units: chromosome}",
               "- {part: calibrate, step: single_variant, script: c.R, folder: steps.single_variant, level: analysis, units: none,",
               "   inputs: [{part: scan}, {part: snv_set.prepare}]}"), f)
  both <- glow_read_parts(f)
  expect_identical(names(both), c("snv_set.prepare", "snv_set.scan", "single_variant.scan", "single_variant.calibrate"))
  expect_identical(vapply(both[["single_variant.calibrate"]]$inputs, `[[`, "", "part"), c("single_variant.scan", "snv_set.prepare"))

  bad <- function(lines, pattern) { writeLines(lines, f); expect_error(glow_read_parts(f), pattern) }
  bad(sub("units: none", "units: weekly", good), "`units` is one of")
  bad(sub("level: region_type", "level: planet", good), "`level` is one of")
  bad(c(good, "  colour: red"), "unknown key\\(s\\) colour")
  bad(sub("part: scan", "part: prepare", good), "the part `prepare` of the step `snv_set` is named twice")
  bad(sub("\\{part: prepare, setting", "{part: nothing, setting", good), "reads the output of `snv_set.nothing`")
  bad(good[-2], "`step` has to be one text value")
  bad(sub("\\{part: prepare, setting: snv_set.calibration.report\\}", "{part: prepare, colour: red}", good), "an input has to be a mapping with `part`")
  expect_error(glow_read_parts(file.path(d, "absent.yml")), "was not found")

  # with a registry, every named setting has to be registered
  s <- gc_setup(d)
  writeLines(c("- {part: a, step: snv_set, script: x.R, folder: steps.snv_set, level: analysis, units: none,",
               "   settings: [snv_set.ld_threshold, shared.chroms, dataset, snv_set.filter_spec.min_mac]}"), f)
  expect_identical(glow_read_parts(f, s$registry)[["snv_set.a"]]$settings[4], "snv_set.filter_spec.min_mac")
  writeLines("- {part: a, step: s, script: x.R, folder: f, level: analysis, units: none, settings: [snv_set.ld_treshold]}", f)
  expect_error(glow_read_parts(f, s$registry), "not registered: snv_set.ld_treshold")
  # a setting named in an input's `unit`, `when` or in `region_settings` is checked too
  one <- function(extra) writeLines(c("- {part: b, step: s, script: x.R, folder: f, level: analysis, units: none}",
                                      paste0("- {part: a, step: s, script: x.R, folder: f, level: analysis, units: none, ", extra, "}")), f)
  one("inputs: [{part: b, unit: \"chr{shared.chroms}\", when: [{setting: snv_set.staar_enabled, is: \"true\"}]}]")
  expect_length(glow_read_parts(f, s$registry)[["s.a"]]$inputs, 1L)
  one("inputs: [{part: b, unit: \"chr{shared.chrom}\"}]"); expect_error(glow_read_parts(f, s$registry), "not registered: shared.chrom")
  one("inputs: [{part: b, when: [{setting: snv_set.staar, is: \"true\"}]}]"); expect_error(glow_read_parts(f, s$registry), "not registered: snv_set.staar")
  one("region_settings: {gene: [snv_set.max_gene]}"); expect_error(glow_read_parts(f, s$registry), "not registered: snv_set.max_gene")
  # the other places that name a setting: units_from, pattern, partial, a file with its pattern and condition, run_when
  one("inputs: [{part: b, units_from: shared.chroms, pattern: shared.label, partial: snv_set.staar_enabled}]")
  inp <- glow_read_parts(f, s$registry)[["s.a"]]$inputs[[1]]
  expect_identical(inp[c("units_from", "pattern", "partial")],
                   list(units_from = "shared.chroms", pattern = "shared.label", partial = "snv_set.staar_enabled"))
  one("inputs: [{part: b, units_from: shared.chrom}]"); expect_error(glow_read_parts(f, s$registry), "not registered: shared.chrom")
  one("inputs: [{part: b, pattern: shared.patern}]"); expect_error(glow_read_parts(f, s$registry), "not registered: shared.patern")
  one("inputs: [{part: b, partial: snv_set.allow}]"); expect_error(glow_read_parts(f, s$registry), "not registered: snv_set.allow")
  one("files: [shared.gds_dir, {setting: snv_set.alpha, pattern: shared.label, when: [{setting: snv_set.staar_enabled, is_not: \"true\"}]}]")
  files <- glow_read_parts(f, s$registry)[["s.a"]]$files
  expect_identical(files[[1]], list(setting = "shared.gds_dir", pattern = NULL, when = list()))
  expect_identical(files[[2]], list(setting = "snv_set.alpha", pattern = "shared.label",
                                    when = list(list(setting = "snv_set.staar_enabled", is_not = "true"))))
  one("files: [{setting: snv_set.alfa}]"); expect_error(glow_read_parts(f, s$registry), "not registered: snv_set.alfa")
  one("files: [{pattern: snv_set.alpha}]"); expect_error(glow_read_parts(f, s$registry), "an entry of `files` is a setting, or a mapping")
  one("run_when: [{setting: snv_set.alpha, is: set}, {region_type: [gene, window]}]")
  expect_identical(glow_read_parts(f, s$registry)[["s.a"]]$run_when,
                   list(list(setting = "snv_set.alpha", is = "set"), list(region_type = c("gene", "window"))))
  one("run_when: [{setting: snv_set.alpha}]"); expect_error(glow_read_parts(f, s$registry), "a condition of `run_when` is a mapping")
  one("inputs: [{part: b, when: {snv_set.alpha: \"1\"}}]"); expect_error(glow_read_parts(f, s$registry), "a condition of `when` is a mapping")
})

test_that("a setting's text separates values that differ in type, digits or order", {
  expect_false(glow_setting_text(5L) == glow_setting_text(5))
  expect_false(glow_setting_text(0.1 + 0.2) == glow_setting_text(0.3))
  expect_false(glow_setting_text(c("a", "b")) == glow_setting_text(c("b", "a")))
  expect_false(glow_setting_text(list(sex = "sex", age = "age")) == glow_setting_text(list(age = "age", sex = "sex")))
  expect_false(glow_setting_text(NULL) == glow_setting_text("NULL"))
  expect_false(glow_setting_text(c("0" = 0L, "1" = 1L)) == glow_setting_text(c(0L, 1L)))
  expect_identical(glow_setting_text(c(21L, 22L)), glow_setting_text(21:22 + 0L))
  expect_identical(glow_setting_text(NULL), "NULL")
  expect_length(glow_setting_text(as.list(1:400)), 1L)             # one line however long the value
  f <- function(maf) sqrt(maf)
  attr(f, "glow_source") <- "function(maf) sqrt(maf)"
  g <- function(maf) sqrt(maf * 2); attr(g, "glow_source") <- "function(maf) sqrt(maf * 2)"
  expect_identical(glow_setting_text(f), "function text: function(maf) sqrt(maf)")
  expect_false(glow_setting_text(f) == glow_setting_text(g))
})

test_that("a fingerprint follows the content of a file and the listing of a folder", {
  d <- withr::local_tempdir()
  f <- file.path(d, "a.txt"); writeLines("one", f)
  first <- glow_file_fingerprint(f)
  expect_match(first, "^md5:[0-9a-f]{32}$")
  writeLines("two", f)
  expect_false(glow_file_fingerprint(f) == first)
  writeLines("one", f)
  expect_identical(glow_file_fingerprint(f), first)                # the content decides, not the time
  expect_match(glow_file_fingerprint(f, md5_limit = 1), "^size@mtime:4@")   # a file over the limit
  expect_identical(glow_file_fingerprint(file.path(d, "absent")), "missing")
  expect_identical(glow_file_fingerprint(NULL), "none")
  sub <- file.path(d, "tree"); dir.create(file.path(sub, "x"), recursive = TRUE)
  expect_identical(glow_file_fingerprint(sub), "dir:0 files")
  writeLines("a", file.path(sub, "x", "a.csv"))
  one <- glow_file_fingerprint(sub)
  expect_match(one, "^dir:1 files, listing md5:")
  writeLines("b", file.path(sub, "b.csv"))
  expect_match(glow_file_fingerprint(sub), "^dir:2 files")
  expect_false(glow_file_fingerprint(sub) == one)
})

test_that("a fingerprint sees a change within one second, and identifies a set of files", {
  d <- withr::local_tempdir()
  tree <- file.path(d, "tree"); dir.create(tree)
  f <- file.path(tree, "a.txt")
  at <- as.POSIXct("2026-11-16 10:00:00", tz = "UTC")
  # a small file of a folder is identified by its content: same size, same time, other bytes
  writeLines("one", f); Sys.setFileTime(f, at)
  first <- glow_file_fingerprint(tree)
  writeLines("two", f); Sys.setFileTime(f, at)
  expect_false(glow_file_fingerprint(tree) == first)
  # a larger file is identified by size and time, and the time counts below the second
  writeLines("one", f); Sys.setFileTime(f, at + 0.1)
  early <- glow_file_fingerprint(tree, small = 1)
  Sys.setFileTime(f, at + 0.9)
  expect_false(glow_file_fingerprint(tree, small = 1) == early)
  expect_false(glow_file_fingerprint(f, md5_limit = 1) == {Sys.setFileTime(f, at + 0.1); glow_file_fingerprint(f, md5_limit = 1)})
  # several paths are one set: the files that are read together from a folder
  g <- file.path(tree, "b.txt"); writeLines("b", g)
  set <- glow_file_fingerprint(c(f, g))
  expect_match(set, "^set:2 files, md5:[0-9a-f]{32}$")
  expect_identical(glow_file_fingerprint(c(g, f)), set)                # the order does not matter
  writeLines("c", file.path(tree, "c.txt"))                            # a file that is not in the set
  expect_identical(glow_file_fingerprint(c(f, g)), set)
  writeLines("B", g)
  expect_false(glow_file_fingerprint(c(f, g)) == set)
  expect_false(glow_file_fingerprint(c(f, file.path(tree, "absent.txt"))) == glow_file_fingerprint(f))
  expect_identical(glow_file_fingerprint(character(0)), "set:0 files")
})

test_that("a finished unit has a record that reads back, and a second run finds it complete", {
  d <- withr::local_tempdir()
  deps <- deps_of(settings = c("snv_set.ld_threshold" = "0.9", "shared.pi_features" = 'c("a b", "007", "N")',
                               "training.b_func" = "function text: function(maf) {\n  sqrt(maf) # \"quoted\" \\ back\n}"),
                  inputs = c("prepare" = "20261116T101010.123-77-abc"), files = c("data.pheno_csv" = "md5:0123"))
  expect_null(glow_unit_record(d, "scan", "chr21"))
  expect_false(glow_unit_state(d, "scan", "chr21")$complete)
  h <- glow_unit_begin(d, "scan", "chr21", deps)
  expect_identical(h$status, "run")
  expect_true(lock_held(d, "scan.chr21"))
  out <- file.path(d, "results", "chr21.csv"); dir.create(dirname(out)); writeLines("x", out)
  written <- glow_unit_complete(h, out)
  expect_false(lock_held(d, "scan.chr21"))                                  # the lock is released
  rec <- glow_unit_record(d, "scan", "chr21")
  expect_identical(rec$part, "scan"); expect_identical(rec$unit, "chr21")
  expect_identical(rec$outputs, "results/chr21.csv")                        # relative to the folder
  expect_false(rec$seeded)
  expect_identical(rec$generation, written$generation)
  # every dependency reads back as the same text, also with quotes, a line break and a backslash
  for (g in names(deps)) expect_identical(rec$dependencies[[g]], deps[[g]][order(names(deps[[g]]))])
  expect_identical(glow_dependency_differences(rec$dependencies, deps), character(0))
  expect_true(glow_unit_state(d, "scan", "chr21", deps)$complete)
  again <- glow_unit_begin(d, "scan", "chr21", deps)
  expect_identical(again$status, "complete")
  expect_null(again$lock)
  expect_false(lock_held(d, "scan.chr21"))
  # a part that is not split has a record without a unit
  whole <- finish_unit(d, part = "aggregate", unit = NULL, file = "aggregated/all.csv")
  expect_null(glow_unit_record(d, "aggregate")$unit)
  expect_true(file.exists(file.path(d, "records", "aggregate.yml")))
  expect_identical(glow_unit_begin(d, "aggregate", NULL, deps_of())$status, "complete")
})

test_that("a deleted output makes its unit, and only its unit, run again", {
  d <- withr::local_tempdir()
  finish_unit(d, unit = "chr21", file = "results/chr21.csv")
  finish_unit(d, unit = "chr22", file = "results/chr22.csv")
  unlink(file.path(d, "results", "chr22.csv"))
  expect_true(glow_unit_state(d, "scan", "chr21", deps_of())$complete)
  state <- glow_unit_state(d, "scan", "chr22", deps_of())
  expect_false(state$complete)
  expect_match(state$reason, "an output that its record lists is missing: results/chr22.csv")
  expect_identical(glow_unit_begin(d, "scan", "chr21", deps_of())$status, "complete")
  first <- glow_unit_record(d, "scan", "chr22")$generation
  h <- finish_unit(d, unit = "chr22", file = "results/chr22.csv")
  expect_identical(h$status, "run")
  expect_true(glow_unit_state(d, "scan", "chr22", deps_of())$complete)
  expect_false(glow_unit_record(d, "scan", "chr22")$generation == first)   # a new generation mark
})

test_that("a record states which process produced the unit when the caller says so", {
  d <- withr::local_tempdir()
  deps <- deps_of()
  out <- file.path(d, "results", "chr21.csv"); dir.create(dirname(out)); writeLines("x", out)
  # (review 2026-11-17) A caller that started a job for a unit could not tell the record of that
  # job from the record of another run. The record can name its producer.
  h <- glow_unit_begin(d, "scan", "chr21", deps)
  glow_unit_complete(h, out)
  expect_identical(length(glow_unit_record(d, "scan", "chr21")$producer), 0L)      # none stated, none recorded
  expect_false(any(grepl("producer", readLines(file.path(d, "records", "scan.chr21.yml")))))
  h <- glow_unit_begin(d, "scan", "chr21", deps, force = TRUE)
  writeLines("x", out)                                                             # a forced begin retires the old output
  written <- glow_unit_complete(h, out, producer = c(slurm_job = "1234567", host = "node \"a\": 1"))
  rec <- glow_unit_record(d, "scan", "chr21")
  expect_identical(rec$producer, c(host = "node \"a\": 1", slurm_job = "1234567"))
  expect_identical(written$producer, rec$producer)
  expect_true(glow_unit_state(d, "scan", "chr21", deps)$complete)                  # the producer is no dependency
  # a producer that is not a named character vector stops the call, and the lock is released
  h <- glow_unit_begin(d, "scan", "chr21", deps, force = TRUE)
  writeLines("x", out)
  expect_error(glow_unit_complete(h, out, producer = "1234567"), "`producer` is NULL or a named character vector")
  expect_false(lock_held(d, "scan.chr21"))
})

test_that("the overwrite guard refuses other dependencies, names them, and yields to --force", {
  d <- withr::local_tempdir()
  finish_unit(d)
  before <- glow_unit_record(d, "scan", "chr21")
  changed <- deps_of(settings = c("snv_set.ld_threshold" = "0.80000000000000004"))
  expect_error(glow_unit_begin(d, "scan", "chr21", changed),
               "made with other dependencies.*setting snv_set.ld_threshold: recorded 0.9, now 0.8")
  # A number is recorded with 17 digits and shown with 15, in a vector as well. Two values that
  # look alike at 15 digits are shown as they are recorded.
  shown <- function(a, b) glow_dependency_differences(list(settings = c(x = a)), list(settings = c(x = b)), seeded = TRUE)
  expect_identical(shown("0.90000000000000002", "0.69999999999999996"), "setting x: recorded 0.9, now 0.7")
  expect_identical(shown("c(0.90000000000000002, 2.4999999999999998e-06)", "c(0.90000000000000002, 5.0000000000000004e-06)"),
                   "setting x: recorded c(0.9, 2.5e-06), now c(0.9, 5e-06)")
  expect_identical(shown("0.90000000000000002", "0.90000000000000013"),
                   "setting x: recorded 0.90000000000000002, now 0.90000000000000013")
  expect_identical(shown("10L", "11L"), "setting x: recorded 10L, now 11L")
  # (review 2026-11-16, third round) The digits inside a quoted text are not a number of the
  # setting: a path or a name with such digits is shown as it is.
  text <- function(x) glow_setting_text(x)
  expect_identical(shown(text("/cohort/0.90000000000000002/file"), text("/cohort/0.80000000000000004/file")),
                   "setting x: recorded \"/cohort/0.90000000000000002/file\", now \"/cohort/0.80000000000000004/file\"")
  expect_identical(shown(text(c("id-1.00000000000000001", "a \"quoted\" 2.50000000000000001")), text("id-1.00000000000000002")),
                   paste0("setting x: recorded ", text(c("id-1.00000000000000001", "a \"quoted\" 2.50000000000000001")),
                          ", now ", text("id-1.00000000000000002")))
  # a number beside a text in one value: the number is shortened, the text is not
  expect_identical(shown("list(0.90000000000000002, \"0.90000000000000002\")", "list(0.69999999999999996, \"0.90000000000000002\")"),
                   "setting x: recorded list(0.9, \"0.90000000000000002\"), now list(0.7, \"0.90000000000000002\")")
  # (review 2026-11-17) The text of a function is shown as the configuration wrote it. Its source
  # can hold a number inside single quotes, which was shortened like a number of the setting.
  fun <- function(src) { f <- function(maf) 1; attr(f, "glow_source") <- src; glow_setting_text(f) }
  expect_identical(shown(fun("function(maf) as.numeric('0.90000000000000002') * 1.00000000000000001"),
                         fun("function(maf) as.numeric('0.80000000000000004') * 1.00000000000000001")),
                   paste0("setting x: recorded function text: function(maf) as.numeric('0.90000000000000002') * 1.00000000000000001, ",
                          "now function text: function(maf) as.numeric('0.80000000000000004') * 1.00000000000000001"))
  # an input or a file is never reformatted
  expect_identical(glow_dependency_differences(list(inputs = c(p = "1.00000000000000001")), list(inputs = c(p = "1.00000000000000002"))),
                   "input p: recorded 1.00000000000000001, now 1.00000000000000002")
  expect_error(glow_unit_begin(d, "scan", "chr21", changed), "a new `analysis` name")
  expect_error(glow_unit_begin(d, "scan", "chr21", changed, level = "dataset"), "a new `dataset` name")
  # the refusal leaves the record, the results and no lock
  expect_identical(glow_unit_record(d, "scan", "chr21"), before)
  expect_false(lock_held(d, "scan.chr21"))
  # each kind of dependency is compared: a new key, a version, an input, a file
  expect_error(glow_unit_begin(d, "scan", "chr21", deps_of(settings = c("snv_set.ld_threshold" = "0.90000000000000002", "snv_set.max_genes" = "5L"))),
               "setting snv_set.max_genes: recorded \\(absent\\), now 5L")
  expect_error(glow_unit_begin(d, "scan", "chr21", deps_of(versions = c(GFisher = "0.3.2", GLOWr = "0.3.7", GLOWpipeline = "0.3.0", GLOWanalyses = "0.4.0"))),
               "version GLOWr: recorded 0.3.6, now 0.3.7")
  expect_error(glow_unit_begin(d, "scan", "chr21", deps_of(files = c("data.pheno_csv" = "md5:1"))), "file data.pheno_csv")
  expect_error(glow_unit_begin(d, "scan", "chr21", deps_of(settings = c("snv_set.ld_threshold" = "0.8"), inputs = c(prepare = "g2"))),
               "setting snv_set.ld_threshold.*input prepare")
  # the refusal comes first, also when an output of the record is gone: results made with other
  # settings are never written over without --force
  second <- file.path(d, "results", "chr21_errors.csv"); writeLines("e", second)
  unlink(glow_unit_record(d, "scan", "chr21")$outputs)
  expect_true(file.exists(second))
  with_two <- glow_unit_begin(d, "scan", "chr21", deps_of(), force = TRUE)
  writeLines("x", file.path(d, "results", "chr21.csv"))
  glow_unit_complete(with_two, c(file.path(d, "results", "chr21.csv"), second))
  before <- glow_unit_record(d, "scan", "chr21")
  unlink(file.path(d, "results", "chr21.csv"))
  gone <- glow_unit_state(d, "scan", "chr21", changed)
  expect_identical(gone$missing, "results/chr21.csv")
  expect_match(gone$differences, "setting snv_set.ld_threshold")
  expect_match(gone$reason, "other dependencies")
  expect_error(glow_unit_begin(d, "scan", "chr21", changed), "made with other dependencies")
  expect_identical(readLines(second), "e")                                 # the other output is untouched
  expect_identical(glow_unit_record(d, "scan", "chr21"), before)
  # --force removes the record and the files it lists, and writes a new record with a new mark
  h <- glow_unit_begin(d, "scan", "chr21", changed, force = TRUE)
  expect_identical(h$status, "run")
  expect_null(glow_unit_record(d, "scan", "chr21"))
  expect_false(file.exists(second))                                        # no file of the old run stays behind
  writeLines("new", file.path(d, "results", "chr21.csv"))
  glow_unit_complete(h, file.path(d, "results", "chr21.csv"))
  after <- glow_unit_record(d, "scan", "chr21")
  expect_identical(after$dependencies$settings, changed$settings)
  expect_identical(after$outputs, "results/chr21.csv")
  expect_false(after$generation == before$generation)
  # a forced run of a unit that is complete runs again too
  expect_identical(glow_unit_begin(d, "scan", "chr21", changed)$status, "complete")
  forced <- glow_unit_begin(d, "scan", "chr21", changed, force = TRUE)
  expect_identical(forced$status, "run")
  glow_unit_release(forced)
})

test_that("a unit whose input was produced again is out of date and runs again without a question", {
  d <- withr::local_tempdir()
  first <- deps_of(inputs = c(prepare = "generation-1", "annotate_favor.chr21" = "seed-abc"))
  finish_unit(d, deps = first)
  before <- glow_unit_record(d, "scan", "chr21")
  # the same settings, a new generation mark of one input: stale, not a conflict
  now <- deps_of(inputs = c(prepare = "generation-2", "annotate_favor.chr21" = "seed-abc"))
  state <- glow_unit_state(d, "scan", "chr21", now)
  expect_false(state$complete)
  expect_true(state$stale)
  expect_match(state$reason, "a part it read was produced again")
  h <- glow_unit_begin(d, "scan", "chr21", now)
  expect_identical(h$status, "run")
  expect_match(h$stale, "input prepare: recorded generation-1, now generation-2")
  expect_null(glow_unit_record(d, "scan", "chr21"))          # the old record is gone before the unit runs
  expect_false(file.exists(file.path(d, "results", "chr21.csv")))   # and so are its results
  writeLines("again", file.path(d, "results", "chr21.csv"))
  glow_unit_complete(h, file.path(d, "results", "chr21.csv"))
  expect_identical(glow_unit_record(d, "scan", "chr21")$dependencies$inputs[["prepare"]], "generation-2")
  expect_false(glow_unit_record(d, "scan", "chr21")$generation == before$generation)
  # an input that appears or disappears is a change of inputs too
  gone <- glow_unit_state(d, "scan", "chr21", deps_of(inputs = c(prepare = "generation-2")))
  expect_true(gone$stale)
  # a changed setting together with a changed input is a conflict, and so is a setting alone
  expect_false(glow_unit_state(d, "scan", "chr21", deps_of(settings = c("snv_set.ld_threshold" = "0.8"), inputs = c(prepare = "generation-3")))$stale)
  expect_false(glow_unit_state(d, "scan", "chr21", deps_of(settings = c("snv_set.ld_threshold" = "0.8"), inputs = now$inputs))$stale)
  # a complete unit is not stale, and a unit without a record is not stale
  expect_false(glow_unit_state(d, "scan", "chr21", now)$stale)
  expect_false(glow_unit_state(d, "scan", "chr22", now)$stale)
})

test_that("a unit that does not finish leaves no record, and an output has to exist", {
  d <- withr::local_tempdir()
  h <- glow_unit_begin(d, "scan", "chr21", deps_of())
  expect_error(glow_unit_complete(h, file.path(d, "results", "absent.csv")), "did not produce")
  expect_null(glow_unit_record(d, "scan", "chr21"))
  expect_false(lock_held(d, "scan.chr21"))                                  # released on the failure too
  h <- glow_unit_begin(d, "scan", "chr21", deps_of())
  glow_unit_release(h)
  expect_false(lock_held(d, "scan.chr21"))
  expect_null(glow_unit_record(d, "scan", "chr21"))
  expect_error(glow_unit_complete(list(status = "complete")), "identical\\(handle\\$status")
})

test_that("an input that is replaced while a unit runs stops its record", {
  d <- withr::local_tempdir()
  deps <- deps_of(inputs = c(prepare = "generation-1", "annotate_favor.chr21" = "seed-abc"))
  h <- glow_unit_begin(d, "scan", "chr21", deps)
  out <- file.path(d, "chr21.csv"); writeLines("x", out)
  expect_error(glow_unit_complete(h, out, inputs_now = c(prepare = "generation-2", "annotate_favor.chr21" = "seed-abc")),
               "was replaced while it ran.*input prepare: recorded generation-1, now generation-2")
  expect_null(glow_unit_record(d, "scan", "chr21"))
  h <- glow_unit_begin(d, "scan", "chr21", deps)
  glow_unit_complete(h, out, inputs_now = deps$inputs)                        # unchanged inputs: recorded
  expect_identical(glow_unit_record(d, "scan", "chr21")$dependencies$inputs, deps$inputs[order(names(deps$inputs))])
  # an outside file is read again at the end too
  with_file <- deps_of(files = c("data.pheno_csv" = "md5:before"))
  h <- glow_unit_begin(d, "scan", "chr22", with_file)
  out22 <- file.path(d, "chr22.csv"); writeLines("x", out22)
  expect_error(glow_unit_complete(h, out22, files_now = c("data.pheno_csv" = "md5:after")),
               "was replaced while it ran.*file data.pheno_csv: recorded md5:before, now md5:after")
  expect_null(glow_unit_record(d, "scan", "chr22"))
  h <- glow_unit_begin(d, "scan", "chr22", with_file)
  glow_unit_complete(h, out22, inputs_now = character(0), files_now = with_file$files)
  expect_false(is.null(glow_unit_record(d, "scan", "chr22")))
})

test_that("a folder among the outputs is recorded file by file", {
  d <- withr::local_tempdir()
  h <- glow_unit_begin(d, "prepare", NULL, deps_of())
  shared <- file.path(d, "shared"); dir.create(file.path(shared, "sub"), recursive = TRUE)
  writeLines("a", file.path(shared, "01-config_used.rds")); writeLines("b", file.path(shared, "sub", "01-chunks.rds"))
  glow_unit_complete(h, shared)
  expect_identical(glow_unit_record(d, "prepare")$outputs, c("shared/01-config_used.rds", "shared/sub/01-chunks.rds"))
  expect_true(glow_unit_state(d, "prepare", NULL, deps_of())$complete)
  # a file that is deleted inside the folder makes the unit incomplete, although the folder exists
  unlink(file.path(shared, "01-config_used.rds"))
  state <- glow_unit_state(d, "prepare", NULL, deps_of())
  expect_false(state$complete)
  expect_identical(state$missing, "shared/01-config_used.rds")
  # the unit runs again, and the file of the old run that the new run does not write is gone
  h <- glow_unit_begin(d, "prepare", NULL, deps_of())
  expect_identical(h$status, "run")
  expect_false(file.exists(file.path(shared, "sub", "01-chunks.rds")))
  expect_true(dir.exists(file.path(shared, "sub")))                        # folders are left
  writeLines("a", file.path(shared, "01-config_used.rds"))
  glow_unit_complete(h, shared)
  expect_identical(glow_unit_record(d, "prepare")$outputs, "shared/01-config_used.rds")
  # a file outside the unit's folder is recorded by its full path and is not deleted at a rerun
  outside <- file.path(withr::local_tempdir(), "ld_scores.rds"); writeLines("x", outside)
  h <- glow_unit_begin(d, "ld", NULL, deps_of()); glow_unit_complete(h, outside)
  expect_identical(glow_unit_record(d, "ld")$outputs, normalizePath(outside))
  h <- glow_unit_begin(d, "ld", NULL, deps_of(), force = TRUE)
  expect_true(file.exists(outside))
  glow_unit_release(h)
})

test_that("the files of an old record are deleted by their literal names and never through a link", {
  skip_on_os("windows")
  # (review 2026-11-16, second round) A recorded name with a `*` deleted what the pattern matched,
  # and a recorded folder that had become a link led the deletion to another place.
  d <- withr::local_tempdir()
  res <- file.path(d, "results"); dir.create(res)
  writeLines("mine", file.path(res, "table[1]*.csv")); writeLines("a neighbor", file.path(res, "table[1]2.csv"))
  writeLines("another", file.path(res, "table1.csv"))
  h <- glow_unit_begin(d, "scan", "chr21", deps_of()); glow_unit_complete(h, file.path(res, "table[1]*.csv"))
  expect_identical(glow_unit_record(d, "scan", "chr21")$outputs, "results/table[1]*.csv")
  h <- glow_unit_begin(d, "scan", "chr21", deps_of(), force = TRUE); glow_unit_release(h)
  expect_identical(sort(list.files(res)), c("table1.csv", "table[1]2.csv"))       # the one recorded file went, no other
  # a recorded folder is replaced by a link to a folder elsewhere: nothing is deleted there
  models <- file.path(d, "models"); dir.create(models); writeLines("m", file.path(models, "model_1.rds"))
  h <- glow_unit_begin(d, "train", NULL, deps_of()); glow_unit_complete(h, models)
  elsewhere <- withr::local_tempdir(); writeLines("not of this run", file.path(elsewhere, "model_1.rds"))
  unlink(models, recursive = TRUE); file.symlink(elsewhere, models)
  h <- glow_unit_begin(d, "train", NULL, deps_of(), force = TRUE); glow_unit_release(h)
  expect_identical(readLines(file.path(elsewhere, "model_1.rds")), "not of this run")
  # an output that is a link to a file elsewhere: the file it points at stays at a rerun
  target <- file.path(elsewhere, "kept.rds"); writeLines("kept", target)
  link <- file.path(d, "linked.rds"); file.symlink(target, link)
  h <- glow_unit_begin(d, "link", NULL, deps_of()); glow_unit_complete(h, link)
  h <- glow_unit_begin(d, "link", NULL, deps_of(), force = TRUE); glow_unit_release(h)
  expect_identical(readLines(target), "kept")
})

# A lock as another process would have left it: the next generation, held by `owner`.
plant_lock <- function(d, stem, owner) {
  base <- lock_base(d, stem)
  lock <- paste0(base, ".", GLOWpipeline:::.rr_lock_state(base)$n + 1L)
  dir.create(lock, recursive = TRUE)
  writeLines(paste(names(owner), owner, sep = ": "), file.path(lock, "owner"))
  lock
}
dead_pid <- function() { p <- 4194000L; while (dir.exists(file.path("/proc", p))) p <- p + 1L; p }
here <- function() unname(Sys.info()[["nodename"]])
# A `squeue` that answers what the test says, in front of the real one.
local_squeue <- function(lines, status = 0L, env = parent.frame()) {
  bin <- withr::local_tempdir(.local_envir = env)
  script <- file.path(bin, "squeue")
  writeLines(c("#!/bin/sh", sprintf("echo '%s'", lines), sprintf("exit %d", status)), script)
  Sys.chmod(script, "755")
  withr::local_envvar(PATH = paste(bin, Sys.getenv("PATH"), sep = ":"), .local_envir = env)
}

test_that("two runs cannot hold one unit, and the lock of a run that is gone is taken over", {
  d <- withr::local_tempdir()
  withr::local_envvar(SLURM_JOB_ID = "")
  h <- glow_unit_begin(d, "scan", "chr21", deps_of())
  expect_identical(basename(h$lock), "scan.chr21.lock.1")                    # the first generation
  owner <- readLines(file.path(h$lock, "owner"))
  expect_identical(sub(": .*", "", owner), c("host", "pid", "started", "slurm_job", "since", "token"))
  expect_identical(GLOWpipeline:::.rr_owner_field(owner, "token"), h$token)
  # the lock is held by this live process: a second begin is refused, and the token is not shown
  err <- tryCatch(glow_unit_begin(d, "scan", "chr21", deps_of()), error = conditionMessage)
  expect_match(err, "Another run holds `scan.chr21`.*delete the folder .*scan.chr21.lock.1 and run again")
  expect_false(grepl("token", err))
  expect_identical(list.files(file.path(d, "records")), "scan.chr21.lock.1")  # the refused try leaves nothing
  expect_identical(list.files(h$lock), "owner")
  # another unit of the same part is free
  other <- glow_unit_begin(d, "scan", "chr22", deps_of())
  expect_identical(other$status, "run")
  glow_unit_release(other); glow_unit_release(h)
  # a released lock keeps its folder, which holds the number of the generation, and is free
  expect_identical(list.files(h$lock), "released")
  expect_false(lock_held(d, "scan.chr21")); expect_false(lock_held(d, "scan.chr22"))
  again <- glow_unit_begin(d, "scan", "chr21", deps_of())
  expect_identical(basename(again$lock), "scan.chr21.lock.2")                # the next generation, and the one before is gone
  expect_identical(lock_folders(d, "scan.chr21"), "scan.chr21.lock.2")
  glow_unit_release(again)
  skip_if_not(dir.exists("/proc"), "needs /proc to tell a dead process")
  # a lock whose process does not exist on this host is taken over, as the next generation
  planted <- plant_lock(d, "scan.chr21", c(host = here(), pid = dead_pid(), started = "1", slurm_job = "", since = "earlier", token = "t0"))
  expect_identical(basename(planted), "scan.chr21.lock.3")
  took <- glow_unit_begin(d, "scan", "chr21", deps_of())
  expect_identical(took$status, "run")
  expect_identical(basename(took$lock), "scan.chr21.lock.4")
  expect_identical(GLOWpipeline:::.rr_owner_field(took$lock, "pid"), as.character(Sys.getpid()))
  expect_false(dir.exists(planted))
  glow_unit_release(took)
  expect_false(lock_held(d, "scan.chr21"))
  # a process of the same number that started at another time is not the owner
  plant_lock(d, "scan.chr21", c(host = here(), pid = Sys.getpid(), started = "1", slurm_job = "", since = "earlier", token = "t1"))
  reused <- glow_unit_begin(d, "scan", "chr21", deps_of())
  expect_identical(reused$status, "run")
  glow_unit_release(reused)
  # a lock of another host without a job cannot be judged, so it holds
  far <- plant_lock(d, "scan.chr21", c(host = "some-other-node", pid = "1", started = "1", slurm_job = "", since = "earlier", token = "t2"))
  expect_error(glow_unit_begin(d, "scan", "chr21", deps_of()), "Another run holds.*some-other-node")
  expect_true(lock_held(d, "scan.chr21"))
  # a generation without an owner cannot be judged either: it holds until a person removes it
  unlink(file.path(far, "owner"))
  expect_error(glow_unit_begin(d, "scan", "chr21", deps_of()), "its owner is not recorded")
  unlink(far, recursive = TRUE)
  free <- glow_unit_begin(d, "scan", "chr21", deps_of())
  expect_identical(free$status, "run")
  glow_unit_release(free)
})

test_that("a named lock is held by one call at a time", {
  d <- file.path(withr::local_tempdir(), "dataset")                          # the folder is created
  first <- glow_lock(d, "slurm-submission")
  expect_identical(basename(first$lock), "slurm-submission.lock.1")
  expect_error(glow_lock(d, "slurm-submission"), "Another run holds `slurm-submission`")
  started <- Sys.time()
  expect_error(glow_lock(d, "slurm-submission", wait = 2), "Another run holds")       # it looked again before it gave up
  expect_gte(as.numeric(difftime(Sys.time(), started, units = "secs")), 2)
  other <- glow_lock(d, "another-name")                                       # a lock of another name is free
  glow_unit_release(other); glow_unit_release(first)
  second <- glow_lock(d, "slurm-submission")
  expect_identical(basename(second$lock), "slurm-submission.lock.2")
  glow_unit_release(second)
})

test_that("a scheduler that does not answer proves nothing about the owner of a lock", {
  skip_on_os("windows")
  d <- withr::local_tempdir()
  elsewhere <- c(host = "some-other-node", pid = "1", started = "1", slurm_job = "4242", since = "earlier", token = "t0")
  stale <- function(owner) GLOWpipeline:::.rr_lock_is_stale(paste(names(owner), owner, sep = ": "))
  # the job is in the queue: the owner is there
  local({ local_squeue("4242"); expect_false(stale(elsewhere)) })
  # the scheduler answers and does not list the job: the owner is gone
  local({ local_squeue(""); expect_true(stale(elsewhere)) })
  local({ local_squeue("slurm_load_jobs error: Invalid job id specified", 1L); expect_true(stale(elsewhere)) })
  # the scheduler fails for another reason: nothing is known, and the lock holds
  local({
    local_squeue("slurm_load_jobs error: Unable to contact slurm controller (connect failure)", 1L)
    expect_identical(stale(elsewhere), NA)
    held <- plant_lock(d, "scan.chr21", elsewhere)
    expect_error(glow_unit_begin(d, "scan", "chr21", deps_of()), "Another run holds")
    expect_identical(GLOWpipeline:::.rr_owner_field(held, "token"), "t0")
    expect_identical(lock_folders(d, "scan.chr21"), basename(held))
  })
  # on the owner's own host the process decides, whatever the queue says
  skip_if_not(dir.exists("/proc"), "needs /proc to tell a dead process")
  local({
    local_squeue("4242")
    expect_true(stale(c(host = here(), pid = dead_pid(), started = "1", slurm_job = "4242", since = "earlier", token = "t1")))
    expect_false(stale(c(host = here(), pid = Sys.getpid(), started = GLOWpipeline:::.rr_process_start(Sys.getpid()),
                         slurm_job = "4242", since = "earlier", token = "t2")))
  })
  expect_identical(GLOWpipeline:::.rr_lock_is_stale(character(0)), NA)
})

test_that("of two runs that find the same dead owner, one takes the lock, and a late one cannot take it from the new owner", {
  skip_if_not(dir.exists("/proc"), "needs /proc to tell a dead process")
  d <- withr::local_tempdir()
  dir.create(file.path(d, "records"))
  dead <- c(host = here(), pid = dead_pid(), started = "1", slurm_job = "", since = "earlier", token = "t0")
  planted <- plant_lock(d, "scan.chr21", dead)
  first <- GLOWpipeline:::.rr_take_lock(d, "scan", "chr21")
  expect_identical(basename(first$dir), "scan.chr21.lock.2")
  # (review 2026-11-16, second round) A second run that judged the dead owner before the first one
  # took over, and acts only now: it still believes that generation 1 is the current one. It
  # tries to create generation 2, which exists, so it gets nothing, and it changes nothing in
  # the folder of the new owner.
  real_state <- GLOWpipeline:::.rr_lock_state
  looks <- 0L
  testthat::local_mocked_bindings(.rr_lock_state = function(base) {
    looks <<- looks + 1L
    if (looks == 1L) list(n = 1L, dir = planted, held = TRUE, owner = paste(names(dead), dead, sep = ": ")) else real_state(base)
  })
  expect_null(GLOWpipeline:::.rr_take_lock(d, "scan", "chr21", stop_if_held = FALSE))
  expect_gte(looks, 2L)                                                       # it looked again and found the live owner
  expect_identical(list.files(first$dir), "owner")
  expect_identical(GLOWpipeline:::.rr_owner_field(first$dir, "token"), first$token)
  expect_true(GLOWpipeline:::.rr_lock_is_mine(first$dir, first$token))
  expect_identical(lock_folders(d, "scan.chr21"), "scan.chr21.lock.2")        # nothing was left beside it
  looks <- 0L
  expect_error(GLOWpipeline:::.rr_take_lock(d, "scan", "chr21"), "Another run holds")
  GLOWpipeline:::.rr_release_lock(first$dir, first$token)
})

test_that("a run that died while it took or held a lock does not block the next run", {
  skip_if_not(dir.exists("/proc"), "needs /proc to tell a dead process")
  d <- withr::local_tempdir()
  dir.create(file.path(d, "records"))
  dead <- c(host = here(), pid = dead_pid(), started = "1", slurm_job = "", since = "earlier", token = "t0")
  # it died after it built its folder and before the folder was moved into place
  stage <- file.path(d, "records", "scan.chr21.lock.new-t0"); dir.create(stage)
  writeLines(paste(names(dead), dead, sep = ": "), file.path(stage, "owner"))
  a <- glow_unit_begin(d, "scan", "chr21", deps_of())
  expect_identical(a$status, "run"); glow_unit_release(a)
  # it died while it held the lock, and the run after it died as well: each is taken over in turn
  plant_lock(d, "scan.chr21", dead)
  plant_lock(d, "scan.chr21", replace(dead, "token", "t1"))
  b <- glow_unit_begin(d, "scan", "chr21", deps_of())
  expect_identical(b$status, "run")
  expect_identical(length(lock_folders(d, "scan.chr21")[!grepl("new-", lock_folders(d, "scan.chr21"))]), 1L)
  # it died between the two steps of a release: the lock counts as released
  file.create(file.path(b$lock, "released"))
  expect_false(lock_held(d, "scan.chr21"))
  c3 <- glow_unit_begin(d, "scan", "chr21", deps_of())
  expect_identical(c3$status, "run")
  # a handle whose lock went to another run neither writes a record nor releases the other's lock
  out <- file.path(d, "chr21.csv"); writeLines("x", out)
  expect_error(glow_unit_complete(b, out), "no longer held by this run")
  expect_null(glow_unit_record(d, "scan", "chr21"))
  glow_unit_release(b)
  expect_true(lock_held(d, "scan.chr21"))
  expect_true(GLOWpipeline:::.rr_lock_is_mine(c3$lock, c3$token))
  glow_unit_complete(c3, out)
  expect_false(lock_held(d, "scan.chr21"))
  expect_identical(glow_unit_record(d, "scan", "chr21")$outputs, "chr21.csv")
})

test_that("a generation mark is new each time and leaves the random numbers alone", {
  set.seed(42); expected <- runif(3)
  set.seed(42)
  marks <- vapply(1:20, function(i) GLOWpipeline:::.rr_generation(), character(1))
  expect_identical(runif(3), expected)
  expect_identical(anyDuplicated(marks), 0L)
  d <- withr::local_tempdir()
  set.seed(42)
  finish_unit(d)
  expect_identical(runif(3), expected)                     # a whole unit does not touch the generator either
})

test_that("supplied files get a seed record that is compared on the settings they were built with", {
  src <- withr::local_tempdir(); d <- withr::local_tempdir()
  from <- file.path(src, c("model_1.rds", "model_2.rds")); writeLines("m1", from[1]); writeLines("m2", from[2])
  to <- file.path(d, "pi_models", basename(from))
  built_with <- c("training.pi_n_models" = "2L", "shared.pi_features" = 'c("cadd_phred", "linsight")')
  expect_identical(glow_seed_unit(d, "pi_train", NULL, from, to, built_with), "seeded")
  expect_identical(readLines(to[2]), "m2")
  rec <- glow_unit_record(d, "pi_train")
  expect_true(rec$seeded)
  expect_match(rec$generation, "^seed-[0-9a-f]{32}$")
  expect_identical(rec$outputs, c("pi_models/model_1.rds", "pi_models/model_2.rds"))
  expect_identical(rec$dependencies$settings, built_with[order(names(built_with))])
  expect_identical(list.files(file.path(d, "records"), pattern = "\\.yml$"), "pi_train.yml")
  expect_false(lock_held(d, "pi_train"))                                      # no lock stays held
  # a second call finds everything in place and writes nothing
  stamp <- file.mtime(file.path(d, "records", "pi_train.yml"))
  expect_identical(glow_seed_unit(d, "pi_train", NULL, from, to, built_with), "present")
  expect_identical(file.mtime(file.path(d, "records", "pi_train.yml")), stamp)
  # a run with the same settings finds the files complete, whatever its versions, inputs and files
  same <- deps_of(settings = built_with, inputs = c(pi_annotate = "g1"), files = c("x" = "md5:1"))
  expect_true(glow_unit_state(d, "pi_train", NULL, same)$complete)
  expect_identical(glow_unit_begin(d, "pi_train", NULL, same)$status, "complete")
  # a run with other settings is refused, and told that the files were supplied
  other <- deps_of(settings = c("shared.pi_features" = 'c("cadd_phred")', "training.pi_n_models" = "2L"))
  expect_error(glow_unit_begin(d, "pi_train", NULL, other), "setting shared.pi_features.*supplied with the dataset")
  # The record may name more settings than the part has (those of the parts before it). The
  # current text of exactly those is asked for, and a difference in any of them counts.
  asked <- NULL
  with_text <- c(deps_of(settings = c("unrelated.setting" = "1L")),
                 list(setting_text = function(keys) { asked <<- keys; built_with[keys] }))
  expect_true(glow_unit_state(d, "pi_train", NULL, with_text)$complete)
  expect_identical(sort(asked), sort(names(built_with)))
  with_text$setting_text <- function(keys) replace(built_with[keys], "training.pi_n_models", "5L")
  expect_match(glow_unit_state(d, "pi_train", NULL, with_text)$differences, "setting training.pi_n_models: recorded 2L, now 5L")
  # a supplied file that is missing or changed is restored under the seed record
  unlink(to[1]); writeLines("edited", to[2])
  expect_identical(glow_seed_unit(d, "pi_train", NULL, from, to, built_with), "seeded")
  expect_identical(c(readLines(to[1]), readLines(to[2])), c("m1", "m2"))
  expect_identical(glow_unit_record(d, "pi_train")$generation, rec$generation)   # the same supplied bytes
})

test_that("supplied files never replace the results of a run, a held unit, or files of unknown origin", {
  src <- withr::local_tempdir()
  from <- file.path(src, "b_model.rds"); writeLines("supplied", from)
  built_with <- c("training.b_method" = '"both"')
  # a unit that a run produced is left alone, also when its file is deleted afterwards
  d <- withr::local_tempdir(); to <- file.path(d, "b_model.rds")
  h <- glow_unit_begin(d, "b_estimate", NULL, deps_of()); writeLines("trained", to); glow_unit_complete(h, to)
  run_record <- glow_unit_record(d, "b_estimate")
  expect_identical(glow_seed_unit(d, "b_estimate", NULL, from, to, built_with), "kept")
  expect_identical(readLines(to), "trained")
  unlink(to)
  expect_identical(glow_seed_unit(d, "b_estimate", NULL, from, to, built_with), "kept")
  expect_false(file.exists(to))                                    # not restored under the run's record
  expect_identical(glow_unit_record(d, "b_estimate"), run_record)
  expect_false(glow_unit_state(d, "b_estimate", NULL, deps_of())$complete)
  # a unit that a run holds is not written to
  d <- withr::local_tempdir(); to <- file.path(d, "b_model.rds")
  h <- glow_unit_begin(d, "b_estimate", NULL, deps_of())
  expect_identical(glow_seed_unit(d, "b_estimate", NULL, from, to, built_with), "busy")
  expect_false(file.exists(to)); expect_null(glow_unit_record(d, "b_estimate"))
  glow_unit_release(h)
  # a file of unknown origin in the place is neither replaced nor given a record
  d <- withr::local_tempdir(); to <- file.path(d, "b_model.rds"); writeLines("who knows", to)
  expect_identical(glow_seed_unit(d, "b_estimate", NULL, from, to, built_with), "foreign")
  expect_identical(readLines(to), "who knows"); expect_null(glow_unit_record(d, "b_estimate"))
  expect_identical(list.files(file.path(d, "records"), pattern = "\\.yml$"), character(0))
  expect_false(lock_held(d, "b_estimate"))
  # a file in the place that holds the supplied bytes is accepted (a copy that stopped half way)
  writeLines("supplied", to)
  expect_identical(glow_seed_unit(d, "b_estimate", NULL, from, to, built_with), "seeded")
  expect_true(glow_unit_record(d, "b_estimate")$seeded)
  # a forced run replaces the supplied files, and then the seed stays away
  h <- glow_unit_begin(d, "b_estimate", NULL, deps_of(), force = TRUE)
  expect_false(file.exists(to))
  expect_identical(glow_seed_unit(d, "b_estimate", NULL, from, to, built_with), "busy")
  writeLines("trained", to); glow_unit_complete(h, to)
  expect_identical(glow_seed_unit(d, "b_estimate", NULL, from, to, built_with), "kept")
  expect_identical(readLines(to), "trained")
})

test_that("the package versions of a run are named", {
  v <- glow_run_versions("0.4.0")
  expect_identical(names(v), c("GFisher", "GLOWr", "GLOWpipeline", "GLOWanalyses"))
  expect_identical(v[["GLOWanalyses"]], "0.4.0")
  expect_identical(v[["GLOWpipeline"]], as.character(utils::packageVersion("GLOWpipeline")))
})

test_that("a folder name with a wildcard reaches no file beside it", {
  # (review 2026-11-16, third round) unlink() expands *, ? and [ ] in a path. The release of a
  # lock in the folder `case*` removed the owner file of the lock in `case-other`, and the
  # deletion of a unit's files removed those of the other folder.
  root <- withr::local_tempdir()
  for (odd in c("case*", "case?other", "case[-]other", "case[!x]other")) {
    plain <- file.path(root, "case-other"); starred <- file.path(root, odd)
    a <- glow_lock(plain, "protect")
    b <- glow_lock(starred, "protect")
    expect_true(file.exists(file.path(a$lock, "owner")))
    glow_unit_release(b)
    expect_true(file.exists(file.path(a$lock, "owner")), info = odd)       # the other lock is still held
    expect_error(glow_lock(plain, "protect"), "protect")
    glow_unit_release(a)
    expect_false(file.exists(file.path(a$lock, "owner")))
    # a unit in the odd folder: begin, complete, begin again with --force, which retires its output
    for (d in c(plain, starred)) { dir.create(file.path(d, "results"), recursive = TRUE, showWarnings = FALSE); writeLines("kept", file.path(d, "results", "t.csv")) }
    h <- glow_unit_begin(starred, "scan", "chr21", deps_of())
    glow_unit_complete(h, file.path(starred, "results", "t.csv"))
    h <- glow_unit_begin(starred, "scan", "chr21", deps_of(), force = TRUE)
    expect_false(file.exists(file.path(starred, "results", "t.csv")))
    expect_identical(readLines(file.path(plain, "results", "t.csv")), "kept")   # the file of the other folder stays
    glow_unit_release(h)
    unlink(c(plain, starred), recursive = TRUE, expand = FALSE)
  }
})

