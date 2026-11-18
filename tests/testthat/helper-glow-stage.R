# ==============================================================================
# File: helper-glow-stage.R
# Purpose: What the tests of the GLOWanalyses stage helper and of the `glow` command
#   share (test-glow_stage_helper.R, test-glow_run_command.R): where GLOWanalyses is,
#   how a stage script starts its unit in a test, and small study configurations.
#   The tests run where GLOWanalyses is next to the package, in the GLOWpipeline
#   monorepo, and are skipped elsewhere.

ga_dir <- function() normalizePath(testthat::test_path("..", "..", "..", "GLOWanalyses"), mustWork = FALSE)
# Whether the package under test is a release tree (the release build writes .release-source
# into it). The folder next to a release tree is the GLOWanalyses tree that was built last,
# which need not be the templates of this code, so the tests that need GLOWanalyses run in
# the monorepo only.
in_release_tree <- function() file.exists(testthat::test_path("..", "..", ".release-source"))
skip_without_helper <- function() {
  testthat::skip_if(in_release_tree(), "a release tree is tested without GLOWanalyses")
  testthat::skip_if_not(file.exists(file.path(ga_dir(), "lib", "glow-stage.R")) &&
                          file.exists(file.path(ga_dir(), "base-config", "parts.yml")),
                        "GLOWanalyses is not next to the package")
}
# The helper's functions, sourced into an environment of their own.
helper_env <- function() {
  env <- new.env(parent = globalenv())
  sys.source(file.path(ga_dir(), "lib", "glow-stage.R"), envir = env)
  env
}
# The options and the part-by-flag rule of each script that the tests start.
script_of <- list(
  "data_prep.plink_to_gds"   = list(options = c(chr = "value")),
  "data_prep.annotate_favor" = list(options = c(chr = "value")),
  "single_variant.scan"      = list(options = c(chr = "value", combine = "flag"), part_if = c(combine = "scan_combine")),
  "single_variant.ld_scores" = list(options = c(chr = "value", combine = "flag"), part_if = c(combine = "ld_scores_combine")),
  "snv_set.scan"             = list(options = c(chr = "value", chunk = "value", "max-genes" = "value", "n-cores" = "value")))
# Run glow_stage_begin() as a stage script would: from the GLOWanalyses directory.
# A unit with nothing to do returns its context with `skipped` in place of quitting.
# `script_file` is the file of the stage script that the call stands for. NA makes no check of
# it (the test process is not a stage script).
begin <- function(h, step, part, args, load = TRUE, script_file = NA_character_) {
  script <- new.env(parent = globalenv())
  withr::local_dir(ga_dir())
  withr::local_options(glow.stage.no_quit = TRUE)
  how <- script_of[[paste(step, part, sep = ".")]]
  st <- h$glow_stage_begin(step, part, options = if (is.null(how)) character(0) else how$options, load = load,
                           part_if = how$part_if, envir = script, args = args, script = script_file)
  list(st = st, env = script)
}
# Start a unit, let it write one file and finish it. Returns the context of the begin.
finish <- function(h, step, part, args, file = "out.txt", script_file = NA_character_) {
  got <- begin(h, step, part, args, load = FALSE, script_file = script_file)
  if (!is.null(got$st$skipped)) return(got)
  folder <- dirname(dirname(got$st$record$handle$record_file))
  out <- file.path(folder, paste0(got$st$part, if (!is.null(got$st$unit)) paste0(".", got$st$unit), ".", file))
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  writeLines("result", out)
  withr::with_dir(ga_dir(), h$glow_stage_end(got$st, outputs = out))
  got
}
# A study of its own: a dataset that is not the example, with raw inputs and outputs in `out`.
# `extra` adds lines to a section, as list(section = lines).
study <- function(out, extra = list(), name = "study.yml") {
  for (d in c("plink", "favor")) dir.create(file.path(out, d), showWarnings = FALSE)
  for (f in c("pheno.csv", "known.csv", file.path("plink", outer(c("chr21", "chr22"), c(".bed", ".bim", ".fam"), paste0))))
    if (!file.exists(file.path(out, f))) writeLines("x", file.path(out, f))
  sections <- list(
    data = c(paste0("  plink_dir: ", file.path(out, "plink")), paste0("  pheno_csv: ", file.path(out, "pheno.csv")),
             paste0("  favor_db: ", file.path(out, "favor"))),
    shared = "  chroms: [21, 22]",
    data_prep = "  pcs: ~",
    training = c(paste0("  b_raw_path: ", file.path(out, "known.csv")), paste0("  pi_raw_path: ", file.path(out, "known.csv")),
                 "  pi_chromosomes: [21, 22]"))
  key_of <- function(lines) sub("^  ([A-Za-z_.]+):.*$", "\\1", lines)
  # a line of `extra` replaces the line of the same setting
  for (section in names(extra))
    sections[[section]] <- c(sections[[section]][!key_of(sections[[section]]) %in% key_of(extra[[section]])], extra[[section]])
  f <- file.path(out, name)
  writeLines(c("dataset: mine", paste0("output_root: ", file.path(out, "runs")),
               unlist(lapply(names(sections), function(k) c(paste0(k, ":"), sections[[k]])))), f)
  f
}
# Give back the lock of a unit that a test started and does not finish.
release <- function(got) { if (!is.null(got$st$record)) glow_unit_release(got$st$record$handle); invisible(got) }
# A configuration that extends another and changes some settings.
derived <- function(out, parent, analysis, lines = character(0)) {
  f <- file.path(out, paste0(analysis, ".yml"))
  writeLines(c(paste0("extends: ", basename(parent)), paste0("analysis: ", analysis), lines), f)
  f
}
# Finish the data preparation and the training of a study, part by part.
finish_upstream <- function(h, cfg) {
  for (k in c("21", "22")) finish(h, "data_prep", "plink_to_gds", c("--config", cfg, "--chr", k))
  for (k in c("21", "22")) finish(h, "data_prep", "annotate_favor", c("--config", cfg, "--chr", k))
  finish(h, "data_prep", "assemble_pheno", c("--config", cfg))
  for (p in c("b_prepare", "b_estimate", "pi_prepare", "pi_annotate", "pi_train", "pi_evaluate"))
    finish(h, "training", p, c("--config", cfg))
}
# Finish the SNV-set step of a gene run of a study, part by part.
finish_gene_run <- function(h, cfg) {
  finish(h, "snv_set", "prepare", c("--config", cfg))
  for (k in c("21", "22")) finish(h, "snv_set", "scan", c("--config", cfg, "--chr", k))
  for (p in c("aggregate", "plots")) finish(h, "snv_set", p, c("--config", cfg))
}
# The state of a run for a configuration, as the helper builds it (for its internal functions).
run_of <- function(h, config, region_type = NULL) {
  withr::local_dir(ga_dir())
  cfg <- glow_absolute_paths(read_glow_config(config, base = "base-config/base_config.yml",
                                              registry = "base-config/settings-registry.yml", root = "."))
  parts <- glow_read_parts("base-config/parts.yml", "base-config/settings-registry.yml")
  h$.gs_run_new(cfg, glow_layout(cfg), parts, region_type)
}

