# ==============================================================================
# File: test-glow_run_command.R
# Purpose: Tests of the dispatcher behind the `glow` command that GLOWanalyses
#   ships (packages/GLOWanalyses/lib/glow-run.R): the command line, which parts
#   belong to a run, the units of a part, the order and the stop at a failure, and
#   what a dry run shows. The tests run where GLOWanalyses is next to the package
#   and are skipped elsewhere (helper-glow-stage.R).

skip_without_command <- function() {
  skip_without_helper()
  testthat::skip_if_not(file.exists(file.path(ga_dir(), "lib", "glow-run.R")), "the glow command is not next to the package")
}
# The dispatcher's functions, sourced from the GLOWanalyses directory without running it.
command_env <- function() {
  env <- new.env(parent = globalenv())
  withr::with_dir(ga_dir(), sys.source("lib/glow-run.R", envir = env))
  env
}
read_study <- function(file) glow_absolute_paths(read_glow_config(
  file, base = file.path(ga_dir(), "base-config", "base_config.yml"),
  registry = file.path(ga_dir(), "base-config", "settings-registry.yml"), root = ga_dir()))
parts_table <- function() glow_read_parts(file.path(ga_dir(), "base-config", "parts.yml"),
                                          file.path(ga_dir(), "base-config", "settings-registry.yml"))
# Run the dispatcher as the launcher does, and keep what it prints.
glow <- function(g, ...) {
  out <- character(0)
  withr::with_dir(ga_dir(), out <- capture.output(g$glow_run_main(c(...))))
  out
}

test_that("the command line names a step and carries no setting", {
  skip_without_command()
  g <- command_env()
  cmd <- g$glow_run_parse(c("run", "snv-set", "--config", "study.yml", "--part", "scan", "--force", "--dry-run"))
  expect_identical(cmd$step, "snv_set")                      # the hyphenated spelling is accepted
  expect_identical(cmd$config, "study.yml")
  expect_identical(cmd$part, "scan")
  expect_true(cmd$force); expect_true(cmd$dry_run); expect_false(cmd$wait)
  plain <- g$glow_run_parse(c("run", "data_prep"))
  expect_true(is.na(plain$config)); expect_null(plain$part); expect_false(plain$force); expect_false(plain$dry_run)
  expect_true(g$glow_run_parse(c("run", "training", "--wait"))$wait)
  expect_error(g$glow_run_parse(c("run", "summary")), "Unknown step `summary`.*Usage: glow run")
  expect_error(g$glow_run_parse(c("go", "training")), "The command is `glow run <step>`")
  expect_error(g$glow_run_parse(c("run")), "The command is `glow run <step>`")
  expect_error(g$glow_run_parse(c("run", "snv_set", "--chr", "21")), "Unknown argument `--chr`.*belongs in the configuration")
  expect_error(g$glow_run_parse(c("run", "snv_set", "--config")), "Unknown argument `--config`")
  expect_error(g$glow_run_parse(c("run", "snv_set", "--n-cores", "4")), "belongs in the configuration")
})

test_that("the control folder of the PI models is asked for only when the run trains PI models", {
  skip_without_command()
  g <- command_env()
  out <- withr::local_tempdir()
  agds <- file.path(out, "agds"); dir.create(agds); models <- file.path(out, "models"); dir.create(models)
  for (k in c(21, 22)) writeLines("x", file.path(agds, sprintf("chr%d_a.gds", k)))
  shared <- c("  chroms: [21, 22]", paste0("  analysis_gds_dir: ", agds), "  analysis_gds_pattern: \"chr{chr}_a.gds\"")
  # supplied annotated GDS files and PI models to train: the folder of the controls has to be named
  supplied <- study(out, list(shared = shared), name = "supplied.yml")
  expect_error(glow(g, "run", "training", "--config", supplied, "--dry-run"),
               "`training.pi_control_source` has to name their folder")
  # (review 2026-11-17) The check ran for every part of the step, also for the PI parts that a
  # study with supplied PI models leaves out. Such a study could not run its B parts.
  both <- study(out, list(shared = shared, snv_set = paste0("  pi_model_dir: ", models)), name = "both.yml")
  shown <- glow(g, "run", "training", "--config", both, "--dry-run")
  expect_true(any(grepl("^  training/b_prepare ", shown)))
  expect_identical(length(grep("training/pi_(prepare|annotate|train|evaluate) +not part of this run", shown)), 4L)
  expect_false(dir.exists(file.path(out, "runs", "mine", "baseline")))            # a dry run wrote no result
})

test_that("a part belongs to a run by its conditions", {
  skip_without_command()
  g <- command_env()
  out <- withr::local_tempdir()
  parts <- parts_table()
  wanted <- function(id, cfg, rt = NULL) g$glow_run_wanted(parts[[id]], cfg, rt)
  plain <- read_study(study(out))
  # a part without a condition, and conditions that hold
  expect_true(wanted("single_variant.scan", plain)$wanted)
  expect_true(wanted("data_prep.plink_to_gds", plain)$wanted)
  expect_true(wanted("training.b_estimate", plain)$wanted)
  expect_true(wanted("single_variant.ld_scores", plain)$wanted)
  # the study has no PCs: the PC part is left out, and the reason is given
  no_pcs <- wanted("data_prep.compute_pcs", plain)
  expect_false(no_pcs$wanted); expect_match(no_pcs$why, "`data_prep.pcs` is null")
  with_pcs <- read_study(study(out, list(data_prep = "  pcs: {id_col: sample.id}"), name = "pcs.yml"))
  expect_true(wanted("data_prep.compute_pcs", with_pcs)$wanted)
  own_pcs <- read_study(study(out, list(data_prep = "  pcs: {path: /data/pcs.rds}"), name = "ownpcs.yml"))
  expect_match(wanted("data_prep.compute_pcs", own_pcs)$why, "`data_prep.pcs.path` is /data/pcs.rds")
  # a path that the user sets leaves out the part that would produce the file
  brought <- read_study(study(out, list(shared = c("  analysis_gds_dir: /data/agds", "  pheno_path: /data/bundle.rds"),
                                        snv_set = c("  pi_model_dir: /data/pi", "  b_func: \"function(maf) 1\"")), name = "brought.yml"))
  for (id in c("data_prep.plink_to_gds", "data_prep.annotate_favor", "data_prep.assemble_pheno", "training.b_prepare",
               "training.b_estimate", "training.pi_prepare", "training.pi_annotate", "training.pi_train", "training.pi_evaluate"))
    expect_false(wanted(id, brought)$wanted, info = id)
  # a value: the LD scores are for the method that reads them
  gc <- read_study(study(out, list(single_variant = "  calibration_method: lambda_gc"), name = "gc.yml"))
  expect_false(wanted("single_variant.ld_scores", gc)$wanted)
  expect_false(wanted("single_variant.ld_scores_combine", gc)$wanted)
  expect_true(wanted("single_variant.calibrate", gc)$wanted)
  # the native comparison: its region type, the comparator and the native switch
  expect_match(wanted("snv_set.staar_coding", plain, "gene")$why, "a part of the coding run only")
  expect_match(wanted("snv_set.staar_coding", plain, "coding")$why, "`snv_set.staar_enabled` is FALSE")
  staar <- read_study(study(out, list(snv_set = c("  staar_enabled: true", "  staar_native: true")), name = "staar.yml"))
  expect_true(wanted("snv_set.staar_coding", staar, "coding")$wanted)
  expect_true(wanted("snv_set.staar_window", staar, "window")$wanted)
  expect_false(wanted("snv_set.staar_window", staar, "coding")$wanted)
  no_modes <- read_study(study(out, list(snv_set = c("  staar_enabled: true", "  staar_modes: []")), name = "nomodes.yml"))
  expect_match(wanted("snv_set.staar_coding", no_modes, "coding")$why, "`snv_set.staar_modes` is null")
  # (review 2026-11-16, second round) A part that its conditions leave out stays in the run while a
  # part of the run reads its output from the layout. With the annotated GDS files supplied, the
  # annotation is left out, and so is the conversion from PLINK. When the PCs are computed from
  # the plain GDS files, the conversion is needed again, and the annotation still is not.
  needed <- function(id, cfg) g$glow_run_wanted(parts[[id]], cfg, NULL, glow_layout(cfg), parts)$wanted
  agds <- c("  analysis_gds_dir: /data/agds")
  supplied <- read_study(study(out, list(shared = agds), name = "agds.yml"))
  expect_false(needed("data_prep.plink_to_gds", supplied)); expect_false(needed("data_prep.annotate_favor", supplied))
  expect_true(needed("data_prep.assemble_pheno", supplied))
  favor_pcs <- read_study(study(out, list(shared = agds, data_prep = "  pcs: {id_col: sample.id}"), name = "agds-pcs.yml"))
  expect_false(needed("data_prep.plink_to_gds", favor_pcs)); expect_true(needed("data_prep.compute_pcs", favor_pcs))
  plain_pcs <- read_study(study(out, list(shared = agds, data_prep = c("  pcs: {id_col: sample.id}", "  pc_source: gds")), name = "agds-plain.yml"))
  expect_true(needed("data_prep.plink_to_gds", plain_pcs)); expect_false(needed("data_prep.annotate_favor", plain_pcs))
  own_gds <- read_study(study(out, list(shared = agds, data_prep = c("  pcs: {id_col: sample.id}", "  pc_source: gds", "  favor_input_gds_dir: /data/gds")),
                              name = "agds-own.yml"))
  expect_false(needed("data_prep.plink_to_gds", own_gds))
  for (id in c("data_prep.plink_to_gds", "data_prep.annotate_favor", "training.b_estimate", "training.pi_train"))
    expect_false(needed(id, brought), info = id)                                 # the closure adds nothing where every reader has a path
})

test_that("the units of a part follow the chromosomes, the chunks and the modes", {
  skip_without_command()
  g <- command_env()
  out <- withr::local_tempdir()
  parts <- parts_table()
  cfg <- read_study(study(out, list(snv_set = c("  region_type: [gene, window]", "  staar_modes: [spa, nospa]"))))
  lay <- glow_layout(cfg)
  units <- function(id, rt = NULL) g$glow_run_units(parts[[id]], cfg, lay, rt)
  expect_identical(units("training.b_estimate"), list(list(unit = NULL, args = character(0))))
  expect_identical(units("data_prep.plink_to_gds"), list(list(unit = "chr21", args = c("--chr", "21")),
                                                         list(unit = "chr22", args = c("--chr", "22"))))
  expect_identical(vapply(units("snv_set.scan", "gene"), `[[`, "", "unit"), c("chr21", "chr22"))
  # a window scan has units only when its preparation has written the chunks
  expect_length(units("snv_set.scan", "window"), 0L)
  shared <- file.path(lay$steps$snv_set[["window"]], "shared"); dir.create(shared, recursive = TRUE)
  saveRDS(data.frame(chunk_id = c(1L, 2L, 3L)), file.path(shared, "01-genome_chunks.rds"))
  win <- units("snv_set.scan", "window")
  expect_identical(vapply(win, `[[`, "", "unit"), c("chunk0001", "chunk0002", "chunk0003"))
  expect_identical(win[[3]]$args, c("--chunk", "3"))
  native <- units("snv_set.staar_window", "window")
  expect_identical(vapply(native, `[[`, "", "unit"), c("spa.chunk0001", "spa.chunk0002", "spa.chunk0003",
                                                       "nospa.chunk0001", "nospa.chunk0002", "nospa.chunk0003"))
  expect_identical(native[[1]]$args, c("--chunk", "1", "--spa", "on"))
  expect_identical(native[[4]]$args, c("--chunk", "1", "--spa", "off"))
  # the unit names are the ones the stage helper gives a script with these arguments
  h <- helper_env()
  for (u in native) {
    opts <- as.list(stats::setNames(u$args[c(FALSE, TRUE)], sub("^--", "", u$args[c(TRUE, FALSE)])))
    expect_identical(h$.gs_unit_name(parts[["snv_set.staar_window"]], opts, "window"), u$unit)
  }
  for (u in win) expect_identical(h$.gs_unit_name(parts[["snv_set.scan"]], list(chunk = u$args[2]), "window"), u$unit)
})

test_that("the units of a part run in order, and the first failure stops the step", {
  skip_without_command()
  g <- command_env()
  out <- withr::local_tempdir()
  script <- file.path(out, "unit.R")
  writeLines(c("a <- commandArgs(trailingOnly = TRUE)", "cat(a, '\\n', file = file.path(dirname(a[1]), 'ran.txt'), append = TRUE)",
               "if ('--chr' %in% a && a[which(a == '--chr') + 1] == '2') quit(status = 3)"), script)
  spec <- list(step = "demo", part = "scan", script = script, call = "--combine")
  units <- lapply(1:3, function(k) list(unit = paste0("chr", k), args = c("--chr", as.character(k))))
  expect_error(capture.output(g$glow_run_local(spec, units, common_args = file.path(out, "cfg.yml"), dry_run = FALSE)),
               "demo/scan chr2 failed \\(exit status 3\\).*run the same command again")
  ran <- readLines(file.path(out, "ran.txt"))
  expect_length(ran, 2L)                                   # chr1 ran, chr2 failed, chr3 was not started
  expect_match(ran[1], "cfg.yml --combine --chr 1")        # the common arguments, the part's own, the unit's
  unlink(file.path(out, "ran.txt"))
  capture.output(n <- g$glow_run_local(spec, units[c(1, 3)], common_args = file.path(out, "cfg.yml"), dry_run = FALSE))
  expect_identical(n, 2L)
  expect_length(readLines(file.path(out, "ran.txt")), 2L)
})

test_that("a dry run shows every unit with its state and starts nothing", {
  skip_without_command()
  g <- command_env(); h <- helper_env()
  out <- withr::local_tempdir()
  cfg <- study(out, list(snv_set = "  region_type: [gene, window]"))
  shown <- glow(g, "run", "data_prep", "--config", cfg, "--dry-run")
  expect_match(shown[1], "glow run data_prep: dataset mine, analysis baseline.*dry run")
  expect_true(any(grepl("data_prep/plink_to_gds chr21 +ready .*Rscript 00-data-prep/plink-to-gds.R --config .* --chr 21", shown)))
  expect_true(any(grepl("data_prep/annotate_favor chr22 +waiting +\\(for plink_to_gds\\)", shown)))
  expect_true(any(grepl("data_prep/compute_pcs +not part of this run: `data_prep.pcs` is null", shown)))
  expect_false(dir.exists(file.path(out, "runs")))          # nothing was created
  # after a part has finished, the dry run says so and shows what is ready next
  for (k in c("21", "22")) finish(h, "data_prep", "plink_to_gds", c("--config", cfg, "--chr", k))
  again <- glow(g, "run", "data_prep", "--config", cfg, "--dry-run")
  expect_true(any(grepl("data_prep/plink_to_gds chr21 +complete", again)))
  expect_true(any(grepl("data_prep/annotate_favor chr22 +ready", again)))
  # one part only, the region types in turn, and a window scan before its preparation
  one <- glow(g, "run", "snv_set", "--config", cfg, "--part", "scan", "--dry-run")
  expect_true(any(grepl("--- region type: gene ---", one, fixed = TRUE)))
  expect_true(any(grepl("snv_set/scan chr21 +waiting .*run-gene.R --config .* --region-type gene --chr 21", one)))
  expect_true(any(grepl("snv_set/scan +its units are known after the preparation has run", one)))
  expect_false(any(grepl("snv_set/prepare", one)))
  # --force is passed on to every unit
  forced <- glow(g, "run", "data_prep", "--config", cfg, "--part", "plink_to_gds", "--force", "--dry-run")
  expect_true(any(grepl("plink_to_gds chr21 +ready +\\(forced\\).* --force --chr 21", forced)))
  # what the command refuses
  expect_error(glow(g, "run", "training", "--config", cfg, "--part", "scan"), "has no part `scan`. Its parts, in order: b_prepare, b_estimate")
  expect_error(glow(g, "run", "training", "--config", file.path(out, "absent.yml")), "Configuration not found")
  r <- file.path(out, "config.R"); writeLines("x <- 1", r)
  expect_error(glow(g, "run", "training", "--config", r), "form of release 0.3.6 \\(a config.R\\).*see MIGRATION.md")
  writeLines("x: 1", file.path(out, "config.txt"))
  expect_error(glow(g, "run", "training", "--config", file.path(out, "config.txt")), "takes a YAML configuration")
  help <- glow(g, "--help")
  expect_identical(help[1], "Usage: glow run <step> [--config <file>] [--part <name>] [--force] [--dry-run] [--wait]")
  # the help says what a dry run writes, and what --wait is for
  expect_true(any(grepl("pre-built files are", help)) && any(grepl("job scripts are written", help)))
  expect_true(any(grepl("--wait +with run.mode slurm", help)))
})

test_that("glow run stops for a configuration that switches the records off", {
  skip_without_helper()
  # (child task C5) run.records false is for the tasks of a workflow engine. The command decides
  # what to run from the records, so it does not start without them.
  out <- withr::local_tempdir()
  cfg <- study(out, list(run = "  records: false"))
  env <- new.env(parent = globalenv())
  withr::with_dir(ga_dir(), sys.source("lib/glow-run.R", envir = env))
  expect_error(withr::with_dir(ga_dir(), env$glow_run_main(c("run", "single_variant", "--config", cfg))),
               "`run.records` is false in this configuration.*`glow run` needs the completion records")
  expect_false(dir.exists(file.path(out, "runs")))
})
