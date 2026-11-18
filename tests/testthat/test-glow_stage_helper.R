# ==============================================================================
# File: test-glow_stage_helper.R
# Purpose: Tests of the stage helper that GLOWanalyses ships
#   (packages/GLOWanalyses/lib/glow-stage.R), which every stage script sources to
#   read its configuration, find its output folder and keep the records of its run.
#   The helper builds on read_glow_config(), glow_layout() and the run records of
#   this package. The tests run where GLOWanalyses is next to the package, in the
#   GLOWpipeline monorepo, and are skipped elsewhere, because the helper is not part
#   of this package.

# Whether the lock of a unit is held now (the locks are generations: glow_run_records.R).
unit_locked <- function(folder, stem) GLOWpipeline:::.rr_lock_state(file.path(folder, "records", paste0(stem, ".lock")))$held

test_that("the command line of a stage script is parsed, and an unknown argument stops it", {
  skip_without_helper()
  h <- helper_env()
  pa <- h$.gs_parse_args(c("--config", "a.yml", "--chr", "21", "--combine", "--region-type", "gene", "--force"),
                         c(chr = "value", combine = "flag", "max-genes" = "value"))
  expect_identical(pa$config, "a.yml")
  expect_identical(pa$region_type, "gene")
  expect_true(pa$force)
  expect_identical(pa$opts, list(chr = "21", combine = TRUE, `max-genes` = NULL))
  none <- h$.gs_parse_args(character(0), c(combine = "flag"))
  expect_true(is.na(none$config))
  expect_null(none$region_type)
  expect_false(none$force)
  expect_identical(none$opts, list(combine = FALSE))
  expect_error(h$.gs_parse_args(c("--chunk", "3"), c(chr = "value")), "Unknown argument: --chunk")
  expect_error(h$.gs_parse_args(c("--chr"), c(chr = "value")), "Unknown argument: --chr")
  expect_error(h$.gs_parse_args(c("stray"), character(0)), "Unknown argument: stray")
})

test_that("annotation tests become the objects of GLOWr, and nothing else changes", {
  skip_without_helper()
  h <- helper_env()
  clauses <- list(list(GENCODE.Category = c("exonic", "splicing"), CADD = list(op = "gt", value = 20)),
                  list(GeneHancer = list(op = "nonempty"), MetaSVM = list(op = "in", value = c("D", "007"))))
  out <- h$.gs_as_predicates(clauses)
  expect_identical(out[[1]]$GENCODE.Category, c("exonic", "splicing"))
  expect_identical(out[[1]]$CADD, GLOWr::annotation_predicate("gt", 20))
  expect_identical(out[[2]]$GeneHancer, GLOWr::annotation_predicate("nonempty"))
  expect_identical(out[[2]]$MetaSVM, GLOWr::annotation_predicate("in", c("D", "007")))
  cats <- h$.gs_as_predicates(list(my_mask = list(list(MetaSVM = list(op = "not_in", value = "T")))))
  expect_identical(cats$my_mask[[1]]$MetaSVM, GLOWr::annotation_predicate("not_in", "T"))
  expect_null(h$.gs_as_predicates(NULL))
  expect_identical(h$.gs_as_predicates(c("plof", "missense")), c("plof", "missense"))
  expect_identical(h$.gs_as_predicates(out), out)
  expect_s3_class(GLOWr::variant_filter(annotation_clauses = out), "glow_variant_filter")
})

test_that("a config.R of release 0.3.6 is refused, with a pointer to the migration guide", {
  skip_without_helper()
  h <- helper_env()
  out <- withr::local_tempdir()
  dir.create(file.path(out, "myrun"))
  cfg <- file.path(out, "myrun", "config.R")
  writeLines(c("gds_dir <- 'somewhere'", "stop('a configuration is never sourced')"), cfg)
  for (call in list(list("single_variant", "scan", c("--config", cfg, "--chr", "2")),
                    list("snv_set", "prepare", c("--config", cfg)),
                    list("training", "b_estimate", c("--config", cfg))))
    expect_error(begin(h, call[[1]], call[[2]], call[[3]]),
                 "is a configuration in the form of release 0.3.6 \\(a config.R\\).*see MIGRATION.md", info = call[[2]])
  expect_false(dir.exists(file.path(out, "myrun", "outputs")))            # nothing was created beside it
  expect_error(begin(h, "training", "b_estimate", c("--config", file.path(out, "absent.yml"))), "Configuration not found")
  writeLines("x: 1", file.path(out, "config.txt"))
  expect_error(begin(h, "training", "b_estimate", c("--config", file.path(out, "config.txt"))), "has to be a YAML file \\(.yml\\)")
  # the context of a run names no form any more, and the single-variant files have one stem
  ok <- begin(h, "data_prep", "plink_to_gds", c("--config", study(out), "--chr", "21"))
  expect_null(ok$st$form)
  expect_identical(ok$st$sv_stem, "single_variant")
  release(ok)
})

test_that("a part that is split needs its unit, and the region type is selected", {
  skip_without_helper()
  h <- helper_env()
  out <- withr::local_tempdir()
  cfg <- study(out, list(snv_set = "  region_type: [gene, window]"))
  expect_error(begin(h, "data_prep", "plink_to_gds", c("--config", cfg)), "runs once per chromosome: pass --chr")
  expect_error(begin(h, "single_variant", "scan", c("--config", cfg)), "runs once per chromosome: pass --chr")
  expect_error(begin(h, "snv_set", "scan", c("--config", cfg, "--region-type", "window")), "runs once per chunk: pass --chunk")
  expect_error(begin(h, "snv_set", "prepare", c("--config", cfg)), "lists 2 region types \\(gene, window\\)")
  expect_error(begin(h, "snv_set", "prepare", c("--config", cfg, "--region-type", "coding")),
               "--region-type coding is not among the region types of the configuration \\(gene, window\\)")
  expect_error(begin(h, "training", "b_prepare", c("--config", cfg, "--region-type", "gene")), "--region-type belongs to the SNV-set step")
  # the unit names
  spec <- function(id) GLOWpipeline::glow_read_parts(file.path(ga_dir(), "base-config", "parts.yml"))[[id]]
  expect_identical(h$.gs_unit_name(spec("data_prep.plink_to_gds"), list(chr = "7"), NULL), "chr7")
  expect_identical(h$.gs_unit_name(spec("snv_set.scan"), list(chr = "21"), "coding"), "chr21")
  expect_identical(h$.gs_unit_name(spec("snv_set.scan"), list(chunk = "12"), "window"), "chunk0012")
  expect_identical(h$.gs_unit_name(spec("snv_set.staar_coding"), list(arrayid = "3", spa = "on"), "coding"), "spa.group0003")
  expect_identical(h$.gs_unit_name(spec("snv_set.staar_window"), list(chunk = "3", spa = "off"), "window"), "nospa.chunk0003")
  expect_null(h$.gs_unit_name(spec("snv_set.prepare"), list(), "gene"))
})

test_that("a script that is not the script of the part and its region type is refused before anything is written", {
  skip_without_helper()
  h <- helper_env()
  out <- withr::local_tempdir()
  # (review 2026-11-16, third round) `Rscript 03-snv-set/run-gene.R --config <a coding run>` wrote
  # gene tables into the coding folder and completed the coding unit, and run-coding.R then
  # found its unit complete and did nothing.
  cfg <- study(out, list(snv_set = "  region_type: [coding]"))
  finish_upstream(h, cfg)
  finish(h, "snv_set", "prepare", c("--config", cfg), script_file = "03-snv-set/prepare.R")
  coding <- file.path(out, "runs", "mine", "baseline", "snv_set", "coding")
  before <- sort(list.files(coding, recursive = TRUE, all.files = TRUE))
  ga <- normalizePath(ga_dir())   # absolute: begin() runs in the GLOWanalyses directory
  refused <- "is not the script of the part `snv_set/scan` for the region type `coding`: that is `03-snv-set/run-coding.R`"
  for (wrong in c("03-snv-set/run-gene.R", file.path(ga, "03-snv-set", "run-window.R"), "03-snv-set/run-noncoding.R",
                  "/no/such/folder/run-coding.R"))
    expect_error(begin(h, "snv_set", "scan", c("--config", cfg, "--chr", "21"), script_file = wrong), refused)
  # (review 2026-11-17) The check compared the two base names, so a link or a copy named
  # run-coding.R that held the gene script completed the coding unit. The files are compared.
  spoof <- withr::local_tempdir()
  file.symlink(file.path(ga, "03-snv-set", "run-gene.R"), file.path(spoof, "run-coding.R"))
  expect_error(begin(h, "snv_set", "scan", c("--config", cfg, "--chr", "21"), script_file = file.path(spoof, "run-coding.R")), refused)
  copies <- withr::local_tempdir()
  file.copy(file.path(ga, "03-snv-set", "run-gene.R"), file.path(copies, "run-coding.R"))
  expect_error(begin(h, "snv_set", "scan", c("--config", cfg, "--chr", "21"), script_file = file.path(copies, "run-coding.R")), refused)
  # a process that does not say which script it runs is refused, and not taken as "no check"
  expect_error(begin(h, "snv_set", "scan", c("--config", cfg, "--chr", "21"), script_file = ""),
               "does not say which script it runs.*whose script is `03-snv-set/run-coding.R`")
  expect_identical(sort(list.files(coding, recursive = TRUE, all.files = TRUE)), before)   # no lock, no record, no file
  expect_false(unit_locked(coding, "scan.chr21"))
  # the right script runs: by a relative path, by an absolute path, and through a link of another name
  got <- begin(h, "snv_set", "scan", c("--config", cfg, "--chr", "21"), script_file = file.path(ga, "03-snv-set", "run-coding.R"))
  expect_null(got$st$skipped); release(got)
  file.symlink(file.path(ga, "03-snv-set", "run-coding.R"), file.path(spoof, "my-coding-scan.R"))
  got <- begin(h, "snv_set", "scan", c("--config", cfg, "--chr", "21"), script_file = file.path(spoof, "my-coding-scan.R"))
  expect_null(got$st$skipped); release(got)
  # with a record of the unit in place, the wrong script is still refused and not told "complete"
  finish(h, "snv_set", "scan", c("--config", cfg, "--chr", "21"), script_file = "03-snv-set/run-coding.R")
  expect_error(begin(h, "snv_set", "scan", c("--config", cfg, "--chr", "21"), script_file = "03-snv-set/run-gene.R"),
               "is not the script of the part")
  # a part with one script for every region type, and a part of another step
  expect_error(begin(h, "snv_set", "aggregate", c("--config", cfg), script_file = "04-summary/plots.R"),
               "`04-summary/plots.R` is not the script of the part `snv_set/aggregate` for the region type `coding`: that is `04-summary/aggregate.R`")
  expect_error(begin(h, "training", "b_estimate", c("--config", cfg), script_file = "01-training/b/prepare-b-training-data.R"),
               "is not the script of the part `training/b_estimate`: that is `01-training/b/estimate-b.R`")
  # a script that serves two parts by a flag is the script of both
  sv <- begin(h, "single_variant", "scan", c("--config", cfg, "--chr", "21"), script_file = "02-single-variant/marginal-scan.R")
  expect_identical(sv$st$unit, "chr21"); release(sv)
})

test_that("the running script is the file that source() reads, or the file R was started with", {
  skip_without_helper()
  h <- helper_env()
  # (review 2026-11-17) A stage script read by source() had no identity, so its check was skipped.
  dir <- withr::local_tempdir()
  inner <- file.path(dir, "inner.R"); outer <- file.path(dir, "outer.R")
  writeLines("seen_inner <- find_script()", inner)
  writeLines(c(sprintf("source(%s, local = TRUE)", deparse(inner)), "seen_outer <- find_script()"), outer)
  env <- new.env(); env$find_script <- h$.gs_running_script
  source(outer, local = env)
  expect_identical(env$seen_inner, inner)     # the innermost file that is being read
  expect_identical(env$seen_outer, outer)
  # outside source(): the file of `Rscript <file>`, or "" when the process has none; never NA
  started_with <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  expect_identical(h$.gs_running_script(), if (length(started_with) == 1L) started_with else "")
})

test_that("a unit runs once: it holds a lock, writes a record and the resolved configuration, and is then complete", {
  skip_without_helper()
  h <- helper_env()
  out <- withr::local_tempdir()
  cfg <- study(out)
  gds <- file.path(out, "runs", "mine", "data", "gds")
  got <- begin(h, "data_prep", "plink_to_gds", c("--config", cfg, "--chr", "21"))
  expect_null(got$st$skipped)
  expect_identical(got$st$part, "plink_to_gds"); expect_identical(got$st$unit, "chr21")
  expect_true(unit_locked(gds, "plink_to_gds.chr21"))
  expect_identical(got$st$paths$gds_dir, gds)
  expect_identical(got$st$paths$plink_dir, file.path(out, "plink"))
  # the resolved configuration of the step: the settings, the derived paths and the provenance
  resolved <- file.path(out, "runs", "mine", "baseline", "data_prep", "config-resolved.yml")
  expect_true(file.exists(resolved))
  expect_true(any(grepl("written_by:.*data_prep/plink_to_gds chr21", readLines(resolved))))
  back <- read_glow_config(resolved, base = file.path(ga_dir(), "base-config", "base_config.yml"),
                           registry = file.path(ga_dir(), "base-config", "settings-registry.yml"), root = ga_dir())
  expect_identical(back$data$plink_dir, file.path(out, "plink"))
  expect_identical(back$shared$chroms, c(21L, 22L))
  # the end: the output has to exist, then the record is written and the lock released
  dir.create(gds, showWarnings = FALSE); writeLines("g", file.path(gds, "chr21.gds"))
  withr::with_dir(ga_dir(), h$glow_stage_end(got$st, outputs = file.path(gds, "chr21.gds")))
  expect_false(unit_locked(gds, "plink_to_gds.chr21"))
  rec <- glow_unit_record(gds, "plink_to_gds", "chr21")
  expect_identical(rec$outputs, "chr21.gds")
  expect_true(all(c("data.plink_pattern", "data_prep.gds_pattern", "data_prep.plink_to_gds_opts") %in% names(rec$dependencies$settings)))
  expect_false("shared.chroms" %in% names(rec$dependencies$settings))       # a per-chromosome unit does not depend on the list
  # the folder of PLINK files is identified by the three files of this chromosome
  expect_match(rec$dependencies$files[["data.plink_dir"]], "^set:3 files")
  expect_identical(names(rec$dependencies$versions), c("GFisher", "GLOWanalyses", "GLOWpipeline", "GLOWr"))
  # a second start has nothing to do, and the other chromosome is a unit of its own
  again <- begin(h, "data_prep", "plink_to_gds", c("--config", cfg, "--chr", "21"))
  expect_identical(again$st$skipped, "complete")
  release(expect_null_skipped <- begin(h, "data_prep", "plink_to_gds", c("--config", cfg, "--chr", "22")))
  expect_null(expect_null_skipped$st$skipped)
  # a file of another chromosome that is added or changed leaves the finished unit complete
  writeLines("new", file.path(out, "plink", "chr22.bed")); writeLines("x", file.path(out, "plink", "chr1.bed"))
  expect_identical(begin(h, "data_prep", "plink_to_gds", c("--config", cfg, "--chr", "21"))$st$skipped, "complete")
  # a changed file of its own chromosome does not: the guard asks before the results are replaced
  writeLines("new", file.path(out, "plink", "chr21.bim"))
  expect_error(begin(h, "data_prep", "plink_to_gds", c("--config", cfg, "--chr", "21")), "file data.plink_dir.*a new `dataset` name")
  writeLines("x", file.path(out, "plink", "chr21.bim"))
  expect_identical(begin(h, "data_prep", "plink_to_gds", c("--config", cfg, "--chr", "21"))$st$skipped, "complete")
  # a deleted output makes the unit run again
  unlink(file.path(gds, "chr21.gds"))
  expect_null(begin(h, "data_prep", "plink_to_gds", c("--config", cfg, "--chr", "21"))$st$skipped)
})

test_that("an input that is not finished stops the run and names the part to run first", {
  skip_without_helper()
  h <- helper_env()
  out <- withr::local_tempdir()
  cfg <- study(out)
  expect_error(begin(h, "data_prep", "annotate_favor", c("--config", cfg, "--chr", "21")),
               "reads the results of 1 part\\(s\\) that are not finished.*data_prep, part plink_to_gds \\(`glow run data_prep`\\).*no completion record")
  finish(h, "data_prep", "plink_to_gds", c("--config", cfg, "--chr", "21"))
  # chromosome 21 can go on, and chromosome 22 still waits for its own unit
  expect_null(begin(h, "data_prep", "annotate_favor", c("--config", cfg, "--chr", "21"))$st$skipped)
  expect_error(begin(h, "data_prep", "annotate_favor", c("--config", cfg, "--chr", "22")), "part plink_to_gds")
  # a part of another step is named with its step
  expect_error(begin(h, "single_variant", "scan", c("--config", cfg, "--chr", "21")),
               "data_prep, part annotate_favor.*data_prep, part assemble_pheno|data_prep, part assemble_pheno.*data_prep, part annotate_favor")
  # a path that the user sets takes the place of the part: no record is asked for
  own <- study(out, list(data_prep = "  favor_input_gds_dir: /data/their/gds"), name = "own.yml")
  got <- begin(h, "data_prep", "annotate_favor", c("--config", own, "--chr", "22"))
  expect_null(got$st$skipped)
  expect_identical(got$st$paths$favor_input_dir, "/data/their/gds")
})

test_that("the settings arrive under the variable names of the scripts, and only the inputs of the part are derived", {
  skip_without_helper()
  h <- helper_env()
  out <- withr::local_tempdir()
  cfg <- study(out, list(single_variant = "  use_spa: false",
                         snv_set = c("  region_type: [coding]",
                                     "  filter_spec:", "    annotation_clauses:", "      - {CADD: {op: ge, value: 10}}",
                                     "  categories:", "    my_mask:", "      - {MetaSVM: {op: in, value: [D]}}")))
  root <- file.path(out, "runs", "mine")
  finish_upstream(h, cfg)
  sv <- begin(h, "single_variant", "scan", c("--config", cfg, "--chr", "21"))$env
  expect_identical(sv$use_SPA, FALSE)                       # the old spelling of the switch
  expect_identical(sv$chroms, c(21L, 22L))
  expect_identical(sv$gds_dir, file.path(root, "data", "gds_favor", "flexible", "gds"))
  expect_identical(sv$gds_pattern, "chr{chr}_hg38_favor.gds")   # the annotated files, not the converted ones
  expect_identical(sv$pheno_path, file.path(root, "baseline", "data_prep", "pheno_covar.rds"))
  expect_null(sv$ld_scores_path)                            # the scan neither writes nor reads the LD scores
  expect_true(exists("symlinked_shared_root", envir = sv, inherits = FALSE))
  expect_null(sv$symlinked_shared_root)
  # The LD scores of a chromosome read that chromosome's scan for the reference variants,
  # and only with a threshold above zero. Without one they do not wait for the scan.
  no_threshold <- release(begin(h, "single_variant", "ld_scores", c("--config", cfg, "--chr", "21")))
  expect_null(no_threshold$st$skipped)
  expect_false(any(grepl("single_variant.scan", names(no_threshold$st$record$handle$dependencies$inputs))))
  with_threshold <- study(out, list(single_variant = c("  use_spa: false", "  ld_ref_maf_min: 0.01")), name = "threshold.yml")
  expect_error(begin(h, "single_variant", "ld_scores", c("--config", with_threshold, "--chr", "21")), "single_variant, part scan")
  ss <- begin(h, "snv_set", "prepare", c("--config", cfg))$env
  expect_identical(ss$region_type, "coding")
  expect_identical(ss$filter_spec$annotation_clauses[[1]]$CADD, GLOWr::annotation_predicate("ge", 10))
  expect_identical(ss$categories$my_mask[[1]]$MetaSVM, GLOWr::annotation_predicate("in", "D"))
  expect_s3_class(do.call(GLOWr::variant_filter, ss$filter_spec), "glow_variant_filter")
  expect_identical(ss$pi_model_dir, file.path(root, "baseline", "training", "pi", "pi_models"))
  expect_identical(ss$b_model_path, file.path(root, "baseline", "training", "b", "b_model.rds"))
  expect_null(ss$b_func)
  expect_null(ss$calibration$report)                        # calibration is off: nothing derived
  dp <- begin(h, "data_prep", "plink_to_gds", c("--config", cfg, "--chr", "22", "--force"))$env
  expect_identical(dp$gds_pattern, "chr{chr}_hg38.gds")     # the converted files
  expect_identical(dp$plink_to_gds_opts, list())
})

test_that("an input is derived only when the run needs it", {
  skip_without_helper()
  h <- helper_env()
  out <- withr::local_tempdir()
  cfg <- study(out)
  finish_upstream(h, cfg)
  root <- file.path(out, "runs", "mine")
  # a B function: the SNV-set preparation neither derives nor depends on the B model
  with_func <- study(out, list(snv_set = "  b_func: \"function(maf) rep(1, length(maf))\""), name = "func.yml")
  got <- release(begin(h, "snv_set", "prepare", c("--config", with_func)))
  expect_null(got$env$b_model_path)
  expect_true(is.function(got$env$b_func))
  expect_false(any(grepl("b_estimate", names(got$st$record$handle$dependencies$inputs))))
  expect_true("training.pi_train" %in% names(got$st$record$handle$dependencies$inputs))
  # calibration on without a source: the report of the single-variant step is an input, and it is missing
  on <- study(out, list(snv_set = "  calibration: {enabled: true}"), name = "on.yml")
  expect_error(begin(h, "snv_set", "prepare", c("--config", on)), "single_variant, part calibrate")
  # the two-file route: the scan is an outside file, and the LD scores are an input for the method that reads them
  writeLines("scan", file.path(out, "scan.csv"))
  two_file <- function(method, name) study(out, list(snv_set = paste0("  calibration: {enabled: true, method: ", method,
                                                     ", single_variant_scan: ", file.path(out, "scan.csv"), "}")), name = name)
  expect_error(begin(h, "snv_set", "prepare", c("--config", two_file("ldsc_intercept", "ldsc.yml"))),
               "single_variant, part ld_scores_combine")
  gc <- release(begin(h, "snv_set", "prepare", c("--config", two_file("lambda_gc", "gc.yml"))))
  expect_null(gc$env$calibration$ld_scores)
  expect_match(gc$st$record$handle$dependencies$files[["snv_set.calibration.single_variant_scan"]], "^md5:")
  # calibration switched off reads no source: a report path that is set is neither checked nor recorded
  in_example <- file.path(ga_dir(), "data-example", "cohort", "known-snps.csv")
  off <- study(out, list(snv_set = paste0("  calibration: {enabled: false, report: ", in_example, "}")), name = "off.yml")
  quiet <- release(begin(h, "snv_set", "prepare", c("--config", off)))
  expect_null(quiet$st$skipped)
  expect_false("snv_set.calibration.report" %in% names(quiet$st$record$handle$dependencies$files))
  # switched on, the same path is an input, and it lies in the example data
  on_example <- study(out, list(snv_set = paste0("  calibration: {enabled: true, report: ", in_example, "}")), name = "on-example.yml")
  expect_error(begin(h, "snv_set", "prepare", c("--config", on_example)), "belong to the shipped example.*snv_set.calibration.report")
  # an explicit path wins over the layout and is a file dependency, not a part
  explicit <- begin(h, "snv_set", "prepare", c("--config", study(out, list(
    snv_set = paste0("  pi_model_dir: ", file.path(out, "plink"))), name = "explicit.yml")))
  expect_identical(explicit$env$pi_model_dir, file.path(out, "plink"))
  expect_match(explicit$st$record$handle$dependencies$files[["snv_set.pi_model_dir"]], "^dir:")
  expect_false(any(grepl("pi_train", names(explicit$st$record$handle$dependencies$inputs))))
})

test_that("the control files of the PI training follow pi_control_format, and the record states the choice", {
  skip_without_helper()
  h <- helper_env()
  out <- withr::local_tempdir()
  favor <- file.path(out, "runs", "mine", "data", "gds_favor", "flexible")
  records_in <- dirname(favor)                       # the records of the annotation are kept above the matching method
  config_of <- function(format) study(out, list(training = c(paste0("  pi_control_format: ", format), "  pi_chromosomes: [21, 22]")),
                                      name = paste0("format-", format, ".yml"))
  source_of <- function(format) { run <- run_of(h, config_of(format)); h$.gs_control_source(run, run$layout) }
  expect_identical(source_of("csv"), file.path(favor, "csv"))
  expect_identical(source_of("gds"), file.path(favor, "gds"))
  # auto: the gds folder, until the annotation of every control chromosome has written and recorded a CSV file
  expect_identical(source_of("auto"), file.path(favor, "gds"))
  annotated <- function(k, with_csv) {
    hd <- glow_unit_begin(records_in, "annotate_favor", paste0("chr", k), list(settings = character(0)), force = TRUE)
    files <- file.path(favor, c("gds", if (with_csv) "csv"), sprintf("chr%d_hg38_favor.%s", k, c("gds", if (with_csv) "csv")))
    for (f in files) { dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE); writeLines("a", f) }
    glow_unit_complete(hd, files)
  }
  annotated(21, with_csv = TRUE); annotated(22, with_csv = FALSE)
  expect_identical(source_of("auto"), file.path(favor, "gds"))
  # a CSV file that was put there by hand is not what "auto" reads: no record lists it
  writeLines("by hand", file.path(favor, "csv", "chr22_hg38_favor.csv"))
  expect_identical(source_of("auto"), file.path(favor, "gds"))
  annotated(22, with_csv = TRUE)
  expect_identical(source_of("auto"), file.path(favor, "csv"))
  # the record of a part that reads the controls states the format that "auto" became
  run <- run_of(h, config_of("auto"))
  expect_identical(h$.gs_effective_value(run, "training.pi_control_format"), "csv")
  expect_identical(h$.gs_unit_dependencies(run, "training.pi_train", NULL)$settings[["training.pi_control_format"]], "\"csv\"")
  expect_identical(h$.gs_derive_inputs(run, run$parts[["training.pi_train"]])$training$pi_control_format, "csv")
  # with control files that the user names, the setting is passed on as written
  own <- run_of(h, study(out, list(training = c("  pi_control_format: auto", paste0("  pi_control_source: ", file.path(out, "plink")))),
                         name = "own-controls.yml"))
  expect_identical(h$.gs_effective_value(own, "training.pi_control_format"), "auto")
})

test_that("the PI parts read the annotation of their control chromosomes, not of the chromosomes of the run", {
  skip_without_helper()
  h <- helper_env()
  out <- withr::local_tempdir()
  # the run scans chromosome 22, and the controls come from chromosomes 21 and 22
  cfg <- study(out, list(shared = "  chroms: [22]"))
  finish(h, "data_prep", "plink_to_gds", c("--config", cfg, "--chr", "22"))
  finish(h, "data_prep", "annotate_favor", c("--config", cfg, "--chr", "22"))
  finish(h, "data_prep", "assemble_pheno", c("--config", cfg))
  for (p in c("pi_prepare", "pi_annotate")) finish(h, "training", p, c("--config", cfg))
  expect_error(begin(h, "training", "pi_train", c("--config", cfg)), "data_prep, part annotate_favor.*unit chr21")
  finish(h, "data_prep", "plink_to_gds", c("--config", cfg, "--chr", "21"))
  finish(h, "data_prep", "annotate_favor", c("--config", cfg, "--chr", "21"))
  trained <- finish(h, "training", "pi_train", c("--config", cfg))
  expect_identical(sort(grep("annotate_favor", names(trained$st$record$handle$dependencies$inputs), value = TRUE)),
                   c("data_prep.annotate_favor.chr21", "data_prep.annotate_favor.chr22"))
  expect_identical(begin(h, "training", "pi_train", c("--config", cfg))$st$skipped, "complete")
  # the annotation of a control chromosome is produced again: the trained models are out of date
  finish(h, "data_prep", "annotate_favor", c("--config", cfg, "--chr", "21", "--force"))
  again <- capture.output(got <- release(begin(h, "training", "pi_train", c("--config", cfg))))
  expect_null(got$st$skipped)
  expect_match(paste(again, collapse = " "), "runs again, because a part it read was produced again.*annotate_favor.chr21")
})

test_that("the record states the file-name pattern that the run uses, also when it is derived", {
  skip_without_helper()
  h <- helper_env()
  out <- withr::local_tempdir()
  gds <- file.path(out, "their-gds"); dir.create(gds)
  writeLines("a", file.path(gds, "a22.gds")); writeLines("b", file.path(gds, "b22.gds"))
  writeLines("p", file.path(out, "bundle.rds"))
  with_pattern <- function(pattern, name) study(out, list(
    shared = c("  chroms: [22]", paste0("  analysis_gds_dir: ", gds), paste0("  pheno_path: ", file.path(out, "bundle.rds"))),
    data_prep = c("  pcs: ~", paste0("  favor_gds_pattern: \"", pattern, "\""))), name = name)
  first <- finish(h, "single_variant", "scan", c("--config", with_pattern("a{chr}.gds", "a.yml"), "--chr", "22"))
  deps <- first$st$record$handle$dependencies
  # `shared.analysis_gds_pattern` is null in the configuration, and the record holds what the scan read
  expect_identical(deps$settings[["shared.analysis_gds_pattern"]], "\"a{chr}.gds\"")
  expect_identical(first$env$gds_pattern, NULL)                # (load = FALSE in `finish`)
  expect_match(deps$files[["shared.analysis_gds_dir"]], "^md5:")           # the one file that is read
  # the other pattern reads the other file: the finished scan is not taken for it
  expect_error(begin(h, "single_variant", "scan", c("--config", with_pattern("b{chr}.gds", "b.yml"), "--chr", "22")),
               "setting shared.analysis_gds_pattern: recorded \"a\\{chr\\}.gds\", now \"b\\{chr\\}.gds\"")
})

test_that("the overwrite guard names the changed setting and yields to --force", {
  skip_without_helper()
  h <- helper_env()
  out <- withr::local_tempdir()
  cfg <- study(out)
  finish_upstream(h, cfg)
  finish(h, "snv_set", "prepare", c("--config", cfg))
  changed <- study(out, list(snv_set = "  ld_threshold: 0.8"), name = "changed.yml")
  expect_error(begin(h, "snv_set", "prepare", c("--config", changed)),
               "made with other dependencies.*setting snv_set.ld_threshold: recorded 0.9.*now 0.8.*a new `analysis` name")
  # a setting of another region type does not touch a finished gene run
  other <- study(out, list(snv_set = c("  window_size: 5000", "  categories: [plof]")), name = "other.yml")
  expect_identical(begin(h, "snv_set", "prepare", c("--config", other))$st$skipped, "complete")
  # a part of the dataset level points at a new dataset name
  moved <- study(out, list(data_prep = "  gds_pattern: \"c{chr}.gds\""), name = "moved.yml")
  expect_error(begin(h, "data_prep", "plink_to_gds", c("--config", moved, "--chr", "21")), "a new `dataset` name")
  # --force runs the unit again, and what read its output is then no longer complete
  finish(h, "snv_set", "scan", c("--config", cfg, "--chr", "21"))
  expect_identical(begin(h, "snv_set", "scan", c("--config", cfg, "--chr", "21"))$st$skipped, "complete")
  forced <- finish(h, "snv_set", "prepare", c("--config", changed, "--force"))
  expect_null(forced$st$skipped)
  expect_error(begin(h, "snv_set", "scan", c("--config", changed, "--chr", "21")),
               "made with other dependencies.*input snv_set.prepare.*setting snv_set.ld_threshold|setting snv_set.ld_threshold.*input snv_set.prepare")
  # under the first configuration the preparation is now the one that is not finished
  expect_error(begin(h, "snv_set", "scan", c("--config", cfg, "--chr", "22")), "snv_set, part prepare.*setting snv_set.ld_threshold")
})

test_that("an input that is replaced while a unit runs stops its record", {
  skip_without_helper()
  h <- helper_env()
  out <- withr::local_tempdir()
  cfg <- study(out)
  finish_upstream(h, cfg)
  finish(h, "snv_set", "prepare", c("--config", cfg))
  running <- begin(h, "snv_set", "scan", c("--config", cfg, "--chr", "21"), load = FALSE)
  finish(h, "snv_set", "prepare", c("--config", cfg, "--force"))             # the preparation is replaced meanwhile
  folder <- file.path(out, "runs", "mine", "baseline", "snv_set", "gene")
  writeLines("r", file.path(folder, "chr21.csv"))
  expect_error(withr::with_dir(ga_dir(), h$glow_stage_end(running$st, outputs = file.path(folder, "chr21.csv"))),
               "was replaced while it ran.*input snv_set.prepare")
  expect_null(glow_unit_record(folder, "scan", "chr21"))
  expect_false(unit_locked(folder, "scan.chr21"))
})

test_that("a derived analysis reuses what did not change and is told what to run again", {
  skip_without_helper()
  h <- helper_env()
  out <- withr::local_tempdir()
  base <- study(out, list(shared = "  pi_features: [cadd_phred, linsight]"))
  finish_upstream(h, base)
  finish(h, "snv_set", "prepare", c("--config", base))
  root <- file.path(out, "runs", "mine")
  # A parent that has finished the preparation of its gene run and no scan: the folder is not
  # finished, so nothing of it is reused. The scan of the derived analysis is not "reused"
  # from results that do not exist, and its preparation runs in its own folder.
  same <- derived(out, base, "same_settings")
  expect_error(begin(h, "snv_set", "scan", c("--config", same, "--chr", "21")), "snv_set, part prepare")
  own_prepare <- finish(h, "snv_set", "prepare", c("--config", same))
  expect_null(own_prepare$st$skipped)
  expect_null(release(begin(h, "snv_set", "scan", c("--config", same, "--chr", "21")))$st$skipped)
  # once the parent has finished the whole gene run, a new derived analysis reuses the folder
  for (k in c("21", "22")) finish(h, "snv_set", "scan", c("--config", base, "--chr", k))
  finish(h, "snv_set", "aggregate", c("--config", base))
  later <- derived(out, base, "later")
  # (the plots of the parent are still missing, so the folder is not finished yet)
  expect_error(begin(h, "snv_set", "aggregate", c("--config", later)), "snv_set, part prepare")
  finish(h, "snv_set", "plots", c("--config", base))
  for (part in c("prepare", "aggregate", "plots"))
    expect_identical(begin(h, "snv_set", part, c("--config", later))$st$skipped, "reused")
  expect_identical(begin(h, "snv_set", "scan", c("--config", later, "--chr", "22"))$st$skipped, "reused")
  expect_false(dir.exists(file.path(root, "later", "snv_set")))
  change_pi <- derived(out, base, "change_pi", c("shared:", "  pi_features: [cadd_phred]"))
  # the SNV-set preparation stops and names the training, with the setting that differs
  expect_error(begin(h, "snv_set", "prepare", c("--config", change_pi)),
               "training, part pi_train \\(`glow run training`\\).*in the analysis `baseline`.*setting shared.pi_features")
  # the B folder did not change: its parts are reused, and nothing runs
  reused <- begin(h, "training", "b_estimate", c("--config", change_pi))
  expect_identical(reused$st$skipped, "reused")
  expect_identical(begin(h, "training", "b_prepare", c("--config", change_pi))$st$skipped, "reused")
  expect_false(dir.exists(file.path(root, "change_pi", "training", "b")))
  # the PI folder did change: every part of it runs in the derived analysis, in its order
  expect_error(begin(h, "training", "pi_train", c("--config", change_pi)), "training, part pi_annotate")
  for (p in c("pi_prepare", "pi_annotate", "pi_train", "pi_evaluate"))
    expect_null(finish(h, "training", p, c("--config", change_pi))$st$skipped)
  expect_true(file.exists(file.path(root, "change_pi", "training", "pi", "records", "pi_train.yml")))
  # now the SNV-set preparation runs: new PI models, the B model and the bundle of the parent
  got <- begin(h, "snv_set", "prepare", c("--config", change_pi))
  expect_null(got$st$skipped)
  expect_identical(got$env$pi_model_dir, file.path(root, "change_pi", "training", "pi", "pi_models"))
  expect_identical(got$env$b_model_path, file.path(root, "baseline", "training", "b", "b_model.rds"))
  expect_identical(got$env$pheno_path, file.path(root, "baseline", "data_prep", "pheno_covar.rds"))
  expect_identical(got$st$output_dir, file.path(root, "change_pi", "snv_set", "gene"))
  resolved <- readLines(file.path(root, "change_pi", "snv_set", "gene", "config-resolved.yml"))
  expect_true(any(grepl("training/b_estimate: reused from baseline", resolved, fixed = TRUE)))
  expect_true(any(grepl("data_prep/assemble_pheno: reused from baseline", resolved, fixed = TRUE)))
  expect_false(any(grepl("pi_train: reused", resolved, fixed = TRUE)))
  # (found on 2026-11-16 while the commands of the README were run as written) A later unit of
  # the folder writes the record again. It still states what the folder reused: the notes were
  # lost, because the later unit met the reused parts while it followed its own dependencies.
  release(got)
  finish(h, "snv_set", "prepare", c("--config", change_pi))
  finish(h, "snv_set", "scan", c("--config", change_pi, "--chr", "21"))
  resolved <- readLines(file.path(root, "change_pi", "snv_set", "gene", "config-resolved.yml"))
  expect_true(any(grepl("written_by: \"snv_set/scan chr21\"", resolved)))
  expect_true(any(grepl("training/b_estimate: reused from baseline", resolved, fixed = TRUE)))
  expect_true(any(grepl("data_prep/assemble_pheno: reused from baseline", resolved, fixed = TRUE)))
  expect_false(any(grepl("pi_train: reused", resolved, fixed = TRUE)))
  # --force runs a reusable part in the derived analysis all the same
  expect_null(begin(h, "training", "b_prepare", c("--config", change_pi, "--force"))$st$skipped)
  # an analysis of another dataset is no parent for reuse
  elsewhere <- derived(out, base, "other_data", c("dataset: theirs"))
  expect_error(begin(h, "training", "b_estimate", c("--config", elsewhere)), "training, part b_prepare")
})

test_that("an input inside the example data stops a run of another dataset", {
  skip_without_helper()
  h <- helper_env()
  out <- withr::local_tempdir()
  f <- file.path(out, "forgot.yml")
  writeLines(c("dataset: mine", paste0("output_root: ", file.path(out, "runs"))), f)       # the raw inputs keep their defaults
  expect_error(begin(h, "data_prep", "plink_to_gds", c("--config", f, "--chr", "21")),
               "The dataset is `mine`, and 1 input\\(s\\) of `data_prep/plink_to_gds` belong to the shipped example.*data.plink_dir")
  expect_error(begin(h, "training", "b_prepare", c("--config", f)), "training.b_raw_path")
  # a link into the example data is caught through its real path
  file.symlink(file.path(ga_dir(), "data-example", "cohort", "known-snps.csv"), file.path(out, "my-snps.csv"))
  linked <- file.path(out, "linked.yml")
  writeLines(c(readLines(f), "training:", paste0("  b_raw_path: ", file.path(out, "my-snps.csv"))), linked)
  expect_error(begin(h, "training", "b_prepare", c("--config", linked)), "training.b_raw_path")
  # a file of one's own passes, and so does the example itself
  expect_null(begin(h, "training", "b_prepare", c("--config", study(out)))$st$skipped)
  # a file inside a dataset folder that is marked as example data is example data too
  marked <- file.path(out, "elsewhere", "example"); dir.create(file.path(marked, "baseline"), recursive = TRUE)
  writeLines("marker", file.path(marked, "EXAMPLE-DATA")); writeLines("m", file.path(marked, "baseline", "known.csv"))
  from_marked <- study(out, list(training = paste0("  pi_raw_path: ", file.path(marked, "baseline", "known.csv"))), name = "marked.yml")
  expect_error(begin(h, "training", "pi_prepare", c("--config", from_marked)), "training.pi_raw_path")
  # A folder of one's own whose files are links into the example data: the folder itself is not
  # example data, the files that are read from it are.
  links <- file.path(out, "my-gds"); dir.create(links)
  example_gds <- file.path(ga_dir(), "data-example", "cohort", "example_gds_favor", "flexible", "gds")
  for (k in c(21, 22)) file.symlink(file.path(example_gds, sprintf("chr%d_favor.gds", k)), file.path(links, sprintf("chr%d_favor.gds", k)))
  writeLines("p", file.path(out, "bundle.rds"))
  through_links <- study(out, list(
    shared = c("  chroms: [21, 22]", paste0("  analysis_gds_dir: ", links), "  analysis_gds_pattern: \"chr{chr}_favor.gds\"",
               paste0("  pheno_path: ", file.path(out, "bundle.rds")))), name = "links.yml")
  expect_error(begin(h, "single_variant", "scan", c("--config", through_links, "--chr", "21")),
               "belong to the shipped example.*shared.analysis_gds_dir")
  # the same folder with files of one's own passes
  unlink(file.path(links, c("chr21_favor.gds", "chr22_favor.gds")))
  for (k in c(21, 22)) writeLines("g", file.path(links, sprintf("chr%d_favor.gds", k)))
  expect_null(release(begin(h, "single_variant", "scan", c("--config", through_links, "--chr", "21")))$st$skipped)
})

# The shipped example with its outputs in a folder of the test.
example_in <- function(out) {
  f <- file.path(out, "example-elsewhere.yml")
  writeLines(c("extends: example.yml", paste0("output_root: ", file.path(out, "runs"))), f)
  file.copy(file.path(ga_dir(), "example.yml"), file.path(out, "example.yml"))
  f
}

test_that("the example's pre-built files are supplied once, marked, linked and recorded", {
  skip_without_helper()
  h <- helper_env()
  out <- withr::local_tempdir()
  f <- example_in(out)
  root <- file.path(out, "runs", "example")
  # a look at the configuration writes nothing at all
  looked <- begin(h, "single_variant", "scan", c("--config", f, "--chr", "21"))
  release(looked)
  unlink(file.path(out, "runs"), recursive = TRUE)
  script <- new.env()
  withr::with_dir(ga_dir(), h$glow_stage_begin("single_variant", "scan", options = c(chr = "value", combine = "flag"),
                                                 records = FALSE, envir = script, args = c("--config", f, "--chr", "21"),
                                                 script = NA_character_))
  expect_false(dir.exists(file.path(out, "runs")))
  expect_identical(script$gds_dir, file.path(root, "data", "gds_favor", "flexible", "gds"))   # the settings are there all the same

  first <- capture.output(got <- begin(h, "single_variant", "scan", c("--config", f, "--chr", "21")))
  expect_match(paste(first, collapse = " "), "supplied the pre-built files of 5 unit")
  expect_true(file.exists(file.path(root, "EXAMPLE-DATA")))
  seeded <- c("data/gds_favor/flexible/gds/chr21_favor.gds", "data/gds_favor/flexible/gds/chr22_favor.gds",
              "baseline/data_prep/pheno_covar.rds", "baseline/training/b/b_model.rds",
              "baseline/training/pi/pi_models/model_1.rds", "baseline/training/pi/pi_models/training_metadata.rds")
  expect_true(all(file.exists(file.path(root, seeded))))
  expect_identical(Sys.readlink(file.path(root, "baseline", "data_prep", "gds_favor")), "../../data/gds_favor")
  # each seeded unit has a seed record, and the scan reads the seeded files as finished results
  expect_true(glow_unit_record(file.path(root, "data", "gds_favor"), "annotate_favor", "chr21")$seeded)
  expect_true(glow_unit_record(file.path(root, "baseline", "data_prep"), "assemble_pheno")$seeded)
  b_record <- glow_unit_record(file.path(root, "baseline", "training", "b"), "b_estimate")
  expect_true(b_record$seeded)
  pi_record <- glow_unit_record(file.path(root, "baseline", "training", "pi"), "pi_train")
  expect_true(all(c("pi_models/model_1.rds", "pi_models/training_metadata.rds") %in% pi_record$outputs))
  # the seed record holds the settings of the part and of the parts before it, and the outside files
  # (the outside files by their content: the record of a clone at another path says the same)
  expect_true(all(c("training.b_method", "training.b_qc_filters", "content of training.b_raw_path") %in% names(b_record$dependencies$settings)))
  expect_identical(b_record$dependencies$settings[["content of training.b_raw_path"]],
                   paste0("md5:", unname(tools::md5sum(file.path(ga_dir(), "data-example", "cohort", "known-snps.csv")))))
  expect_identical(b_record$dependencies$settings[["content of training.b_training_rds"]], "none")
  expect_null(got$st$skipped)
  expect_match(got$st$record$handle$dependencies$inputs[["data_prep.annotate_favor.chr21"]], "^seed-")
  expect_identical(got$env$pheno_path, file.path(root, "baseline", "data_prep", "pheno_covar.rds"))
  # a second run copies nothing and keeps the generation marks
  second <- capture.output(release(begin(h, "single_variant", "scan", c("--config", f, "--chr", "22"))))
  expect_false(any(grepl("supplied", second)))
  expect_identical(glow_unit_record(file.path(root, "baseline", "training", "b"), "b_estimate")$generation, b_record$generation)
  release(got)
})

test_that("the pre-built files stand for the settings they were built with, whatever the configuration says later", {
  skip_without_helper()
  h <- helper_env()
  out <- withr::local_tempdir()
  f <- example_in(out)
  gene <- function(name, lines) { g <- file.path(out, name); writeLines(c("extends: example-elsewhere.yml", lines, "snv_set:", "  region_type: [gene]"), g); g }
  # the shipped settings: the SNV-set preparation of the example starts from the supplied files
  expect_null(release(begin(h, "snv_set", "prepare", c("--config", gene("gene.yml", character(0)))))$st$skipped)
  # a derived analysis with other PI features is told to train its own models
  change_pi <- gene("change_pi.yml", c("analysis: change_pi", "shared:", "  pi_features: [cadd_phred]"))
  expect_error(begin(h, "snv_set", "prepare", c("--config", change_pi)), "training, part pi_train.*setting shared.pi_features")
  # a setting of a part before the seeded one counts too: the B model was built from other training data
  earlier <- gene("earlier.yml", c("analysis: earlier", "training:", "  b_filter_autosomes: false"))
  expect_error(begin(h, "snv_set", "prepare", c("--config", earlier)), "training, part b_estimate.*setting training.b_filter_autosomes")
  # and so does an outside file of such a part
  writeLines("other", file.path(out, "other-snps.csv"))
  other_raw <- gene("other_raw.yml", c("analysis: other_raw", "training:", paste0("  b_raw_path: ", file.path(out, "other-snps.csv"))))
  expect_error(begin(h, "snv_set", "prepare", c("--config", other_raw)), "training, part b_estimate.*content of training.b_raw_path")
  # (review 2026-11-16, second round) The same bytes under another path are the same input,
  # as in a clone of the example at another place.
  file.copy(file.path(ga_dir(), "data-example", "cohort", "known-snps.csv"), file.path(out, "same-snps.csv"))
  same_raw <- gene("same_raw.yml", c("analysis: same_raw", "training:", paste0("  b_raw_path: ", file.path(out, "same-snps.csv")),
                                     paste0("  pi_raw_path: ", file.path(out, "same-snps.csv"))))
  expect_null(release(begin(h, "snv_set", "prepare", c("--config", same_raw)))$st$skipped)
  # A file in the folder of the supplied models that their record does not list would be read
  # as a model of the ensemble. The run stops and names it.
  stray <- file.path(out, "runs", "example", "baseline", "training", "pi", "pi_models", "model_999.rds")
  writeLines("not a supplied model", stray)
  expect_error(begin(h, "snv_set", "prepare", c("--config", file.path(out, "gene.yml"))),
               "holds 1 file\\(s\\) that the completion record of `training/pi_train` does not list: model_999.rds")
  unlink(stray)
  expect_null(release(begin(h, "snv_set", "prepare", c("--config", file.path(out, "gene.yml"))))$st$skipped)
  # The copy of example.yml is edited: five PI models become two. The supplied files are the five
  # models still, and their record still says so, so the run is told to train.
  edited <- sub("pi_n_models: 5", "pi_n_models: 2", readLines(file.path(out, "example.yml")))
  expect_false(identical(edited, readLines(file.path(out, "example.yml"))))
  writeLines(edited, file.path(out, "example.yml"))
  expect_error(begin(h, "snv_set", "prepare", c("--config", file.path(out, "gene.yml"))),
               "training, part pi_train.*setting training.pi_n_models: recorded 5L, now 2L")
  pi_dir <- file.path(out, "runs", "example", "baseline", "training", "pi")
  expect_identical(glow_unit_record(pi_dir, "pi_train")$dependencies$settings[["training.pi_n_models"]], "5L")
  # A run trains the two models. From then on the unit is the run's: a deleted model is not
  # restored from the supplied files under the run's record.
  for (p in c("pi_prepare", "pi_annotate")) finish(h, "training", p, c("--config", file.path(out, "gene.yml")))
  trained <- begin(h, "training", "pi_train", c("--config", file.path(out, "gene.yml"), "--force"))
  models <- file.path(pi_dir, "pi_models")
  expect_identical(list.files(models), character(0))           # the supplied models went with their record
  dir.create(models, showWarnings = FALSE); writeLines("m", file.path(models, "model_1.rds")); writeLines("m", file.path(models, "model_2.rds"))
  withr::with_dir(ga_dir(), h$glow_stage_end(trained$st, outputs = models))
  mark <- glow_unit_record(pi_dir, "pi_train")$generation
  unlink(file.path(models, "model_2.rds"))
  expect_error(begin(h, "training", "pi_evaluate", c("--config", file.path(out, "gene.yml"))), "training, part pi_train")
  expect_false(file.exists(file.path(models, "model_2.rds")))
  expect_identical(glow_unit_record(pi_dir, "pi_train")$generation, mark)
  expect_false(glow_unit_state(pi_dir, "pi_train")$complete)
})

test_that("pre-built files that are not the collected set are not used", {
  skip_without_helper()
  h <- helper_env()
  # (review 2026-11-17) A collection that stopped between two files left a new model beside the
  # settings file of the collection before, and a run of the example took the pair as one set.
  root <- withr::local_tempdir(); withr::local_dir(root)
  dir.create("data-example/cohort/models", recursive = TRUE)
  writeLines("b", "data-example/cohort/b.rds"); writeLines("m1", "data-example/cohort/models/m1.rds")
  manifest <- list(source = "data-example/cohort", files = list(list(from = "b.rds"), list(from = "models")))
  files <- c("b.rds", "models/m1.rds")
  expect_identical(h$.gs_seed_files(manifest), files)
  frozen <- list(artifacts = as.list(stats::setNames(unname(tools::md5sum(file.path("data-example/cohort", files))), files)))
  expect_identical(h$.gs_seed_problems(manifest, frozen), character(0))
  # a file of another collection, a file more, a file less, and a settings file without checksums
  writeLines("another model", "data-example/cohort/b.rds")
  expect_match(h$.gs_seed_problems(manifest, frozen), "not the file that was collected: b.rds")
  writeLines("b", "data-example/cohort/b.rds"); writeLines("m2", "data-example/cohort/models/m2.rds")
  expect_match(h$.gs_seed_problems(manifest, frozen), "not of the collected set: models/m2.rds")
  unlink("data-example/cohort/models", recursive = TRUE); dir.create("data-example/cohort/models")
  expect_match(h$.gs_seed_problems(manifest, frozen), "missing: models/m1.rds")
  expect_match(h$.gs_seed_problems(manifest, list()), "states no checksums of the pre-built files")
  # the marker of a collection that began and did not end outweighs everything else
  writeLines("m1", "data-example/cohort/models/m1.rds")
  expect_identical(h$.gs_seed_problems(manifest, frozen), character(0))
  writeLines("x", h$GLOW_SEED_MARKER)
  expect_match(h$.gs_seed_problems(manifest, frozen), "began and did not end")
})

test_that("the shipped seed settings are those that the shipped example gives", {
  skip_without_helper()
  h <- helper_env()
  withr::local_dir(ga_dir())
  built <- h$.gs_seed_settings("example.yml")
  shipped <- yaml::read_yaml("data-example/seed-settings.yml")
  expect_identical(sort(setdiff(names(shipped), "artifacts")), sort(names(built)))
  # (review 2026-11-17) The settings file states the checksum of every pre-built file, and the
  # shipped files are that set: none missing, none more, none replaced.
  manifest <- yaml::read_yaml("data-example/seed-manifest.yml")
  in_place <- h$.gs_seed_files(manifest)
  expect_identical(sort(names(shipped$artifacts)), sort(in_place))
  expect_identical(unname(unlist(shipped$artifacts)[in_place]), unname(tools::md5sum(file.path(manifest$source, in_place))))
  expect_identical(h$.gs_seed_problems(manifest, shipped), character(0))
  expect_false(file.exists(h$GLOW_SEED_MARKER))
  for (id in names(built)) {
    expect_identical(unlist(shipped[[id]]$settings), unlist(built[[id]]$settings), info = id)
    expect_identical(names(shipped[[id]]$inputs), names(built[[id]]$inputs), info = id)
    expect_identical(unlist(shipped[[id]]$inputs), unlist(built[[id]]$inputs), info = id)
    # an outside file is fixed by its content, never by a path or a time: a clone lies elsewhere
    expect_true(all(grepl("^(none|md5:[0-9a-f]{32}|dir:[0-9]+ files, listing md5:[0-9a-f]{32})$", unlist(built[[id]]$inputs))), info = id)
  }
  # (review 2026-11-16, second round) The shipped files are those that the settings describe: the
  # folder of the PI models holds the number of models that the seed settings state, which is
  # also the number that the training recorded for itself.
  shipped_models <- list.files(file.path("data-example", "cohort", "pi_models"), pattern = "^model_[0-9]+\\.rds$")
  stated <- shipped[["training.pi_train"]]$settings[["training.pi_n_models"]]
  expect_identical(paste0(length(shipped_models), "L"), stated)
  expect_identical(as.integer(readRDS(file.path("data-example", "cohort", "pi_models", "training_metadata.rds"))$n_models), length(shipped_models))
  # the closure of a seeded part reaches the parts before it
  parts <- glow_read_parts("base-config/parts.yml")
  closure <- h$.gs_seed_closure(parts, "training.pi_train")
  expect_true(all(c("training.pi_n_models", "training.pi_qc_filters", "training.favor_features", "data_prep.favor_features",
                    "data.plink_pattern", "content of training.pi_raw_path", "content of data.favor_db", "content of data.plink_dir",
                    "content of training.pi_control_source") %in% closure))
  # every setting of a seed is a value, never a path as text: a path would differ between installations
  registry <- attr(read_glow_config("example.yml", base = "base-config/base_config.yml",
                                    registry = "base-config/settings-registry.yml", root = "."), "glow_registry")$entries
  type_of <- vapply(registry, function(e) e$spec$type, character(1))
  names(type_of) <- vapply(registry, function(e) if (e$section == "top") e$name else paste(e$section, e$name, sep = "."), character(1))
  for (id in names(built)) expect_false(any(type_of[intersect(names(built[[id]]$settings), names(type_of))] == "path"), info = id)
})

test_that("a setting is read from the configuration, and passing it on the command line of a YAML run stops", {
  skip_without_helper()
  h <- helper_env()
  out <- withr::local_tempdir()
  cfg <- study(out, list(shared = c("  chroms: [21, 22]", "  n_cores: 3"), snv_set = c("  max_genes: 7", "  allow_incomplete: true")))
  finish_upstream(h, cfg)
  finish(h, "snv_set", "prepare", c("--config", cfg))
  st <- release(begin(h, "snv_set", "scan", c("--config", cfg, "--chr", "21"), load = FALSE))$st
  expect_identical(h$glow_stage_option(st, "max-genes", c("snv_set", "max_genes"), NA_integer_), 7L)
  expect_null(h$glow_stage_option(st, "max-genes", c("snv_set", "max_windows")))     # a null setting: the default
  # the workers are the current setting, so a changed number reaches a new unit without a new preparation
  expect_identical(h$glow_stage_option(st, "n-cores", c("shared", "n_cores"), NULL), 3L)
  # the limit is a dependency of the scan, so a limited run is not taken for a full one
  expect_identical(st$record$handle$dependencies$settings[["snv_set.max_genes"]], "7L")
  # The option of release 0.3.6 with a YAML configuration: what the unit does would differ from
  # what its record states. The script stops before it locks or writes anything, and names the setting.
  folder <- file.path(out, "runs", "mine", "baseline", "snv_set", "gene", "records")
  before <- list.files(folder)
  expect_error(begin(h, "snv_set", "scan", c("--config", cfg, "--chr", "22", "--max-genes", "3"), load = FALSE),
               "--max-genes is not an option of a YAML configuration. Set `snv_set.max_genes` in the configuration")
  expect_error(begin(h, "snv_set", "scan", c("--config", cfg, "--chr", "22", "--n-cores", "4"), load = FALSE),
               "--n-cores is not an option of a YAML configuration. Set `shared.n_cores`")
  expect_identical(list.files(folder), before)
  # the same stop where a script asks for the option's value
  as_given <- st; as_given$opts[["max-genes"]] <- "3"
  expect_error(h$glow_stage_option(as_given, "max-genes", c("snv_set", "max_genes"), NA_integer_),
               "--max-genes is not an option of a YAML configuration. Set `snv_set.max_genes` in the configuration")
  # the selectors of a unit are not settings
  expect_null(release(begin(h, "snv_set", "scan", c("--config", cfg, "--chr", "22"), load = FALSE))$st$skipped)
  # a limited scan and an unlimited one are different units of work: the finished limited scan is not complete for the other
  finish(h, "snv_set", "scan", c("--config", cfg, "--chr", "21"))
  unlimited <- study(out, list(shared = c("  chroms: [21, 22]", "  n_cores: 3"), snv_set = "  allow_incomplete: true"), name = "unlimited.yml")
  finish(h, "snv_set", "prepare", c("--config", unlimited))          # same settings: complete, nothing runs
  expect_error(begin(h, "snv_set", "scan", c("--config", unlimited, "--chr", "21"), load = FALSE),
               "setting snv_set.max_genes: recorded 7L, now NULL")
})

test_that("the units of a part follow from the configuration, and a unit's name gives its arguments", {
  skip_without_helper()
  h <- helper_env()
  out <- withr::local_tempdir()
  cfg <- study(out, list(snv_set = c("  region_type: [window, coding]", "  staar_enabled: true", "  staar_native: true",
                                     "  staar_modes: [spa, nospa]", "  gene_num_in_array: 50")))
  folder <- file.path(out, "runs", "mine", "baseline", "snv_set", "window")
  window <- run_of(h, cfg, "window")
  units <- function(run, id, dir = folder) h$.gs_expected_units(run, run$parts[[id]], dir)
  expect_identical(units(window, "snv_set.prepare"), list(NULL))
  expect_identical(units(window, "data_prep.annotate_favor"), c("chr21", "chr22"))
  # a window part has units once the preparation wrote the table of chunks, whatever records exist
  expect_identical(units(window, "snv_set.scan"), character(0))
  expect_identical(units(window, "snv_set.staar_window"), character(0))
  dir.create(file.path(folder, "shared"), recursive = TRUE)
  saveRDS(data.frame(chunk_id = c(1L, 2L, 7L)), file.path(folder, "shared", "01-genome_chunks.rds"))
  expect_identical(units(window, "snv_set.scan"), c("chunk0001", "chunk0002", "chunk0007"))
  # the native comparison runs every chunk once per mode: one finished unit does not make the part finished
  expect_identical(units(window, "snv_set.staar_window"),
                   c("spa.chunk0001", "spa.chunk0002", "spa.chunk0007", "nospa.chunk0001", "nospa.chunk0002", "nospa.chunk0007"))
  dir.create(file.path(folder, "records")); writeLines("part: staar_window", file.path(folder, "records", "staar_window.nospa.chunk0001.yml"))
  expect_length(units(window, "snv_set.staar_window"), 6L)
  # the arguments that select a unit
  expect_identical(h$.gs_unit_args(NULL), character(0))
  expect_identical(h$.gs_unit_args("chr21"), c("--chr", "21"))
  expect_identical(h$.gs_unit_args("chunk0007"), c("--chunk", "7"))
  expect_identical(h$.gs_unit_args("spa.group0012"), c("--arrayid", "12", "--spa", "on"))
  expect_identical(h$.gs_unit_args("nospa.chunk0003"), c("--chunk", "3", "--spa", "off"))
  expect_error(h$.gs_unit_args("monday"), "Not the name of a unit")
  # the gene groups of the native coding comparison, as STAARpipeline numbers them
  skip_if_not_installed("STAARpipeline")
  coding <- run_of(h, cfg, "coding")
  groups <- h$.gs_gene_groups(coding$cfg)
  genes_info <- get("genes_info", envir = asNamespace("STAARpipeline"))
  per_chr <- ceiling(table(factor(genes_info[, 2], levels = 1:22)) / 50)
  expect_length(groups, per_chr[["21"]] + per_chr[["22"]])
  expect_identical(max(groups), as.integer(sum(per_chr)))                    # chromosome 22 holds the last groups
  expect_identical(min(groups), as.integer(sum(per_chr[1:20]) + 1L))
  native <- units(coding, "snv_set.staar_coding", file.path(out, "runs", "mine", "baseline", "snv_set", "coding"))
  expect_length(native, 2L * length(groups))
  expect_identical(native[1L], sprintf("spa.group%04d", groups[1L]))
  # the coding aggregation reads the native tables, so it waits for them: all of them, in both modes
  deps <- h$.gs_unit_dependencies(coding, "snv_set.aggregate", NULL)
  expect_true("snv_set.staar_coding" %in% names(deps$inputs))
  expect_identical(deps$inputs[["snv_set.staar_coding"]], "missing")
  # a gene run has no native coding tables to wait for, and a coding run without the comparison neither
  gene <- run_of(h, study(out, name = "gene.yml"), "gene")
  expect_false(any(grepl("staar_coding", names(h$.gs_unit_dependencies(gene, "snv_set.aggregate", NULL)$inputs))))
  plain <- run_of(h, study(out, list(snv_set = "  region_type: [coding]"), name = "plain.yml"), "coding")
  expect_false(any(grepl("staar_coding", names(h$.gs_unit_dependencies(plain, "snv_set.aggregate", NULL)$inputs))))
})

test_that("an aggregation may run on the finished units of a scan when the configuration allows it", {
  skip_without_helper()
  h <- helper_env()
  out <- withr::local_tempdir()
  cfg <- study(out)
  finish_upstream(h, cfg)
  finish(h, "snv_set", "prepare", c("--config", cfg))
  finish(h, "snv_set", "scan", c("--config", cfg, "--chr", "21"))
  # by default the aggregation waits for every unit of the scan
  expect_error(begin(h, "snv_set", "aggregate", c("--config", cfg), load = FALSE), "snv_set, part scan.*unit chr22")
  # with allow_incomplete the unit that is missing is recorded as absent
  partial <- study(out, list(snv_set = "  allow_incomplete: true"), name = "partial.yml")
  got <- finish(h, "snv_set", "aggregate", c("--config", partial))
  inputs <- got$st$record$handle$dependencies$inputs
  expect_identical(inputs[["snv_set.scan.chr22"]], "absent")
  expect_false(inputs[["snv_set.scan.chr21"]] %in% c("absent", "missing"))
  # when the missing unit finishes, the aggregation is out of date and runs again
  finish(h, "snv_set", "scan", c("--config", partial, "--chr", "22"))
  expect_null(release(begin(h, "snv_set", "aggregate", c("--config", partial), load = FALSE))$st$skipped)
  # a setting that only another region type reads does not touch the finished aggregation
  finish(h, "snv_set", "aggregate", c("--config", partial))
  other <- study(out, list(snv_set = c("  allow_incomplete: true", "  merge_gap: 9", "  cmac_cutoff: 3", "  gene_num_in_array: 10")), name = "other.yml")
  expect_identical(begin(h, "snv_set", "aggregate", c("--config", other), load = FALSE)$st$skipped, "complete")
})

test_that("a part records the files it wrote: in a shared folder, in a folder of its own, and for one unit of a scan", {
  skip_without_helper()
  h <- helper_env()
  out <- withr::local_tempdir()
  # two parts write into one folder, one after the other
  shared <- file.path(out, "aggregated"); dir.create(shared)
  before <- h$glow_stage_watch(shared)
  writeLines("t", file.path(shared, "table.csv"))
  expect_identical(h$glow_stage_written(shared, before), file.path(shared, "table.csv"))
  before <- h$glow_stage_watch(shared)
  writeLines("p", file.path(shared, "plot.pdf"))
  expect_identical(h$glow_stage_written(shared, before), file.path(shared, "plot.pdf"))     # not the table of the other part
  before <- h$glow_stage_watch(shared)
  Sys.sleep(0.05); writeLines("t2", file.path(shared, "table.csv"))
  expect_identical(h$glow_stage_written(shared, before), file.path(shared, "table.csv"))     # a file that is written again
  expect_identical(h$glow_stage_written(shared, h$glow_stage_watch(shared)), character(0))
  # (review 2026-11-16, second round) a file whose name starts with a dot is a file the part wrote
  before <- h$glow_stage_watch(shared)
  writeLines("l", file.path(shared, ".family_lambdas_all.rds"))
  expect_identical(h$glow_stage_written(shared, before), file.path(shared, ".family_lambdas_all.rds"))
  # a folder of the part's own starts empty, and only inside the part's output folder
  yaml <- list(output_dir = out)
  models <- file.path(out, "pi_models"); dir.create(file.path(models, "sub"), recursive = TRUE)
  writeLines("old", file.path(models, "model_9.rds")); writeLines("old", file.path(models, "sub", "x.rds"))
  h$glow_stage_clear_dir(yaml, models)
  expect_identical(list.files(models, recursive = TRUE), character(0))
  expect_error(h$glow_stage_clear_dir(list(output_dir = models), out), "is not inside the output folder")
  # the files of one unit of a scan: its tables and its side files, and not those of another unit
  results <- file.path(out, "results"); dir.create(results)
  for (f in c("glow_chr2.csv", "glow_chr2_errors.csv", "glow_chr21.csv", "glow_chr21_summary.csv")) writeLines("x", file.path(results, f))
  for (k in c("evidence", "result_objects")) { dir.create(file.path(out, k)); writeLines("x", file.path(out, k, sprintf("%s_chr2.rds", k))) }
  writeLines("x", file.path(out, "evidence", "evidence_chr21.rds"))
  expect_identical(sort(basename(h$glow_scan_unit_files(yaml, results, "glow_chr2", "chr2"))),
                   c("evidence_chr2.rds", "glow_chr2.csv", "glow_chr2_errors.csv", "result_objects_chr2.rds"))
  h$glow_scan_unit_files(yaml, results, "glow_chr2", "chr2", clear = TRUE)
  expect_identical(h$glow_scan_unit_files(yaml, results, "glow_chr2", "chr2"), character(0))
  expect_identical(sort(basename(h$glow_scan_unit_files(yaml, results, "glow_chr21", "chr21"))),
                   c("evidence_chr21.rds", "glow_chr21.csv", "glow_chr21_summary.csv"))
})

test_that("an aggregation reads the files of the finished units of this run, and no other file", {
  skip_without_helper(); skip_on_os("windows")
  # (review 2026-11-16, second round) The aggregation took every table it found in the results
  # folder: the table of a unit that had not finished, and of a chromosome that an earlier
  # configuration of the folder had scanned.
  h <- helper_env()
  out <- withr::local_tempdir()
  cfg <- study(out, list(snv_set = "  allow_incomplete: true"))
  finish_upstream(h, cfg)
  finish(h, "snv_set", "prepare", c("--config", cfg))
  finish(h, "snv_set", "scan", c("--config", cfg, "--chr", "21"))
  folder <- file.path(out, "runs", "mine", "baseline", "snv_set", "gene")
  writeLines("of a unit that stopped", file.path(folder, "scan.chr22.out.txt"))
  writeLines("of a chromosome that is not in the run any more", file.path(folder, "scan.chr20.out.txt"))
  got <- begin(h, "snv_set", "aggregate", c("--config", cfg), load = FALSE)
  view <- withr::with_dir(ga_dir(), h$glow_stage_input_dir(got$st))
  expect_identical(view, file.path(folder, ".inputs-aggregate"))
  expect_identical(sort(list.files(view, recursive = TRUE)), c("prepare.out.txt", "scan.chr21.out.txt"))
  expect_identical(Sys.readlink(file.path(view, "scan.chr21.out.txt")), file.path(folder, "scan.chr21.out.txt"))
  expect_identical(readLines(file.path(view, "scan.chr21.out.txt")), "result")
  # the view is gone at the end of the unit, and it held links only: the results stay
  agg <- file.path(folder, "aggregate.out.txt"); writeLines("a", agg)
  withr::with_dir(ga_dir(), h$glow_stage_end(got$st, outputs = agg))
  expect_false(dir.exists(view))
  expect_true(all(file.exists(file.path(folder, c("prepare.out.txt", "scan.chr21.out.txt", "scan.chr22.out.txt")))))
  # once the second unit has finished, its table is in the view of the next aggregation
  unlink(file.path(folder, "scan.chr22.out.txt"))
  finish(h, "snv_set", "scan", c("--config", cfg, "--chr", "22"))
  again <- begin(h, "snv_set", "aggregate", c("--config", cfg), load = FALSE)
  view <- withr::with_dir(ga_dir(), h$glow_stage_input_dir(again$st))
  expect_identical(sort(list.files(view, recursive = TRUE)), c("prepare.out.txt", "scan.chr21.out.txt", "scan.chr22.out.txt"))
  release(again)
  # a run without records reads the run folder itself
  expect_identical(h$glow_stage_input_dir(list(output_dir = "/some/run", record = NULL)), "/some/run")
})

test_that("a folder is emptied without following a link, and only inside the output folder of its part", {
  skip_without_helper(); skip_on_os("windows")
  # (review 2026-11-16, second round) The folder of the PI models was emptied through any link in
  # it, and the subfolder could be named `records`.
  h <- helper_env()
  out <- withr::local_tempdir()
  st <- list(output_dir = file.path(out, "pi"))
  models <- file.path(st$output_dir, "pi_models"); dir.create(file.path(models, "sub"), recursive = TRUE)
  writeLines("m", file.path(models, "model_1.rds")); writeLines("s", file.path(models, "sub", "old.rds")); writeLines("h", file.path(models, ".hidden"))
  elsewhere <- file.path(out, "elsewhere"); dir.create(elsewhere); writeLines("keep", file.path(elsewhere, "precious.txt"))
  file.symlink(elsewhere, file.path(models, "linked-folder"))
  file.symlink(file.path(elsewhere, "precious.txt"), file.path(models, "linked-file"))
  h$glow_stage_clear_dir(st, models)
  expect_identical(list.files(models, recursive = TRUE, all.files = TRUE, no.. = TRUE), character(0))
  expect_identical(readLines(file.path(elsewhere, "precious.txt")), "keep")       # what the links pointed at is untouched
  # a folder that is a link to another place, a folder outside, and a folder of the workflow itself
  file.symlink(elsewhere, file.path(st$output_dir, "alias"))
  expect_error(h$glow_stage_clear_dir(st, file.path(st$output_dir, "alias")), "is not inside the output folder")
  expect_error(h$glow_stage_clear_dir(st, elsewhere), "is not inside the output folder")
  dir.create(file.path(st$output_dir, "records")); writeLines("r", file.path(st$output_dir, "records", "pi_train.yml"))
  expect_error(h$glow_stage_clear_dir(st, file.path(st$output_dir, "records")), "a folder of the workflow itself")
  expect_true(file.exists(file.path(st$output_dir, "records", "pi_train.yml")))
  expect_identical(readLines(file.path(elsewhere, "precious.txt")), "keep")
})

test_that("names that would make two units write one file, or a part write into the records, stop the run", {
  skip_without_helper()
  # (review 2026-11-16, second round) `gds_pattern: all.gds` made both chromosomes write one file,
  # and `pi_models_subdir: records` had the training empty its own records folder.
  h <- helper_env()
  out <- withr::local_tempdir()
  one_file <- study(out, list(data_prep = "  gds_pattern: all.gds"), name = "one.yml")
  expect_error(begin(h, "data_prep", "plink_to_gds", c("--config", one_file, "--chr", "21")),
               "`data_prep.gds_pattern` is `all.gds`.*has to contain \\{chr\\}")
  expect_false(dir.exists(file.path(out, "runs")))                               # before any folder or lock
  for (key in c("favor_gds_pattern", "favor_csv_pattern"))
    expect_error(begin(h, "data_prep", "annotate_favor",
                       c("--config", study(out, list(data_prep = paste0("  ", key, ": annotated.out")), name = paste0(key, ".yml")), "--chr", "21")),
                 paste0("`data_prep.", key, "` is `annotated.out`"))
  for (name in c("records", "jobs", "logs"))
    expect_error(begin(h, "training", "pi_train",
                       c("--config", study(out, list(training = paste0("  pi_models_subdir: ", name)), name = paste0(name, ".yml")))),
                 paste0("`training.pi_models_subdir` is `", name, "`, a folder that the workflow keeps for itself"))
  # the shipped names pass
  expect_null(h$.gs_check_names(run_of(h, study(out, name = "plain.yml"))$cfg))
})

test_that("an annotated GDS folder of the user stands in for the annotation for every part that reads it", {
  skip_without_helper()
  # (review 2026-11-16, second round) With shared.analysis_gds_dir set, the annotation was left
  # out of the run while the phenotype assembly still waited for its record and read its folder.
  h <- helper_env()
  out <- withr::local_tempdir()
  mine <- file.path(out, "my-agds"); dir.create(mine)
  for (k in c("21", "22")) writeLines("g", file.path(mine, paste0("study_chr", k, ".gds")))
  supplied <- c("  chroms: [21, 22]", paste0("  analysis_gds_dir: ", mine), "  analysis_gds_pattern: study_chr{chr}.gds")
  cfg <- study(out, list(shared = supplied))
  got <- begin(h, "data_prep", "assemble_pheno", c("--config", cfg))
  expect_null(got$st$skipped)                                                    # no record of the annotation is asked for
  expect_identical(got$st$paths$read_favor_gds_dir, mine)                        # and the script reads the user's folder and pattern
  expect_identical(got$st$paths$read_favor_gds_pattern, "study_chr{chr}.gds")
  deps <- got$st$record$handle$dependencies
  expect_false(any(grepl("annotate_favor", names(deps$inputs))))
  # the one file it reads, the chromosome of sample_gds_chr, identifies the input
  expect_identical(deps$files[["shared.analysis_gds_dir"]], glow_file_fingerprint(file.path(mine, "study_chr22.gds")))
  release(got)
  # the PCs from the annotated files read the same folder, and from the plain files the folder of the conversion
  pcs <- study(out, list(shared = supplied, data_prep = "  pcs: {id_col: sample.id}"), name = "pcs.yml")
  from_favor <- begin(h, "data_prep", "compute_pcs", c("--config", pcs))
  expect_identical(from_favor$st$paths$read_favor_gds_dir, mine)
  expect_null(from_favor$st$paths$read_gds_dir)
  release(from_favor)
  plain <- study(out, list(shared = supplied, data_prep = c("  pcs: {id_col: sample.id}", "  pc_source: gds")), name = "plain.yml")
  expect_error(begin(h, "data_prep", "compute_pcs", c("--config", plain)), "data_prep, part plink_to_gds")
  for (k in c("21", "22")) finish(h, "data_prep", "plink_to_gds", c("--config", plain, "--chr", k))
  from_plain <- begin(h, "data_prep", "compute_pcs", c("--config", plain))
  expect_identical(from_plain$st$paths$read_gds_dir, file.path(out, "runs", "mine", "data", "gds"))
  release(from_plain)
  # a folder of plain GDS files of the user stands in for the conversion in the same way
  own_gds <- file.path(out, "my-gds"); dir.create(own_gds)
  for (k in c("21", "22")) writeLines("g", file.path(own_gds, paste0("chr", k, "_hg38.gds")))
  brought <- study(out, list(shared = supplied, data_prep = c("  pcs: {id_col: sample.id}", "  pc_source: gds",
                                                             paste0("  favor_input_gds_dir: ", own_gds))), name = "brought.yml")
  from_own <- begin(h, "data_prep", "compute_pcs", c("--config", brought))
  expect_identical(from_own$st$paths$read_gds_dir, own_gds)
  expect_false(any(grepl("plink_to_gds", names(from_own$st$record$handle$dependencies$inputs))))
  expect_true("data_prep.favor_input_gds_dir" %in% names(from_own$st$record$handle$dependencies$files))
  release(from_own)
  # without a folder of the user the scripts read what the run wrote
  own <- study(out, name = "own.yml")
  for (k in c("21", "22")) finish(h, "data_prep", "annotate_favor", c("--config", own, "--chr", k))
  default <- begin(h, "data_prep", "assemble_pheno", c("--config", own))
  expect_identical(default$st$paths$read_favor_gds_dir, file.path(out, "runs", "mine", "data", "gds_favor", "flexible", "gds"))
  expect_identical(default$st$paths$read_favor_gds_pattern, default$env$favor_gds_pattern)
  release(default)
})

test_that("with run.records false a unit gets its settings and the resolved configuration, and keeps no record", {
  skip_without_helper()
  # (child task C5) A task of a workflow engine runs in a container of its own: there is no
  # record to compare with, and the engine keeps track of what is finished.
  h <- helper_env()
  out <- withr::local_tempdir()
  cfg <- study(out, list(run = "  records: false",
                         shared = c("  chroms: [21, 22]", paste0("  analysis_gds_dir: ", file.path(out, "plink")),
                                    "  analysis_gds_pattern: \"chr{chr}.bed\"", paste0("  pheno_path: ", file.path(out, "pheno.csv")))))
  step_dir <- file.path(out, "runs", "mine", "baseline", "single_variant")
  got <- begin(h, "single_variant", "scan", c("--config", cfg, "--chr", "21"))
  expect_null(got$st$skipped)
  expect_null(got$st$record)                                   # no lock, no record to write
  expect_identical(got$st$unit, "chr21")
  expect_identical(got$env$gds_dir, file.path(out, "plink"))
  expect_identical(got$env$pheno_path, file.path(out, "pheno.csv"))
  expect_true(file.exists(file.path(step_dir, "config-resolved.yml")))
  expect_true(any(grepl("^\\s*records: false", readLines(file.path(step_dir, "config-resolved.yml")))))
  expect_false(dir.exists(file.path(step_dir, "records")))
  expect_null(withr::with_dir(ga_dir(), h$glow_stage_end(got$st, outputs = "anything")))
  expect_false(dir.exists(file.path(step_dir, "records")))
  # nothing is skipped at a second start, and a part is not stopped for an input without a record
  expect_null(begin(h, "single_variant", "scan", c("--config", cfg, "--chr", "21"))$st$skipped)
  expect_null(begin(h, "single_variant", "calibrate", c("--config", cfg))$st$skipped)
  # a derived input that the configuration leaves null is placed in the task's own layout
  ld <- begin(h, "single_variant", "ld_scores", c("--config", cfg, "--chr", "22"))
  expect_identical(ld$env$ld_scores_path, file.path(step_dir, "ld_scores", "ld_scores.rds"))
  # a setting on the command line is refused all the same
  expect_error(begin(h, "snv_set", "scan", c("--config", cfg, "--chr", "21", "--n-cores", "2"), load = FALSE),
               "--n-cores is not an option of a YAML configuration")
})

test_that("the files of a scan unit are cleared in their own folder only, also under a folder name with a wildcard", {
  skip_without_helper()
  h <- helper_env()
  # (review 2026-11-16, third round) glow_scan_unit_files(clear = TRUE) under the output root
  # `case*` removed the unit's table of the root `case-other` as well.
  root <- withr::local_tempdir()
  st_of <- function(dir) list(output_dir = dir)
  for (d in c("case-other", "case*")) {
    dir.create(file.path(root, d, "results"), recursive = TRUE); dir.create(file.path(root, d, "evidence"))
    writeLines("x", file.path(root, d, "results", "glow_chr21.csv")); writeLines("x", file.path(root, d, "results", "glow_chr21_errors.csv"))
    writeLines("x", file.path(root, d, "evidence", "evidence_chr21.rds"))
  }
  starred <- file.path(root, "case*")
  listed <- h$glow_scan_unit_files(st_of(starred), file.path(starred, "results"), "glow_chr21", "chr21")
  expect_identical(length(listed), 3L)
  h$glow_scan_unit_files(st_of(starred), file.path(starred, "results"), "glow_chr21", "chr21", clear = TRUE)
  expect_identical(length(list.files(starred, recursive = TRUE)), 0L)
  expect_identical(length(list.files(file.path(root, "case-other"), recursive = TRUE)), 3L)
})

test_that("supplied annotated GDS files need a named folder of PI controls", {
  skip_without_helper()
  h <- helper_env()
  out <- withr::local_tempdir()
  agds <- file.path(out, "agds"); dir.create(agds)
  for (k in c(21, 22)) writeLines("x", file.path(agds, sprintf("chr%d_a.gds", k)))
  # (review 2026-11-16, third round) With the annotated GDS files supplied, the annotation is
  # left out of the run, and the PI training waited for it: "run data_prep annotate_favor
  # first", which that run could not do.
  supplied <- study(out, list(shared = c("  chroms: [21, 22]", paste0("  analysis_gds_dir: ", agds), "  analysis_gds_pattern: \"chr{chr}_a.gds\"")),
                    name = "supplied.yml")
  for (part in c("pi_train", "pi_evaluate"))
    expect_error(begin(h, "training", part, c("--config", supplied)),
                 "supplies the annotated GDS files \\(`shared.analysis_gds_dir`\\).*`training.pi_control_source` has to name their folder")
  expect_false(dir.exists(file.path(out, "runs")))                       # nothing was written
  # the parts that do not read the controls are not concerned
  got <- begin(h, "training", "b_prepare", c("--config", supplied)); expect_null(got$st$skipped); release(got)
  # with the folder named, the training goes on to what it really waits for
  named <- study(out, list(shared = c("  chroms: [21, 22]", paste0("  analysis_gds_dir: ", agds), "  analysis_gds_pattern: \"chr{chr}_a.gds\""),
                           training = c(paste0("  pi_control_source: ", agds), "  pi_control_format: gds")), name = "named.yml")
  expect_error(begin(h, "training", "pi_train", c("--config", named)), "training, part pi_annotate")
  # without supplied GDS files nothing changes: the training waits for the annotation of the dataset
  expect_error(begin(h, "training", "pi_train", c("--config", study(out, name = "plain.yml"))), "pi_annotate|annotate_favor")
})

