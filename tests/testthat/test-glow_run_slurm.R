# ==============================================================================
# File: test-glow_run_slurm.R
# Purpose: Tests of the SLURM mode of the `glow` command that GLOWanalyses ships
#   (packages/GLOWanalyses/lib/glow-slurm.R): the resources of a job, the job
#   script, the plan of a call as a dry run shows it (which units would be
#   submitted, after which jobs, and which have to wait), and calls that submit to a
#   scheduler of the test's own. Nothing reaches a real scheduler. The
#   tests run where GLOWanalyses is next to the package and are skipped elsewhere
#   (helper-glow-stage.R).

skip_without_slurm_mode <- function() {
  skip_without_helper()
  testthat::skip_if_not(file.exists(file.path(ga_dir(), "lib", "glow-slurm.R")), "the SLURM mode is not next to the package")
}
slurm_env <- function() {
  env <- new.env(parent = globalenv())
  withr::with_dir(ga_dir(), { sys.source("lib/glow-run.R", envir = env); sys.source("lib/glow-slurm.R", envir = env) })
  env
}
slurm_study <- function(out, extra = list(), name = "slurm.yml", run = character(0)) {
  extra$run <- c("  mode: slurm", run)
  read <- study(out, extra, name = name)
  list(file = read, cfg = glow_absolute_paths(read_glow_config(
    read, base = file.path(ga_dir(), "base-config", "base_config.yml"),
    registry = file.path(ga_dir(), "base-config", "settings-registry.yml"), root = ga_dir())))
}
slurm_parts <- function() glow_read_parts(file.path(ga_dir(), "base-config", "parts.yml"),
                                          file.path(ga_dir(), "base-config", "settings-registry.yml"))

test_that("the resources of a job come from the table, the configuration and the workers", {
  skip_without_slurm_mode()
  g <- slurm_env(); parts <- slurm_parts()
  out <- withr::local_tempdir()
  plain <- slurm_study(out)$cfg
  expect_identical(g$glow_slurm_resources(parts[["training.b_estimate"]], plain), list(mem = "16G", time = "06:00:00", cpus = 1L))
  # the values that the job templates of release 0.3.6 used for the heavy parts
  expect_identical(g$glow_slurm_resources(parts[["data_prep.annotate_favor"]], plain), list(mem = "40G", time = "04:00:00", cpus = 1L))
  expect_identical(g$glow_slurm_resources(parts[["snv_set.staar_coding"]], plain)$time, "08:00:00")
  # the configuration overrides them: `default` for every part, then the part's own
  over <- slurm_study(out, name = "over.yml", run = c("  slurm:", "    resources:",
    "      default: {time: \"12:00:00\", cpus: 2}", "      annotate_favor: {mem: 64G}"))$cfg
  expect_identical(g$glow_slurm_resources(parts[["data_prep.annotate_favor"]], over), list(mem = "64G", time = "12:00:00", cpus = 2L))
  expect_identical(g$glow_slurm_resources(parts[["training.b_estimate"]], over), list(mem = "16G", time = "12:00:00", cpus = 2L))
  # a part that starts workers asks for as many cores as it starts workers, whatever `cpus` says
  cores <- slurm_study(out, list(shared = "  n_cores: 8"), name = "cores.yml",
                       run = c("  slurm:", "    resources:", "      default: {cpus: 2}"))$cfg
  expect_identical(g$glow_slurm_resources(parts[["snv_set.scan"]], cores)$cpus, 8L)
  expect_identical(g$glow_slurm_resources(parts[["single_variant.scan"]], cores)$cpus, 8L)
  expect_identical(g$glow_slurm_resources(parts[["snv_set.aggregate"]], cores)$cpus, 2L)
  # a value that sbatch would not read as one word stops the call
  blank <- slurm_study(out, name = "blank.yml", run = c("  slurm:", "    resources:", "      default: {mem: \"16 G\"}"))$cfg
  expect_error(g$glow_slurm_resources(parts[["training.b_estimate"]], blank), "are not ones sbatch reads: mem `16 G`")
})

# A scheduler of the test's own, in front of the real one: `sbatch` gives the next job id and,
# unless `finish_at_once`, puts the job in a queue file, and `squeue` prints that file. `state`
# holds the queue (`queue`: one line `<id> <name>` per task), the arguments of every sbatch call
# (`args-<id>`) and the number of squeue calls (`looks`). `squeue_status` and `sbatch_status`
# make the program fail, `free_after` empties the queue of the lines of other jobs at that
# squeue call, and `empty_after` empties it of every line. Like sbatch, the program takes the job name from its command line first, then
# from SBATCH_JOB_NAME, then from the script, and `env-<id>` holds the SBATCH_* variables it
# was started with.
local_scheduler <- function(finish_at_once = FALSE, squeue_status = 0L, sbatch_status = 0L, sbatch_says = NULL,
                            free_after = NA, empty_after = NA, env = parent.frame()) {
  bin <- withr::local_tempdir(.local_envir = env)
  state <- file.path(bin, "state"); dir.create(state)
  writeLines(character(0), file.path(state, "queue")); writeLines("7000", file.path(state, "next")); writeLines("0", file.path(state, "looks"))
  writeLines(c("#!/bin/bash", sprintf("S=%s", shQuote(state)),
               "n=$(( $(cat $S/looks) + 1 )); echo $n > $S/looks",
               if (!is.na(free_after)) sprintf("if [ $n -ge %d ]; then grep -v ' other' $S/queue > $S/queue.new; mv $S/queue.new $S/queue; fi", free_after),
               if (!is.na(empty_after)) sprintf("if [ $n -ge %d ]; then : > $S/queue; fi", empty_after),
               if (squeue_status != 0L) c("echo 'slurm_load_jobs error: Unable to contact slurm controller' >&2", sprintf("exit %d", squeue_status)),
               "cat $S/queue"), file.path(bin, "squeue"))
  writeLines(c("#!/bin/bash", sprintf("S=%s", shQuote(state)),
               "script=\"${@: -1}\"",
               "name=''; for a in \"$@\"; do case \"$a\" in --job-name=*) name=\"${a#--job-name=}\";; esac; done",
               "if [ -z \"$name\" ]; then name=\"${SBATCH_JOB_NAME:-}\"; fi",
               "if [ -z \"$name\" ]; then name=$(grep -m1 '^#SBATCH --job-name=' \"$script\" | sed 's/^#SBATCH --job-name=//'); fi",
               "id=$(( $(cat $S/next) + 1 )); echo $id > $S/next",
               "printf '%s\\n' \"$@\" > $S/args-$id",
               "env | grep '^SBATCH_' > $S/env-$id || true",
               if (!finish_at_once) "echo \"$id $name\" >> $S/queue",
               if (is.null(sbatch_says)) "echo $id" else sprintf("echo %s", shQuote(sbatch_says)),
               sprintf("exit %d", sbatch_status)), file.path(bin, "sbatch"))
  Sys.chmod(file.path(bin, c("squeue", "sbatch")), "755")
  withr::local_envvar(PATH = paste(bin, Sys.getenv("PATH"), sep = ":"), .local_envir = env)
  withr::local_options(glow.slurm.reserve_wait = 0, .local_envir = env)
  state
}
submissions <- function(state) length(list.files(state, pattern = "^args-"))
# Run the command and keep what it prints.
slurm_run <- function(g, ...) { o <- character(0); withr::with_dir(ga_dir(), o <- capture.output(g$glow_run_main(c("run", ...)))); o }
job_name <- function(g, st, parts, id, unit = NULL, rt = NULL)
  withr::with_dir(ga_dir(), g$glow_slurm_job_name(st$cfg, glow_layout(st$cfg), parts[[id]], rt, unit))


test_that("a job script checks the configuration, sets up R, runs one stage script and starts no other job", {
  skip_without_slurm_mode(); skip_on_os("windows")
  g <- slurm_env(); parts <- slurm_parts()
  out <- withr::local_tempdir()
  st <- slurm_study(out, run = c("  slurm:", "    setup: [\"module load R/4.3\"]", "    conda_env: my_env"))
  cfg <- st$cfg
  jobs <- file.path(out, "jobs"); logs <- file.path(out, "logs", "slurm")
  chain <- withr::with_dir(ga_dir(), g$glow_slurm_chain(cfg, attr(parts, "glow_md5")))
  name <- "glow-mine-data_prep-plink_to_gds-chr21-0a1b2c3d"
  path <- withr::with_dir(ga_dir(), g$glow_slurm_job_script(name,
    c("00-data-prep/plink-to-gds.R", "--config", file.path(out, "my study.yml"), "--chr", "21"),
    parts[["data_prep.plink_to_gds"]], cfg, jobs, logs, chain))
  expect_identical(path, file.path(jobs, paste0(name, ".sh")))
  expect_true(dir.exists(logs))                                    # the log folder exists before submission
  s <- readLines(path)
  expect_identical(s[1], "#!/usr/bin/env bash")
  expect_true(paste0("#SBATCH --job-name=", name) %in% s)
  expect_true(paste0("#SBATCH --output=", file.path(logs, paste0(name, "-%j.log"))) %in% s)
  expect_true(all(c("#SBATCH --cpus-per-task=1", "#SBATCH --mem=16G", "#SBATCH --time=04:00:00") %in% s))
  # the user's setup lines come before the R environment, and the activation is not under nounset
  expect_lt(which(s == "module load R/4.3"), which(grepl("^conda activate 'my_env'", s)))
  expect_lt(which(s == "set +u"), which(grepl("^conda activate", s)))
  expect_gt(which(s == "set -u"), which(grepl("^conda activate", s)))
  expect_true(paste0("cd '", ga_dir(), "'") %in% s)
  # the one command, with an argument that has a space kept as one argument
  expect_true(paste0("\"$RSCRIPT\" '00-data-prep/plink-to-gds.R' '--config' '", file.path(out, "my study.yml"), "' '--chr' '21'") %in% s)
  expect_true("RSCRIPT=Rscript" %in% s)
  expect_false(any(grepl("sbatch|srun", s[!startsWith(s, "#")])))
  # the checksums of the configuration: the files of the chain, the registry and the table of parts,
  # checked before anything else runs
  sums <- file.path(jobs, paste0(name, ".config.md5"))
  expect_identical(readLines(sums), chain)
  expect_true(all(vapply(c(st$file, file.path(ga_dir(), "base-config", c("base_config.yml", "settings-registry.yml", "parts.yml"))),
                         function(f) any(endsWith(chain, paste0("  ", normalizePath(f)))), logical(1))))
  expect_lt(which(grepl("^md5sum --quiet -c ", s)), which(s == "module load R/4.3"))
  expect_identical(system2("md5sum", c("--quiet", "-c", shQuote(sums)), stdout = FALSE, stderr = FALSE), 0L)
  # (review 2026-11-16, second round) The configuration is edited after the submission: the job
  # asked for its resources under the earlier settings, so it stops before it sets up R.
  writeLines(c(readLines(st$file), "single_variant:", "  mac_cutoff: 3"), st$file)
  err <- file.path(out, "job.err")
  expect_identical(suppressWarnings(system2("bash", shQuote(path), stdout = FALSE, stderr = err)), 1L)
  expect_match(paste(readLines(err), collapse = " "), "a file of the configuration changed after this job was submitted")
  expect_false(any(grepl("conda", readLines(err))))                 # it did not get as far as the R setup
  # a log folder with a blank: the path is not written as a directive, and sbatch gets it as an argument
  spaced <- withr::with_dir(ga_dir(), g$glow_slurm_job_script("j", c("x.R"), parts[["snv_set.scan"]], cfg, jobs, file.path(out, "my logs")))
  expect_false(any(grepl("^#SBATCH --output", readLines(spaced))))
  # a part that starts workers states the cores it asks for and the workers it starts
  cores <- slurm_study(out, list(shared = "  n_cores: 4"), name = "cores.yml")$cfg
  s3 <- readLines(withr::with_dir(ga_dir(), g$glow_slurm_job_script("k", c("x.R"), parts[["snv_set.scan"]], cores, jobs, logs)))
  expect_true("#SBATCH --cpus-per-task=4" %in% s3)
  expect_true(any(grepl("asks for 4 core\\(s\\), and the unit starts shared.n_cores = 4 worker\\(s\\)", s3)))
})

test_that("a job calls the Rscript that the configuration names, and no text of a setting is run by the shell", {
  skip_without_slurm_mode(); skip_on_os("windows")
  g <- slurm_env(); parts <- slurm_parts()
  out <- withr::local_tempdir()
  # (review 2026-11-16, second round) With `custom-r` and `Rscript` in one folder the job ran
  # `Rscript`, whatever run.slurm.rscript said.
  bin <- file.path(out, "bin"); dir.create(bin)
  for (prog in c("custom-r", "Rscript")) {
    writeLines(c("#!/bin/sh", sprintf("echo %s \"$@\" > %s", prog, shQuote(file.path(out, "ran.txt")))), file.path(bin, prog))
    Sys.chmod(file.path(bin, prog), "755")
  }
  own <- slurm_study(out, name = "own.yml", run = c("  slurm:", paste0("    rscript: ", file.path(bin, "custom-r"))))$cfg
  path <- withr::with_dir(ga_dir(), g$glow_slurm_job_script("own-r", c("x.R", "--chr", "21"), parts[["snv_set.scan"]], own,
                                                             file.path(out, "jobs"), file.path(out, "logs")))
  s <- readLines(path)
  expect_false(any(grepl("conda", s[!startsWith(s, "#")])))
  expect_true(paste0("RSCRIPT='", file.path(bin, "custom-r"), "'") %in% s)
  expect_true(paste0("export PATH='", bin, "':\"$PATH\"") %in% s)
  expect_identical(system2("bash", shQuote(path), stdout = FALSE, stderr = FALSE), 0L)
  expect_identical(readLines(file.path(out, "ran.txt")), "custom-r x.R --chr 21")
  # A path that holds shell syntax is a name and nothing else, also inside the message of a job
  # whose program does not exist.
  sentinel <- file.path(out, "sentinel")
  odd <- slurm_study(out, name = "odd.yml", run = c("  slurm:", sprintf("    rscript: \"%s/$(touch %s)/Rscript\"", out, sentinel)))$cfg
  path <- withr::with_dir(ga_dir(), g$glow_slurm_job_script("odd-r", c("x.R"), parts[["snv_set.scan"]], odd,
                                                             file.path(out, "jobs"), file.path(out, "logs")))
  err <- file.path(out, "odd.err")
  expect_identical(suppressWarnings(system2("bash", shQuote(path), stdout = FALSE, stderr = err)), 1L)
  expect_match(paste(readLines(err), collapse = " "), "run.slurm.rscript is not an executable file: .*\\$\\(touch ", perl = TRUE)
  expect_false(file.exists(sentinel))
  env_odd <- slurm_study(out, name = "envodd.yml", run = c("  slurm:", sprintf("    conda_env: \"x$(touch %s)\"", sentinel)))$cfg
  s <- readLines(withr::with_dir(ga_dir(), g$glow_slurm_job_script("odd-env", c("x.R"), parts[["snv_set.scan"]], env_odd,
                                                                    file.path(out, "jobs"), file.path(out, "logs"))))
  activate <- grep("^conda activate", s, value = TRUE)
  expect_false(grepl("\\$\\(touch", gsub("'[^']*'", "", activate)))            # the text appears inside single quotes only
})

test_that("run.slurm.options is split by quotes, and an option that the command sets itself is refused", {
  skip_without_slurm_mode()
  g <- slurm_env()
  expect_identical(g$glow_slurm_options(""), character(0))
  expect_identical(g$glow_slurm_options(NULL), character(0))
  expect_identical(g$glow_slurm_options("--partition=short  --account=lab"), c("--partition=short", "--account=lab"))
  # (review 2026-11-16, second round) --comment="my study" was passed as two arguments.
  expect_identical(g$glow_slurm_options("--partition=short --comment=\"my study\" -A 'my lab' --constraint=a\\ b"),
                   c("--partition=short", "--comment=my study", "-A", "my lab", "--constraint=a b"))
  expect_identical(g$glow_slurm_options("--qos=long -p short --exclusive --nodes=1 --mem-bind=local --time-min=10 --export=ALL --comment=''"),
                   c("--qos=long", "-p", "short", "--exclusive", "--nodes=1", "--mem-bind=local", "--time-min=10", "--export=ALL", "--comment="))
  expect_error(g$glow_slurm_options("--comment=\"my study"), "a quote that is not closed")
  # what the command sets itself: the long form, a shortened long form, the short form, with the value attached or apart
  for (bad in c("--job-name=renamed", "--job-name renamed", "--job-n=renamed", "-J renamed", "-Jrenamed", "--cpus-per-task=4", "-c 4", "-c4",
                "--output=o.log", "-o o.log", "--error=e.log", "-e e.log", "--chdir=/tmp", "-D /tmp", "--array=1-3", "-a 1-3",
                "--dependency=afterok:1", "-d afterok:1", "--kill-on-invalid-dep=no", "--parsable", "--wait", "-W", "--wrap=hostname",
                "--mem=4G", "--mem-per-cpu=2G", "--time=10", "-t 10", "--test-only", "-Hc4", "--partition=short --job-name=x"))
    expect_error(g$glow_slurm_options(bad), "which the command sets itself", info = bad)
  expect_error(g$glow_slurm_options("--mem=4G"), "`run.slurm.options` carries `--mem=4G` \\(--mem\\).*run.slurm.resources")
})

test_that("a job name belongs to one unit: names that run together and two output roots do not share it", {
  skip_without_slurm_mode()
  g <- slurm_env(); parts <- slurm_parts()
  out <- withr::local_tempdir()
  named <- function(dataset, analysis, root = file.path(out, "runs"), file) {
    f <- file.path(out, file)
    writeLines(c(sub("^dataset: .*", paste0("dataset: ", dataset), sub("^output_root: .*", paste0("output_root: ", root), readLines(study(out)))),
                 paste0("analysis: ", analysis), "run:", "  mode: slurm"), f)
    list(file = f, cfg = read_study_slurm(f))
  }
  read_study_slurm <- function(f) glow_absolute_paths(read_glow_config(f, base = file.path(ga_dir(), "base-config", "base_config.yml"),
    registry = file.path(ga_dir(), "base-config", "settings-registry.yml"), root = ga_dir()))
  one <- named("a-b", "c", file = "one.yml"); two <- named("a", "b-c", file = "two.yml")
  n1 <- job_name(g, one, parts, "single_variant.calibrate"); n2 <- job_name(g, two, parts, "single_variant.calibrate")
  expect_match(n1, "^glow-a-b-c-single_variant-calibrate-[0-9a-f]{8}$")
  # (review 2026-11-16, second round) These two studies had one job name, so a job of the one hid a unit of the other.
  expect_identical(sub("-[0-9a-f]{8}$", "", n1), sub("-[0-9a-f]{8}$", "", n2))
  expect_false(identical(n1, n2))
  elsewhere <- named("a-b", "c", root = file.path(out, "other-root"), file = "three.yml")
  expect_false(identical(n1, job_name(g, elsewhere, parts, "single_variant.calibrate")))
  expect_identical(n1, job_name(g, one, parts, "single_variant.calibrate"))            # the same unit, the same name
  expect_false(identical(job_name(g, one, parts, "single_variant.scan", "chr21"), job_name(g, one, parts, "single_variant.scan", "chr22")))
  # a unit of the dataset level is one unit for every analysis of the dataset, and has one job name
  sibling <- named("a-b", "d", file = "four.yml")
  expect_identical(job_name(g, one, parts, "data_prep.plink_to_gds", "chr21"), job_name(g, sibling, parts, "data_prep.plink_to_gds", "chr21"))
  expect_match(job_name(g, one, parts, "data_prep.plink_to_gds", "chr21"), "^glow-a-b-data_prep-plink_to_gds-chr21-[0-9a-f]{8}$")
  expect_false(identical(job_name(g, one, parts, "single_variant.calibrate"), job_name(g, sibling, parts, "single_variant.calibrate")))
  # a real path: the same folder through a link has the same name
  skip_on_os("windows")
  dir.create(file.path(out, "runs")); file.symlink(file.path(out, "runs"), file.path(out, "runs-link"))
  linked <- named("a-b", "c", root = file.path(out, "runs-link"), file = "five.yml")
  expect_identical(job_name(g, linked, parts, "single_variant.calibrate"), n1)
})

test_that("the queue and a submission are read only from a scheduler that answers", {
  skip_without_slurm_mode(); skip_on_os("windows")
  g <- slurm_env()
  local({
    state <- local_scheduler()
    writeLines(c("101 glow-a", "102_1 other-array", "102_2 other-array", "102_3 other-array", "103 glow-a"), file.path(state, "queue"))
    q <- g$glow_slurm_queue()
    expect_identical(q$tasks, 5L)                                              # an array task counts as one task each
    expect_identical(q$ids[["glow-a"]], c("101", "103"))                       # two jobs of one name are both known
    expect_identical(length(q$ids[["other-array"]]), 3L)
    writeLines(character(0), file.path(state, "queue"))
    expect_identical(g$glow_slurm_queue()$tasks, 0L)
    writeLines("squeue: something that is not a queue line", file.path(state, "queue"))
    expect_error(g$glow_slurm_queue(), "squeue printed a line that is not a job id and a job name")
  })
  # (review 2026-11-16, second round) A squeue that failed was read as an empty queue.
  local({
    local_scheduler(squeue_status = 2L)
    expect_error(g$glow_slurm_queue(), "squeue failed \\(exit status 2\\): .*Unable to contact slurm controller.*submits nothing")
  })
  script <- file.path(withr::local_tempdir(), "job.sh"); writeLines("#SBATCH --job-name=glow-x", script)
  local({ local_scheduler(); expect_identical(g$glow_slurm_submit(c("--parsable", script), "glow-x"), "7001") })
  local({ local_scheduler(sbatch_says = "7001;cluster"); expect_identical(g$glow_slurm_submit(c("--parsable", script), "glow-x"), "7001") })
  # an sbatch that printed a number and failed was read as a submitted job
  local({ local_scheduler(sbatch_status = 7L); expect_error(g$glow_slurm_submit(c("--parsable", script), "glow-x"), "sbatch failed for glow-x \\(exit status 7\\)") })
  local({ local_scheduler(sbatch_says = "Submitted batch job 7001")
          expect_error(g$glow_slurm_submit(c("--parsable", script), "glow-x"), "sbatch did not return a job id for glow-x") })
  local({ local_scheduler(sbatch_says = "7001abc"); expect_error(g$glow_slurm_submit(c("--parsable", script), "glow-x"), "did not return a job id") })
})

test_that("a dry run of the SLURM mode shows what would be submitted, after what, and what waits", {
  skip_without_slurm_mode()
  g <- slurm_env(); h <- helper_env(); parts <- slurm_parts()
  out <- withr::local_tempdir()
  st <- slurm_study(out, list(snv_set = "  region_type: [gene, window]"),
                    run = c("  slurm:", "    options: \"--partition=short --account=lab\""))
  local_scheduler()                                    # an empty queue of the test's own, not the queue of the person who runs the tests
  dry <- function(...) slurm_run(g, ..., "--config", st$file, "--dry-run")
  hex <- "[0-9a-f]{8}"
  shown <- dry("data_prep")
  line <- function(pattern) grep(pattern, shown, value = TRUE)
  # the first part is ready; a unit that reads it would be submitted after its job, chromosome by chromosome
  expect_match(line("^  plink_to_gds chr21 "), paste0("would submit: sbatch '--parsable' '--job-name=glow-mine-data_prep-plink_to_gds-chr21-", hex,
                                                      "' '--output=.*/logs/slurm/glow-mine-data_prep-plink_to_gds-chr21-", hex,
                                                      "-%j.log' '--cpus-per-task=1' '--mem=[^']*' '--time=[^']*' '--partition=short' '--account=lab' ",
                                                      "'.*/jobs/glow-mine-data_prep-plink_to_gds-chr21-", hex, ".sh'"))
  expect_match(line("^  annotate_favor chr21 "), paste0("'--kill-on-invalid-dep=yes' '--dependency=afterok:<glow-mine-data_prep-plink_to_gds-chr21-", hex, ">' '--partition=short'"))
  expect_false(grepl("plink_to_gds-chr22", line("^  annotate_favor chr21 ")))
  # (review 2026-11-16, second round) The phenotype assembly reads the annotated file of one
  # chromosome (sample_gds_chr, 22). It was made to wait for every chromosome.
  expect_match(line("^  assemble_pheno "), paste0("'--dependency=afterok:<glow-mine-data_prep-annotate_favor-chr22-", hex, ">'"))
  expect_false(grepl("annotate_favor-chr21", line("^  assemble_pheno ")))
  expect_true(any(grepl("dry run: the job scripts are written, nothing is submitted", shown)))
  # the job scripts are written, each in the jobs folder of its step, with the checksums of the configuration beside it
  jobs <- list.files(file.path(out, "runs", "mine", "baseline", "data_prep", "jobs"))
  expect_true(any(grepl(paste0("^glow-mine-baseline-data_prep-assemble_pheno-", hex, "\\.sh$"), jobs)))
  expect_true(any(grepl(paste0("^glow-mine-baseline-data_prep-assemble_pheno-", hex, "\\.config\\.md5$"), jobs)))
  # (review 2026-11-16, second round) A part of another step that is not finished cannot be waited
  # for by a job of this call. pi_train was submitted after pi_annotate alone.
  shown <- dry("training")
  expect_match(line("^  pi_annotate "), "afterok:<glow-mine-baseline-training-pi_prepare-")     # it reads the step's own part only
  expect_match(line("^  pi_train "), "waits for data_prep/annotate_favor, which this call does not run")
  expect_false(any(grepl("would submit", line("^  pi_(train|evaluate) "))))
  expect_match(line("^  pi_evaluate "), "waits")
  expect_match(line("^  b_estimate "), "afterok:<glow-mine-baseline-training-b_prepare-")       # what does not read it goes on
  shown <- dry("training", "--part", "pi_train")
  expect_match(line("^  pi_train "), "waits for .*training/pi_annotate.*, which this call does not run")
  expect_match(line("^  pi_train "), "waits for .*data_prep/annotate_favor.*, which this call does not run")
  # a window scan before its preparation: its units are not known, and what reads them waits
  shown <- dry("snv_set")
  expect_match(line("^  scan \\(window\\) "), "its units are known after the part before it has run")
  expect_match(line("^  aggregate window "), "waits for the units of scan")
  expect_match(line("^  plots window "), "waits")
  # after the data preparation and the training have finished, the gene run is one chain
  local_cfg <- study(out, list(snv_set = "  region_type: [gene, window]"), name = "local.yml")
  finish_upstream(h, local_cfg)
  shown <- dry("snv_set")
  expect_match(line("^  prepare gene "), "would submit: sbatch '--parsable' '--job-name=[^']*' '--output=[^']*' '--cpus-per-task=[0-9]+' '--mem=[^']*' '--time=[^']*' '--partition=short'")
  expect_match(line("^  scan gene chr22 "), paste0("afterok:<glow-mine-baseline-snv_set-prepare-gene-", hex, ">'"))
  expect_match(line("^  aggregate gene "), paste0("afterok:<glow-mine-baseline-snv_set-prepare-gene-", hex, ">:<glow-mine-baseline-snv_set-scan-gene-chr21-", hex,
                                                  ">:<glow-mine-baseline-snv_set-scan-gene-chr22-", hex, ">'"))
  expect_match(line("^  plots gene "), "afterok:<glow-mine-baseline-snv_set-aggregate-gene-")
  script <- list.files(file.path(out, "runs", "mine", "baseline", "snv_set", "gene", "jobs"), pattern = "scan-gene-chr21-.*\\.sh$", full.names = TRUE)
  expect_match(readLines(script), "\"\\$RSCRIPT\" '03-snv-set/run-gene.R' '--config' .* '--region-type' 'gene' '--chr' '21'", all = FALSE)
  # an aggregation that may run on an unfinished scan runs after the scan jobs, whatever their end
  lenient <- slurm_study(out, list(snv_set = c("  region_type: [gene]", "  allow_incomplete: true")), name = "lenient.yml")
  shown <- slurm_run(g, "snv_set", "--config", lenient$file, "--dry-run")
  expect_match(line("^  aggregate gene "), paste0("afterok:<glow-mine-baseline-snv_set-prepare-gene-", hex, ">,afterany:<glow-mine-baseline-snv_set-scan-gene-chr21-", hex,
                                                  ">:<glow-mine-baseline-snv_set-scan-gene-chr22-", hex, ">'"))
})

test_that("the units of a table that is about to be written again are not planned from the old table", {
  skip_without_slurm_mode()
  g <- slurm_env(); h <- helper_env()
  out <- withr::local_tempdir()
  st <- slurm_study(out, list(snv_set = "  region_type: [window]"))
  local_cfg <- study(out, list(snv_set = "  region_type: [window]"), name = "local.yml")
  finish_upstream(h, local_cfg)
  finish(h, "snv_set", "prepare", c("--config", local_cfg, "--region-type", "window"))
  folder <- file.path(out, "runs", "mine", "baseline", "snv_set", "window")
  dir.create(file.path(folder, "shared"), showWarnings = FALSE)
  saveRDS(data.frame(chunk_id = 1:2), file.path(folder, "shared", "01-genome_chunks.rds"))
  local_scheduler()
  shown <- slurm_run(g, "snv_set", "--config", st$file, "--dry-run")
  line <- function(pattern) grep(pattern, shown, value = TRUE)
  # the preparation is finished: its table is the one of this run, and its two chunks are planned
  expect_match(line("^  prepare window "), "complete")
  expect_match(line("^  scan window chunk0001 "), "would submit")
  expect_match(line("^  scan window chunk0002 "), "would submit")
  expect_match(line("^  aggregate window "), "afterok:<[^>]*chunk0001-[0-9a-f]{8}>:<[^>]*chunk0002-")
  # (review 2026-11-16, second round) A preparation that runs again may write another table, so
  # the chunks of the old one are not submitted with it.
  shown <- slurm_run(g, "snv_set", "--config", st$file, "--dry-run", "--force")
  expect_match(line("^  prepare window "), "would submit")
  expect_false(any(grepl("chunk000", shown)))
  expect_match(line("^  scan \\(window\\) "), "its units are known after the part before it has run")
  expect_match(line("^  aggregate window "), "waits for the units of scan")
})

test_that("a unit in the queue is left alone and followed, whatever its record says", {
  skip_without_slurm_mode(); skip_on_os("windows")
  g <- slurm_env(); h <- helper_env(); parts <- slurm_parts()
  out <- withr::local_tempdir()
  st <- slurm_study(out, run = c("  slurm:", "    poll_seconds: 1"))
  local_cfg <- study(out, name = "local.yml")
  for (k in c("21", "22")) finish(h, "data_prep", "plink_to_gds", c("--config", local_cfg, "--chr", k))
  state <- local_scheduler()
  # everything of the part is finished and nothing is queued
  shown <- slurm_run(g, "data_prep", "--part", "plink_to_gds", "--config", st$file)
  expect_true(any(grepl("every unit is finished", shown))); expect_identical(submissions(state), 0L)
  # (review 2026-11-16, second round) A forced rerun of chr21 is in the queue. Its old record is
  # still there, and the call said that every unit was finished.
  name21 <- job_name(g, st, parts, "data_prep.plink_to_gds", "chr21")
  writeLines(paste("6500", name21), file.path(state, "queue"))
  shown <- slurm_run(g, "data_prep", "--part", "plink_to_gds", "--config", st$file)
  expect_match(grep("^  plink_to_gds chr21 ", shown, value = TRUE), "in the queue as job 6500")
  expect_false(any(grepl("every unit is finished", shown)))
  expect_true(any(grepl("0 job\\(s\\) submitted, 1 already in the queue", shown)))
  expect_identical(submissions(state), 0L)                                    # and it is not submitted a second time
  # a unit that reads it would be out of date when it has run, so it follows the queued job
  for (k in c("21", "22")) finish(h, "data_prep", "annotate_favor", c("--config", local_cfg, "--chr", k))
  shown <- slurm_run(g, "data_prep", "--config", st$file)
  expect_match(grep("^  annotate_favor chr21 ", shown, value = TRUE), "submitted as job 7001, afterok:6500")
  expect_match(grep("^  annotate_favor chr22 ", shown, value = TRUE), "complete")
  args <- readLines(file.path(state, "args-7001"))
  expect_true(all(c("--parsable", "--kill-on-invalid-dep=yes", "--dependency=afterok:6500") %in% args))
  expect_match(args[2L], "^--job-name=glow-mine-data_prep-annotate_favor-chr21-[0-9a-f]{8}$")
  expect_match(args[3L], "^--output=.*/logs/slurm/glow-mine-data_prep-annotate_favor-chr21-[0-9a-f]{8}-%j.log$")
  # a second call finds both in the queue and submits nothing
  before <- submissions(state)
  shown <- slurm_run(g, "data_prep", "--config", st$file)
  expect_match(grep("^  annotate_favor chr21 ", shown, value = TRUE), "in the queue as job 7001")
  expect_identical(submissions(state), before)
})

test_that("a scheduler that fails stops the call, and a second call at the same time waits its turn", {
  skip_without_slurm_mode(); skip_on_os("windows")
  g <- slurm_env(); parts <- slurm_parts()
  out <- withr::local_tempdir()
  st <- slurm_study(out)
  # (review 2026-11-16, second round) A squeue that failed was read as an empty queue, and work was submitted.
  local({
    state <- local_scheduler(squeue_status = 2L)
    expect_error(slurm_run(g, "data_prep", "--part", "plink_to_gds", "--config", st$file), "squeue failed \\(exit status 2\\)")
    expect_identical(submissions(state), 0L)
  })
  # an sbatch that fails: the call stops, says what is queued already, and no later unit depends on the failed one
  local({
    state <- local_scheduler(sbatch_status = 7L)
    err <- tryCatch(slurm_run(g, "data_prep", "--config", st$file), error = conditionMessage)
    expect_match(err, "sbatch failed for glow-mine-data_prep-plink_to_gds-chr21-[0-9a-f]{8} \\(exit status 7\\)")
    expect_match(err, "0 job\\(s\\) were submitted before it")
    expect_identical(submissions(state), 1L)                                  # it stopped at the first failure
  })
  # the reservation of the dataset is released after a call that stopped, and held while another call submits
  lay <- glow_layout(st$cfg)
  expect_false(GLOWpipeline:::.rr_lock_state(file.path(lay$dataset_dir, "records", "slurm-submission.lock"))$held)
  local({
    state <- local_scheduler()
    other_call <- glow_lock(lay$dataset_dir, "slurm-submission")
    expect_error(slurm_run(g, "data_prep", "--part", "plink_to_gds", "--config", st$file), "Another run holds `slurm-submission`")
    expect_identical(submissions(state), 0L)
    glow_unit_release(other_call)
    shown <- slurm_run(g, "data_prep", "--part", "plink_to_gds", "--config", st$file)
    expect_identical(submissions(state), 2L)
    expect_true(any(grepl("2 job\\(s\\) submitted, 0 already in the queue", shown)))
    # a dry run reads and writes without the reservation
    held <- glow_lock(lay$dataset_dir, "slurm-submission")
    expect_true(any(grepl("dry run", slurm_run(g, "data_prep", "--config", st$file, "--dry-run"))))
    glow_unit_release(held)
  })
  # an option that the command sets itself stops the call before anything is submitted
  renamed <- slurm_study(out, name = "renamed.yml", run = c("  slurm:", "    options: \"--job-name=renamed\""))
  local({
    state <- local_scheduler()
    expect_error(slurm_run(g, "data_prep", "--config", renamed$file), "which the command sets itself")
    expect_identical(submissions(state), 0L)
  })
})

test_that("the cap counts every task of the user and a dry run, and --wait waits for room and judges its jobs", {
  skip_without_slurm_mode(); skip_on_os("windows")
  g <- slurm_env(); parts <- slurm_parts()
  out <- withr::local_tempdir()
  capped <- slurm_study(out, run = c("  slurm:", "    max_queued: 1", "    poll_seconds: 1"))
  # (review 2026-11-16, second round) A dry run showed two submissions under a cap of one.
  local({
    local_scheduler()
    shown <- slurm_run(g, "data_prep", "--part", "plink_to_gds", "--config", capped$file, "--dry-run")
    expect_match(grep("^  plink_to_gds chr21 ", shown, value = TRUE), "would submit")
    expect_match(grep("^  plink_to_gds chr22 ", shown, value = TRUE), "held back by run.slurm.max_queued")
  })
  # an array of another study fills the cap task by task
  wide <- slurm_study(out, name = "wide.yml", run = c("  slurm:", "    max_queued: 20", "    poll_seconds: 1"))
  local({
    state <- local_scheduler()
    writeLines(sprintf("900_%d other", 1:20), file.path(state, "queue"))
    shown <- slurm_run(g, "data_prep", "--part", "plink_to_gds", "--config", wide$file)
    expect_identical(length(grep("held back by run.slurm.max_queued", shown)), 2L)
    expect_identical(submissions(state), 0L)
    writeLines(sprintf("900_%d other", 1:19), file.path(state, "queue"))
    shown <- slurm_run(g, "data_prep", "--part", "plink_to_gds", "--config", wide$file)
    expect_identical(submissions(state), 1L)
  })
  # (review 2026-11-16, second round) With the cap taken by another job, --wait stopped with
  # "nothing is queued". It waits for room. Here the other job leaves at the third look, the unit
  # is submitted, and its job ends without a finished unit, which --wait reports.
  local({
    state <- local_scheduler(finish_at_once = TRUE, free_after = 3L)
    writeLines("555 other", file.path(state, "queue"))
    shown <- character(0)
    err <- tryCatch(withr::with_dir(ga_dir(), shown <- capture.output(
      g$glow_run_main(c("run", "data_prep", "--part", "plink_to_gds", "--config", capped$file, "--wait")), type = "output")), error = conditionMessage)
    # (review 2026-11-16, second round) and a job that left the queue without its unit is an
    # error at once: the failed unit was submitted again and again by the same call.
    expect_match(err, "glow run --wait: 1 of 1 unit\\(s\\) left the queue without finishing: plink_to_gds chr21")
    expect_identical(submissions(state), 1L)
    expect_gte(as.integer(readLines(file.path(state, "looks"))), 3L)
  })
  # without a cap: one round of two jobs, which both end without a record, and no second round
  open <- slurm_study(out, name = "open.yml", run = c("  slurm:", "    poll_seconds: 1"))
  local({
    state <- local_scheduler(finish_at_once = TRUE)
    expect_error(slurm_run(g, "data_prep", "--part", "plink_to_gds", "--config", open$file, "--wait"),
                 "2 of 2 unit\\(s\\) left the queue without finishing")
    expect_identical(submissions(state), 2L)
  })
  # nothing can be submitted and nothing is queued: --wait says so instead of waiting for ever
  local({
    state <- local_scheduler()
    expect_error(slurm_run(g, "training", "--part", "pi_train", "--config", open$file, "--wait"),
                 "nothing is queued and 1 unit\\(s\\) cannot be submitted")
  })
})

test_that("--wait takes a forced unit as finished only by the record that its own job wrote", {
  skip_without_slurm_mode(); skip_on_os("windows")
  g <- slurm_env(); h <- helper_env(); parts <- slurm_parts()
  withr::local_envvar(SLURM_JOB_ID = NA)            # the tests may run inside a job of their own
  out <- withr::local_tempdir()
  st <- slurm_study(out, run = c("  slurm:", "    poll_seconds: 1"))
  call <- function(..., config = st$file) slurm_run(g, "training", "--part", "b_prepare", "--config", config, "--wait", ...)
  refused <- "glow run --wait: 1 of 1 unit\\(s\\) left the queue without finishing: b_prepare"
  # (review 2026-11-16, third round) The unit is finished. Its forced job leaves the queue
  # without having run (the scheduler of the test runs nothing), so the record is the old one.
  # --wait reported success, because a valid record was there.
  state <- local_scheduler(finish_at_once = TRUE)
  finish(h, "training", "b_prepare", c("--config", st$file))
  folder <- file.path(out, "runs", "mine", "baseline", "training", "b")
  record <- file.path(folder, "records", "b_prepare.yml")
  before <- readLines(record)
  expect_identical(length(glow_unit_record(folder, "b_prepare")$producer), 0L)   # a run outside a job names none
  expect_error(call("--force"), refused)                                         # job 7001
  expect_identical(readLines(record), before)
  # (review 2026-11-17) A record that another run wrote in the meantime changed the file, and
  # the changed file passed for the work of the job. The record has to name the job.
  withr::with_envvar(c(SLURM_JOB_ID = "6999"), finish(h, "training", "b_prepare", c("--config", st$file, "--force")))
  expect_identical(glow_unit_record(folder, "b_prepare")$producer, c(slurm_job = "6999"))
  expect_error(call("--force"), refused)                                         # job 7002
  # (review 2026-11-17) The record of the job itself counts whenever it was written, also before
  # the call looked at the queue: here it is in place when the call starts.
  withr::with_envvar(c(SLURM_JOB_ID = "7003"), finish(h, "training", "b_prepare", c("--config", st$file, "--force")))
  shown <- call("--force")                                                       # job 7003
  expect_true(any(grepl("every unit is finished", shown)))
  expect_identical(submissions(state), 3L)
  # a unit that is not forced is finished by any valid record: nothing is submitted for it
  expect_true(any(grepl("every unit is finished", call())))
  expect_identical(submissions(state), 3L)
  # (review 2026-11-17) A forced unit of a derived analysis is not finished by the folder that its
  # parent finished: with --force the unit does not reuse, so its job has to write a record.
  finish(h, "training", "b_estimate", c("--config", st$file))
  child <- derived(out, st$file, "child")
  shown <- call(config = child)                                                  # not forced: reused, no job
  expect_true(any(grepl("every unit is finished", shown))); expect_identical(submissions(state), 3L)
  expect_error(call("--force", config = child), refused)                         # job 7004 ran nothing
  expect_false(file.exists(file.path(out, "runs", "mine", "child", "training", "b", "records", "b_prepare.yml")))
})

test_that("a job script says whether it was forced, and a queued job of an earlier call is judged by that", {
  skip_without_slurm_mode(); skip_on_os("windows")
  g <- slurm_env(); h <- helper_env(); parts <- slurm_parts()
  withr::local_envvar(SLURM_JOB_ID = NA)
  out <- withr::local_tempdir()
  st <- slurm_study(out, run = c("  slurm:", "    poll_seconds: 1"))
  lay <- glow_layout(st$cfg); spec <- parts[["training.b_prepare"]]
  name <- job_name(g, st, parts, "training.b_prepare")
  p <- list(spec = spec, rt = NULL, unit = NULL, name = name)
  script_of <- function(argv) withr::with_dir(ga_dir(), g$glow_slurm_job_script(name, argv, spec, st$cfg,
    file.path(out, "runs", "mine", "baseline", "training", "b", "jobs"), file.path(out, "logs")))
  forced <- function() withr::with_dir(ga_dir(), g$glow_slurm_job_forced(lay, p))
  expect_true(forced())                                             # no script: the stricter reading
  path <- script_of(c("01-training/b/prepare-b-training-data.R", "--config", st$file))
  expect_true("# forced: no" %in% readLines(path)); expect_false(forced())
  path <- script_of(c("01-training/b/prepare-b-training-data.R", "--config", st$file, "--force"))
  expect_true("# forced: yes" %in% readLines(path)); expect_true(forced())
  # The forced job of an earlier call is in the queue (6500) and its unit has an old, valid record.
  # A call without --force follows it. The job leaves the queue without a record of its own.
  finish(h, "training", "b_prepare", c("--config", st$file))
  local({
    state <- local_scheduler(empty_after = 2L)
    writeLines(paste("6500", name), file.path(state, "queue"))                 # leaves at the second look
    expect_error(slurm_run(g, "training", "--part", "b_prepare", "--config", st$file, "--wait"),
                 "glow run --wait: 1 of 1 unit\\(s\\) left the queue without finishing: b_prepare")
    expect_identical(submissions(state), 0L)
  })
  # The same with the record that job 6500 wrote: finished, although the record was there first.
  withr::with_envvar(c(SLURM_JOB_ID = "6500"), finish(h, "training", "b_prepare", c("--config", st$file, "--force")))
  local({
    state <- local_scheduler(empty_after = 2L)
    writeLines(paste("6500", name), file.path(state, "queue"))
    shown <- slurm_run(g, "training", "--part", "b_prepare", "--config", st$file, "--wait")
    expect_true(any(grepl("every unit is finished", shown)))
    expect_identical(submissions(state), 0L)
  })
})

test_that("the checksums of a job are those of the files as they were read", {
  skip_without_slurm_mode(); skip_on_os("windows")
  g <- slurm_env()
  out <- withr::local_tempdir()
  st <- slurm_study(out, list(shared = "  n_cores: 1"))
  read_sum <- unname(tools::md5sum(st$file))
  # (review 2026-11-16, third round) The file is edited after it was read. The checksums were
  # taken at the submission, so the job passed its check and ran with settings other than
  # those its resources were asked for.
  writeLines(sub("n_cores: 1", "n_cores: 8", readLines(st$file), fixed = TRUE), st$file)
  parts <- slurm_parts()
  chain <- withr::with_dir(ga_dir(), g$glow_slurm_chain(st$cfg, attr(parts, "glow_md5")))
  mine <- grep(paste0("  ", normalizePath(st$file), "$"), chain, value = TRUE)
  expect_identical(sub("  .*$", "", mine), read_sum)
  expect_false(identical(read_sum, unname(tools::md5sum(st$file))))
  sums <- file.path(out, "chain.md5"); writeLines(chain, sums)
  expect_false(identical(suppressWarnings(system2("md5sum", c("--quiet", "-c", shQuote(sums)), stdout = FALSE, stderr = FALSE)), 0L))
  # the registry and the table of parts are in the list, each once
  expect_identical(length(grep("settings-registry.yml$|parts.yml$", chain)), 2L)
  expect_identical(anyDuplicated(sub("^[0-9a-f]+  ", "", chain)), 0L)
  # and a call stops before it submits anything
  expect_error(withr::with_dir(ga_dir(), g$glow_slurm_chain_current(chain)), "changed while the command was running: .*slurm.yml")
  local({
    state <- local_scheduler()
    stale <- st$cfg                                   # the configuration as read before the edit
    lay <- glow_layout(stale); parts <- slurm_parts()
    cmd <- list(step = "data_prep", part = "plink_to_gds", dry_run = FALSE, wait = FALSE, force = FALSE)
    of_step <- Filter(function(p) identical(p$id, "data_prep.plink_to_gds"), parts)
    expect_error(withr::with_dir(ga_dir(), g$glow_run_slurm(cmd, st$file, stale, lay, of_step, parts)), "changed while the command was running")
    expect_identical(submissions(state), 0L)
  })
  # (review 2026-11-17) A file was parsed first and hashed by its path afterwards, so an edit
  # between the two gave the settings of one file with the checksum of another. The checksum
  # is that of the bytes that were parsed. Each file here is edited while the parser works on it.
  edited_while_parsed <- function(file, read) {
    before <- unname(tools::md5sum(file))
    text <- rawToChar(readBin(file, "raw", file.size(file)))
    # the parser gets the text of this file (and of no other file of the call): the edit comes then
    suppressMessages(trace("yaml.load", where = asNamespace("yaml"), print = FALSE,
                           tracer = bquote(if (identical(as.character(string), .(text)))
                             cat("# edited while it was parsed\n", file = .(file), append = TRUE))))
    on.exit(suppressMessages(untrace("yaml.load", where = asNamespace("yaml"))), add = TRUE)
    got <- read(file)
    expect_false(identical(unname(tools::md5sum(file)), before))     # the file on disk is another one now
    expect_identical(got, before)                                    # and the checksum is that of what was read
  }
  copy <- function(name) { f <- file.path(out, name); file.copy(file.path(ga_dir(), "base-config", name), f); f }
  mine <- file.path(out, "edited.yml"); file.copy(st$file, mine)
  edited_while_parsed(mine, function(f) {
    chain <- attr(read_glow_config(f, base = file.path(ga_dir(), "base-config", "base_config.yml"),
                                   registry = file.path(ga_dir(), "base-config", "settings-registry.yml"), root = ga_dir()), "glow_chain")
    chain$md5[chain$file == normalizePath(f)]
  })
  registry <- copy("settings-registry.yml")
  edited_while_parsed(registry, function(f) .gc_read_registry(f)$md5)
  table <- copy("parts.yml")
  edited_while_parsed(table, function(f) attr(glow_read_parts(f), "glow_md5"))
  # The table of parts is read once, by the command, against the registry of the configuration.
  # The SLURM mode gets that table and its checksum and does not read the file again.
  expect_identical(attr(parts, "glow_md5"), unname(tools::md5sum(file.path(ga_dir(), "base-config", "parts.yml"))))
  expect_identical(names(glow_read_parts(file.path(ga_dir(), "base-config", "parts.yml"), attr(st$cfg, "glow_registry"))), names(parts))
  expect_false(any(grepl("glow_read_parts(\"", readLines(file.path(ga_dir(), "lib", "glow-slurm.R")), fixed = TRUE)))
  expect_error(withr::with_dir(ga_dir(), g$glow_slurm_chain(st$cfg, NULL)), "the checksum of the table of parts as it was read")
})

test_that("sbatch gets what the command sets on its command line and none of it from the environment", {
  skip_without_slurm_mode(); skip_on_os("windows")
  g <- slurm_env(); parts <- slurm_parts()
  out <- withr::local_tempdir()
  st <- slurm_study(out, run = c("  slurm:", "    resources:", "      plink_to_gds: {mem: 3G, time: \"00:20:00\"}"))
  # (review 2026-11-16, third round) sbatch reads an option from the environment before the
  # #SBATCH lines. With SBATCH_JOB_NAME set, every job had that name, and a second call did
  # not find its jobs in the queue and submitted them again.
  withr::local_envvar(SBATCH_JOB_NAME = "site-default", SBATCH_ARRAY_INX = "1-3", SBATCH_TIMELIMIT = "1", SBATCH_PARTITION = "short")
  state <- local_scheduler()
  slurm_run(g, "data_prep", "--part", "plink_to_gds", "--config", st$file)
  expect_identical(submissions(state), 2L)
  slurm_run(g, "data_prep", "--part", "plink_to_gds", "--config", st$file)
  expect_identical(submissions(state), 2L)                                 # the second call found both in the queue
  name <- job_name(g, st, parts, "data_prep.plink_to_gds", unit = "chr21")
  expect_true(any(grepl(paste0(" ", name, "$"), readLines(file.path(state, "queue")))))
  expect_false(any(grepl("site-default", readLines(file.path(state, "queue")))))
  args <- readLines(file.path(state, "args-7001"))
  expect_true(all(c(paste0("--job-name=", name), "--cpus-per-task=1", "--mem=3G", "--time=00:20:00") %in% args))
  seen <- readLines(file.path(state, "env-7001"))
  expect_identical(seen, "SBATCH_PARTITION=short")                         # a variable of the site stays, the others are gone
  # the alias of the cores is among the options that the command sets itself
  expect_error(g$glow_slurm_options("--tres-per-task=cpu=8"), "which the command sets itself")
})

test_that("a scheduler that does not answer in time stops the call, and settings of the mode are checked", {
  skip_without_slurm_mode(); skip_on_os("windows")
  g <- slurm_env()
  out <- withr::local_tempdir()
  st <- slurm_study(out)
  local({
    state <- local_scheduler()
    bin <- dirname(state)
    writeLines(c("#!/bin/bash", "sleep 20"), file.path(bin, "squeue")); Sys.chmod(file.path(bin, "squeue"), "755")
    withr::local_options(glow.slurm.command_timeout = 1)
    started <- Sys.time()
    expect_error(slurm_run(g, "data_prep", "--part", "plink_to_gds", "--config", st$file), "squeue failed .*no answer within 1 seconds")
    expect_lt(as.numeric(difftime(Sys.time(), started, units = "secs")), 15)
    expect_identical(submissions(state), 0L)
    withr::local_options(glow.slurm.command_timeout = 0)
    expect_error(g$glow_slurm_queue(), "glow.slurm.command_timeout has to be a positive number")
  })
  local({
    state <- local_scheduler()
    for (bad in list(c("    max_queued: 0", "`run.slurm.max_queued` has to be 1 or more"),
                     c("    max_queued: -1", "`run.slurm.max_queued` has to be 1 or more"),
                     c("    poll_seconds: 0", "`run.slurm.poll_seconds` has to be a positive number"),
                     c("    resources: {plink_to_gdss: {mem: 4G}}", "names 1 part\\(s\\) that the workflow does not have: plink_to_gdss"))) {
      wrong <- slurm_study(out, name = "wrong.yml", run = c("  slurm:", bad[1L]))
      expect_error(slurm_run(g, "data_prep", "--part", "plink_to_gds", "--config", wrong$file), bad[2L])
    }
    expect_identical(submissions(state), 0L)
  })
})

