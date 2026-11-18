# ==============================================================================
# File: glow_run_records.R
# Purpose: The records of a run on the derived layout. A unit of work (one run of
#   a part of a step, or one chromosome or chunk of a part that is split) writes a
#   small completion record when it has finished: what it depended on, which files
#   it produced, and a generation mark that is new at every run. The records make
#   three things possible: a rerun continues from the first unit that is not
#   complete, a run that would overwrite results made with other settings is
#   refused, and an analysis can reuse a finished part of the analysis it extends.
#
# EXPORTED FUNCTIONS:
#   - glow_read_parts()                The table of parts and what each depends on
#   - glow_setting_text()              A setting's value as one line of text
#   - glow_file_fingerprint()          What identifies an outside file or folder
#   - glow_run_versions()              The package versions a run depends on
#   - glow_unit_record()               The completion record of a unit, or NULL
#   - glow_unit_state()                Whether a unit is complete, and if not, why
#   - glow_dependency_differences()    Recorded against current dependencies
#   - glow_unit_begin()                Lock a unit and decide: run, or complete
#   - glow_unit_complete()             Write the record and release the lock
#   - glow_unit_release()              Release the lock of a unit that did not finish
#   - glow_lock()                      A named lock, for the submission of the SLURM mode
#   - glow_seed_unit()                 Supply the files of a unit instead of running it
#
# INTERNAL HELPERS:
#   - .rr_stem(), .rr_record_file(), .rr_lock_dir()   Names in the records folder
#   - .rr_write_record(), .rr_yaml_lines()            The record as YAML text
#   - .rr_named_text()                 A map of the record as a named character vector
#   - .rr_take_lock(), .rr_release_lock(), .rr_lock_is_stale()   The lock of a unit
#   - .rr_lock_generations(), .rr_lock_state(), .rr_lock_is_mine()   The generations of a lock
#   - .rr_read_owner(), .rr_owner_field(), .rr_process_start()   The owner of a lock
#   - .rr_retire_outputs(), .rr_relative(), .rr_text_md5()       Recorded files
#   - .rr_generation()                 A mark that is new at every run
#

#################### EXPORTED MAIN FUNCTIONS ####################

#' Read the table of parts
#'
#' A step of a workflow is a sequence of parts, and each part declares what it
#' depends on: the settings it reads, the settings that reach it through a
#' forwarded configuration, the outputs of other parts, and the outside files. The
#' table is a YAML file with one entry per part. This function reads and checks it.
#'
#' @param file Path of the table (a YAML list of mappings).
#' @param registry `NULL`, the path of the settings registry, or the registry as
#'   [read_glow_config()] read it (the attribute `glow_registry` of a
#'   configuration). When given, every setting that the table names has to be
#'   registered. The registry of a configuration that was already read is the
#'   one to pass, so that the file is not read a second time.
#'
#' @return A list named by `step.part`, which identifies a part, because two
#'   steps may each have a part of the same name. The attribute `glow_md5` holds
#'   the MD5 checksum of the bytes of `file` that were parsed. Each element has `id` (that
#'   name), `part`, `step`, `script`,
#'   `folder` (a key of the layout, such as `steps.training_b`), `level`
#'   (`dataset`, `analysis` or `region_type`), `units` (`none`, `chromosome`,
#'   `chunk`, `region`, `group_mode` or `chunk_mode`), the character vectors
#'   `settings` and `forwarded`, the list `region_settings` (a character vector of
#'   settings per region type) and the lists `inputs`, `files` and `run_when`. An
#'   input is a list with `part` (the `step.part` of the part that is read), and
#'   `setting`, `unit`, `units_from`, `pattern`, `partial` and `when`, each of
#'   which may be `NULL` (`when` is a list of conditions, empty when there is none). An entry of
#'   `files` is a list with `setting`, `pattern` and `when`.
#'
#' @details
#' An entry may use the keys `part`, `step`, `script`, `call`, `folder`, `level`,
#' `units`, `settings`, `region_settings`, `forwarded`, `files`, `inputs`,
#' `run_when`, `slurm` and `description`.
#' `units: region` is the unit of an SNV-set scan, a chromosome or a chunk by the
#' region type. `region_settings` maps a region type to the settings that only a
#' run of that type reads. `inputs` is a list of mappings with `part` (the part
#' whose output is read: a part of the same step by its name, a part of another
#' step as `step.part`) and, optionally, `setting` (the setting that holds the
#' path, which a user may set to override the derivation), `unit` (the one unit
#' that is read, where a setting in braces is replaced by its value, as in
#' `chr{data_prep.sample_gds_chr}`), `units_from` (a setting that lists the
#' chromosomes whose units are read, when they are not the chromosomes of the
#' run), `pattern` (the setting that holds the file-name pattern of the files that
#' are read from a folder the user set, so that only those files identify the
#' input), `partial` (a logical setting that, when true, lets the part run on the
#' units of the input that are finished) and `when` (conditions that all have to
#' hold for the input to be read).
#' An entry of `files` is a setting that names an outside file or folder, or a
#' mapping with `setting` and, optionally, `pattern` and `when`. A condition, of
#' `when` and of `run_when`, is a mapping with `setting` and either `is` or
#' `is_not`, whose value is `null`, `set` or the text of a value, or a mapping
#' with `region_type` and the region types it holds for. A setting is named
#' `section.name`, with a dot and a field for a value inside a map
#' (`snv_set.calibration.report`). The function stops on an unknown key, a part
#' named twice, an input that names no part of the table, and, with a registry, a
#' setting that is not registered. Computational complexity: linear in the size of
#' the table.
#'
#' @examples
#' \dontrun{
#' parts <- glow_read_parts("base-config/parts.yml", "base-config/settings-registry.yml")
#' parts$b_estimate$settings
#' }
#'
#' @seealso [glow_unit_begin()]
#' @export
glow_read_parts <- function(file, registry = NULL) {
  stopifnot(is.character(file), length(file) == 1L)
  if (!file.exists(file)) stop("The table of parts was not found: ", file, call. = FALSE)
  tree <- .gc_parse_file(file, top = "seq")
  raw <- .gc_plain(tree)
  keys <- c("part", "step", "script", "call", "folder", "level", "units", "settings", "region_settings",
            "forwarded", "files", "inputs", "run_when", "slurm", "description")
  levels_ok <- c("dataset", "analysis", "region_type")
  units_ok  <- c("none", "chromosome", "chunk", "region", "group_mode", "chunk_mode")
  text_vec <- function(x) if (is.null(x) || length(x) == 0L) character(0) else as.character(unlist(x))
  one_text <- function(x) is.character(x) && length(x) == 1L
  # a list of conditions: each names a setting and what it has to be (`is`) or must not
  # be (`is_not`), or the region types it holds for
  conditions <- function(x, where, what) {
    out <- list()
    for (cond in x) {
      ok <- is.list(cond) && !is.null(names(cond)) &&
        (identical(names(cond), "region_type") ||
           (length(cond) == 2L && "setting" %in% names(cond) && any(c("is", "is_not") %in% names(cond)) &&
              all(vapply(cond, one_text, logical(1)))))
      if (!ok) stop(where, ": a condition of `", what, "` is a mapping with `setting` and `is` (or `is_not`), ",
                    "or with `region_type`.", call. = FALSE)
      out[[length(out) + 1L]] <- if (!is.null(cond$region_type)) list(region_type = text_vec(cond$region_type))
                                 else cond[c("setting", intersect(c("is", "is_not"), names(cond)))]
    }
    out
  }
  parts <- list()
  for (i in seq_along(raw)) {
    e <- raw[[i]]
    where <- sprintf("%s, entry %d", file, i)
    if (!is.list(e) || is.null(names(e))) stop(where, ": an entry has to be a mapping.", call. = FALSE)
    unknown <- setdiff(names(e), keys)
    if (length(unknown)) stop(where, ": unknown key(s) ", paste(unknown, collapse = ", "), ".", call. = FALSE)
    for (need in c("part", "step", "script", "folder", "level", "units"))
      if (!one_text(e[[need]])) stop(where, ": `", need, "` has to be one text value.", call. = FALSE)
    id <- paste(e$step, e$part, sep = ".")
    if (id %in% names(parts)) stop(where, ": the part `", e$part, "` of the step `", e$step, "` is named twice.", call. = FALSE)
    if (!e$level %in% levels_ok) stop(where, ": `level` is one of ", paste(levels_ok, collapse = ", "), ".", call. = FALSE)
    if (!e$units %in% units_ok) stop(where, ": `units` is one of ", paste(units_ok, collapse = ", "), ".", call. = FALSE)
    inputs <- list()
    for (inp in e$inputs) {
      if (!is.list(inp) || !one_text(inp[["part"]]) ||
          length(setdiff(names(inp), c("part", "setting", "unit", "units_from", "pattern", "partial", "when"))))
        stop(where, ": an input has to be a mapping with `part` and, optionally, `setting`, `unit`, `units_from`, ",
             "`pattern`, `partial` and `when`.", call. = FALSE)
      # a part of the same step is named without the step
      inputs[[length(inputs) + 1L]] <- list(
        part = if (grepl(".", inp[["part"]], fixed = TRUE)) inp[["part"]] else paste(e$step, inp[["part"]], sep = "."),
        # by exact name: `$unit` would also match `units_from`, and `$part` `partial`
        setting = inp[["setting"]], unit = inp[["unit"]], units_from = inp[["units_from"]],
        pattern = inp[["pattern"]], partial = inp[["partial"]], when = conditions(inp[["when"]], where, "when"))
    }
    files <- list()
    for (f in e$files) {
      if (one_text(f)) f <- list(setting = f)
      if (!is.list(f) || !one_text(f[["setting"]]) || length(setdiff(names(f), c("setting", "pattern", "when"))))
        stop(where, ": an entry of `files` is a setting, or a mapping with `setting` and, optionally, `pattern` ",
             "and `when`.", call. = FALSE)
      files[[length(files) + 1L]] <- list(setting = f[["setting"]], pattern = f[["pattern"]],
                                          when = conditions(f[["when"]], where, "when"))
    }
    region_settings <- list()
    if (!is.null(e$region_settings)) {
      if (!is.list(e$region_settings) || is.null(names(e$region_settings)))
        stop(where, ": `region_settings` has to be a mapping from a region type to settings.", call. = FALSE)
      region_settings <- lapply(e$region_settings, text_vec)
    }
    if (!is.null(e$slurm) && (!is.list(e$slurm) || length(setdiff(names(e$slurm), c("mem", "time", "cpus")))))
      stop(where, ": `slurm` is a mapping with `mem`, `time` and `cpus`.", call. = FALSE)
    parts[[id]] <- list(id = id, part = e$part, step = e$step, script = e$script, call = text_vec(e$call),
                        run_when = conditions(e$run_when, where, "run_when"), slurm = e$slurm,
                        folder = e$folder, level = e$level, units = e$units, settings = text_vec(e$settings),
                        region_settings = region_settings, forwarded = text_vec(e$forwarded),
                        files = files, inputs = inputs)
  }
  input_parts <- function(p) vapply(p$inputs, `[[`, character(1), "part")
  for (p in parts) {
    absent <- setdiff(input_parts(p), names(parts))
    if (length(absent))
      stop(file, ": the part `", p$id, "` reads the output of `", paste(absent, collapse = "`, `"),
           "`, which the table does not hold.", call. = FALSE)
  }
  if (!is.null(registry)) {
    registered <- names((if (is.list(registry)) registry else .gc_read_registry(registry))$entries)
    # `top.dataset` in the registry is the top-level setting `dataset`
    registered <- c(registered, sub("^top\\.", "", registered[startsWith(registered, "top.")]))
    of_conditions <- function(conds) unlist(lapply(conds, `[[`, "setting"))
    in_braces <- function(x) if (is.null(x)) NULL else regmatches(x, gregexpr("(?<=\\{)[^}]+(?=\\})", x, perl = TRUE))[[1]]
    for (p in parts) {
      named <- unique(c(p$settings, unlist(p$region_settings), p$forwarded, of_conditions(p$run_when),
                        unlist(lapply(p$files, function(f) c(f$setting, f$pattern, of_conditions(f$when)))),
                        unlist(lapply(p$inputs, function(i) c(i$setting, i$units_from, i$pattern, i$partial, of_conditions(i$when),
                                                              in_braces(i$unit))))))
      # a field inside a map: the setting is the first two parts of the name
      setting <- vapply(strsplit(named, ".", fixed = TRUE), function(k) paste(utils::head(k, 2L), collapse = "."), character(1))
      bad <- named[!(setting %in% registered | named %in% registered)]
      if (length(bad))
        stop(file, ": the part `", p$id, "` names setting(s) that are not registered: ",
             paste(bad, collapse = ", "), ".", call. = FALSE)
    }
  }
  attr(parts, "glow_md5") <- attr(tree, "glow_md5")   # the checksum of the bytes that were parsed
  parts
}

#' A setting's value as one line of text
#'
#' The dependencies of a unit are stored and compared as text. This function gives
#' the text of one value: the R expression that rebuilds it, with every digit of a
#' number and with integers marked as such, so that two values have the same text
#' exactly when they are the same value of the same type, in the same order.
#'
#' @param value The value of a setting, as [read_glow_config()] returns it.
#'
#' @return One string. A function gives its source text as the configuration wrote
#'   it. A null gives `NULL`.
#'
#' @details Settings are compared after the conversion by their registered type,
#'   never as the YAML text a user wrote, so `1e6` and `1000000` are the same value
#'   of a numeric setting. Computational complexity: linear in the size of the value.
#'
#' @examples
#' glow_setting_text(c(21L, 22L))
#' glow_setting_text(list(sex = "sex", age = "age"))
#'
#' @export
glow_setting_text <- function(value) {
  if (is.function(value)) {
    src <- attr(value, "glow_source")
    if (is.null(src)) src <- paste(deparse(value, control = c("keepNA", "keepInteger")), collapse = " ")
    return(paste("function text:", paste(src, collapse = " ")))
  }
  paste(deparse(value, width.cutoff = 500L,
                control = c("keepNA", "keepInteger", "niceNames", "showAttributes", "digits17")),
        collapse = " ")
}

#' What identifies an outside file or folder
#'
#' A unit depends on files that no part of the workflow produced, such as a
#' phenotype table. The record keeps a fingerprint of each one, and a changed
#' fingerprint means changed input.
#'
#' @param path One path, several paths (a set of files that are read together,
#'   such as the per-chromosome files of a folder), or `NULL`.
#' @param md5_limit Numeric. The largest file, in bytes, whose content is read for
#'   a checksum. Default: 64 MB.
#' @param small Numeric. Inside a folder, the largest file, in bytes, whose
#'   content is read. Default: 1 MB.
#'
#' @return One string. For a file up to `md5_limit`: `md5:` and its checksum. For
#'   a larger file: its size and modification time, because reading whole-genome
#'   files for every unit would cost more than the unit. For a folder: the number
#'   of files and a checksum of their names with, per file, the checksum of a file
#'   up to `small` and the size and modification time of a larger one. For several
#'   paths: the number of files and a checksum of their names and fingerprints.
#'   For a path that does not exist: `missing`. For `NULL`: `none`.
#'
#' @details Size and modification time do not prove equal content. They detect a
#'   file that was replaced, which is the case the record has to catch, and they
#'   miss a replacement that keeps both. The modification time is taken to the
#'   microsecond. A folder is listed recursively. Computational complexity: linear
#'   in the size of the files that are read, and in the number of files of a folder.
#'
#' @examples
#' f <- tempfile(); writeLines("a", f)
#' glow_file_fingerprint(f)
#'
#' @export
glow_file_fingerprint <- function(path, md5_limit = 64 * 1024^2, small = 1024^2) {
  if (is.null(path)) return("none")
  stopifnot(is.character(path))
  if (length(path) != 1L) {
    if (!length(path)) return("set:0 files")
    path <- sort(path, method = "radix")                # the order of the C locale, in every session
    each <- vapply(path, glow_file_fingerprint, character(1), md5_limit = md5_limit, small = small)
    return(sprintf("set:%d files, md5:%s", length(path), .rr_text_md5(paste(basename(path), each, sep = "="))))
  }
  if (!file.exists(path)) return("missing")
  stamp <- function(info) sprintf("%.0f@%.6f", info$size, as.numeric(info$mtime))
  if (dir.exists(path)) {
    rel <- sort(list.files(path, recursive = TRUE, all.files = TRUE, no.. = TRUE), method = "radix")
    if (!length(rel)) return("dir:0 files")
    info <- file.info(file.path(path, rel))
    id <- stamp(info)
    read <- !is.na(info$size) & info$size <= small
    id[read] <- paste0("md5:", unname(tools::md5sum(file.path(path, rel[read]))))
    return(sprintf("dir:%d files, listing md5:%s", length(rel), .rr_text_md5(paste(rel, id, sep = "="))))
  }
  info <- file.info(path)
  if (info$size <= md5_limit) return(paste0("md5:", unname(tools::md5sum(path))))
  paste0("size@mtime:", stamp(info))
}

#' The package versions a run depends on
#'
#' A default of a package can change between versions, so the versions are part of
#' the dependencies of every unit.
#'
#' @param glowanalyses Character. The version of the workflow templates, which is
#'   not an R package and is passed by the caller.
#'
#' @return A named character vector: `GFisher`, `GLOWr`, `GLOWpipeline` and
#'   `GLOWanalyses`. A package that is not installed gives `not installed`.
#'
#' @examples
#' glow_run_versions("0.4.0")
#'
#' @export
glow_run_versions <- function(glowanalyses) {
  stopifnot(is.character(glowanalyses), length(glowanalyses) == 1L)
  pkg <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed")
  c(GFisher = pkg("GFisher"), GLOWr = pkg("GLOWr"), GLOWpipeline = pkg("GLOWpipeline"),
    GLOWanalyses = glowanalyses)
}

#' The completion record of a unit
#'
#' @param folder The folder of the part (its records are in `records/` there).
#' @param part The name of the part.
#' @param unit The name of the unit, such as `chr21`, or `NULL` for a part that is
#'   not split.
#'
#' @return `NULL` when the unit has no record. Otherwise a list with `part`,
#'   `unit`, `generation`, `completed`, `seeded` (logical), `outputs` (paths
#'   relative to `folder`, or absolute for a file outside it), `dependencies`,
#'   a list of the four named character vectors `settings`, `inputs`, `files` and
#'   `versions`, and `producer`, a named character vector that says which process
#'   wrote the record (empty when the record states none; see
#'   [glow_unit_complete()]).
#'
#' @examples
#' glow_unit_record(tempdir(), "scan", "chr21")   # NULL: no record
#'
#' @seealso [glow_unit_complete()]
#' @export
glow_unit_record <- function(folder, part, unit = NULL) {
  file <- .rr_record_file(folder, part, unit)
  if (!file.exists(file)) return(NULL)
  raw <- .gc_plain(.gc_parse_file(file))
  deps <- raw$dependencies
  list(part = raw$part, unit = raw$unit, generation = raw$generation, completed = raw$completed,
       seeded = identical(raw$seeded, "true"),
       outputs = if (length(raw$outputs)) as.character(unlist(raw$outputs)) else character(0),
       dependencies = list(settings = .rr_named_text(deps$settings), inputs = .rr_named_text(deps$inputs),
                           files = .rr_named_text(deps$files), versions = .rr_named_text(deps$versions)),
       producer = .rr_named_text(raw$producer))
}

#' Whether a unit is complete
#'
#' A unit is complete when its record exists and every output the record lists
#' exists. With `dependencies`, the record also has to state the same
#' dependencies.
#'
#' @inheritParams glow_unit_record
#' @param dependencies The current dependencies of the unit (see
#'   [glow_unit_begin()]), or `NULL` to check the record and its outputs only.
#'
#' @return A list: `complete` (logical), `reason` (a sentence when not complete,
#'   else an empty string), `record` (the record, or `NULL`), `differences`
#'   (see [glow_dependency_differences()]), `stale` (logical) and `missing` (the
#'   outputs of the record that do not exist). A unit is stale when its record
#'   differs from the current dependencies in its inputs only: a part it read was
#'   produced again since, with the same settings. Its results belong to an input
#'   that is gone, so it runs again without a question.
#'
#' @details The dependencies are compared whether or not the outputs exist, so a
#'   record with other settings is reported as such also when one of its files
#'   was deleted. For a seed record (see [glow_seed_unit()]) the settings that
#'   the record names are compared: `dependencies$setting_text`, when it is a
#'   function, is called with those names and returns their current text.
#'
#' @examples
#' glow_unit_state(tempdir(), "scan", "chr21")$reason
#'
#' @export
glow_unit_state <- function(folder, part, unit = NULL, dependencies = NULL) {
  none <- character(0)
  rec <- glow_unit_record(folder, part, unit)
  if (is.null(rec)) return(list(complete = FALSE, reason = "it has no completion record", record = NULL,
                                differences = none, stale = FALSE, missing = none))
  out <- ifelse(grepl("^/", rec$outputs), rec$outputs, file.path(folder, rec$outputs))
  gone <- rec$outputs[!file.exists(out)]
  diffs <- none
  if (!is.null(dependencies)) {
    current <- dependencies
    if (rec$seeded && is.function(dependencies$setting_text))
      current$settings <- dependencies$setting_text(names(rec$dependencies$settings))
    diffs <- glow_dependency_differences(rec$dependencies, current, seeded = rec$seeded)
  }
  stale <- length(diffs) > 0L && all(startsWith(diffs, "input "))
  reason <- if (length(diffs) && !stale) "its record states other dependencies than this run"
            else if (stale) "a part it read was produced again since it ran"
            else if (length(gone)) paste0("an output that its record lists is missing: ", paste(gone, collapse = ", "))
            else ""
  list(complete = !length(diffs) && !length(gone), reason = reason, record = rec, differences = diffs,
       stale = stale, missing = gone)
}

#' Recorded against current dependencies
#'
#' @param recorded The `dependencies` of a completion record.
#' @param current The dependencies of the current run, in the same form: a list
#'   of the named character vectors `settings`, `inputs`, `files` and `versions`.
#' @param seeded Logical. `TRUE` for the record of files that were supplied and
#'   not produced by a run (see [glow_seed_unit()]). Such a record states the
#'   settings the files were built with and nothing else, so only the settings
#'   are compared.
#'
#' @return A character vector with one line per difference, naming the group and
#'   the key with the recorded and the current value. Empty when the two are the
#'   same.
#'
#' @details Two dependencies are the same when they have the same keys and every
#'   value is the same text (see [glow_setting_text()]). A key on one side only is
#'   a difference. Long values are shortened in the message, not in the
#'   comparison.
#'
#' @examples
#' a <- list(settings = c(x = "1L"), inputs = character(0), files = character(0),
#'           versions = c(GLOWr = "0.3.6"))
#' b <- a; b$settings[["x"]] <- "2L"
#' glow_dependency_differences(a, b)
#'
#' @export
glow_dependency_differences <- function(recorded, current, seeded = FALSE) {
  short <- function(x) if (is.na(x)) "(absent)" else if (nchar(x) > 140L) paste0(substr(x, 1L, 137L), "...") else x
  # A number of a setting is kept with 17 significant digits, which is exact and reads
  # badly (0.90000000000000002). The message shows 15 digits (0.9), unless the two
  # values then look alike. A character value stands in double quotes in the text of a
  # setting (glow_setting_text()), and what is inside the quotes is shown as it is: a
  # path or a name can hold such digits too.
  readable <- function(a, b) {
    shorter <- function(x) {
      # The text of a function is shown as the configuration wrote it: its source can hold a
      # number inside single quotes, which the rule for double quotes below does not see.
      if (is.na(x) || startsWith(x, "function text:")) return(x)
      quotes <- gregexpr('"(?:[^"\\\\]|\\\\.)*"', x, perl = TRUE)
      quoted <- regmatches(x, quotes)[[1L]]
      outside <- regmatches(x, quotes, invert = TRUE)[[1L]]
      outside <- vapply(outside, function(part) {
        m <- gregexpr("-?[0-9]+\\.[0-9]{15,}(e[-+]?[0-9]+)?", part)
        regmatches(part, m) <- lapply(regmatches(part, m), function(tok)
          vapply(tok, function(t) format(as.numeric(t), digits = 15), character(1)))
        part
      }, character(1), USE.NAMES = FALSE)
      paste0(c(rbind(outside, c(quoted, ""))), collapse = "")
    }
    sa <- shorter(a); sb <- shorter(b)
    if (identical(sa, sb)) c(a, b) else c(sa, sb)
  }
  out <- character(0)
  for (group in if (seeded) "settings" else c("settings", "inputs", "files", "versions")) {
    r <- recorded[[group]]; cur <- current[[group]]
    if (is.null(r)) r <- character(0)
    if (is.null(cur)) cur <- character(0)
    for (k in union(names(r), names(cur))) {
      a <- if (k %in% names(r)) r[[k]] else NA_character_
      b <- if (k %in% names(cur)) cur[[k]] else NA_character_
      if (!identical(a, b)) {
        shown <- if (group == "settings") readable(a, b) else c(a, b)
        out <- c(out, sprintf("%s %s: recorded %s, now %s", sub("s$", "", group), k, short(shown[1L]), short(shown[2L])))
      }
    }
  }
  out
}

#' Begin a unit of work
#'
#' Takes the lock of the unit and compares the dependencies of the current run
#' with the completion record in the folder. The result says whether the unit has
#' to run.
#'
#' @inheritParams glow_unit_record
#' @param dependencies The dependencies of the current run: a list of the four
#'   named character vectors `settings`, `inputs`, `files` and `versions`.
#' @param force Logical. `TRUE` removes an existing record and runs the unit
#'   whatever the record states.
#' @param level `"analysis"` or `"dataset"`. It decides which way forward the
#'   refusal names: a new analysis name, or a new dataset name for an asset that
#'   belongs to the dataset.
#'
#' @return A handle, a list with `status` and what [glow_unit_complete()] needs.
#'   `status` is `"run"` (the unit has to run, and its lock is held) or
#'   `"complete"` (the record states the same dependencies and its outputs exist,
#'   and no lock is held).
#'
#' @details
#' The overwrite guard: when a record exists and states other settings, outside
#' files or versions, the function stops, prints each difference with the recorded
#' and the current value, and gives the two ways forward. It does so whether or
#' not the recorded outputs still exist. A record that differs in its inputs only
#' is out of date, because a part it read was produced again with the same
#' settings: the unit runs again, and the handle carries the differences in
#' `stale`. A record whose dependencies are the same and whose outputs are partly
#' gone makes the unit run again.
#'
#' A unit that runs again first loses its record and then the files that the
#' record lists inside the folder. So the folder never holds a record that
#' describes other results than its files, and a file that the new run does not
#' write (an optional table that was switched off) does not stay behind.
#'
#' The lock is a folder in `records/`, `<unit>.lock.<n>`, that holds a file `owner`
#' with the host, the process, the SLURM job and a token of the owner. Each
#' acquisition creates the next generation `n`, complete, and moves it into place
#' in one step, so two runs cannot hold the lock. A released generation keeps its
#' folder, with a file `released`, so the numbers only grow. A held lock is taken
#' over only when its owner is known to be gone: the process does not exist on
#' this host (or its number was given to another process), or the SLURM job is not
#' in the queue. A scheduler that does not answer proves nothing, and the lock
#' stays. Of several runs that find the same dead owner, one creates the next
#' generation. Any other held lock stops the function with the owner named.
#' Assumption: moving a folder to a name that is free is atomic, and fails when
#' the name is taken by a folder that is not empty, which holds for POSIX file
#' systems. Computational complexity: constant, plus the number of recorded files
#' of a unit that runs again.
#'
#' @examples
#' \dontrun{
#' h <- glow_unit_begin(folder, "scan", "chr21", dependencies)
#' if (h$status == "run") {
#'   # ... the work ...
#'   glow_unit_complete(h, outputs = file.path(folder, "results", "chr21.csv"))
#' }
#' }
#'
#' @seealso [glow_unit_complete()], [glow_unit_release()]
#' @export
glow_unit_begin <- function(folder, part, unit = NULL, dependencies, force = FALSE,
                            level = c("analysis", "dataset")) {
  level <- match.arg(level)
  stopifnot(is.character(folder), length(folder) == 1L, is.list(dependencies))
  dir.create(file.path(folder, "records"), recursive = TRUE, showWarnings = FALSE)
  lock <- .rr_take_lock(folder, part, unit)
  handle <- list(status = "run", folder = folder, part = part, unit = unit, lock = lock$dir, token = lock$token,
                 dependencies = dependencies[c("settings", "inputs", "files", "versions")],
                 record_file = .rr_record_file(folder, part, unit))
  ok <- FALSE
  on.exit(if (!ok) .rr_release_lock(lock$dir, lock$token), add = TRUE)
  state <- glow_unit_state(folder, part, unit, dependencies)
  if (!is.null(state$record)) {
    if (!isTRUE(force)) {
      if (state$complete) {
        .rr_release_lock(lock$dir, lock$token)
        handle$status <- "complete"; handle$lock <- NULL; handle$token <- NULL
        ok <- TRUE
        return(handle)
      }
      if (length(state$differences) && !state$stale) {
        what <- .rr_stem(part, unit)
        stop("The folder ", folder, " holds finished results of `", what, "` that were made with other ",
             "dependencies than this run:\n", paste0("  - ", state$differences, collapse = "\n"), "\n",
             if (isTRUE(state$record$seeded)) "They are files that were supplied with the dataset, not produced by a run. " else "",
             "Two ways forward: give this run ", if (level == "dataset") "a new `dataset` name" else "a new `analysis` name",
             ", which keeps the finished results, or run again with --force, which replaces them.", call. = FALSE)
      }
      # an input was produced again: the results are out of date
      if (state$stale) handle$stale <- state$differences
    }
    # The unit runs again. Its record goes first and its recorded files after it.
    unlink(handle$record_file, expand = FALSE)
    .rr_retire_outputs(folder, state$record$outputs)
  }
  ok <- TRUE
  handle
}

#' Complete a unit of work
#'
#' Writes the completion record of a unit and releases its lock. The record is
#' written last and in one move, so a unit that stopped half way has no record.
#'
#' @param handle The handle that [glow_unit_begin()] returned with the status
#'   `"run"`.
#' @param outputs Character vector. The files and folders the unit produced. Each
#'   has to exist. A folder stands for every file in it at this moment, and the
#'   record lists those files, so a file that is deleted later is noticed. Paths
#'   inside the unit's folder are recorded relative to it.
#' @param inputs_now `NULL`, or the `inputs` dependencies read again at this
#'   moment. When they differ from those at the begin, an input was replaced while
#'   the unit ran. The function then stops and writes no record.
#' @param files_now `NULL`, or the `files` dependencies read again at this moment,
#'   compared in the same way.
#' @param producer `NULL`, or a named character vector that says which process
#'   produced the unit, such as `c(slurm_job = "1234567")`. It is written into the
#'   record as it is and is no dependency: a caller that started a job for the
#'   unit reads it to tell whether the record is the one that job wrote.
#'
#' @return The record that was written, invisibly.
#'
#' @details The record gets a generation mark that is new at every run. A part
#'   that reads this unit's output stores the mark among its own dependencies, so
#'   it finds out when the output is replaced. The two checks at the end catch an
#'   input that differs between the begin and the end of the unit. They cannot
#'   catch one that was replaced and put back in between. The function also stops
#'   when the lock is no longer held under the handle's token.
#'
#' @seealso [glow_unit_begin()]
#' @export
glow_unit_complete <- function(handle, outputs = character(0), inputs_now = NULL, files_now = NULL, producer = NULL) {
  stopifnot(is.list(handle), identical(handle$status, "run"), is.character(outputs))
  on.exit(.rr_release_lock(handle$lock, handle$token), add = TRUE)
  if (!is.null(producer) && !(is.character(producer) && length(producer) > 0L && !anyNA(producer) &&
                              !is.null(names(producer)) && all(nzchar(names(producer)))))
    stop("glow_unit_complete(): `producer` is NULL or a named character vector.", call. = FALSE)
  what <- .rr_stem(handle$part, handle$unit)
  now <- list(inputs = inputs_now, files = files_now)
  diffs <- character(0)
  for (group in names(now)) if (!is.null(now[[group]]))
    diffs <- c(diffs, glow_dependency_differences(stats::setNames(list(handle$dependencies[[group]]), group),
                                                  stats::setNames(list(now[[group]]), group)))
  if (length(diffs))
    stop("An input of `", what, "` was replaced while it ran, so its results are not recorded as complete:\n",
         paste0("  - ", diffs, collapse = "\n"), "\nRun it again.", call. = FALSE)
  gone <- outputs[!file.exists(outputs)]
  if (length(gone)) stop("`", what, "` did not produce: ", paste(gone, collapse = ", "), call. = FALSE)
  if (!.rr_lock_is_mine(handle$lock, handle$token))
    stop("The lock of `", what, "` in ", handle$folder, " is no longer held by this run, so its results are not ",
         "recorded as complete. Run it again.", call. = FALSE)
  files <- unlist(lapply(outputs, function(p) if (dir.exists(p))
    list.files(p, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE) else p))
  record <- list(part = handle$part, unit = handle$unit, generation = .rr_generation(),
                 completed = format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"),
                 seeded = FALSE, outputs = .rr_relative(handle$folder, files), dependencies = handle$dependencies,
                 producer = .rr_named_text(as.list(producer)))
  .rr_write_record(record, handle$record_file)
  invisible(record)
}

#' Release the lock of a unit that did not finish
#'
#' @param handle A handle of [glow_unit_begin()].
#' @return `NULL`, invisibly. No record is written. A lock that another run took
#'   over in the meantime is left alone.
#' @seealso [glow_unit_complete()]
#' @export
glow_unit_release <- function(handle) {
  if (is.list(handle) && !is.null(handle$lock)) .rr_release_lock(handle$lock, handle$token)
  invisible(NULL)
}

#' Hold a named lock in a folder
#'
#' The SLURM mode of the `glow` command reads the queue and submits jobs. Two calls
#' that did this at the same time would each find a unit missing from the queue and
#' submit it. So a call holds a lock while it reads the queue and submits. The lock
#' is the lock of the units (see [glow_unit_begin()]): its owner is recorded, and a
#' lock whose owner is known to be gone is taken over.
#'
#' @param folder Character. The folder whose `records/` holds the lock. It is
#'   created when it does not exist.
#' @param name Character. The name of the lock.
#' @param wait Numeric. For how many seconds to wait for a lock that another
#'   process holds, looking again every two seconds. 0 stops at once.
#'
#' @return A handle for [glow_unit_release()].
#'
#' @details
#' Computational complexity: constant per look. Assumption: the same as for the
#' lock of a unit, that moving a folder to a free name is atomic.
#'
#' @examples
#' folder <- tempfile("dataset-")
#' reservation <- glow_lock(folder, "slurm-submission")
#' glow_unit_release(reservation)
#'
#' @seealso [glow_unit_begin()], [glow_unit_release()]
#' @export
glow_lock <- function(folder, name, wait = 0) {
  stopifnot(is.character(folder), length(folder) == 1L, is.character(name), length(name) == 1L,
            is.numeric(wait), length(wait) == 1L, wait >= 0)
  dir.create(file.path(folder, "records"), recursive = TRUE, showWarnings = FALSE)
  until <- Sys.time() + wait
  repeat {
    lock <- .rr_take_lock(folder, name, NULL, stop_if_held = Sys.time() >= until)
    if (!is.null(lock)) return(list(folder = folder, part = name, lock = lock$dir, token = lock$token))
    Sys.sleep(2)
  }
}

#' Supply the files of a unit instead of running it
#'
#' An example dataset ships files that a full run would produce, so that a new
#' user can start from a later step. This function puts the supplied files of one
#' unit in place and writes its seed record, which states the settings the files
#' were built with. A run whose settings are the same reads them as finished
#' results, and a run with other settings is told to produce its own.
#'
#' @inheritParams glow_unit_record
#' @param from Character vector. The supplied files, which have to exist.
#' @param to Character vector of the same length. Where each file belongs.
#' @param settings Named character vector. The settings the files were built
#'   with, as [glow_setting_text()] gives them. They are fixed when the files are
#'   built and do not follow a later edit of a configuration.
#'
#' @return One string, invisibly: `"present"` (the files and the seed record are
#'   in place, nothing was written), `"seeded"` (files were copied or the record
#'   was written), `"kept"` (a run produced this unit, so nothing is supplied),
#'   `"busy"` (another run holds the unit) or `"foreign"` (a destination holds a
#'   file that is not the supplied one and no seed record explains it, so nothing
#'   is supplied and the unit stays without a record).
#'
#' @details
#' The function never writes over the results of a run: a unit with a record
#' that a run produced is left alone, also when one of its files is missing. It
#' takes the lock of the unit, so it does not write while the unit runs. A file
#' is copied under a temporary name and moved into place, so a copy that stopped
#' half way leaves no file under the final name. A destination that already holds
#' the supplied bytes is accepted. Under a seed record the supplied files are
#' restored when they are missing or differ, because a seed record stands for the
#' supplied bytes. The generation mark of a seed record is a checksum of the
#' supplied files, so it stays the same as long as they do. Computational
#' complexity: linear in the size of the supplied files, which are read for their
#' checksums at every call.
#'
#' @seealso [glow_unit_state()], [glow_dependency_differences()]
#' @export
glow_seed_unit <- function(folder, part, unit = NULL, from, to, settings) {
  stopifnot(is.character(from), length(from) >= 1L, all(file.exists(from)), is.character(to),
            length(to) == length(from), is.character(settings), !is.null(names(settings)) || !length(settings))
  o <- order(to, method = "radix"); from <- from[o]; to <- to[o]
  settings <- settings[order(names(settings), method = "radix")]
  sums <- unname(tools::md5sum(from))
  generation <- paste0("seed-", .rr_text_md5(sums))
  same <- function() { m <- unname(tools::md5sum(to)); !is.na(m) & m == sums }   # a missing file gives NA
  intact <- function(rec) !is.null(rec) && rec$seeded && identical(rec$generation, generation) &&
    identical(rec$dependencies$settings, settings)
  rec <- glow_unit_record(folder, part, unit)
  if (!is.null(rec) && !rec$seeded) return(invisible("kept"))
  if (intact(rec) && all(same())) return(invisible("present"))
  dir.create(file.path(folder, "records"), recursive = TRUE, showWarnings = FALSE)
  lock <- .rr_take_lock(folder, part, unit, stop_if_held = FALSE)
  if (is.null(lock)) return(invisible("busy"))
  on.exit(.rr_release_lock(lock$dir, lock$token), add = TRUE)
  rec <- glow_unit_record(folder, part, unit)        # read again, now that the unit is held
  if (!is.null(rec) && !rec$seeded) return(invisible("kept"))
  ok <- same()
  if (is.null(rec) && any(file.exists(to) & !ok)) return(invisible("foreign"))
  for (k in which(!ok)) {
    dir.create(dirname(to[k]), recursive = TRUE, showWarnings = FALSE)
    tmp <- paste0(to[k], ".seed-", Sys.getpid())
    if (!file.copy(from[k], tmp, overwrite = TRUE) || !file.rename(tmp, to[k])) {
      unlink(tmp, expand = FALSE)
      stop("Could not copy ", from[k], " to ", to[k], call. = FALSE)
    }
  }
  if (!intact(rec))
    .rr_write_record(list(part = part, unit = unit, generation = generation,
                          completed = format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"), seeded = TRUE,
                          outputs = .rr_relative(folder, to),
                          dependencies = list(settings = settings, inputs = character(0), files = character(0),
                                              versions = character(0))),
                     .rr_record_file(folder, part, unit))
  invisible("seeded")
}

#################### INTERNAL HELPER FUNCTIONS ####################

#' The name of a unit in file names and messages
#' @param part The part.
#' @param unit The unit, or NULL.
#' @return `part`, or `part.unit`.
#' @keywords internal
#' @noRd
.rr_stem <- function(part, unit) if (is.null(unit)) part else paste(part, unit, sep = ".")

# The record file and the lock folder of a unit.
.rr_record_file <- function(folder, part, unit) file.path(folder, "records", paste0(.rr_stem(part, unit), ".yml"))
# (the base name of the lock: its generations are the folders <base>.1, <base>.2, ...)
.rr_lock_dir    <- function(folder, part, unit) file.path(folder, "records", paste0(.rr_stem(part, unit), ".lock"))

#' A map of a record as a named character vector
#' @param x The map as `.gc_plain()` returns it: a list of strings, or NULL.
#' @return A named character vector, sorted by name. Empty for an empty map.
#' @keywords internal
#' @noRd
.rr_named_text <- function(x) {
  if (is.null(x) || length(x) == 0L) return(stats::setNames(character(0), character(0)))
  v <- vapply(x, function(e) paste(as.character(unlist(e)), collapse = " "), character(1))
  v[order(names(v), method = "radix")]
}

#' The lines of a record as YAML
#'
#' @param x A named list whose values are strings, character vectors (written as
#'   a list) or named character vectors and lists (written as a mapping).
#' @param indent The indentation of the lines.
#' @return Character vector of lines. Every scalar and every key is in double
#'   quotes, so no text is read as a number, a logical or a null.
#' @keywords internal
#' @noRd
.rr_yaml_lines <- function(x, indent = "") {
  out <- character(0)
  for (k in names(x)) {
    v <- x[[k]]
    key <- paste0(indent, .gc_scalar_text(k, force_quote = TRUE), ":")
    if (is.null(v)) next
    if (is.list(v) || !is.null(names(v))) {
      v <- as.list(v)
      if (!length(v)) out <- c(out, paste(key, "{}"))
      else out <- c(out, key, .rr_yaml_lines(v, paste0(indent, "  ")))
    } else if (length(v) == 1L) {
      out <- c(out, paste(key, .gc_scalar_text(as.character(v), force_quote = TRUE)))
    } else if (!length(v)) {
      out <- c(out, paste(key, "[]"))
    } else {
      out <- c(out, key, paste0(indent, "  - ", vapply(as.character(v), .gc_scalar_text, character(1), force_quote = TRUE)))
    }
  }
  out
}

#' Write a completion record in one move
#'
#' @param record The record (see `glow_unit_record()`).
#' @param file The record file.
#' @return The file, invisibly. The text is written to a temporary file in the
#'   same folder and renamed, so a reader never sees half a record.
#' @keywords internal
#' @noRd
.rr_write_record <- function(record, file) {
  deps <- record$dependencies
  sorted <- function(v) if (length(v)) v[order(names(v), method = "radix")] else v
  body <- list(part = record$part, unit = record$unit, generation = record$generation,
               completed = record$completed, seeded = if (isTRUE(record$seeded)) "true" else NULL,
               outputs = as.character(record$outputs),
               dependencies = list(settings = sorted(deps$settings), inputs = sorted(deps$inputs),
                                   files = sorted(deps$files), versions = sorted(deps$versions)),
               producer = if (length(record$producer)) sorted(record$producer) else NULL)
  # `outputs` of length one has to stay a list, and an empty group an empty map
  lines <- c("# The completion record of one unit of work, written by GLOWpipeline when the unit had",
             "# finished. It states what the unit depended on and which files it produced. Do not edit.")
  for (k in names(body)) {
    v <- body[[k]]
    if (is.null(v)) next
    if (k == "outputs") {
      lines <- c(lines, if (length(v)) c("\"outputs\":", paste0("  - ", vapply(v, .gc_scalar_text, character(1), force_quote = TRUE)))
                        else "\"outputs\": []")
    } else if (k == "dependencies") {
      lines <- c(lines, "\"dependencies\":")
      for (g in names(v)) {
        gv <- v[[g]]
        lines <- c(lines, if (length(gv)) c(paste0("  \"", g, "\":"), .rr_yaml_lines(as.list(gv), "    "))
                          else paste0("  \"", g, "\": {}"))
      }
    } else lines <- c(lines, .rr_yaml_lines(stats::setNames(list(v), k)))
  }
  tmp <- tempfile(pattern = paste0(basename(file), "."), tmpdir = dirname(file), fileext = ".tmp")
  writeLines(enc2utf8(lines), tmp, useBytes = TRUE)
  if (!file.rename(tmp, file)) { unlink(tmp, expand = FALSE); stop("Could not write the record ", file, call. = FALSE) }
  invisible(file)
}

#' Take the lock of a unit
#'
#' @inheritParams glow_unit_record
#' @param stop_if_held Logical. `TRUE`: stop when another run holds the lock.
#'   `FALSE`: return `NULL` then.
#' @return `list(dir, token)`: the lock folder of this acquisition and the token
#'   that says the lock is this process's.
#' @details
#' The lock of a unit is a sequence of generations, the folders
#' `records/<unit>.lock.1`, `.lock.2` and so on. The highest generation is the
#' current one. It is held while it has a file `owner` and no file `released`.
#' Taking the lock is creating the next generation, whether the current one was
#' released or its owner is known to be gone: the folder is built beside its place,
#' with the `owner` file in it, and moved into place. The move fails when the
#' folder exists, so of any number of runs one creates a given generation, and a
#' generation is never seen without its owner. Generations only grow and a folder
#' is never replaced, so a run that judged an earlier owner cannot act on a later
#' one: it either loses the move, or it finds a higher generation right after it
#' and gives its own up. No step leaves a state that blocks the next run, so a run
#' that dies at any point is recovered from. The generations below the new one are
#' deleted by the new owner.
#' @keywords internal
#' @noRd
.rr_take_lock <- function(folder, part, unit, stop_if_held = TRUE) {
  base <- .rr_lock_dir(folder, part, unit)
  token <- .rr_generation()
  pid <- Sys.getpid()
  me <- c(host = unname(Sys.info()[["nodename"]]), pid = as.character(pid), started = .rr_process_start(pid),
          slurm_job = Sys.getenv("SLURM_JOB_ID", ""), since = format(Sys.time(), "%Y-%m-%d %H:%M:%S UTC", tz = "UTC"),
          token = token)
  state <- NULL
  for (attempt in 1:5) {
    state <- .rr_lock_state(base)
    if (state$held && !isTRUE(.rr_lock_is_stale(state$owner))) break       # a live owner, or one that cannot be judged
    mine <- paste0(base, ".", state$n + 1L)
    stage <- paste0(base, ".new-", token)
    dir.create(stage, showWarnings = FALSE)
    writeLines(paste(names(me), me, sep = ": "), file.path(stage, "owner"))
    won <- suppressWarnings(file.rename(stage, mine))
    if (!won) { unlink(stage, recursive = TRUE, expand = FALSE); next }                    # another run created this generation: look again
    newer <- .rr_lock_generations(base)
    if (max(newer) > state$n + 1L) { unlink(mine, recursive = TRUE, expand = FALSE); next }  # the view was out of date: this generation is void
    for (n in newer[newer <= state$n]) unlink(paste0(base, ".", n), recursive = TRUE, expand = FALSE)
    return(list(dir = mine, token = token))
  }
  if (!stop_if_held) return(NULL)
  owner <- state$owner
  stop("Another run holds `", .rr_stem(part, unit), "` in ", folder, " (",
       if (length(owner)) paste(owner[!startsWith(owner, "token: ")], collapse = "; ") else "its owner is not recorded",
       "). Wait for it to finish. If that run is gone, delete the folder ", paste0(base, ".", state$n),
       " and run again.", call. = FALSE)
}

#' Release a lock that is still this process's
#'
#' @param lock The lock folder of the acquisition.
#' @param token The token of the handle.
#' @return `TRUE` when the lock was released, `FALSE` when it was not this
#'   process's any more (or did not exist). The folder stays, with a file
#'   `released` in place of `owner`: it keeps the number of the generation, which
#'   the next acquisition counts on from.
#' @keywords internal
#' @noRd
.rr_release_lock <- function(lock, token) {
  if (is.null(lock) || is.null(token) || !dir.exists(lock)) return(invisible(FALSE))
  if (!identical(.rr_owner_field(lock, "token"), token)) return(invisible(FALSE))
  file.create(file.path(lock, "released"))
  unlink(file.path(lock, "owner"), expand = FALSE)
  invisible(TRUE)
}

# The numbers of the generations of a lock that exist, in rising order.
.rr_lock_generations <- function(base) {
  prefix <- paste0(basename(base), ".")
  found <- list.files(dirname(base))
  rest <- substring(found[startsWith(found, prefix)], nchar(prefix) + 1L)
  sort(as.integer(rest[grepl("^[0-9]+$", rest)]))
}

# The current state of a lock: the highest generation `n` (0 when there is none), its
# folder, whether it is held, and the lines of its owner.
#
# Other runs release and take the lock while this one looks. A generation that shows
# neither an owner nor a release was caught between two looks: it was released in the
# meantime (the look is taken again), or a newer generation took its place and its
# folder is gone (the generations are listed again). A generation that still stands
# without an owner and without a release is one that nobody can judge, and it holds.
.rr_lock_state <- function(base) {
  for (look in 1:10) {
    gens <- .rr_lock_generations(base)
    if (!length(gens)) return(list(n = 0L, dir = NULL, held = FALSE, owner = character(0)))
    n <- max(gens); dir <- paste0(base, ".", n)
    released <- file.exists(file.path(dir, "released"))
    owner <- if (released) character(0) else .rr_read_owner(dir)
    if (!released && !length(owner)) {
      if (file.exists(file.path(dir, "released"))) released <- TRUE
      else if (!dir.exists(dir)) next
    }
    break
  }
  list(n = n, dir = dir, held = !released, owner = owner)
}

# Whether a lock folder is still the current generation and still held under a token.
.rr_lock_is_mine <- function(lock, token) {
  if (is.null(lock) || is.null(token) || !dir.exists(lock)) return(FALSE)
  base <- sub("\\.[0-9]+$", "", lock)
  identical(.rr_owner_field(lock, "token"), token) && !file.exists(file.path(lock, "released")) &&
    identical(max(.rr_lock_generations(base)), as.integer(sub("^.*\\.", "", lock)))
}

# The lines of the `owner` file of a lock folder (empty when there is none), and one
# field of them. `.rr_owner_field()` takes the lines or the lock folder.
.rr_read_owner <- function(lock) {
  f <- file.path(lock, "owner")
  if (!file.exists(f)) return(character(0))
  # the file can go between the two lines, when its owner releases the lock
  tryCatch(suppressWarnings(readLines(f, warn = FALSE)), error = function(e) character(0))
}
.rr_owner_field <- function(owner, key) {
  if (length(owner) == 1L && dir.exists(owner)) owner <- .rr_read_owner(owner)
  v <- sub(paste0("^", key, ": "), "", owner[startsWith(owner, paste0(key, ": "))])
  if (length(v)) v[1L] else ""
}

#' What tells two processes of the same number apart
#'
#' @param pid A process number.
#' @return The start time of the process as the system counts it (the 22nd field
#'   of `/proc/<pid>/stat`), or an empty string where that file does not exist.
#' @keywords internal
#' @noRd
.rr_process_start <- function(pid) {
  f <- file.path("/proc", pid, "stat")
  if (!file.exists(f)) return("")
  line <- tryCatch(readLines(f, n = 1L, warn = FALSE), error = function(e) "")
  # the command name, in parentheses, may hold spaces: count from behind it
  rest <- strsplit(sub("^.*\\) ", "", line), " ", fixed = TRUE)[[1]]
  if (length(rest) >= 20L) rest[20L] else ""
}

#' Whether the owner of a lock is gone
#'
#' @param owner The lines of the lock's `owner` file.
#' @return `TRUE` when the owner is known to be gone, `FALSE` when it is known to
#'   be there, `NA` when that cannot be told.
#' @details On the owner's host the process tells: it does not exist, or the
#'   process of that number started at another time. This comes first, because a
#'   session that is a SLURM job outlives the processes started in it. For an
#'   owner on another host, its SLURM job tells: `squeue` answers and does not list
#'   the job, or says that the job number is not valid any more. A scheduler that
#'   fails for another reason gives `NA`.
#' @keywords internal
#' @noRd
.rr_lock_is_stale <- function(owner) {
  if (!length(owner)) return(NA)
  host <- .rr_owner_field(owner, "host"); job <- .rr_owner_field(owner, "slurm_job")
  pid <- suppressWarnings(as.integer(.rr_owner_field(owner, "pid")))
  if (identical(host, unname(Sys.info()[["nodename"]])) && !is.na(pid) && dir.exists("/proc")) {
    if (!dir.exists(file.path("/proc", pid))) return(TRUE)
    started <- .rr_owner_field(owner, "started"); now <- .rr_process_start(pid)
    return(nzchar(started) && nzchar(now) && !identical(started, now))
  }
  if (nzchar(job) && nzchar(Sys.which("squeue"))) {
    answer <- tryCatch(suppressWarnings(system2("squeue", c("-h", "-j", job, "-o", "%A"), stdout = TRUE, stderr = TRUE)),
                       error = function(e) NULL)
    if (is.null(answer)) return(NA)
    status <- attr(answer, "status")
    if (is.null(status) || status == 0L) return(!any(trimws(answer) == job))
    if (any(grepl("Invalid job id", answer, fixed = TRUE))) return(TRUE)
  }
  NA
}

# The checksum of some lines of text.
.rr_text_md5 <- function(lines) {
  tmp <- tempfile(); on.exit(unlink(tmp, expand = FALSE), add = TRUE)
  writeLines(lines, tmp)
  unname(tools::md5sum(tmp))
}

# Paths as a record lists them: relative to the unit's folder when inside it.
.rr_relative <- function(folder, paths) {
  if (!length(paths)) return(character(0))
  prefix <- paste0(normalizePath(folder), "/")
  full <- normalizePath(paths)
  sort(unique(ifelse(startsWith(full, prefix), substring(full, nchar(prefix) + 1L), full)), method = "radix")
}

# Delete the files that a record lists inside its folder, and nothing else.
#
# A recorded name is taken literally, so a name with a `*` in it deletes that file and
# not what the pattern would match. A file is deleted only when the folder it lies in
# is, in its real form, the unit's folder or a folder inside it: a recorded folder that
# was replaced by a link to another place is not followed. A link is removed as a link.
# Folders, and files outside the unit's folder, are left: a run writes over the latter.
.rr_retire_outputs <- function(folder, outputs) {
  root <- normalizePath(folder)
  for (rel in outputs[!grepl("^/", outputs)]) {
    path <- file.path(folder, rel)
    is_link <- nzchar(Sys.readlink(path))
    if (!is_link && (!file.exists(path) || dir.exists(path))) next
    parent <- normalizePath(dirname(path), mustWork = FALSE)
    if (!identical(parent, root) && !startsWith(parent, paste0(root, "/"))) next
    unlink(path, expand = FALSE)
  }
  invisible(NULL)
}

#' A mark that is new at every run
#'
#' @return One string: the time to the millisecond, the process, and a random
#'   tail. The tail comes from `tempfile()`, which does not touch the state of R's
#'   random number generator, so a run with a seed stays reproducible.
#' @keywords internal
#' @noRd
.rr_generation <- function() {
  sprintf("%s-%d-%s", format(Sys.time(), "%Y%m%dT%H%M%OS3", tz = "UTC"), Sys.getpid(),
          sub("^file", "", basename(tempfile())))
}
