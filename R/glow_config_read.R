# ==============================================================================
# File: glow_config_read.R
# Purpose: Read a GLOW configuration: follow the chain of YAML files from a user
#   file to the base, merge them, and check every value against the settings
#   registry.
#
# EXPORTED FUNCTIONS:
#   - read_glow_config()    A user file (or none) -> the resolved, checked settings
#
# INTERNAL HELPERS:
#   - .gc_chain()           The files of the chain, base first, each parsed
#   - .gc_check_names()     Unknown sections and settings, with the nearest name
#   - .gc_check_values()    The values one file sets, each by its registered type
#   - .gc_merge()           One file's settings over the settings so far
#   - .gc_check_study()     The names in `dataset` and `analysis`
#
# The reading of a single file, the registry and the conversion by type are in

#################### EXPORTED MAIN FUNCTIONS ####################

#' Read and check a GLOW configuration
#'
#' Resolves a YAML configuration against the base configuration and the settings
#' registry. A user file holds only the settings that differ from the base. Its
#' `extends` line names the file it builds on, and a file without that line builds
#' on the base, so configurations form a chain that ends at the base.
#'
#' @param config Path of the user's YAML configuration, or `NULL` to resolve the
#'   base alone.
#' @param base Path of the base configuration (`base_config.yml`), which sets every
#'   registered setting.
#' @param registry Path of the settings registry (`settings-registry.yml`), which
#'   gives every setting its section, type and default.
#' @param root Directory from which a relative data path of the configuration is
#'   meant. It is recorded with the result and not used to read the configuration
#'   files themselves. Default: the working directory.
#'
#' @return An object of class `glow_config`: a list with the top-level settings
#'   (such as `dataset`) and one list per section (such as `snv_set`), every value
#'   converted to its registered type. A setting that is null is present with the
#'   value `NULL`. A setting of type `r_function` is an R function whose source
#'   text is in the attribute `glow_source`. The attributes `glow_chain` (the files
#'   read, base first, with their checksums), `glow_registry` (the registry as
#'   read, with the type, the `derived` mark and the old name of every setting) and
#'   `glow_root` record where the settings came from.
#'
#' @details
#' Algorithm:
#' \enumerate{
#'   \item Read the registry and check that every default is valid under its type.
#'   \item Follow the chain: the user file, the file its `extends` names (read
#'     relative to the file that contains the line), and so on to the base. Files
#'     are compared by their real paths, so a cycle is found. A file may name the
#'     base by its bare file name from any folder.
#'   \item Check each file on its own: a section or a setting that the registry
#'     does not hold is an error, with the nearest registered name, and so is a
#'     value that cannot be read as the type of its setting, even when a later file
#'     of the chain replaces it.
#'   \item Merge from the base toward the user file. Sections merge setting by
#'     setting, and a setting that is a map with declared fields merges field by
#'     field, unless the registry marks it `replace_whole`. Every other value (a
#'     scalar, a list, a map with free keys, a value that takes one of several
#'     forms) replaces the inherited value whole. An explicit null replaces the
#'     inherited value. A section with nothing after its colon sets nothing.
#'   \item Convert each value by its registered type.
#' }
#' Every scalar is read as the text the user wrote and converted by the registered
#' type of its setting, so `no` stays text for a text setting, `1e6` is a number
#' for a numeric one, and `010` is ten. A logical takes `true` or `false`. A list
#' of integers may be written as a range, `1-22`.
#'
#' Problems of names and types are collected and reported in one error, each with
#' the file it is in. A file that cannot be found or parsed, and a cycle, stop the
#' function at once.
#'
#' Assumptions: the function checks the structure and the types of a configuration
#' and nothing else. It does not check that a path exists, and it does not fill a
#' null that a later step derives. A setting of type `r_function` is parsed and
#' never run. The function it returns has the global environment, so, when a step
#' calls it later, it sees base R and the attached packages and nothing of the
#' configuration; a function that depends on other objects of the session is not
#' reproducible from the configuration alone.
#'
#' Computational complexity: linear in the number of settings and files.
#'
#' @examples
#' \dontrun{
#' cfg <- read_glow_config("cohort0203.yml",
#'                         base = "base-config/base_config.yml",
#'                         registry = "base-config/settings-registry.yml")
#' cfg$dataset
#' cfg$snv_set$ld_threshold
#' }
#'
#' @seealso \code{\link{write_resolved_config}}, \code{\link{write_base_config}}
#' @export
read_glow_config <- function(config = NULL, base, registry, root = ".") {
  stopifnot(is.null(config) || (is.character(config) && length(config) == 1L),
            is.character(base), length(base) == 1L,
            is.character(registry), length(registry) == 1L,
            is.character(root), length(root) == 1L)
  if (!dir.exists(root)) stop("`root` is not a directory: ", root, call. = FALSE)
  reg   <- .gc_read_registry(registry)
  chain <- .gc_chain(config, base)

  acc <- new.env(); acc$problems <- character()
  merged <- .gc_empty_map()
  for (link in chain) {
    tree <- link$tree
    tree$extends <- NULL        # the chain is already followed
    tree$provenance <- NULL     # the reserved key of a resolved configuration
    .gc_check_names(tree, reg, link$file, acc)
    .gc_check_values(tree, reg, acc)
    merged <- .gc_merge(merged, tree, reg)
  }

  # Convert every registered setting. The base sets all of them, so a setting that
  # is absent after the merge means a base that is older than the registry.
  cfg <- list()
  for (e in reg$entries) {
    node <- if (e$section == "top") merged[[e$name]] else if (.gc_is_map(merged[[e$section]])) merged[[e$section]][[e$name]]
    if (is.null(node)) {
      acc$problems <- c(acc$problems, paste0(
        "`", e$key, "` is set by no file of the chain. The base configuration ", base,
        " has to set every registered setting; regenerate it from the registry."))
      next
    }
    value <- .gc_convert(node, e$spec, e$key, acc)
    if (e$section == "top") cfg[e$name] <- list(value)
    else {
      if (is.null(cfg[[e$section]])) cfg[[e$section]] <- list()
      cfg[[e$section]][e$name] <- list(value)
    }
  }
  .gc_check_study(cfg, merged, acc)

  # the per-file check and the conversion of the merged values report the same problem twice
  problems <- unique(acc$problems)
  if (length(problems))
    stop("The configuration has ", length(problems), " problem(s):\n",
         paste0("  - ", problems, collapse = "\n"), call. = FALSE)

  structure(cfg, class = "glow_config",
            glow_chain = data.frame(file = vapply(chain, `[[`, character(1), "file"),
                                    md5 = vapply(chain, `[[`, character(1), "md5"),
                                    stringsAsFactors = FALSE),
            glow_registry = reg,
            glow_root = normalizePath(root))
}

#################### INTERNAL HELPER FUNCTIONS ####################

#' The files of a configuration chain, base first
#'
#' @param config The user file, or NULL.
#' @param base The base configuration.
#' @return A list, base first and the user file last, of `list(file, md5, tree)`.
#' @details `extends` is read from the parsed mapping, wherever it stands in the
#'   file, and its path is relative to the file that contains it. A file may also
#'   name the base by its bare file name, wherever the base is installed. A path
#'   with a folder in it that does not exist is an error, even when its last part
#'   is the base's file name. Files are compared by real path.
#' @keywords internal
#' @noRd
.gc_chain <- function(config, base) {
  if (!file.exists(base)) stop("Base configuration not found: ", base, call. = FALSE)
  base_real <- normalizePath(base)
  links <- list(); seen <- character()
  cur <- config
  while (!is.null(cur)) {
    if (!file.exists(cur))
      stop("Configuration file not found: ", cur,
           if (length(seen)) paste0("\n  named by `extends` in the chain: ", paste(seen, collapse = " -> ")) else "",
           call. = FALSE)
    real <- normalizePath(cur)
    if (identical(real, base_real)) break                     # the chain reached the base
    if (real %in% seen)
      stop("The configuration files extend each other in a cycle: ",
           paste(c(seen, real), collapse = " -> "), call. = FALSE)
    seen <- c(seen, real)
    tree <- .gc_parse_file(real, top = "map")
    links <- c(list(list(file = real, md5 = attr(tree, "glow_md5"), tree = tree)), links)   # of the bytes parsed
    ext <- tree$extends
    if (is.null(ext) || .gc_is_null(ext)) break               # no `extends`: builds on the base
    if (!.gc_is_scalar(ext) || !nzchar(trimws(ext$text)))
      stop(real, ": `extends` has to name one configuration file.", call. = FALSE)
    nxt <- trimws(ext$text)
    is_abs <- grepl("^(/|~|[A-Za-z]:[/\\\\])", nxt)
    cand <- if (is_abs) path.expand(nxt) else file.path(dirname(real), nxt)
    # the bare file name of the base means the base, wherever it is installed
    bare_base <- !is_abs && identical(nxt, basename(base_real))
    cur <- if (!file.exists(cand) && bare_base) base_real else cand
  }
  base_tree <- .gc_parse_file(base_real, top = "map")
  if (!is.null(base_tree$extends) && !.gc_is_null(base_tree$extends))
    stop(base_real, ": the base configuration cannot extend another file.", call. = FALSE)
  c(list(list(file = base_real, md5 = attr(base_tree, "glow_md5"), tree = base_tree)), links)
}

#' Check the names of one file against the registry
#'
#' Records a problem for a top-level name that is neither a top-level setting nor a
#' section, for a section that is not a mapping, and for a setting that its section
#' does not hold. The message names the nearest registered name, and the section
#' that does hold the name when another one does. A section with nothing after its
#' colon is accepted as one that sets nothing, which is what a user file looks like
#' when every setting of a section is commented out.
#' @param tree One file's parsed tree.
#' @param reg The registry.
#' @param file The file, for the messages.
#' @param acc The problem collector.
#' @return `NULL`, invisibly.
#' @keywords internal
#' @noRd
.gc_check_names <- function(tree, reg, file, acc) {
  sec_of  <- vapply(reg$entries, `[[`, character(1), "section")
  name_of <- vapply(reg$entries, `[[`, character(1), "name")
  top_names <- name_of[sec_of == "top"]
  sections  <- setdiff(reg$sections, "top")
  add <- function(msg) acc$problems <- c(acc$problems, paste0(file, ": ", sub("?.", "?", msg, fixed = TRUE)))
  elsewhere <- function(n) {
    hit <- unique(sec_of[name_of == n & sec_of != "top"])
    if (length(hit)) paste0("; it is a setting of the section `", paste(hit, collapse = "`, `"), "`") else ""
  }
  for (k in names(tree)) {
    if (k %in% top_names) next
    if (!k %in% sections) {
      hint <- elsewhere(k)
      if (!nzchar(hint)) hint <- .gc_nearest(k, c(top_names, sections))
      add(paste0("`", k, "` is not a setting or a section", hint, "."))
      next
    }
    sec <- tree[[k]]
    if (.gc_is_null(sec)) next                                             # a section that sets nothing
    if (!.gc_is_map(sec)) { add(paste0("the section `", k, "` has to be a mapping of settings.")); next }
    known <- name_of[sec_of == k]
    for (n in setdiff(names(sec), known)) {
      hint <- .gc_nearest(n, known)
      if (!nzchar(hint)) {
        hit <- unique(sec_of[name_of == n])
        if (length(hit)) hint <- paste0("; it is a setting of ",
          paste(ifelse(hit == "top", "the top level", paste0("the section `", hit, "`")), collapse = " and "))
      }
      add(paste0("`", k, ".", n, "` is not a registered setting", hint, "."))
    }
  }
  invisible(NULL)
}

#' Check the values one file sets, each by its registered type
#'
#' Every file of a chain is held to the registry on its own, so that a value that
#' cannot be read is reported in the file that holds it, even when a file later in
#' the chain replaces it. A map with declared fields may be partial here: the
#' fields a file does not name come from the files before it.
#' @param tree One file's parsed tree.
#' @param reg The registry.
#' @param acc The problem collector.
#' @return `NULL`, invisibly.
#' @keywords internal
#' @noRd
.gc_check_values <- function(tree, reg, acc) {
  for (e in reg$entries) {
    node <- if (e$section == "top") tree[[e$name]] else if (.gc_is_map(tree[[e$section]])) tree[[e$section]][[e$name]]
    # a map that merges field by field may be partial in one file
    partial <- identical(e$spec$type, "map") && !is.null(e$spec$fields) && !isTRUE(e$replace_whole)
    if (!is.null(node)) .gc_convert(node, e$spec, e$key, acc, partial = partial)
  }
  invisible(NULL)
}

#' Merge one file's settings over the settings so far
#'
#' Works on the parsed trees, before any conversion, so an explicit null is still a
#' token and replaces the inherited value like any other value.
#' @param old The merged tree so far.
#' @param new One file's tree.
#' @param reg The registry, which says which settings are maps with declared fields
#'   and which of those are replaced whole.
#' @return The merged tree.
#' @keywords internal
#' @noRd
.gc_merge <- function(old, new, reg) {
  sections <- setdiff(reg$sections, "top")
  one <- function(old_value, new_value, e) {
    if (is.null(e) || isTRUE(e$replace_whole)) new_value else .gc_merge_value(old_value, new_value, e$spec)
  }
  for (k in names(new)) {
    if (!k %in% sections) {                                   # a top-level setting (an unknown name is reported apart)
      old[[k]] <- one(old[[k]], new[[k]], reg$entries[[k]])
      next
    }
    sec <- new[[k]]
    if (!.gc_is_map(sec)) next                                # a null or malformed section sets nothing
    if (is.null(old[[k]])) old[[k]] <- .gc_empty_map()
    for (n in names(sec)) old[[k]][[n]] <- one(old[[k]][[n]], sec[[n]], reg$entries[[paste0(k, ".", n)]])
  }
  old
}

#' Merge one value over the inherited one
#'
#' A map with declared fields merges field by field; everything else replaces.
#' @param old The inherited node.
#' @param new The node of the later file.
#' @param spec The type spec of the value.
#' @return The merged node.
#' @keywords internal
#' @noRd
.gc_merge_value <- function(old, new, spec) {
  if (identical(spec$type, "map") && !is.null(spec$fields) && .gc_is_map(old) && .gc_is_map(new)) {
    for (f in names(new)) {
      fs <- spec$fields[[f]]
      old[[f]] <- if (is.null(fs)) new[[f]] else .gc_merge_value(old[[f]], new[[f]], fs)
    }
    return(old)
  }
  new
}

#' Check the names of the study
#'
#' `dataset` and `analysis` become folder names, so each has to be one name of
#' letters, digits, `_`, `-` and `.`, other than `.` and `..`. `data` is the folder
#' of the dataset-level assets and cannot be an analysis.
#' @param cfg The converted settings.
#' @param merged The merged tree, whose nodes name the file each value came from.
#' @param acc The problem collector.
#' @return `NULL`, invisibly.
#' @keywords internal
#' @noRd
.gc_check_study <- function(cfg, merged, acc) {
  for (k in intersect(c("dataset", "analysis"), names(cfg))) {
    v <- cfg[[k]]
    if (!is.character(v) || length(v) != 1L) next             # a type problem is already recorded
    if (!grepl("^[A-Za-z0-9._-]+$", v) || v %in% c(".", ".."))
      .gc_problem(acc, k, merged[[k]], "is used as a folder name: letters, digits, `_`, `-` and `.` only")
    else if (k == "analysis" && v == "data")
      .gc_problem(acc, k, merged[[k]], "cannot be `data`, the folder of the dataset-level assets")
  }
  invisible(NULL)
}
