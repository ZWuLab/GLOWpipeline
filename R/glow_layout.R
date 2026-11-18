# ==============================================================================
# File: glow_layout.R
# Purpose: The derived run layout of a GLOW configuration: every output folder
#   follows from the three names `output_root`, `dataset` and `analysis`, so no
#   configuration declares an output folder.
#
# EXPORTED FUNCTIONS:
#   - glow_layout()            A resolved configuration -> the folders of its run
#   - glow_absolute_paths()    A resolved configuration with every path made absolute
#
# INTERNAL HELPERS:
#   - .gl_absolute()           One path against a root
#   - .gl_absolute_value()     The paths inside one value, by its registered type
#

#################### EXPORTED MAIN FUNCTIONS ####################

#' The folders of a run, derived from its configuration
#'
#' Turns a resolved configuration into the folders its steps read and write. The
#' layout has two levels. The assets that depend on the genotypes alone (the GDS
#' files, the annotated GDS files and the principal components) belong to the
#' dataset and sit under `data/`. Everything that depends on the outcome, the
#' covariates or an analysis setting sits under the analysis.
#'
#' ```
#' <output_root>/<dataset>/
#'   data/gds/                       the converted GDS files
#'   data/gds_favor/<match>/gds/     the annotated GDS files
#'   data/gds_favor/<match>/csv/     the annotation tables
#'   data/pcs/                       the principal components
#'   <analysis>/data_prep/           the phenotype bundle
#'   <analysis>/training/b/          the B training table and model
#'   <analysis>/training/pi/         the PI case tables and models
#'   <analysis>/single_variant/      the single-variant scan
#'   <analysis>/snv_set/<region_type>/   one folder per region type
#' ```
#'
#' @param cfg A `glow_config` object, as [read_glow_config()] returns it. It has
#'   to hold the top-level settings `output_root`, `dataset` and `analysis`.
#' @param match Character. The name of the level under `data/gds_favor/`, which
#'   is the method that matched the variants to the annotation database. Default:
#'   the setting `shared.favor_match_method`.
#' @param region_types Character vector. The region types that get a folder
#'   under `snv_set/`. Default: the setting `snv_set.region_type`.
#'
#' @return A list of absolute paths: `output_root`, `dataset_dir`, `analysis_dir`,
#'   the list `data` (`dir`, `gds`, `gds_favor`, `gds_favor_gds`, `gds_favor_csv`,
#'   `pcs`) and the list `steps` (`data_prep`, `training`, `training_b`,
#'   `training_pi`, `single_variant`, `snv_set_dir`, and `snv_set`, a character
#'   vector named by region type). The elements `dataset` and `analysis` repeat
#'   the two names.
#'
#' @details
#' A relative `output_root` is read from the root the configuration was loaded
#' with (the attribute `glow_root`). The function builds paths and nothing else:
#' it creates no folder and reads no file. Computational complexity: constant.
#'
#' Assumptions: `dataset` and `analysis` are single names without a path
#' separator, which [read_glow_config()] checks.
#'
#' @examples
#' \dontrun{
#' cfg <- read_glow_config("cohort0203.yml", base = "base-config/base_config.yml",
#'                         registry = "base-config/settings-registry.yml")
#' lay <- glow_layout(cfg)
#' lay$steps$snv_set[["gene"]]
#' }
#'
#' @seealso [read_glow_config()], [glow_absolute_paths()]
#' @export
glow_layout <- function(cfg, match = cfg$shared$favor_match_method,
                        region_types = cfg$snv_set$region_type) {
  if (!inherits(cfg, "glow_config"))
    stop("`cfg` has to be a configuration read by read_glow_config().", call. = FALSE)
  one_name <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  for (nm in c("output_root", "dataset", "analysis"))
    if (!one_name(cfg[[nm]]))
      stop("The configuration has no usable `", nm, "` (one non-empty text value).", call. = FALSE)
  if (!one_name(match))
    stop("`match` has to be one non-empty text value (the setting `shared.favor_match_method`).",
         call. = FALSE)
  if (is.null(region_types)) region_types <- character(0)
  if (!is.character(region_types) || anyNA(region_types) || any(!nzchar(region_types)) ||
      anyDuplicated(region_types))
    stop("`region_types` has to be a vector of distinct, non-empty names.", call. = FALSE)
  # A name becomes a folder, so it may not leave its parent.
  for (nm in c(cfg$dataset, cfg$analysis, match, region_types))
    if (grepl("[/\\\\]", nm) || nm %in% c(".", ".."))
      stop("`", nm, "` cannot name a folder of the layout.", call. = FALSE)

  root        <- .gl_absolute(cfg$output_root, attr(cfg, "glow_root"))
  dataset_dir <- file.path(root, cfg$dataset)
  data_dir    <- file.path(dataset_dir, "data")
  favor       <- file.path(data_dir, "gds_favor")
  analysis    <- file.path(dataset_dir, cfg$analysis)
  training    <- file.path(analysis, "training")
  snv_set_dir <- file.path(analysis, "snv_set")

  list(
    output_root  = root,
    dataset      = cfg$dataset,
    analysis     = cfg$analysis,
    dataset_dir  = dataset_dir,
    data = list(
      dir           = data_dir,
      gds           = file.path(data_dir, "gds"),
      gds_favor     = favor,
      gds_favor_gds = file.path(favor, match, "gds"),
      gds_favor_csv = file.path(favor, match, "csv"),
      pcs           = file.path(data_dir, "pcs")),
    analysis_dir = analysis,
    steps = list(
      data_prep      = file.path(analysis, "data_prep"),
      training       = training,
      training_b     = file.path(training, "b"),
      training_pi    = file.path(training, "pi"),
      single_variant = file.path(analysis, "single_variant"),
      snv_set_dir    = snv_set_dir,
      snv_set        = stats::setNames(file.path(snv_set_dir, region_types), region_types)))
}

#' A configuration with every path in absolute form
#'
#' Returns the configuration with each value of the registered type `path` made
#' absolute, so that the record of a run does not depend on the directory it is
#' read in. A relative path is read from the root the configuration was loaded
#' with. Paths inside a map, a list or one of several forms are found through the
#' registered type of the setting.
#'
#' @param cfg A `glow_config` object, as [read_glow_config()] returns it.
#'
#' @return The configuration, of the same class and with the same attributes. A
#'   path that was absolute, and a null, are unchanged.
#'
#' @details
#' A path is made absolute as text: the root is put in front of it, and `.` and
#' `..` are left as written. The file does not have to exist, and a symbolic link
#' is not followed. Computational complexity: linear in the number of settings.
#'
#' @examples
#' \dontrun{
#' cfg <- glow_absolute_paths(read_glow_config(NULL, base = "base-config/base_config.yml",
#'                                             registry = "base-config/settings-registry.yml"))
#' cfg$data$pheno_csv
#' }
#'
#' @seealso [read_glow_config()], [write_resolved_config()]
#' @export
glow_absolute_paths <- function(cfg) {
  if (!inherits(cfg, "glow_config"))
    stop("`cfg` has to be a configuration read by read_glow_config().", call. = FALSE)
  root <- attr(cfg, "glow_root")
  for (e in attr(cfg, "glow_registry")$entries) {
    if (e$section == "top") {
      cfg[e$name] <- list(.gl_absolute_value(cfg[[e$name]], e$spec, root))
    } else if (!is.null(cfg[[e$section]])) {
      cfg[[e$section]][e$name] <- list(.gl_absolute_value(cfg[[e$section]][[e$name]], e$spec, root))
    }
  }
  cfg
}

#################### INTERNAL HELPER FUNCTIONS ####################

#' One path against a root
#'
#' @param path Character. One path, absolute or relative.
#' @param root Character. The directory a relative path is read from.
#' @return The path in absolute form, as text. `~` is expanded.
#' @keywords internal
#' @noRd
.gl_absolute <- function(path, root) {
  path <- path.expand(path)
  windows <- identical(.Platform$OS.type, "windows")
  if (!windows && grepl("^[A-Za-z]:[/\\\\]", path))
    stop("The path `", path, "` starts with a drive letter, which is a Windows path. On this system write ",
         "an absolute path with a leading slash, or a relative path.", call. = FALSE)
  if (startsWith(path, "/") || (windows && grepl("^[A-Za-z]:[/\\\\]", path))) path else file.path(root, path)
}

#' The paths inside one value, by its registered type
#'
#' @param value The value of a setting, or of a part of one.
#' @param spec The registered type of that value (a list with `type` and, for a
#'   value with parts, `fields`, `values`, `items` or `one_of`).
#' @param root Character. The directory a relative path is read from.
#' @return The value with every part of type `path` made absolute.
#' @details A map with declared fields is walked field by field, a map with one
#'   type for all values and a sequence item by item. For a value that takes one
#'   of several forms, the form is the first one the value is valid for.
#' @keywords internal
#' @noRd
.gl_absolute_value <- function(value, spec, root) {
  if (is.null(value)) return(NULL)
  type <- spec$type
  if (identical(type, "path")) return(.gl_absolute(value, root))
  if (identical(type, "map") && is.list(value)) {
    if (!is.null(spec$fields)) {
      for (f in intersect(names(value), names(spec$fields)))
        value[f] <- list(.gl_absolute_value(value[[f]], spec$fields[[f]], root))
    } else if (!is.null(spec$values)) {
      for (i in seq_along(value)) value[i] <- list(.gl_absolute_value(value[[i]], spec$values, root))
    }
    return(value)
  }
  if (identical(type, "sequence") && is.list(value) && !is.null(spec$items)) {
    for (i in seq_along(value)) value[i] <- list(.gl_absolute_value(value[[i]], spec$items, root))
    return(value)
  }
  if (identical(type, "one_of")) {
    for (form in spec$one_of)
      if (!nzchar(.gc_valid(value, form))) return(.gl_absolute_value(value, form, root))
  }
  value
}
