# ==============================================================================
# Load a Run's Saved Configuration
# ==============================================================================
#
# Reads the resolved configuration a scan run saved during its preparation step
# (<output_dir>/shared/01-config_used.rds) and returns it as a list, for the
# summary and reporting steps to use.
#
# EXPORTED FUNCTIONS:
#   - load_config_snapshot()   Read <output_dir>/shared/01-config_used.rds


#################### EXPORTED MAIN FUNCTIONS ####################

#' Load a Run's Saved Configuration
#'
#' Reads the resolved configuration a scan run saved during its preparation step
#' (\code{<output_dir>/shared/01-config_used.rds}) and returns it as a list, for
#' the summary and reporting steps to use.
#'
#' @details
#' Errors if the file is absent (the run's preparation step has not been run
#' yet). Cost is that of a single \code{readRDS}.
#'
#' @param output_dir Character. The run's \code{outputs/} directory. Required.
#'
#' @return The resolved-configuration list saved by the run's preparation step.
#'
#' @export
load_config_snapshot <- function(output_dir) {
  snap_path <- file.path(output_dir, "shared", "01-config_used.rds")
  if (!file.exists(snap_path)) {
    stop("Pair-1 config snapshot not found at ", snap_path,
         ". Run 01-prepare-shared.R first.")
  }
  readRDS(snap_path)
}
