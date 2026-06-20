# ==============================================================================
# Load Per-Chromosome Evidence Packets
# ==============================================================================
#
# Loads the per-chr evidence packets for a vector of chromosomes and returns a
# flat named list keyed by region label, so a drill-down stage need not
# materialise all 22 packets when only a handful of chrs contain top hits.
#
# EXPORTED FUNCTIONS:
#   - load_evidence_for_chrs()   Load per-chr evidence packets -> flat named list


#################### EXPORTED MAIN FUNCTIONS ####################

#' Load Per-Chromosome Evidence Packets into a Flat Named List
#'
#' Loads the per-chr evidence packets (\code{evidence_chr\{N\}.rds}) for a vector
#' of chromosomes and returns a flat named list keyed by region label. Used by a
#' drill-down stage to avoid materialising all 22 packets when only a handful of
#' chrs contain top hits.
#'
#' @details
#' Packets are flat-merged across chromosomes; genome-wide region labels are
#' unique, so no label collides. Errors loudly if any requested chr's packet is
#' missing.
#'
#' @param evidence_dir Character. Directory holding the per-chr evidence packets.
#'   Required.
#' @param chrs Integer vector of chromosomes whose packets to load.
#'
#' @return A flat named list of region evidence-packet entries, keyed by region
#'   label.
#'
#' @export
load_evidence_for_chrs <- function(evidence_dir, chrs) {
  out <- list()
  for (c in chrs) {
    f <- file.path(evidence_dir, sprintf("evidence_chr%d.rds", c))
    if (!file.exists(f)) {
      stop("Missing evidence packet: ", f)
    }
    ev <- readRDS(f)
    # Flat-merge; if any region label collides across chrs we'd lose data.
    # Genome-wide region labels are unique (verified in the 2026-05-14
    # closeout), so this is safe.
    out <- c(out, ev)
  }
  out
}
