# ==============================================================================
# Genome-Wide Assembly of Gene-Centric Coding-Mask Results
# ==============================================================================
#
# Combines the per-chromosome results of a gene-centric coding-mask scan into one
# genome-wide table and computes its significance. A coding-mask scan tests, for
# each gene, several functional variant categories (masks); this function gathers
# those results from up to three sources and lines them up for comparison:
#   - GLOW          : GLOW's own per-gene, per-category results.
#   - native STAAR  : STAARpipeline's STAAR-O results, run two ways -- without
#                     and with the saddlepoint approximation (SPA).
# It row-binds each source across chromosomes, joins the native STAAR results
# onto the GLOW results by (gene, category) to form a side-by-side method
# comparison, then derives per-category and genome-wide significance from the
# p-values (multiple-testing thresholds and counts are computed from the data,
# never hardcoded). Output rows are ordered gene-first: one block per gene at its
# genomic position, with that gene's variant categories beneath it.
#
# This function only READS a run's result files and RETURNS tables; it writes
# nothing to disk. The companion summary scripts persist the returned tables.
#
# EXPORTED FUNCTIONS:
#   - aggregate_coding_results()  3-source assembly + merge + significance
#
# INTERNAL HELPER FUNCTIONS:
#   - .acr_assert_same_schema()   Refuse to rbind tables whose column sets drift
#   - .acr_read_glow_chr()        rbind the GLOW per-chr CSVs + rename *_glowG
#   - .acr_read_native()          rbind one native mode's per array-task CSVs
#   - .acr_expected_arrayids()    Expected native array tasks for config$chroms
#   - .acr_report_native_only()   Category-aware native-not-in-GLOW reporting
#   - .acr_null_or()              a if non-NULL else b


#################### EXPORTED FUNCTIONS ####################

#' Assemble Gene-Centric Coding-Mask Results Genome-Wide
#'
#' Combines the per-chromosome output of a gene-centric coding-mask scan into one
#' genome-wide table and computes its significance. It reads up to three result
#' sources a coding run produces -- GLOW's own results, and (optionally)
#' STAARpipeline's STAAR-O results run without and with the saddlepoint
#' approximation (SPA) -- joins them by gene and variant category into a
#' side-by-side method comparison, and derives per-category and genome-wide
#' significance from the p-values. This is the coding-mask counterpart of
#' \code{\link{aggregate_scan_results}}.
#'
#' @details
#' \strong{Three sources.}
#' \enumerate{
#'   \item \emph{GLOW} -- the per-chromosome \code{glow_coding_chr*.csv}
#'     row-bound genome-wide. The two baked-in STAAR omnibus columns on GLOW's
#'     exact post-filter genotype matrix (\code{STAAR_O}, \code{ACAT_O}) are
#'     renamed \code{STAAR_O_glowG} / \code{ACAT_O_glowG} so the merge can tell
#'     them from the native engine's \code{STAAR_O_native}.
#'   \item \emph{native non-SPA} -- the \code{results-staar-nospa/staar_coding_*.csv}
#'     per array-task tables row-bound (columns suffixed \code{_native}).
#'   \item \emph{native SPA} -- the \code{results-staar-spa/staar_coding_*.csv}
#'     per array-task tables row-bound (columns suffixed \code{_spa}).
#' }
#'
#' \strong{Merge.} GLOW is the base: it emits a cell for every \eqn{\ge 1}-variant
#' gene; native only for \eqn{\ge 2}-variant cells, so native \eqn{\subseteq}
#' GLOW. A left join from GLOW captures every real cell; native cells with no GLOW
#' row are reported category-aware -- a cell in a GLOW-tested category is WARNED
#' (a genuine extraction mismatch), one in a standard category this run's GLOW set
#' does not test is INFO (expected: native always emits all 7 standard masks).
#'
#' \strong{Significance (computed at runtime, never hardcoded).} For each
#' \code{(category, method)} with \eqn{n_{tested} > 0} finite p-values, the
#' per-category Bonferroni cutoff is \eqn{\alpha / n_{tested}} and the genome-wide
#' cutoff is \code{gw_alpha}; counts of cells clearing each are tallied. The
#' significant-cells table lists every \code{(gene, category)} cell that any
#' present method makes significant (its per-\code{(category, method)} Bonferroni
#' cutoff, looked up from the counts table, or \code{gw_alpha}).
#'
#' \strong{Completeness.} The number of expected vs present GLOW chromosomes and
#' native array tasks is reconciled; \code{allow_partial = FALSE} (default) stops on
#' any gap (a missing chr/task distorts the significance denominators),
#' \code{TRUE} downgrades it to a warning and proceeds on what is present.
#'
#' \strong{Gene-first output.} \code{coding_merged} and
#' \code{coding_sig_by_category_cells} are sorted by
#' \code{(chr, gene_start, category)} -- one block per gene at its genomic
#' position with its categories underneath. \code{glow_coding_all},
#' \code{staar_coding_nospa_all} and \code{staar_coding_spa_all} keep their
#' \code{(chr, start[/gene], category)} order.
#'
#' @param output_dir Character. The run's \code{outputs/} directory; the function
#'   reads \code{results-glow/}, \code{results-staar-<mode>/} under it.
#' @param chroms Integer vector of the chromosomes the run covered (the GLOW
#'   completeness denominator, e.g. \code{config$chroms}).
#' @param gene_num_in_array Integer. Genes per native STAAR array task (the
#'   native array-task mapping, e.g. \code{config$gene_num_in_array}); used only
#'   to reconcile native completeness.
#' @param staar_modes Character vector of native modes to read (default
#'   \code{c("spa", "nospa")}; pass \code{character(0)} for a GLOW-only run).
#' @param genes_info Optional data frame whose second column gives each gene's
#'   chromosome (the \code{STAARpipeline::genes_info} layout). Used ONLY to derive
#'   the native array-task completeness denominator when \code{staar_modes} is
#'   non-empty; pass \code{STAARpipeline::genes_info} (or an equivalent table) to
#'   enable native reconciliation. Not needed for GLOW-only aggregation
#'   (\code{staar_modes = character(0)}) -- the package has no hard dependency on
#'   STAARpipeline.
#' @param alpha Family-wise error rate for the per-category Bonferroni cutoff
#'   (default \code{0.05}).
#' @param gw_alpha Genome-wide significance threshold (default \code{2.5e-6},
#'   the STAARpipeline convention).
#' @param primary_test Column ranking the significant-cells table and driving the
#'   within-category order (default \code{"GLOW_Omni"}).
#' @param allow_partial Logical (default \code{FALSE}). When \code{FALSE}, a
#'   missing GLOW chromosome file or native array task is an error; when \code{TRUE} it is a
#'   warning and the assembly proceeds on the present files.
#'
#' @return A named list:
#' \describe{
#'   \item{glow_coding_all}{Genome-wide GLOW table (gene x category), with the
#'     \code{*_glowG} rename and \code{label, chr, category} fronted.}
#'   \item{staar_coding_nospa_all, staar_coding_spa_all}{Genome-wide native
#'     tables per mode (0-row schema-only data.frame if a mode found nothing or is
#'     not requested).}
#'   \item{coding_merged}{One row per \code{(gene, category)}: the four method
#'     columns plus each source's sub-tests and count columns, sorted gene-first.}
#'   \item{coding_sig_by_category}{Per \code{(category, method)}: \code{n_tested},
#'     thresholds, and significant-cell counts.}
#'   \item{coding_sig_by_category_cells}{One row per significant
#'     \code{(gene, category)} cell, gene-first, with every method's p-value, the
#'     count chain, a global GLOW rank, and the significance flags.}
#'   \item{methods}{The method columns present in \code{coding_merged}.}
#'   \item{completeness}{Character vector reconciling expected vs present GLOW
#'     chromosomes and native array tasks (the completeness-report lines).}
#' }
#'
#' @seealso [aggregate_scan_results()], [run_scan_unit()], [build_scan_regions()]
#'
#' @references
#' Li, X., Li, Z., Zhou, H., et al. (2020). Dynamic incorporation of multiple in
#' silico functional annotations empowers rare variant association analysis of large
#' whole-genome sequencing studies at scale. Nature Genetics, 52, 969-983.
#' doi:10.1038/s41588-020-0676-4
#'
#' Li, Z., Li, X., Zhou, H., et al. (2022). A framework for detecting noncoding
#' rare-variant associations of large-scale whole-genome sequencing studies.
#' Nature Methods, 19, 1599-1611. doi:10.1038/s41592-022-01640-x
#'
#' @examples
#' \dontrun{
#' agg <- aggregate_coding_results(
#'   output_dir       = "runs/maf0.5_snv_7cat/outputs",
#'   chroms           = 1:22,
#'   gene_num_in_array = 500,
#'   staar_modes      = c("spa", "nospa"))
#' head(agg$coding_merged)
#' agg$coding_sig_by_category
#' }
#'
#' @export
aggregate_coding_results <- function(output_dir,
                                     chroms,
                                     gene_num_in_array,
                                     staar_modes = c("spa", "nospa"),
                                     genes_info = NULL,
                                     alpha = 0.05,
                                     gw_alpha = 2.5e-6,
                                     primary_test = "GLOW_Omni",
                                     allow_partial = FALSE) {

  stopifnot(dir.exists(output_dir))
  chroms <- as.integer(chroms)
  prov_lines <- character(0L)        # completeness provenance -> coding_completeness.txt

  # ===========================================================================
  # (1) GLOW: rbind per-chr glow_coding_chr*.csv, rename *_glowG
  # ===========================================================================
  glow_dir   <- file.path(output_dir, "results-glow")
  glow_files <- file.path(glow_dir, sprintf("glow_coding_chr%d.csv", chroms))
  present_glow <- file.exists(glow_files)
  missing_glow <- glow_files[!present_glow]
  prov_lines <- c(prov_lines,
    sprintf("GLOW per-chr: expected %d, present %d, missing %d%s",
            length(glow_files), sum(present_glow), length(missing_glow),
            if (length(missing_glow))
              paste0(" (chr ", paste(chroms[!present_glow], collapse = ","), ")")
            else ""))
  if (length(missing_glow) > 0L) {
    msg <- paste0("Missing GLOW per-chr file(s):\n  ",
                  paste(missing_glow, collapse = "\n  "))
    if (allow_partial) message("WARNING (allow_partial): ", msg)
    else stop(msg, call. = FALSE)
    glow_files <- glow_files[present_glow]
  }

  glow_all <- .acr_read_glow_chr(glow_files)

  # ===========================================================================
  # (2) native STAAR, both modes
  # ===========================================================================
  # The expected native array-task ids are only needed when reconciling native
  # STAAR results (non-empty staar_modes); they are derived from the caller-
  # supplied genes_info. GLOW-only aggregation (staar_modes = character(0)) skips
  # this entirely, so neither genes_info nor STAARpipeline is required for it.
  exp_ids <- if (length(staar_modes) > 0L)
    .acr_expected_arrayids(chroms, gene_num_in_array, genes_info) else integer(0L)
  staar_tables <- list()
  for (mode in staar_modes) {
    mdir <- file.path(output_dir, paste0("results-staar-", mode))
    if (!dir.exists(mdir)) {
      msg <- sprintf("Native STAAR dir missing for mode '%s': %s", mode, mdir)
      if (allow_partial) { message("WARNING (allow_partial): ", msg); next }
      else stop(msg, call. = FALSE)
    }
    present_ids <- as.integer(sub("^staar_coding_(\\d+)\\.csv$", "\\1",
                  list.files(mdir, pattern = "^staar_coding_\\d+\\.csv$")))
    missing_ids <- setdiff(exp_ids, present_ids)
    prov_lines <- c(prov_lines,
      sprintf("native STAAR [%s]: expected %d arrayid(s), present %d, missing %d%s",
              mode, length(exp_ids), length(intersect(exp_ids, present_ids)),
              length(missing_ids),
              if (length(missing_ids))
                paste0(" (e.g. ", paste(utils::head(missing_ids, 8), collapse = ", "), ")")
              else ""))
    if (length(missing_ids) > 0L) {
      msg <- sprintf(
        "Native STAAR mode '%s': %d/%d expected arrayid file(s) missing (e.g. %s)",
        mode, length(missing_ids), length(exp_ids),
        paste(utils::head(missing_ids, 8), collapse = ", "))
      if (allow_partial) message("WARNING (allow_partial): ", msg)
      else stop(msg, call. = FALSE)
    }
    flat <- .acr_read_native(mdir)
    if (is.null(flat)) {
      message(sprintf("  native STAAR mode '%s': 0 (gene, category) cells with results.",
                      mode))
      flat <- data.frame(gene = character(0L), chr = character(0L),
                         category = character(0L), stringsAsFactors = FALSE)
    }
    staar_tables[[mode]] <- flat
  }

  nospa_all <- .acr_null_or(staar_tables[["nospa"]],
                            data.frame(gene = character(0L), chr = character(0L),
                                       category = character(0L),
                                       stringsAsFactors = FALSE))
  spa_all   <- .acr_null_or(staar_tables[["spa"]],
                            data.frame(gene = character(0L), chr = character(0L),
                                       category = character(0L),
                                       stringsAsFactors = FALSE))

  # ===========================================================================
  # (3) Merge (GLOW base; left join native; native subset-of GLOW)  [from 04]
  # ===========================================================================
  glow <- glow_all
  names(glow)[names(glow) == "label"] <- "gene"   # uniform (gene, category) key

  # The native sources used in the merge: drop the duplicated chr column (GLOW's
  # is identical), and treat a 0-row table as "not present" (04's read_native).
  prep_native <- function(d) {
    if (is.null(d) || nrow(d) == 0L) return(NULL)
    d$chr <- NULL
    d
  }
  nospa <- if ("nospa" %in% staar_modes) prep_native(nospa_all) else NULL
  spa   <- if ("spa"   %in% staar_modes) prep_native(spa_all)   else NULL

  glow_keys <- paste(glow$gene, glow$category)
  .acr_report_native_only(nospa, "nospa", glow, glow_keys)
  .acr_report_native_only(spa,   "spa",   glow, glow_keys)

  merged <- glow
  if (!is.null(nospa))
    merged <- merge(merged, nospa, by = c("gene", "category"), all.x = TRUE)
  if (!is.null(spa))
    merged <- merge(merged, spa,   by = c("gene", "category"), all.x = TRUE)
  # Gene-first OUTPUT (design Q6): (chr, gene_start, category) -- one block per
  # gene at its genomic position with categories underneath. `start` is the gene
  # start carried in from GLOW. (Prior 04 sorted (chr, gene-alphabetical,
  # category); content is identical, only row order changes.)
  merged <- merged[order(.chr_rank(merged$chr), merged$start, merged$category),
                   , drop = FALSE]
  rownames(merged) <- NULL

  # ===========================================================================
  # (4) Per-category significance (computed from data; no hardcoded results)
  # ===========================================================================
  methods <- intersect(
    c("GLOW_Omni", "STAAR_O_glowG", "STAAR_O_native", "STAAR_B_native"),
    names(merged))
  sig_rows <- list()
  for (cat0 in sort(unique(merged$category))) {
    sub <- merged[merged$category == cat0, , drop = FALSE]
    for (meth in methods) {
      p <- sub[[meth]]
      n_tested <- sum(!is.na(p) & is.finite(p))
      if (n_tested == 0L) next
      bonf <- alpha / n_tested
      sig_rows[[length(sig_rows) + 1L]] <- data.frame(
        category = cat0, method = meth, n_tested = n_tested,
        min_p = min(p, na.rm = TRUE),
        bonferroni_threshold = bonf,
        n_sig_bonferroni = sum(p <= bonf, na.rm = TRUE),
        genomewide_threshold = gw_alpha,
        n_sig_genomewide = sum(p <= gw_alpha, na.rm = TRUE),
        stringsAsFactors = FALSE)
    }
  }
  sig <- if (length(sig_rows) > 0L) do.call(rbind, sig_rows) else
    data.frame(category = character(0L), method = character(0L),
               n_tested = integer(0L), stringsAsFactors = FALSE)

  # ===========================================================================
  # (5) Significant (gene, category) cells: the cells behind (4)        [from 04]
  # ===========================================================================
  if (!primary_test %in% names(merged))
    stop("coding_merged lacks the primary-test column '", primary_test,
         "' (primary_test).", call. = FALSE)
  bonf_thr <- stats::setNames(sig$bonferroni_threshold,
                              paste(sig$category, sig$method))

  flag_bonf <- matrix(FALSE, nrow(merged), length(methods),
                      dimnames = list(NULL, methods))
  flag_gw   <- flag_bonf
  for (meth in methods) {
    p   <- merged[[meth]]
    thr <- bonf_thr[paste(merged$category, meth)]
    ok  <- !is.na(p) & is.finite(p)
    flag_bonf[, meth] <- ok & !is.na(thr) & p <= thr
    flag_gw[,   meth] <- ok & p <= gw_alpha
  }
  join_methods <- function(flag_row) paste(methods[flag_row], collapse = ",")
  cells <- merged
  is_significant         <- rowSums(flag_bonf) > 0L | rowSums(flag_gw) > 0L
  cells$sig_methods_bonf <- if (length(methods))
    apply(flag_bonf, 1L, join_methods) else ""
  cells$sig_methods_gw   <- if (length(methods))
    apply(flag_gw,   1L, join_methods) else ""

  # Global rank by the primary GLOW test (ascending p; non-finite last).
  cells$rank <- rank(
    ifelse(is.finite(cells[[primary_test]]), cells[[primary_test]], Inf),
    ties.method = "first")

  # Keep only the significant cells. Gene-first OUTPUT (design Q6): order by
  # (chr, gene_start, category) so the cells sit at their genomic position, one
  # gene block at a time. (Prior 04 grouped by category then ascending p; content
  # is the same set of cells, only the row order changes.)
  cells <- cells[is_significant, , drop = FALSE]
  cells <- cells[order(.chr_rank(cells$chr), cells$start, cells$category),
                 , drop = FALSE]
  rownames(cells) <- NULL

  # Column order: cell key/info, the present method p-values, then annotations.
  info_cols <- intersect(
    c("gene", "chr", "start", "end", "category",
      "n_variants", "n_after_annotation", "n_after_filter",
      "n_after_ld", "n_after_collapse", "cMAC", "is_single"), names(cells))
  ann_cols  <- c("rank", "sig_methods_bonf", "sig_methods_gw")
  front     <- c(info_cols, methods, ann_cols)
  cells <- cells[, c(front, setdiff(names(cells), front)), drop = FALSE]

  list(
    glow_coding_all              = glow_all,
    staar_coding_nospa_all       = nospa_all,
    staar_coding_spa_all         = spa_all,
    coding_merged                = merged,
    coding_sig_by_category       = sig,
    coding_sig_by_category_cells = cells,
    methods                      = methods,
    completeness                 = prov_lines)
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Refuse to rbind per array-task tables whose column SETS drift
#'
#' rbind matches columns by name, so only the column SET must agree (a pure order
#' difference is harmless); a differing set is real drift. Names the offending
#' file.
#'
#' @param parts List of data.frames.
#' @param files Their file paths (for the error message).
#' @param what A label for the source in the error message.
#' @return Invisibly NULL; stops on drift.
#' @keywords internal
#' @noRd
.acr_assert_same_schema <- function(parts, files, what) {
  if (length(parts) <= 1L) return(invisible())
  cols1 <- colnames(parts[[1L]])
  for (k in seq_along(parts)[-1L]) {
    if (!setequal(colnames(parts[[k]]), cols1)) {
      extra   <- setdiff(colnames(parts[[k]]), cols1)
      missing <- setdiff(cols1, colnames(parts[[k]]))
      stop(sprintf("%s: per-file column SETS differ -- %s vs %s (extra: %s; missing: %s).",
                   what, basename(files[1L]), basename(files[k]),
                   if (length(extra))   paste(extra,   collapse = ",") else "none",
                   if (length(missing)) paste(missing, collapse = ",") else "none"),
           call. = FALSE)
    }
  }
  invisible()
}


#' rbind the GLOW per-chr CSVs + rename the baked-in STAAR omnibi to *_glowG
#'
#' @param glow_files Character vector of \code{glow_coding_chr*.csv} paths
#'   (already filtered to present files).
#' @return The genome-wide GLOW table with \code{STAAR_O}/\code{ACAT_O} renamed
#'   to \code{STAAR_O_glowG}/\code{ACAT_O_glowG}, ordered \code{(chr, start)},
#'   leading with \code{label, chr, category}.
#' @keywords internal
#' @noRd
.acr_read_glow_chr <- function(glow_files) {
  glow_parts <- lapply(glow_files, utils::read.csv, stringsAsFactors = FALSE,
                       colClasses = c(chr = "character"))
  keep <- vapply(glow_parts, nrow, integer(1)) > 0L
  glow_parts <- glow_parts[keep]; glow_kept_files <- glow_files[keep]
  if (length(glow_parts) == 0L)
    stop("No non-empty GLOW per-chr tables found.", call. = FALSE)
  .acr_assert_same_schema(glow_parts, glow_kept_files, "GLOW per-chr tables")
  glow_all <- do.call(rbind, glow_parts)

  # Rename the baked-in STAAR omnibi on GLOW's G so the merge can tell them from
  # the native engine's STAAR_O_native (design §6.5).
  ren <- c(STAAR_O = "STAAR_O_glowG", ACAT_O = "ACAT_O_glowG")
  for (old in names(ren)) {
    j <- which(names(glow_all) == old)
    if (length(j) == 1L) names(glow_all)[j] <- unname(ren[[old]])
  }
  glow_all <- glow_all[order(as.integer(glow_all$chr), glow_all$start),
                       , drop = FALSE]
  # Front the (gene, chr, category) key -- gene first, consistent with native.
  front <- c("label", "chr", "category")
  front <- front[front %in% names(glow_all)]
  glow_all <- glow_all[, c(front, setdiff(names(glow_all), front)), drop = FALSE]
  rownames(glow_all) <- NULL
  glow_all
}


#' rbind one native mode's per array-task CSVs (curated mode-keyed schema)
#'
#' The coding scan writes one flat per-array-task CSV per mode, so the genome assembler
#' is a plain rbind (no nested-list walk). 0-row CSVs are dropped.
#'
#' @param results_dir A \code{results-staar-<mode>/} directory.
#' @return The mode's genome-wide table ordered \code{(chr, gene, category)}, or
#'   \code{NULL} if no non-empty per array-task CSV exists.
#' @keywords internal
#' @noRd
.acr_read_native <- function(results_dir) {
  files <- list.files(results_dir, pattern = "^staar_coding_\\d+\\.csv$",
                      full.names = TRUE)
  parts <- lapply(files, function(f)
    utils::read.csv(f, stringsAsFactors = FALSE, colClasses = c(chr = "character")))
  keep  <- vapply(parts, nrow, integer(1)) > 0L
  parts <- parts[keep]; files <- files[keep]
  if (length(parts) == 0L) return(NULL)
  .acr_assert_same_schema(parts, files,
                          paste0("native STAAR per array-task (", basename(results_dir), ")"))
  out <- do.call(rbind, parts)
  out <- out[order(as.integer(out$chr), out$gene, out$category), , drop = FALSE]
  rownames(out) <- NULL
  out
}


#' Expected native STAAR array tasks for config$chroms
#'
#' Mirrors the native array-task mapping: genes per chromosome (from the supplied
#' \code{genes_info}) chunked by \code{gene_num_in_array}, cumulated, then the
#' run's chromosomes' id ranges concatenated.
#'
#' @param chroms Integer chromosomes the run covered.
#' @param gene_num_in_array Genes per array task.
#' @param genes_info Data frame whose 2nd column is each gene's chromosome (the
#'   STAARpipeline::genes_info layout). Passed in by the caller so this package
#'   carries no hard STAARpipeline dependency.
#' @return Integer vector of expected array-task ids.
#' @keywords internal
#' @noRd
.acr_expected_arrayids <- function(chroms, gene_num_in_array, genes_info) {
  # genes_info supplies the native array-task layout (2nd column = each gene's
  # chromosome). It is passed in by the caller (e.g. STAARpipeline::genes_info)
  # so this package needs no hard STAARpipeline dependency; reached only for
  # native STAAR reconciliation (non-empty staar_modes).
  if (is.null(genes_info)) {
    stop("Native STAAR completeness (non-empty `staar_modes`) needs a ",
         "`genes_info` table (each gene's chromosome in column 2), e.g. ",
         "genes_info = STAARpipeline::genes_info. Use staar_modes = character(0) ",
         "for GLOW-only aggregation.", call. = FALSE)
  }
  g  <- ceiling(table(genes_info[, 2]) / gene_num_in_array)
  cs <- cumsum(g)
  ids <- integer(0L)
  for (chr in chroms) {
    lo <- if (chr == 1L) 1L else cs[chr - 1L] + 1L
    ids <- c(ids, lo:cs[chr])
  }
  as.integer(ids)
}


#' Category-aware reporting of native (gene, category) cells with no GLOW row
#'
#' Native subset-of GLOW under matched masks, so a native cell with no GLOW row
#' splits two ways (design §6.1): (1) a category GLOW ALSO tests -> a genuine
#' mismatch (WARN; check extraction); (2) a standard category GLOW does NOT test
#' in this run (subset/custom-category run) -> expected (INFO; native always
#' emits all 7 standard masks).
#'
#' @param native A prepped native table (chr dropped) or NULL.
#' @param label "nospa" / "spa", for the message.
#' @param glow The GLOW base table (for the tested-category set).
#' @param glow_keys \code{paste(glow$gene, glow$category)}.
#' @return Invisibly NULL; emits messages only.
#' @keywords internal
#' @noRd
.acr_report_native_only <- function(native, label, glow, glow_keys) {
  if (is.null(native)) return(invisible())
  glow_cats <- unique(glow$category)
  nk        <- paste(native$gene, native$category)
  extra     <- native[!(nk %in% glow_keys), c("gene", "category"), drop = FALSE]
  if (nrow(extra) == 0L) return(invisible())
  mism   <- extra[  extra$category %in% glow_cats, , drop = FALSE]   # (1)
  noncat <- extra[!(extra$category %in% glow_cats), , drop = FALSE]  # (2)
  if (nrow(mism) > 0L) {
    message(sprintf("NOTE: %d native-%s cell(s) in a GLOW-tested category have no GLOW row ",
                    nrow(mism), label),
            "(unexpected under matched masks; check extraction): ",
            paste(utils::head(paste(mism$gene, mism$category), 5), collapse = "; "))
  }
  if (nrow(noncat) > 0L) {
    cats <- sort(unique(noncat$category))
    message(sprintf("  native-%s also tested %d cell(s) in %d standard categor%s NOT in this run's GLOW set ",
                    label, nrow(noncat), length(cats),
                    if (length(cats) == 1L) "y" else "ies"),
            sprintf("(expected; native emits all 7 standard masks): %s",
                    paste(cats, collapse = ", ")))
  }
}


#' Return a if Non-NULL, Else b
#' @keywords internal
#' @noRd
.acr_null_or <- function(a, b) if (is.null(a)) b else a
