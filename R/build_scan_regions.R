# ==============================================================================
# Region Dispatch for the Unified GLOW Scan (build_scan_regions)
# ==============================================================================
#
# The one structural fork of the unified scan: gene/window run ONE
# filter over many regions; coding runs MANY filters (one per category) over the
# SAME regions.
# build_scan_regions() returns a per-chromosome `regions` table PLUS either a
# single `filters` object (gene/window) or a NAMED category -> filter map
# (coding). run_scan_unit() then loops glow_scan_chunk() once (single filter) or
# once per category (coding). The per-type region/filter construction:
#   - gene  : define_regions_gene(chr) (or config$resolved_gene_list subset for
#             a custom gene list) + do.call(variant_filter, config$filter_spec)
#   - window: define_regions_window(...) + the boundary trim (start <= chunk_end)
#   - coding: define_regions_gene(chr) + a category -> filter map via
#             resolve_category_filter()
#
# EXPORTED FUNCTIONS:
#   - build_scan_regions()   regions table + filter(s) for one chr/chunk
#
# INTERNAL HELPER FUNCTIONS:
#   - .bsr_get()                  Read one config field from an env/list.
#   - .resolve_category_filter()  One category spec -> a glow_variant_filter
#                                 (builtin name -> coding_filter; DNF clauses ->
#                                 variant_filter), applying the run's filter_spec.


#################### EXPORTED MAIN FUNCTIONS ####################

#' Build the Regions Table and Filter(s) for One Scan Unit
#'
#' Region dispatch for the unified GLOW scan. Returns the
#' per-chromosome \code{regions} table plus the variant filter(s) for a single
#' scan unit (a chromosome for gene/coding, a genome chunk for window). This is
#' the one structural fork of the scan: \code{gene}/\code{window} run a single
#' filter over many regions, while \code{coding} runs many filters (one per
#' category) over the same regions, so the return carries either one filter
#' object or a named \code{category -> filter} map.
#'
#' @details
#' \strong{\code{region_type = "gene"}.} \code{regions} is
#' \code{GLOWr::define_regions_gene(chr = chr)} — the bundled STAARpipeline gene
#' table for the chromosome — unless \code{config$resolved_gene_list} is present
#' (a user-supplied custom gene list), in which case that list is subset to
#' \code{chr}. \code{filters} is the single filter
#' \code{do.call(GLOWr::variant_filter, config$filter_spec)}.
#'
#' \strong{\code{region_type = "window"}.} Requires a \code{chunk} row (one row
#' of the genome-chunk table from \code{\link{define_genome_chunks}}, with
#' \code{chunk_start}, \code{chunk_end}, \code{span_end}). \code{regions} is
#' \code{GLOWr::define_regions_window(chr, chunk$chunk_start, chunk$span_end,
#' config$window_size, config$step_size)} trimmed by the sliding boundary rule
#' (keep \code{start <= chunk$chunk_end}): windows are generated to the
#' chromosome \code{span_end} so boundary-spanning windows are kept, then only
#' those owned by this chunk (start within \code{[chunk_start, chunk_end]}) are
#' retained, so the union of chunks equals one full-span slide. \code{filters} is
#' the single \code{variant_filter} from \code{config$filter_spec}.
#'
#' \strong{\code{region_type = "coding"}.} \code{regions} is the chromosome gene
#' table (as for \code{gene}). \code{filters} is a NAMED list mapping each
#' category label to its resolved \code{glow_variant_filter}:
#' \code{config$categories} is either a character vector of builtin coding
#' category names (each resolved via \code{GLOWr::coding_filter}) or a named list
#' mapping labels to a builtin name OR a list of custom DNF annotation clauses
#' (resolved via \code{GLOWr::variant_filter} with the run's
#' \code{filter_spec}). The MAF/MAC/min-variant axes come from
#' \code{config$filter_spec} so a category never hardcodes them.
#'
#' @param region_type Character: one of \code{"gene"}, \code{"window"},
#'   \code{"coding"}.
#' @param chr Integer or character chromosome.
#' @param config An environment or named list of resolved config values. Fields
#'   used per type: \code{filter_spec} (all); \code{resolved_gene_list}
#'   (optional, gene); \code{window_size}, \code{step_size} (window);
#'   \code{categories} (coding).
#' @param chunk A single-row data.frame (one row of the genome-chunk table) for
#'   \code{region_type = "window"}; ignored otherwise. Must carry
#'   \code{chunk_start}, \code{chunk_end}, \code{span_end}. Default NULL.
#'
#' @return A named list:
#'   \describe{
#'     \item{regions}{A \code{data.frame} with at least \code{chr}, \code{start},
#'       \code{end}, \code{label} (the \code{define_regions_*} shape).}
#'     \item{filters}{A single \code{glow_variant_filter} (gene/window) OR a
#'       named list mapping each category label to a \code{glow_variant_filter}
#'       (coding).}
#'   }
#'
#' @examples
#' \dontrun{
#' env <- new.env(); source("runs/myrun/config.R", local = env)
#' # Gene: one filter, the chr gene table.
#' g <- build_scan_regions("gene", chr = 22, config = env)
#' # Window: trimmed to the chunk boundary.
#' chunk <- define_genome_chunks(22, window_size = env$window_size,
#'                               step_size = env$step_size)[1, ]
#' w <- build_scan_regions("window", chr = 22, config = env, chunk = chunk)
#' # Coding: a category -> filter map.
#' co <- build_scan_regions("coding", chr = 22, config = env)
#' names(co$filters)
#' }
#'
#' @seealso \code{GLOWr::define_regions_gene}, \code{GLOWr::define_regions_window},
#'   \code{GLOWr::variant_filter}, \code{GLOWr::coding_filter},
#'   \code{\link{run_scan_unit}}
#' @export
build_scan_regions <- function(region_type, chr, config, chunk = NULL) {

  region_type <- match.arg(region_type, c("gene", "window", "coding"))
  filter_spec <- .bsr_get(config, "filter_spec")
  if (is.null(filter_spec)) {
    stop("`config$filter_spec` is required (the list of variant_filter() args).",
         call. = FALSE)
  }

  if (region_type == "gene") {
    # Gene table: bundled STAARpipeline genes, OR the custom gene list
    # subset to this chr (a build_scan_regions plug-in). Subsetting on chr keeps
    # the resolved-list shape identical to define_regions_gene(chr).
    resolved_gene_list <- .bsr_get(config, "resolved_gene_list")
    if (!is.null(resolved_gene_list)) {
      regions <- resolved_gene_list[
        as.character(resolved_gene_list$chr) == as.character(chr), ,
        drop = FALSE]
      rownames(regions) <- NULL
    } else {
      regions <- GLOWr::define_regions_gene(chr = as.character(chr))
    }
    filters <- do.call(GLOWr::variant_filter, filter_spec)
    return(list(regions = regions, filters = filters))
  }

  if (region_type == "window") {
    if (is.null(chunk) || !is.data.frame(chunk) || nrow(chunk) != 1L) {
      stop("`region_type = \"window\"` requires a single-row `chunk` (a row of ",
           "the genome-chunk table).", call. = FALSE)
    }
    needed <- c("chunk_start", "chunk_end", "span_end")
    miss <- setdiff(needed, names(chunk))
    if (length(miss) > 0L) {
      stop("`chunk` is missing column(s): ", paste(miss, collapse = ", "),
           " (expected the define_genome_chunks() shape).", call. = FALSE)
    }
    window_size <- .bsr_get(config, "window_size")
    step_size   <- .bsr_get(config, "step_size")
    if (is.null(window_size) || is.null(step_size)) {
      stop("`config$window_size` and `config$step_size` are required for ",
           "region_type = \"window\".", call. = FALSE)
    }
    # Boundary rule: generate to
    # the chromosome span_end so boundary-spanning windows are kept, then keep
    # only windows whose start is owned by this chunk (start <= chunk_end).
    regions <- GLOWr::define_regions_window(chr, chunk$chunk_start,
                                            chunk$span_end,
                                            window_size, step_size)
    regions <- regions[regions$start <= chunk$chunk_end, , drop = FALSE]
    rownames(regions) <- NULL
    filters <- do.call(GLOWr::variant_filter, filter_spec)
    return(list(regions = regions, filters = filters))
  }

  # region_type == "coding": same regions, a category -> filter map.
  regions <- GLOWr::define_regions_gene(chr = as.character(chr))
  category_specs  <- .bsr_get(config, "categories")
  if (is.null(category_specs)) {
    stop("`config$categories` is required for region_type = \"coding\".",
         call. = FALSE)
  }
  category_labels <- if (is.list(category_specs)) names(category_specs) else
    category_specs
  if (is.list(category_specs) &&
      (is.null(category_labels) || any(!nzchar(category_labels)))) {
    stop("config$categories given as a list must be NAMED (each entry needs a ",
         "category label).", call. = FALSE)
  }
  filters <- stats::setNames(
    lapply(category_labels, function(category) {
      spec <- if (is.list(category_specs)) category_specs[[category]] else category
      .resolve_category_filter(spec, filter_spec)
    }),
    category_labels)
  list(regions = regions, filters = filters)
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Read One Config Field from an Environment or List
#'
#' @param config Environment or named list of config values.
#' @param name Character. Field name.
#' @return The field value, or NULL if unset/absent.
#' @keywords internal
#' @noRd
.bsr_get <- function(config, name) {
  if (is.environment(config)) {
    get0(name, envir = config, inherits = FALSE)
  } else if (is.list(config)) {
    config[[name]]
  } else {
    stop("`config` must be an environment or a named list.")
  }
}

#' Resolve One Category Spec to a glow_variant_filter
#'
#' Resolve one category spec to a filter. A length-1 string
#' is a builtin coding-category name (\code{GLOWr::coding_filter}); a non-empty
#' list is a set of custom DNF annotation clauses
#' (\code{GLOWr::variant_filter}). Either way the run's \code{filter_spec}
#' supplies the MAF axis, \code{min_mac}, \code{min_variants} (and
#' \code{variant_type} for the custom path), so a category never hardcodes them.
#'
#' @param spec A length-1 string (builtin name) or a non-empty list (DNF clauses).
#' @param fs The run's \code{filter_spec} (a list of \code{variant_filter} args).
#' @return A \code{glow_variant_filter}.
#' @keywords internal
#' @noRd
.resolve_category_filter <- function(spec, fs) {
  if (is.character(spec) && length(spec) == 1L) {
    GLOWr::coding_filter(spec, rare_maf_cutoff = fs$rare_maf_cutoff,
                         min_mac = fs$min_mac, min_variants = fs$min_variants)
  } else if (is.list(spec) && length(spec) > 0L) {
    GLOWr::variant_filter(variant_type = fs$variant_type,
                          annotation_clauses = spec,
                          rare_maf_cutoff = fs$rare_maf_cutoff,
                          min_mac = fs$min_mac, min_variants = fs$min_variants)
  } else {
    stop("Each category must be a builtin coding-category name (length-1 ",
         "string) or a non-empty list of DNF annotation clauses; got class '",
         class(spec)[1L], "'.", call. = FALSE)
  }
}
