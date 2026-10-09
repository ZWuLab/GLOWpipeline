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
#   - .segment_large_genes()      Tile genes above gene_max_records into
#                                 gene_segment_bp segments (region_type gene).
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
#' @param variant_index A \code{glow_variant_index} of the unit's GDS from
#'   \code{GLOWr::build_variant_index()}, or NULL (default). Used only by
#'   \code{region_type = "gene"}, for the large-gene segmentation below.
#'
#' @section Large-gene segmentation (gene):
#' At whole-genome-sequencing density a gene body can hold tens of thousands
#' of records, and the per-gene cost and memory grow with the square of the
#' variant count. When \code{variant_index} is supplied and
#' \code{config$gene_max_records} is a finite number, every gene's GDS records
#' inside its span are counted from the index (before any filter, so the count
#' depends only on the file), and a gene with more records than that is
#' replaced by non-overlapping segments of \code{config$gene_segment_bp} base
#' pairs (default 10,000) that start at the gene start; the last segment ends
#' at the gene end. Each segment is its own region row with the gene's
#' \code{region_id}, its own \code{start} and \code{end}, and
#' \code{label = <gene label>_seg<k>}, and it is tested and reported like a
#' gene; there is no gene-level combination of the segments. A gene whose span
#' fits one segment is left as it is. The returned table then also carries
#' \code{n_records} (the index count of each row) and the attributes
#' \code{n_genes_segmented} and \code{n_segments}. With \code{variant_index}
#' NULL, or \code{gene_max_records} NULL or infinite, nothing is segmented.
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
build_scan_regions <- function(region_type, chr, config, chunk = NULL,
                               variant_index = NULL) {

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
    # Large-gene segmentation (see the roxygen section): needs the index and a
    # finite threshold; otherwise the gene table is returned as built.
    gene_max_records <- .bsr_get(config, "gene_max_records")
    if (!is.null(variant_index) && !is.null(gene_max_records) &&
        is.finite(gene_max_records)) {
      gene_segment_bp <- .bsr_get(config, "gene_segment_bp")
      if (is.null(gene_segment_bp)) gene_segment_bp <- 10000L
      regions <- .segment_large_genes(regions, variant_index,
                                      gene_max_records, gene_segment_bp)
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

#' Tile Genes Above a Record Count Into Fixed-Length Segments
#'
#' Counts the GDS records in each gene's span from the variant index and
#' replaces every gene above \code{gene_max_records} by non-overlapping
#' segments of \code{gene_segment_bp} from the gene start (the last one ends
#' at the gene end). Segment rows keep the gene's \code{region_id} and carry
#' \code{label = <label>_seg<k>}. The result carries \code{n_records} per row
#' and the attributes \code{n_genes_segmented} and \code{n_segments}.
#'
#' @param regions The gene table (\code{region_id, chr, start, end, label}).
#' @param variant_index A \code{glow_variant_index}.
#' @param gene_max_records Numeric threshold on the record count.
#' @param gene_segment_bp Integer segment length in base pairs.
#' @return The regions table, segmented where needed, in genomic order.
#' @keywords internal
#' @noRd
.segment_large_genes <- function(regions, variant_index, gene_max_records,
                                 gene_segment_bp) {
  if (!inherits(variant_index, "glow_variant_index")) {
    stop("`variant_index` must be a glow_variant_index from ",
         "GLOWr::build_variant_index().", call. = FALSE)
  }
  gene_segment_bp <- as.integer(gene_segment_bp)
  if (is.na(gene_segment_bp) || gene_segment_bp < 1L) {
    stop("`config$gene_segment_bp` must be a positive integer (base pairs).",
         call. = FALSE)
  }
  if (nrow(regions) == 0L) {
    regions$n_records <- integer(0)
    attr(regions, "n_genes_segmented") <- 0L
    attr(regions, "n_segments") <- 0L
    return(regions)
  }
  n_rec <- GLOWr::count_index_records(variant_index, regions$chr,
                                      regions$start, regions$end)
  span  <- regions$end - regions$start + 1L
  split_gene <- n_rec > gene_max_records & span > gene_segment_bp

  rows <- lapply(seq_len(nrow(regions)), function(i) {
    r <- regions[i, , drop = FALSE]
    if (!split_gene[i]) {
      r$n_records <- n_rec[i]
      return(r)
    }
    starts <- seq.int(r$start, r$end, by = gene_segment_bp)
    ends   <- pmin(starts + gene_segment_bp - 1L, r$end)
    seg <- r[rep(1L, length(starts)), , drop = FALSE]
    seg$start <- as.integer(starts)
    seg$end   <- as.integer(ends)
    seg$label <- sprintf("%s_seg%d", r$label, seq_along(starts))
    seg$n_records <- GLOWr::count_index_records(variant_index, seg$chr,
                                                seg$start, seg$end)
    seg
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "n_genes_segmented") <- sum(split_gene)
  attr(out, "n_segments") <- sum(vapply(rows[split_gene], nrow, integer(1)))
  out
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
