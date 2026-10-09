# ==============================================================================
# Shared Preparation Step for Multi-Region GLOW Scans (prepare_scan_run)
# ==============================================================================
#
# Consolidates the near-identical preparation that the SNV-set scan
# types (gene, sliding-window, coding) had each forked into their own
# `01-prepare-shared.R`. The SPINE captured here is: output-dir creation,
# GDS-presence check, phenotype load + sample-alignment + GLOWr::fit_null_model,
# PI-model load + PI-feature resolve, B-source dispatch (`b_func` XOR
# `b_model_path` + probe), input checksums, and the BASE config_used snapshot
# fields. Each module's TYPE-SPECIFIC extras (STAAR nulls, gene-list resolver,
# annotation-median cache, chunk table, GC calibration, ...) stay in that
# module's own preparation script and run AFTER prepare_scan_run; each script also keeps
# WRITING its own artifacts, so every module's exact on-disk layout/snapshot is
# preserved byte-identical.
#
# EXPORTED FUNCTIONS:
#   - prepare_scan_run()   Resolve + validate the shared preparation inputs; fit the
#                          null model; return a bundle (and write shared
#                          01-null_model.rds + 01-input_checksums.rds).
#
# ALSO EXPORTED:
#   - ensure_shared_output_dir()  Best-effort, portable hook that makes a run's
#                          output_dir a shared-backed symlink (the design's
#                          sanctioned integration point for opt-in born-sharing).
#
# INTERNAL HELPERS:
#   - .psr_get()           Read one config field from an env/list (NULL if unset).
#   - .psr_has()           TRUE if a config field is set (exists and non-NULL).
#   - .psr_null_or()       a if non-NULL else b.


#################### EXPORTED MAIN FUNCTIONS ####################

#' Resolve the Shared Preparation Step of a Multi-Region GLOW Scan
#'
#' Captures the preparation logic common to every GLOW SNV-set scan module
#' (gene, sliding-window, or coding): create the output directory tree,
#' validate that every per-chromosome GDS file is present, load the prepared
#' phenotype/covariate object and align its samples to GDS sample order, fit the
#' GLOWr null model, load the PI models and resolve the PI feature set, dispatch
#' the B source (\code{b_func} \emph{xor} \code{b_model_path}, with a probe/load
#' sanity check), compute input checksums, and assemble the \emph{base} fields of
#' the resolved-config snapshot.
#'
#' @details
#' \strong{What it does (the spine).} The steps below capture the preparation
#' shared by every scan type, so a run's prepared outputs are stable and
#' reproducible:
#' \enumerate{
#'   \item Create \code{<output_dir>/shared} (always) plus any directories named
#'     in \code{make_dirs}.
#'   \item Validate \code{gds_dir}, \code{pheno_path}, \code{pi_model_dir} exist;
#'     build per-chromosome GDS paths for \code{chroms} from \code{gds_pattern}
#'     and error if any are missing.
#'   \item Dispatch the B source: exactly one of \code{b_func} (a closed-form
#'     function, probed at MAF 0.1 for a finite scalar) and \code{b_model_path}
#'     (a trained model, loaded via \code{GLOWr::load_B_model}) must be set.
#'   \item Read \code{pheno_path}; intersect its \code{sample_id} with the first
#'     requested chromosome's GDS sample set; align \code{Y}/\code{X} to GDS
#'     order.
#'   \item Fit the null model with \code{GLOWr::fit_null_model} (trait from
#'     \code{prep$trait}, default \code{"binary"}).
#'   \item Load the PI models (\code{GLOWr::load_PI_models}); resolve PI features
#'     (config \code{pi_features}, else \code{GLOWr:::.default_PI_features()}).
#'   \item Compute MD5 checksums of every input file (passive audit trail).
#' }
#'
#' \strong{What it does NOT do.} Type-specific preparation work stays in the calling
#' script and runs after \code{prepare_scan_run}: STAAR null models, the custom
#' gene-list resolver, the per-chromosome annotation-median cache + chromosome
#' spans, the genome chunk table, the essentialdb annotation catalog, and GC
#' calibration. Those are folded into each module's own snapshot and artifact
#' writes so each scan type's outputs stay reproducible.
#'
#' \strong{Side effects (writes).} Writes \code{<output_dir>/shared/01-null_model.rds}
#' and \code{<output_dir>/shared/01-input_checksums.rds}. The caller still writes
#' \code{01-config_used.rds} (base fields from \code{base_snapshot} + its extras)
#' and any module-specific artifacts. \strong{Born-shared:} when run inside a repo
#' that ships \code{tools/shared-storage/scaffold_shared_output.sh}, a not-yet-existing
#' \code{output_dir} is first created as a symlink into the lab-shared storage (so the
#' run's outputs live off \code{/home} and are registered for clone reconstruction);
#' this is a portable no-op otherwise (plain local dir), and the output \emph{content}
#' is unchanged either way. \strong{Controlled per run by the config field}
#' \code{symlinked_shared_root}: a directory path symlinks \code{output_dir} into
#' \code{<path>/<repo-relative-leaf>}; \code{NULL} keeps it a plain local dir; an
#' \emph{absent} field uses the shared-storage tooling's default root
#' (\code{tools/shared-storage/shared-root.conf}). Set it in the run's \code{config.R} /
#' base-config.
#'
#' \strong{Complexity.} Dominated by \code{fit_null_model} (one GLM fit on the
#' aligned cohort) and one GDS open per requested chromosome for presence/sample
#' checks. Memory is O(n samples \eqn{\times} covariates) for the design matrix.
#'
#' @param config An environment or named list of resolved configuration values
#'   (typically the environment populated by \code{source(config_path,
#'   local = TRUE)}). Required fields: \code{gds_dir}, \code{gds_pattern},
#'   \code{pheno_path}, \code{pi_model_dir}, \code{filter_spec}. Optional:
#'   \code{pi_features} (NULL -> GLOWr defaults), \code{b_func},
#'   \code{b_model_path} (exactly one of these two), and \code{symlinked_shared_root}
#'   (where a new \code{output_dir} physically lives — see Born-shared below: a path
#'   symlinks outputs into it; \code{NULL} keeps them local; absent uses the default
#'   shared root). The resolved value is recorded in the snapshot.
#' @param output_dir Character. Run output directory; \code{shared/} and any
#'   \code{make_dirs} are created beneath it.
#' @param chroms Integer vector of chromosomes whose GDS files must be present
#'   (default \code{1:22}). The first element's GDS is opened for sample
#'   alignment.
#' @param make_dirs Character vector of additional directory names (relative to
#'   \code{output_dir}) to create alongside \code{shared/}. Default
#'   \code{character(0)} (only \code{shared/}).
#' @param config_path Character or NULL. Path to the sourced config file, stored
#'   in \code{base_snapshot$config_path} and read line-by-line into
#'   \code{base_snapshot$config_source}. If NULL both are \code{NULL}.
#' @param verbose Integer. 0 = silent (default); >= 1 prints progress, matching
#'   the calling scripts' \code{cat()} chatter.
#'
#' @return A named list (the "bundle"):
#'   \describe{
#'     \item{null_model}{Fitted \code{GLOWr::fit_null_model} object (GDS order).}
#'     \item{pi_models}{Loaded PI models (\code{GLOWr::load_PI_models}).}
#'     \item{pi_features}{Resolved character vector of PI feature names.}
#'     \item{b_func}{The closed-form B function, or NULL.}
#'     \item{b_model}{The loaded trained B model, or NULL.}
#'     \item{filter_spec}{The config's \code{filter_spec} (echoed back).}
#'     \item{sample_ids}{Aligned sample IDs in GDS order.}
#'     \item{checksums}{Named list of MD5 sums for every input file.}
#'     \item{base_snapshot}{Named list of the base config_used fields shared by
#'       all modules (paths, resolved filter/feature settings, trait, provenance
#'       git_sha/session_info/prepared_at). The caller merges these with its
#'       module-specific snapshot fields.}
#'     \item{gds_files}{Named (by chr) character vector of per-chr GDS paths.}
#'   }
#'
#' @examples
#' \dontrun{
#' env <- new.env(); source("runs/myrun/config.R", local = env)
#' bundle <- prepare_scan_run(env, output_dir = "runs/myrun/outputs",
#'                            chroms = 22L, verbose = 1)
#' str(bundle$base_snapshot)
#' }
#'
#' @seealso \code{GLOWr::fit_null_model}, \code{GLOWr::load_PI_models},
#'   \code{GLOWr::load_B_model}
#' @export
prepare_scan_run <- function(config, output_dir, chroms = 1:22,
                             make_dirs = character(0),
                             config_path = NULL, verbose = 0) {

  # ---- Read config fields (env or list, NULL if unset) ----
  gds_dir         <- .psr_get(config, "gds_dir")
  gds_pattern     <- .psr_get(config, "gds_pattern")
  pheno_path      <- .psr_get(config, "pheno_path")
  pi_model_dir    <- .psr_get(config, "pi_model_dir")
  filter_spec     <- .psr_get(config, "filter_spec")
  ld_threshold    <- .psr_get(config, "ld_threshold")
  mac_threshold   <- .psr_get(config, "mac_threshold")
  collapse_method <- .psr_get(config, "collapse_method")
  pi_features     <- .psr_get(config, "pi_features")
  b_func          <- .psr_get(config, "b_func")
  b_model_path    <- .psr_get(config, "b_model_path")
  # `symlinked_shared_root` decides WHERE outputs/ lives (resolved into 3 states in the
  # Output-dirs block). Distinguish "absent" from an explicit NULL: both read back as NULL,
  # so check key existence too.
  ssr_set               <- .psr_key_exists(config, "symlinked_shared_root")
  symlinked_shared_root <- if (ssr_set) .psr_get(config, "symlinked_shared_root") else NULL

  chroms <- as.integer(chroms)
  stopifnot(length(chroms) >= 1L, all(chroms >= 1L & chroms <= 22L))

  # ---- Output dirs ----
  # `symlinked_shared_root` (config) decides WHERE outputs/ physically lives; the repo always
  # shows outputs/ in place. Three states:
  #   absent  -> born-shared into the shared-storage tooling's DEFAULT root (back-compat:
  #              a config that never mentions the field is shared, not silently local).
  #   NULL    -> OFF: outputs/ is a plain local directory under the repo.
  #   "<dir>" -> born-shared into <dir>/<repo-relative-path>.
  # Born-shared makes output_dir a symlink into shared storage BEFORE the dir.create below
  # populates it (so outputs land off /home and are registered for clone reconstruction). It is
  # a portable no-op outside a repo that ships tools/shared-storage/ / when the volume is
  # unavailable -> the dir.create then makes a plain local dir. See ensure_shared_output_dir
  born_shared <- if (ssr_set) !is.null(symlinked_shared_root) else TRUE
  if (born_shared)
    ensure_shared_output_dir(output_dir, share_root = symlinked_shared_root, verbose = verbose)
  shared_dir <- file.path(output_dir, "shared")
  for (d in c(shared_dir, file.path(output_dir, make_dirs))) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }

  # ---- B-source dispatch (mutually exclusive) ----
  # Either b_func (closed form, no training) or b_model_path (trained B model)
  # is non-NULL; never both. Matches GLOWr::prepare_glow_input's contract.
  b_func_set       <- !is.null(b_func)
  b_model_path_set <- !is.null(b_model_path)
  if (b_func_set == b_model_path_set) {
    stop("Config error: exactly one of `b_func` and `b_model_path` must be ",
         "non-NULL (got b_func=", b_func_set, ", b_model_path=",
         b_model_path_set, ").")
  }
  if (b_func_set) {
    stopifnot(is.function(b_func))
    # Sanity: produces a single finite numeric for a single MAF.
    probe <- b_func(0.1)
    if (!is.numeric(probe) || length(probe) != 1L || !is.finite(probe)) {
      stop("b_func(0.1) did not return a finite scalar numeric.")
    }
    if (verbose >= 1L)
      cat(sprintf("B source: closed-form b_func; b_func(0.1) = %.6f\n", probe))
  }
  if (b_model_path_set) {
    stopifnot(file.exists(b_model_path))
    if (verbose >= 1L)
      cat(sprintf("B source: trained model at %s\n", b_model_path))
  }

  # ---- Validate inputs ----
  stopifnot(dir.exists(gds_dir))
  stopifnot(file.exists(pheno_path))
  stopifnot(dir.exists(pi_model_dir))

  # Verify per-chr GDS files exist for the requested chromosomes.
  gds_files <- vapply(chroms, function(chr) {
    file.path(gds_dir, gsub("\\{chr\\}", chr, gds_pattern))
  }, character(1))
  names(gds_files) <- as.character(chroms)
  missing_gds <- gds_files[!file.exists(gds_files)]
  if (length(missing_gds) > 0L) {
    stop("Missing per-chr GDS files:\n  ", paste(missing_gds, collapse = "\n  "))
  }

  # ---- QC-node check: missing QC values (a VCF FILTER of ".") count as not passing ----
  # GLOWr >= 0.2.1 excludes a variant whose QC value is NA instead of stopping the
  # region. Count those variants up front, so that a run built from a VCF never
  # loses them silently: a warning with the per-chromosome counts, and a stop when
  # a chromosome would contribute no variant at all.
  if (!is.null(filter_spec)) {
    .psr_check_qc_node(gds_files,
                       qc_label = .psr_null_or(filter_spec$qc_label, "annotation/filter"),
                       verbose = verbose)
  }

  # ---- Load prepared phenotype/covariate data ----
  prep <- readRDS(pheno_path)
  if (verbose >= 1L)
    cat(sprintf("Pheno: n=%d (%d cases, %d controls), %d covariates\n",
                length(prep$Y), sum(prep$Y == 1L), sum(prep$Y == 0L),
                ncol(prep$X)))

  # ---- Sample alignment using the first requested chr's GDS (set is chr-independent) ----
  gds <- SeqArray::seqOpen(gds_files[1L], readonly = TRUE)
  gds_sample_ids <- SeqArray::seqGetData(gds, "sample.id")
  SeqArray::seqClose(gds)

  common_ids <- intersect(prep$sample_id, gds_sample_ids)
  if (length(common_ids) == 0L) stop("No overlapping samples between pheno and GDS")
  gds_common_ids <- gds_sample_ids[gds_sample_ids %in% common_ids]
  prep_idx <- match(gds_common_ids, prep$sample_id)

  Y_aligned  <- prep$Y[prep_idx]
  X_aligned  <- prep$X[prep_idx, , drop = FALSE]
  sample_ids <- prep$sample_id[prep_idx]
  stopifnot(all(sample_ids == gds_common_ids))
  if (verbose >= 1L)
    cat(sprintf("Aligned: %d common samples (in GDS order)\n", length(sample_ids)))

  # ---- Fit null model ----
  trait_resolved <- .psr_null_or(prep$trait, "binary")
  if (verbose >= 1L)
    cat("Fitting null model (trait =", trait_resolved, ")...\n")
  null_model <- GLOWr::fit_null_model(X = X_aligned, Y = Y_aligned,
                                      trait = trait_resolved,
                                      sample_id = sample_ids)
  if (verbose >= 1L)
    cat(sprintf("  n = %d, trait = %s\n", null_model$n, null_model$trait))

  # ---- Load PI (sanity) + resolve features; load B model if a path was given ----
  pi_models <- GLOWr::load_PI_models(pi_model_dir)
  if (verbose >= 1L)
    cat(sprintf("PI: %d %s models loaded\n",
                pi_models$n_models, pi_models$model_type))
  b_model <- NULL
  if (b_model_path_set) {
    b_model <- GLOWr::load_B_model(b_model_path)
    if (verbose >= 1L)
      cat(sprintf("B model loaded from %s\n", b_model_path))
  }
  pi_feats <- if (is.null(pi_features)) GLOWr:::.default_PI_features() else pi_features

  # ---- Checksums (input provenance) ----
  # MD5 fingerprints of every input file. Same content -> same hash; one byte
  # changed -> entirely different hash. A small passive audit trail (not actively
  # verified by Stage 2).
  input_checksums <- list(
    pheno_path     = tools::md5sum(pheno_path),
    b_model_path   = if (b_model_path_set) tools::md5sum(b_model_path) else NA_character_,
    pi_model_files = tools::md5sum(list.files(pi_model_dir, full.names = TRUE,
                                              pattern = "\\.rds$")),
    gds_files      = tools::md5sum(gds_files)
  )

  # ---- Base config_used snapshot fields (shared across all modules) ----
  base_snapshot <- list(
    config_path        = config_path,
    config_source      = if (!is.null(config_path)) readLines(config_path) else NULL,
    gds_dir            = gds_dir,
    gds_pattern        = gds_pattern,
    pheno_path         = pheno_path,
    pi_model_dir       = pi_model_dir,
    b_model_path       = if (b_model_path_set) b_model_path else NULL,
    b_func             = if (b_func_set) b_func else NULL,
    filter_spec        = filter_spec,
    ld_threshold       = ld_threshold,
    mac_threshold      = mac_threshold,
    collapse_method    = collapse_method,
    pi_features        = pi_feats,
    output_dir         = output_dir,
    symlinked_shared_root = if (!ssr_set) "(default; tools/shared-storage/shared-root.conf)"
                            else if (is.null(symlinked_shared_root)) "(none; local output_dir)"
                            else as.character(symlinked_shared_root),
    trait              = null_model$trait,
    glow_test_used_spa = (null_model$trait == "binary"),
    git_sha            = tryCatch(
      system("git rev-parse HEAD", intern = TRUE, ignore.stderr = TRUE)[1L],
      error = function(e) NA_character_),
    session_info       = utils::capture.output(utils::sessionInfo()),
    prepared_at        = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  )

  # ---- Persist the shared artifacts the spine owns ----
  saveRDS(null_model,      file.path(shared_dir, "01-null_model.rds"))
  saveRDS(input_checksums, file.path(shared_dir, "01-input_checksums.rds"))

  list(
    null_model    = null_model,
    pi_models     = pi_models,
    pi_features   = pi_feats,
    b_func        = if (b_func_set) b_func else NULL,
    b_model       = b_model,
    filter_spec   = filter_spec,
    sample_ids    = sample_ids,
    checksums     = input_checksums,
    base_snapshot = base_snapshot,
    gds_files     = gds_files
  )
}


#' Make a Run's Output Dir a Shared-Backed Symlink (Best-Effort, Portable)
#'
#' The integration point for born-shared outputs: a portable, guarded hook that
#' turns a run's
#' \code{output_dir} into a symlink into lab-shared storage \emph{before} the run
#' populates it, so the bytes live off \code{/home} and the leaf is registered for
#' clone reconstruction. \code{prepare_scan_run()} calls it for SNV-set runs; the
#' GLOWanalyses upstream templates call it directly when a run \emph{opts in} (config
#' field \code{symlinked_shared_root} non-NULL). This function does \emph{not} read
#' config — the caller decides whether to call it and passes the resolved
#' \code{share_root}.
#'
#' @details
#' If \code{output_dir} does not yet exist and lies inside a git repository that ships
#' \code{tools/shared-storage/scaffold_shared_output.sh}, it invokes that helper to
#' create \code{output_dir} as a symlink into shared storage and register it in the
#' manifest. With \code{share_root = NULL} the shared location is the helper's default
#' (\code{tools/shared-storage/shared-root.conf}, overridable by \code{GLOW_SHARE_ROOT});
#' with an explicit \code{share_root} the symlink points into
#' \code{<share_root>/<repo-relative-leaf>} and that exact target is recorded in the
#' manifest (so a clone reconstructs it regardless of any default).
#'
#' It is a deliberate \strong{no-op} (returning \code{FALSE}, after which the caller's own
#' \code{dir.create} makes a plain local dir) when: the dir already exists (real or
#' symlink); the parent does not exist yet; the path is not inside a git repo (e.g. a
#' \code{tempfile()} under \code{/tmp}, so unit tests are unaffected); the path is outside
#' the enclosing repo; the leaf is empty (the repo root); the helper is absent; or the
#' helper fails (e.g. the shared volume is unavailable). No shared path is hardcoded here,
#' so the package stays portable.
#'
#' @param output_dir Character. The run output directory (relative or absolute).
#' @param share_root Character or NULL. Explicit shared root to symlink into; \code{NULL}
#'   uses the shared-storage tooling's default root.
#' @param verbose Integer. \code{>= 1} prints a one-line note when a symlink is made.
#' @return Logical, invisibly: \code{TRUE} iff a shared-backed symlink was created.
#'
#' @examples
#' \dontrun{
#' # Opt-in born-sharing from a workflow template: symlink this run's outputs/ into a
#' # shared root before creating it (no-op outside a repo that ships the tooling).
#' ensure_shared_output_dir("runs/myrun/outputs", share_root = "/project/lab/shared")
#' }
#'
#' @seealso \code{\link{prepare_scan_run}}
#' @export
ensure_shared_output_dir <- function(output_dir, share_root = NULL, verbose = 0) {
  # Already present (real dir or existing symlink)? Nothing to do.
  # (Sys.readlink() returns NA for a non-existent path, and nzchar(NA) is TRUE, so
  # guard with !is.na() — otherwise a not-yet-created output_dir would falsely match.)
  rl0 <- Sys.readlink(output_dir)
  if (dir.exists(output_dir) || (!is.na(rl0) && nzchar(rl0))) return(invisible(FALSE))
  parent <- dirname(output_dir)
  if (!dir.exists(parent)) return(invisible(FALSE))  # nowhere to host the symlink yet
  # Locate the enclosing git repo and its shared-storage helper.
  repo <- tryCatch(
    suppressWarnings(system2("git", c("-C", parent, "rev-parse", "--show-toplevel"),
                             stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0))
  if (length(repo) != 1L || !nzchar(repo) || !dir.exists(repo)) return(invisible(FALSE))
  helper <- file.path(repo, "tools", "shared-storage", "scaffold_shared_output.sh")
  if (!file.exists(helper)) return(invisible(FALSE))
  # Require output_dir to be INSIDE the repo (guards against /tmp test dirs etc.).
  # Build the absolute path from the EXISTING parent: normalizePath(mustWork=FALSE) can
  # leave a path with a not-yet-existing leaf *relative* on Linux, which would fail the
  # prefix test. The parent exists (checked above), so it resolves to a real absolute path.
  parent_abs <- normalizePath(parent, mustWork = TRUE)
  repo_n     <- normalizePath(repo,   mustWork = TRUE)
  abs_out    <- file.path(parent_abs, basename(output_dir))
  prefix     <- paste0(repo_n, .Platform$file.sep)
  if (substr(abs_out, 1L, nchar(prefix)) != prefix) return(invisible(FALSE))
  leaf <- substring(abs_out, nchar(prefix) + 1L)
  if (!nzchar(leaf)) return(invisible(FALSE))  # never scaffold the repo root itself
  # A specific root -> pass it as an explicit --target so the bytes go there AND the manifest
  # records the exact target (clone reconstruction then resolves regardless of any default).
  # NULL -> let the scaffold helper use its default root (tools/shared-storage/shared-root.conf).
  args <- if (is.null(share_root)) c(leaf, "output")
          else c(leaf, "--target", file.path(share_root, leaf), "output")
  status <- tryCatch(
    system2(helper, args, stdout = FALSE, stderr = FALSE),
    error = function(e) 1L)
  rl1 <- Sys.readlink(output_dir)
  made <- identical(as.integer(status), 0L) && !is.na(rl1) && nzchar(rl1)
  if (made && verbose >= 1L)
    cat(sprintf("Output dir born-shared: %s -> %s\n", leaf, rl1))
  invisible(made)
}


#################### INTERNAL HELPER FUNCTIONS ####################

#' Read One Config Field from an Environment or List
#'
#' @param config Environment or named list of config values.
#' @param name Character. Field name.
#' @return The field value, or NULL if unset/absent.
#' @keywords internal
#' @noRd
.psr_get <- function(config, name) {
  if (is.environment(config)) {
    get0(name, envir = config, inherits = FALSE)
  } else if (is.list(config)) {
    config[[name]]
  } else {
    stop("`config` must be an environment or a named list.")
  }
}

#' Is a Config Field Set (Exists and Non-NULL)?
#'
#' @param config Environment or named list of config values.
#' @param name Character. Field name.
#' @return Logical scalar.
#' @keywords internal
#' @noRd
.psr_has <- function(config, name) {
  !is.null(.psr_get(config, name))
}

#' Does a Config Field KEY Exist (Even If Set to NULL)?
#'
#' Distinguishes an absent field from one explicitly set to \code{NULL} (both read
#' back as \code{NULL} via \code{.psr_get}); needed for the 3-state
#' \code{symlinked_shared_root} (absent vs NULL vs path).
#'
#' @param config Environment or named list of config values.
#' @param name Character. Field name.
#' @return Logical scalar.
#' @keywords internal
#' @noRd
.psr_key_exists <- function(config, name) {
  if (is.environment(config)) {
    exists(name, envir = config, inherits = FALSE)
  } else if (is.list(config)) {
    name %in% names(config)
  } else {
    stop("`config` must be an environment or a named list.")
  }
}

#' Return a if Non-NULL, Else b
#'
#' @param a,b Any values.
#' @return \code{a} when non-NULL, else \code{b}.
#' @keywords internal
#' @noRd
.psr_null_or <- function(a, b) if (is.null(a)) b else a



#' Count the variants whose QC value is missing, per chromosome
#'
#' \code{SeqArray::seqVCF2GDS()} stores a VCF \code{FILTER} of \code{"."}
#' (filters not applied) as \code{NA}. GLOWr's variant filter (0.2.1 and later)
#' treats such a value as not passing, so a run silently loses those variants
#' unless it is told. This helper reads the QC node of every per-chromosome GDS
#' before anything is fitted. A GDS without the node passes every variant (the
#' GLOWr fallback) and counts as 0. A chromosome whose QC values are all missing
#' would contribute no variant at all, so it stops the run. Any other missing
#' count is reported as a warning.
#'
#' @param gds_files Named (by chromosome) character vector of per-chr GDS paths.
#' @param qc_label Character GDS node path of the QC status (the filter spec's
#'   \code{qc_label}, default \code{"annotation/filter"}).
#' @param verbose Integer. At 1 or more, a one-line confirmation when nothing is missing.
#'
#' @return Invisibly, a \code{data.frame} with \code{chr}, \code{n_variants} and
#'   \code{n_missing_qc}, one row per chromosome.
#' @keywords internal
#' @noRd
.psr_check_qc_node <- function(gds_files, qc_label = "annotation/filter", verbose = 0) {
  counts <- lapply(names(gds_files), function(chr) {
    g <- SeqArray::seqOpen(gds_files[[chr]], readonly = TRUE)
    on.exit(SeqArray::seqClose(g), add = TRUE)
    n  <- length(SeqArray::seqGetData(g, "variant.id"))
    qc <- tryCatch(SeqArray::seqGetData(g, qc_label), error = function(e) NULL)
    data.frame(chr = chr, n_variants = n,
               n_missing_qc = if (is.null(qc)) 0L else sum(is.na(qc)),
               stringsAsFactors = FALSE)
  })
  counts <- do.call(rbind, counts)
  all_missing <- counts$n_variants > 0L & counts$n_missing_qc == counts$n_variants
  if (any(all_missing)) {
    stop("Every variant on chr ", paste(counts$chr[all_missing], collapse = ", "),
         " has a missing value in `", qc_label, "` (a VCF FILTER of \".\"), so no variant ",
         "would pass QC. If the VCF was not filtered, set FILTER to PASS when converting it ",
         "to GDS, or point `filter_spec$qc_label` at a node that carries the QC status.",
         call. = FALSE)
  }
  some <- counts$n_missing_qc > 0L
  if (any(some)) {
    warning(sprintf(paste0("%d variant(s) have a missing value in `%s` (a VCF FILTER of \".\") ",
                           "and are excluded from every set: %s"),
                    sum(counts$n_missing_qc), qc_label,
                    paste(sprintf("chr%s %d of %d", counts$chr[some], counts$n_missing_qc[some],
                                  counts$n_variants[some]), collapse = "; ")),
            call. = FALSE)
  } else if (verbose >= 1L) {
    cat("QC node `", qc_label, "`: no missing values\n", sep = "")
  }
  invisible(counts)
}
