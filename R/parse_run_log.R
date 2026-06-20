# ==============================================================================
# Run-Log Parsers
# ==============================================================================
#
# Parse a per-chr run log (gene-based scan) or a per-chunk run log
# (sliding-window scan) into a one-row summary data.frame plus an errors
# data.frame, for the post-analysis diagnostics stages. The errors data.frame
# uses a shared `label` column so write_error_drilldowns() consumes both
# identically.
#
# EXPORTED FUNCTIONS:
#   - parse_chr_log()     Parse one log_chrN.txt -> $summary + $errors
#   - parse_chunk_log()   Parse one log_chunk_NNNN.txt -> $summary + $errors


#################### EXPORTED MAIN FUNCTIONS ####################

#' Parse a Per-Chromosome Run Log
#'
#' Parses one \code{log_chrN.txt} (written by the per-chr GLOW runner) into a
#' one-row summary data.frame of run statistics plus a sub-data.frame of error
#' rows for regions that errored.
#'
#' @details
#' The chromosome id is parsed from the filename. The gene-list count regex
#' accepts both \code{"(\\d+) coding genes"} and \code{"(\\d+) genes"}. Each
#' error line's optional leading word (\code{"STAAR ERROR in ..."} vs the bare
#' \code{"ERROR in ..."} main-pipeline catch-all) is parsed into the
#' \code{source} column; an absent prefix is labelled \code{"GLOW"}.
#'
#' @param log_path Character. Path to a \code{log_chrN.txt} file.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{$summary}}{one-row data.frame: \code{chr}, \code{n_samples},
#'       \code{trait}, \code{gene_list_n}, \code{filter}, \code{attempted},
#'       \code{results}, \code{skipped}, \code{errors}, \code{runtime_sec}.}
#'     \item{\code{$errors}}{data.frame with columns \code{chr}, \code{label},
#'       \code{source}, \code{message} (0 rows if no errors).}
#'   }
#'
#' @export
parse_chr_log <- function(log_path) {
  if (!file.exists(log_path)) stop("Log file not found: ", log_path)
  lines <- readLines(log_path)

  # Chromosome ID is encoded in the filename
  chr <- as.integer(sub(".*log_chr(\\d+)\\.txt$", "\\1", log_path))

  #This function is for pulling single values from the log by regex pattern. It returns the first match group of the first line that matches the pattern, coerced to integer or numeric if requested.
  pull <- function(pattern, default = NA_character_, integer = FALSE,
                   numeric = FALSE) {
    h <- grep(pattern, lines, value = TRUE)
    if (length(h) == 0L) return(default)
    v <- sub(paste0(".*", pattern, ".*"), "\\1", h[1L])
    if (integer) return(as.integer(v))
    if (numeric) return(as.numeric(v))
    v
  }

  n_samples   <- pull("n_samples = (\\d+)",           integer = TRUE)
  trait       <- pull("trait = ([a-zA-Z]+)")
  # Accept either "(\\d+) coding genes" (the coding-gene-list convention) or
  # "(\\d+) genes" (runs that source a non-coding-only gene list, e.g. a
  # custom gene list). The regex is anchored
  # before the genes count, so both wordings parse identically.
  gene_list_n <- pull("Gene list \\(chr\\d+\\): (\\d+) (?:coding )?genes",
                       integer = TRUE)
  filter_line <- pull("(Filter:.*)$")

  done_line <- grep("Done\\. genes_attempted=", lines, value = TRUE)
  if (length(done_line) == 0L) {
    # Run did not complete; surface NAs so callers can detect
    attempted <- NA_integer_; results <- NA_integer_
    skipped   <- NA_integer_; errors  <- NA_integer_
    runtime   <- NA_real_
  } else {
    attempted <- as.integer(sub(".*genes_attempted=(\\d+).*", "\\1",
                                 done_line[1L]))
    results   <- as.integer(sub(".*results=(\\d+).*",   "\\1", done_line[1L]))
    skipped   <- as.integer(sub(".*skipped=(\\d+).*",   "\\1", done_line[1L]))
    errors    <- as.integer(sub(".*errors=(\\d+).*",    "\\1", done_line[1L]))
    runtime   <- as.numeric(sub(".*runtime=([0-9.]+)s.*", "\\1", done_line[1L]))
  }

  err_lines <- grep("ERROR in", lines, value = TRUE)
  if (length(err_lines) == 0L) {
    err_df <- data.frame(chr        = integer(0L),
                         label      = character(0L),
                         source     = character(0L),
                         message    = character(0L),
                         stringsAsFactors = FALSE)
  } else {
    # Two log formats from 02-run-glow-chr.R: "[ts]   STAAR ERROR in GENE: msg"
    # (the STAAR sub-step) and "[ts]   ERROR in GENE: msg" (the main per-gene
    # catch-all). Capture the optional leading word as the error `source` so the
    # diagnostics report can say WHICH component failed; an absent prefix is the
    # main pipeline, labelled "GLOW". perl = TRUE for the non-capturing group.
    rx <- "(?:(\\w+) )?ERROR in ([^:]+): (.*)$"
    g  <- regmatches(err_lines, regexec(rx, err_lines, perl = TRUE))
    g  <- vapply(g, function(x) if (length(x) >= 4L) x[2:4]
                                 else rep(NA_character_, 3L),
                 character(3L))
    src <- trimws(g[1L, ])
    src[is.na(src) | src == ""] <- "GLOW"
    err_df <- data.frame(chr        = rep(chr, ncol(g)),
                         label      = trimws(g[2L, ]),
                         source     = src,
                         message    = trimws(g[3L, ]),
                         stringsAsFactors = FALSE)
  }

  summary <- data.frame(
    chr         = chr,
    n_samples   = n_samples,
    trait       = trait,
    gene_list_n = gene_list_n,
    filter      = filter_line,
    attempted   = attempted,
    results     = results,
    skipped     = skipped,
    errors      = errors,
    runtime_sec = runtime,
    stringsAsFactors = FALSE
  )

  list(summary = summary, errors = err_df)
}


#' Parse a Per-Chunk Sliding-Window Run Log
#'
#' Sliding-window analog of \code{\link{parse_chr_log}}. Parses one
#' \code{log_chunk_NNNN.txt} (written by the per-chunk scan runner) into a
#' one-row summary plus a per-window error data.frame.
#'
#' @details
#' Accepts both the older \code{"[ts] Starting chunk N: chrX [a, b], ..."}
#' header and the newer \code{"[ts] Starting chunk N [run R]: chrX [a, b], ..."}
#' header (\code{run_name} is \code{NA} when absent), so existing per-chunk logs
#' parse without regeneration. The window label always matches
#' \code{"chrN:start-end"}, which lets the parser locate the colon separating
#' label from message even though the label itself contains a colon.
#'
#' @param log_path Character. Path to a \code{log_chunk_NNNN.txt} file.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{$summary}}{one-row data.frame: \code{chunk_id},
#'       \code{run_name}, \code{chr}, \code{chunk_start}, \code{chunk_end},
#'       \code{est_n_windows}, \code{n_windows}, \code{n_ok},
#'       \code{n_skip_empty}, \code{n_staar_skip_single}, \code{n_error},
#'       \code{engine_runtime_secs}, \code{wall_secs}. NA-filled fields for an
#'       incomplete run.}
#'     \item{\code{$errors}}{data.frame with columns \code{chunk_id},
#'       \code{chr}, \code{label}, \code{message} (0 rows if no errors).}
#'   }
#'
#' @export
parse_chunk_log <- function(log_path) {
  if (!file.exists(log_path)) stop("Log file not found: ", log_path)
  lines <- readLines(log_path)
  chunk_id <- as.integer(sub(".*log_chunk_(\\d+)\\.txt$", "\\1", log_path))

  start_line <- grep("Starting chunk", lines, value = TRUE)
  start_line <- if (length(start_line) > 0L) start_line[1L] else NA_character_

  if (is.na(start_line)) {
    run_name <- NA_character_; chr <- NA_character_
    chunk_start <- NA_integer_; chunk_end <- NA_integer_
    est_n_windows <- NA_integer_
  } else {
    chr <- sub(".*: chr([^ ]+) \\[.*", "\\1", start_line)
    chunk_start <- as.integer(sub(".*\\[(\\d+),\\s*\\d+\\].*", "\\1",
                                  start_line))
    chunk_end   <- as.integer(sub(".*\\[\\d+,\\s*(\\d+)\\].*", "\\1",
                                  start_line))
    est_n_windows <- as.integer(sub(".*est_n_windows=(\\d+).*", "\\1",
                                    start_line))
    # Optional "[run R]" between "chunk N" and ":". perl = TRUE is required
    # for the non-capturing group (?:...); without it, regexec uses POSIX
    # ERE and treats (?: literally. An empty captured group is what
    # regexec returns when the optional group is missing; map to NA.
    m <- regmatches(start_line,
                    regexec("Starting chunk \\d+(?: \\[run ([^\\]]+)\\])?:",
                            start_line, perl = TRUE))[[1L]]
    run_name <- if (length(m) >= 2L && nzchar(m[2L])) m[2L] else NA_character_
  }

  done_line <- grep("Done chunk", lines, value = TRUE)
  done_line <- if (length(done_line) > 0L) done_line[1L] else NA_character_

  if (is.na(done_line)) {
    n_windows <- NA_integer_; n_ok <- NA_integer_
    n_skip_empty <- NA_integer_; n_staar_skip_single <- NA_integer_
    n_error <- NA_integer_
    engine_runtime_secs <- NA_real_; wall_secs <- NA_real_
  } else {
    n_windows <- as.integer(sub(".*n_windows=(\\d+).*", "\\1", done_line))
    n_ok      <- as.integer(sub(".*n_ok=(\\d+).*", "\\1", done_line))
    n_skip_empty <- as.integer(sub(".*n_skip_empty=(\\d+).*", "\\1",
                                   done_line))
    n_staar_skip_single <- as.integer(sub(".*n_staar_skip_single=(\\d+).*",
                                          "\\1", done_line))
    n_error <- as.integer(sub(".*n_error=(\\d+).*", "\\1", done_line))
    engine_runtime_secs <- as.numeric(sub(".*engine_runtime=([0-9.]+)s.*",
                                          "\\1", done_line))
    wall_secs <- as.numeric(sub(".*wall=([0-9.]+)s.*", "\\1", done_line))
  }

  # Errors: parse the "N window error(s):" header and the following N
  # "[ts]   <label>: <message>" lines. The label always matches
  # "chrN:start-end" (set by glow_scan_chunk()), which lets us locate the
  # colon separating label from message even though the label itself
  # contains a colon.
  err_df <- data.frame(chunk_id = integer(0L),
                       chr      = character(0L),
                       label    = character(0L),
                       message  = character(0L),
                       stringsAsFactors = FALSE)
  err_header <- grep("window error\\(s\\):", lines)
  if (length(err_header) > 0L && !is.na(n_error) && n_error > 0L) {
    hdr <- err_header[1L]
    tail_lines <- lines[seq(hdr + 1L,
                            min(hdr + n_error, length(lines)))]
    # perl = TRUE so \\s and \\d work; POSIX ERE doesn't support them.
    rx <- "^\\[[^\\]]+\\]\\s+(chr[^:\\s]+:\\d+-\\d+):\\s+(.*)$"
    m  <- regmatches(tail_lines, regexec(rx, tail_lines, perl = TRUE))
    g  <- vapply(m, function(x) if (length(x) >= 3L) x[2:3]
                                 else rep(NA_character_, 2L),
                 character(2L))
    labels   <- trimws(g[1L, ])
    messages <- trimws(g[2L, ])
    keep <- !is.na(labels) & nzchar(labels)
    if (any(keep)) {
      err_df <- data.frame(
        chunk_id = rep(chunk_id, sum(keep)),
        chr      = sub("^chr([^:]+):.*", "\\1", labels[keep]),
        label    = labels[keep],
        message  = messages[keep],
        stringsAsFactors = FALSE)
    }
  }

  summary <- data.frame(
    chunk_id            = chunk_id,
    run_name            = run_name,
    chr                 = chr,
    chunk_start         = chunk_start,
    chunk_end           = chunk_end,
    est_n_windows       = est_n_windows,
    n_windows           = n_windows,
    n_ok                = n_ok,
    n_skip_empty        = n_skip_empty,
    n_staar_skip_single = n_staar_skip_single,
    n_error             = n_error,
    engine_runtime_secs = engine_runtime_secs,
    wall_secs           = wall_secs,
    stringsAsFactors    = FALSE
  )
  list(summary = summary, errors = err_df)
}
