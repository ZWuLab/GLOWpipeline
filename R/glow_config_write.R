# ==============================================================================
# File: glow_config_write.R
# Purpose: Write GLOW configurations as YAML that reads back to the same settings:
#   the resolved configuration of a run, and the base configuration generated from
#   the settings registry.
#
# EXPORTED FUNCTIONS:
#   - write_resolved_config()   The settings a run used -> <dir>/config-resolved.yml
#   - write_base_config()       The registry -> base_config.yml, every default with
#                               its description
#
# INTERNAL HELPERS:
#   - .gc_settings_lines()      Every setting of a registry as YAML lines
#   - .gc_setting_lines()       `name: value` lines of one setting
#   - .gc_render()              A typed value -> YAML, by its registered type
#   - .gc_form()                The form of a `one_of` setting that a value has
#   - .gc_num_text()            A double with the fewest digits that read back
#   - .gc_scalar_text()         Text as a YAML scalar, quoted when it has to be
#   - .gc_needs_escape()        Whether text holds a character YAML cannot hold as is
#
# Why not yaml::as.yaml(): it writes a logical as yes/no, loses digits of a double
# at its default precision, and writes a named vector as a plain list
# here serializes by the registered type, so that the reader of glow_config_read.R
# returns identical settings.

#################### EXPORTED MAIN FUNCTIONS ####################

#' Write the resolved configuration of a run
#'
#' Writes every setting of a resolved configuration, with the value in force, to
#' `config-resolved.yml` in a folder. The file is the record of what a run used,
#' and its settings read back, through \code{\link{read_glow_config}}, to the same
#' values.
#'
#' @param cfg A `glow_config` object from \code{\link{read_glow_config}}.
#' @param dir The folder to write into. It is created when it does not exist.
#' @param provenance A named list of what only the caller knows, for example
#'   `list(glowanalyses = "0.4.0", command = "glow run snv_set --config my.yml")`.
#'   Each element is a character, numeric or logical vector, and each name is used
#'   once. The names `files`, `registry`, `root`, `packages` and `written_at` are
#'   taken by the function.
#'
#' @return The path of the written file, invisibly.
#'
#' @details
#' The settings are written section by section in the order of the registry. A
#' logical is written as `true` or `false`, a number with the fewest significant
#' digits (15 to 17) that read back to the same double, a named vector as a
#' mapping, a function as its source text, and text that could be taken for a
#' number, a boolean or a null in quotes. A character that YAML cannot hold as it
#' is, such as a control character, is written as an escape.
#'
#' Seventeen significant digits identify a double, so the text of a number is exact.
#' Reading it back is R's conversion of text to a number. That conversion returns the
#' same double for a number of up to 15 significant digits on every build of R, which
#' covers a setting as a user writes it. For a number that needs 17 digits, or one
#' near the largest or the smallest double, it is exact where R computes in long
#' double (\code{capabilities("long.double")}) and can differ by one unit in the
#' last place elsewhere, for example on Apple silicon.
#'
#' Every value is checked against the registered type of its setting before
#' anything is written: the R type, whole numbers and the range of an integer,
#' infinity, the allowed words, and the fields of a map. A value that was changed
#' to something its setting does not take stops the function with the setting
#' named, so the function never writes a file that cannot be read back.
#'
#' After the settings comes one reserved key, `provenance`, which
#' \code{read_glow_config()} skips: the configuration files of the chain with their
#' checksums, the registry with its checksum, the root for relative paths, the
#' versions of GFisher, GLOWr and GLOWpipeline, the time of writing, and the
#' caller's entries. The provenance describes the configuration that was written.
#' A resolved file that is read and written again keeps its settings and gets a new
#' provenance, in which the resolved file itself is the last file of the chain.
#'
#' Computational complexity: linear in the size of the configuration.
#'
#' @examples
#' \dontrun{
#' cfg <- read_glow_config("cohort0203.yml",
#'                         base = "base-config/base_config.yml",
#'                         registry = "base-config/settings-registry.yml")
#' write_resolved_config(cfg, "runs/cohort0203/baseline/training",
#'                       provenance = list(command = "glow run training"))
#' }
#'
#' @seealso \code{\link{read_glow_config}}
#' @export
write_resolved_config <- function(cfg, dir, provenance = list()) {
  if (!inherits(cfg, "glow_config")) stop("`cfg` has to come from read_glow_config().", call. = FALSE)
  stopifnot(is.character(dir), length(dir) == 1L, is.list(provenance))
  taken <- c("files", "registry", "root", "packages", "written_at")
  if (length(provenance)) {
    nm <- names(provenance)
    if (is.null(nm) || any(!nzchar(nm))) stop("`provenance` has to be a named list.", call. = FALSE)
    if (anyDuplicated(nm)) stop("`provenance` uses the name `", nm[anyDuplicated(nm)], "` twice.", call. = FALSE)
    if (any(nm %in% taken)) stop("`provenance` cannot use the names ", paste(taken, collapse = ", "), ".", call. = FALSE)
    ok <- vapply(provenance, function(v) is.atomic(v) && !is.null(v) && is.null(names(v)) && !anyNA(v), logical(1))
    if (!all(ok)) stop("Each element of `provenance` has to be a plain vector without missing values.", call. = FALSE)
  }
  reg <- attr(cfg, "glow_registry")
  chain <- attr(cfg, "glow_chain")

  lines <- c("# Resolved GLOW configuration: every setting with the value this run used.",
             "# Written by GLOWpipeline::write_resolved_config(). Do not edit: it is a record.",
             "", .gc_settings_lines(reg, function(e) if (e$section == "top") cfg[[e$name]] else cfg[[e$section]][[e$name]],
                                    with_comments = FALSE))

  # The provenance block. `read_glow_config()` skips this key.
  pkg <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed")
  q <- function(x) .gc_scalar_text(x, force_quote = TRUE)
  prov <- c("", "provenance:", "  files:")
  for (i in seq_len(nrow(chain)))
    prov <- c(prov, paste0("    - file: ", q(chain$file[i])), paste0("      md5: ", q(chain$md5[i])))
  prov <- c(prov,
            "  registry:", paste0("    file: ", q(reg$file)), paste0("    md5: ", q(reg$md5)),
            paste0("  root: ", q(attr(cfg, "glow_root"))),
            "  packages:",
            paste0("    ", c("GFisher", "GLOWr", "GLOWpipeline"), ": ",
                   vapply(c(pkg("GFisher"), pkg("GLOWr"), pkg("GLOWpipeline")), q, character(1), USE.NAMES = FALSE)),
            paste0("  written_at: ", q(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))))
  for (i in seq_along(provenance)) {
    v <- provenance[[i]]; k <- paste0("  ", .gc_yaml_key(names(provenance)[i]), ":")
    txt <- vapply(as.character(v), q, character(1), USE.NAMES = FALSE)
    prov <- c(prov, if (!length(v)) paste(k, "[]") else if (length(v) == 1L) paste(k, txt)
                    else c(k, paste0("    - ", txt)))
  }

  if (!dir.exists(dir) && !dir.create(dir, recursive = TRUE))
    stop("Cannot create the folder ", dir, call. = FALSE)
  out <- file.path(dir, "config-resolved.yml")
  writeLines(enc2utf8(c(lines, prov)), out, useBytes = TRUE)
  invisible(out)
}

#' Write the base configuration from the settings registry
#'
#' Generates `base_config.yml`: every registered setting with its default, under
#' its description as comment lines, section by section. The file is what a user
#' reads to learn the settings and what \code{\link{read_glow_config}} reads for
#' the defaults, so the two cannot differ.
#'
#' @param registry Path of the settings registry (`settings-registry.yml`).
#' @param file Path of the file to write, or `NULL` to write nothing.
#'
#' @return The lines of the base configuration, invisibly. With `file = NULL` the
#'   caller can compare them with a committed file, to find out whether that file
#'   is stale.
#'
#' @details
#' The output is a function of the registry alone: it holds no date and no path of
#' the machine, so writing it twice gives the same bytes.
#'
#' Assumptions: the registry is valid, which \code{read_glow_config()} also checks.
#'
#' Computational complexity: linear in the number of settings.
#'
#' @examples
#' \dontrun{
#' write_base_config("base-config/settings-registry.yml", "base-config/base_config.yml")
#' }
#'
#' @seealso \code{\link{read_glow_config}}
#' @export
write_base_config <- function(registry, file = NULL) {
  stopifnot(is.character(registry), length(registry) == 1L,
            is.null(file) || (is.character(file) && length(file) == 1L))
  reg <- .gc_read_registry(registry)
  lines <- c(
    "# ==============================================================================",
    "# base_config.yml - every GLOWanalyses setting with its default",
    "#",
    "# GENERATED from settings-registry.yml by GLOWpipeline::write_base_config().",
    "# Do not edit this file: change the registry and regenerate it.",
    "#",
    "# A configuration of your own holds only the settings that differ. It builds on",
    "# this file unless its `extends` line names another configuration to build on.",
    "# ==============================================================================",
    "", .gc_settings_lines(reg, function(e) e$default, with_comments = TRUE))
  if (!is.null(file)) writeLines(enc2utf8(lines), file, useBytes = TRUE)
  invisible(lines)
}

#################### INTERNAL HELPER FUNCTIONS ####################

#' The settings of a registry as YAML lines, top level first, then the sections
#'
#' Every value is checked against its registered type before it is rendered, so a
#' value outside its type stops the writer with the setting named.
#' @param reg A registry from `.gc_read_registry()`.
#' @param get A function of an entry that returns the value to write.
#' @param with_comments Whether each setting gets its description as comment lines.
#' @return A character vector of lines.
#' @keywords internal
#' @noRd
.gc_settings_lines <- function(reg, get, with_comments) {
  lines <- character()
  for (sec in reg$sections) {
    es <- Filter(function(e) e$section == sec, reg$entries)
    if (!length(es)) next
    indent <- if (sec == "top") 0L else 2L
    if (sec != "top") {
      if (with_comments) lines <- c(lines, "", paste0("# ---- ", sec, " ", strrep("-", max(3L, 70L - nchar(sec)))))
      else if (length(lines)) lines <- c(lines, "")
      lines <- c(lines, paste0(sec, ":"))
    }
    for (e in es) {
      value <- get(e)
      msg <- .gc_valid(value, e$spec)
      if (nzchar(msg)) stop("`", e$key, "` ", msg, ", so it cannot be written.", call. = FALSE)
      if (with_comments) {
        pad <- strrep(" ", indent)
        note <- c(e$description,
                  if (!is.null(e$spec$choices)) paste0("Allowed: ", paste(e$spec$choices, collapse = ", "), "."),
                  if (isTRUE(e$derived)) "Null: found through the layout when the run needs it.")
        lines <- c(lines, if (length(lines) && nzchar(lines[length(lines)]) && !grepl(":$", lines[length(lines)])) "",
                   unlist(lapply(note, function(t) strwrap(t, width = 94L - indent, prefix = paste0(pad, "# ")))))
      }
      lines <- c(lines, .gc_setting_lines(e$name, value, e$spec, indent, e$key))
    }
  }
  lines
}

#' `key: value` lines of one setting at an indentation
#' @param key The key.
#' @param value The typed value.
#' @param spec Its type spec.
#' @param indent The number of blanks before the key.
#' @param where The qualified name of the setting, for messages.
#' @return A character vector of lines.
#' @keywords internal
#' @noRd
.gc_setting_lines <- function(key, value, spec, indent, where) {
  pad <- strrep(" ", indent)
  r <- .gc_render(value, spec, where)
  k <- paste0(pad, .gc_yaml_key(key), ":")
  switch(r$kind,
    inline  = paste(k, r$text),
    literal = c(paste(k, r$header), paste0(pad, "  ", r$lines)),
    block   = c(k, paste0(pad, "  ", r$lines)))
}

#' Render a typed value as YAML, by its registered type
#'
#' The value has passed `.gc_valid()` for the spec, so the function renders and
#' does not check again.
#' @param value The typed value.
#' @param spec Its type spec.
#' @param where The qualified name of the setting, for messages.
#' @return `list(kind = "inline", text)` for a value that fits on the key's line,
#'   `list(kind = "block", lines)` for mapping or sequence lines to put under the
#'   key, or `list(kind = "literal", header, lines)` for a block scalar.
#' @keywords internal
#' @noRd
.gc_render <- function(value, spec, where) {
  inline <- function(t) list(kind = "inline", text = t)
  if (is.null(value)) return(inline("~"))
  spec <- .gc_form(value, spec, where)
  # a list of scalars: on the key's line when it is short, else one item per line
  flow <- function(texts) {
    if (!length(texts)) return(inline("[]"))
    one <- paste0("[", paste(texts, collapse = ", "), "]")
    if (nchar(one) <= 80L) inline(one) else list(kind = "block", lines = paste0("- ", texts))
  }
  switch(spec$type,
    logical = inline(if (value) "true" else "false"),
    integer = , number = inline(.gc_num_text(value)),
    string = , path = , file_name = {
      if (.gc_block_ok(value)) {
        body <- strsplit(sub("\n$", "", value), "\n", fixed = TRUE)[[1L]]
        list(kind = "literal", header = if (endsWith(value, "\n")) "|" else "|-", lines = body)
      } else inline(.gc_scalar_text(value))
    },
    integer_list = {
      v <- as.numeric(value)
      # Three or more consecutive integers are written as the range a-b, the form the
      # reader takes: both ends finite and not negative.
      if (length(v) >= 3L && all(is.finite(v)) && v[1L] >= 0 && all(diff(v) == 1))
        inline(paste0(.gc_num_text(v[1L]), "-", .gc_num_text(v[length(v)])))
      else flow(vapply(v, .gc_num_text, character(1)))
    },
    number_list = flow(vapply(as.numeric(value), .gc_num_text, character(1))),
    string_list = flow(vapply(value, .gc_scalar_text, character(1), USE.NAMES = FALSE)),
    named_vector = {
      if (!length(value)) return(inline("{}"))
      txt <- vapply(seq_along(value), function(i) {
        v <- value[[i]]
        switch(spec$values$type, logical = if (v) "true" else "false",
               integer = , number = .gc_num_text(v), string = .gc_scalar_text(v)) }, character(1))
      list(kind = "block", lines = paste0(vapply(names(value), .gc_yaml_key, character(1)), ": ", txt))
    },
    map = {
      if (!length(value)) return(inline("{}"))
      # a map with declared fields is written in the declared order, the order the reader returns
      keys <- if (!is.null(spec$fields)) intersect(names(spec$fields), names(value)) else names(value)
      lines <- character()
      for (k in keys) {
        fs <- if (!is.null(spec$fields)) spec$fields[[k]] else spec$values
        lines <- c(lines, .gc_setting_lines(k, value[[k]], fs, 0L, paste0(where, ".", k)))
      }
      list(kind = "block", lines = lines)
    },
    sequence = {
      if (!length(value)) return(inline("[]"))
      lines <- character()
      for (i in seq_along(value)) {
        r <- .gc_render(value[[i]], spec$items, sprintf("%s[%d]", where, i))
        lines <- c(lines, switch(r$kind,
          inline  = paste0("- ", r$text),
          literal = c(paste("-", r$header), paste0("  ", r$lines)),
          block   = c(paste0("- ", r$lines[1L]), if (length(r$lines) > 1L) paste0("  ", r$lines[-1L]))))
      }
      list(kind = "block", lines = lines)
    },
    r_function = {
      src <- attr(value, "glow_source")
      if (is.null(src)) src <- paste(deparse(value), collapse = "\n")
      if (.gc_block_ok(src)) list(kind = "literal", header = "|-", lines = strsplit(src, "\n", fixed = TRUE)[[1L]])
      else inline(.gc_scalar_text(src, force_quote = TRUE))
    })
}

#' The form of a `one_of` setting that a value has
#'
#' For a spec that is not `one_of` the spec itself is returned. For a `one_of`
#' spec, the first form under which the value is fully valid, by the same check the
#' writer applies to every value, so two forms that are both maps are told apart by
#' their fields, and an empty list from an empty map by its names.
#' @param value The typed value.
#' @param spec Its type spec.
#' @param where The qualified name of the setting, for the message.
#' @return The spec of the form.
#' @keywords internal
#' @noRd
.gc_form <- function(value, spec, where) {
  if (!identical(spec$type, "one_of")) return(spec)
  for (alt in spec$one_of) if (!nzchar(.gc_valid(value, alt, exact = TRUE))) return(.gc_form(value, alt, where))
  stop("`", where, "` holds a value that fits none of the forms of its registered type.", call. = FALSE)
}

# ---- Scalars as YAML text -------------------------------------------------------

#' A double as text, with the fewest significant digits that read back to it
#'
#' Tries 15, 16 and 17 significant digits. Seventeen always identify an IEEE
#' double; fewer give `0.1` in place of `0.10000000000000001`. A whole number below
#' 1e15 is written without an exponent, and infinity as `.inf`.
#' @param x One finite or infinite double.
#' @return The text.
#' @keywords internal
#' @noRd
.gc_num_text <- function(x) {
  if (is.infinite(x)) return(if (x > 0) ".inf" else "-.inf")
  if (x == round(x) && abs(x) < 1e15) return(sprintf("%.0f", x))
  for (d in 15:17) {
    t <- sprintf("%.*g", d, x)
    if (identical(as.numeric(t), as.numeric(x))) break
  }
  t
}

#' Whether text may stand unquoted
#'
#' It starts with a letter, `/` or `_`, holds only characters that mean nothing to
#' YAML, and cannot be taken for a boolean or a null. Text that starts with a digit,
#' a sign or a dot is always quoted, which covers every number and the parser's
#' reserved words such as `.inf` and `.na.character`.
#' @param x One string.
#' @return `TRUE` or `FALSE`.
#' @keywords internal
#' @noRd
.gc_plain_ok <- function(x) {
  grepl("^[A-Za-z/_][A-Za-z0-9/._+-]*$", x) &&
    !tolower(x) %in% c("true", "false", "yes", "no", "on", "off", "y", "n", "null")
}

#' Whether text holds a character that YAML cannot hold as it is
#'
#' The control characters (except tab, line feed and carriage return, which have
#' escapes of their own), the C1 controls, among them NEL, which YAML reads as a
#' line break, the line and paragraph separators, and three code points YAML
#' forbids. Such text is written in double quotes with the character escaped.
#' @param x One string.
#' @return `TRUE` or `FALSE`.
#' @keywords internal
#' @noRd
.gc_needs_escape <- function(x) {
  cp <- utf8ToInt(enc2utf8(x))
  if (anyNA(cp)) return(TRUE)
  any((cp < 32L & !cp %in% c(9L, 10L, 13L)) | (cp >= 127L & cp <= 159L) |
        cp %in% c(8232L, 8233L, 65279L, 65534L, 65535L))
}

#' Text as a YAML scalar on one line
#'
#' Unquoted when that is safe, else in double quotes with backslash escapes.
#' @param x One string.
#' @param force_quote Whether to quote text that could stand unquoted.
#' @return The scalar as text.
#' @keywords internal
#' @noRd
.gc_scalar_text <- function(x, force_quote = FALSE) {
  x <- enc2utf8(x)
  if (!force_quote && .gc_plain_ok(x)) return(x)
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub("\"", "\\\"", x, fixed = TRUE)
  x <- gsub("\n", "\\n", x, fixed = TRUE)
  x <- gsub("\t", "\\t", x, fixed = TRUE)
  x <- gsub("\r", "\\r", x, fixed = TRUE)
  if (.gc_needs_escape(x)) {
    cp <- utf8ToInt(x)
    bad <- (cp < 32L) | (cp >= 127L & cp <= 159L) | cp %in% c(8232L, 8233L, 65279L, 65534L, 65535L)
    ch <- vapply(cp, intToUtf8, character(1))
    ch[bad] <- ifelse(cp[bad] < 256L, sprintf("\\x%02X", cp[bad]), sprintf("\\u%04X", cp[bad]))
    x <- paste(ch, collapse = "")
  }
  paste0("\"", x, "\"")
}

# A key of a mapping, quoted when it has to be.
.gc_yaml_key <- function(k) .gc_scalar_text(k)

#' Whether multi-line text can be written as a block scalar
#'
#' Several lines, none with a carriage return, a tab or a character that needs an
#' escape, the first not starting with a blank, no line ending in a blank (a block
#' scalar would keep it, an editor might not), and at most one final newline.
#' @param x One string.
#' @return `TRUE` or `FALSE`.
#' @keywords internal
#' @noRd
.gc_block_ok <- function(x) {
  grepl("\n", x, fixed = TRUE) && !grepl("[\r\t]", x) && !grepl("^[ \n]", x) &&
    !grepl(" (\n|$)", x) && !grepl("\n\n$", x) && !.gc_needs_escape(x)
}
