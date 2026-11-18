# ==============================================================================
# File: glow_config_registry.R
# Purpose: The internal machinery of the YAML configuration loader: reading a YAML
#   file without losing what the user wrote, reading the settings registry,
#   converting a value by the registered type of its setting, and checking a typed
#   value against its type.
#
# INTERNAL HELPERS (no exported function; the exported readers and writers are in
# glow_config_read.R and glow_config_write.R):
#   - .gc_parse_file()      a YAML file -> a tree of scalar tokens, lists and maps
#   - .gc_read_registry()   the settings registry -> entries with type specs and defaults
#   - .gc_spec()            a type spec written in YAML -> a type spec as an R list
#   - .gc_convert()         a token tree + a type spec -> a typed R value, problems collected
#   - .gc_valid()           a typed R value + a type spec -> "" or what is wrong with it
#   - .gc_nearest()         the nearest registered name, for "did you mean"
#
# yaml package converts plain scalars before a caller sees them, and for some the
# text is then lost (`010` arrives as 8, `3000000000` as a missing value, `no` as
# FALSE). The loader therefore replaces the handler of every scalar kind by one
# that returns the text as written together with its kind, and converts that text
# itself, by the type the registry gives the setting.

#################### INTERNAL HELPER FUNCTIONS ####################

# ---- Reading a YAML file as written -------------------------------------------

# The scalar kinds of the yaml package (2.3.12) whose handler the loader replaces.
# Each handler receives the scalar's text. `expr` is the kind of a value tagged
# `!expr`; it is caught here so that it can be refused, never evaluated.
.gc_scalar_kinds <- c("null", "bool#yes", "bool#no", "bool", "bool#na", "int", "int#hex",
                      "int#oct", "int#base60", "int#na", "float", "float#fix", "float#exp",
                      "float#base60", "float#inf", "float#neginf", "float#nan", "float#na",
                      "str", "str#na", "timestamp", "timestamp#iso8601", "timestamp#spaced",
                      "timestamp#ymd", "timestamp#hms", "expr")

#' The handlers that keep a YAML file as written
#'
#' One handler per scalar kind keeps the text and remembers the kind. A scalar
#' handler that receives anything but one string was called for a list or a map
#' that carries a scalar tag (`!!str [a, b]`), and it returns a token of kind `bad`.
#' The handlers of plain lists and maps mark them, so that a list or a map without
#' a mark is known to carry a tag of its own.
#' @return A named list of functions for `yaml::yaml.load(handlers = )`.
#' @keywords internal
#' @noRd
.gc_handlers <- function() {
  make <- function(kind) {
    force(kind)
    function(x) {
      if (!is.character(x) || length(x) != 1L || is.na(x))
        return(structure(list(text = NA_character_), kind = "bad", class = "glow_yaml_scalar"))
      structure(list(text = x), kind = kind, class = "glow_yaml_scalar")
    }
  }
  mark <- function(what) { force(what); function(x) { attr(x, "glow_node") <- what; x } }
  c(stats::setNames(lapply(.gc_scalar_kinds, make), .gc_scalar_kinds),
    list(seq = mark("seq"), map = mark("map")))
}

# The three node kinds of a parsed tree. A node from the parser carries the mark of
# its handler. A node built inside the loader carries none and is told apart by its
# names: an empty map is a list with names character(0), an empty list has no names.
.gc_is_scalar <- function(x) inherits(x, "glow_yaml_scalar")
.gc_is_map <- function(x) {
  if (!is.list(x) || .gc_is_scalar(x)) return(FALSE)
  m <- attr(x, "glow_node")
  if (is.null(m)) !is.null(names(x)) else identical(m, "map")
}
.gc_is_seq <- function(x) {
  if (!is.list(x) || .gc_is_scalar(x)) return(FALSE)
  m <- attr(x, "glow_node")
  if (is.null(m)) is.null(names(x)) else identical(m, "seq")
}
.gc_is_null <- function(x) .gc_is_scalar(x) && identical(attr(x, "kind"), "null")

# An empty map node, as the loader builds one.
.gc_empty_map <- function() structure(stats::setNames(list(), character(0)), glow_node = "map")

#' Read one YAML file into a tree of scalar tokens
#'
#' @param path File to read.
#' @param top Either "map" (a configuration) or "seq" (the registry): what the top
#'   level of the document has to be.
#' @return The parsed tree. Every scalar is a `glow_yaml_scalar` token holding the
#'   text as written; every node carries the attribute `glow_file`. An empty file,
#'   or one with comments only, gives an empty map (or an empty list).
#' @details Refused, with the file named: a file that is not a readable regular
#'   file or not valid UTF-8, a second YAML document, a value tagged as an R
#'   expression, any other tag on a value, a list or a map, a merge key that would
#'   replace a key, a duplicated or empty key, nesting deeper than 64 levels, and
#'   any warning of the parser. Nothing in the file is evaluated.
#' @keywords internal
#' @noRd
.gc_parse_file <- function(path, top = c("map", "seq")) {
  top <- match.arg(top)
  if (!file.exists(path)) stop("Configuration file not found: ", path, call. = FALSE)
  if (dir.exists(path)) stop("Cannot read ", path, ": it is a folder, not a file.", call. = FALSE)
  fail <- function(...) stop("Cannot read ", path, ": ", ..., call. = FALSE)

  # Everything that can go wrong while reading is reported with the file's name.
  tree <- tryCatch({
    # Read the bytes and drop a UTF-8 byte-order mark, so that the result does not
    # depend on the session's locale.
    raw <- readBin(path, what = "raw", n = file.size(path))
    md5 <- .gc_md5_bytes(raw)        # of the bytes that are parsed below, not of the file later
    if (length(raw) >= 3L && identical(as.integer(raw[1:3]), c(239L, 187L, 191L))) raw <- raw[-(1:3)]
    if (any(raw == as.raw(0L))) stop("the file holds a NUL byte; a configuration is a text file")
    text <- rawToChar(raw)
    if (!validUTF8(text)) stop("the file is not valid UTF-8")
    Encoding(text) <- "UTF-8"

    # One document only. The parser returns the first document and drops the rest
    # without a word. A line that starts with `---` is always a document marker in
    # YAML (the parser refuses one inside a quoted scalar), so the lines are split at
    # every line break YAML knows (LF, CR, CRLF, NEL, LS, PS) and the markers counted.
    # One marker may open the document, after directives (`%YAML ...`) and comments.
    lines <- strsplit(text, "\r\n|[\r\n\u0085\u2028\u2029]", perl = TRUE)[[1L]]
    content <- which(!grepl("^\\s*(#.*)?$", lines) & !grepl("^(%|---(\\s|$))", lines))
    marks   <- grep("^---(\\s|$)", lines)
    if (length(marks) > 1L || (length(marks) == 1L && length(content) && marks > content[1L]))
      stop("the file holds more than one YAML document (a line starting with `---` ",
           "after the first setting, or two such lines). A configuration is one document")

    withCallingHandlers(
      yaml::yaml.load(text, handlers = .gc_handlers(), eval.expr = FALSE, merge.warning = TRUE),
      # every warning of the parser is a failed or ambiguous reading: make it an error
      warning = function(w) stop(conditionMessage(w), call. = FALSE))
  }, error = function(e) fail(sub("\\.$", "", trimws(conditionMessage(e))), "."))

  # An empty document (no content at all) is no overrides. A document that consists
  # of a null, a scalar or the wrong kind of collection is not a configuration.
  if (is.null(tree)) return(structure(.gc_tag_file(if (top == "map") .gc_empty_map() else list(), path), glow_md5 = md5))
  ok_top <- if (top == "map") identical(attr(tree, "glow_node"), "map") else identical(attr(tree, "glow_node"), "seq")
  if (!ok_top)
    fail("the top level has to be a ", if (top == "map") "mapping (`name: value` lines)" else "list of entries", ".")
  tryCatch(.gc_check_tree(tree, path, where = "", depth = 0L),
           error = function(e) {
             msg <- conditionMessage(e)
             if (startsWith(msg, "Cannot read ")) stop(msg, call. = FALSE)
             fail(msg)                                        # for example a C stack error on absurd nesting
           })
  structure(.gc_tag_file(tree, path), glow_md5 = md5)
}

#' The MD5 checksum of bytes in memory
#'
#' @param raw A raw vector, the content of a file as it was read.
#' @return The checksum that `tools::md5sum()` gives for a file with these bytes.
#' @details A file is parsed from one read of its bytes, and its checksum has to be
#'   that of the same bytes: a checksum taken from the path afterwards can be that
#'   of an edited file. Base R hashes files only, so the bytes are written to a
#'   temporary file and that file is hashed.
#' @keywords internal
#' @noRd
.gc_md5_bytes <- function(raw) {
  tmp <- tempfile("glow-md5-")
  on.exit(unlink(tmp, expand = FALSE), add = TRUE)
  writeBin(raw, tmp)
  unname(tools::md5sum(tmp))
}

#' Check that a parsed tree holds only what a configuration may hold
#'
#' Every leaf has to be a scalar token of a known kind, every list and map has to
#' carry the mark of the plain list or map handler (a missing mark means a tag),
#' and no map may have an empty key.
#' @param node A node of the tree.
#' @param path The file, for the message.
#' @param where The position of the node, for the message.
#' @param depth Its depth.
#' @return `TRUE`, invisibly; an error otherwise.
#' @keywords internal
#' @noRd
.gc_check_tree <- function(node, path, where, depth) {
  at <- if (nzchar(where)) paste0("`", where, "`") else "the top level"
  fail <- function(...) stop("Cannot read ", path, ": ", ..., call. = FALSE)
  if (depth > 64L) fail(at, " is nested more than 64 levels deep.")
  if (.gc_is_scalar(node)) {
    kind <- attr(node, "kind")
    if (identical(kind, "expr"))
      fail(at, " is tagged as an R expression (`!expr`). A configuration holds values, ",
           "and the loader evaluates nothing.")
    if (identical(kind, "bad"))
      fail(at, " carries a YAML tag that does not fit its value. Write the value without a tag.")
    return(invisible(TRUE))
  }
  if (!is.list(node) || is.null(attr(node, "glow_node")))
    fail(at, " carries a YAML tag this loader does not read. Write the value without a tag.")
  nms <- names(node)
  if (identical(attr(node, "glow_node"), "map") && length(node) && (is.null(nms) || any(!nzchar(nms))))
    fail(at, " has an empty key. Every key of a mapping needs a name.")
  for (i in seq_along(node)) {
    child <- if (identical(attr(node, "glow_node"), "seq")) sprintf("%s[%d]", where, i)
             else if (nzchar(where)) paste0(where, ".", nms[i]) else nms[i]
    .gc_check_tree(node[[i]], path, child, depth + 1L)
  }
  invisible(TRUE)
}

#' Mark every node of a tree with the file it came from
#'
#' The mark is the attribute `glow_file`, which the messages after the merge use.
#' @param node A node of the tree.
#' @param path The file.
#' @return The node, marked, with its children marked.
#' @keywords internal
#' @noRd
.gc_tag_file <- function(node, path) {
  if (is.list(node) && !.gc_is_scalar(node)) {
    for (i in seq_along(node)) node[[i]] <- .gc_tag_file(node[[i]], path)
  }
  attr(node, "glow_file") <- path
  node
}

# ---- Type specs ---------------------------------------------------------------

.gc_types <- c("logical", "integer", "number", "string", "path", "file_name", "integer_list", "number_list",
               "string_list", "named_vector", "map", "sequence", "one_of", "r_function")
.gc_spec_keys <- c("type", "null_ok", "inf_ok", "choices", "fields", "required", "values", "items", "one_of")

#' The text of a registry field that has to be plain text
#' @param node The token.
#' @param what What it is, for the message.
#' @param path The registry.
#' @return The text; an error when the node is not a text value.
#' @keywords internal
#' @noRd
.gc_meta_text <- function(node, what, path) {
  if (!.gc_is_scalar(node) || .gc_is_null(node))
    stop("Registry ", path, ": ", what, " has to be a text value.", call. = FALSE)
  node$text
}

#' The value of a registry field that has to be true or false
#' @param node The token.
#' @param what What it is, for the message.
#' @param path The registry.
#' @return `TRUE` or `FALSE`; an error otherwise.
#' @keywords internal
#' @noRd
.gc_meta_flag <- function(node, what, path) {
  t <- if (.gc_is_scalar(node)) tolower(node$text) else ""
  if (!t %in% c("true", "false"))
    stop("Registry ", path, ": ", what, " has to be true or false.", call. = FALSE)
  t == "true"
}

#' Turn a type spec written in the registry into an R list
#'
#' A spec is a mapping with `type` and, by type, `null_ok`, `inf_ok`, `choices`
#' (allowed words of a string or of the items of a string_list), `fields` (a map
#' with declared keys), `required` (the fields a map cannot do without), `values`
#' (the type of every value of a free-key map or of a named_vector), `items` (the
#' type of every item of a sequence) and `one_of` (a list of specs, tried in order).
#' @param node The spec as parsed.
#' @param what What it is, for the message.
#' @param path The registry.
#' @param inside The type that contains this spec, or "".
#' @return The spec as an R list.
#' @details A null is allowed where the spec that receives the value says so. A
#'   form of a `one_of` and the values of a `named_vector` never receive a null
#'   (the `one_of` decides, and a named vector holds no null), so `null_ok` on them
#'   would promise something the loader does not do, and the registry is refused.
#' @keywords internal
#' @noRd
.gc_spec <- function(node, what, path, inside = "") {
  if (!.gc_is_map(node)) stop("Registry ", path, ": ", what, " has to be a mapping with `type`.", call. = FALSE)
  bad <- setdiff(names(node), .gc_spec_keys)
  if (length(bad))
    stop("Registry ", path, ": ", what, " has the unknown key(s) ", paste(bad, collapse = ", "), ".", call. = FALSE)
  if (is.null(node$type)) stop("Registry ", path, ": ", what, " has no `type`.", call. = FALSE)
  spec <- list(type = .gc_meta_text(node$type, paste0(what, ": type"), path))
  if (!spec$type %in% .gc_types)
    stop("Registry ", path, ": ", what, " has the unknown type `", spec$type, "`. Known types: ",
         paste(.gc_types, collapse = ", "), ".", call. = FALSE)
  if (!is.null(node$null_ok)) spec$null_ok <- .gc_meta_flag(node$null_ok, paste0(what, ": null_ok"), path)
  if (!is.null(node$inf_ok))  spec$inf_ok  <- .gc_meta_flag(node$inf_ok,  paste0(what, ": inf_ok"),  path)
  if (isTRUE(spec$null_ok) && inside %in% c("one_of", "named_vector"))
    stop("Registry ", path, ": ", what, " has `null_ok`, which has no effect ",
         if (inside == "one_of") "on a form of a `one_of`; put it on the `one_of` itself"
         else "on the values of a `named_vector`", ".", call. = FALSE)
  if (!is.null(node$choices)) {
    if (!.gc_is_seq(node$choices)) stop("Registry ", path, ": ", what, ": choices has to be a list.", call. = FALSE)
    spec$choices <- vapply(node$choices, .gc_meta_text, character(1), what = paste0(what, ": choices"), path = path)
  }
  sub <- function(n, key, inside = "") .gc_spec(n, paste0(what, ": ", key), path, inside)
  if (!is.null(node$fields)) {
    if (!.gc_is_map(node$fields)) stop("Registry ", path, ": ", what, ": fields has to be a mapping.", call. = FALSE)
    spec$fields <- stats::setNames(lapply(names(node$fields), function(k) sub(node$fields[[k]], paste0("fields.", k))),
                                   names(node$fields))
  }
  if (!is.null(node$required)) {
    if (!.gc_is_seq(node$required)) stop("Registry ", path, ": ", what, ": required has to be a list of field names.", call. = FALSE)
    spec$required <- vapply(node$required, .gc_meta_text, character(1), what = paste0(what, ": required"), path = path)
    if (is.null(spec$fields) || !all(spec$required %in% names(spec$fields)))
      stop("Registry ", path, ": ", what, ": every name in `required` has to be one of its `fields`.", call. = FALSE)
  }
  if (!is.null(node$values)) spec$values <- sub(node$values, "values", if (spec$type == "named_vector") "named_vector" else "")
  if (!is.null(node$items))  spec$items  <- sub(node$items, "items")
  if (!is.null(node$one_of)) {
    if (!.gc_is_seq(node$one_of)) stop("Registry ", path, ": ", what, ": one_of has to be a list of specs.", call. = FALSE)
    spec$one_of <- lapply(seq_along(node$one_of), function(i) sub(node$one_of[[i]], sprintf("one_of[%d]", i), "one_of"))
  }
  # what each type needs
  need <- switch(spec$type,
    named_vector = if (is.null(spec$values) || !spec$values$type %in% c("logical", "integer", "number", "string"))
                     "`values` with a scalar type (logical, integer, number or string)",
    map          = if (is.null(spec$fields) == is.null(spec$values)) "exactly one of `fields` and `values`",
    sequence     = if (is.null(spec$items)) "`items`",
    one_of       = if (!length(spec$one_of)) "a non-empty `one_of` list",
    NULL)
  if (!is.null(need)) stop("Registry ", path, ": ", what, " (type ", spec$type, ") needs ", need, ".", call. = FALSE)
  spec
}

# ---- The registry -------------------------------------------------------------

.gc_entry_keys <- c("name", "section", "default", "derived", "replace_whole", "was", "description",
                    .gc_spec_keys)

#' Read the settings registry
#'
#' @param path The registry, a YAML list with one entry per setting.
#' @return A list with `entries` (a named list keyed `section.name`, or `name` for a
#'   top-level setting; each entry has `name`, `section`, `key`, `spec`, `default`
#'   (the typed default), `default_node`, `derived`, `replace_whole`, `was`, `description`),
#'   `sections` (in the order of first appearance, "top" first), `file` and `md5`.
#' @details An entry's section is `top` for a setting written at the top level of a
#'   configuration. The registry is read with the same reader as a configuration,
#'   so a default is held to the rules a user's value is held to. The names
#'   `extends` and `provenance` are reserved and cannot be settings or sections.
#' @keywords internal
#' @noRd
.gc_read_registry <- function(path) {
  tree <- .gc_parse_file(path, top = "seq")
  if (!length(tree)) stop("Registry ", path, " holds no setting.", call. = FALSE)
  entries <- list()
  for (i in seq_along(tree)) {
    node <- tree[[i]]
    what <- sprintf("entry %d", i)
    if (!.gc_is_map(node)) stop("Registry ", path, ": ", what, " has to be a mapping.", call. = FALSE)
    bad <- setdiff(names(node), .gc_entry_keys)
    if (length(bad))
      stop("Registry ", path, ": ", what, " has the unknown key(s) ", paste(bad, collapse = ", "), ".", call. = FALSE)
    for (k in c("name", "section", "type", "default", "description"))
      if (!k %in% names(node)) stop("Registry ", path, ": ", what, " has no `", k, "`.", call. = FALSE)
    name    <- .gc_meta_text(node$name, paste0(what, ": name"), path)
    section <- .gc_meta_text(node$section, paste0(what, " (", name, "): section"), path)
    what    <- if (section == "top") name else paste0(section, ".", name)
    ok_id   <- function(x) grepl("^[A-Za-z][A-Za-z0-9_]*$", x)
    if (!ok_id(name) || !ok_id(section))
      stop("Registry ", path, ": `", what, "`: a name and a section are letters, digits and `_`.", call. = FALSE)
    if (name %in% c("extends", "provenance") || section %in% c("extends", "provenance"))
      stop("Registry ", path, ": `", what, "`: `extends` and `provenance` are reserved.", call. = FALSE)
    spec_node <- structure(node[intersect(names(node), .gc_spec_keys)], glow_node = "map")
    spec <- .gc_spec(spec_node, paste0("`", what, "`"), path)
    if (what %in% names(entries)) stop("Registry ", path, ": `", what, "` is defined twice.", call. = FALSE)
    entries[[what]] <- list(
      name = name, section = section, key = what, spec = spec, default_node = node$default,
      derived = if (is.null(node$derived)) FALSE else .gc_meta_flag(node$derived, paste0("`", what, "`: derived"), path),
      # a map with declared fields merges field by field over the chain, unless the entry says
      # that a user's value replaces the inherited one whole (fields that are alternatives)
      replace_whole = if (is.null(node$replace_whole)) FALSE
                      else .gc_meta_flag(node$replace_whole, paste0("`", what, "`: replace_whole"), path),
      was = if (is.null(node$was)) NULL else .gc_plain(node$was),
      description = gsub("\\s+", " ", trimws(.gc_meta_text(node$description, paste0("`", what, "`: description"), path))))
  }
  # a top-level setting and a section cannot share a name
  sections <- unique(vapply(entries, `[[`, character(1), "section"))
  top_names <- vapply(Filter(function(e) e$section == "top", entries), `[[`, character(1), "name")
  clash <- intersect(top_names, setdiff(sections, "top"))
  if (length(clash)) stop("Registry ", path, ": `", clash[1L], "` is both a top-level setting and a section.", call. = FALSE)
  # every default has to be valid under its own type
  acc <- new.env(); acc$problems <- character()
  for (k in names(entries)) entries[[k]]["default"] <- list(.gc_convert(entries[[k]]$default_node, entries[[k]]$spec, k, acc))
  if (length(acc$problems))
    stop("Registry ", path, " has ", length(acc$problems), " invalid default(s):\n",
         paste0("  - ", acc$problems, collapse = "\n"), call. = FALSE)
  list(entries = entries, sections = c("top", setdiff(sections, "top")),
       file = normalizePath(path), md5 = attr(tree, "glow_md5"))        # the checksum of the bytes that were parsed
}

#' A token tree as plain R, for registry fields that are not settings
#' @param node A node of a parsed tree.
#' @return `NULL` for a null, the text of a scalar, a character vector for a list of
#'   texts, and a list otherwise.
#' @keywords internal
#' @noRd
.gc_plain <- function(node) {
  if (.gc_is_null(node)) return(NULL)
  if (.gc_is_scalar(node)) return(node$text)
  out <- lapply(node, .gc_plain)
  attributes(out) <- list(names = names(node))
  if (.gc_is_seq(node) && all(vapply(out, function(v) is.character(v) && length(v) == 1L, logical(1))))
    return(unlist(out))
  out
}

# ---- Converting a value by its registered type ---------------------------------

# Numbers are written in decimal or scientific notation. The loader converts the
# text itself, so `010` is ten and never octal eight.
.gc_num_pattern <- "^[+-]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][+-]?[0-9]+)?$"
.gc_inf_pattern <- "^[+-]?\\.?(inf|Inf|INF)$"

#' Record one problem of a configuration
#' @param acc The environment whose `problems` vector collects the messages.
#' @param where The qualified setting.
#' @param node The offending node, which names its file.
#' @param what What is wrong.
#' @return `NULL`, invisibly, so that a converter can return the call.
#' @keywords internal
#' @noRd
.gc_problem <- function(acc, where, node, what) {
  file <- attr(node, "glow_file")
  written <- if (.gc_is_scalar(node) && !.gc_is_null(node)) sprintf(" (written: %s)", encodeString(node$text, quote = "`")) else ""
  # a message that ends in a question ("did you mean ...?") takes no full stop after it
  msg <- paste0(if (!is.null(file)) paste0(file, ": ") else "", "`", where, "` ", what, written, ".")
  acc$problems <- c(acc$problems, sub("?.", "?", msg, fixed = TRUE))
  invisible(NULL)
}

# What a node is, in the words of a message.
.gc_shape <- function(node) if (.gc_is_scalar(node)) "a single value" else if (.gc_is_map(node)) "a mapping" else "a list"

#' Whether a number written in decimal or scientific notation is whole
#'
#' Decided from the text, because a double cannot tell `1.0000000000000001` from 1.
#' The value is the digits of the mantissa times a power of ten, and it is whole
#' when no nonzero digit is left behind the decimal point.
#' @param t Text that matches `.gc_num_pattern`.
#' @return `TRUE` or `FALSE`.
#' @keywords internal
#' @noRd
.gc_text_is_whole <- function(t) {
  m <- regmatches(t, regexec("^[+-]?([0-9]*)\\.?([0-9]*)([eE]([+-]?[0-9]+))?$", t))[[1L]]
  ex <- if (nzchar(m[5L])) as.numeric(m[5L]) else 0
  frac <- sub("0+$", "", m[3L])                       # the digits behind the point, without trailing zeros
  if (nzchar(frac)) return(ex >= nchar(frac))         # the exponent has to move every one of them left of the point
  # no fractional digit: whole unless a negative exponent cuts into the nonzero integer digits
  int <- sub("^0+", "", m[2L])
  if (ex >= 0 || !nzchar(int)) return(TRUE)
  -ex <= nchar(int) - nchar(sub("0+$", "", int))      # as many trailing zeros as the exponent removes
}

#' The text of a scalar read as a number
#' @param node The token.
#' @param spec Its type spec (`inf_ok`).
#' @param where The setting.
#' @param acc The problem collector.
#' @return A double, or `NULL` with a problem recorded. A nonzero number too small
#'   for a double is refused, because it would arrive as zero.
#' @keywords internal
#' @noRd
.gc_number <- function(node, spec, where, acc) {
  if (!.gc_is_scalar(node)) return(.gc_problem(acc, where, node, paste0("expects a number, not ", .gc_shape(node))))
  t <- trimws(node$text)
  if (grepl(.gc_inf_pattern, t)) {
    if (isTRUE(spec$inf_ok)) return(if (startsWith(t, "-")) -Inf else Inf)
    return(.gc_problem(acc, where, node, "expects a finite number"))
  }
  if (!grepl(.gc_num_pattern, t))
    return(.gc_problem(acc, where, node, "expects a number in decimal or scientific notation, such as 10, 0.5 or 1e6"))
  v <- as.numeric(t)
  if (!is.finite(v)) return(.gc_problem(acc, where, node, "is too large to be read as a number"))
  if (v == 0 && grepl("[1-9]", sub("[eE].*$", "", t)))
    return(.gc_problem(acc, where, node, "is too small to be read as a number other than zero"))
  v
}

#' The text of a scalar read as an integer
#' @inheritParams .gc_number
#' @return An integer (or an infinity when `inf_ok`), or `NULL` with a problem recorded.
#' @keywords internal
#' @noRd
.gc_integer <- function(node, spec, where, acc) {
  if (.gc_is_scalar(node)) {
    t <- trimws(node$text)
    if (grepl(.gc_num_pattern, t) && !.gc_text_is_whole(t))
      return(.gc_problem(acc, where, node, "expects a whole number"))
  }
  v <- .gc_number(node, spec, where, acc)
  if (is.null(v) || is.infinite(v)) return(v)
  if (abs(v) > .Machine$integer.max) return(.gc_problem(acc, where, node, "is outside the range of an integer"))
  as.integer(v)
}

#' The text of a scalar read as true or false
#' @inheritParams .gc_number
#' @return `TRUE` or `FALSE`, or `NULL` with a problem recorded.
#' @keywords internal
#' @noRd
.gc_logical <- function(node, spec, where, acc) {
  if (!.gc_is_scalar(node)) return(.gc_problem(acc, where, node, paste0("expects true or false, not ", .gc_shape(node))))
  t <- tolower(trimws(node$text))
  if (!t %in% c("true", "false")) return(.gc_problem(acc, where, node, "expects true or false"))
  t == "true"
}

#' The text of a scalar, as written
#' @inheritParams .gc_number
#' @return The text, or `NULL` with a problem recorded (a list, a map, a null, or a
#'   word outside the `choices` of the spec).
#' @keywords internal
#' @noRd
.gc_string <- function(node, spec, where, acc) {
  if (!.gc_is_scalar(node)) return(.gc_problem(acc, where, node, paste0("expects a single text value, not ", .gc_shape(node))))
  if (.gc_is_null(node)) return(.gc_problem(acc, where, node, "holds a null where text is expected; quote the word if it is meant as text"))
  if (!is.null(spec$choices) && !node$text %in% spec$choices)
    return(.gc_problem(acc, where, node, paste0("expects one of: ", paste(spec$choices, collapse = ", "))))
  if (identical(spec$type, "file_name") && nzchar(.gc_file_name_problem(node$text)))
    return(.gc_problem(acc, where, node, .gc_file_name_problem(node$text)))
  node$text
}

#' Why a text is not the name of one file or folder
#'
#' A setting of the type `file_name` names a file or a folder inside a place that
#' the workflow decides, so it must not lead out of that place.
#'
#' @param text One string.
#' @return An empty string for a name, else the problem as a phrase.
#' @keywords internal
#' @noRd
.gc_file_name_problem <- function(text) {
  if (!nzchar(text) || grepl("[/\\\\]", text) || text %in% c(".", "..") || startsWith(text, "~"))
    return("expects the name of one file or folder, without a folder in front of it")
  ""
}

#' A list of scalars of one kind: a YAML list, or one item
#' @inheritParams .gc_number
#' @param one The converter of one item.
#' @param empty The value of an empty list.
#' @return A vector, or `NULL` with a problem recorded.
#' @keywords internal
#' @noRd
.gc_scalar_list <- function(node, spec, where, acc, one, empty) {
  if (.gc_is_map(node)) return(.gc_problem(acc, where, node, "expects a list, not a mapping"))
  items <- if (.gc_is_scalar(node)) list(node) else node
  out <- vector("list", length(items))
  for (i in seq_along(items)) {
    w <- if (.gc_is_scalar(node)) where else sprintf("%s[%d]", where, i)
    if (!.gc_is_scalar(items[[i]]) || .gc_is_null(items[[i]])) { .gc_problem(acc, w, items[[i]], "has to be a single value"); next }
    v <- one(items[[i]], spec, w, acc)
    if (!is.null(v)) out[[i]] <- v
  }
  if (any(vapply(out, is.null, logical(1)))) return(NULL)   # a problem was recorded
  if (!length(out)) return(empty)
  unlist(out)
}

#' Convert a parsed value by a type spec
#'
#' @param node A token tree (from `.gc_parse_file()`).
#' @param spec A type spec (from `.gc_spec()`).
#' @param where The qualified name of the setting, for messages.
#' @param acc An environment whose `problems` vector collects the messages.
#' @param partial Whether a map with declared fields may lack its `required` fields.
#'   True when one file of a chain is checked on its own: the fields it does not
#'   name come from the files before it. It applies to the value itself, not to the
#'   values inside it.
#' @return The typed R value: a logical, an integer (or `Inf`), a double, a string,
#'   a vector, a named vector, a list (a map with declared fields in the declared
#'   order, a map with free keys in the order written), an R function with the
#'   attribute `glow_source`, or `NULL`. On a problem the function records a
#'   message and returns `NULL`.
#' @keywords internal
#' @noRd
.gc_convert <- function(node, spec, where, acc, partial = FALSE) {
  if (.gc_is_null(node)) {
    if (isTRUE(spec$null_ok)) return(NULL)
    return(.gc_problem(acc, where, node, "is null, and this setting needs a value"))
  }
  switch(spec$type,
    logical = .gc_logical(node, spec, where, acc),
    integer = .gc_integer(node, spec, where, acc),
    number  = .gc_number(node, spec, where, acc),
    string  = ,
    file_name = ,
    path    = .gc_string(node, spec, where, acc),
    integer_list = {
      t <- if (.gc_is_scalar(node)) trimws(node$text) else ""
      if (grepl("^[0-9]+\\s*-\\s*[0-9]+$", t)) {                 # the range form a-b
        ab <- suppressWarnings(as.numeric(strsplit(t, "\\s*-\\s*")[[1L]]))
        if (any(ab > .Machine$integer.max) || ab[1L] > ab[2L])
          return(.gc_problem(acc, where, node, "is not a range a-b with a not above b"))
        return(seq.int(as.integer(ab[1L]), as.integer(ab[2L])))
      }
      if (grepl("^[0-9]+\\s*:\\s*[0-9]+$", t))
        return(.gc_problem(acc, where, node, "is an R range; write it as a-b, for example 1-22"))
      .gc_scalar_list(node, spec, where, acc, .gc_integer, integer(0))
    },
    number_list = .gc_scalar_list(node, spec, where, acc, .gc_number, numeric(0)),
    string_list = .gc_scalar_list(node, spec, where, acc, .gc_string, character(0)),
    named_vector = {
      if (!.gc_is_map(node)) return(.gc_problem(acc, where, node, paste0("expects a mapping from names to values, not ", .gc_shape(node))))
      one <- switch(spec$values$type, logical = .gc_logical, integer = .gc_integer, number = .gc_number, string = .gc_string)
      empty <- switch(spec$values$type, logical = logical(0), integer = integer(0), number = numeric(0), string = character(0))
      out <- vector("list", length(node))
      for (i in seq_along(node)) {
        w <- paste0(where, ".", names(node)[i])
        if (!.gc_is_scalar(node[[i]]) || .gc_is_null(node[[i]])) { .gc_problem(acc, w, node[[i]], "has to be a single value"); next }
        v <- one(node[[i]], spec$values, w, acc)
        if (!is.null(v)) out[[i]] <- v
      }
      if (any(vapply(out, is.null, logical(1)))) return(NULL)
      stats::setNames(if (length(out)) unlist(out) else empty, names(node))
    },
    map = {
      if (!.gc_is_map(node)) return(.gc_problem(acc, where, node, paste0("expects a mapping, not ", .gc_shape(node))))
      out <- stats::setNames(vector("list", length(node)), names(node))
      known <- rep(TRUE, length(node))
      for (i in seq_along(node)) {
        k <- names(node)[i]; w <- paste0(where, ".", k)
        fs <- if (!is.null(spec$fields)) spec$fields[[k]] else spec$values
        if (is.null(fs)) {
          .gc_problem(acc, w, node[[i]], paste0("is not a field of `", where, "`", .gc_nearest(k, names(spec$fields))))
          known[i] <- FALSE
          next
        }
        out[i] <- list(.gc_convert(node[[i]], fs, w, acc))
      }
      out <- out[known]
      if (!partial) for (f in setdiff(spec$required, names(node)))
        .gc_problem(acc, where, node, paste0("needs the field `", f, "`"))
      # a map with declared fields comes back in the declared order, whatever order the file
      # used, so that a written file reads back to an identical value
      if (!is.null(spec$fields)) out <- out[intersect(names(spec$fields), names(out))]
      out
    },
    sequence = {
      if (!.gc_is_seq(node)) return(.gc_problem(acc, where, node, paste0("expects a list, not ", .gc_shape(node))))
      lapply(seq_along(node), function(i) .gc_convert(node[[i]], spec$items, sprintf("%s[%d]", where, i), acc))
    },
    one_of = {
      tried <- vector("list", length(spec$one_of))
      for (j in seq_along(spec$one_of)) {
        try_acc <- new.env(); try_acc$problems <- character()
        v <- .gc_convert(node, spec$one_of[[j]], where, try_acc)
        if (!length(try_acc$problems)) return(v)
        tried[[j]] <- try_acc$problems
      }
      # When the value has the shape of exactly one form, that form's own problems say what
      # is wrong (a bad field of a grid, say). Otherwise the forms are listed.
      shaped <- which(vapply(spec$one_of, function(a) .gc_shape_fits(node, a), logical(1)))
      if (length(shaped) == 1L) {
        acc$problems <- c(acc$problems, tried[[shaped]])
        return(invisible(NULL))
      }
      forms <- vapply(spec$one_of, .gc_spec_words, character(1))
      .gc_problem(acc, where, node, paste0("fits none of the forms this setting takes: ", paste(forms, collapse = "; or ")))
    },
    r_function = {
      if (!.gc_is_scalar(node)) return(.gc_problem(acc, where, node, "expects the text of one R function"))
      exprs <- tryCatch(parse(text = node$text, keep.source = FALSE), error = function(e) NULL)
      if (is.null(exprs) || length(exprs) != 1L || !is.call(exprs[[1L]]) || !identical(exprs[[1L]][[1L]], as.name("function")))
        return(.gc_problem(acc, where, node, "expects the text of exactly one R function definition, such as function(maf) sqrt(maf * (1 - maf))"))
      # Evaluating a `function(...) body` expression only builds the closure; neither
      # the body nor a default argument runs. The function gets the global environment,
      # so it sees base R and the attached packages (stats among them) when it is called
      # later, like a function typed at the prompt. The source text is kept for the writer.
      fn <- eval(exprs[[1L]], envir = globalenv())
      attr(fn, "glow_source") <- trimws(node$text)
      fn
    })
}

#' Whether a node has the shape a type takes
#' @param node A node of a parsed tree.
#' @param spec A type spec.
#' @return `TRUE` when the node is a scalar, a list or a map as the type expects.
#' @keywords internal
#' @noRd
.gc_shape_fits <- function(node, spec) {
  switch(spec$type,
    logical = , integer = , number = , string = , path = , file_name = , r_function = .gc_is_scalar(node),
    integer_list = , number_list = , string_list = .gc_is_scalar(node) || .gc_is_seq(node),
    sequence = .gc_is_seq(node),
    named_vector = , map = .gc_is_map(node),
    one_of = TRUE)
}

#' A type spec in a few words, for the message of a `one_of`
#' @param spec A type spec.
#' @return A phrase such as "a list of numbers".
#' @keywords internal
#' @noRd
.gc_spec_words <- function(spec) {
  switch(spec$type,
    logical = "true or false", integer = "a whole number", number = "a number",
    string = if (is.null(spec$choices)) "text" else paste0("one of ", paste(spec$choices, collapse = ", ")),
    path = "a path", file_name = "the name of one file or folder", integer_list = "a list of whole numbers", number_list = "a list of numbers",
    string_list = "a list of text items", named_vector = "a mapping from names to values",
    map = if (!is.null(spec$fields))
            paste0("a mapping with ", paste(vapply(names(spec$fields), function(k)
              paste0(k, " (", .gc_spec_words(spec$fields[[k]]), ")"), character(1)), collapse = ", "))
          else "a mapping",
    sequence = "a list", r_function = "the text of an R function", one_of = "one of several forms")
}

#' The nearest registered name, as a question
#'
#' Near means an edit distance of at most half the length of the name (and at least
#' 2), or one name being the beginning of the other, as `min_vars` is of
#' `min_variants` up to its last letter.
#' @param name The name the user wrote.
#' @param candidates The registered names.
#' @return "; did you mean `x`?", or "" when no candidate is near.
#' @keywords internal
#' @noRd
.gc_nearest <- function(name, candidates) {
  if (!length(candidates)) return("")
  d <- utils::adist(name, candidates, ignore.case = TRUE)[1L, ]
  best <- which.min(d)
  stem <- substr(name, 1L, max(3L, nchar(name) - 1L))
  near <- d[best] <= max(2L, nchar(name) %/% 2L) ||
    startsWith(tolower(candidates[best]), tolower(stem)) || startsWith(tolower(name), tolower(candidates[best]))
  if (near) paste0("; did you mean `", candidates[best], "`?") else ""
}

# ---- Checking a typed value against its type ------------------------------------

#' What is wrong with a typed R value under a type spec
#'
#' The counterpart of `.gc_convert()` for a value that is already in R: the writer
#' uses it before it writes, and to find the form of a `one_of` that a value has.
#' It checks everything the reader checks: the R type, whole numbers and the range
#' of an integer, infinity against `inf_ok`, words against `choices`, the keys of a
#' map against its fields, and every part of a value with parts.
#' @param value The R value.
#' @param spec Its type spec.
#' @param exact Whether an empty list has to say what it is: an empty map then needs
#'   names (`setNames(list(), character(0))`) and an empty sequence none. Used to
#'   tell the forms of a `one_of` apart.
#' @return "" when the value is valid, else a few words on what is wrong.
#' @keywords internal
#' @noRd
.gc_valid <- function(value, spec, exact = FALSE) {
  if (is.null(value)) return(if (isTRUE(spec$null_ok)) "" else "is NULL, and its type does not allow a null")
  num_ok <- function(v, whole, inf_ok) {
    if (!is.numeric(v) || anyNA(v)) return("is not a number")
    if (any(is.infinite(v)) && !isTRUE(inf_ok)) return("is infinite, and its type does not allow infinity")
    fin <- v[is.finite(v)]
    if (whole && any(fin != round(fin))) return("is not a whole number")
    if (whole && any(abs(fin) > .Machine$integer.max)) return("is outside the range of an integer")
    ""
  }
  words_ok <- function(v, choices) {
    if (!is.character(v) || anyNA(v)) return("is not text")
    if (!is.null(choices) && !all(v %in% choices)) return(paste0("is not one of: ", paste(choices, collapse = ", ")))
    ""
  }
  one <- function(v) length(v) == 1L && is.null(names(v))
  switch(spec$type,
    logical = if (is.logical(value) && one(value) && !is.na(value)) "" else "is not TRUE or FALSE",
    integer = if (!one(value)) "is not a single number" else num_ok(value, TRUE, spec$inf_ok),
    number  = if (!one(value)) "is not a single number" else num_ok(value, FALSE, spec$inf_ok),
    string  = , path = if (!is.character(value) || !one(value)) "is not a single text value" else words_ok(value, spec$choices),
    file_name = if (!is.character(value) || !one(value) || is.na(value)) "is not a single text value"
                else sub("^expects", "is not", .gc_file_name_problem(value)),
    integer_list = if (!is.null(names(value))) "is a named vector" else num_ok(value, TRUE, spec$inf_ok),
    number_list  = if (!is.null(names(value))) "is a named vector" else num_ok(value, FALSE, spec$inf_ok),
    string_list  = if (!is.null(names(value))) "is a named vector" else words_ok(value, spec$choices),
    named_vector = {
      if (!is.atomic(value)) return("is not a named vector")
      nm <- names(value)
      if (length(value) && (is.null(nm) || any(!nzchar(nm)) || anyDuplicated(nm))) return("needs a distinct name for every value")
      if (exact && is.null(nm)) return("has no names")
      vs <- spec$values
      switch(vs$type,
        logical = if (is.logical(value) && !anyNA(value)) "" else "holds a value that is not TRUE or FALSE",
        integer = num_ok(unname(value), TRUE, vs$inf_ok),
        number  = num_ok(unname(value), FALSE, vs$inf_ok),
        string  = words_ok(unname(value), vs$choices))
    },
    map = {
      if (!is.list(value) || is.function(value)) return("is not a mapping")
      nm <- names(value)
      if (!length(value)) return(if (exact && is.null(nm)) "is an empty list without names"
                                 else if (length(spec$required)) paste0("lacks the field `", spec$required[1L], "`") else "")
      if (is.null(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) return("needs a distinct name for every entry")
      lacking <- setdiff(spec$required, nm)
      if (length(lacking)) return(paste0("lacks the field `", lacking[1L], "`"))
      for (k in nm) {
        fs <- if (!is.null(spec$fields)) spec$fields[[k]] else spec$values
        if (is.null(fs)) return(paste0("has the entry `", k, "`, which is not one of its fields"))
        msg <- .gc_valid(value[[k]], fs)
        if (nzchar(msg)) return(paste0("has the entry `", k, "`, which ", msg))
      }
      ""
    },
    sequence = {
      if (!is.list(value) || is.function(value) || !is.null(names(value))) return("is not a list without names")
      for (i in seq_along(value)) {
        msg <- .gc_valid(value[[i]], spec$items)
        if (nzchar(msg)) return(sprintf("has the item %d, which %s", i, msg))
      }
      ""
    },
    one_of = if (any(vapply(spec$one_of, function(a) !nzchar(.gc_valid(value, a, exact = TRUE)), logical(1)))) ""
             else "fits none of the forms of its type",
    r_function = if (is.function(value)) "" else "is not a function")
}
