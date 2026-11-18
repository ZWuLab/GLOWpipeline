# ==============================================================================
# File: test-glow_config_registry.R
# Purpose: Tests of the internals in R/glow_config_registry.R: reading a YAML file
#   as written, the registry, and the conversion of a value by its registered type

# Write bytes or lines to a file and return the message of the refusal, or "" when it reads.
gc_refusal <- function(f, content) {
  if (is.raw(content)) writeBin(content, f) else writeLines(content, f)
  tryCatch({ .gc_parse_file(f); "" }, error = function(e) conditionMessage(e))
}

# ---- Reading a file as written -------------------------------------------------

test_that("every scalar arrives as the text the user wrote", {
  d <- withr::local_tempdir()
  f <- file.path(d, "a.yml")
  writeLines(c("a: Y", "b: no", "c: 007", "d: 010", "e: 0x10", "f: 3000000000", "g: 1e6",
               "h: 2026-11-16", "i: ~", "j: \"null\"", "k: [1]", "l: 1", "m: []", "n: {}",
               "o: .na", "p: .na.character", "q: .na.integer", "r: .na.real"), f)
  t <- .gc_parse_file(f)
  txt <- function(k) t[[k]]$text
  expect_identical(txt("a"), "Y")
  expect_identical(txt("b"), "no")
  expect_identical(txt("c"), "007")
  expect_identical(txt("d"), "010")
  expect_identical(txt("e"), "0x10")
  expect_identical(txt("f"), "3000000000")
  expect_identical(txt("g"), "1e6")
  expect_identical(txt("h"), "2026-11-16")
  expect_true(.gc_is_null(t$i))
  expect_false(.gc_is_null(t$j))            # a quoted "null" is text
  expect_true(.gc_is_seq(t$k))              # a list of one item is still a list
  expect_true(.gc_is_scalar(t$l))
  expect_true(.gc_is_seq(t$m))              # an empty list
  expect_true(.gc_is_map(t$n))              # differs from an empty map
  # the parser's four words for a missing value are text like any other
  expect_identical(c(txt("o"), txt("p"), txt("q"), txt("r")), c(".na", ".na.character", ".na.integer", ".na.real"))
})

test_that("an empty file, a comment-only file, a byte-order mark and a directive are accepted", {
  d <- withr::local_tempdir()
  f <- file.path(d, "empty.yml"); file.create(f)
  expect_length(.gc_parse_file(f), 0L)
  expect_true(.gc_is_map(.gc_parse_file(f)))
  writeLines("# nothing here", f)
  expect_length(.gc_parse_file(f), 0L)
  writeLines(c("   ", "", "  # only blanks and a comment"), f)
  expect_length(.gc_parse_file(f), 0L)
  writeBin(c(as.raw(c(0xEF, 0xBB, 0xBF)), charToRaw("a: 1\n")), f)
  expect_identical(.gc_parse_file(f)$a$text, "1")
  writeLines(c("%YAML 1.1", "---", "a: only"), f)                 # one document with a directive
  expect_identical(.gc_parse_file(f)$a$text, "only")
  writeLines(c("---", "a: 1"), f)                                 # a leading document marker
  expect_identical(.gc_parse_file(f)$a$text, "1")
  writeLines(c("a: |", "  text", "  --- not a marker", "b: 1 # --- in a comment"), f)
  expect_identical(.gc_parse_file(f)$b$text, "1")
})

test_that("an anchor and its alias are plain YAML", {
  d <- withr::local_tempdir()
  f <- file.path(d, "a.yml")
  writeLines(c("a: &v [x, y]", "b: *v", "c: {<<: {p: 1}, q: 2}"), f)
  t <- .gc_parse_file(f)
  expect_identical(vapply(t$b, function(z) z$text, ""), c("x", "y"))
  expect_identical(names(t$c), c("p", "q"))                       # a merge that replaces nothing
})

test_that("the YAML a configuration may not use is refused, with the file named", {
  d <- withr::local_tempdir()
  f <- file.path(d, "bad.yml")
  bad <- function(content, pattern) {
    msg <- gc_refusal(f, content)
    expect_match(msg, pattern, info = paste(pattern, "|", msg))
    expect_match(msg, "bad.yml", fixed = TRUE, info = msg)        # every refusal names the file
  }
  bad(c("a: 1", "---", "a: 2"), "more than one YAML document")
  bad(c("---", "---", "a: 2"), "more than one YAML document")
  bad(charToRaw("a: first\r---\ra: second\r"), "more than one YAML document")          # lone CR line ends
  bad(c(charToRaw("a: first"), as.raw(c(0xC2, 0x85)), charToRaw("---"), as.raw(c(0xC2, 0x85)), charToRaw("a: second\n")),
      "more than one YAML document")                                                     # NEL line ends
  bad("a: !expr 1 + 1", "R expression")
  bad("a: !custom abc", "tag this loader does not read")
  bad("a: !custom {b: 1}", "tag this loader does not read")       # a tag on a map
  bad("a: !custom [1, 2]", "tag this loader does not read")       # a tag on a list
  bad("a: !!str [x, y]", "tag that does not fit its value")       # a scalar tag on a list
  bad("a: !!str {}", "tag that does not fit its value")
  bad("a: {<<: {n: 1}, n: 8}", "Duplicate map key")               # a merge that replaces a key
  bad(c("a: 1", "a: 2"), "Duplicate map key")
  bad("a: {\"\": 1}", "empty key")
  bad(c("- a", "- b"), "top level has to be a mapping")
  bad("just text", "top level has to be a mapping")
  bad("~", "top level has to be a mapping")                       # a document that is only a null
  bad("null", "top level has to be a mapping")
  bad("a:\n\tb: 1", "Cannot read")                                # a tab as indentation
  bad(as.raw(c(120, 58, 32, 255, 10)), "not valid UTF-8")
  bad(as.raw(c(120, 58, 32, 0, 10)), "NUL byte")
  bad(paste0("a: ", strrep("[", 200), strrep("]", 200)), "nested more than 64 levels")
  expect_error(.gc_parse_file(file.path(d, "absent.yml")), "not found.*absent.yml")
  expect_error(.gc_parse_file(d), "is a folder")
})

# ---- The type table ------------------------------------------------------------

test_that("logical takes true or false only", {
  expect_true(gc_conv("true", "{type: logical}")$value)
  expect_false(gc_conv("FALSE", "{type: logical}")$value)
  expect_length(gc_conv("True", "{type: logical}")$problems, 0L)
  for (w in c("yes", "no", "on", "Y", "n", "1", "!!bool no", "[true]"))
    expect_match(gc_conv(w, "{type: logical}")$problems, "expects true or false", info = w)
})

test_that("integer reads decimal and scientific notation and refuses the rest", {
  int <- function(x) gc_conv(x, "{type: integer}")
  expect_identical(int("22")$value, 22L)
  expect_identical(int("007")$value, 7L)        # decimal, and never octal
  expect_identical(int("010")$value, 10L)
  expect_identical(int("1e6")$value, 1000000L)
  expect_identical(int("1.0")$value, 1L)
  expect_identical(int("5.")$value, 5L)
  expect_identical(int("100e-2")$value, 1L)
  expect_identical(int("1.5e1")$value, 15L)
  expect_identical(int("-3")$value, -3L)
  expect_identical(int("+5")$value, 5L)
  expect_identical(int("2147483647")$value, .Machine$integer.max)
  expect_identical(int("-2147483647")$value, -.Machine$integer.max)
  expect_match(int("1.5")$problems, "whole number")
  # whole is decided from the text: a double cannot tell these from a whole number
  for (w in c("1e-400", "1.0000000000000001", "2147483647.00000001", "10e-2"))
    expect_match(int(w)$problems, "whole number", info = w)
  for (w in c("3000000000", "2147483648", "-2147483648", "1e400"))
    expect_match(int(w)$problems, "outside the range|too large", info = w)
  for (w in c("0x10", "1_000", "1:22", "abc", ".nan"))
    expect_match(int(w)$problems, "decimal or scientific", info = w)
  expect_match(int(".inf")$problems, "finite")
  expect_identical(gc_conv(".inf", "{type: integer, inf_ok: true}")$value, Inf)
  expect_identical(gc_conv("Inf", "{type: integer, inf_ok: true}")$value, Inf)
  expect_identical(gc_conv("-.inf", "{type: integer, inf_ok: true}")$value, -Inf)
})

test_that("number reads what the default parser loses, and refuses what a double cannot hold", {
  num <- function(x) gc_conv(x, "{type: number}")
  expect_identical(num("0.9")$value, 0.9)
  expect_identical(num(".5")$value, 0.5)
  expect_identical(num("2.5e-6")$value, 2.5e-6)
  expect_identical(num("1e6")$value, 1e6)                   # the default parser reads text
  expect_identical(num("3000000000")$value, 3e9)            # the default parser reads NA
  expect_identical(num("0")$value, 0)
  expect_identical(num("0e-5")$value, 0)
  expect_identical(num(" 5 ")$value, 5)
  expect_match(num("abc")$problems, "expects a number")
  expect_match(num("[1, 2]")$problems, "not a list")
  expect_match(num("1e400")$problems, "too large")
  expect_match(num("1e-400")$problems, "too small")         # it would arrive as zero
  expect_identical(gc_conv("-.inf", "{type: number, inf_ok: true}")$value, -Inf)
  expect_match(num(".inf")$problems, "finite")
})

test_that("a string takes any scalar as it was written", {
  str <- function(x) gc_conv(x, "{type: string}")$value
  expect_identical(str("Y"), "Y")
  expect_identical(str("no"), "no")
  expect_identical(str("007"), "007")
  expect_identical(str("2026-11-16"), "2026-11-16")
  expect_identical(str("1e6"), "1e6")
  expect_identical(str("\"null\""), "null")
  expect_identical(str("chr{chr}.gds"), "chr{chr}.gds")
  expect_identical(str(".na.character"), ".na.character")
  for (type in c("string", "path")) {
    expect_match(gc_conv("[a, b]", sprintf("{type: %s}", type))$problems, "single text value, not a list", info = type)
    expect_match(gc_conv("{a: b}", sprintf("{type: %s}", type))$problems, "single text value, not a mapping", info = type)
  }
  expect_match(gc_conv("~", "{type: string}")$problems, "is null")
  expect_identical(gc_conv("fst", "{type: string, choices: [csv, fst]}")$value, "fst")
  expect_match(gc_conv("tsv", "{type: string, choices: [csv, fst]}")$problems, "one of: csv, fst")
})

test_that("the null spellings are nulls where a null is allowed, and a problem where it is not", {
  for (w in c("~", "null", "Null", "NULL", "")) {
    r <- gc_conv(w, "{type: integer, null_ok: true}")
    expect_null(r$value, info = w)
    expect_length(r$problems, 0L)                # a null, not a failed conversion
    expect_match(gc_conv(w, "{type: integer}")$problems, "is null", info = w)
  }
  r <- gc_conv("~", "{type: string, null_ok: true}")
  expect_null(r$value); expect_length(r$problems, 0L)
})

test_that("lists of scalars: a list, one item, and the range form", {
  il <- function(x) gc_conv(x, "{type: integer_list}")
  expect_identical(il("1-22")$value, 1:22)
  expect_identical(il("1 - 3")$value, 1:3)
  expect_identical(il("[21, 22]")$value, c(21L, 22L))
  expect_identical(il("[-2, -1, 0]")$value, c(-2L, -1L, 0L))
  expect_identical(il("22")$value, 22L)
  expect_identical(il("[]")$value, integer(0))
  expect_match(il("1:22")$problems, "write it as a-b")
  expect_match(il("22-1")$problems, "not a range")
  expect_match(il("[1, x]")$problems, "decimal or scientific")
  expect_match(il("[1, ~]")$problems, "single value")
  expect_match(il("{a: 1}")$problems, "not a mapping")
  expect_identical(gc_conv("[0.5, 1e-3]", "{type: number_list}")$value, c(0.5, 0.001))
  expect_identical(gc_conv("0.5", "{type: number_list}")$value, 0.5)
  expect_match(gc_conv("[0.5, x]", "{type: number_list}")$problems, "expects a number")
  # words that the default parser would turn into logicals stay text
  expect_identical(gc_conv("[a, n, y, no]", "{type: string_list}")$value, c("a", "n", "y", "no"))
  expect_identical(gc_conv("cadd", "{type: string_list}")$value, "cadd")
  expect_identical(gc_conv("[]", "{type: string_list}")$value, character(0))
  expect_match(gc_conv("[a, [b]]", "{type: string_list}")$problems, "single value")
  # choices hold for every item
  expect_identical(gc_conv("[spa, nospa]", "{type: string_list, choices: [spa, nospa]}")$value, c("spa", "nospa"))
  expect_match(gc_conv("[spa, both]", "{type: string_list, choices: [spa, nospa]}")$problems, "one of: spa, nospa")
})

test_that("a named vector is read into a named R vector that keeps its order", {
  nv <- gc_conv("{\"1\": 0, \"2\": 1}", "{type: named_vector, values: {type: integer}}")$value
  expect_identical(nv, c("1" = 0L, "2" = 1L))
  expect_true(is.numeric(unname(nv[as.character(c(2, 1, 2))])))   # what the phenotype step does
  expect_identical(gc_conv("{b: 0.5, a: 2}", "{type: named_vector, values: {type: number}}")$value,
                   c(b = 0.5, a = 2))
  expect_match(gc_conv("{a: [1]}", "{type: named_vector, values: {type: integer}}")$problems, "single value")
  expect_match(gc_conv("[1, 2]", "{type: named_vector, values: {type: integer}}")$problems, "mapping")
})

test_that("maps with declared fields, maps with free keys and sequences", {
  fields <- "{type: map, fields: {a: {type: integer}, b: {type: string, null_ok: true}}}"
  expect_identical(gc_conv("{a: 1, b: ~}", fields)$value, list(a = 1L, b = NULL))
  expect_identical(names(gc_conv("{b: x, a: 1}", fields)$value), c("a", "b"))          # the declared order
  expect_match(gc_conv("{a: 1, bb: x}", fields)$problems, "not a field .* did you mean `b`")
  expect_match(gc_conv("[1, 2]", fields)$problems, "expects a mapping, not a list")
  req <- "{type: map, required: [a], fields: {a: {type: integer}, b: {type: string}}}"
  expect_match(gc_conv("{b: x}", req)$problems, "needs the field `a`")
  expect_length(gc_conv("{a: 1}", req)$problems, 0L)
  free <- "{type: map, values: {type: string}}"
  expect_identical(gc_conv("{z: 1, a: no}", free)$value, list(z = "1", a = "no"))     # the order written
  expect_identical(gc_conv("{}", free)$value, stats::setNames(list(), character(0)))
  expect_match(gc_conv("{a: [x]}", free)$problems, "`x.a` expects a single text value")
  seqs <- "{type: sequence, items: {type: integer}}"
  expect_identical(gc_conv("[1, 2]", seqs)$value, list(1L, 2L))
  expect_identical(gc_conv("[]", seqs)$value, list())
  expect_match(gc_conv("1", seqs)$problems, "expects a list")
  expect_match(gc_conv("{a: 1}", seqs)$problems, "expects a list, not a mapping")
  expect_match(gc_conv("[1, x]", seqs)$problems, "`x\\[2\\]` expects a number")
})

test_that("one_of takes the first form that fits, and says what is wrong when none does", {
  either <- paste0("{type: one_of, one_of: [{type: number_list}, {type: map, required: [from, to], ",
                   "fields: {from: {type: number}, to: {type: number}, n: {type: integer, null_ok: true}}}]}")
  expect_identical(gc_conv("[0.1, 0.2]", either)$value, c(0.1, 0.2))
  expect_identical(gc_conv("{from: 0, to: 0.5}", either)$value, list(from = 0, to = 0.5))
  expect_identical(gc_conv("{n: ~, to: 0.5, from: 0}", either)$value, list(from = 0, to = 0.5, n = NULL))
  # a value shaped like exactly one form gets that form's own message
  expect_match(gc_conv("abc", either)$problems, "expects a number")
  expect_match(gc_conv("{from: 0, to: 0.5, n: many}", either)$problems, "`x.n` expects a number")
  expect_match(gc_conv("{to: 0.5}", either)$problems, "needs the field `from`")
  # two forms of the same shape: the forms are listed
  two <- "{type: one_of, one_of: [{type: map, fields: {a: {type: integer}}}, {type: map, fields: {b: {type: integer}}}]}"
  expect_identical(gc_conv("{b: 2}", two)$value, list(b = 2L))
  expect_match(gc_conv("{c: 2}", two)$problems, "fits none of the forms .* a \\(a whole number\\).* b \\(a whole number\\)")
})

test_that("an R function is parsed, never run, and keeps its source text", {
  f <- gc_conv("\"function(maf) sqrt(maf * (1 - maf))\"", "{type: r_function}")$value
  expect_true(is.function(f))
  expect_equal(f(0.5), 0.5)
  expect_identical(attr(f, "glow_source"), "function(maf) sqrt(maf * (1 - maf))")
  expect_equal(gc_conv("\"function(p) qnorm(p)\"", "{type: r_function}")$value(0.5), 0)   # stats is seen
  # text that is not exactly one function definition is refused without being evaluated
  # forward slashes only: the path is written inside a double-quoted YAML text, where a
  # backslash starts an escape
  marker <- gsub("\\", "/", file.path(withr::local_tempdir(), "ran"), fixed = TRUE)
  for (code in c(sprintf("\"file.create('%s')\"", marker),
                 sprintf("\"local({ file.create('%s'); function(x) x })\"", marker),
                 sprintf("\"(function() { file.create('%s'); function(x) x })()\"", marker),
                 "\"function(x) x; 1\"", "\"function(x\""))
    expect_match(gc_conv(code, "{type: r_function}")$problems, "exactly one R function", info = code)
  # a default argument and a body with a side effect do not run when the function is built
  g <- gc_conv(sprintf("\"function(x = file.create('%s')) file.create('%s')\"", marker, marker), "{type: r_function}")$value
  expect_true(is.function(g))
  expect_false(file.exists(marker))
})

# ---- The registry --------------------------------------------------------------

test_that("the test registry reads, with typed defaults", {
  s <- gc_setup(withr::local_tempdir())
  reg <- .gc_read_registry(s$registry)
  expect_identical(reg$sections[1L], "top")
  expect_identical(reg$entries[["shared.chroms"]]$default, 1:22)
  expect_null(reg$entries[["shared.pi_features"]]$default)
  expect_true(reg$entries[["shared.gds_dir"]]$derived)
  expect_true(reg$entries[["data_prep.source"]]$replace_whole)
  expect_identical(reg$entries[["data_prep.outcome"]]$default$map, c("0" = 0L, "1" = 1L))
  expect_identical(reg$entries[["snv_set.alpha"]]$default, 2.5e-6)
})

test_that("a registry with a mistake is refused", {
  d <- withr::local_tempdir()
  f <- file.path(d, "r.yml")
  bad <- function(lines, pattern) { writeLines(lines, f); expect_error(.gc_read_registry(f), pattern) }
  bad("- {name: a, section: top, type: whatever, default: 1, description: d}", "unknown type")
  bad("- {name: a, section: top, type: integer, description: d}", "has no `default`")
  bad("- {name: a, section: top, type: integer, default: x, description: d}", "invalid default")
  bad(c("- {name: a, section: top, type: integer, default: 1, description: d}",
        "- {name: a, section: top, type: integer, default: 2, description: d}"), "defined twice")
  bad("- {name: a, section: top, type: integer, default: 1, description: d, typo: 1}", "unknown key")
  bad("- {name: extends, section: top, type: string, default: x, description: d}", "reserved")
  bad("- {name: a, section: top, type: map, default: {}, description: d}", "exactly one of `fields` and `values`")
  bad(c("- {name: s, section: top, type: string, default: x, description: d}",
        "- {name: a, section: s, type: integer, default: 1, description: d}"), "both a top-level setting and a section")
  # a null_ok that the loader would never honor
  bad("- {name: a, section: top, type: one_of, one_of: [{type: integer, null_ok: true}], default: 1, description: d}",
      "`null_ok`, which has no effect on a form of a `one_of`")
  bad("- {name: a, section: top, type: named_vector, values: {type: integer, null_ok: true}, default: {}, description: d}",
      "`null_ok`, which has no effect on the values of a `named_vector`")
  bad("- {name: a, section: top, type: map, required: [z], fields: {b: {type: integer}}, default: {}, description: d}",
      "every name in `required` has to be one of its `fields`")
  bad("a: 1", "list of entries")
})
