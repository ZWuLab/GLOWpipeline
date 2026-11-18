# ==============================================================================
# File: test-glow_config_write.R
# Purpose: Tests of write_resolved_config() and write_base_config()
#   (R/glow_config_write.R): a written file reads back to identical settings, and
#   the generated base is a function of the registry alone

# The settings of a configuration without its attributes, functions as their source.
gc_plain_settings <- function(cfg) {
  strip <- function(x) {
    if (is.function(x)) return(attr(x, "glow_source"))
    if (is.list(x)) return(lapply(x, strip))
    x
  }
  strip(unclass(cfg)[names(cfg)])
}

# Write a configuration and read the written file back.
gc_round_trip <- function(s, cfg) {
  f <- write_resolved_config(cfg, withr::local_tempdir(.local_envir = parent.frame()))
  read_glow_config(f, base = s$base, registry = s$registry, root = s$dir)
}

test_that("a resolved configuration reads back to identical settings, for every type", {
  s <- gc_setup(withr::local_tempdir())
  cfg <- gc_read(s, c(
    "dataset: cohort0203",
    "analysis: \"007\"",
    "limits: {b: 7}",
    "shared:",
    "  chroms: [1, 2, 22]",
    "  pi_features: [cadd, \"no\", n, \"1e6\", linsight, .na.character, \".inf\"]",
    "  label: \"a: b # not a comment\"",
    "  gds_dir: /data/with space/gds",
    "data_prep:",
    "  outcome: {col: status, map: {\"1\": 0, \"2\": 1}}",
    "  covariates:",
    "    sex: {col: SEX, map: {M: 0, F: 1.5}}",
    "    age: AGE",
    "snv_set:",
    "  ld_threshold: 0.1",
    "  gene_max_records: .inf",
    "  staar_enabled: true",
    "  output_format: fst",
    "  use_spa: false",
    "  alpha: 2.5e-8",
    "  weights: [0.1, 3.141592653589793, 1e300, 5e-324]",
    "  categories: {lof: plof, mis: missense}",
    "  b_func: |",
    "    function(maf) {",
    "      sqrt(maf * (1 - maf))   # the closed form, with a \"quote\"",
    "    }",
    "  filter_spec:",
    "    variant_type: Indel",
    "    annotation_clauses:",
    "      - {cadd_phred: {op: gt, value: 20}, category: [exonic, splicing]}",
    "      - {genehancer: {op: nonempty}}",
    "training:",
    "  b_target_maf: {n: 50, to: 0.5, from: 0.001}"))
  out <- withr::local_tempdir()
  f <- write_resolved_config(cfg, out, provenance = list(glowanalyses = "0.4.0", command = "glow run snv_set"))
  expect_identical(basename(f), "config-resolved.yml")
  back <- read_glow_config(f, base = s$base, registry = s$registry, root = s$dir)
  expect_identical(gc_plain_settings(back), gc_plain_settings(cfg))
  expect_equal(back$snv_set$b_func(0.5), 0.5)
  # the record shows the spelling the reader takes
  txt <- readLines(f)
  expect_true(any(grepl("^  staar_enabled: true$", txt)))
  expect_true(any(grepl("^  use_spa: false$", txt)))
  expect_true(any(grepl("^analysis: \"007\"$", txt)))
  expect_true(any(grepl("^  gene_max_records: \\.inf$", txt)))
  expect_false(any(grepl(": (yes|no)$", txt)))
  # a second write of what was read keeps the settings (the provenance describes the new chain)
  again <- gc_round_trip(s, back)
  expect_identical(gc_plain_settings(again), gc_plain_settings(cfg))
})

test_that("the base alone round-trips, nulls included", {
  s <- gc_setup(withr::local_tempdir())
  cfg <- read_glow_config(NULL, base = s$base, registry = s$registry, root = s$dir)
  f <- write_resolved_config(cfg, withr::local_tempdir())
  back <- read_glow_config(f, base = s$base, registry = s$registry, root = s$dir)
  expect_identical(gc_plain_settings(back), gc_plain_settings(cfg))
  expect_true(any(grepl("^  chroms: 1-22$", readLines(f))))            # the range form
})

test_that("every tested double reads back unchanged through a written file", {
  set.seed(20261116)
  x <- c(.Machine$double.xmin, .Machine$double.xmax, 5e-324, 1 + .Machine$double.eps, 1 - .Machine$double.eps / 2,
         0.1, 1 / 3, pi, 2.5e-6, 1e6, 3e9, 1e15, 123456789012345678, -7.25e-12, 0, -0,
         2^stats::runif(3000, -1070, 1020) * sample(c(-1, 1), 3000, replace = TRUE), stats::runif(1000))
  txt <- vapply(x, .gc_num_text, character(1))
  expect_true(all(grepl(.gc_num_pattern, txt)))             # the writer stays inside the reader's grammar
  expect_identical(as.numeric(txt), x)
  expect_identical(.gc_num_text(0.1), "0.1")                # and uses no more digits than needed
  expect_identical(.gc_num_text(1e6), "1000000")
  expect_identical(.gc_num_text(Inf), ".inf")
  # the complete round trip: all of them written to a file as a setting and read back
  s <- gc_setup(withr::local_tempdir())
  cfg <- read_glow_config(NULL, base = s$base, registry = s$registry, root = s$dir)
  cfg$snv_set$weights <- x
  expect_identical(gc_round_trip(s, cfg)$snv_set$weights, x)
})

test_that("text that needs quotes gets them, and plain text stays plain", {
  expect_identical(.gc_scalar_text("cadd_phred"), "cadd_phred")
  expect_identical(.gc_scalar_text("data-example/cohort/pheno.csv"), "data-example/cohort/pheno.csv")
  for (w in c("true", "No", "y", "null", "007", "1e6", ".5", ".inf", ".na", ".na.character", ".na.real", "2026-11-16",
              "chr{chr}.gds", "a b", "a: b", "#x", "-x", "", "~", "[a]", "a,b", "+1", "1-22", "0x10"))
    expect_match(.gc_scalar_text(w), "^\".*\"$", info = w)
  expect_identical(.gc_scalar_text("say \"hi\"\\now"), "\"say \\\"hi\\\"\\\\now\"")
})

test_that("every string a setting can hold reads back unchanged", {
  s <- gc_setup(withr::local_tempdir())
  cfg <- read_glow_config(NULL, base = s$base, registry = s$registry, root = s$dir)
  hard <- c(" leading", "trailing ", "tab\there", "cr\rhere", "two\nlines", "\nstarts with a line break", "ends with two\n\n",
            "trailing blank \nin a line", "say \"hi\"", "back\\slash", "caf\u00e9 \u4e2d\u6587", "# not a comment", "a: b",
            "-dash", "?q", "@at", "`tick`", "%pct", "&amp", "*star", "!bang", "|pipe", ">gt", "a, b", "[x]", "{y}",
            "true", "no", "null", "~", "Inf", ".inf", "NaN", "1-22", "0x10", "1e6", "007", ".na.character",
            # control characters and the line breaks only YAML knows
            "bell\u0007", "back\u0008space", "esc\u001b", "del\u007f", "nel\u0085x", "c1\u009f", "ls\u2028x", "ps\u2029x",
            "bom\ufeffx")
  for (v in hard) {
    cfg$shared$label <- v
    expect_identical(gc_round_trip(s, cfg)$shared$label, v, info = encodeString(v))
  }
  cfg$shared$label <- "x"
  cfg$shared$pi_features <- hard[!grepl("\n", hard)]            # the same texts as the items of a list
  expect_identical(gc_round_trip(s, cfg)$shared$pi_features, cfg$shared$pi_features)
  cfg$shared["pi_features"] <- list(NULL)
  cfg$data_prep$covariates <- stats::setNames(as.list(rep("col", 6)), c("a b", "007", "true", "x: y", "#k", "nel\u0085k"))
  expect_identical(gc_round_trip(s, cfg)$data_prep$covariates, cfg$data_prep$covariates)     # keys too
})

test_that("lists of integers are written in a form the reader takes", {
  s <- gc_setup(withr::local_tempdir())
  cfg <- read_glow_config(NULL, base = s$base, registry = s$registry, root = s$dir)
  for (v in list(c(-2L, -1L, 0L), c(3L, 2L, 1L), 22L, c(1L, 2L), c(5L, 7L, 9L), integer(0), 0:3, c(21L, 22L, 23L))) {
    cfg$shared$chroms <- v
    expect_identical(gc_round_trip(s, cfg)$shared$chroms, v, info = paste(v, collapse = ","))
  }
  cfg$shared$chroms <- 0:3
  expect_true(any(grepl("^  chroms: 0-3$", readLines(write_resolved_config(cfg, withr::local_tempdir())))))
  cfg$shared$chroms <- c(-2L, -1L, 0L)                          # a negative run is not written as a range
  expect_true(any(grepl("^  chroms: \\[-2, -1, 0\\]$", readLines(write_resolved_config(cfg, withr::local_tempdir())))))
})

test_that("the form of a one_of is found by its content, not by its shape", {
  d <- withr::local_tempdir()
  s <- gc_setup(d, c(gc_test_registry_text(),
    "- name: two_maps",
    "  section: shared",
    "  type: one_of",
    "  one_of:",
    "    - {type: map, fields: {a: {type: number}}}",
    "    - {type: map, fields: {b: {type: number}}}",
    "  default: {a: 1}",
    "  description: Two forms that are both maps.",
    "- name: map_or_list",
    "  section: shared",
    "  type: one_of",
    "  one_of:",
    "    - {type: map, values: {type: string}}",
    "    - {type: sequence, items: {type: string}}",
    "  default: []",
    "  description: A map with free keys, or a list."))
  cfg <- gc_read(s, c("shared:", "  two_maps: {b: 2}"))
  expect_identical(gc_round_trip(s, cfg)$shared$two_maps, list(b = 2))
  expect_identical(cfg$shared$map_or_list, list())               # an empty list stays an empty list
  expect_identical(gc_round_trip(s, cfg)$shared$map_or_list, list())
  cfg2 <- gc_read(s, c("shared:", "  map_or_list: {}"), name = "u2.yml")
  expect_identical(gc_round_trip(s, cfg2)$shared$map_or_list, stats::setNames(list(), character(0)))
  # every form of the settings of the test registry
  for (txt in list(c("snv_set:", "  categories: [plof]"), c("snv_set:", "  categories: {lof: plof}"),
                   c("training:", "  b_target_maf: [0.1]"), c("training:", "  b_target_maf: {from: 0, to: 1}"),
                   c("data_prep:", "  outcome: {col: y, map: {a: 0.5}}"),
                   c("data_prep:", "  outcome: {col: y, map: \"function(v) as.numeric(v == 2)\"}"),
                   c("data_prep:", "  outcome: {col: y, map: ~}"))) {
    cf <- gc_read(s, txt, name = "form.yml")
    expect_identical(gc_plain_settings(gc_round_trip(s, cf)), gc_plain_settings(cf), info = paste(txt, collapse = " "))
  }
})

test_that("the provenance block is written, and skipped when the file is read", {
  s <- gc_setup(withr::local_tempdir())
  cfg <- gc_read(s, "dataset: d1")
  f <- write_resolved_config(cfg, withr::local_tempdir(),
                             provenance = list(command = "glow run \"x\"\nsecond line", parts = c("a", "b"), none = character(0),
                                               "odd key: #" = 1L))
  txt <- readLines(f)
  i <- which(txt == "provenance:")
  expect_length(i, 1L)
  prov <- txt[i:length(txt)]
  expect_true(any(grepl("base_config.yml", prov, fixed = TRUE)))
  expect_true(any(grepl("user.yml", prov, fixed = TRUE)))
  expect_true(any(grepl(unname(tools::md5sum(s$registry)), prov, fixed = TRUE)))
  expect_true(any(grepl("GLOWpipeline:", prov, fixed = TRUE)))
  expect_true(any(grepl("^  none: \\[\\]$", prov)))               # an empty vector is an empty list
  # the block parses, with the caller's entries as given, and the reader skips it
  tree <- .gc_parse_file(f)
  expect_identical(tree$provenance$command$text, "glow run \"x\"\nsecond line")
  expect_identical(vapply(tree$provenance$parts, function(z) z$text, ""), c("a", "b"))
  expect_length(tree$provenance$none, 0L)
  expect_identical(tree$provenance[["odd key: #"]]$text, "1")
  expect_identical(read_glow_config(f, base = s$base, registry = s$registry, root = s$dir)$dataset, "d1")
  # what the caller may not pass
  expect_error(write_resolved_config(cfg, tempdir(), provenance = list(files = "x")), "cannot use the names")
  expect_error(write_resolved_config(cfg, tempdir(), provenance = list("x")), "named list")
  expect_error(write_resolved_config(cfg, tempdir(), provenance = list(x = "a", x = "b")), "uses the name `x` twice")
  expect_error(write_resolved_config(cfg, tempdir(), provenance = list(x = NA)), "plain vector")
  expect_error(write_resolved_config(list(), tempdir()), "read_glow_config")
})

test_that("a value outside its registered type stops the writer, with the setting named", {
  s <- gc_setup(withr::local_tempdir())
  base <- gc_read(s, "dataset: d1")
  stops <- function(change, pattern) {
    cfg <- base; cfg <- change(cfg)
    expect_error(write_resolved_config(cfg, withr::local_tempdir()), pattern)
  }
  stops(function(c) { c$snv_set$ld_threshold <- "high"; c }, "`snv_set.ld_threshold` is not a number")
  stops(function(c) { c$snv_set["ld_threshold"] <- list(NULL); c }, "`snv_set.ld_threshold` is NULL")
  stops(function(c) { c$snv_set$mac_threshold <- 1.5; c }, "`snv_set.mac_threshold` is not a whole number")
  stops(function(c) { c$snv_set$mac_threshold <- 3e9; c }, "`snv_set.mac_threshold` is outside the range of an integer")
  stops(function(c) { c$snv_set$alpha <- Inf; c }, "`snv_set.alpha` is infinite")
  stops(function(c) { c$snv_set$output_format <- "xml"; c }, "`snv_set.output_format` is not one of: csv, fst")
  stops(function(c) { c$snv_set$staar_enabled <- NA; c }, "`snv_set.staar_enabled` is not TRUE or FALSE")
  stops(function(c) { c$snv_set$filter_spec$min_vars <- 2L; c }, "`snv_set.filter_spec` has the entry `min_vars`, which is not one of its fields")
  stops(function(c) { c$snv_set$filter_spec$min_mac <- 0.5; c }, "has the entry `min_mac`, which is not a whole number")
  stops(function(c) { c$training$b_target_maf <- list(n = 5L); c }, "`training.b_target_maf` fits none of the forms")
  stops(function(c) { c$shared$chroms <- c(1.5, 2); c }, "`shared.chroms` is not a whole number")
  stops(function(c) { c$snv_set$b_func <- "function(x) x"; c }, "`snv_set.b_func` is not a function")
  # a value changed to another valid one is written
  ok <- base; ok$snv_set$gene_max_records <- Inf; ok$snv_set$mac_threshold <- 12
  expect_identical(gc_round_trip(s, ok)$snv_set$mac_threshold, 12L)
})

test_that("the generated base is a function of the registry alone", {
  s <- gc_setup(withr::local_tempdir())
  first <- readLines(s$base)
  expect_identical(write_base_config(s$registry, NULL), first)        # a second run gives the same lines
  expect_identical(write_base_config(s$registry, s$base), first)
  expect_identical(readLines(s$base), first)
  # every setting is there with its description, and the file says it is generated
  expect_true(any(grepl("GENERATED", first)))
  expect_true(any(grepl("^  # The LD threshold\\.$", first)))
  expect_true(any(grepl("^  ld_threshold: 0.9$", first)))
  expect_true(any(grepl("^  # Allowed: csv, fst\\.$", first)))
  expect_true(any(grepl("^  # Null: found through the layout", first)))
  # a changed registry changes the lines, which is what a staleness check compares
  writeLines(sub("default: 0.9", "default: 0.8", gc_test_registry_text(), fixed = TRUE), s$registry)
  expect_false(identical(write_base_config(s$registry, NULL), first))
})
