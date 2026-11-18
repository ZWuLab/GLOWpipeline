# ==============================================================================
# File: helper-glow-config.R
# Purpose: A small settings registry and helpers shared by the tests of the YAML
#   configuration loader (test-glow_config_registry.R, test-glow_config_read.R,
#   test-glow_config_write.R). The registry holds one setting of every registered
#   type, so the tests need no GLOWanalyses tree.

gc_test_registry_text <- function() c(
  "- {name: dataset, section: top, type: string, default: example, description: The dataset.}",
  "- {name: analysis, section: top, type: string, default: baseline, description: The analysis.}",
  "- name: limits",
  "  section: top",
  "  type: map",
  "  fields: {a: {type: integer}, b: {type: integer}}",
  "  default: {a: 1, b: 2}",
  "  description: A top-level map with declared fields.",
  "- {name: chroms, section: shared, type: integer_list, default: 1-22, description: Chromosomes.}",
  "- {name: pi_features, section: shared, type: string_list, null_ok: true, default: ~, description: Features.}",
  "- {name: label, section: shared, type: string, default: x, description: A label.}",
  "- {name: gds_dir, section: shared, type: path, null_ok: true, derived: true, default: ~, description: The GDS folder.}",
  "- {name: ld_threshold, section: snv_set, type: number, default: 0.9, description: The LD threshold.}",
  "- {name: mac_threshold, section: snv_set, type: integer, default: 10, description: The MAC threshold.}",
  "- {name: gene_max_records, section: snv_set, type: integer, inf_ok: true, default: 10000, description: Records.}",
  "- {name: staar_enabled, section: snv_set, type: logical, default: false, description: The comparator.}",
  "- {name: output_format, section: snv_set, type: string, choices: [csv, fst], default: csv, description: The format.}",
  "- {name: use_spa, section: snv_set, type: logical, null_ok: true, default: ~, description: SPA.}",
  "- {name: b_func, section: snv_set, type: r_function, null_ok: true, default: ~, description: The B function.}",
  "- {name: alpha, section: snv_set, type: number, default: 2.5e-6, description: The level.}",
  "- {name: weights, section: snv_set, type: number_list, default: [1, 25], description: Weights.}",
  "- name: filter_spec",
  "  section: snv_set",
  "  type: map",
  "  fields:",
  "    rare_maf_cutoff: {type: number}",
  "    variant_type: {type: string, choices: [SNV, Indel, variant]}",
  "    min_mac: {type: integer}",
  "    annotation_clauses:",
  "      type: sequence",
  "      null_ok: true",
  "      items:",
  "        type: map",
  "        values:",
  "          type: one_of",
  "          one_of:",
  "            - {type: string_list}",
  "            - type: map",
  "              fields:",
  "                op: {type: string}",
  "                value: {type: number, null_ok: true}",
  "  default: {rare_maf_cutoff: 0.5, variant_type: SNV, min_mac: 1, annotation_clauses: ~}",
  "  description: The variant filter.",
  "- name: categories",
  "  section: snv_set",
  "  type: one_of",
  "  one_of:",
  "    - {type: string_list}",
  "    - {type: map, values: {type: string}}",
  "  default: [plof, missense]",
  "  description: The coding categories.",
  "- name: outcome",
  "  section: data_prep",
  "  type: map",
  "  fields:",
  "    col: {type: string}",
  "    map:",
  "      type: one_of",
  "      null_ok: true",
  "      one_of:",
  "        - {type: named_vector, values: {type: integer}}",
  "        - {type: named_vector, values: {type: number}}",
  "        - {type: r_function}",
  "  default: {col: outcome, map: {\"0\": 0, \"1\": 1}}",
  "  description: The outcome column and its recoding.",
  "- name: source",
  "  section: data_prep",
  "  type: map",
  "  replace_whole: true",
  "  fields: {node: {type: string, null_ok: true}, col: {type: string, null_ok: true}}",
  "  default: {col: outcome}",
  "  description: A map whose fields are alternatives.",
  "- name: covariates",
  "  section: data_prep",
  "  type: map",
  "  values:",
  "    type: one_of",
  "    one_of:",
  "      - {type: string}",
  "      - type: map",
  "        fields:",
  "          col: {type: string}",
  "          map: {type: named_vector, values: {type: number}}",
  "  default: {sex: sex, age: age}",
  "  description: The covariates, in the order of the design matrix.",
  "- name: b_target_maf",
  "  section: training",
  "  type: one_of",
  "  null_ok: true",
  "  one_of:",
  "    - {type: number_list}",
  "    - type: map",
  "      required: [from, to]",
  "      fields: {from: {type: number}, to: {type: number}, n: {type: integer, null_ok: true}}",
  "  default: ~",
  "  description: Target MAFs for the B table.")

# Write the test registry and its generated base into `dir`; return their paths.
gc_setup <- function(dir, registry_text = gc_test_registry_text()) {
  reg <- file.path(dir, "settings-registry.yml")
  base <- file.path(dir, "base_config.yml")
  writeLines(registry_text, reg)
  write_base_config(reg, base)
  list(dir = dir, registry = reg, base = base)
}

# Write a user configuration into the test folder and return its path.
gc_file <- function(s, name, text) {
  f <- file.path(s$dir, name)
  writeLines(text, f)
  f
}

# Read a user configuration given as text against the test registry.
gc_read <- function(s, text, name = "user.yml") {
  read_glow_config(gc_file(s, name, text), base = s$base, registry = s$registry, root = s$dir)
}

# Convert one YAML value by a type spec written as YAML; returns the value and the problems.
gc_conv <- function(value_text, spec_text) {
  d <- withr::local_tempdir()
  vf <- file.path(d, "v.yml"); sf <- file.path(d, "s.yml")
  writeLines(paste0("x: ", value_text), vf)
  writeLines(paste0("x: ", spec_text), sf)
  spec <- .gc_spec(.gc_parse_file(sf)$x, "x", sf)
  acc <- new.env(); acc$problems <- character()
  v <- .gc_convert(.gc_parse_file(vf)$x, spec, "x", acc)
  list(value = v, problems = acc$problems)
}
