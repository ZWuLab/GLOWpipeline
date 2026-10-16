# ==============================================================================
# Tests for aggregate_coding_results() — 3-source assembly + merge + significance
# ==============================================================================
#
# Synthetic per-chr GLOW + per-task native CSVs written to a tempdir laid out
# like a coding run's outputs/ (results-glow/, results-staar-nospa/,
# results-staar-spa/); no GDS needed. The load-bearing cases: the STAAR_O ->
# STAAR_O_glowG rename, the native-subset-of-GLOW (gene, category) left join,
# the per-category Bonferroni + genome-wide significance computed from the data,
# the gene-first OUTPUT sort (chr, gene_start, category), the completeness gate,
# and the category-aware native-only reporting.
#
# Native-completeness tests pass a minimal in-test genes_info (.mk_genes_info)
# to aggregate_coding_results(), so they need neither STAARpipeline nor a GDS.


# A synthetic per-chr GLOW table: gene-level identity + counts + GLOW + the
# baked-in STAAR omnibi (STAAR_O / ACAT_O -> renamed *_glowG by the function).
# `cats` is a per-row category vector. Built category-grouped on purpose (the
# input is NOT pre-sorted gene-first), so the gene-first OUTPUT sort is exercised.
.mk_glow_chr <- function(chr, gene, start, category, glow_omni,
                         staar_o = 0.5, acat_o = 0.5) {
  n <- length(gene)
  data.frame(
    label            = gene,
    chr              = as.character(chr),
    category         = category,
    start            = as.integer(start),
    end              = as.integer(start + 1000L),
    midpoint         = as.integer(start + 500L),
    n_variants       = rep(10L, n),
    n_after_annotation = rep(8L, n),
    n_after_filter   = rep(6L, n),
    n_after_ld       = rep(6L, n),
    n_after_collapse = rep(5L, n),
    cMAC             = rep(50, n),
    is_single        = rep(FALSE, n),
    GLOW_Omni        = as.numeric(glow_omni),
    GLOW_SKAT_equ    = as.numeric(glow_omni),
    STAAR_O          = as.numeric(staar_o),
    ACAT_O           = as.numeric(acat_o),
    stringsAsFactors = FALSE, check.names = FALSE)
}

# A synthetic native per-task table (the curated mode-keyed schema; the function
# drops `chr` in the merge and keeps GLOW's). `suffix` = "_native" or "_spa".
.mk_native_task <- function(chr, gene, category, pcol, p) {
  n <- length(gene)
  df <- data.frame(gene = gene, chr = as.character(chr), category = category,
                   stringsAsFactors = FALSE, check.names = FALSE)
  df[[pcol]] <- as.numeric(p)
  df
}

# Lay out a minimal coding-run outputs/ tree in a tempdir and return its path.
.mk_run_dir <- function(glow_by_chr, nospa_by_task = list(), spa_by_task = list()) {
  od <- tempfile("acr_run_"); dir.create(od)
  gdir <- file.path(od, "results-glow"); dir.create(gdir)
  for (chr in names(glow_by_chr)) {
    utils::write.csv(glow_by_chr[[chr]],
                     file.path(gdir, sprintf("glow_coding_chr%s.csv", chr)),
                     row.names = FALSE)
  }
  if (length(nospa_by_task)) {
    nd <- file.path(od, "results-staar-nospa"); dir.create(nd)
    for (id in names(nospa_by_task))
      utils::write.csv(nospa_by_task[[id]],
                       file.path(nd, sprintf("staar_coding_%s.csv", id)),
                       row.names = FALSE)
  }
  if (length(spa_by_task)) {
    sd <- file.path(od, "results-staar-spa"); dir.create(sd)
    for (id in names(spa_by_task))
      utils::write.csv(spa_by_task[[id]],
                       file.path(sd, sprintf("staar_coding_%s.csv", id)),
                       row.names = FALSE)
  }
  od
}

# A minimal genes_info: its 2nd column is each gene's chromosome (mirroring
# STAARpipeline::genes_info), one gene per autosome so .acr_expected_arrayids()'s
# chr-indexed cumsum is well-defined for any chroms in 1:22. Lets the native-
# completeness path run without the STAARpipeline package.
.mk_genes_info <- function() data.frame(gene = paste0("G", 1:22), chr = 1:22)


test_that("3-source assembly renames *_glowG and fronts the gene key", {
  # One GLOW chr (gene-grouped input not sorted gene-first), no native.
  glow22 <- .mk_glow_chr("22",
    gene = c("BBB", "AAA"), start = c(200, 100),
    category = c("missense", "missense"), glow_omni = c(0.3, 0.4))
  od <- .mk_run_dir(list("22" = glow22))

  agg <- aggregate_coding_results(od, chroms = 22L, gene_num_in_array = 500L,
                                  staar_modes = character(0))

  # STAAR_O / ACAT_O are renamed; the originals are gone.
  expect_true(all(c("STAAR_O_glowG", "ACAT_O_glowG") %in% names(agg$glow_coding_all)))
  expect_false(any(c("STAAR_O", "ACAT_O") %in% names(agg$glow_coding_all)))
  # glow_coding_all leads with label, chr, category.
  expect_equal(names(agg$glow_coding_all)[1:3], c("label", "chr", "category"))
  # merged uses `gene` (renamed from label), gene-first by start (AAA at 100).
  expect_equal(agg$coding_merged$gene[1], "AAA")
  expect_equal(agg$coding_merged$gene[2], "BBB")
})


test_that("gene-first OUTPUT sort orders coding_merged by (chr, gene_start, category)", {
  # Two chrs; within chr1 two genes at different starts, each with 2 categories.
  # Input is category-grouped (all missense, then all synonymous) -> the function
  # must re-sort to (chr, start, category).
  glow1 <- .mk_glow_chr("1",
    gene = c("GENEB", "GENEA", "GENEB", "GENEA"),
    start = c(500, 100, 500, 100),
    category = c("missense", "missense", "synonymous", "synonymous"),
    glow_omni = c(0.2, 0.3, 0.4, 0.5))
  glow2 <- .mk_glow_chr("2",
    gene = "GENEC", start = 50, category = "missense", glow_omni = 0.6)
  od <- .mk_run_dir(list("1" = glow1, "2" = glow2))

  agg <- aggregate_coding_results(od, chroms = c(1L, 2L), gene_num_in_array = 500L,
                                  staar_modes = character(0))
  m <- agg$coding_merged
  # chr1 GENEA (start 100) block first (its categories alpha-sorted), then GENEB
  # (start 500), then chr2.
  expect_equal(paste(m$chr, m$gene, m$category),
               c("1 GENEA missense", "1 GENEA synonymous",
                 "1 GENEB missense", "1 GENEB synonymous",
                 "2 GENEC missense"))
})


test_that("native left join keeps every GLOW cell; native subset-of GLOW", {
  glow22 <- .mk_glow_chr("22",
    gene = c("AAA", "BBB"), start = c(100, 200),
    category = c("missense", "missense"), glow_omni = c(1e-7, 0.5))
  # native (>= 2-variant) only for AAA; suffix columns _native / _spa.
  nospa <- list("0001" = .mk_native_task("22", "AAA", "missense",
                                         "STAAR_O_native", 1e-6))
  spa   <- list("0001" = .mk_native_task("22", "AAA", "missense",
                                         "STAAR_B_native", 1e-5))
  od <- .mk_run_dir(list("22" = glow22), nospa, spa)

  agg <- aggregate_coding_results(od, chroms = 22L, gene_num_in_array = 500L,
                                  staar_modes = c("spa", "nospa"),
                                  genes_info = .mk_genes_info(),
                                  allow_partial = TRUE)
  m <- agg$coding_merged
  expect_equal(nrow(m), 2L)                       # both GLOW cells kept
  # native columns present; AAA has values, BBB is NA (native subset-of GLOW).
  expect_true(all(c("STAAR_O_native", "STAAR_B_native") %in% names(m)))
  expect_equal(m$STAAR_O_native[m$gene == "AAA"], 1e-6)
  expect_true(is.na(m$STAAR_O_native[m$gene == "BBB"]))
  expect_equal(m$STAAR_B_native[m$gene == "AAA"], 1e-5)
  # the four-method method set is present.
  expect_setequal(agg$methods,
    c("GLOW_Omni", "STAAR_O_glowG", "STAAR_O_native", "STAAR_B_native"))
})


test_that("per-category significance is computed from the data (Bonferroni + GW)", {
  # 4 missense cells; one GLOW_Omni below alpha/4 and below the GW threshold.
  glow22 <- .mk_glow_chr("22",
    gene = c("G1", "G2", "G3", "G4"), start = c(100, 200, 300, 400),
    category = rep("missense", 4),
    glow_omni = c(1e-8, 0.5, 0.5, 0.5))
  od <- .mk_run_dir(list("22" = glow22))

  agg <- aggregate_coding_results(od, chroms = 22L, gene_num_in_array = 500L,
                                  staar_modes = character(0),
                                  alpha = 0.05, gw_alpha = 2.5e-6)
  sig <- agg$coding_sig_by_category
  row <- sig[sig$category == "missense" & sig$method == "GLOW_Omni", ]
  expect_equal(row$n_tested, 4L)
  expect_equal(row$bonferroni_threshold, 0.05 / 4)
  expect_equal(row$n_sig_bonferroni, 1L)
  expect_equal(row$genomewide_threshold, 2.5e-6)
  expect_equal(row$n_sig_genomewide, 1L)
  # the significant-cells table has exactly that one cell, flagged by GLOW_Omni.
  cells <- agg$coding_sig_by_category_cells
  expect_equal(nrow(cells), 1L)
  expect_equal(cells$gene, "G1")
  expect_true(grepl("GLOW_Omni", cells$sig_methods_bonf))
  expect_true(grepl("GLOW_Omni", cells$sig_methods_gw))
})


test_that("significant-cells table is gene-first and carries the count chain + rank", {
  # Two significant genes on different chrs at different starts.
  glow <- rbind(
    .mk_glow_chr("2", gene = "LATER", start = 9000, category = "missense",
                 glow_omni = 1e-9),
    .mk_glow_chr("1", gene = "EARLY", start = 100,  category = "missense",
                 glow_omni = 1e-9))
  od <- .mk_run_dir(list("1" = glow[glow$chr == "1", ],
                         "2" = glow[glow$chr == "2", ]))
  agg <- aggregate_coding_results(od, chroms = c(1L, 2L),
                                  gene_num_in_array = 500L,
                                  staar_modes = character(0))
  cells <- agg$coding_sig_by_category_cells
  expect_equal(nrow(cells), 2L)
  # gene-first: chr1 EARLY before chr2 LATER.
  expect_equal(cells$gene, c("EARLY", "LATER"))
  # count chain + rank columns are fronted.
  expect_true(all(c("gene", "chr", "start", "category", "n_variants",
                    "n_after_annotation", "rank", "sig_methods_bonf",
                    "sig_methods_gw") %in% names(cells)))
})


test_that("completeness reconciles GLOW chrs; missing chr stops unless allow_partial", {
  glow1 <- .mk_glow_chr("1", gene = "AAA", start = 100, category = "missense",
                        glow_omni = 0.5)
  od <- .mk_run_dir(list("1" = glow1))     # chr2 absent

  # chroms = c(1, 2): chr2 GLOW file missing -> stop by default.
  expect_error(
    aggregate_coding_results(od, chroms = c(1L, 2L), gene_num_in_array = 500L,
                             staar_modes = character(0)),
    "Missing GLOW per-chr")
  # allow_partial downgrades to a message and proceeds on chr1.
  msgs <- character(0)
  agg <- withCallingHandlers(
    aggregate_coding_results(od, chroms = c(1L, 2L), gene_num_in_array = 500L,
                             staar_modes = character(0), allow_partial = TRUE),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_true(any(grepl("allow_partial", msgs)))
  expect_equal(nrow(agg$coding_merged), 1L)
  expect_true(any(grepl("missing 1", agg$completeness)))
})


test_that("native-only cell in a NON-GLOW category is reported as info, not a mismatch", {
  # GLOW tests only `missense`; native (allow_partial) emits a `synonymous` cell
  # with no GLOW row -> category-aware INFO (expected; native emits all masks).
  glow22 <- .mk_glow_chr("22", gene = "AAA", start = 100, category = "missense",
                         glow_omni = 0.5)
  nospa <- list("0001" = rbind(
    .mk_native_task("22", "AAA", "missense",   "STAAR_O_native", 0.4),
    .mk_native_task("22", "ZZZ", "synonymous", "STAAR_O_native", 0.3)))
  od <- .mk_run_dir(list("22" = glow22), nospa)

  msgs <- character(0)
  withCallingHandlers(
    agg <- aggregate_coding_results(od, chroms = 22L, gene_num_in_array = 500L,
                                    staar_modes = "nospa",
                                    genes_info = .mk_genes_info(), allow_partial = TRUE),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  joined <- paste(msgs, collapse = " ")
  # The synonymous-only native cell is INFO (expected), not "check extraction".
  expect_true(grepl("NOT in this run's GLOW set", joined))
  expect_false(grepl("check extraction", joined))
  # The GLOW missense cell still merges; the synonymous native cell is dropped
  # (native subset-of GLOW; left join from GLOW).
  expect_equal(nrow(agg$coding_merged), 1L)
  expect_equal(agg$coding_merged$gene, "AAA")
})


test_that("missing output_dir / unknown primary_test error", {
  expect_error(
    aggregate_coding_results("/no/such/dir", chroms = 1L, gene_num_in_array = 500L),
    "dir.exists")
  glow22 <- .mk_glow_chr("22", gene = "AAA", start = 100, category = "missense",
                         glow_omni = 0.5)
  od <- .mk_run_dir(list("22" = glow22))
  expect_error(
    aggregate_coding_results(od, chroms = 22L, gene_num_in_array = 500L,
                             staar_modes = character(0), primary_test = "NOPE"),
    "primary-test column")
})


test_that("result_prefix and type_label route the non-coding files and name the tables", {
  g22 <- .mk_glow_chr(22, c("A", "B", "A"), c(100L, 500L, 10000L),
                      c("upstream", "upstream", "enhancer_CAGE"), c(1e-3, 0.2, 0.05))
  od <- tempfile("acr_nc_"); dir.create(file.path(od, "results-glow"), recursive = TRUE)
  utils::write.csv(g22, file.path(od, "results-glow", "glow_noncoding_chr22.csv"), row.names = FALSE)
  on.exit(unlink(od, recursive = TRUE), add = TRUE)
  agg <- aggregate_coding_results(od, chroms = 22L, gene_num_in_array = 500L,
                                  staar_modes = character(0), result_prefix = "glow_noncoding",
                                  type_label = "noncoding")
  expect_true(all(c("glow_noncoding_all", "noncoding_merged", "noncoding_sig_by_category",
                    "noncoding_sig_by_category_cells", "thresholds", "completeness_table")
                  %in% names(agg)))
  expect_false("coding_merged" %in% names(agg))
  expect_equal(nrow(agg$noncoding_merged), 3L)
  # The coding default still reads glow_coding_chr*.csv and names its tables as before.
  expect_error(aggregate_coding_results(od, chroms = 22L, gene_num_in_array = 500L,
                                        staar_modes = character(0)), "glow_coding_chr22")
})

test_that("the gene anchor keeps a gene's masks together when their spans differ", {
  # The review's counterexample: A upstream at 100, B upstream at 500, A
  # enhancer at 10,000. A (start, category) sort separates A's two masks.
  g22 <- .mk_glow_chr(22, c("A", "B", "A"), c(100L, 500L, 10000L),
                      c("upstream", "upstream", "enhancer_CAGE"), c(1e-3, 0.2, 0.05))
  od <- tempfile("acr_nc_"); dir.create(file.path(od, "results-glow"), recursive = TRUE)
  utils::write.csv(g22, file.path(od, "results-glow", "glow_noncoding_chr22.csv"), row.names = FALSE)
  on.exit(unlink(od, recursive = TRUE), add = TRUE)
  anchor <- data.frame(gene = c("A", "B"), chr = "22", anchor = c(100L, 500L), stringsAsFactors = FALSE)
  agg <- aggregate_coding_results(od, chroms = 22L, gene_num_in_array = 500L,
                                  staar_modes = character(0), result_prefix = "glow_noncoding",
                                  type_label = "noncoding", gene_anchor = anchor)
  m <- agg$noncoding_merged
  expect_equal(m$gene, c("A", "A", "B"))
  expect_equal(m$category, c("enhancer_CAGE", "upstream", "upstream"))
  # Without the anchor the old (start, category) order separates them.
  agg0 <- aggregate_coding_results(od, chroms = 22L, gene_num_in_array = 500L,
                                   staar_modes = character(0), result_prefix = "glow_noncoding",
                                   type_label = "noncoding")
  expect_equal(agg0$noncoding_merged$gene, c("A", "B", "A"))
  # noncoding_gene_anchor() on a table: the bundled table's start for a
  # protein-coding gene, the smallest start otherwise.
  tbl <- data.frame(label = c("ZZZ_not_a_gene", "ZZZ_not_a_gene", "RABGGTB"), chr = c("1", "1", "1"),
                    start = c(500L, 100L, 999L), stringsAsFactors = FALSE)
  an <- noncoding_gene_anchor(tbl)
  expect_equal(an$anchor[an$gene == "ZZZ_not_a_gene"], 100L)
  expect_equal(an$family[an$gene == "ZZZ_not_a_gene"], "other")
  expect_equal(an$family[an$gene == "RABGGTB"], "protein_coding")
  expect_equal(an$anchor[an$gene == "RABGGTB"], 75786197L)
})

test_that("the family threshold, its denominators, the reference constants and the completeness table", {
  # Unequal mask sizes; one cell with a missing comparator value (STAAR_O NA).
  g21 <- .mk_glow_chr(21, c("A", "B", "C", "A"), c(1L, 2L, 3L, 4L),
                      c("upstream", "upstream", "upstream", "UTR"),
                      c(1e-9, 0.5, 0.5, 0.5), staar_o = c(0.1, NA, 0.3, 0.4))
  g22 <- .mk_glow_chr(22, c("D", "D"), c(10L, 20L), c("upstream", "ncRNA"), c(0.5, 1e-8))
  od <- tempfile("acr_nc_"); dir.create(file.path(od, "results-glow"), recursive = TRUE)
  utils::write.csv(g21, file.path(od, "results-glow", "glow_noncoding_chr21.csv"), row.names = FALSE)
  utils::write.csv(g22, file.path(od, "results-glow", "glow_noncoding_chr22.csv"), row.names = FALSE)
  utils::write.csv(data.frame(category = "UTR", label = "D", message = "planted"),
                   file.path(od, "results-glow", "glow_noncoding_chr22_errors.csv"), row.names = FALSE)
  on.exit(unlink(od, recursive = TRUE), add = TRUE)
  refs <- data.frame(category = c("upstream", "ncRNA", "all"), cutoff = c(3.57e-7, 2.5e-6, 3.57e-7),
                     derivation = c("0.05 / (20,000 x 7)", "0.05 / 20,000", "0.05 / (20,000 x 7)"),
                     stringsAsFactors = FALSE)
  agg <- aggregate_coding_results(od, chroms = 21:22, gene_num_in_array = 500L,
                                  staar_modes = character(0), result_prefix = "glow_noncoding",
                                  type_label = "noncoding", family_alpha = 0.05,
                                  reference_cutoffs = refs)
  thr <- agg$thresholds
  # Per mask: alpha / n_tested; across all: family_alpha / n_all = 0.05 / 6 for GLOW,
  # 0.05 / 5 for the comparator (one NA).
  expect_equal(thr$cutoff[thr$method == "GLOW_Omni" & thr$category == "upstream"], 0.05 / 4)
  expect_equal(thr$n_tested[thr$method == "GLOW_Omni" & thr$category == "all"], 6L)
  expect_equal(thr$cutoff[thr$method == "GLOW_Omni" & thr$category == "all"], 0.05 / 6)
  expect_equal(thr$n_tested[thr$method == "STAAR_O_glowG" & thr$category == "all"], 5L)
  expect_equal(thr$cutoff[thr$method == "STAAR_O_glowG" & thr$category == "all"], 0.05 / 5)
  expect_match(thr$cutoff_rule[thr$method == "GLOW_Omni" & thr$category == "all"], "family_alpha / n_all")
  expect_equal(thr$reference_cutoff[thr$method == "GLOW_Omni" & thr$category == "upstream"], 3.57e-7)
  expect_equal(thr$reference_derivation[thr$method == "GLOW_Omni" & thr$category == "ncRNA"], "0.05 / 20,000")
  expect_true(is.na(thr$reference_cutoff[thr$method == "GLOW_Omni" & thr$category == "UTR"]))
  # The computed across-family cutoff is what flags the cells: both 1e-9 and 1e-8 clear 0.05 / 6.
  sig <- agg$noncoding_sig_by_category
  expect_equal(sig$genomewide_threshold[sig$method == "GLOW_Omni" & sig$category == "upstream"], 0.05 / 6)
  expect_equal(sum(sig$n_sig_genomewide[sig$method == "GLOW_Omni"]), 2L)
  # Without family_alpha the fixed gw_alpha stays, as for coding.
  agg0 <- aggregate_coding_results(od, chroms = 21:22, gene_num_in_array = 500L,
                                   staar_modes = character(0), result_prefix = "glow_noncoding",
                                   type_label = "noncoding")
  thr0 <- agg0$thresholds
  expect_equal(thr0$cutoff[thr0$method == "GLOW_Omni" & thr0$category == "all"], 2.5e-6)
  expect_match(thr0$cutoff_rule[thr0$method == "GLOW_Omni" & thr0$category == "all"], "fixed gw_alpha")
  # Completeness: per chromosome, the cells, genes, tested and missing primary
  # p-values, and the error rows (NA where no errors file exists).
  ct <- agg$completeness_table
  expect_equal(ct$chr, c("21", "22"))
  expect_equal(ct$n_cells, c(4L, 2L))
  expect_equal(ct$n_genes, c(3L, 1L))
  expect_equal(ct$n_tested_primary, c(4L, 2L))
  expect_equal(ct$n_error_rows, c(0L, 1L))
  # A missing chromosome still stops (the family is never defined on a partial run).
  expect_error(aggregate_coding_results(od, chroms = 20:22, gene_num_in_array = 500L,
                                        staar_modes = character(0), result_prefix = "glow_noncoding",
                                        type_label = "noncoding", family_alpha = 0.05),
               "Missing GLOW")
})

# ---- 2026-10-16: the implementation review round (D3, D6) ----

test_that("without a gene anchor the coding order is (chr, start, category), tied starts included", {
  g1 <- .mk_glow_chr("1", gene = c("A", "A", "B", "B"), start = c(100L, 100L, 100L, 100L),
                     category = c("missense", "plof", "missense", "plof"),
                     glow_omni = c(1e-9, 1e-9, 1e-9, 1e-9))
  od <- .mk_run_dir(list("1" = g1))
  on.exit(unlink(od, recursive = TRUE), add = TRUE)
  agg <- aggregate_coding_results(od, chroms = 1L, gene_num_in_array = 500L,
                                  staar_modes = character(0))
  key <- function(d) paste(d$gene, d$category)
  expect_equal(key(agg$coding_merged), c("A missense", "B missense", "A plof", "B plof"))
  expect_equal(key(agg$coding_sig_by_category_cells),
               c("A missense", "B missense", "A plof", "B plof"))
  # With a gene anchor the gene comes before the category.
  anchor <- data.frame(gene = c("A", "B"), chr = "1", anchor = c(100L, 100L),
                       stringsAsFactors = FALSE)
  agg2 <- aggregate_coding_results(od, chroms = 1L, gene_num_in_array = 500L,
                                   staar_modes = character(0), gene_anchor = anchor)
  expect_equal(key(agg2$coding_merged), c("A missense", "A plof", "B missense", "B plof"))
})

test_that("the completeness table lists every requested chromosome with its status and the sidecar tallies", {
  g21 <- .mk_glow_chr(21, c("A", "B"), c(1L, 2L), c("upstream", "UTR"), c(0.5, 0.4))
  g22 <- .mk_glow_chr(22, "X", 1L, "upstream", 0.5)[0, ]                      # all skipped or errored
  od <- tempfile("acr_cmp_"); gdir <- file.path(od, "results-glow")
  dir.create(gdir, recursive = TRUE)
  on.exit(unlink(od, recursive = TRUE), add = TRUE)
  utils::write.csv(g21, file.path(gdir, "glow_noncoding_chr21.csv"), row.names = FALSE)
  utils::write.csv(g22, file.path(gdir, "glow_noncoding_chr22.csv"), row.names = FALSE)
  # A zero-row errors sidecar on 21 (a clean unit), two error rows on 22.
  utils::write.csv(data.frame(category = character(0), label = character(0), message = character(0)),
                   file.path(gdir, "glow_noncoding_chr21_errors.csv"), row.names = FALSE)
  utils::write.csv(data.frame(category = c("upstream", "UTR"), label = "C", message = "planted"),
                   file.path(gdir, "glow_noncoding_chr22_errors.csv"), row.names = FALSE)
  sm <- function(cat, att, ok, skip, err)
    data.frame(category = cat, n_attempted = att, n_ok = ok, n_skip_empty = skip,
               n_staar_skip_single = 0L, n_error = err, runtime_secs = 1)
  utils::write.csv(rbind(sm("upstream", 3L, 1L, 2L, 0L), sm("UTR", 1L, 1L, 0L, 0L)),
                   file.path(gdir, "glow_noncoding_chr21_summary.csv"), row.names = FALSE)
  utils::write.csv(rbind(sm("upstream", 2L, 0L, 1L, 1L), sm("UTR", 1L, 0L, 0L, 1L)),
                   file.path(gdir, "glow_noncoding_chr22_summary.csv"), row.names = FALSE)
  msgs <- character(0)
  agg <- withCallingHandlers(
    aggregate_coding_results(od, chroms = 20:22, gene_num_in_array = 500L,
                             staar_modes = character(0), result_prefix = "glow_noncoding",
                             type_label = "noncoding", family_alpha = 0.05, allow_partial = TRUE),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  ct <- agg$completeness_table
  expect_equal(ct$chr, c("20", "21", "22"))
  expect_equal(ct$status, c("no results file", "present", "empty results file"))
  expect_equal(ct$n_cells, c(0L, 2L, 0L))
  expect_equal(ct$n_attempted, c(NA_integer_, 4L, 3L))
  expect_equal(ct$n_ok, c(NA_integer_, 2L, 0L))
  expect_equal(ct$n_skip_empty, c(NA_integer_, 2L, 1L))
  expect_equal(ct$n_errored, c(NA_integer_, 0L, 2L))
  expect_equal(ct$n_error_rows, c(0L, 0L, 2L))
  # Without any sidecar (the coding route) the tallies and the error rows are NA.
  od2 <- .mk_run_dir(list("1" = .mk_glow_chr("1", "A", 100L, "missense", 0.5)))
  on.exit(unlink(od2, recursive = TRUE), add = TRUE)
  ct2 <- aggregate_coding_results(od2, chroms = 1L, gene_num_in_array = 500L,
                                  staar_modes = character(0))$completeness_table
  expect_equal(ct2$status, "present")
  expect_true(is.na(ct2$n_attempted) && is.na(ct2$n_error_rows))
})
