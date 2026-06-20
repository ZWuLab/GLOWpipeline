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
