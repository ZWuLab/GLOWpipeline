#' Helper: Create Synthetic Test aGDS for Component 3 Testing
#'
#' Builds a small annotated GDS (aGDS) file with hand-coded annotations
#' so that every filtering outcome is deterministic and verifiable.
#'
#' File Log (reverse chronological order):
#' - 2026-10-15: Updated by Claude Code (Fable 5.1), prompted by ZWu. Planted
#'               the gene-linkage strings of the gene-centric non-coding masks
#'               v1 shapes of upstream/downstream records (two genes; the
#'               two-category "upstream;downstream" on 29:30), UTR records with
#'               transcript details (a multi-gene record on 76:78, "UTR5;UTR3"
#'               on 79:80), ncRNA records of bundled chromosome 22 ncRNA genes
#'               on 36:40, GeneHancer strings in the v1 layout (a header
#'               without targets on 69), and CAGE/DHS presence on records that
#'               also carry GeneHancer. Both copies changed identically.
#' - 2026-06-15: Updated by Claude Code (Opus 4.8 / 1M ctx), prompted by ZWu.
#'               create_test_agds() now SKIPs (was stop()) when a GDS
#'               Bioconductor package (gdsfmt/SeqArray/SNPRelate) is unavailable,
#'               so a CI runner without those binaries skips the GDS tests rather
#'               than failing R CMD check. Mirrors the same fix in the GLOWr copy.
#' - 2026-05-26: Copied verbatim from GLOWr (packages/GLOWr/tests/testthat/
#'               helper-create-test-agds.R) by Claude Code (Opus 4.7), prompted
#'               by ZWu, so GLOWpipeline's engine integration test is
#'               self-contained. Keep in sync with the GLOWr original; it is
#'               self-contained (only gdsfmt/SeqArray/SNPRelate, no GLOWr
#'               internals).
#' - 2026-02-20: Created by r-developer agent - Synthetic aGDS with hand-coded
#'               annotations for 125 variants across 5 genes on chr22

# ==============================================================================
# MAIN FUNCTION
# ==============================================================================

#' Create a Synthetic Test aGDS File
#'
#' Generates a small aGDS with 125 variants across 5 genes on chromosome 22,
#' 50 samples with MAF-based random genotypes, and fully specified FAVOR-style
#' annotations under \code{annotation/info/FunctionalAnnotation/}.
#'
#' The file is compatible with SeqArray and follows the node structure produced
#' by FAVORannotator (individual typed sub-nodes, not a single matrix).
#'
#' @param file_path Character. Output file path for the aGDS file.
#'
#' @return Invisible \code{file_path}. The file is created as a side effect.
#'
#' @details
#' \strong{Variant layout (125 total, all chr22):}
#' \itemize{
#'   \item GENE_A (40 variants, pos 1000-5000): mixed coding/noncoding
#'   \item GENE_B (40 variants, pos 10000-15000): mixed coding/noncoding
#'   \item GENE_C (30 variants, pos 20000-25000): all synonymous (negative ctrl)
#'   \item GENE_D (10 variants, pos 30000-32000): small gene, min-variant test
#'   \item GENE_E (5 variants, pos 40000-41000): rv_num_cutoff skip test
#' }
#'
#' \strong{Annotation nodes} are stored under
#' \code{annotation/info/FunctionalAnnotation/} as individual sub-nodes with
#' \code{compress = "LZMA_ra"}.
#'
#' \strong{QC filter:} \code{annotation/filter} is a character vector of
#' "PASS" or "" (empty string). Approximately 10 variants are non-PASS.
#'
#' Requires: SNPRelate, SeqArray, gdsfmt.
create_test_agds <- function(file_path) {

  # -- Check required packages. SKIP (not error) the calling test if a GDS
  #    package is unavailable -- e.g. a CI runner that could not install the
  #    Bioconductor binaries -- so a missing optional dependency never fails the
  #    suite. skip() called here propagates to the calling test (or skips the
  #    whole file when create_test_agds() is invoked at file top level). --
  for (pkg in c("SNPRelate", "SeqArray", "gdsfmt")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      testthat::skip(paste0("Package '", pkg, "' is required but not installed."))
    }
  }

  # ============================================================================
  # Step 1: Define all variant properties
  # ============================================================================

  set.seed(42)

  n_samples  <- 50L
  n_variants <- 125L

  # --- Positions ---
  # GENE_A: 40 variants in [1000, 5000]
  pos_A <- as.integer(seq(1000, 5000, length.out = 40))
  # GENE_B: 40 variants in [10000, 15000]
  pos_B <- as.integer(seq(10000, 15000, length.out = 40))
  # GENE_C: 30 variants in [20000, 25000]
  pos_C <- as.integer(seq(20000, 25000, length.out = 30))
  # GENE_D: 10 variants in [30000, 32000]
  pos_D <- as.integer(seq(30000, 32000, length.out = 10))
  # GENE_E: 5 variants in [40000, 41000]
  pos_E <- as.integer(seq(40000, 41000, length.out = 5))

  positions <- c(pos_A, pos_B, pos_C, pos_D, pos_E)
  stopifnot(length(positions) == n_variants)

  # --- MAF design ---
  # Most variants are rare (MAF < 0.01). ~10 are common (MAF > 0.05).
  # GENE_E: only 1 rare variant (rest are common, testing rv_num_cutoff skip).
  maf <- numeric(n_variants)

  # GENE_A (indices 1:40) -- 4 common, 36 rare
  maf[1:4]   <- runif(4, 0.06, 0.15)
  maf[5:40]  <- runif(36, 0.001, 0.009)

  # GENE_B (indices 41:80) -- 4 common, 36 rare
  maf[41:44] <- runif(4, 0.06, 0.15)
  maf[45:80] <- runif(36, 0.001, 0.009)

  # GENE_C (indices 81:110) -- 2 common, 28 rare
  maf[81:82]  <- runif(2, 0.06, 0.15)
  maf[83:110] <- runif(28, 0.001, 0.009)

  # GENE_D (indices 111:120) -- all rare

  maf[111:120] <- runif(10, 0.001, 0.009)

  # GENE_E (indices 121:125) -- 1 rare variant, 4 common
  maf[121]     <- 0.005
  maf[122:125] <- runif(4, 0.06, 0.15)

  # --- Annotation vectors (length = n_variants) ---
  # Initialize all as empty strings or zeros
  genecode_category  <- rep("", n_variants)
  exonic_category    <- rep("", n_variants)
  genecode_info      <- rep("", n_variants)
  exonic_info        <- rep("", n_variants)
  metasvm_pred       <- rep("", n_variants)
  genehancer         <- rep("", n_variants)
  cage_tc            <- rep("", n_variants)
  cage_promoter      <- rep("", n_variants)
  rdhs               <- rep("", n_variants)
  rsid               <- rep("", n_variants)

  cadd_phred              <- runif(n_variants, 0, 40)
  linsight                <- runif(n_variants, 0, 1)
  fathmm_xf               <- runif(n_variants, 0, 1)
  apc_conservation         <- runif(n_variants, 0, 50)
  apc_conservation_v2      <- runif(n_variants, 0, 50)
  apc_epigenetics          <- runif(n_variants, 0, 50)
  apc_epigenetics_active   <- runif(n_variants, 0, 50)
  apc_epigenetics_repressed <- runif(n_variants, 0, 50)
  apc_epigenetics_transcription <- runif(n_variants, 0, 50)
  apc_local_nucleotide_diversity <- runif(n_variants, 0, 50)
  apc_local_nucleotide_diversity_v3 <- runif(n_variants, 0, 50)
  apc_mappability          <- runif(n_variants, 0, 50)
  apc_mutation_density     <- runif(n_variants, 0, 50)
  apc_protein_function     <- runif(n_variants, 0, 50)
  apc_protein_function_v3  <- runif(n_variants, 0, 50)
  apc_proximity_to_tsstes  <- runif(n_variants, 0, 50)
  apc_transcription_factor <- runif(n_variants, 0, 50)

  # --- rsid for all variants ---
  rsid <- paste0("rs", seq(100001, length.out = n_variants))

  # ==========================================================================
  # GENE_A (indices 1:40): 5 stopgain, 5 nonsynonymous(D), 5 nonsynonymous(T),
  #   5 synonymous, 10 upstream, 10 intronic
  # ==========================================================================

  # Stopgain (indices 1:5)
  genecode_category[1:5]  <- "exonic"
  exonic_category[1:5]    <- "stopgain"
  genecode_info[1:5]      <- "GENE_A"
  exonic_info[1:5]        <- "GENE_A"

  # Nonsynonymous damaging (indices 6:10)
  genecode_category[6:10]  <- "exonic"
  exonic_category[6:10]    <- "nonsynonymous SNV"
  genecode_info[6:10]      <- "GENE_A"
  exonic_info[6:10]        <- "GENE_A"
  metasvm_pred[6:10]       <- "D"

  # Nonsynonymous tolerated (indices 11:15)
  genecode_category[11:15]  <- "exonic"
  exonic_category[11:15]    <- "nonsynonymous SNV"
  genecode_info[11:15]      <- "GENE_A"
  exonic_info[11:15]        <- "GENE_A"
  metasvm_pred[11:15]       <- "T"

  # Synonymous (indices 16:20)
  genecode_category[16:20]  <- "exonic"
  exonic_category[16:20]    <- "synonymous SNV"
  genecode_info[16:20]      <- "GENE_A"
  exonic_info[16:20]        <- "GENE_A"

  # Upstream (indices 21:30): the FAVOR v1 shape of a record upstream of two
  # genes ("G,G"); 29:30 are the two-category shape "G;G" with the category
  # "upstream;downstream" (upstream of GENE_A, downstream of GENE_B).
  genecode_category[21:30]  <- "upstream"
  genecode_info[21:30]      <- "GENE_A,GENE_B"
  genecode_category[29:30]  <- "upstream;downstream"
  genecode_info[29:30]      <- "GENE_A;GENE_B"

  # Intronic (indices 31:35); 36:40 are ncRNA records of genes present in the
  # bundled ncRNA table of chromosome 22 (POTEH-AS1, LINC01297), 40 with the
  # two-category shape "G;P" of an ncRNA_exonic;splicing record.
  genecode_category[31:40]  <- "intronic"
  genecode_info[31:40]      <- "GENE_A"
  genecode_category[36:39]  <- "ncRNA_exonic"
  genecode_info[36:39]      <- "POTEH-AS1"
  genecode_category[40]     <- "ncRNA_exonic;splicing"
  genecode_info[40]         <- "POTEH-AS1;LINC01297(ENST00000000003.1:exon2:c.10-2A>G)"

  # ==========================================================================
  # GENE_B (indices 41:80): 5 stoploss, 5 splicing, 5 frameshift, 5 missense,
  #   10 downstream, 10 UTR3
  # ==========================================================================

  # Stoploss (indices 41:45)
  genecode_category[41:45]  <- "exonic"
  exonic_category[41:45]    <- "stoploss"
  genecode_info[41:45]      <- "GENE_B"
  exonic_info[41:45]        <- "GENE_B"

  # Splicing (indices 46:50)
  genecode_category[46:50]  <- "splicing"
  genecode_info[46:50]      <- "GENE_B"

  # Frameshift deletion (indices 51:55)
  genecode_category[51:55]  <- "exonic"
  exonic_category[51:55]    <- "frameshift deletion"
  genecode_info[51:55]      <- "GENE_B"
  exonic_info[51:55]        <- "GENE_B"

  # Missense / nonsynonymous damaging (indices 56:60)
  genecode_category[56:60]  <- "exonic"
  exonic_category[56:60]    <- "nonsynonymous SNV"
  genecode_info[56:60]      <- "GENE_B"
  exonic_info[56:60]        <- "GENE_B"
  metasvm_pred[56:60]       <- "D"

  # Downstream (indices 61:70): downstream of GENE_B and GENE_C.
  genecode_category[61:70]  <- "downstream"
  genecode_info[61:70]      <- "GENE_B,GENE_C"

  # UTR (indices 71:80): the v1 shape "P" (gene with a transcript detail that
  # contains a comma) for 71:75; a single-category multi-gene UTR record
  # ("P,P") for 76:78; the two-category "UTR5;UTR3" shape ("P;P") for 79:80.
  genecode_category[71:80]  <- "UTR3"
  genecode_info[71:80]      <- "GENE_B(ENST00000000001.1:c.*100A>G,ENST00000000002.1:c.*95A>G)"
  genecode_info[76:78]      <- "GENE_B(ENST00000000001.1:c.*16G>A),GENE_C(ENST00000000004.1:c.*5C>T)"
  genecode_category[79:80]  <- "UTR5;UTR3"
  genecode_info[79:80]      <- "GENE_B(ENST00000000001.1:c.-10A>G);GENE_C(ENST00000000004.1:c.*7T>C)"

  # ==========================================================================
  # GENE_C (indices 81:110): All synonymous (negative control)
  # ==========================================================================
  genecode_category[81:110]  <- "exonic"
  exonic_category[81:110]    <- "synonymous SNV"
  genecode_info[81:110]      <- "GENE_C"
  exonic_info[81:110]        <- "GENE_C"

  # ==========================================================================
  # GENE_D (indices 111:120): 2 stopgain, 8 nonsynonymous
  # ==========================================================================

  # Stopgain (indices 111:112)
  genecode_category[111:112]  <- "exonic"
  exonic_category[111:112]    <- "stopgain"
  genecode_info[111:112]      <- "GENE_D"
  exonic_info[111:112]        <- "GENE_D"

  # Nonsynonymous damaging (indices 113:120)
  genecode_category[113:120]  <- "exonic"
  exonic_category[113:120]    <- "nonsynonymous SNV"
  genecode_info[113:120]      <- "GENE_D"
  exonic_info[113:120]        <- "GENE_D"
  metasvm_pred[113:120]       <- "D"

  # ==========================================================================
  # GENE_E (indices 121:125): 1 rare variant -- tests rv_num_cutoff skip
  # All exonic/nonsynonymous for simplicity
  # ==========================================================================
  genecode_category[121:125]  <- "exonic"
  exonic_category[121:125]    <- "nonsynonymous SNV"
  genecode_info[121:125]      <- "GENE_E"
  exonic_info[121:125]        <- "GENE_E"
  metasvm_pred[121:125]       <- "D"

  # --- Populate a few noncoding regulatory annotations for testing ---
  # Some upstream/downstream/intronic variants get enhancer/promoter/RDHS
  # values. GeneHancer strings have the v1 layout (header, then the connected
  # genes in the element's order); 69 is a header without targets.
  genehancer[21:23]  <- "Name=0.50;genehancer_id=GH22E000001;connected_gene=GENE_A;score=10.50;connected_gene=GENE_B;score=1.20"
  genehancer[61:63]  <- "Name=0.30;genehancer_id=GH22E000002;connected_gene=GENE_C;score=5.00;connected_gene=GENE_B;score=0.40"
  genehancer[69]     <- "Name=0.10;genehancer_id=GH22E000003"
  cage_tc[22:26]     <- "TC_chr22_1000"
  cage_tc[64:66]     <- "TC_chr22_2000"
  cage_tc[69]        <- "TC_chr22_2001"
  cage_promoter[27:29] <- "promoter_GENE_A"
  cage_promoter[67:69] <- "promoter_GENE_B"
  rdhs[21]           <- "RDHS_000"
  rdhs[30]           <- "RDHS_001"
  rdhs[61]           <- "RDHS_002"
  rdhs[70]           <- "RDHS_002"
  rdhs[31:35]        <- "RDHS_003"

  # --- QC Filter: ~10 variants set to non-PASS (empty string) ---
  filter_vals <- rep("PASS", n_variants)
  # Mark specific variants as non-PASS for QC testing
  # Spread across multiple genes: indices 3, 8, 18, 43, 53, 63, 73, 85, 115, 123
  non_pass_idx <- c(3L, 8L, 18L, 43L, 53L, 63L, 73L, 85L, 115L, 123L)
  filter_vals[non_pass_idx] <- ""

  # ============================================================================
  # Step 2: Generate genotypes from binomial(2, maf) per variant
  # ============================================================================

  # Genotype matrix: rows = variants, cols = samples (snpgdsCreateGeno format)
  geno <- matrix(0L, nrow = n_variants, ncol = n_samples)
  for (v in seq_len(n_variants)) {
    geno[v, ] <- rbinom(n_samples, size = 2, prob = maf[v])
  }

  # ============================================================================
  # Step 3: Create SNP GDS via snpgdsCreateGeno
  # ============================================================================

  tmp_snp <- tempfile(fileext = ".gds")
  on.exit(unlink(tmp_snp), add = TRUE)

  sample_ids <- sprintf("sample_%03d", seq_len(n_samples))
  snp_ids    <- seq_len(n_variants)

  # Generate REF/ALT allele pairs (A/G for odd indices, C/T for even)
  alleles <- ifelse(seq_len(n_variants) %% 2 == 1, "A/G", "C/T")

  SNPRelate::snpgdsCreateGeno(
    gds.fn         = tmp_snp,
    genmat         = geno,
    sample.id      = sample_ids,
    snp.id         = snp_ids,
    snp.chromosome = rep(22L, n_variants),
    snp.position   = positions,
    snp.allele     = alleles,
    snpfirstdim    = TRUE
  )

  # ============================================================================
  # Step 4: Convert to SeqArray GDS via seqSNP2GDS
  # ============================================================================

  suppressMessages({
    SeqArray::seqSNP2GDS(tmp_snp, file_path, verbose = FALSE)
  })

  # ============================================================================
  # Step 5: Open with gdsfmt write access and add annotation nodes
  # ============================================================================

  gds <- gdsfmt::openfn.gds(file_path, readonly = FALSE)
  on.exit(gdsfmt::closefn.gds(gds), add = TRUE)

  # -- Replace annotation/filter with character PASS/empty values --
  # seqSNP2GDS creates an Int32 factor filter; we need character for STAARpipeline
  filt_node <- gdsfmt::index.gdsn(gds, "annotation/filter")
  gdsfmt::delete.gdsn(filt_node)
  gdsfmt::add.gdsn(
    gdsfmt::index.gdsn(gds, "annotation"), "filter",
    val = filter_vals, compress = "LZMA_ra", closezip = TRUE
  )

  # -- Create annotation/info/FunctionalAnnotation folder --
  info_node <- gdsfmt::index.gdsn(gds, "annotation/info")
  fa_folder <- gdsfmt::addfolder.gdsn(info_node, "FunctionalAnnotation")

  # -- Helper to add a single annotation sub-node --
  .add_anno <- function(name, val) {
    gdsfmt::add.gdsn(fa_folder, name, val = val,
                      compress = "LZMA_ra", closezip = TRUE)
  }

  # -- String annotations --
  .add_anno("genecode_comprehensive_category",        genecode_category)
  .add_anno("genecode_comprehensive_exonic_category",  exonic_category)
  .add_anno("genecode_comprehensive_info",             genecode_info)
  .add_anno("genecode_comprehensive_exonic_info",      exonic_info)
  .add_anno("metasvm_pred",                            metasvm_pred)
  .add_anno("genehancer",                              genehancer)
  .add_anno("cage_tc",                                 cage_tc)
  .add_anno("cage_promoter",                           cage_promoter)
  .add_anno("rdhs",                                    rdhs)
  .add_anno("rsid",                                    rsid)

  # -- Numeric annotations --
  .add_anno("cadd_phred",                       cadd_phred)
  .add_anno("linsight",                         linsight)
  .add_anno("fathmm_xf",                        fathmm_xf)
  .add_anno("apc_conservation",                    apc_conservation)
  .add_anno("apc_conservation_v2",                 apc_conservation_v2)
  .add_anno("apc_epigenetics",                     apc_epigenetics)
  .add_anno("apc_epigenetics_active",              apc_epigenetics_active)
  .add_anno("apc_epigenetics_repressed",           apc_epigenetics_repressed)
  .add_anno("apc_epigenetics_transcription",       apc_epigenetics_transcription)
  .add_anno("apc_local_nucleotide_diversity",      apc_local_nucleotide_diversity)
  .add_anno("apc_local_nucleotide_diversity_v3",   apc_local_nucleotide_diversity_v3)
  .add_anno("apc_mappability",                     apc_mappability)
  .add_anno("apc_mutation_density",                apc_mutation_density)
  .add_anno("apc_protein_function",                apc_protein_function)
  .add_anno("apc_protein_function_v3",             apc_protein_function_v3)
  .add_anno("apc_proximity_to_tsstes",             apc_proximity_to_tsstes)
  .add_anno("apc_transcription_factor",            apc_transcription_factor)

  # ============================================================================
  # Step 6: Close GDS (handled by on.exit) and return path
  # ============================================================================

  invisible(file_path)
}
