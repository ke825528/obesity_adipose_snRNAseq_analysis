#############################################
# Script:    00_config.R
# Purpose:   Single source of truth for the whole pipeline: directory layout,
#            random seed, sample sheet, group definitions, QC thresholds,
#            processing parameters, differential expression settings,
#            enrichment settings, CellChat settings and the genes of interest.
#            Sourced by every other script. Contains no analysis code.
# Related manuscript figures: n/a (standard workflow, no figure-specific code)
# Required R packages: here
#############################################
#
# Every parameter that affects a result is defined here and nowhere else, so the
# whole pipeline can be re-configured from this one file.
#
# Path layout:
#   data/     input data (not distributed with this repository)
#   results/  intermediate objects, result tables and marker gene tables
#   figures/  figures
#
# Scripts must be run in numerical order; each one sources this file and can be
# started from the repository root or from any of its sub-directories. All paths
# are resolved from the repository root, so the working directory is never
# changed anywhere in this repository.

suppressPackageStartupMessages(library(here))

project_dir <- here::here()

data_dir    <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(project_dir, "figures")

results_deg_dir      <- file.path(results_dir, "deg")
results_enrich_dir   <- file.path(results_dir, "enrichment")
results_cellchat_dir <- file.path(results_dir, "cellchat")
results_tables_dir   <- file.path(results_dir, "tables")

for (d in c(data_dir, results_dir, figures_dir, results_deg_dir,
            results_enrich_dir, results_cellchat_dir, results_tables_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ===========================================================================
# 1. Organism
# ===========================================================================
# The single switch for everything that depends on the species. The annotation
# database, the KEGG organism code, the MSigDB species, the mitochondrial gene
# prefix, the GSEA gene set file and the CellChat database are all derived from
# it below, so they cannot disagree with each other.
organism <- "human"   # "human" or "mouse"

if (organism == "human") {
  enrichment_orgdb  <- "org.Hs.eg.db"
  enrichment_kegg   <- "hsa"
  msigdb_species    <- "Homo sapiens"
  mito_pattern      <- "^MT-"
  gmt_filename      <- "c2.cp.kegg.v7.5.1.symbols.gmt"
  cellchat_db       <- "CellChatDB.human"
  cellchat_ppi      <- "PPI.human"
} else if (organism == "mouse") {
  enrichment_orgdb  <- "org.Mm.eg.db"
  enrichment_kegg   <- "mmu"
  msigdb_species    <- "Mus musculus"
  mito_pattern      <- "^mt-"
  gmt_filename      <- "m2.cp.kegg.v7.5.1.symbols.gmt"
  cellchat_db       <- "CellChatDB.mouse"
  cellchat_ppi      <- "PPI.mouse"
} else {
  stop("`organism` must be 'human' or 'mouse'")
}

# ===========================================================================
# 2. Reproducibility and resources
# ===========================================================================
seed <- 12345
set.seed(seed)

# Number of workers for the parallelised CellChat steps, and the maximum size
# allowed for objects exported to those workers. The default of 500 MiB is too
# small for a full expression matrix, so it is raised here.
n_workers        <- 4
future_max_bytes <- 12 * 1024^3

# ===========================================================================
# 3. Samples and groups
# ===========================================================================
# data/sample_sheet.csv maps an anonymised sample identifier to the 10x Genomics
# directory holding that sample and to the study group it belongs to. Edit this
# file to run the pipeline on other data; nothing else needs to change.
sample_sheet <- read.csv(file.path(data_dir, "sample_sheet.csv"),
                         stringsAsFactors = FALSE)

stopifnot(all(c("sample_id", "data_dir", "group") %in% colnames(sample_sheet)))
stopifnot(!anyDuplicated(sample_sheet$sample_id))

# The two groups compared throughout the pipeline. The first entry is the
# experimental group, the second the control group.
group_levels <- c("group_A", "group_B")

stopifnot(all(unique(sample_sheet$group) %in% group_levels))

# Absolute path of each sample's 10x directory, named by sample identifier.
sample_dirs <- setNames(file.path(data_dir, sample_sheet$data_dir),
                        sample_sheet$sample_id)

# Derive the group label of every cell from `orig.ident`. `orig.ident` is set to
# the sample identifier when the Seurat objects are created (01), so this is a
# well-defined 1:1 mapping between sample and group.
assign_sample_type <- function(orig_ident) {
  group <- sample_sheet$group[match(orig_ident, sample_sheet$sample_id)]
  ifelse(is.na(group), NA_character_, group)
}

# ===========================================================================
# 4. Quality control  (ONE set of thresholds for the whole repository)
# ===========================================================================
# These values are the only QC parameters in the repository. They are used by
# 01_import_and_qc.R and by no other script; downstream scripts never filter
# cells again.
qc_min_cells_per_gene    <- 3      # a gene must be detected in >= 3 cells
qc_min_features_per_cell <- 300    # a cell must be detected with >= 300 genes
qc_max_features_per_cell <- 7500   # ... and with < 7500 genes
qc_max_percent_mt        <- 30     # and have < 30% mitochondrial counts

# ===========================================================================
# 5. Standard processing parameters
# ===========================================================================
hvg_method         <- "vst"
hvg_n              <- 2000       # variable features used for scaling and PCA
pca_dims           <- 1:20       # PCs used for PCA, Harmony, the neighbourhood
                                 # graph, UMAP and the clustering
harmony_group_by   <- "batch"    # metadata column corrected by Harmony
cluster_resolution <- 1          # resolution of the main clustering

# Sub-clustering of one compartment (05). These values are independent of the
# main clustering above.
subcluster_celltypes    <- c("monocyte", "macrophage")  # cell types to subset
subcluster_hvg_n        <- 3000
subcluster_pca_dims     <- 1:20
subcluster_resolution   <- 0.5
subcluster_label_prefix <- "myeloid_"   # neutral sub-cluster labels, see 05

# Marker gene detection
marker_test            <- "wilcox"
marker_only_pos        <- TRUE
marker_logfc_threshold <- 0.25
marker_pval_cutoff     <- 0.05
marker_top_n           <- 10

# ===========================================================================
# 6. Differential expression between the two groups
# ===========================================================================
# `sample_type` is the metadata column holding the group of every cell.
de_group_column        <- "sample_type"
de_test                <- "wilcox"
de_logfc_threshold     <- 0.25
de_min_pct             <- 0.1
de_min_cells_per_group <- 3      # a comparison is skipped below this group size
# A gene is reported as differentially expressed if both hold:
de_padj_cutoff         <- 0.05
de_logfc_cutoff        <- 0.5

# ===========================================================================
# 7. Functional enrichment
# ===========================================================================
enrichment_go_ontologies     <- c("BP", "CC", "MF")
enrichment_pvalue_cutoff     <- 0.05
enrichment_qvalue_cutoff     <- 0.05
enrichment_padjust_method    <- "BH"
enrichment_show_category     <- 10
enrichment_gene_padj_cutoff  <- 0.05
enrichment_gene_logfc_cutoff <- 0.5
enrichment_min_genes         <- 10   # warn below this number of input genes

# GSEA / fgsea
gsea_min_gs_size    <- 10
gsea_pvalue_cutoff  <- 0.05
gsea_top_pathways   <- 10        # pathways plotted individually
gsea_min_overlap    <- 100       # warn below this overlap with the gene sets
fgsea_padj_cutoff   <- 0.05
fgsea_top_pathways  <- 20

# Cell type whose differential expression result is used for the enrichment
# analysis. Set to NULL to use the cell type with the most significant genes.
enrichment_celltype <- NULL

# MSigDB KEGG gene sets in GMT format (used by clusterProfiler::GSEA). The file
# name follows `organism`; download the matching file from
#   https://www.gsea-msigdb.org/gsea/msigdb/
# and place it in data/.
gsea_gmt_file <- file.path(data_dir, gmt_filename)

# ===========================================================================
# 8. CellChat
# ===========================================================================
# The database and the PPI network follow `organism`; see section 1.
cellchat_group_by  <- "celltype"   # metadata column with the cell type labels
cellchat_min_cells <- 10

# ===========================================================================
# 9. Genes of interest (descriptive visualisation only, see 08)
# ===========================================================================
# These genes are only plotted (FeaturePlot / VlnPlot / DotPlot); no analysis is
# performed on them. Change the vector to plot your own genes. Note that PGRN is
# the protein product of the GRN gene, hence the gene symbol "GRN" below.
genes_of_interest <- c("LGMN", "GRN", "SORT1")

# ===========================================================================
# 10. Helpers
# ===========================================================================
# Make a string safe to use as a file name (cell type names contain spaces).
sanitise_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  gsub("^_+|_+$", "", x)
}
