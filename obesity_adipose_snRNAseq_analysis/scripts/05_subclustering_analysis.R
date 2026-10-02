#############################################
# Script:    05_subclustering_analysis.R
# Purpose:   Standard sub-clustering of one cell compartment: subset the cells of
#            interest, re-normalise and re-cluster them independently of the main
#            clustering, and characterise the resulting sub-clusters.
# Related manuscript figures: n/a (standard workflow, no figure-specific code)
# Required R packages: Seurat, dplyr, ggplot2
#############################################
#
# Input
#   results/03_annotated.rds   annotated object written by 03
#
# Main steps
#   1. Subset the compartment named in `subcluster_celltypes` (00_config.R).
#   2. Normalise, select variable features, scale.
#   3. PCA, ElbowPlot, UMAP, neighbourhood graph and clustering with the
#      sub-clustering parameters from 00_config.R.
#   4. Sub-cluster markers and a marker panel for the compartment.
#   5. Store neutral sub-cluster labels ("myeloid_0", "myeloid_1", ...) in
#      `subcluster`. Use the marker tables to replace them with biological names.
#
# Output
#   results/05_subclustered.rds
#   results/tables/05_subcluster_markers_full.csv
#   results/tables/05_subcluster_markers_top.csv
#   figures/05_subcluster_*.pdf

source(here::here("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

annotated <- readRDS(file.path(results_dir, "03_annotated.rds"))

# ---------------------------------------------------------------------------
# 1. Subset the compartment
# ---------------------------------------------------------------------------
cells <- colnames(annotated)[annotated$celltype %in% subcluster_celltypes]
if (!length(cells)) {
  stop("None of the cell types in `subcluster_celltypes` (",
       paste(subcluster_celltypes, collapse = ", "), ") is present in ",
       "results/03_annotated.rds. Adjust `subcluster_celltypes` in ",
       "scripts/00_config.R.")
}

sub <- subset(annotated, cells = cells)
Idents(sub) <- "celltype"
message("05: sub-clustering ", length(cells), " cells from ",
        paste(subcluster_celltypes, collapse = ", "))

# ---------------------------------------------------------------------------
# 2. Normalisation
# ---------------------------------------------------------------------------
sub <- NormalizeData(sub)
sub <- FindVariableFeatures(sub, selection.method = hvg_method,
                            nfeatures = subcluster_hvg_n)
sub <- ScaleData(sub)

# ---------------------------------------------------------------------------
# 3. PCA, UMAP and clustering
# ---------------------------------------------------------------------------
sub <- RunPCA(sub, npcs = max(subcluster_pca_dims))

ggsave(file.path(figures_dir, "05_subcluster_elbow_plot.pdf"),
       ElbowPlot(sub, ndims = max(subcluster_pca_dims), reduction = "pca"),
       width = 7, height = 5)

sub <- RunUMAP(sub, dims = subcluster_pca_dims)
sub <- FindNeighbors(sub, dims = subcluster_pca_dims)
sub <- FindClusters(sub, resolution = subcluster_resolution)

ggsave(file.path(figures_dir, "05_subcluster_umap_clusters.pdf"),
       DimPlot(sub, reduction = "umap", label = TRUE) + ggtitle("Sub-clusters"),
       width = 8, height = 6)
ggsave(file.path(figures_dir, "05_subcluster_umap_by_sample.pdf"),
       DimPlot(sub, reduction = "umap", group.by = "orig.ident"),
       width = 8, height = 6)
ggsave(file.path(figures_dir, "05_subcluster_umap_by_group.pdf"),
       DimPlot(sub, reduction = "umap", group.by = "sample_type"),
       width = 8, height = 6)

# ---------------------------------------------------------------------------
# 4. Sub-cluster markers
# ---------------------------------------------------------------------------
markers <- FindAllMarkers(
  object          = sub,
  test.use        = marker_test,
  only.pos        = marker_only_pos,
  logfc.threshold = marker_logfc_threshold
)

write.csv(markers, file = file.path(results_tables_dir, "05_subcluster_markers_full.csv"),
          row.names = FALSE)

top_markers <- markers %>%
  filter(p_val < marker_pval_cutoff) %>%
  select(gene, everything()) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = marker_top_n) %>%
  as.data.frame()

write.csv(top_markers, file = file.path(results_tables_dir, "05_subcluster_markers_top.csv"),
          row.names = FALSE)

# Marker panel for the compartment (myeloid markers; gene symbols, HGNC)
myeloid_markers <- c("CD14", "CD68", "FCGR3A", "MRC1", "LYVE1", "C1QA", "C1QB",
                     "APOE", "SPP1", "S100A8", "S100A9", "FCN1", "VCAN", "LYZ",
                     "CD163", "HLA-DRA", "CD1C", "FCER1A", "CLEC9A", "MARCO")

present <- intersect(myeloid_markers, rownames(sub))
missing <- setdiff(myeloid_markers, rownames(sub))
if (length(missing)) {
  message("  marker panel: gene(s) not detected and skipped: ",
          paste(missing, collapse = ", "))
}

if (length(present)) {
  pdf(file.path(figures_dir, "05_subcluster_marker_dotplot.pdf"), width = 14, height = 6)
  print(DotPlot(sub, features = present) + RotatedAxis() + ggtitle("Sub-cluster markers"))
  dev.off()
} else {
  message("  marker panel: no gene present, panel skipped")
}

# ---------------------------------------------------------------------------
# 5. Neutral sub-cluster labels
# ---------------------------------------------------------------------------
# Cluster labels are deliberately neutral: they identify the sub-clusters, but
# they do not claim a biological identity. Replace them with the names you
# assign from results/tables/05_subcluster_markers_top.csv.
sub$subcluster <- paste0(subcluster_label_prefix, as.character(sub$seurat_clusters))
Idents(sub) <- "subcluster"

write.csv(as.data.frame(table(subcluster = sub$subcluster,
                              sample = sub$orig.ident)),
          file = file.path(results_tables_dir, "05_subcluster_cell_counts.csv"),
          row.names = FALSE)

saveRDS(sub, file = file.path(results_dir, "05_subclustered.rds"))

message("05: ", length(levels(sub$seurat_clusters)),
        " sub-clusters at resolution ", subcluster_resolution)
