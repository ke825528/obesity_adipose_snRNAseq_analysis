#############################################
# Script:    02_integration_and_clustering.R
# Purpose:   Standard pre-processing of the merged object: normalisation,
#            variable feature selection, scaling, PCA, Harmony batch correction,
#            UMAP embedding and graph-based clustering.
# Related manuscript figures: n/a (standard workflow, no figure-specific code)
# Required R packages: Seurat, harmony, ggplot2
#############################################
#
# Input
#   results/01_merged_qc.rds   merged, QC-filtered object written by 01
#
# Main steps
#   1. NormalizeData -> FindVariableFeatures -> ScaleData
#   2. RunPCA, ElbowPlot (diagnostic for the number of PCs)
#   3. RunHarmony on the `batch` column created in 01
#   4. RunUMAP / FindNeighbors / FindClusters on the Harmony embedding
#   5. JoinLayers, so that the marker analysis in 03 sees a single layer
#
# Output
#   results/02_clustered.rds                  normalised, integrated, clustered
#   figures/02_elbow_plot.pdf
#   figures/02_umap_by_sample.pdf
#   figures/02_umap_by_cluster.pdf

source(here::here("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
})

combined <- readRDS(file.path(results_dir, "01_merged_qc.rds"))

# ---------------------------------------------------------------------------
# 1. Normalisation, variable features, scaling
# ---------------------------------------------------------------------------
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, selection.method = hvg_method, nfeatures = hvg_n)
combined <- ScaleData(combined)

# ---------------------------------------------------------------------------
# 2. PCA
# ---------------------------------------------------------------------------
combined <- RunPCA(combined, npcs = max(pca_dims))

p_elbow <- ElbowPlot(combined, ndims = max(pca_dims), reduction = "pca")
ggsave(file.path(figures_dir, "02_elbow_plot.pdf"), p_elbow, width = 7, height = 5)

# ---------------------------------------------------------------------------
# 3. Harmony batch correction
# ---------------------------------------------------------------------------
# dims.use is passed explicitly so that the embedding used downstream is the one
# selected by `pca_dims` and not Harmony's own default.
combined <- RunHarmony(combined, group.by.vars = harmony_group_by,
                       dims.use = pca_dims)

# ---------------------------------------------------------------------------
# 4. UMAP and clustering on the Harmony embedding
# ---------------------------------------------------------------------------
combined <- RunUMAP(combined, reduction = "harmony", dims = pca_dims)
combined <- FindNeighbors(combined, reduction = "harmony", dims = pca_dims)
combined <- FindClusters(combined, resolution = cluster_resolution)

ggsave(file.path(figures_dir, "02_umap_by_sample.pdf"),
       DimPlot(combined, reduction = "umap", group.by = "orig.ident"),
       width = 8, height = 6)
ggsave(file.path(figures_dir, "02_umap_by_cluster.pdf"),
       DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE),
       width = 8, height = 6)

write.csv(data.frame(cell = colnames(combined),
                     cluster = as.character(combined$seurat_clusters),
                     sample = combined$orig.ident,
                     group = combined$sample_type),
          file = file.path(results_tables_dir, "02_cell_cluster_assignment.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 5. Collapse the per-sample layers before the marker analysis
# ---------------------------------------------------------------------------
combined <- JoinLayers(combined)

saveRDS(combined, file = file.path(results_dir, "02_clustered.rds"))

message("02: ", length(levels(combined$seurat_clusters)), " clusters at resolution ",
        cluster_resolution)
