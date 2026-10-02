#############################################
# Script:    03_marker_analysis_and_annotation.R
# Purpose:   Standard cluster marker analysis and cell type annotation: cluster
#            markers, per-cell-type marker gene panels, and assignment of a
#            `celltype` label to every cluster.
# Related manuscript figures: n/a (standard workflow, no figure-specific code)
# Required R packages: Seurat, dplyr, ggplot2
#############################################
#
# Input
#   results/02_clustered.rds   clustered object written by 02
#
# Main steps
#   1. FindAllMarkers per cluster (parameters from 00_config.R).
#   2. Top marker genes per cluster -> results/tables/03_cluster_markers_top.csv
#   3. Marker gene panels -> figures/03_annotation_marker_dotplots.pdf
#   4. Apply the cluster -> cell type map below and write the `celltype` metadata
#      column. Clusters that are not in the map are labelled "unassigned".
#
# Output
#   results/03_annotated.rds                     object with `celltype`
#   results/tables/03_cluster_markers_full.csv
#   results/tables/03_cluster_markers_top.csv
#   results/tables/03_cluster_celltype_assignment.csv
#   figures/03_annotation_marker_dotplots.pdf
#   figures/03_umap_celltype.pdf
#
# NOTE on the annotation map: the cluster -> cell type assignment below refers
# to the clustering produced by 02 on THIS dataset. It is dataset-specific
# content, not a pipeline parameter. When the pipeline is run on other data,
# inspect results/tables/03_cluster_markers_top.csv and replace the map - the
# rest of the script does not need to change. Clusters that are absent from the
# map are labelled "unassigned" and are excluded from the differential expression
# in 04.

source(here::here("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

clustered <- readRDS(file.path(results_dir, "02_clustered.rds"))

# ---------------------------------------------------------------------------
# 1.-2. Cluster markers
# ---------------------------------------------------------------------------
markers <- FindAllMarkers(
  object          = clustered,
  test.use        = marker_test,
  only.pos        = marker_only_pos,
  logfc.threshold = marker_logfc_threshold
)

write.csv(markers, file = file.path(results_tables_dir, "03_cluster_markers_full.csv"),
          row.names = FALSE)

top_markers <- markers %>%
  filter(p_val < marker_pval_cutoff) %>%
  select(gene, everything()) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = marker_top_n) %>%
  as.data.frame()

write.csv(top_markers, file = file.path(results_tables_dir, "03_cluster_markers_top.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 3. Marker gene panels
# ---------------------------------------------------------------------------
# Gene symbols follow HGNC nomenclature: Seurat resolves features by gene symbol,
# so protein/CD designations such as "CD45" or "CD29" cannot be used here (CD45 is
# PTPRC, CD29 is ITGB1, CD73 is NT5E, CD146 is MCAM, ...).
marker_panels <- list(
  adipocyte      = c("ADIPOQ", "LEP", "PLIN1", "PLIN4", "RBP4", "FABP4"),
  adsc           = c("PDGFRA", "THY1", "ENG", "NT5E", "ITGB1", "CD44"),
  preadipocyte   = c("PDGFRA", "CD36", "PPARG", "CEBPA"),
  fibroblast     = c("COL1A1", "COL3A1", "DCN", "LUM"),
  endothelial    = c("PECAM1", "VWF", "CDH5", "CLDN5", "RAMP2"),
  lymphatic      = c("PROX1", "PDPN", "LYVE1"),
  pericyte       = c("RGS5", "PDGFRB", "MCAM", "CSPG4", "NOTCH3"),
  smooth_muscle  = c("ACTA2", "MYH11", "TAGLN"),
  epithelial     = c("EPCAM", "KRT8", "KRT18", "KRT19", "CLDN3"),
  mesothelial    = c("MSLN", "UPK3B", "KRT19", "WT1"),
  monocyte       = c("S100A8", "S100A9", "FCN1", "VCAN", "LYZ"),
  macrophage     = c("CD68", "CD14", "FCGR3A", "MRC1", "C1QA", "APOE"),
  dendritic      = c("HLA-DRA", "HLA-DRB1", "CD1C", "FCER1A", "CLEC9A"),
  mast           = c("CPA3", "KIT", "TPSAB1", "MS4A2"),
  nk             = c("NCAM1", "KLRD1", "KLRK1", "NKG7", "GNLY"),
  t_cell         = c("CD3D", "CD3E", "CD2", "IL7R", "CD8A"),
  b_cell         = c("MS4A1", "CD79A", "CD79B", "CD19", "CD22"),
  plasma         = c("JCHAIN", "MZB1", "IGHG3", "DERL3"),
  neuron         = c("PLP1", "S100B", "SOX10")
)

# Plot only the genes that are actually present in the object, so that a missing
# marker gene produces a message instead of an error.
plot_panel <- function(obj, genes, title) {
  present <- intersect(genes, rownames(obj))
  missing <- setdiff(genes, rownames(obj))
  if (length(missing)) {
    message("  panel '", title, "': gene(s) not detected and skipped: ",
            paste(missing, collapse = ", "))
  }
  if (!length(present)) {
    message("  panel '", title, "': no gene present, panel skipped")
    return(NULL)
  }
  DotPlot(obj, features = present) + RotatedAxis() + ggtitle(title)
}

pdf(file.path(figures_dir, "03_annotation_marker_dotplots.pdf"), width = 12, height = 6)
for (nm in names(marker_panels)) {
  p <- plot_panel(clustered, marker_panels[[nm]], nm)
  if (!is.null(p)) print(p)
}
dev.off()

# ---------------------------------------------------------------------------
# 4. Cluster -> cell type assignment
# ---------------------------------------------------------------------------
# Dataset-specific cell-type annotation.
# Re-evaluate this mapping when applying the pipeline to another dataset.
# It refers to the clusters produced by 02 on this dataset only; inspect
# results/tables/03_cluster_markers_top.csv and edit the table below. Clusters
# that are absent from the map are labelled "unassigned" and are excluded from
# the differential expression in 04.
celltype_map <- data.frame(
  cluster = c("0", "1", "2", "3", "4", "6", "7", "8", "9", "10",
              "12", "13", "14", "15", "16", "17", "18", "19", "20",
              "21", "22", "23", "24"),
  celltype = c("mesothelial cell", "adipocyte", "native B cell", "adipocyte",
               "fibroblast cell", "T cell", "monocyte", "Neuron", "adipocyte",
               "epithelial cell", "macrophage", "endothelial cell", "B cell",
               "mast cell", "preadipocyte", "adsc", "Vascular smooth muscle cell",
               "T cell", "NK cell", "adipocyte", "adipocyte", "fibroblast cell",
               "macrophage"),
  stringsAsFactors = FALSE
)

clustered$celltype <- celltype_map$celltype[
  match(as.character(clustered$seurat_clusters), celltype_map$cluster)
]
clustered$celltype[is.na(clustered$celltype)] <- "unassigned"

unassigned <- setdiff(as.character(clustered$seurat_clusters),
                      celltype_map$cluster)
if (length(unassigned)) {
  message("  cluster(s) not present in celltype_map, labelled 'unassigned': ",
          paste(unassigned, collapse = ", "))
}

write.csv(as.data.frame(table(cluster = clustered$seurat_clusters,
                              celltype = clustered$celltype)),
          file = file.path(results_tables_dir, "03_cluster_celltype_assignment.csv"),
          row.names = FALSE)

Idents(clustered) <- "celltype"

ggsave(file.path(figures_dir, "03_umap_celltype.pdf"),
       DimPlot(clustered, reduction = "umap", group.by = "celltype", label = TRUE) +
         ggtitle("Annotated cell types"),
       width = 10, height = 7)

saveRDS(clustered, file = file.path(results_dir, "03_annotated.rds"))

message("03: ", length(unique(clustered$celltype)), " cell type labels assigned")
