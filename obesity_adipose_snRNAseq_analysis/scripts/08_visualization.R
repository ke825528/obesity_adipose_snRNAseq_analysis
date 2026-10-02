#############################################
# Script:    08_visualization.R
# Purpose:   Standard descriptive figures: cell type composition, UMAP embeddings
#            and expression of the genes named in `genes_of_interest`.
#            Descriptive only - no statistical analysis is performed here.
# Related manuscript figures: n/a (standard workflow, no figure-specific code)
# Required R packages: Seurat, ggplot2, ggpointdensity, viridis, patchwork
#############################################
#
# Input
#   results/03_annotated.rds   annotated object written by 03
#
# Main steps
#   1. Cell type composition per sample and per group.
#   2. UMAP coloured by cell type, by group and by sample.
#   3. A UMAP density rendering.
#   4. FeaturePlot / VlnPlot / DotPlot of the genes in `genes_of_interest`
#      (00_config.R). This is a descriptive panel only: it adds no analysis and
#      no interpretation.
#
# Output
#   figures/08_*.pdf
#   results/tables/08_celltype_composition_per_sample.csv
#   results/tables/08_celltype_composition_per_group.csv

source(here::here("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(ggpointdensity)
  library(viridis)
  library(patchwork)
})

annotated <- readRDS(file.path(results_dir, "03_annotated.rds"))
Idents(annotated) <- "celltype"

# A palette with exactly one colour per cell type, so that a change in the number
# of cell types cannot make the colour scale fail.
celltype_levels  <- sort(unique(as.character(annotated$celltype)))
celltype_colours <- setNames(grDevices::hcl.colors(length(celltype_levels), "Dark 3"),
                             celltype_levels)

# ---------------------------------------------------------------------------
# 1. Cell type composition
# ---------------------------------------------------------------------------
comp_sample <- as.data.frame(prop.table(table(celltype = annotated$celltype,
                                              sample   = annotated$orig.ident),
                                       margin = 2))
write.csv(comp_sample,
          file = file.path(results_tables_dir, "08_celltype_composition_per_sample.csv"),
          row.names = FALSE)

p_comp_sample <- ggplot(comp_sample, aes(x = sample, y = Freq, fill = celltype)) +
  geom_col(width = 0.8, colour = "black", linewidth = 0.2) +
  scale_fill_manual(values = celltype_colours) +
  labs(x = "Sample", y = "Proportion of cells", fill = "Cell type") +
  theme_classic()
ggsave(file.path(figures_dir, "08_celltype_composition_per_sample.pdf"), p_comp_sample,
       width = 10, height = 6)

comp_group <- as.data.frame(prop.table(table(celltype = annotated$celltype,
                                             group    = annotated$sample_type),
                                       margin = 2))
colnames(comp_group) <- c("celltype", "group", "proportion")
write.csv(comp_group,
          file = file.path(results_tables_dir, "08_celltype_composition_per_group.csv"),
          row.names = FALSE)

p_comp_group <- ggplot(comp_group, aes(x = celltype, y = proportion, fill = group)) +
  geom_col(position = "dodge", colour = "black", linewidth = 0.2) +
  labs(x = "Cell type", y = "Proportion of cells", fill = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(figures_dir, "08_celltype_composition_per_group.pdf"), p_comp_group,
       width = 10, height = 5)

# ---------------------------------------------------------------------------
# 2. UMAP embeddings
# ---------------------------------------------------------------------------
p_umap_celltype <- DimPlot(annotated, reduction = "umap", group.by = "celltype",
                           label = TRUE, label.size = 3) +
  ggtitle("Cell types")
ggsave(file.path(figures_dir, "08_umap_celltype.pdf"), p_umap_celltype, width = 10, height = 7)

p_umap_group <- DimPlot(annotated, reduction = "umap", group.by = "sample_type") +
  ggtitle("Study group")
ggsave(file.path(figures_dir, "08_umap_group.pdf"), p_umap_group, width = 8, height = 6)

p_umap_sample <- DimPlot(annotated, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Sample")
ggsave(file.path(figures_dir, "08_umap_sample.pdf"), p_umap_sample, width = 9, height = 6)

# ---------------------------------------------------------------------------
# 3. UMAP density
# ---------------------------------------------------------------------------
emb <- as.data.frame(Embeddings(annotated, "umap"))
colnames(emb) <- c("umap_1", "umap_2")

p_density <- ggplot(emb, aes(x = umap_1, y = umap_2)) +
  geom_pointdensity(size = 0.3, shape = 16) +
  viridis::scale_color_viridis(option = "A") +
  labs(x = "UMAP 1", y = "UMAP 2", colour = "Density") +
  theme_bw()
ggsave(file.path(figures_dir, "08_umap_density.pdf"), p_density, width = 8, height = 6)

# ---------------------------------------------------------------------------
# 4. Genes of interest (descriptive panel)
# ---------------------------------------------------------------------------
# `genes_of_interest` is set in 00_config.R. Only genes that are present in the
# object are plotted, so a missing gene produces a message instead of an error.
genes_present <- intersect(genes_of_interest, rownames(annotated))
genes_missing <- setdiff(genes_of_interest, rownames(annotated))

if (length(genes_missing)) {
  message("08: gene(s) of interest not detected and skipped: ",
          paste(genes_missing, collapse = ", "))
}

if (length(genes_present)) {
  p_feature <- FeaturePlot(annotated, features = genes_present, reduction = "umap",
                           pt.size = 0.1, ncol = min(3, length(genes_present)))
  ggsave(file.path(figures_dir, "08_genes_of_interest_featureplot.pdf"), p_feature,
         width = 5 * min(3, length(genes_present)), height = 5 * ceiling(length(genes_present) / 3))

  p_violin <- VlnPlot(annotated, features = genes_present, pt.size = 0,
                      ncol = min(3, length(genes_present))) &
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(figures_dir, "08_genes_of_interest_vlnplot.pdf"), p_violin,
         width = 6 * min(3, length(genes_present)), height = 5)

  p_dot <- DotPlot(annotated, features = genes_present) + RotatedAxis() +
    ggtitle("Genes of interest")
  ggsave(file.path(figures_dir, "08_genes_of_interest_dotplot.pdf"), p_dot,
         width = 8, height = 6)
} else {
  message("08: none of the genes in `genes_of_interest` is present; panel skipped")
}

message("08: figures written to ", figures_dir)
