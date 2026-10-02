#############################################
# Script:    06_cellchat_analysis.R
# Purpose:   Standard CellChat workflow: infer cell-cell communication between
#            the annotated cell types and visualise the resulting network.
#            This script is a generic template - it contains no comparison
#            between study groups and no cell-type-specific settings, so that it
#            can be pointed at any annotated Seurat object.
# Related manuscript figures: n/a (standard workflow, no figure-specific code)
# Required R packages: CellChat, Seurat, ComplexHeatmap, future
#############################################
#
# Input
#   results/03_annotated.rds   annotated object written by 03 (any annotated
#                              Seurat object with a cell type column will do;
#                              point `cellchat_group_by` at that column)
#
# Main steps
#   1. createCellChat on the cell type annotation.
#   2. Load the ligand-receptor database and the PPI network configured in
#      00_config.R.
#   3. identifyOverExpressedGenes / identifyOverExpressedInteractions.
#   4. projectData, computeCommunProb and computeCommunProbPathway.
#   5. filterCommunication and aggregateNet.
#   6. Write the communication tables and draw the aggregated network.
#
# Output
#   results/cellchat/06_cellchat.rds
#   results/cellchat/06_communication_LR.csv
#   results/cellchat/06_communication_pathway.csv
#   figures/06_*.pdf
#
# This is a standard CellChat workflow, not a reproduction of a specific figure.
# Only the aggregated network is drawn by default so that a full pipeline run
# stays fast; any further figure can be produced from the saved CellChat object.

source(here::here("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(ComplexHeatmap)
})

options(stringsAsFactors = FALSE)

# CellChat functions return either a ggplot or a ComplexHeatmap depending on the
# function and on the installed version; this writes both to a file.
save_plot <- function(x, file, width = 8, height = 6) {
  if (is.null(x)) {
    message("  nothing returned for ", basename(file), ", skipped")
    return(invisible(NULL))
  }
  if (inherits(x, c("Heatmap", "HeatmapList"))) {
    pdf(file, width = width, height = height)
    ComplexHeatmap::draw(x)
    dev.off()
  } else {
    ggplot2::ggsave(filename = file, plot = x, width = width, height = height)
  }
}

annotated <- readRDS(file.path(results_dir, "03_annotated.rds"))

if (!cellchat_group_by %in% colnames(annotated@meta.data)) {
  stop("Column '", cellchat_group_by, "' is not present in the object. Set ",
       "`cellchat_group_by` in scripts/00_config.R to the metadata column ",
       "holding the cell type annotation.")
}

Idents(annotated) <- cellchat_group_by

# ---------------------------------------------------------------------------
# 1.-2. Create the object and load the database
# ---------------------------------------------------------------------------
cellchat <- createCellChat(object = annotated, group.by = cellchat_group_by)

# The ligand-receptor database and the PPI network are selected in 00_config.R.
if (!exists(cellchat_db)) {
  stop("Object '", cellchat_db, "' is not available. Check `cellchat_db` in ",
       "scripts/00_config.R (use CellChatDB.mouse for mouse data).")
}
cellchat@DB <- get(cellchat_db)

levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))
message("06: ", length(levels(cellchat@idents)), " cell groups, ",
        ncol(annotated), " cells")

p_db <- showDatabaseCategory(get(cellchat_db))
save_plot(p_db, file.path(figures_dir, "06_cellchat_database_category.pdf"),
          width = 8, height = 6)

# ---------------------------------------------------------------------------
# 3.-5. Communication probabilities and aggregated network
# ---------------------------------------------------------------------------
cellchat <- subsetData(cellchat)

# The default limit for objects exported to the workers (500 MiB) is too small
# for a full expression matrix, so it is raised before the parallel step.
future::plan("multisession", workers = n_workers)
options(future.globals.maxSize = future_max_bytes)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Project the averaged expression onto the PPI network before computing the
# communication probabilities (CellChat tutorial workflow). Set
# `cellchat_ppi <- NULL` in 00_config.R to skip the projection and use the
# raw averaged expression instead.
if (is.null(cellchat_ppi)) {
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
} else {
  if (!exists(cellchat_ppi)) {
    stop("Object '", cellchat_ppi, "' is not available. Check `cellchat_ppi` in ",
         "scripts/00_config.R (use PPI.mouse for mouse data).")
  }
  cellchat <- projectData(cellchat, get(cellchat_ppi))
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE)
}

# Drop communication involving cell groups with too few cells
cellchat <- filterCommunication(cellchat, min.cells = cellchat_min_cells)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

df_net <- subsetCommunication(cellchat)
write.csv(df_net, file = file.path(results_cellchat_dir, "06_communication_LR.csv"))

df_netP <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df_netP, file = file.path(results_cellchat_dir, "06_communication_pathway.csv"))

# ---------------------------------------------------------------------------
# 6. Basic network figures
# ---------------------------------------------------------------------------
# Only the aggregated network is drawn by default, to keep `run_all.R` fast.
# The full communication tables are written to results/cellchat/ above, and the
# CellChat object is saved below, so any further figure can be produced on demand
# from the saved object.

# Number and strength of the interactions between all cell groups
pdf(file.path(figures_dir, "06_network_number_and_strength.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE,
                 label.edge = FALSE, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE,
                 label.edge = FALSE, title.name = "Interaction strength")
dev.off()
par(mfrow = c(1, 1))

# Heatmaps of the aggregated network
save_plot(netVisual_heatmap(cellchat),
          file.path(figures_dir, "06_network_heatmap_count.pdf"), width = 8, height = 6)
save_plot(netVisual_heatmap(cellchat, measure = "weight"),
          file.path(figures_dir, "06_network_heatmap_strength.pdf"), width = 8, height = 6)

saveRDS(cellchat, file = file.path(results_cellchat_dir, "06_cellchat.rds"))

message("06: ", length(cellchat@netP$pathways), " signalling pathways detected; ",
        "communication tables and the CellChat object written to ", results_cellchat_dir)
