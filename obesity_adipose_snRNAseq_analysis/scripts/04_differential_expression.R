#############################################
# Script:    04_differential_expression.R
# Purpose:   Differential expression between the two study groups, computed
#            separately within every annotated cell type.
# Related manuscript figures: n/a (standard workflow, no figure-specific code)
# Required R packages: Seurat, dplyr
#############################################
#
# Input
#   results/03_annotated.rds   annotated object written by 03
#
# Main steps
#   1. For every annotated cell type, test group_A versus group_B with the
#      parameters from 00_config.R (Wilcoxon by default).
#   2. Write one table per cell type and one combined table.
#   3. Write a summary table of the number of differentially expressed genes per
#      cell type; this summary also drives which cell type the enrichment
#      analysis in 07 uses (see `enrichment_celltype` in 00_config.R).
#
# Output
#   results/deg/04_DE_<celltype>_<groupA>_vs_<groupB>.csv   one table per cell type
#   results/deg/04_DE_all_celltypes.csv
#   results/tables/04_DE_summary.csv
#
# Cells labelled "unassigned" by 03 are excluded.

source(here::here("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

annotated <- readRDS(file.path(results_dir, "03_annotated.rds"))
Idents(annotated) <- "celltype"

celltypes <- setdiff(sort(unique(as.character(annotated$celltype))), "unassigned")
message("04: comparing ", group_levels[1], " vs ", group_levels[2],
        " within ", length(celltypes), " cell types")

de_all    <- list()
de_summary <- data.frame()

for (ct in celltypes) {
  cells <- colnames(annotated)[annotated$celltype == ct]
  sub   <- subset(annotated, cells = cells)

  Idents(sub) <- de_group_column
  n_per_group <- table(Idents(sub))

  # A group comparison needs enough cells in both groups
  if (length(n_per_group) < 2 || any(n_per_group < de_min_cells_per_group)) {
    message("  ", ct, ": skipped (",
            paste(names(n_per_group), n_per_group, sep = "=", collapse = ", "), ")")
    next
  }

  de <- FindMarkers(
    object          = sub,
    ident.1         = group_levels[1],
    ident.2         = group_levels[2],
    test.use        = de_test,
    logfc.threshold = de_logfc_threshold,
    min.pct         = de_min_pct
  )

  if (nrow(de) == 0) {
    message("  ", ct, ": no gene passed the logFC threshold")
    next
  }

  de$gene     <- rownames(de)
  de$celltype <- ct

  fname <- file.path(
    results_deg_dir,
    sprintf("04_DE_%s_%s_vs_%s.csv", sanitise_filename(ct),
            sanitise_filename(group_levels[1]), sanitise_filename(group_levels[2]))
  )
  write.csv(de, file = fname, row.names = FALSE)

  significant <- de$p_val_adj < de_padj_cutoff & abs(de$avg_log2FC) > de_logfc_cutoff

  de_summary <- rbind(de_summary, data.frame(
    celltype        = ct,
    n_cells         = length(cells),
    n_tested        = nrow(de),
    n_significant   = sum(significant),
    n_up            = sum(significant & de$avg_log2FC > 0),
    n_down          = sum(significant & de$avg_log2FC < 0)
  ))

  de_all[[ct]] <- de
  message("  ", ct, ": ", nrow(de), " genes tested, ", sum(significant), " significant")
}

if (!length(de_all)) {
  stop("No cell type produced a differential expression result. Check the ",
       "annotation in 03 and the group assignments in data/sample_sheet.csv.")
}

de_combined <- bind_rows(de_all)
write.csv(de_combined, file = file.path(results_deg_dir, "04_DE_all_celltypes.csv"),
          row.names = FALSE)

de_summary <- de_summary[order(-de_summary$n_significant), ]
write.csv(de_summary, file = file.path(results_tables_dir, "04_DE_summary.csv"),
          row.names = FALSE)

saveRDS(de_summary, file = file.path(results_dir, "04_DE_summary.rds"))

message("04: results written for ", length(de_all), " cell types; summary in ",
        "results/tables/04_DE_summary.csv")
