#############################################
# Script:    01_import_and_qc.R
# Purpose:   Read the per-sample 10x Genomics matrices, create one Seurat object
#            per sample, compute QC metrics, apply the single set of QC
#            thresholds defined in 00_config.R, merge the samples and establish
#            the sample/group/batch metadata columns.
# Related manuscript figures: n/a (standard workflow, no figure-specific code)
# Required R packages: Seurat, ggplot2
#############################################
#
# Input
#   data/<sample_dir>/  10x Genomics output (barcodes.tsv.gz, features.tsv.gz,
#                       matrix.mtx.gz) for every sample listed in
#                       data/sample_sheet.csv
#
# Main steps
#   1. Read one 10x matrix per sample, create a Seurat object
#      (min.cells = qc_min_cells_per_gene, min.features = qc_min_features_per_cell)
#      and set `project` to the sample identifier, so that `orig.ident` holds a
#      well-defined sample label from the start.
#   2. Mitochondrial fraction per cell (pattern = mito_pattern).
#   3. QC violin plots before filtering.
#   4. Cell filtering with the single set of thresholds from 00_config.R.
#   5. QC violin plots after filtering and a table of retained cells per sample.
#   6. Merge the samples and create the metadata columns used downstream:
#      `sample_type` (group) and `batch` (sample), both derived from `orig.ident`.
#
# Output
#   results/01_merged_qc.rds                       merged, QC-filtered object
#   results/tables/01_cells_per_sample_after_QC.csv
#   figures/01_qc_violin_before_filtering.pdf
#   figures/01_qc_violin_after_filtering.pdf
#
# No other script in this repository filters cells or re-creates these metadata
# columns.

source(here::here("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# ---------------------------------------------------------------------------
# 1. One Seurat object per sample
# ---------------------------------------------------------------------------
scRNAlist <- vector("list", nrow(sample_sheet))
names(scRNAlist) <- sample_sheet$sample_id

for (i in seq_len(nrow(sample_sheet))) {
  counts <- Read10X(data.dir = sample_dirs[i])
  scRNAlist[[i]] <- CreateSeuratObject(
    counts,
    project      = sample_sheet$sample_id[i],
    min.cells    = qc_min_cells_per_gene,
    min.features = qc_min_features_per_cell
  )
}

# ---------------------------------------------------------------------------
# 2.-5. QC metrics, filtering and diagnostic plots
# ---------------------------------------------------------------------------
qc_before <- list()
qc_after  <- list()

for (i in seq_along(scRNAlist)) {
  obj <- scRNAlist[[i]]

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mito_pattern)

  qc_before[[names(scRNAlist)[i]]] <-
    VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

  # The QC thresholds come from 00_config.R only; the boolean vector is built
  # explicitly so that the values that are applied are auditable.
  keep <- obj$nFeature_RNA > qc_min_features_per_cell &
          obj$nFeature_RNA < qc_max_features_per_cell &
          obj$percent.mt   < qc_max_percent_mt

  if (!any(keep)) {
    stop("No cell of sample '", names(scRNAlist)[i], "' passed the QC thresholds ",
         "in scripts/00_config.R. Check the thresholds and the input data.")
  }

  obj <- subset(obj, cells = colnames(obj)[keep])
  scRNAlist[[i]] <- obj

  qc_after[[names(scRNAlist)[i]]] <-
    VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
}

pdf(file.path(figures_dir, "01_qc_violin_before_filtering.pdf"), width = 12, height = 5)
for (p in qc_before) print(p)
dev.off()

pdf(file.path(figures_dir, "01_qc_violin_after_filtering.pdf"), width = 12, height = 5)
for (p in qc_after) print(p)
dev.off()

write.csv(data.frame(sample = names(scRNAlist),
                     cells  = vapply(scRNAlist, ncol, integer(1))),
          file = file.path(results_tables_dir, "01_cells_per_sample_after_QC.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 6. Merge and create the metadata used downstream
# ---------------------------------------------------------------------------
# `y` holds all samples but the first; `add.cell.ids` therefore has one
# identifier per object, i.e. length(y) + 1 entries. With a single sample there
# is nothing to merge.
if (length(scRNAlist) > 1) {
  combined <- merge(scRNAlist[[1]],
                    y = scRNAlist[-1],
                    add.cell.ids = names(scRNAlist),
                    project = "combined")
} else {
  combined <- scRNAlist[[1]]
}

# Metadata contract for the whole pipeline, created once and only here:
#   orig.ident   sample identifier (set by `project` above)
#   sample_type  study group, derived from orig.ident via the sample sheet
#   batch        the variable Harmony corrects for (one batch per sample)
combined$sample_type <- assign_sample_type(combined$orig.ident)
combined$batch       <- combined$orig.ident

if (anyNA(combined$sample_type)) {
  stop("Some cells could not be assigned to a group. Check that every entry of ",
       "`data_dir` in data/sample_sheet.csv exists and belongs to a known group.")
}

saveRDS(combined, file = file.path(results_dir, "01_merged_qc.rds"))

message("01: ", ncol(combined), " cells retained across ",
        length(unique(combined$orig.ident)), " samples; per-sample counts in ",
        "results/tables/01_cells_per_sample_after_QC.csv")
