# Obesity adipose snRNA-seq analysis

## Overview

This repository provides a standardized single-nucleus RNA-seq analysis workflow
for the adipose tissue dataset associated with the manuscript. The workflow
includes quality control, clustering, cell-type annotation, differential
expression, functional enrichment, CellChat analysis, and descriptive
visualization.

This repository provides a reproducible analysis workflow and is not intended to
reproduce every manuscript-specific exploratory analysis.

## Workflow

1. Data import and quality control
2. Normalization and Harmony integration
3. Clustering and cell-type annotation
4. Differential expression
5. Subclustering
6. CellChat analysis
7. GO / KEGG / GSEA
8. Visualization

## Repository structure

```text
obesity_adipose_snRNAseq_analysis/
├── README.md
├── R_packages.txt
├── run_all.R
├── data/
├── scripts/
├── results/
└── figures/
```

`scripts/` contains:

```text
00_config.R                         paths and all analysis parameters
01_import_and_qc.R                  import, QC, sample metadata
02_integration_and_clustering.R     normalization, Harmony, UMAP, clustering
03_marker_analysis_and_annotation.R cluster markers, cell-type annotation
04_differential_expression.R        group comparison within each cell type
05_subclustering_analysis.R         subclustering of one compartment
06_cellchat_analysis.R              standard CellChat workflow
07_pathway_analysis.R               GO / KEGG / GSEA
08_visualization.R                  composition, UMAP, gene expression panel
```

## Input data

* 10x Genomics matrices, one directory per sample.
* Sample information in `data/sample_sheet.csv` (sample identifier, directory,
  group).
* Raw patient-level data are not included.

Data accessions: **OMIX016079** (snRNA-seq), **PRJNA763819** (bulk RNA-seq).

## Run

```bash
Rscript run_all.R
```

The scripts can also be run individually in numerical order.

## Software

R >= 4.2
Seurat >= 5.0

Main dependencies:

* Seurat
* harmony
* CellChat
* clusterProfiler
* fgsea
* msigdbr
* ggplot2

The complete list is in `R_packages.txt`.
