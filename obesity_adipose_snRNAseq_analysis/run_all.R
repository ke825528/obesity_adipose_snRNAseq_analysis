#############################################
# Script:    run_all.R
# Purpose:   Run the complete pipeline in order. This is the entry point for a
#            reproducible run of the whole repository.
# Related manuscript figures: n/a (standard workflow, no figure-specific code)
# Required R packages: here
#############################################
#
# Usage, from the repository root:
#
#   Rscript run_all.R
#
# Input
#   data/sample_sheet.csv  and the 10x directories it points to
#   data/c2.cp.kegg.v7.5.1.symbols.gmt  (only needed by 07)
#
# The pipeline is a sequence of scripts that pass objects to each other through
# results/*.rds:
#
#   01_import_and_qc.R                   -> results/01_merged_qc.rds
#   02_integration_and_clustering.R      -> results/02_clustered.rds
#   03_marker_analysis_and_annotation.R  -> results/03_annotated.rds
#   04_differential_expression.R         -> results/deg/*.csv
#   05_subclustering_analysis.R          -> results/05_subclustered.rds
#   06_cellchat_analysis.R               -> results/cellchat/*.rds
#   07_pathway_analysis.R                -> results/enrichment/*.csv
#   08_visualization.R                   -> figures/*.pdf
#
# Each script is also runnable on its own once the objects it needs exist.

source(here::here("scripts", "00_config.R"))

step_scripts <- c(
  "01_import_and_qc.R",
  "02_integration_and_clustering.R",
  "03_marker_analysis_and_annotation.R",
  "04_differential_expression.R",
  "05_subclustering_analysis.R",
  "06_cellchat_analysis.R",
  "07_pathway_analysis.R",
  "08_visualization.R"
)

# Fail early with a clear message rather than part-way through the pipeline.
missing_samples <- sample_sheet$data_dir[
  !dir.exists(file.path(data_dir, sample_sheet$data_dir))
]
if (length(missing_samples)) {
  stop("These sample directories are missing from data/: ",
       paste(missing_samples, collapse = ", "), "\n",
       "See README.md ('Input data') for how to provide the input data.")
}

# The GMT file is only needed by 07, which skips its GSEA section if it is
# absent, so this is a warning rather than an error.
if (!file.exists(gsea_gmt_file)) {
  warning("KEGG GMT file not found: ", gsea_gmt_file, "\n",
          "  The GSEA section of 07_pathway_analysis.R will be skipped. See the ",
          "download note for `gsea_gmt_file` in scripts/00_config.R.")
}

for (s in step_scripts) {
  path <- file.path(project_dir, "scripts", s)
  if (!file.exists(path)) stop("Missing pipeline script: ", path)

  message("\n", strrep("=", 70))
  message("Running ", s)
  message(strrep("=", 70))

  source(path)
}

message("\nPipeline finished. Figures: ", figures_dir, "  Results: ", results_dir)
