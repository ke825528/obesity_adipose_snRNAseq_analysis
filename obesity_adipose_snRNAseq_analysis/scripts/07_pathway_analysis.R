#############################################
# Script:    07_pathway_analysis.R
# Purpose:   Standard functional enrichment of a differential expression result:
#            GO over-representation, KEGG over-representation, GSEA against the
#            MSigDB KEGG gene sets and fgsea against the MSigDB C2:KEGG
#            collection.
# Related manuscript figures: n/a (standard workflow, no figure-specific code)
# Required R packages: clusterProfiler, enrichplot, fgsea, msigdbr, ggplot2
#############################################
#
# Input
#   results/deg/04_DE_all_celltypes.csv   differential expression written by 04
#   results/tables/04_DE_summary.csv      per-cell-type summary written by 04
#   data/<c2.cp.kegg.v7.5.1.symbols.gmt>  MSigDB KEGG gene sets (see 00_config.R)
#
# Main steps
#   1. Pick the differential expression result to analyse: `enrichment_celltype`
#      from 00_config.R, or the cell type with the most significant genes if it
#      is NULL.
#   2. GO over-representation for every ontology in `enrichment_go_ontologies`.
#   3. KEGG over-representation.
#   4. GSEA of the ranked log2 fold changes against the KEGG GMT file.
#   5. fgsea of the same ranking against MSigDB C2:KEGG.
#
# The organism used by all four analyses is derived from `organism` in
# 00_config.R (annotation database, KEGG organism code and MSigDB species are set
# together), so a human/mouse mismatch cannot occur.
#
# Output
#   results/enrichment/07_GO_<ontology>.csv
#   results/enrichment/07_KEGG.csv
#   results/enrichment/07_GSEA.csv
#   results/enrichment/07_fgsea.csv
#   results/enrichment/gsea/*.pdf
#   figures/07_*.pdf

source(here::here("scripts", "00_config.R"))

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(enrichment_orgdb, character.only = TRUE)
})

gsea_fig_dir <- file.path(results_enrich_dir, "gsea")
if (!dir.exists(gsea_fig_dir)) dir.create(gsea_fig_dir, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Select the differential expression result to analyse
# ---------------------------------------------------------------------------
de_all     <- read.csv(file.path(results_deg_dir, "04_DE_all_celltypes.csv"),
                       stringsAsFactors = FALSE)
de_summary <- read.csv(file.path(results_tables_dir, "04_DE_summary.csv"),
                       stringsAsFactors = FALSE)

if (is.null(enrichment_celltype)) {
  enrichment_celltype <- de_summary$celltype[which.max(de_summary$n_significant)]
  message("07: `enrichment_celltype` is NULL, using the cell type with the most ",
          "significant genes: ", enrichment_celltype)
}

if (!enrichment_celltype %in% de_summary$celltype) {
  stop("Cell type '", enrichment_celltype, "' has no differential expression ",
       "result. Available: ", paste(de_summary$celltype, collapse = ", "))
}

de <- de_all[de_all$celltype == enrichment_celltype, ]
if (!nrow(de)) stop("No differential expression rows for '", enrichment_celltype, "'")

message("07: enrichment of ", enrichment_celltype, " (", nrow(de), " genes tested, ",
        "organism = ", organism, ")")

# Ranked list for GSEA / fgsea: log2 fold change named by gene symbol
gene_list <- de$avg_log2FC
names(gene_list) <- de$gene
gene_list <- sort(gene_list[!is.na(gene_list)], decreasing = TRUE)

# Gene set for over-representation
is_significant <- !is.na(de$p_val_adj) & !is.na(de$avg_log2FC) &
                  de$p_val_adj < enrichment_gene_padj_cutoff &
                  abs(de$avg_log2FC) > enrichment_gene_logfc_cutoff
genes_significant <- de$gene[is_significant]
message("07: ", length(genes_significant), " significant genes used for GO / KEGG")

if (length(genes_significant) < enrichment_min_genes) {
  warning("Fewer than ", enrichment_min_genes, " significant genes; ",
          "over-representation results will be unstable. Consider relaxing ",
          "`enrichment_gene_padj_cutoff` or `enrichment_gene_logfc_cutoff` in ",
          "00_config.R.")
}

if (length(genes_significant) < 1) {
  stop("No significant gene for '", enrichment_celltype, "'. Nothing to enrich. ",
       "Relax `enrichment_gene_padj_cutoff` / `enrichment_gene_logfc_cutoff` in ",
       "00_config.R, or choose another cell type with `enrichment_celltype`.")
}

# ---------------------------------------------------------------------------
# 2. GO over-representation
# ---------------------------------------------------------------------------
# enrichGO() returns NULL (not an empty table) when no term passes the cut-offs,
# so the NULL case is handled before the result is converted.
for (ont in enrichment_go_ontologies) {
  ego <- enrichGO(gene          = genes_significant,
                  OrgDb         = enrichment_orgdb,
                  keyType       = "SYMBOL",
                  ont           = ont,
                  pAdjustMethod = enrichment_padjust_method,
                  pvalueCutoff  = enrichment_pvalue_cutoff,
                  qvalueCutoff  = enrichment_qvalue_cutoff)

  if (is.null(ego) || !nrow(as.data.frame(ego))) {
    message("  GO ", ont, ": no term passed the cut-offs")
    next
  }

  ego_df <- as.data.frame(ego)
  write.csv(ego_df, file = file.path(results_enrich_dir, paste0("07_GO_", ont, ".csv")),
            row.names = FALSE)

  p <- barplot(ego, showCategory = enrichment_show_category) +
    ggtitle(paste("GO", ont, "-", enrichment_celltype))
  ggsave(file.path(figures_dir, paste0("07_GO_", ont, "_barplot.pdf")), p,
         width = 9, height = 6)
}

# ---------------------------------------------------------------------------
# 3. KEGG over-representation
# ---------------------------------------------------------------------------
# KEGG queries the KEGG REST API; a failure (e.g. no internet access) is reported
# instead of stopping the pipeline.
entrez <- tryCatch(
  bitr(genes_significant, fromType = "SYMBOL", toType = "ENTREZID",
       OrgDb = enrichment_orgdb),
  error = function(e) {
    message("  gene symbol -> ENTREZID conversion failed: ", conditionMessage(e))
    NULL
  }
)

ekegg <- NULL
if (!is.null(entrez) && nrow(entrez)) {
  ekegg <- tryCatch(
    enrichKEGG(gene          = entrez$ENTREZID,
               organism      = enrichment_kegg,
               pAdjustMethod = enrichment_padjust_method,
               pvalueCutoff  = enrichment_pvalue_cutoff,
               qvalueCutoff  = enrichment_qvalue_cutoff),
    error = function(e) {
      message("  KEGG over-representation failed (network access?): ",
              conditionMessage(e))
      NULL
    }
  )
}

if (!is.null(ekegg) && nrow(as.data.frame(ekegg))) {
  write.csv(as.data.frame(ekegg), file = file.path(results_enrich_dir, "07_KEGG.csv"),
            row.names = FALSE)
  p_kegg <- dotplot(ekegg, showCategory = enrichment_show_category) +
    ggtitle(paste("KEGG -", enrichment_celltype))
  ggsave(file.path(figures_dir, "07_KEGG_dotplot.pdf"), p_kegg, width = 9, height = 7)
}

# ---------------------------------------------------------------------------
# 4. GSEA against the KEGG GMT file
# ---------------------------------------------------------------------------
# The GMT file is an external download. Its absence is reported and the GSEA
# part is skipped rather than stopping the pipeline, so that sections 1-3 and 08
# still run.
geneset <- NULL

if (!file.exists(gsea_gmt_file)) {
  warning("KEGG GMT file not found: ", gsea_gmt_file, "\n",
          "  GSEA (section 4) is skipped. See the download note for ",
          "`gsea_gmt_file` in scripts/00_config.R.")
} else {
  geneset <- read.gmt(gsea_gmt_file)
  n_overlap <- length(intersect(names(gene_list), unique(geneset$gene)))
  message("07: ", n_overlap, " of ", length(gene_list),
          " genes in the ranked list are present in the GMT gene sets")
  if (n_overlap < gsea_min_overlap) {
    warning("Very little overlap between the ranked list and the GMT file. ",
            "Check that the gene sets match the species of the data ",
            "(organism is '", organism, "').")
  }

  egmt <- GSEA(gene_list, TERM2GENE = geneset,
               minGSSize = gsea_min_gs_size, pvalueCutoff = gsea_pvalue_cutoff,
               verbose = FALSE)

  if (is.null(egmt) || !nrow(as.data.frame(egmt))) {
    message("  GSEA: no gene set passed the cut-offs")
  } else {
    gsea_df <- as.data.frame(egmt)
    write.csv(gsea_df, file = file.path(results_enrich_dir, "07_GSEA.csv"),
              row.names = FALSE)

    gsea_show <- head(gsea_df[order(gsea_df$p.adjust), ], gsea_top_pathways)
    for (i in seq_len(nrow(gsea_show))) {
      pdf(file.path(gsea_fig_dir, paste0("07_GSEA_", sanitise_filename(gsea_show$ID[i]), ".pdf")),
          width = 11, height = 5)
      print(gseaplot2(egmt, geneSetID = gsea_show$ID[i],
                      title = gsea_show$Description[i], pvalue_table = TRUE))
      dev.off()
    }
  }
}

# ---------------------------------------------------------------------------
# 5. fgsea against MSigDB C2:KEGG
# ---------------------------------------------------------------------------
# msigdbr >= 10 renamed `category`/`subcategory` to `collection`/`subcollection`;
# if the call fails for any reason the gene sets are built from the GMT file
# instead, so the pipeline does not stop here.
fetch_msigdb_kegg <- function(species) {
  if (utils::packageVersion("msigdbr") >= "10.0.0") {
    m <- tryCatch(
      msigdbr(species = species, collection = "C2", subcollection = "CP:KEGG_LEGACY"),
      error = function(e) msigdbr(species = species, collection = "C2",
                                  subcollection = "CP:KEGG")
    )
  } else {
    m <- msigdbr(species = species, category = "C2", subcategory = "KEGG")
  }
  split(x = m$gene_symbol, f = m$gs_name)
}

fgsea_sets <- tryCatch(
  fetch_msigdb_kegg(msigdb_species),
  error = function(e) {
    if (is.null(geneset)) stop(e)
    message("  msigdbr failed (", conditionMessage(e),
            "); falling back to the GMT gene sets")
    split(x = geneset$gene, f = geneset$term)
  }
)

fgsea_res <- fgsea(pathways = fgsea_sets, stats = gene_list)
fgsea_res <- fgsea_res[order(fgsea_res$padj), ]

# `leadingEdge` is a list column; collapse it so the CSV stays readable
fgsea_out <- as.data.frame(fgsea_res)
fgsea_out$leadingEdge <- vapply(fgsea_out$leadingEdge, paste,
                                character(1), collapse = ",")
write.csv(fgsea_out, file = file.path(results_enrich_dir, "07_fgsea.csv"),
          row.names = FALSE)

fgsea_sig <- fgsea_out[!is.na(fgsea_out$padj) & fgsea_out$padj < fgsea_padj_cutoff, ]

if (nrow(fgsea_sig)) {
  fgsea_top <- head(fgsea_sig[order(fgsea_sig$NES), ], fgsea_top_pathways)
  p_nes <- ggplot(fgsea_top, aes(x = reorder(pathway, NES), y = NES, fill = NES > 0)) +
    geom_col() +
    coord_flip() +
    labs(x = "Pathway", y = "Normalised enrichment score",
         title = paste("fgsea (MSigDB C2:KEGG) -", enrichment_celltype),
         fill = "Activated") +
    theme_minimal()
  ggsave(file.path(figures_dir, "07_fgsea_NES.pdf"), p_nes, width = 9, height = 7)

  top_pathway <- fgsea_sig$pathway[which.max(abs(fgsea_sig$NES))]
  p_enrich <- plotEnrichment(fgsea_sets[[top_pathway]], gene_list) +
    labs(title = top_pathway)
  ggsave(file.path(gsea_fig_dir, "07_fgsea_top_pathway.pdf"), p_enrich, width = 8, height = 5)
}

message("07: done (organism = ", organism, ", MSigDB species = ",
        msigdb_species, ")")
