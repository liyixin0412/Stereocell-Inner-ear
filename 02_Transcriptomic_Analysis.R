# Key transcriptomic analyses. Select one independent task.
# Tasks: p1, adult, subtype, hearing
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
args <- commandArgs(trailingOnly=TRUE)
task <- if (length(args)) args[1] else "subtype"
cell_type <- if (length(args) >= 2) args[2] else "Epi"
out_dir <- "results/02"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

if (task == "p1") {
library(Seurat)
library(dplyr)
library(Matrix)
obj <- readRDS("data/P1_annotated.rds")
DefaultAssay(obj) <- "RNA"
stopifnot(all(c("celltype", "source") %in% colnames(obj@meta.data)))
qc_summary <- obj@meta.data %>%
  group_by(celltype, source) %>%
  summarise(n_cells = n(), mean_genes = mean(nFeature_RNA),
            mean_UMI = mean(nCount_RNA), .groups = "drop")
composition <- obj@meta.data %>% count(celltype, source) %>%
  group_by(celltype) %>% mutate(proportion = n / sum(n)) %>% ungroup()
write.csv(qc_summary, file.path(out_dir, "P1_QC_summary.csv"), row.names = FALSE)
write.csv(composition, file.path(out_dir, "P1_source_composition.csv"), row.names = FALSE)


counts <- GetAssayData(
  obj,
  assay = "RNA",
  slot = "counts"
)

meta <- obj@meta.data

celltypes <- unique(meta$celltype)

detected_list <- lapply(celltypes, function(ct) {

  cells_sc <- rownames(meta)[
    meta$celltype == ct &
    meta$source == "Single-Cell"
  ]

  cells_sn <- rownames(meta)[
    meta$celltype == ct &
    meta$source == "Single-Nuclear"
  ]

  if(length(cells_sc) == 0 | length(cells_sn) == 0) {
    return(NULL)
  }

  n_gene_sc <- sum(
    Matrix::rowSums(counts[, cells_sc, drop = FALSE] > 0) > 0
  )

  n_gene_sn <- sum(
    Matrix::rowSums(counts[, cells_sn, drop = FALSE] > 0) > 0
  )

  data.frame(
    celltype = ct,
    Single_Cell = n_gene_sc,
    Single_Nuclear = n_gene_sn
  )
})

df_detected <- bind_rows(detected_list)

df_detected
write.csv(df_detected, file.path(out_dir, "P1_detected_genes.csv"), row.names=FALSE)



} else if (task == "adult") {
obj <- readRDS('data/Merge-Cochlea-Final-V3.rds')
rna <- GetAssayData(obj,assay="RNA",slot="counts")

sn_cells <- colnames(obj)[which(obj$source=="snRNA")]
st_cells <- colnames(obj)[which(obj$source=="StereoCell")]

prop_sn <- Matrix::rowMeans(
  rna[,sn_cells,drop=FALSE] > 0,
  na.rm=TRUE
)

prop_st <- Matrix::rowMeans(
  rna[,st_cells,drop=FALSE] > 0,
  na.rm=TRUE
)

genes_sn <- rownames(rna)[
  !is.na(prop_sn) & prop_sn >= 0.01
]

genes_st <- rownames(rna)[
  !is.na(prop_st) & prop_st >= 0.01
]

genes_sn <- unique(genes_sn[!is.na(genes_sn) & genes_sn!=""])
genes_st <- unique(genes_st[!is.na(genes_st) & genes_st!=""])


shared_genes <- intersect(genes_sn, genes_st)
avg_expr <- AverageExpression(obj, assays="RNA", slot="data", group.by="source")$RNA
comparison <- data.frame(x=log1p(avg_expr[, "snRNA"]), y=log1p(avg_expr[, "StereoCell"]))
pearson_r <- cor(comparison$x, comparison$y, method="pearson")
fit <- lm(y ~ x, data=comparison)
r_squared <- summary(fit)$r.squared

write.csv(data.frame(shared=length(shared_genes),
  only_10x=length(setdiff(genes_sn, genes_st)),
  only_StereoCell=length(setdiff(genes_st, genes_sn)),
  Pearson_r=pearson_r, R_squared=r_squared), file.path(out_dir, "Adult_method_comparison.csv"), row.names=FALSE)
write.csv(comparison, file.path(out_dir, "Adult_expression_comparison.csv"), row.names=TRUE)


} else if (task == "subtype") {
  paths <- c(SGN="data/Merge-SGN-Final.rds",
             Glia="data/Merge-Glia-Final.rds",
             Epi="data/Merge-Epi-Final-V2.rds")
  stopifnot(cell_type %in% names(paths))
  obj <- readRDS(paths[[cell_type]])
  group_col <- if (cell_type == "Glia") "Glia_subtype" else "subtype_Final"
  stopifnot(group_col %in% colnames(obj@meta.data))
  DefaultAssay(obj) <- "RNA"
  Idents(obj) <- group_col
  # Thresholds retained from the supplied notebooks.
  if (cell_type == "SGN") {
    markers <- FindAllMarkers(obj, only.pos=TRUE) %>%
      filter(avg_log2FC > 0.25)
  } else {
    markers <- FindAllMarkers(obj, only.pos=TRUE,
                              logfc.threshold=0.5, min.pct=0.2)
  }
  write.csv(markers, file.path(out_dir, paste0(cell_type, "_markers.csv")), row.names=FALSE)
  p <- DimPlot(obj, group.by=group_col, label=TRUE)
  ggsave(file.path(out_dir, paste0(cell_type, "_UMAP.pdf")), p, width=7, height=5)
  avg <- AverageExpression(obj, assays="RNA", slot="data", group.by=group_col)$RNA
  correlation <- cor(log1p(avg), method="pearson")
  write.csv(correlation, file.path(out_dir, paste0(cell_type, "_Pearson.csv")))
  # Unfiltered signed ranking for GO-BP GSEA.
  library(msigdbr)
  library(fgsea)
table(obj[[group_col, drop=TRUE]])
Idents(obj) <- group_col

rank_table <- FindAllMarkers(
  object = obj,
  assay = "RNA",
  slot = "data",
  test.use = "wilcox",
  only.pos = FALSE,
  logfc.threshold = 0,
  min.pct = 0.01,
  return.thresh = Inf
)


m_df <- msigdbr(
  species = "Mus musculus",
  category = "C5"
) %>%
  filter(grepl("^GOBP_", gs_name)) %>%
  distinct(gs_name, gene_symbol)

stopifnot(nrow(m_df) > 0)

fgsea_sets <- split(m_df$gene_symbol, m_df$gs_name)


set.seed(123)
gsea_list <- list()
rank_list <- list()

for (cl in unique(as.character(rank_table$cluster))) {

  deg <- rank_table %>%
    filter(cluster == cl, is.finite(avg_log2FC)) %>%
    distinct(gene, .keep_all = TRUE)

  ranks <- setNames(deg$avg_log2FC, deg$gene)
  ranks <- sort(ranks, decreasing = TRUE)
  rank_list[[cl]] <- ranks

  res <- fgseaMultilevel(
    pathways = fgsea_sets,
    stats = ranks,
    minSize = 15,
    maxSize = 500,
    scoreType = "std"
  ) %>%
    as.data.frame()

  if (nrow(res) == 0) next

  res$cluster <- cl
  res$leadingEdge <- vapply(
    res$leadingEdge, paste, character(1), collapse = ","
  )
  gsea_list[[cl]] <- res
}

gsea_all <- bind_rows(gsea_list) %>%
  arrange(cluster, padj)

write.csv(
  gsea_all,
  file.path(out_dir, paste0(cell_type, "_GSEA_GOBP.csv")),
  row.names = FALSE
)

} else if (task == "hearing") {
  # Use independently curated lists and independently computed stage matrices.
  stage <- if (length(args) >= 2) args[2] else "P1"
  stopifnot(stage %in% c("P1", "Adult"))
  path <- if (stage == "P1") "data/P1_SN_annotated.rds" else "data/Adult_StereoCell_annotated.rds"
  obj <- readRDS(path)
  gene_df <- read.csv("data/Hearing_loss_genes_HL-adult-P1.csv", stringsAsFactors=FALSE)
  disease_col <- paste0(stage, "_Disease")
  genes <- gene_df[[stage]]
  disease <- gene_df[[disease_col]]
  use <- !is.na(genes) & nzchar(genes) & !is.na(disease) & nzchar(disease)
  gene_anno <- unique(data.frame(Gene_human=genes[use], Disease=disease[use]))
  # Retains the symbol conversion used in the figure code; not an orthology lookup.
  gene_anno$Gene <- paste0(toupper(substr(gene_anno$Gene_human,1,1)),
                          tolower(substring(gene_anno$Gene_human,2)))
  gene_anno <- gene_anno[gene_anno$Gene %in% rownames(obj[["RNA"]]), ]
  genes_use <- unique(gene_anno$Gene)
  stopifnot(length(genes_use) > 0)
  avg <- AverageExpression(obj, assays="RNA", features=genes_use,
                           group.by="celltype", slot="data")$RNA
  z <- t(scale(t(as.matrix(avg[genes_use, , drop=FALSE]))))
  z[is.na(z)] <- 0
  write.csv(z, file.path(out_dir, paste0(stage, "_hearing_genes_Zscore.csv")))
  write.csv(gene_anno, file.path(out_dir, paste0(stage, "_hearing_gene_annotation.csv")), row.names=FALSE)
  pheatmap::pheatmap(z, cluster_rows=FALSE, cluster_cols=FALSE,
    color=colorRampPalette(c("#00A6E8","white","#FF5A36"))(100),
    breaks=seq(-2,2,length.out=101),
    filename=file.path(out_dir, paste0(stage, "_hearing_genes.pdf")))
} else stop("task must be p1, adult, subtype or hearing")
