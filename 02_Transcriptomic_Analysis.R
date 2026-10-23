# Transcriptomic comparisons, subtype markers and GO-BP GSEA.
# Computational results are exported in full; plotting sections are examples.
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)

save_plot <- function(plot, name, out_dir, width = 7, height = 5) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), plot,
         width = width, height = height, bg = "white")
}

# 1. Compare expression profiles and cell-type composition ---------------------
compare_methods <- function(obj, celltype_col, sources, out_dir,
                            cumulative_genes = FALSE) {
  DefaultAssay(obj) <- "RNA"
  stopifnot(length(sources) == 2, celltype_col %in% names(obj@meta.data),
            all(sources %in% as.character(obj$source)),
            all(c("nFeature_RNA", "nCount_RNA", "percent.mt") %in% names(obj@meta.data)),
            !anyNA(obj@meta.data[[celltype_col]]))
  obj <- subset(obj, cells = colnames(obj)[obj$source %in% sources])
  obj$source <- factor(obj$source, levels = sources)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  summary <- obj@meta.data %>%
    group_by(celltype = .data[[celltype_col]], source) %>%
    summarise(n_cells = n(), mean_genes = mean(nFeature_RNA),
              mean_UMI = mean(nCount_RNA), .groups = "drop")
  composition <- summary %>% group_by(celltype) %>%
    mutate(proportion = n_cells / sum(n_cells)) %>% ungroup()
  write.csv(summary, file.path(out_dir, "QC_summary.csv"), row.names = FALSE)
  write.csv(composition, file.path(out_dir, "source_composition.csv"), row.names = FALSE)

  # Example: QC distributions with mean values.
  plots <- lapply(c("nFeature_RNA", "nCount_RNA", "percent.mt"), function(feature) {
    means <- obj@meta.data %>% group_by(source) %>%
      summarise(value = mean(.data[[feature]], na.rm = TRUE), .groups = "drop")
    VlnPlot(obj, features = feature, group.by = "source", pt.size = 0) +
      geom_point(data = means, aes(source, value), inherit.aes = FALSE,
                 shape = 23, fill = "white", size = 3) +
      geom_text(data = means, aes(source, value, label = round(value, 1)),
                inherit.aes = FALSE, vjust = -1, size = 3) +
      theme_classic() + theme(legend.position = "none") + labs(x = NULL)
  })
  save_plot(patchwork::wrap_plots(plots, ncol = 3), "QC_violin", out_dir, 12, 4)
  p <- ggplot(composition, aes(celltype, proportion * 100, fill = source)) +
    geom_col(width = 0.75) + labs(x = NULL, y = "Proportion (%)", fill = NULL) +
    theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p, "source_composition", out_dir, 10, 5)

  counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
  meta <- obj@meta.data[colnames(counts), , drop = FALSE]
  if (cumulative_genes) {
    # P1: unique genes detected in at least one cell/nucleus of each cell type.
    detected <- bind_rows(lapply(unique(as.character(meta[[celltype_col]])), function(ct) {
      a <- rownames(meta)[meta[[celltype_col]] == ct & meta$source == sources[1]]
      b <- rownames(meta)[meta[[celltype_col]] == ct & meta$source == sources[2]]
      if (!length(a) || !length(b)) return(NULL)
      data.frame(celltype = ct, n_A = length(a), n_B = length(b),
        genes_A = sum(Matrix::rowSums(counts[, a, drop = FALSE] > 0) > 0),
        genes_B = sum(Matrix::rowSums(counts[, b, drop = FALSE] > 0) > 0))
    }))
    stopifnot(nrow(detected) > 0)
    write.csv(detected, file.path(out_dir, "detected_genes_by_celltype.csv"), row.names = FALSE)
    limit <- max(detected$genes_A, detected$genes_B) * 1.08
    p <- ggplot(detected, aes(genes_A, genes_B)) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(size = 3, color = "#2E86DE") +
      ggrepel::geom_text_repel(aes(label = celltype), seed = 123, max.overlaps = Inf) +
      coord_equal(xlim = c(0, limit), ylim = c(0, limit)) +
      labs(x = paste("Detected genes:", sources[1]), y = paste("Detected genes:", sources[2])) +
      theme_classic()
    save_plot(p, "detected_genes_by_celltype", out_dir, 6, 6)
  } else {
    # Adult: genes detected in at least 1% of nuclei in each source.
    gene_sets <- setNames(lapply(sources, function(s) {
      cells <- rownames(meta)[meta$source == s]
      proportion <- Matrix::rowMeans(counts[, cells, drop = FALSE] > 0)
      rownames(counts)[is.finite(proportion) & proportion >= 0.01]
    }), sources)
    write.csv(bind_rows(lapply(sources, function(s)
      data.frame(source = rep(s, length(gene_sets[[s]])), gene = gene_sets[[s]]))),
      file.path(out_dir, "detected_genes_1pct.csv"), row.names = FALSE)
    write.csv(data.frame(shared = length(intersect(gene_sets[[1]], gene_sets[[2]])),
      only_A = length(setdiff(gene_sets[[1]], gene_sets[[2]])),
      only_B = length(setdiff(gene_sets[[2]], gene_sets[[1]]))),
      file.path(out_dir, "gene_overlap.csv"), row.names = FALSE)
    venn <- VennDiagram::venn.diagram(gene_sets, filename = NULL,
      fill = c("#95e1d3", "#fce38a"), alpha = 0.6)
    pdf(file.path(out_dir, "gene_overlap_Venn.pdf"), width = 5, height = 5)
    grid::grid.draw(venn)
    dev.off()
  }
  avg <- AverageExpression(obj, assays = "RNA", slot = "data", group.by = "source")$RNA
  df <- data.frame(gene = rownames(avg), x = log1p(avg[, sources[1]]),
                   y = log1p(avg[, sources[2]]))
  fit <- lm(y ~ x, data = df)
  r <- cor(df$x, df$y, method = "pearson")
  r2 <- summary(fit)$r.squared
  write.csv(df, file.path(out_dir, "mean_expression.csv"), row.names = FALSE)
  write.csv(data.frame(Pearson_r = r, R_squared = r2,
    intercept = coef(fit)[1], slope = coef(fit)[2]),
    file.path(out_dir, "expression_correlation.csv"), row.names = FALSE)
  p <- ggplot(df, aes(x, y)) + geom_point(size = 0.5, alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE, color = "#D94B4B") +
    annotate("text", x = min(df$x), y = max(df$y), hjust = 0, vjust = 1,
             label = sprintf("Pearson r = %.3f\nR² = %.3f", r, r2)) +
    labs(x = paste("log1p(mean expression):", sources[1]),
         y = paste("log1p(mean expression):", sources[2])) + theme_classic()
  save_plot(p, "expression_correlation", out_dir, 6, 6)
}

# 2. Subtype markers and expression similarity --------------------------------
analyse_subtypes <- function(obj, group_col, cell_type, out_dir) {
  DefaultAssay(obj) <- "RNA"
  stopifnot(group_col %in% names(obj@meta.data),
            !anyNA(obj@meta.data[[group_col]]))
  Idents(obj) <- group_col
  stopifnot(nlevels(droplevels(Idents(obj))) >= 2)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # Preserve the different marker thresholds in the study notebooks.
  if (cell_type == "SGN") {
    markers <- FindAllMarkers(obj, only.pos = TRUE)
    if (nrow(markers)) markers <- filter(markers, avg_log2FC > 0.25)
  } else {
    markers <- FindAllMarkers(obj, only.pos = TRUE, logfc.threshold = 0.5, min.pct = 0.2)
  }
  write.csv(markers, file.path(out_dir, "markers.csv"), row.names = FALSE)
  save_plot(DimPlot(obj, group.by = group_col, label = TRUE), "subtype_UMAP", out_dir)

  # Data-driven plotting example; no predefined marker genes.
  if (nrow(markers)) {
    display <- markers %>% group_by(cluster) %>%
      slice_max(avg_log2FC, n = 3, with_ties = FALSE) %>% ungroup()
    p <- DotPlot(obj, features = unique(display$gene), group.by = group_col) +
      RotatedAxis() + labs(x = NULL, y = NULL)
    save_plot(p, "marker_dotplot_example", out_dir, 10, 5)
  }

  avg <- AverageExpression(obj, assays = "RNA", slot = "data", group.by = group_col)$RNA
  correlation <- cor(log1p(as.matrix(avg)), method = "pearson")
  if (any(!is.finite(correlation))) stop("Non-finite correlations: check subtype expression.")
  write.csv(correlation, file.path(out_dir, "expression_Pearson.csv"))
  # Same standardization and clustering conventions as the original notebooks.
  if (cell_type == "SGN") {
    scaled <- t(scale(t(correlation)))
    hc <- TRUE
  } else {
    scaled <- scale(correlation)
    hc <- hclust(as.dist(1 - correlation))
  }
  scaled[!is.finite(scaled)] <- 0
  pheatmap::pheatmap(scaled, scale = "none", cluster_rows = hc, cluster_cols = hc,
    color = viridisLite::viridis(100), breaks = seq(-1, 1, length.out = 101),
    filename = file.path(out_dir, "subtype_expression_correlation.pdf"))
  if (cell_type %in% c("Glia", "Epi")) run_gsea(obj, out_dir)
}

# 3. GO-BP GSEA: complete signed gene rankings ---------------------------------
run_gsea <- function(obj, out_dir) {
  rank_table <- FindAllMarkers(obj, assay = "RNA", slot = "data", test.use = "wilcox",
    only.pos = FALSE, logfc.threshold = 0, min.pct = 0.01, return.thresh = Inf)
  write.csv(rank_table, file.path(out_dir, "GSEA_gene_rankings.csv"), row.names = FALSE)
  if (!nrow(rank_table)) {
    message("No ranked genes available; GSEA skipped.")
    return(invisible(NULL))
  }
  sets <- msigdbr::msigdbr(species = "Mus musculus", category = "C5") %>%
    filter(grepl("^GOBP_", gs_name)) %>% distinct(gs_name, gene_symbol)
  if (!nrow(sets)) stop("No mouse GO-BP gene sets found.")
  pathways <- split(sets$gene_symbol, sets$gs_name)
  set.seed(123)
  results <- lapply(unique(as.character(rank_table$cluster)), function(cl) {
    genes <- rank_table %>% filter(cluster == cl, is.finite(avg_log2FC)) %>%
      distinct(gene, .keep_all = TRUE)
    ranks <- sort(setNames(genes$avg_log2FC, genes$gene), decreasing = TRUE)
    res <- as.data.frame(fgsea::fgseaMultilevel(
      pathways = pathways, stats = ranks, minSize = 15, maxSize = 500, scoreType = "std"))
    if (!nrow(res)) return(NULL)
    res$cluster <- cl
    res$leadingEdge <- vapply(res$leadingEdge, paste, character(1), collapse = ",")
    res
  })
  results <- bind_rows(results)
  if (!nrow(results)) {
    message("No gene sets passed the size/overlap filters; no enrichment plot generated.")
    return(invisible(NULL))
  }
  results <- arrange(results, cluster, padj)
  write.csv(results, file.path(out_dir, "GSEA_GOBP.csv"), row.names = FALSE)
  plot_gsea(results, out_dir)
}

plot_gsea <- function(results, out_dir, top_n = 5) {
  # Example display only: lowest adjusted P values; full results are retained.
  display <- results %>% filter(is.finite(pval), pval > 0, is.finite(padj)) %>%
    group_by(cluster) %>% slice_min(padj, n = top_n, with_ties = FALSE) %>%
    ungroup() %>% arrange(cluster, pval)
  if (!nrow(display)) return(invisible(NULL))
  display$id <- factor(seq_len(nrow(display)), levels = rev(seq_len(nrow(display))))
  labels <- setNames(gsub("_", " ", sub("^GOBP_", "", display$pathway)), display$id)
  p <- ggplot(display, aes(-log10(pval), id, fill = cluster)) +
    geom_col(width = 0.8) +
    facet_grid(cluster ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = labels) + labs(x = expression(-log[10](P)), y = NULL) +
    theme_minimal() + theme(legend.position = "none", panel.grid.major.y = element_blank())
  save_plot(p, "GSEA_barplot_example", out_dir, 10, 8)
}

# Example execution -----------------------------------------------------------
# A. Method comparisons: final annotated Seurat objects are required.
p1 <- readRDS("data/P1_annotated.rds")
compare_methods(p1, "celltype", c("Single-Cell", "Single-Nuclear"),
                "results/02/P1", cumulative_genes = TRUE)

adult <- readRDS("data/Merge-Cochlea-Final-V3.rds")
compare_methods(adult, "celltype_1013", c("snRNA", "StereoCell"),
                "results/02/Adult")

# B. Subtype analysis: Glia example; use the corresponding object/column for SGN/Epi.
# SGN: Merge-SGN-Final.rds, subtype_Final
# Epi: Merge-Epi-Final-V2.rds, subtype_Final
cell_type <- "Glia"
obj <- readRDS("data/Merge-Glia-Final.rds")
group_col <- "Glia_subtype"
analyse_subtypes(obj, group_col, cell_type, file.path("results/02", cell_type))
