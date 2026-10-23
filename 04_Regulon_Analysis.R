# Regulon activity across annotated subtypes.
# No predefined transcription factors or subtype-specific display lists.
library(Seurat)
library(ggplot2)
library(tidyr)

row_zscore <- function(x) {
  s <- apply(x, 1, sd)
  s[!is.finite(s) | s == 0] <- 1
  sweep(sweep(x, 1, rowMeans(x), "-"), 1, s, "/")
}

analyse_regulons <- function(obj, auc, group_col, out_dir) {
  stopifnot(is.matrix(auc), nrow(auc) > 0, ncol(auc) > 0,
            !is.null(rownames(auc)), !is.null(colnames(auc)),
            is.numeric(auc), all(is.finite(auc)), !anyDuplicated(rownames(auc)),
            !anyDuplicated(colnames(auc)), all(rownames(auc) %in% colnames(obj)),
            group_col %in% names(obj@meta.data))
  # Match cell barcodes before attaching subtype labels.
  meta <- obj@meta.data[rownames(auc), , drop = FALSE]
  stopifnot(identical(rownames(meta), rownames(auc)), !anyNA(meta[[group_col]]))
  groups <- if (is.factor(meta[[group_col]])) levels(droplevels(meta[[group_col]])) else
    unique(as.character(meta[[group_col]]))
  stopifnot(length(groups) >= 2)
  mean_auc <- do.call(cbind, lapply(groups, function(g)
    colMeans(auc[meta[[group_col]] == g, , drop = FALSE])))
  dimnames(mean_auc) <- list(colnames(auc), groups)
  mean_z <- row_zscore(mean_auc)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(mean_auc, file.path(out_dir, "subtype_mean_AUC.csv"))
  write.csv(mean_z, file.path(out_dir, "subtype_relative_activity.csv"))

  # Example 1: full mean-AUC and relative-activity heatmaps.
  annotation <- data.frame(Subtype = factor(groups, levels = groups), row.names = groups)
  colors <- setNames(grDevices::hcl.colors(length(groups), "Dark 3"), groups)
  height <- max(5, 2 + 0.15 * nrow(mean_auc))
  auc_range <- range(mean_auc)
  if (diff(auc_range) == 0) auc_range <- auc_range + c(-1, 1) * 1e-6
  pheatmap::pheatmap(mean_auc, cluster_rows = FALSE, cluster_cols = FALSE,
    breaks = seq(auc_range[1], auc_range[2], length.out = 101),
    annotation_col = annotation, annotation_colors = list(Subtype = colors),
    border_color = NA, color = colorRampPalette(c("#FFF7EC", "#FC8D59", "#7F0000"))(100),
    filename = file.path(out_dir, "mean_AUC_heatmap.pdf"), width = 8, height = height)
  pheatmap::pheatmap(mean_z, cluster_rows = FALSE, cluster_cols = FALSE,
    annotation_col = annotation, annotation_colors = list(Subtype = colors),
    border_color = NA, color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(100),
    breaks = seq(-1, 1, length.out = 101),
    filename = file.path(out_dir, "relative_activity_heatmap.pdf"), width = 8, height = height)

  # Example 2: four variable regulons; this is a plotting example, not TF selection.
  variance <- apply(mean_auc, 1, var)
  variable <- colnames(auc)[apply(auc, 2, var) > 0]
  features <- head(names(sort(variance[variable], decreasing = TRUE)), 4)
  if (!length(features)) {
    message("AUC values are constant; heatmaps saved and ridge plot skipped.")
    return(invisible(mean_auc))
  }
  df <- as.data.frame(auc[, features, drop = FALSE])
  df$subtype <- factor(meta[[group_col]], levels = groups)
  df <- pivot_longer(df, cols = all_of(features), names_to = "Regulon", values_to = "AUC")
  df$Regulon <- factor(df$Regulon, levels = features)
  p <- ggplot(df, aes(AUC, subtype, fill = subtype)) +
    ggridges::geom_density_ridges(alpha = 0.8, scale = 1.2, color = "white") +
    facet_wrap(~Regulon, ncol = 2, scales = "free_x") +
    scale_fill_manual(values = colors) + scale_y_discrete(limits = rev(groups)) +
    labs(x = "Regulon AUC", y = NULL) + theme_classic() +
    theme(legend.position = "none", strip.text = element_text(face = "bold"))
  ggsave(file.path(out_dir, "regulon_ridges_example.pdf"), p, width = 8, height = 7)
  invisible(mean_auc)
}

# Example execution -----------------------------------------------------------
# Use the same lineage and cell barcodes as in module 03.
cell_type <- "Glia"
obj <- readRDS("data/Merge-Glia-Final.rds")
group_col <- "Glia_subtype"
auc_file <- file.path("results/03", cell_type, paste0(cell_type, "_auc_matrix.csv"))
auc <- as.matrix(read.csv(auc_file, row.names = 1, check.names = FALSE))
analyse_regulons(obj, auc, group_col, file.path("results/04", cell_type))
# For SGN/Epi, use their annotated RDS and subtype_Final column.
