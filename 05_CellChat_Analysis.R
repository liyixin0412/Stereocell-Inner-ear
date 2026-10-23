# CellChat inference and basic network/pathway visualization.
# Run SGN-Glia and SGN-Epi as separate input objects.
library(ggplot2)

infer_communication <- function(obj, group_col) {
  stopifnot(group_col %in% names(obj@meta.data), !anyNA(obj[[group_col, drop = TRUE]]))
  cellchat <- CellChat::createCellChat(
    object = Seurat::GetAssayData(obj, assay = "RNA", slot = "data"),
    meta = obj@meta.data, group.by = group_col)
  cellchat <- CellChat::setIdent(cellchat, ident.use = group_col)
  cellchat@DB <- CellChat::subsetDB(CellChat::CellChatDB.mouse,
    search = c("ECM-Receptor", "Secreted Signaling", "Cell-Cell Contact"))
  future::plan("sequential")
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  # raw.use means unprojected signaling expression, not raw UMI counts.
  cellchat <- CellChat::computeCommunProb(cellchat, raw.use = TRUE)
  cellchat <- CellChat::filterCommunication(cellchat, min.cells = 10)
  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  CellChat::aggregateNet(cellchat)
}

plot_lr <- function(lr, groups, out_dir) {
  # Use the same raw-probability and P-value definitions as the study plots.
  display <- lr[is.finite(lr$prob) & lr$prob > 0 &
                  is.finite(lr$pval) & lr$pval < 0.05, , drop = FALSE]
  if (!nrow(display)) return(invisible(NULL))
  # Illustrate one inferred pathway; no manuscript-specific pathway is fixed.
  display <- display[display$pathway_name == unique(display$pathway_name)[1], , drop = FALSE]
  probability_max <- max(display$prob)
  combinations <- unlist(lapply(groups, function(x) paste(x, groups, sep = " -> ")),
                         use.names = FALSE)
  pdf(file.path(out_dir, "ligand_receptor_bubble_example.pdf"), width = 14, height = 8,
      onefile = TRUE)
  for (pathway in unique(display$pathway_name)) {
    df <- display[display$pathway_name == pathway, , drop = FALSE]
    df$cell_pair <- factor(paste(df$source, df$target, sep = " -> "), levels = combinations)
    df$pair <- factor(df$interaction_name_2, levels = rev(unique(df$interaction_name_2)))
    df$p_group <- factor(ifelse(df$pval < 0.01, "p < 0.01", "0.01 <= p < 0.05"),
                        levels = c("p < 0.01", "0.01 <= p < 0.05"))
    p <- ggplot(df, aes(cell_pair, pair)) +
      geom_point(aes(color = prob, size = p_group)) +
      scale_color_gradientn(colors = c("#29409E", "#1A99B8", "#8FCCA3", "#FCD462", "#C71A29"),
                            limits = c(0, probability_max), name = "Raw probability") +
      scale_size_manual(values = c("p < 0.01" = 3.2, "0.01 <= p < 0.05" = 1.9),
                        drop = FALSE, name = "p-value") +
      scale_x_discrete(drop = FALSE) + labs(x = NULL, y = NULL, title = pathway) +
      theme_bw() + theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 7),
                         panel.grid.minor = element_blank())
    print(p)
  }
  dev.off()
}

export_and_plot <- function(cellchat, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # Export matrices even if no communication was retained.
  write.csv(cellchat@net$count, file.path(out_dir, "interaction_counts.csv"))
  write.csv(cellchat@net$weight, file.path(out_dir, "interaction_weights.csv"))
  if (!any(is.finite(cellchat@net$weight) & cellchat@net$weight > 0)) {
    message("No positive communication weights; matrices exported.")
    return(invisible(NULL))
  }
  lr <- CellChat::subsetCommunication(cellchat)
  write.csv(lr, file.path(out_dir, "ligand_receptor_interactions.csv"), row.names = FALSE)
  write.csv(CellChat::subsetCommunication(cellchat, slot.name = "netP"),
            file.path(out_dir, "pathway_interactions.csv"), row.names = FALSE)

  # Example 1: strongest approximately 50% of aggregate edges (counts and weights).
  weights <- cellchat@net$weight
  counts <- cellchat@net$count
  positive <- weights[is.finite(weights) & weights > 0]
  if (!length(positive)) {
    message("No positive network weights; tables exported.")
    return(invisible(NULL))
  }
  keep <- is.finite(weights) & weights > 0 & weights >= unname(quantile(positive, 0.5))
  size <- as.numeric(table(cellchat@idents)[rownames(weights)])
  pdf(file.path(out_dir, "communication_networks.pdf"), width = 12, height = 6,
      onefile = TRUE)
  weights[!keep] <- 0
  counts[!keep] <- 0
  par(mfrow = c(1, 2))
  CellChat::netVisual_circle(counts, vertex.weight = size, weight.scale = TRUE,
    label.edge = FALSE, title.name = "Top 50% edges - interaction number")
  CellChat::netVisual_circle(weights, vertex.weight = size, weight.scale = TRUE,
    label.edge = FALSE, title.name = "Top 50% edges - interaction strength")
  dev.off()

  # Example 2: one inferred pathway, without a predefined biological selection.
  if (length(cellchat@netP$pathways)) {
    pdf(file.path(out_dir, "pathway_heatmap_example.pdf"), width = 8, height = 7,
        onefile = TRUE)
    pathway <- cellchat@netP$pathways[1]
    ht <- CellChat::netVisual_heatmap(cellchat, signaling = pathway, color.heatmap = "Reds")
    ComplexHeatmap::draw(ht)
    dev.off()
  }
  # Example 3: an LR bubble plot; no aggregate-edge filter is applied here.
  plot_lr(lr, levels(cellchat@idents), out_dir)
}

# Example execution -----------------------------------------------------------
# SGN-Glia example. Run SGN-Epi separately using its own object and annotation.
obj <- readRDS("data/Merge-SGN_GLia.rds")
group_col <- "subtype_1004"
out_dir <- "results/05/SGN_Glia"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cellchat <- infer_communication(obj, group_col)
saveRDS(cellchat, file.path(out_dir, "cellchat.rds"))
export_and_plot(cellchat, out_dir)

# To plot existing results, run these lines instead of the inference block:
# cellchat <- readRDS("data/Stereocell-SGN_Glia_cellchat.rds")
# export_and_plot(cellchat, out_dir)
