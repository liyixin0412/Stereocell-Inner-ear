# Mean regulon activity, representative figures and optional LR target overlap.
library(Seurat)
library(ggplot2)
library(ggridges)
library(tidyr)
args <- commandArgs(trailingOnly=TRUE)
cell_type <- if (length(args)) args[1] else "Epi"
paths <- c(SGN="data/Merge-SGN-Final.rds", Glia="data/Merge-Glia-Final.rds",
           Epi="data/Merge-Epi-Final-V2.rds")
stopifnot(cell_type %in% names(paths))
obj <- readRDS(paths[[cell_type]])
group_col <- if (cell_type == "Glia") "Glia_subtype" else "subtype_Final"
# Set group_col="subtype_Final_new" if that is the final SGN annotation in your object.
out_dir <- file.path("results/04", cell_type)
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
auc <- as.matrix(read.csv(file.path("results/03", cell_type,
  paste0(cell_type, "_auc_matrix.csv")), row.names=1, check.names=FALSE))
stopifnot(!anyDuplicated(rownames(auc)), all(rownames(auc) %in% colnames(obj)),
          is.numeric(auc), all(is.finite(auc)))
meta <- obj@meta.data[rownames(auc), , drop=FALSE]
stopifnot(group_col %in% names(meta), !anyNA(meta[[group_col]]))
if (cell_type == "SGN") {
group_colors <- c(
  "SGN-Ia"       = "#FDB45B",
  "SGN-Ib"      = "#68C1D9",
  "SGN-Ic" = "#B8DFAF",
  "SGN-II"            = "#D92E92"
)

regulon_order <- c(
    "Atf1(+)", "Crem(+)", "Mafb(+)", "Rxrg(+)",
    "Etv5(+)", "Ing4(+)", "Ikzf2(+)",
    "Creb1(+)", "Rest(+)", "Stat3(+)", "Nfatc2(+)",
    "Cux1(+)", "Irx1(+)", 
    "Foxp2(+)", "Mafg(+)", "Irf3(+)","Gata3(+)", 
    "Meis2(+)", "Tcf7l2(+)","Runx2(+)"
)

} else if (cell_type == "Glia") {
group_colors <- c(
  "mSC1"         = "#FFE66D",
  "nmSC1"       = "#FDB45B",
  "nmSC2"      = "#68C1D9",
  "Oligodendrocytes" = "#B8DFAF",
  "SGC"            = "#D92E92"
)

# 手动设置 regulon 顺序
regulon_order <- c(
     "Nr2f2(+)", "Zfp558(+)", "Polr3a(+)",
    "Sox17(+)",
    "Bach2(+)", "Kdm5a(+)",  "Sox6(+)", 
    "Sox8(+)", "Srf(+)", "Zfp108(+)", "Klf15(+)", "Mzf1(+)",
    "Prdm16(+)","Mybl1(+)"
)


} else if (cell_type == "Epi") {
group_colors <- c(
  "SV-MC"         = "#FFE66D",
  "ISC_OSC"       = "#FDB45B",
  "SV-Basal"      = "#68C1D9",
  "Spindle cells" = "#B8DFAF",
  "RM"            = "#D92E92",
  "IPhC_IBC"      = "#9FC5E8",
  "Hair Cells"    = "#C7B5E3",
  "TBC"           = "#A9D18E"
)

# 手动设置 regulon 顺序
regulon_order <- c(
  "Hmx2(+)",
  "Arnt2(+)",
  "Gtf2i(+)",
  "Pbx3(+)",
  "Pou3f1(+)",
  "Irf8(+)",
  "Rorb(+)",
  "Zmat4(+)",
  "6720489N17Rik(+)",
  "Brf1(+)"
)


}
group_order <- names(group_colors)
stopifnot(setequal(as.character(meta[[group_col]]), group_order))
mean_auc <- sapply(group_order, function(g) colMeans(auc[meta[[group_col]] == g, , drop=FALSE]))
rownames(mean_auc) <- colnames(auc)
row_zscore <- function(x) {
  s <- apply(x, 1, sd)
  s[!is.finite(s) | s == 0] <- 1
  sweep(sweep(x, 1, rowMeans(x), "-"), 1, s, "/")
}
write.csv(mean_auc, file.path(out_dir, "subtype_mean_AUC.csv"))
write.csv(row_zscore(mean_auc), file.path(out_dir, "subtype_relative_activity.csv"))
# Lists copied from the supplied visualization notebook; edit for final panel selection.
stopifnot(all(regulon_order %in% rownames(mean_auc)))
pheatmap::pheatmap(row_zscore(mean_auc[regulon_order, , drop=FALSE]),
  cluster_rows=FALSE, cluster_cols=FALSE, border_color=NA,
  color=colorRampPalette(c("#2166AC","#F7F7F7","#B2182B"))(100),
  breaks=seq(-1,1,length.out=101),
  filename=file.path(out_dir, "regulon_relative_activity.pdf"), width=7, height=8)

ridge_regulons <- switch(cell_type,
  SGN=c("Rxrg(+)","Zfp708(+)","Nfatc2(+)","Creb1(+)","Etv5(+)","Gata3(+)",
        "Cux1(+)","Grhl1(+)","Zbtb7c(+)","Nr2c2(+)"),
  Glia=c("Mybl1(+)","Zfp558(+)","Sox8(+)","Sox6(+)","Bach2(+)","Kdm5a(+)"),
  Epi=c("Arnt2(+)","Pbx3(+)","Rorb(+)","Zmat4(+)"))
stopifnot(all(ridge_regulons %in% colnames(auc)))
df <- as.data.frame(auc[, ridge_regulons, drop=FALSE])
df$subtype <- factor(meta[[group_col]], levels=group_order)
df <- pivot_longer(df, cols=all_of(ridge_regulons), names_to="Regulon", values_to="AUC")
df$Regulon <- factor(df$Regulon, levels=ridge_regulons)
p <- ggplot(df, aes(AUC, subtype, fill=subtype)) +
  geom_density_ridges(alpha=0.8, scale=1.2, color="white") +
  facet_wrap(~Regulon, ncol=2, scales="free_x") +
  scale_fill_manual(values=group_colors) +
  scale_y_discrete(limits=rev(group_order)) +
  theme_classic() + theme(legend.position="none") + labs(x="Regulon AUC", y=NULL)
ggsave(file.path(out_dir, "regulon_ridge.pdf"), p, width=8,
       height=max(5, 2*ceiling(length(ridge_regulons)/2)))

# Optional RSS: only enable if SCENIC is installed; not required for the mean-AUC panels.
run_rss <- FALSE
if (run_rss) {
  annotation <- setNames(as.character(meta[[group_col]]), rownames(meta))
  rss <- SCENIC::calcRSS(AUC=t(auc), cellAnnotation=annotation)
  write.csv(rss, file.path(out_dir, "RSS_matrix.csv"))
}

# Optional connection to existing CellChat results: exact gene-set overlap.
# Supply the corresponding files in data/; SGN–Glia and SGN–Epi stay separate.
target_file <- file.path("results/03", cell_type, paste0(cell_type, "_regulon_targets.csv"))
lr_files <- c(SGN_Glia="data/SGN_Glia_bidirectional_LR.csv",
              SGN_Epi="data/SGN_Epi_bidirectional_LR.csv")
if (cell_type == "Glia") lr_files <- lr_files["SGN_Glia"]
if (cell_type == "Epi") lr_files <- lr_files["SGN_Epi"]
if (file.exists(target_file)) {
  targets <- read.csv(target_file, stringsAsFactors=FALSE)
  # Edit if your annotation uses different labels. Matching is case-sensitive.
  groups <- switch(cell_type,
    SGN=c("SGN-Ia","SGN-Ib","SGN-Ic","SGN-II"),
    Glia=c("nmSC1","nmSC2","mSC1","SGC","Oligodendrocytes",
           "2.nmSC2","3.nmSC1","6.SGC","7.mSC1","9.Oligodendrocytes"),
    Epi=c("Hair Cells","IPhC_IBC","ISC_OSC"))
  for (context in names(lr_files)) {
    if (!file.exists(lr_files[[context]])) next
    lr <- read.csv(lr_files[[context]], stringsAsFactors=FALSE)
    stopifnot(all(c("source","target","ligand","receptor","prob","pval","pathway_name") %in% names(lr)))
    lr <- lr[is.finite(lr$prob) & is.finite(lr$pval) & lr$prob > 0 & lr$pval < 0.05, ]
    hits <- list()
    for (role in c("ligand","receptor")) {
      side <- if (role == "ligand") "source" else "target"
      for (i in which(lr[[side]] %in% groups)) {
        # Split multimeric complexes into gene components.
        genes <- strsplit(lr[[role]][i], "[_+&(),[:space:]]+")[[1]]
        h <- targets[targets$Target %in% genes, , drop=FALSE]
        if (!nrow(h)) next
        h$role <- role
        h$source <- lr$source[i]; h$receiver <- lr$target[i]
        h$pathway <- lr$pathway_name[i]
        h$LR_ligand <- lr$ligand[i]; h$LR_receptor <- lr$receptor[i]
        hits[[length(hits)+1]] <- h
      }
    }
    if (length(hits)) write.csv(unique(do.call(rbind, hits)),
      file.path(out_dir, paste0(context, "_TF_LR_overlap.csv")), row.names=FALSE)
  }
}
# Overlap indicates candidate association, not TF causality or subtype-specific regulation.
