# Infer SGN–Glia and SGN–Epi networks in separate runs.
args <- commandArgs(trailingOnly=TRUE)
context <- if (length(args)) args[1] else "SGN_Glia"
stopifnot(context %in% c("SGN_Glia","SGN_Epi"))
input_file <- if (length(args) >= 2) args[2] else
  if (context == "SGN_Glia") "data/Merge-SGN_GLia.rds" else "data/Merge-SGN_Epi.rds"
group_col <- if (length(args) >= 3) args[3] else "subtype_1004"
out_dir <- file.path("results/05", context)
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
library(Seurat)
library(CellChat)
library(dplyr)
obj <- readRDS(input_file)
options(future.globals.maxSize=10000000000)
cellchat <- createCellChat(object = obj@assays$RNA@data,
                           meta = obj@meta.data,
                           group.by = group_col)
cellchat <- setIdent(cellchat, ident.use = group_col)

cellchat@DB <- subsetDB(CellChatDB.mouse,
  search=c("ECM-Receptor", "Secreted Signaling", "Cell-Cell Contact"))


cellchat <- subsetData(cellchat)
future::plan('sequential')

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, file.path(out_dir, paste0(context, "_cellchat.rds")))


sgn.use <- c(
  "SGN-Ic","SGN-Ib",
  "SGN-Ia","SGN-II"
)

partners <- if (context == "SGN_Glia") {
  c("2.nmSC2","3.nmSC1","6.SGC","7.mSC1","9.Oligodendrocytes")
} else c("Hair Cells","IPhC_IBC","ISC_OSC")

stopifnot(
  all(c(sgn.use,partners) %in% levels(cellchat@idents))
)


lr <- subsetCommunication(cellchat)

lr.cross <- lr[
  (lr$source %in% sgn.use & lr$target %in% partners) |
  (lr$source %in% partners & lr$target %in% sgn.use),
]

write.csv(
  lr.cross,
  file.path(out_dir, paste0(context, "_bidirectional_LR.csv")),
  row.names = FALSE
)

write.csv(
  cellchat@net$weight,
  file.path(out_dir, "interaction_weights.csv")
)


net_weight <- cellchat@net$weight
net_count  <- cellchat@net$count

keep_fraction <- 0.5

positive_weights <- net_weight[
  is.finite(net_weight) & net_weight > 0
]

if (!length(positive_weights)) stop("没有可绘制的互作。")

cutoff <- unname(
  quantile(positive_weights, probs = 1 - keep_fraction)
)

keep <- is.finite(net_weight) &
        net_weight > 0 &
        net_weight >= cutoff

net_weight[!keep] <- 0
net_count[!keep]  <- 0


lr_display <- lr.cross %>% filter(is.finite(prob), is.finite(pval), prob > 0, pval < 0.05)
write.csv(lr_display, file.path(out_dir, paste0(context, "_LR_display.csv")), row.names=FALSE)
# Network display only: keep edges at or above the median positive weight.
# The LR table above is not filtered by this 50% network-display rule.
group_size <- as.numeric(table(cellchat@idents)[rownames(net_weight)])
pdf(file.path(out_dir, "communication_networks.pdf"), width=12, height=6)
par(mfrow=c(1,2))
netVisual_circle(net_count, vertex.weight=group_size, weight.scale=TRUE,
                 label.edge=FALSE, title.name="Interaction number")
netVisual_circle(net_weight, vertex.weight=group_size, weight.scale=TRUE,
                 label.edge=FALSE, title.name="Interaction strength")
dev.off()
pathways <- if (context == "SGN_Glia") c("FGF","LAMININ") else c("NRG","NRXN")
for (pw in intersect(pathways, cellchat@netP$pathways)) {
  pdf(file.path(out_dir, paste0(pw, "_heatmap.pdf")), width=7, height=6)
  ComplexHeatmap::draw(netVisual_heatmap(cellchat, signaling=pw))
  dev.off()
}
