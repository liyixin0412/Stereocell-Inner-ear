library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

gene_df <- read.csv(
  "Hearing_loss_genes_HL-adult-P1.csv",
  stringsAsFactors = FALSE
)

to_mouse <- function(x) {
  paste0(
    toupper(substr(x, 1, 1)),
    tolower(substr(x, 2, nchar(x)))
  )
}

disease_colors <- c(
  "Age-related HLG (ARHL)" = "#F0F39A",
  "Auditory Neuropathy" = "#FB8072",
  "Autosomal Dominant Nonsyndromic HL" = "#CCEBC5",
  "Autosomal Recessive Nonsyndromic HL" = "#BC80BD",
  "Cisplatin-related HL (CRHL)" = "#BEBADA",
  "Mitochondrial Nonsyndromic HL" = "#FDB462",
  "Sex-Linked Nonsyndromic HL" = "#80B1D3",
  "Syndromic HL" = "#8DD3C7"
)

get_gene_info <- function(obj, gene_column, disease_column) {
  gene_df %>%
    filter(
      !is.na(.data[[gene_column]]),
      .data[[gene_column]] != "",
      !is.na(.data[[disease_column]]),
      .data[[disease_column]] != ""
    ) %>%
    distinct(
      .data[[disease_column]],
      .data[[gene_column]],
      .keep_all = TRUE
    ) %>%
    transmute(
      Disease = .data[[disease_column]],
      Gene = to_mouse(.data[[gene_column]])
    ) %>%
    filter(Gene %in% rownames(obj[["RNA"]]))
}

plot_hearing_loss_heatmap <- function(
  obj,
  gene_info,
  group_by,
  markers,
  title,
  output
) {
  genes <- unique(gene_info$Gene)

  DefaultAssay(obj) <- "RNA"

  avg_mat <- AverageExpression(
    obj,
    assays = "RNA",
    features = genes,
    group.by = group_by,
    slot = "data"
  )$RNA

  avg_mat <- avg_mat[genes, , drop = FALSE]

  mat_scaled <- t(scale(t(avg_mat)))
  mat_scaled[is.na(mat_scaled)] <- 0
  mat_scaled <- mat_scaled[genes, , drop = FALSE]

  gene_anno <- gene_info %>%
    distinct(Gene, Disease, .keep_all = TRUE)

  gene_anno <- gene_anno[
    match(rownames(mat_scaled), gene_anno$Gene),
    ,
    drop = FALSE
  ]

  row_anno <- data.frame(
    Disease = gene_anno$Disease
  )
  rownames(row_anno) <- gene_anno$Gene

  ha_left <- rowAnnotation(
    Disease = row_anno$Disease,
    col = list(Disease = disease_colors),
    show_annotation_name = FALSE
  )

  genes_label <- unique(c(
    markers$marker1,
    markers$marker2
  ))

  mark_idx <- which(
    rownames(mat_scaled) %in% genes_label
  )

  mark_genes <- rownames(mat_scaled)[mark_idx]

  mark_cols <- ifelse(
    mark_genes %in% markers$marker2,
    "red",
    "black"
  )

  ha_right <- rowAnnotation(
    mark = anno_mark(
      at = mark_idx,
      labels = mark_genes,
      labels_gp = gpar(
        fontsize = 10,
        col = mark_cols
      ),
      link_gp = gpar(
        col = mark_cols,
        lwd = 0.8
      ),
      link_width = unit(5, "mm"),
      padding = unit(2, "mm")
    )
  )

  ht <- Heatmap(
    mat_scaled,
    name = "Z-score",
    col = colorRamp2(
      c(-2, 0, 2),
      c("#00A6E8", "#FFFFFF", "#FF5A36")
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = TRUE,
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 10),
    rect_gp = gpar(
      col = "grey60",
      lwd = 0.5
    ),
    left_annotation = ha_left,
    right_annotation = ha_right,
    column_title = title
  )

  pdf(output, width = 13, height = 10)

  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right"
  )

  dev.off()
}

# Adult cochlea
adult_obj <- readRDS("Merge-Cochlea-Final-V3.rds")

adult_cell_order <- c(
  "SGN",
  "SL Fibrocytes",
  "Vascular Endo",
  "NA MSC",
  "SV-IC",
  "Glia",
  "Fibrocytes",
  "Pericytes",
  "Myeloid cells",
  "Osteogenic cells",
  "SV-MC",
  "ISC_OSC",
  "SV-Basal",
  "Spindle cells",
  "RM",
  "IPhC_IBC",
  "Hair Cells",
  "TBC"
)

adult_obj$celltype_1013_Epi <- factor(
  adult_obj$celltype_1013_Epi,
  levels = adult_cell_order
)

adult_markers <- list(
  marker1 = c(
    "Ccbe1", "Cntnap2", "Fat2", "Arhgap35", "Sh2d4b",
    "Usp48", "Ncoa3", "Pde1c", "Pi4kb", "Slc22a4",
    "Gas2", "Gosr2", "Bdp1", "Msrb3", "Eps8l2",
    "Col4a6", "Aifm1", "Smpx", "Mepe", "Foxl1",
    "Gjb2", "Gjb6"
  ),
  marker2 = "Ndp"
)

adult_gene_info <- get_gene_info(
  adult_obj,
  gene_column = "Adult",
  disease_column = "Adult_Disease"
)

plot_hearing_loss_heatmap(
  obj = adult_obj,
  gene_info = adult_gene_info,
  group_by = "celltype_1013_Epi",
  markers = adult_markers,
  title = "Adult Hearing Loss Genes",
  output = "results/Merge-Adult_DiseaseGene_mark.pdf"
)

# P1 cochlea
p1_obj <- readRDS("P1_去双胞_260907.rds")

p1_cell_order <- c(
  "Osteogenic cells",
  "Immature SGN",
  "Mesenchyme",
  "Chondrogenic MSC",
  "SL Fibrocytes",
  "Glia",
  "Myeloid cells",
  "Vascular Endo",
  "Ery",
  "SGN",
  "SV-IC",
  "ISC_OSC",
  "SV-MC",
  "RM",
  "Hair Cells",
  "TBC"
)

p1_obj$celltype <- factor(
  p1_obj$celltype,
  levels = p1_cell_order
)

p1_markers <- list(
  marker1 = c(
    "Gjb2", "Slc12a2", "Kcnq4", "Slc26a5", "Otof",
    "Tmc1", "Myo7a", "Pcdh15", "Tecta", "Pou4f3"
  ),
  marker2 = "Ndp"
)

p1_gene_info <- get_gene_info(
  p1_obj,
  gene_column = "P1",
  disease_column = "P1_Disease"
)

plot_hearing_loss_heatmap(
  obj = p1_obj,
  gene_info = p1_gene_info,
  group_by = "celltype",
  markers = p1_markers,
  title = "P1 Hearing Loss Genes",
  output = "results/Merge-P1_DiseaseGene_mark.pdf"
)
