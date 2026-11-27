# Stereo-cell preprocessing, CCA integration and doublet detection.
# Edit the example input paths below and run the required section.
# No custom random seeds were set in the original analysis.
# Default Seurat v4.3.0 seeds are explicitly documented where applicable.

library(Seurat)
library(Matrix)
library(ggplot2)

# 1. Sample preprocessing ------------------------------------------------------
make_counts <- function(gem, cell_col) {
  genes <- unique(as.character(gem$geneID))
  cells <- unique(as.character(gem[[cell_col]]))
  sparseMatrix(i = match(gem$geneID, genes), j = match(gem[[cell_col]], cells), x = gem$MIDCount, dims = c(length(genes),
    length(cells)), dimnames = list(genes, cells))
}

sample_reduction <- function(obj) {
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, verbose = FALSE)
  obj <- FindClusters(obj, verbose = FALSE)
  RunUMAP(obj, dims = 1:30, verbose = FALSE)
}

preprocess_sample <- function(input_file, sample_id) {
  library(SoupX)
  gem <- data.table::fread(input_file)
  stopifnot(all(c("geneID", "label", "x", "y", "MIDCount") %in% names(gem)))
  background <- gem[gem$label == 0, ]
  foreground <- gem[gem$label != 0, ]
  stopifnot(nrow(background) > 0, nrow(foreground) > 0)
  # Background aggregation follows the original Stereo-cell notebook.
  background$cell <- paste(round(background$x/50), round(background$y/50), sep = "_")
  bg <- make_counts(background, "cell")
  counts <- make_counts(foreground, "label")
  obj <- CreateSeuratObject(counts, min.features = 100)
  obj$SN <- sample_id
  obj <- sample_reduction(obj)
  out_dir <- file.path("results/01", sample_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, file.path(out_dir, paste0(sample_id, "_raw.rds")))

  counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
  genes <- union(rownames(bg), rownames(counts))
  pad_genes <- function(x) {
    missing <- setdiff(genes, rownames(x))
    if (length(missing)) {
      zero <- sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(length(missing),
        ncol(x)), dimnames = list(missing, colnames(x)))
      x <- rbind(x, zero)
    }
    x[genes, , drop = FALSE]
  }
  bg <- pad_genes(bg)
  counts <- pad_genes(counts)
  soup <- SoupChannel(tod = cbind(counts, bg), toc = counts)
  soup <- setClusters(soup, obj@meta.data[rownames(soup$metaData), "seurat_clusters"])
  soup <- autoEstCont(soup)
  corrected <- adjustCounts(soup)
  obj <- CreateSeuratObject(corrected[Matrix::rowSums(corrected) > 0, , drop = FALSE])
  obj$SN <- sample_id
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-|^mt-")
  p <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, ncol = 3)
  ggsave(file.path(out_dir, "QC_before_filter.pdf"), p, width = 12, height = 4)
  obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA < 20000 & percent.mt <
    10)
  obj <- sample_reduction(obj)
  saveRDS(obj, file.path(out_dir, paste0(sample_id, "-QC.rds")))
  p <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, ncol = 3)
  ggsave(file.path(out_dir, "QC_after_filter.pdf"), p, width = 12, height = 4)
}

# 2. P1 integration ------------------------------------------------------------

integrate_P1 <- function() {
  files <- c(
    "results/01/Y40235F9/Y40235F9-QC.rds",
    "results/01/Y40236L9/Y40236L9-QC.rds"
  )
  samples <- c("SN", "SC")
  sources <- c("Single-Nuclear", "Single-Cell")
  resolution <- 0.3

  stopifnot(all(file.exists(files)))

  objects <- lapply(seq_along(files), function(i) {
    x <- readRDS(files[i])
    x$sample <- samples[i]
    x$source <- sources[i]
    DefaultAssay(x) <- "RNA"
    x <- NormalizeData(x, verbose = FALSE)
    FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  })

  names(objects) <- samples
  features <- SelectIntegrationFeatures(objects, nfeatures = 3000)

  # CCA integration
  # RunCCA uses seed.use = 42 by default in Seurat 4.3.0.
  anchors <- FindIntegrationAnchors(
    object.list = objects,
    anchor.features = features,
    reference = NULL,
    normalization.method = "LogNormalize",
    scale = TRUE,
    reduction = "cca",
    l2.norm = TRUE,
    dims = 1:30,
    k.anchor = 5,
    k.filter = 200,
    k.score = 30,
    max.features = 200,
    nn.method = "annoy",
    n.trees = 50,
    eps = 0
  )

  obj <- IntegrateData(
    anchorset = anchors,
    normalization.method = "LogNormalize",
    dims = 1:30,
    k.weight = 100,
    sd.weight = 1,
    eps = 0
  )

  DefaultAssay(obj) <- "integrated"
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, seed.use = 42, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:30, seed.use = 42)
  obj <- FindNeighbors(obj, dims = 1:30)
  obj <- FindClusters(obj, resolution = resolution, random.seed = 0)

  # Manual cell-type annotation
  # Cell types were manually annotated based on
  # marker genes listed in the supplementary tables.
  celltype <- c(
    "0"  = "Osteogenic cells",
    "1"  = "Immature SGN",
    "2"  = "Mesenchyme",
    "3"  = "TBC",
    "4"  = "Chondrogenic MSC",
    "5"  = "SL Fibrocytes",
    "6"  = "Glia",
    "7"  = "SV-MC",
    "8"  = "ISC_OSC",
    "9"  = "Myeloid cells",
    "10" = "Vascular Endo",
    "11" = "Ery",
    "12" = "SGN",
    "13" = "RM",
    "14" = "SV-IC",
    "15" = "Hair Cells"
  )

  obj$celltype <- unname(
    celltype[as.character(obj$seurat_clusters)]
  )

  table(obj$seurat_clusters, obj$celltype)

  out_dir <- file.path("results/01", "P1")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, file.path(out_dir, "integrated.rds"))

  for (field in c("sample", "source", "seurat_clusters")) {
    p <- DimPlot(obj, group.by = field, label = field == "seurat_clusters")
    ggsave(file.path(out_dir, paste0("UMAP_", field, ".pdf")), p, width = 7, height = 5)
  }

  # Doublet detection
  library(scDblFinder)
  library(SingleCellExperiment)

  sce <- as.SingleCellExperiment(obj, assay = "RNA")
  sce$sample <- obj$source
  # No explicit random seed was set for doublet detection.
  sce <- scDblFinder(sce, samples = "sample")

  stopifnot(identical(colnames(obj), colnames(sce)))

  obj$scDblFinder.class <- sce$scDblFinder.class
  obj$scDblFinder.score <- sce$scDblFinder.score

  write.csv(
    as.data.frame(table(obj$source, obj$scDblFinder.class)),
    file.path(out_dir, "doublet_summary.csv"),
    row.names = FALSE
  )

  p <- DimPlot(obj, group.by = "scDblFinder.class", raster = FALSE)
  ggsave(file.path(out_dir, "UMAP_doublets.pdf"), p, width = 7, height = 5)

  obj <- subset(obj, scDblFinder.class == "singlet")
  saveRDS(obj, file.path(out_dir, "singlets.rds"))
}


# 3. Adult integration ---------------------------------------------------------

integrate_Adult <- function() {
  files <- c(
    "results/02/Y01818K9/Y01818K9-QC.rds",
    "results/02/Y40237E4/Y40237E4-QC.rds",
    "results/02/3M_Cochlea_snRNA/3M_Cochlea_snRNA-QC.rds"
  )
  samples <- c("Y01818K9", "Y40237E4", "3M_sn_Cochlea")
  sources <- c("StereoCell", "StereoCell", "snRNA")
  resolution <- 0.5

  stopifnot(all(file.exists(files)))

  objects <- lapply(seq_along(files), function(i) {
    x <- readRDS(files[i])
    x$sample <- samples[i]
    x$source <- sources[i]
    DefaultAssay(x) <- "RNA"
    x <- NormalizeData(x, verbose = FALSE)
    FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  })

  names(objects) <- samples
  features <- SelectIntegrationFeatures(objects, nfeatures = 3000)

  # CCA integration
  # RunCCA uses seed.use = 42 by default in Seurat 4.3.0.
  anchors <- FindIntegrationAnchors(
    object.list = objects,
    anchor.features = features,
    reference = NULL,
    normalization.method = "LogNormalize",
    scale = TRUE,
    reduction = "cca",
    l2.norm = TRUE,
    dims = 1:30,
    k.anchor = 5,
    k.filter = 200,
    k.score = 30,
    max.features = 200,
    nn.method = "annoy",
    n.trees = 50,
    eps = 0
  )

  obj <- IntegrateData(
    anchorset = anchors,
    normalization.method = "LogNormalize",
    dims = 1:30,
    k.weight = 100,
    sd.weight = 1,
    eps = 0
  )

  DefaultAssay(obj) <- "integrated"
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, seed.use = 42, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:30, seed.use = 42)
  obj <- FindNeighbors(obj, dims = 1:30)
  obj <- FindClusters(obj, resolution = resolution, random.seed = 0)

  # Manual cell-type annotation
  # Cell types were manually annotated based on
  # marker genes listed in the supplementary tables.
  celltype <- c(
    "0"  = "TBC",
    "1"  = "Glia",
    "2"  = "SV-MC",
    "3"  = "SV-IC",
    "4"  = "SL Fibrocytes",
    "5"  = "NA MSC",
    "6"  = "RM",
    "7"  = "Vascular Endo",
    "8"  = "SCs",
    "9"  = "Glia",
    "10" = "NA MSC",
    "11" = "SV-Basal",
    "12" = "SL Fibrocytes",
    "13" = "Hair Cells",
    "14" = "SCs",
    "15" = "Osteogenic cells",
    "16" = "Fibrocytes",
    "17" = "SGN",
    "18" = "Myeloid cells",
    "19" = "Pericytes",
    "20" = "Glia",
    "21" = "RM",
    "22" = "SCs",
    "23" = "SV-IC",
    "24" = "SGN"
  )

  obj$celltype <- unname(
    celltype[as.character(obj$seurat_clusters)]
  )

  table(obj$seurat_clusters, obj$celltype)

  out_dir <- file.path("results/02", "Adult")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, file.path(out_dir, "integrated.rds"))

  for (field in c("sample", "source", "seurat_clusters")) {
    p <- DimPlot(obj, group.by = field, label = field == "seurat_clusters")
    ggsave(file.path(out_dir, paste0("UMAP_", field, ".pdf")), p, width = 7, height = 5)
  }

  # Doublet detection
  library(scDblFinder)
  library(SingleCellExperiment)

  sce <- as.SingleCellExperiment(obj, assay = "RNA")
  sce$sample <- obj$source
  # No explicit random seed was set for doublet detection.
  sce <- scDblFinder(sce, samples = "sample")

  stopifnot(identical(colnames(obj), colnames(sce)))

  obj$scDblFinder.class <- sce$scDblFinder.class
  obj$scDblFinder.score <- sce$scDblFinder.score

  write.csv(
    as.data.frame(table(obj$source, obj$scDblFinder.class)),
    file.path(out_dir, "doublet_summary.csv"),
    row.names = FALSE
  )

  p <- DimPlot(obj, group.by = "scDblFinder.class", raster = FALSE)
  ggsave(file.path(out_dir, "UMAP_doublets.pdf"), p, width = 7, height = 5)

  obj <- subset(obj, scDblFinder.class == "singlet")
  saveRDS(obj, file.path(out_dir, "singlets.rds"))
}


# 4. Example execution ---------------------------------------------------------

# A. Process each Stereo-cell sample separately by changing these two settings.
input_file <- "data/Y40235F9.gem.gz"
sample_id <- "Y40235F9"
preprocess_sample(input_file, sample_id)

# B. Run the required integration after all QC objects are available.
integrate_P1()
# integrate_Adult()
