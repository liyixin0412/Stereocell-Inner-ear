# Stereo-cell preprocessing, CCA integration and doublet detection.
# Edit the example input paths below and run the required section.

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

# 2. P1 or adult integration ---------------------------------------------------
integrate_samples <- function(stage) {
  if (stage == "P1") {
    files <- c("results/01/Y40235F9/Y40235F9-QC.rds", "results/01/Y40236L9/Y40236L9-QC.rds")
    samples <- c("SN", "SC")
    sources <- c("Single-Nuclear", "Single-Cell")
    resolution <- 0.5
  } else if (stage == "Adult") {
    files <- c("data/3M_Cochlea_snRNA-QC.rds", "results/01/Y01818K9/Y01818K9-QC.rds", "results/01/Y40237E4/Y40237E4-QC.rds")
    samples <- c("3M_sn_Cochlea", "Y01818K9", "Y40237E4")
    sources <- c("snRNA", "StereoCell", "StereoCell")
    resolution <- 0.3
  } else stop("stage must be P1 or Adult")
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
  anchors <- FindIntegrationAnchors(object.list = objects, anchor.features = features, reduction = "cca",
    dims = 1:30)
  obj <- IntegrateData(anchors, dims = 1:30)
  DefaultAssay(obj) <- "integrated"
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- RunUMAP(obj, dims = 1:30)
  obj <- FindNeighbors(obj, dims = 1:30)
  obj <- FindClusters(obj, resolution = resolution)
  out_dir <- file.path("results/01", stage)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(obj, file.path(out_dir, "integrated.rds"))
  for (field in c("sample", "source", "seurat_clusters")) {
    p <- DimPlot(obj, group.by = field, label = field == "seurat_clusters")
    ggsave(file.path(out_dir, paste0("UMAP_", field, ".pdf")), p, width = 7, height = 5)
  }

  # 3. Doublet detection: preserve the original source grouping.
  library(scDblFinder)
  library(SingleCellExperiment)
  sce <- as.SingleCellExperiment(obj, assay = "RNA")
  sce$sample <- obj$source
  sce <- scDblFinder(sce, samples = "sample")
  stopifnot(identical(colnames(obj), colnames(sce)))
  obj$scDblFinder.class <- sce$scDblFinder.class
  obj$scDblFinder.score <- sce$scDblFinder.score
  write.csv(as.data.frame(table(obj$source, obj$scDblFinder.class)), file.path(out_dir, "doublet_summary.csv"),
    row.names = FALSE)
  p <- DimPlot(obj, group.by = "scDblFinder.class", raster = FALSE)
  ggsave(file.path(out_dir, "UMAP_doublets.pdf"), p, width = 7, height = 5)
  obj <- subset(obj, scDblFinder.class == "singlet")
  saveRDS(obj, file.path(out_dir, "singlets.rds"))
  # Final biological annotation is supplied separately for downstream modules.
}

# Example execution -----------------------------------------------------------
# A. Process each Stereo-cell sample separately by changing these two settings.
input_file <- "data/Y40235F9.gem.gz"
sample_id <- "Y40235F9"
preprocess_sample(input_file, sample_id)

# B. Run after ALL required QC objects are available (paths in integrate_samples).
# P1: Y40235F9, Y40236L9; Adult: public 10x, Y01818K9, Y40237E4.
stage <- "P1"  # change to "Adult" for the adult comparison
integrate_samples(stage)
