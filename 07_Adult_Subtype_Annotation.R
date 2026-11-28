library(Seurat)

input_file <- "Adult_Fianl_metadata.rds"

# 1. Adult SGN subtype annotation ------------------------------------------------

obj <- readRDS(input_file)
obj <- subset(obj, celltype == "SGN" & source == "StereoCell")

DefaultAssay(obj) <- "integrated"
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, seed.use = 123, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, seed.use = 123)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.3, random.seed = 0)

p <- DimPlot(obj, group.by = "seurat_clusters", raster = FALSE, label = TRUE)
print(p)

# Manual subtype annotation
# Cell types were manually annotated based on
# marker genes listed in the supplementary tables.
celltype <- c(
  "0" = "SGN-Ib",
  "1" = "SGN-Ic",
  "2" = "SGN-Ia",
  "3" = "SGN-II"
)

obj$celltype <- unname(celltype[as.character(obj$seurat_clusters)])
print(table(obj$seurat_clusters, obj$celltype))

saveRDS(obj, "Adult_SGN.rds")

# 2. Adult glial subtype annotation ----------------------------------------------

obj <- readRDS(input_file)
obj <- subset(obj, celltype == "Glia" & source == "StereoCell")

DefaultAssay(obj) <- "integrated"
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, seed.use = 123, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, seed.use = 123)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.8, random.seed = 0)

p <- DimPlot(obj, group.by = "seurat_clusters", raster = FALSE, label = TRUE)
print(p)

# Manual subtype annotation
# Cell types were manually annotated based on
# marker genes listed in the supplementary tables.
celltype <- c(
  "0" = "nmSchwann2",
  "1" = "nmSchwann1",
  "2" = "SGC (satellite glial cells)",
  "3" = "mSchwann1",
  "4" = "SGC (satellite glial cells)",
  "5" = "nmSchwann2",
  "6" = "nmSchwann1",
  "7" = "Osteogenic cells",
  "8" = "nmSchwann1"
)

obj$celltype <- unname(celltype[as.character(obj$seurat_clusters)])
print(table(obj$seurat_clusters, obj$celltype))

saveRDS(obj, "Adult_Glia.rds")

# 3. Adult epithelial cell annotation --------------------------------------------

obj <- readRDS(input_file)
obj <- subset(obj, celltype %in% c("Hair Cells", "RM", "SCs", "SV-Basal", "SV-MC", "TBC") & source == "StereoCell")

DefaultAssay(obj) <- "integrated"
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, seed.use = 123, verbose = FALSE)
obj <- RunUMAP(obj, dims = 1:30, seed.use = 123)
obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = 0.3, random.seed = 0)

p <- DimPlot(obj, group.by = "seurat_clusters", raster = FALSE, label = TRUE)
print(p)

# Manual subtype annotation
# Cell types were manually annotated based on
# marker genes listed in the supplementary tables.
celltype <- c(
  "0" = "ISC_OSC",
  "1" = "SV_MC",
  "2" = "SV_MC",
  "3" = "TBC",
  "4" = "SV-Basal",
  "5" = "Spindle cells",
  "6" = "RM",
  "7" = "IPhc_IBC",
  "8" = "Hair Cells",
  "9" = "RM",
  "10" = "ISC_OSC",
  "11" = "SV_MC",
  "12" = "SV_MC"
)

obj$celltype <- unname(celltype[as.character(obj$seurat_clusters)])
print(table(obj$seurat_clusters, obj$celltype))

saveRDS(obj, "Adult_Epi.rds")
