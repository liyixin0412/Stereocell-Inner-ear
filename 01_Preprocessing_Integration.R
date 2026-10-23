# Representative Stereo-cell preprocessing and P1 SC/SN integration.
# Run from the repository root. Annotate clusters separately before downstream analysis.
args <- commandArgs(trailingOnly = TRUE)
task <- if (length(args)) args[1] else "integrate"
dir.create("results/01", recursive=TRUE, showWarnings=FALSE)
dir.create("data", showWarnings=FALSE)
library(Seurat)
library(Matrix)
if (task == "preprocess") {
  input_file <- if (length(args) >= 2) args[2] else "data/Cell_GetExp_gene_with_background.txt"
  sample_id <- if (length(args) >= 3) args[3] else "Y40235F9"
library(Seurat)
library(data.table)
library(Matrix)
library(SoupX)
gem=fread(input_file)
gem_nocell=gem[gem$label==0,]
gem_nocell$cell= paste0( round(gem_nocell$x/50,digits = 0),"_",round(gem_nocell$y/50,digits = 0) )

gene=1:length(unique(gem_nocell$geneID))
names(gene)=unique(gem_nocell$geneID)
cell=1:length(unique(gem_nocell$cell))
names(cell)=unique(gem_nocell$cell)
mat=sparseMatrix(i = gene[gem_nocell$geneID],j = cell[ gem_nocell$cell ], x = gem_nocell$MIDCount)
rownames(mat)=names(gene)
colnames(mat)=names(cell)
gem=gem[gem$label!=0,]
gene=1:length(unique(gem$geneID))
names(gene)=unique(gem$geneID)
cell=1:length(unique(gem$label))
names(cell)=unique(gem$label)
mat1=sparseMatrix(i = gene[gem$geneID],j=cell[ as.character(gem$label) ], x=gem$MIDCount)
rownames(mat1)=names(gene)
colnames(mat1)=names(cell)

SN <- sample_id

obj1=CreateSeuratObject(counts = mat1,min.features=100)#,min.cells=10)
obj1$SN <-SN

obj1 <- NormalizeData(object = obj1,verbose = FALSE)
obj1 <- FindVariableFeatures(object = obj1,verbose = FALSE)
obj1 <- ScaleData(object = obj1,verbose = FALSE)
obj1 <- RunPCA(object = obj1,verbose = FALSE)
obj1 <- FindNeighbors(object = obj1,verbose = FALSE)
obj1 <- FindClusters(object = obj1,verbose = FALSE)
obj1 <- RunUMAP(object = obj1,dims = 1:30,verbose = FALSE)

library(SoupX)
mat1=obj1@assays$RNA@counts
genes=unique(c(rownames(mat),rownames(mat1)))
restgene=genes[!genes %in% rownames(mat)]

tmp_mat=Matrix(rep(0,length(restgene)*ncol(mat)),ncol = ncol(mat))
rownames(tmp_mat)=restgene
mat=rbind(mat,tmp_mat)
restgene=genes[!genes %in% rownames(mat1)]

tmp_mat=Matrix(rep(0,length(restgene)*ncol(mat1)),ncol = ncol(mat1))
rownames(tmp_mat)=restgene
mat1=rbind(mat1,tmp_mat)

tod=cbind(mat1[genes,],mat[genes,])

sc=SoupChannel(tod = tod,toc = mat1[genes,])
sc = setClusters(sc,obj1@meta.data[rownames( sc$metaData),"seurat_clusters"])
sc = autoEstCont(sc)
out = adjustCounts(sc)

obj2=CreateSeuratObject(counts=out[rowSums(out)>0,])
DefaultAssay(obj2) <- "RNA"
obj2[["percent.mt"]] <- PercentageFeatureSet(obj2, pattern = "^MT-|^mt-")
obj2 <- subset(obj2, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & 
                                nCount_RNA < 20000 & 
                                percent.mt < 10)
obj2 <- NormalizeData(object = obj2,verbose = FALSE)
obj2 <- FindVariableFeatures(object = obj2,verbose = FALSE)
obj2 <- ScaleData(object = obj2,verbose = FALSE)
obj2 <- RunPCA(object = obj2,verbose = FALSE)
obj2 <- FindNeighbors(object = obj2,verbose = FALSE)
obj2 <- FindClusters(object = obj2,verbose = FALSE)
obj2 <- RunUMAP(object = obj2,dims = 1:30,verbose = FALSE)

saveRDS(obj2, file.path("data", paste0(sample_id, "-QC.rds")))



} else if (task == "integrate") {
obj1 <- readRDS('data/Y40235F9-QC.rds')
obj2 <- readRDS('data/Y40236L9-QC.rds')

obj1$sample <- "SN"
obj2$sample <- "SC"

obj1$source <- "Single-Nuclear"
obj2$source <- "Single-Cell"

obj.list <- list(
  `SN` = obj1,
  `SC` = obj2
)

obj.list <- lapply(obj.list, function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  x
})

features <- SelectIntegrationFeatures(
  object.list = obj.list,
  nfeatures = 3000
)

anchors <- FindIntegrationAnchors(
  object.list = obj.list,
  anchor.features = features,
  reduction = "cca",
  dims = 1:30
)

obj.integrated <- IntegrateData(
  anchorset = anchors,
  dims = 1:30
)

DefaultAssay(obj.integrated) <- "integrated"
obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = FALSE)
obj.integrated <- RunUMAP(obj.integrated, dims = 1:30)
obj.integrated <- FindNeighbors(obj.integrated, dims = 1:30)
obj.integrated <- FindClusters(obj.integrated, resolution = 0.5)


obj <- obj.integrated
library(scDblFinder)
library(SingleCellExperiment)
sce <- as.SingleCellExperiment(obj, assay = "RNA")
sce$sample <- obj$source
sce <- scDblFinder(sce, samples = "sample")
obj$scDblFinder.class <- sce$scDblFinder.class
obj$scDblFinder.score <- sce$scDblFinder.score
obj <- subset(obj, scDblFinder.class == "singlet")
saveRDS(obj, "results/01/P1_singlets.rds")

} else stop("task must be preprocess or integrate")
