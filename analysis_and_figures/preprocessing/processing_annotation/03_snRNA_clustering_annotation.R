library(Seurat)
library(harmony)
library(future)
library(BiocParallel)
library(scDblFinder)
library(igraph)
library(leidenAlg)
library(reticulate)
Sys.setenv(RETICULATE_PYTHON="/usr/bin/python3")
library(reticulate)
sc<-import("scanpy")
set.seed(1234)

snRNA<-readRDS(file="objects/snRNA/snRNA_merged_only.RDS")

snRNA.list <- SplitObject(snRNA, split.by = "sample")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

snRNA.list <- lapply(X = snRNA.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = TRUE)
    x <- FindVariableFeatures(x, verbose = TRUE, nfeatures = 3000)
    x <- CellCycleScoring(x, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
})

features <- SelectIntegrationFeatures(object.list = snRNA.list, nfeatures = 3000)
snRNA.list <- lapply(X = snRNA.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

merged<- Merge_Seurat_List(
  snRNA.list,
  merge.data = TRUE
)

merged <- ScaleData(merged, verbose = FALSE)
merged <- FindVariableFeatures(merged, nfeatures = 3000)
merged <- RunPCA(merged, verbose = FALSE)
merged <- RunUMAP(merged, dims = 1:50)

sce<-as.SingleCellExperiment(snRNA)
sce

#by random
sce <- scDblFinder(sce, samples="sample", BPPARAM=MulticoreParam(3))
table(sce$scDblFinder.class)
snRNA2<-as.Seurat(sce)
snRNA2
table(sce$scDblFinder.class)
snRNA3<-subset(snRNA2,subset = scDblFinder.class == "singlet")
saveRDS(snRNA3, file="objects/snRNA/snRNA_scDblFinder_subset1.RDS")

snRNA3 <- ScaleData(snRNA3, verbose = TRUE, vars.to.regress = c("percent.mt","S.Score","G2M.Score"))
snRNA3 <- RunPCA(snRNA3, verbose = TRUE)
snRNA3 <- RunHarmony(snRNA3, group.by.vars = "sample", kmeans_init_nstart = 20, kmeans_init_iter_max = 100)
snRNA3 <- RunUMAP(snRNA3, dims = 1:50, umap.method = "umap-learn", metric = "correlation", reduction = "harmony", min.dist = 0.2)
snRNA3 <- FindNeighbors(snRNA3, reduction = "harmony", dims = 1:50)
snRNA1 <- FindClusters(snRNA3, resolution = 0.6, algorithm = 2)

Idents(snRNA2_sub)<-'seurat_clusters'

snRNA1_sub<-subset(snRNA1, idents=c("0","1","2","3","4","5","6","7","8","9","10",
                                    "11","12","13","14","15","16","17","18","19","20",
                                    "21","22","23","24","25","26","27","28"))

#random doublet removal clusters look the best

snRNA1_sub <- RunPCA(snRNA1_sub, verbose = TRUE)
snRNA1_sub <- RunHarmony(snRNA1_sub, group.by.vars = "sample", kmeans_init_nstart = 20, kmeans_init_iter_max = 100)
snRNA1_sub <- RunUMAP(snRNA1_sub, dims = 1:50, umap.method = "umap-learn", metric = "correlation", reduction = "harmony", min.dist = 0.2)

gc()

gc()

adata_snRNA <- sc$AnnData(
  X   = Matrix::t(GetAssayData(snRNAsub)),
  obs = snRNAsub[[]],
  var = GetAssay(snRNAsub)[[]]
)
adata_snRNA$obsm$update(umap = Embeddings(snRNAsub, "umap"))
adata_snRNA$obsm$update(harmony = Embeddings(snRNAsub, "harmony"))

sp <- import("scipy.sparse")
adata_snRNA$X  <- sp$csc_matrix(adata_snRNA$X)
adata_snRNA

adata_snRNA$write_h5ad("snRNA_NSForest1.h5ad")

#noticed 25 and 29 have both oligo and neuronal markers with low total feature/counts so remove:

snRNAsub2<-subset(snRNAsub, idents=c("1","2","3","4","5","6","7","8","9","10",
                                     "11","12","13","14","15","16","17","18","19",
                                     "20","21","22","23","24","26","27","28","29"))
snRNAsub2
rm(snRNAsub)
gc()
snRNAsub2 <- FindVariableFeatures(snRNAsub2, nfeatures=3000)
snRNAsub2 <- ScaleData(snRNAsub2, verbose = TRUE, vars.to.regress = c("percent.mt","S.Score","G2M.Score"))
gc()
snRNAsub2 <- RunPCA(snRNAsub2, verbose = TRUE)
snRNAsub2 <- RunHarmony(snRNAsub2, group.by.vars = "sample", kmeans_init_nstart = 20, kmeans_init_iter_max = 100)
snRNAsub2 <- RunUMAP(snRNAsub2, dims = 1:50, umap.method = "umap-learn", metric = "correlation", reduction = "harmony", min.dist = 0.3)
gc()
snRNAsub2 <- FindNeighbors(snRNAsub2, reduction = "harmony", dims = 1:50)
snRNAsub2 <- FindClusters(snRNAsub2, resolution = 0.8, algorithm = 4, method = "igraph")

saveRDS(snRNAsub2, file="objects/snRNA/snRNA_NSForest1.RDS")

adata_snRNA <- sc$AnnData(
  X   = Matrix::t(GetAssayData(snRNAsub2)),
  obs = snRNAsub2[[]],
  var = GetAssay(snRNAsub2)[[]]
)
adata_snRNA$obsm$update(umap = Embeddings(snRNAsub2, "umap"))
adata_snRNA$obsm$update(harmony = Embeddings(snRNAsub2, "harmony"))

sp <- import("scipy.sparse")
adata_snRNA$X  <- sp$csc_matrix(adata_snRNA$X)
adata_snRNA

adata_snRNA$write_h5ad("snRNA_NSForest1.h5ad")

##end session and restart
library(Seurat)
library(harmony)
library(future)
library(BiocParallel)
library(scDblFinder)
library(igraph)
library(leidenAlg)
library(reticulate)
Sys.setenv(RETICULATE_PYTHON="/usr/bin/python3")
library(reticulate)
sc<-import("scanpy")
set.seed(1234)

snRNA<-readRDS(file="/mnt/WORKHORSE/C9ALSFTLD_multiome/objects/snRNA/round2_NSForest.RDS")

snRNA<-subset(snRNA, idents=c("1","2","3","4","5","6","7","8","9","10",
                              "11","12","13","14","15","16","17","18","19",
                              "20","21","22","23","24","25","26","27"))
snRNA

snRNA <- RenameIdents(snRNA, `1` = "oligo1", `2` = "oligo2",  `3` = "oligo3", `4` = "astro1", 
                      `5` = "opc", `6` = "oligo4", `7` = "ex1", `8` = "micro", `9` = "ex2",
                      `10` = "in1", `11` = "astro2", `12` = "in2", `13` = "in3", `14` = "ex3", 
                      `15` = "ex4", `16` = "oligo5", `17` = "ex5", `18` = "in4", `19` = "ex6", 
                      `20` = "in5", `21` = "endo_peri", `22` = "oligo6", `23` = "ex7", `24` = "ex8", 
                      `25` = "ex9", `26` = "immune", `27` = "ex10")

Idents(snRNA) <- factor(Idents(snRNA),
                        levels = c("oligo1", "oligo2", "oligo3", "oligo4", "oligo5", "oligo6",
                                   "opc", "astro1", "astro2", "endo_peri", "micro", "immune",
                                   "ex1", "ex2", "ex3", "ex4", "ex5", "ex6", "ex7", "ex8", "ex9", "ex10",
                                   "in1", "in2", "in3", "in4", "in5"))

snRNA$cell_clusters <- Idents(snRNA)

Idents(snRNA) <- "cell_clusters"
cell_clusters2 <- c("Oligo", "Oligo", "Oligo", "Oligo", "Oligo", "Oligo",
                    "OPC", "Astro1", "Astro2", "Endo_Peri", "Micro", "Tcells",
                    "Ex1", "Ex2", "Ex3", "Ex4", "Ex5", "Ex6", "Ex7", "Ex8", "Ex9", "Ex10",
                    "In1", "In2", "In3", "In4", "In5")
names(cell_clusters2) <- levels(snRNA)
snRNA <- RenameIdents(snRNA, cell_clusters2)
snRNA$celltype_clusters <- Idents(snRNA)

celltypes <- c("Oligodendrocytes","Oligodendrocytes","Oligodendrocytes","Oligodendrocytes","Oligodendrocytes","Oligodendrocytes",
               "OPC","Astrocytes","Astrocytes","Endothelial-Pericytes","Microglia","T-cells",
               "Excitatory","Excitatory","Excitatory","Excitatory","Excitatory","Excitatory",
               "Excitatory","Excitatory","Excitatory","Excitatory",
               "Inhibitory","Inhibitory","Inhibitory","Inhibitory","Inhibitory")
names(celltypes) <- levels(snRNA)
snRNA < -RenameIdents(snRNA,celltypes)
snRNA$major_celltypes <- Idents(snRNA)
snRNA$diagnosis_celltypes <- paste0(snRNA$diagnosis, "_", snRNA$major_celltypes)

diagnosis <- c("control","C9ALSFTLD","C9ALSnoFTLD","sALSnoFTLD")
snRNA$diagnosis <- factor(snRNA$diagnosis, levels=diagnosis)

controls <- c("CTRL1","CTRL2","CTRL3","CTRL4","CTRL5","CTRL6",
              "C9ALSFTLD1","C9ALSFTLD2","C9ALSFTLD3","C9ALSFTLD4","C9ALSFTLD5","C9ALSFTLD6",
              "C9ALSnoFTLD1","C9ALSnoFTLD2","C9ALSnoFTLD3",
              "sALSnoFTLD1","sALSnoFTLD2","sALSnoFTLD3","sALSnoFTLD4","sALSnoFTLD5","sALSnoFTLD6","sALSnoFTLD7","sALSnoFTLD8")
snRNA$sample <- factor(snRNA$sample, levels = controls)

snRNA$diagnosis_celltypes <- factor(snRNA$diagnosis_celltypes, c(
  "control_Oligodendrocytes","C9ALSFTLD_Oligodendrocytes","C9ALSnoFTLD_Oligodendrocytes","sALSnoFTLD_Oligodendrocytes",
  "control_OPC","C9ALSFTLD_OPC","C9ALSnoFTLD_OPC","sALSnoFTLD_OPC",
  "control_Astrocytes","C9ALSFTLD_Astrocytes","C9ALSnoFTLD_Astrocytes","sALSnoFTLD_Astrocytes",
  "control_Endothelial-Pericytes","C9ALSFTLD_Endothelial-Pericytes","C9ALSnoFTLD_Endothelial-Pericytes","sALSnoFTLD_Endothelial-Pericytes",
  "control_Microglia","C9ALSFTLD_Microglia","C9ALSnoFTLD_Microglia","sALSnoFTLD_Microglia",
  "control_T-cells","C9ALSFTLD_T-cells","C9ALSnoFTLD_T-cells","sALSnoFTLD_T-cells",
  "control_Excitatory","C9ALSFTLD_Excitatory","C9ALSnoFTLD_Excitatory","sALSnoFTLD_Excitatory",
  "control_Inhibitory","C9ALSFTLD_Inhibitory","C9ALSnoFTLD_Inhibitory","sALSnoFTLD_Inhibitory"
))

Idents(snRNA) <- "diagnosis"
pseudobulk <- c("pseudobulk","pseudobulk","pseudobulk","pseudobulk")
names(pseudobulk) <- levels(snRNA)
snRNA <- RenameIdents(snRNA, pseudobulk)
snRNA$pseudobulk <- Idents(snRNA)

saveRDS(snRNA,file="objects/snRNA/annotated_snRNA.RDS")

adata_snRNA <- sc$AnnData(
  X   = Matrix::t(GetAssayData(snRNA)),
  obs = snRNA[[]],
  var = GetAssay(snRNA)[[]]
)
adata_snRNA$obsm$update(umap = Embeddings(snRNA, "umap"))
adata_snRNA$obsm$update(harmony = Embeddings(snRNA, "harmony"))

sp <- import("scipy.sparse")
adata_snRNA$X  <- sp$csc_matrix(adata_snRNA$X)
adata_snRNA

adata_snRNA$write_h5ad("snRNAannotated.h5ad")