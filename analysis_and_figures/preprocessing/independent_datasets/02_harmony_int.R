## 02_harmony_integration.R
library(Seurat)
library(harmony)
library(scCustomize)
set.seed(1234)

# Process Li2023 dataset
li_obj <- readRDS("li2023_merged.RDS")

li_obj <- NormalizeData(li_obj)
li_obj <- FindVariableFeatures(li_obj, nfeatures = 3000)
li_obj <- CellCycleScoring(li_obj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
li_obj <- ScaleData(li_obj, vars.to.regress = c("G2M.Score", "S.Score"))
li_obj <- RunPCA(li_obj)

li_integrated <- RunHarmony(li_obj, group.by.vars = "sample")
li_integrated <- RunUMAP(li_integrated, dims = 1:50, reduction = "harmony")
li_integrated <- FindNeighbors(li_integrated, reduction = "harmony", dims = 1:50)
li_integrated <- FindClusters(li_integrated, resolution = 0.6)

# Process Pineda2024 dataset
pineda_obj <- readRDS("pineda2024_merged.RDS")

pineda_obj <- NormalizeData(pineda_obj)
pineda_obj <- FindVariableFeatures(pineda_obj, nfeatures = 3000)
pineda_obj <- CellCycleScoring(pineda_obj, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
pineda_obj <- ScaleData(pineda_obj, vars.to.regress = c("G2M.Score", "S.Score"))
pineda_obj <- RunPCA(pineda_obj)

pineda_integrated <- RunHarmony(pineda_obj, group.by.vars = "Sample_ID")
pineda_integrated <- RunUMAP(pineda_integrated, dims = 1:50, reduction = "harmony")
pineda_integrated <- FindNeighbors(pineda_integrated, reduction = "harmony", dims = 1:50)
pineda_integrated <- FindClusters(pineda_integrated, resolution = 0.6)

saveRDS(li_integrated, "li2023_integrated.RDS")
saveRDS(pineda_integrated, "pineda2024_integrated.RDS")