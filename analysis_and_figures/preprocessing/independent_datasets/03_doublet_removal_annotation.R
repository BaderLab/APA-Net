## 03_doublet_removal.R
library(Seurat)
library(scDblFinder)
library(BiocParallel)
set.seed(1234)

# Process Li2023 dataset
li_obj <- readRDS("li2023_integrated.RDS")
li_sce <- as.SingleCellExperiment(li_obj)

li_sce_rand <- scDblFinder(li_sce, samples="sample", BPPARAM=MulticoreParam(3))
li_rand <- as.Seurat(li_sce_rand)
li_rand_clean <- subset(li_rand, subset = scDblFinder.class == "singlet")

li_sce_clust <- scDblFinder(li_sce, samples="sample", BPPARAM=MulticoreParam(3), clusters=TRUE)
li_clust <- as.Seurat(li_sce_clust)
li_clust_clean <- subset(li_clust, subset = scDblFinder.class == "singlet")

# Process Pineda2024 dataset  
pineda_obj <- readRDS("pineda2024_integrated.RDS")
pineda_sce <- as.SingleCellExperiment(pineda_obj)

pineda_sce_rand <- scDblFinder(pineda_sce, samples="Sample_ID", BPPARAM=MulticoreParam(3))
pineda_rand <- as.Seurat(pineda_sce_rand)
pineda_rand_clean <- subset(pineda_rand, subset = scDblFinder.class == "singlet")

pineda_sce_clust <- scDblFinder(pineda_sce, samples="Sample_ID", BPPARAM=MulticoreParam(3), clusters=TRUE)  
pineda_clust <- as.Seurat(pineda_sce_clust)
pineda_clust_clean <- subset(pineda_clust, subset = scDblFinder.class == "singlet")

saveRDS(li_rand_clean, "li2023_doublets_removed_random.RDS")
saveRDS(li_clust_clean, "li2023_doublets_removed_clusters.RDS")
saveRDS(pineda_rand_clean, "pineda2024_doublets_removed_random.RDS")
saveRDS(pineda_clust_clean, "pineda2024_doublets_removed_clusters.RDS")

## 04_cell_type_annotation.R
library(Seurat)
library(scCustomize)
set.seed(1234)

# Annotate Li2023 dataset
li_obj <- readRDS("li2023_doublets_removed_random.RDS")

li_celltypes <- c("Oligo1", "Oligo2", "L23_IT_1", "Astro1", "L23_IT_2", 
                  "Oligo3", "OPC", "Micro_PVM", "IN_VIP", "L5_IT_1",
                  "L4_IT", "IN_SST", "L6_IT", "IN_PVALB", "Astro2",
                  "Endothelial", "L6_CT", "L6b", "IN_LAMP5", "L56_NP")

names(li_celltypes) <- levels(li_obj)
li_obj <- RenameIdents(li_obj, li_celltypes)
li_obj$celltype <- Idents(li_obj)

li_major_types <- c("Oligodendrocytes", "Oligodendrocytes", "Excitatory", "Astrocytes",
                    "Excitatory", "Oligodendrocytes", "OPC", "Microglia", "Inhibitory",
                    "Excitatory", "Excitatory", "Inhibitory", "Excitatory", "Inhibitory", 
                    "Astrocytes", "Endothelial", "Excitatory", "Excitatory", "Inhibitory",
                    "Excitatory")

names(li_major_types) <- levels(li_obj)
li_obj <- RenameIdents(li_obj, li_major_types)
li_obj$major_celltype <- Idents(li_obj)

# Annotate Pineda2024 dataset
pineda_obj <- readRDS("pineda2024_doublets_removed_random.RDS")

pineda_celltypes <- c("Oligo1", "Oligo2", "L23_IT_1", "Astro1", "L23_IT_2",
                      "Oligo3", "OPC", "Micro_PVM", "IN_VIP", "L5_IT_1",
                      "L4_IT", "IN_SST", "L6_IT", "IN_PVALB", "Astro2",
                      "Endothelial", "L6_CT", "L6b", "IN_LAMP5", "L56_NP")

names(pineda_celltypes) <- levels(pineda_obj)
pineda_obj <- RenameIdents(pineda_obj, pineda_celltypes)
pineda_obj$celltype <- Idents(pineda_obj)

pineda_major_types <- c("Oligodendrocytes", "Oligodendrocytes", "Excitatory", "Astrocytes",
                        "Excitatory", "Oligodendrocytes", "OPC", "Microglia", "Inhibitory",
                        "Excitatory", "Excitatory", "Inhibitory", "Excitatory", "Inhibitory",
                        "Astrocytes", "Endothelial", "Excitatory", "Excitatory", "Inhibitory",
                        "Excitatory")

names(pineda_major_types) <- levels(pineda_obj)
pineda_obj <- RenameIdents(pineda_obj, pineda_major_types)
pineda_obj$major_celltype <- Idents(pineda_obj)

saveRDS(li_obj, "li2023_annotated.RDS")
saveRDS(pineda_obj, "pineda2024_annotated.RDS")