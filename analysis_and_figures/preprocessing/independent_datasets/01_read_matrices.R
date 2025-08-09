## 01_read_matrices.R
library(Seurat)
library(scCustomize)
library(Matrix)
set.seed(1234)

# Read Li2023 dataset (frontal cortex)
li_ID <- c("A1","A2","A3","A4","A5","A6",
           "F1","F3","F4","F5",
           "C1","C2","C3","C4","C5","C6")
li_sample <- c(paste0("C9ALS", 1:6), "C9FTD1","C9FTD3","C9FTD4","C9FTD5", paste0("CTRL", 1:6))
li_diagnosis <- c(rep("C9ALS", 6), rep("C9FTD", 4), rep("CTRL", 6))
li_sex <- c("male", "female", "female", "female", "male", "male",
            "female", "male", "male", "male", 
            "female", "male", "male", "male", "female", "female")
li_age <- c(72,66,72,59,63,61,
            56,64,66,71,
            67,64,70,56,73,56)

li_files <- list.files(path = "GSE219280_RAW/frontal", pattern = "*.h5", full.names = T)
li_data <- list()

for(i in seq_len(length(li_files))) {
  counts <- Read10X_h5(li_files[i])
  li_data[[i]] <- CreateSeuratObject(counts, project = li_sample[i])
  
  li_data[[i]][["percent.mt"]] <- PercentageFeatureSet(li_data[[i]], pattern = "^MT-")
  li_data[[i]][["percent.rpl"]] <- PercentageFeatureSet(li_data[[i]], pattern = "^RPL")
  li_data[[i]][["percent.rps"]] <- PercentageFeatureSet(li_data[[i]], pattern = "^RPS")
  
  li_data[[i]] <- subset(li_data[[i]], subset = nFeature_RNA > 200 & nFeature_RNA < 12000)
  li_data[[i]]$sample_ID <- li_ID[i]
  li_data[[i]]$sex <- li_sex[i]
  li_data[[i]]$sample <- li_sample[i] 
  li_data[[i]]$diagnosis <- li_diagnosis[i]
  li_data[[i]]$age <- li_age[i]
  li_data[[i]]$study <- "Li2023"
}

li_merged <- Merge_Seurat_List(li_data, add.cell.ids = NULL, merge.data = TRUE)

# Read Pineda2024 dataset 
pineda_files <- list.dirs(path = "GSE174332_RAW", full.names = FALSE, recursive = FALSE)
pineda_data <- list()

for(i in seq_along(pineda_files)) {
  mat <- readMM(file.path(pineda_files[i], "counts_fil.mtx"))
  row_meta <- read.delim(file.path(pineda_files[i], "row_metadata.tsv"), header=TRUE, row.names=1)
  col_meta <- read.delim(file.path(pineda_files[i], "col_metadata.tsv"), header=TRUE, row.names=1)
  
  rownames(mat) <- row_meta$Gene
  colnames(mat) <- rownames(col_meta)
  
  pineda_data[[i]] <- CreateSeuratObject(mat, meta.data = col_meta)
  pineda_data[[i]]$study <- "Pineda2024"
}

pineda_merged <- Merge_Seurat_List(pineda_data, add.cell.ids = NULL, merge.data = TRUE)

saveRDS(li_merged, "li2023_merged.RDS")
saveRDS(pineda_merged, "pineda2024_merged.RDS")


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

## 05_differential_expression.R
library(Seurat)
library(tidyverse)
library(DESeq2)
library(pheatmap)
set.seed(1234)

process_dataset <- function(seurat_obj, dataset_name) {
  # Find markers for each cell type
  Idents(seurat_obj) <- "major_celltype"
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25)
  write.csv(markers, paste0(dataset_name, "_celltype_markers.csv"))
  
  # Pseudobulk analysis
  counts <- AggregateExpression(seurat_obj,
                               group.by = c("major_celltype", "sample"),
                               slot = "counts")$RNA
  
  metadata <- seurat_obj[[]]
  
  dds <- DESeqDataSetFromMatrix(countData = counts,
                               colData = metadata,
                               design = ~ diagnosis + sex)
  
  dds <- DESeq(dds)
  res <- results(dds, name="diagnosis_C9ALS_vs_control")
  write.csv(res, paste0(dataset_name, "_diagnosis_DEGs.csv"))
  
  return(list(markers = markers, degs = res))
}

# Process both datasets
li_obj <- readRDS("li2023_annotated.RDS")
li_results <- process_dataset(li_obj, "li2023")

pineda_obj <- readRDS("pineda2024_annotated.RDS")
pineda_results <- process_dataset(pineda_obj, "pineda2024")

# Compare results between datasets
common_markers <- inner_join(
  li_results$markers %>% as_tibble() %>% filter(p_val_adj < 0.05),
  pineda_results$markers %>% as_tibble() %>% filter(p_val_adj < 0.05),
  by = c("gene", "cluster")
)

common_degs <- inner_join(
  as_tibble(li_results$degs) %>% filter(padj < 0.05),
  as_tibble(pineda_results$degs) %>% filter(padj < 0.05),
  by = "gene"
)

write.csv(common_markers, "common_celltype_markers.csv")
write.csv(common_degs, "common_disease_DEGs.csv")