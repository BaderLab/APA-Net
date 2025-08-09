library(Seurat)
library(future)
set.seed(1234)

setwd("/ALS_multiome")

# Define sample names and their corresponding metadata
sample_names <- c(paste0("CTRL", 1:6), paste0("C9ALSFTLD", 1:6), paste0("C9ALSnoFTLD", 1:3), paste0("sALSnoFTLD", 1:8))
snRNA_diagnosis <- c(rep("control", 6), rep("C9ALS", 9), rep("sALS", 8))
snRNA_sex <- c("female", "male", "female", "male", "female", "female",
                "male", "male", "male", "female", "female", "male",
                "female", "female", "female",
                "male", "female", "male", "female", "male", "female", "male", "male")
chemistry_labels <- c(rep("V3", 14), "V2", "V2", "V3", "V2", "V3", "V2", "V2", "V2", "V3", "V3", "V2", "V2", "V3")
snRNA_age<-c(50,59,72,72,53,48,60,60,58,47,59,72,71,66,88,59,57,65,88,43,72,37,50)

cellranger_path <- "/ALS_snRNA"
seurat_objects <- list()

for (i in seq_along(sample_names)) {
  sample_name <- sample_names[i]
  input_path <- file.path(cellranger_path, sample_name, "outs/filtered_feature_bc_matrix")

  # % mitochondrial, RPL, and RPS genes
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.rpl"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPL")
  seurat_obj[["percent.rps"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RPS")
  
  # remove high mitochondrial reads
  mito_thresh <- median(seurat_obj$percent.mt) + mad(seurat_obj$percent.mt) * 3
  drop_mito <- seurat_obj$percent.mt > mito_thresh | seurat_obj$percent.mt > 50
  seurat_obj <- seurat_obj[, !drop_mito]
  
  # subset object based on minimum and maximum number of genes
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 12000)
  
  # add metadata to subsetted object
  seurat_obj$sample <- sample_name
  seurat_obj$diagnosis <- snRNA_diagnosis[i]
  seurat_obj$sex <- snRNA_sex[i]
  seurat_obj$chemistry <- chemistry_labels[i]
  seurat_obj$age <- snRNA_age[i]
  
  # Store in the list
  seurat_objects[[sample_name]] <- seurat_obj
}

# Merge all samples for further processing
merger <- Reduce(function(x, y) merge(x, y, add.cell.ids = sample_names, project = "ALS_multiome"), seurat_objects)
DefaultAssay(merger)<-"RNA"

# save
saveRDS(merger,file="objects/snRNA/snRNA_merged_only.RDS")

