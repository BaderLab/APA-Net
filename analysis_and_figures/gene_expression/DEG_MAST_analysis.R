library(tidyverse)
library(MAST)
library(Seurat)
library(magrittr)
library(future)
set.seed(1234)

# Function to run MAST analysis
run_MAST <- function(snRNA, selected_cluster, selected_disease_1, selected_disease_2,
                     output_dir = "results/MAST") {
    message(sprintf("Processing %s: %s vs %s", selected_cluster, selected_disease_1, selected_disease_2))
    
    # Subset data
    sub_idx <- which(snRNA$celltypes == selected_cluster & 
                     snRNA$disease %in% c(selected_disease_1, selected_disease_2))
    subset_data <- snRNA[, sub_idx]
    
    # Prepare data for MAST
    data_used <- GetAssayData(subset_data[["RNA"]], slot = "counts")
    data_used@x <- log2((10^6 * data_used@x / rep.int(colSums(data_used), diff(data_used@p))) + 1)
    
    # Create SCA object
    sca <- FromMatrix(
        exprsArray = Matrix::as.matrix(data_used),
        cData = data.frame(subset_data@meta.data),
        fData = data.frame(gene_name = rownames(data_used))
    )
    
    # Scale covariates
    colData(sca)$cn_genes_on <- scale(colData(sca)$nFeature_RNA)
    colData(sca)$cn_age <- scale(colData(sca)$age)
    colData(sca)$cn_nCount_RNA <- scale(colData(sca)$nCount_RNA)
    colData(sca)$cn_percent.mt <- scale(colData(sca)$percent.mt)
    
    # Filter for expressed genes
    expressed_genes <- freq(sca) > 0.1
    sca <- sca[expressed_genes, ]
    
    # Set up conditions
    cond <- factor(colData(sca)$disease)
    cond <- relevel(cond, selected_disease_2)
    colData(sca)$disease <- cond
    
    # Run MAST
    zlm_fit <- zlm(~disease + cn_genes_on + cn_nCount_RNA + cn_percent.mt + 
                    cn_age + sex + (1|sample),
                   sca,
                   method = 'glmer',
                   ebayes = FALSE,
                   strictConvergence = FALSE,
                   fitArgsD = list(nAGQ = 0))
    
    # Calculate differential expression
    lrt_term <- paste0("disease", selected_disease_1)
    summary_obj <- summary(zlm_fit, doLRT = lrt_term)
    
    # Process results
    result_dt <- summary_obj$datatable
    fcHurdle <- merge(
        result_dt[contrast==lrt_term & component=='H',.(primerid, `Pr(>Chisq)`)],
        result_dt[contrast==lrt_term & component=='logFC', .(primerid, coef, ci.hi, ci.lo)],
        by='primerid'
    )
    fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    
    # Save results
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    filename <- file.path(output_dir, 
                         sprintf("MAST_%s_%s_vs_%s.txt", 
                                selected_cluster, selected_disease_1, selected_disease_2))
    write_tsv(as_tibble(fcHurdle), filename)
    
    return(fcHurdle)
}

# Run analysis for each combination
run_celltype_analyses <- function(snRNA) {
    celltypes <- c("Oligo", "OPC", "Astro", "Micro", "Exc_upper", "Exc_int", "Exc_deep", "Inh")
    conditions <- list(
        c("C9ALS", "control"),
        c("sALS", "control"),
        c("C9ALS", "sALS")
    )
    
    for (celltype in celltypes) {
        for (cond in conditions) {
            run_MAST(snRNA, celltype, cond[1], cond[2])
        }
    }
}

run_subclass_analyses <- function(snRNA) {
    subclasses <- c("L23-IT", "L4-IT", "In-LAMP5", "In-SNCG", "In-SST", "In-VIP")
    conditions <- list(
        c("ALSFTLD", "control"),
        c("ALSnoFTLD", "control"),
        c("ALSFTLD", "ALSnoFTLD")
    )
    
    for (subclass in subclasses) {
        for (cond in conditions) {
            run_MAST(snRNA, subclass, cond[1], cond[2])
        }
    }
}

# Read processed data
# Source: Gittings et al. 2023
snRNA <- readRDS("processed/Gittings2023_processed.rds")

# Run analyses
run_celltype_analyses(snRNA)
run_subclass_analyses(snRNA)