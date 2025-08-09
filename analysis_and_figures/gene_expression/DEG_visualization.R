library(tidyverse)
library(pheatmap)
library(EnhancedVolcano)
library(circlize)
set.seed(1234)

#' Process DEG results for comparison
#' @param deg_data List of DEG results from MAST analysis
#' @param external_data List of external dataset results
#' @param z_score_cutoff Cutoff for z-score filtering
process_comparison_data <- function(deg_data, external_data, z_score_cutoff = 1.5) {
    # Calculate z-scores for internal data
    Z <- 1.96
    
    internal_data <- deg_data %>%
        mutate(SE = (ci.hi - ci.lo) / (2 * Z),
               z_score = model_log2FC / SE)
    
    # Process external data
    external_processed <- map(external_data, function(dataset) {
        if("ci.hi" %in% names(dataset) & "ci.lo" %in% names(dataset)) {
            dataset %>%
                mutate(SE = (ci.hi - ci.lo) / (2 * Z),
                       z_score = model_log2FC / SE)
        } else {
            dataset
        }
    })
    
    return(list(internal = internal_data, 
                external = external_processed))
}

#' Generate Heatmap for Cross-Study Comparison
#' @param data_list List of processed DEG results
#' @param cell_type Cell type to analyze
#' @param comparison_type Type of comparison (e.g., "C9ALS_vs_control")
generate_comparison_heatmap <- function(data_list, cell_type, comparison_type) {
    # Combine data and find common genes
    common_genes <- Reduce(intersect, 
                          lapply(data_list, function(x) x$gene))
    
    # Prepare heatmap data
    heatmap_data <- map_dfr(names(data_list), function(study) {
        data_list[[study]] %>%
            filter(gene %in% common_genes) %>%
            select(gene, z_score) %>%
            spread(key = study, value = z_score)
    })
    
    # Create heatmap matrix
    heatmap_matrix <- as.matrix(heatmap_data[,-1])
    rownames(heatmap_matrix) <- heatmap_data$gene
    
    # Define color scheme
    colors <- colorRamp2(
        c(min(heatmap_matrix), 0, max(heatmap_matrix)),
        c("blue", "white", "red")
    )
    
    # Generate heatmap
    pheatmap(heatmap_matrix,
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             scale = "none",
             color = colors(100),
             main = paste("DEG Comparison -", cell_type, comparison_type),
             fontsize = 10,
             filename = paste0("results/heatmaps/", 
                             cell_type, "_", 
                             comparison_type, 
                             "_comparison.pdf"))
}

#' Generate Volcano Plots
#' @param deg_results DEG analysis results
#' @param cell_types Vector of cell types to analyze
#' @param comparisons Vector of comparisons to analyze
generate_volcano_plots <- function(deg_results, cell_types, comparisons) {
    # Create output directory
    dir.create("results/volcanos", recursive = TRUE, showWarnings = FALSE)
    
    # Generate plots for each cell type and comparison
    walk(cell_types, function(cell_type) {
        walk(comparisons, function(comparison) {
            data <- deg_results %>%
                filter(celltype == cell_type, 
                       comparison == !!comparison)
            
            if(nrow(data) > 0) {
                pdf(paste0("results/volcanos/", 
                          cell_type, "_", 
                          comparison, 
                          "_volcano.pdf"),
                    width = 10, height = 10)
                
                plot <- EnhancedVolcano(data,
                    lab = data$gene,
                    x = 'avg_log2FC',
                    y = 'fdr',
                    pCutoff = 0.01,
                    FCcutoff = 0.5,
                    title = paste(comparison, "in", cell_type),
                    subtitle = '',
                    caption = '',
                    labSize = 0,
                    pointSize = 3.0,
                    colAlpha = 1,
                    legendPosition = 'top')
                
                print(plot)
                dev.off()
            }
        })
    })
}

#' Compare DEGs Across Studies
#' @param internal_data Internal DEG results
#' @param external_datasets List of external datasets
#' @param cell_type Cell type to analyze
compare_deg_lists <- function(internal_data, external_datasets, cell_type) {
    # Get DEGs from each dataset
    deg_lists <- list(
        internal = internal_data %>%
            filter(celltype == cell_type, fdr < 0.05) %>%
            pull(gene),
        external = map(external_datasets, function(dataset) {
            dataset %>%
                filter(celltype == cell_type, fdr < 0.05) %>%
                pull(gene)
        })
    )
    
    # Find overlapping genes
    overlap <- Reduce(intersect, deg_lists)
    
    # Calculate statistics
    stats <- list(
        overlap_size = length(overlap),
        overlap_genes = overlap,
        individual_sizes = map_int(deg_lists, length)
    )
    
    return(stats)
}

#' Main Visualization Workflow
#' @param deg_results Internal DEG results
#' @param external_data List of external datasets
#' @param cell_types Vector of cell types to analyze
#' @param comparisons Vector of comparisons to analyze
main_visualization <- function(deg_results, external_data, cell_types, comparisons) {
    # Process data
    processed_data <- process_comparison_data(deg_results, external_data)
    
    # Generate heatmaps
    walk(cell_types, function(ct) {
        walk(comparisons, function(comp) {
            generate_comparison_heatmap(processed_data, ct, comp)
        })
    })
    
    # Generate volcano plots
    generate_volcano_plots(processed_data$internal, cell_types, comparisons)
    
    # Compare DEGs across studies
    comparison_results <- map(cell_types, function(ct) {
        compare_deg_lists(processed_data$internal, 
                         processed_data$external, 
                         ct)
    })
    names(comparison_results) <- cell_types
    
    return(comparison_results)
}