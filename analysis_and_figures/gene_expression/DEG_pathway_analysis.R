library(fgsea)
library(gprofiler2)
library(tidyverse)
library(org.Hs.eg.db)
library(GO.db)
library(reactome.db)

# Prepare GMT files for pathway analysis
prepare_gmt_data <- function(min_size = 15, max_size = 500) {
    # Filter gene sets by size
    adjust_gene_set_size <- function(gene_set, min_size, max_size) {
        if (length(gene_set) < min_size || length(gene_set) > max_size) {
            return(NULL)
        }
        return(gene_set)
    }
    
    # Convert Ensembl IDs to symbols
    convert_ensembl_to_symbols <- function(ensembl_ids) {
        symbols <- select(org.Hs.eg.db, 
                         keys = ensembl_ids, 
                         columns = "SYMBOL", 
                         keytype = "ENSEMBL")
        unique(na.omit(symbols$SYMBOL))
    }
    
    # Filter for GO and Reactome terms
    filter_gmt_data <- function(gmt_list, prefixes = c("GO:", "REAC:")) {
        pattern <- paste0("^(", paste(prefixes, collapse = "|"), ")")
        gmt_list[names(gmt_list) %in% names(gmt_list)[grepl(pattern, names(gmt_list))]]
    }
    
    # Get pathway descriptions
    get_pathway_descriptions <- function(pathway_ids) {
        go_ids <- pathway_ids[grepl("^GO:", pathway_ids)]
        reactome_ids <- pathway_ids[grepl("^REAC:", pathway_ids)]
        
        descriptions <- rep(NA, length(pathway_ids))
        
        if (length(go_ids) > 0) {
            go_terms <- AnnotationDbi::select(GO.db, 
                                            keys = go_ids,
                                            columns = "TERM",
                                            keytype = "GOID")$TERM
            descriptions[match(go_ids, pathway_ids)] <- go_terms
        }
        
        if (length(reactome_ids) > 0) {
            reactome_clean <- gsub("REAC:R-HSA-", "", reactome_ids)
            reactome_terms <- AnnotationDbi::select(reactome.db,
                                                  keys = reactome_clean,
                                                  columns = "DEFINITION",
                                                  keytype = "PATHID")$DEFINITION
            descriptions[match(reactome_ids, pathway_ids)] <- reactome_terms
        }
        
        return(descriptions)
    }
    
    # Process and save pathway data
    pathways <- gprofiler2::get_gmt_from_url()
    pathways_filtered <- lapply(pathways, adjust_gene_set_size, min_size, max_size)
    pathways_filtered <- Filter(Negate(is.null), pathways_filtered)
    pathways_filtered <- filter_gmt_data(pathways_filtered)
    
    # Add descriptions
    descriptions <- get_pathway_descriptions(names(pathways_filtered))
    names(pathways_filtered) <- descriptions[match(names(pathways_filtered), names(descriptions))]
    
    return(pathways_filtered)
}

# Run GSEA analysis
run_gsea_analysis <- function(deg_results, pathways, output_dir = "results/pathway_analysis") {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    gsea_results <- lapply(names(deg_results), function(comparison) {
        # Prepare ranked gene list
        deg_data <- deg_results[[comparison]]
        ranks <- setNames(deg_data$coef, deg_data$primerid)
        
        # Run GSEA
        fgsea_res <- fgsea(pathways = pathways,
                          stats = ranks,
                          minSize = 15,
                          maxSize = 500)
        
        # Add comparison name and sort by NES
        fgsea_res$comparison <- comparison
        fgsea_res <- fgsea_res[order(abs(fgsea_res$NES), decreasing = TRUE), ]
        
        return(fgsea_res)
    })
    names(gsea_results) <- names(deg_results)
    
    # Save results
    saveRDS(gsea_results, file.path(output_dir, "gsea_results.rds"))
    
    # Create summary tables
    lapply(names(gsea_results), function(name) {
        write_tsv(gsea_results[[name]], 
                 file.path(output_dir, paste0("gsea_", name, ".txt")))
    })
    
    return(gsea_results)
}

# Run gProfiler analysis
run_gprofiler_analysis <- function(deg_results, output_dir = "results/pathway_analysis") {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Run analysis for each comparison
    gprofiler_results <- lapply(names(deg_results), function(comparison) {
        deg_data <- deg_results[[comparison]]
        significant_genes <- deg_data$primerid[deg_data$fdr < 0.05]
        
        gost_res <- gost(
            query = significant_genes,
            organism = "hsapiens",
            correction_method = "g_SCS",
            source = c("GO:BP", "GO:MF", "GO:CC", "REAC"),
            evcodes = TRUE
        )
        
        if (!is.null(gost_res$result)) {
            gost_res$result$comparison <- comparison
        }
        
        return(gost_res)
    })
    names(gprofiler_results) <- names(deg_results)
    
    # Save results
    saveRDS(gprofiler_results, file.path(output_dir, "gprofiler_results.rds"))
    
    # Create summary tables
    lapply(names(gprofiler_results), function(name) {
        if (!is.null(gprofiler_results[[name]]$result)) {
            write_tsv(gprofiler_results[[name]]$result,
                     file.path(output_dir, paste0("gprofiler_", name, ".txt")))
        }
    })
    
    return(gprofiler_results)
}

# Main workflow
main <- function() {
    # Read processed DEG results
    deg_results <- readRDS("results/processed_DEG/DEG_results_filtered.rds")
    
    # Prepare pathway data
    pathways <- prepare_gmt_data()
    
    # Run analyses
    gsea_results <- run_gsea_analysis(deg_results, pathways)
    gprofiler_results <- run_gprofiler_analysis(deg_results)
    
    # Return results list
    return(list(
        gsea = gsea_results,
        gprofiler = gprofiler_results
    ))
}

# Run analysis
results <- main()