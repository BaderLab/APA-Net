library(tidyverse)
library(gprofiler2)
library(data.table)
library(org.Hs.eg.db)
library(GO.db)
library(reactome.db)
set.seed(1234)

process_apa_data <- function(directory = "data/APA", output_dir = "results/APA") {
    # Read APA result files
    files <- list.files(directory, pattern = "\\.txt$", full.names = TRUE)
    apa_list <- lapply(files, read_table)
    names(apa_list) <- basename(files) %>% sub("\\.txt$", "", .)
    
    # Remove REDu and REDi columns if present
    apa_list <- lapply(apa_list, function(df) {
        df %>% select(-any_of(c("REDu", "REDi")))
    })
    
    # Organize by cell type
    cell_types <- c("Astro", "Micro", "Oligo", "OPC", "Exc_upper", "Exc_int", "Exc_deep", "Inh")
    apa_by_celltype <- lapply(cell_types, function(ct) {
        apa_list[grep(ct, names(apa_list))]
    })
    names(apa_by_celltype) <- cell_types
    
    # Save processed data
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(apa_by_celltype, file.path(output_dir, "processed_apa_data.rds"))
    
    return(apa_by_celltype)
}

run_apa_enrichment <- function(apa_data, output_dir = "results/APA/enrichment") {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Function to run gProfiler for one dataset
    run_gprofiler_set <- function(gene_set, name) {
        gostres <- gost(
            gene_set,
            organism = "hsapiens",
            ordered_query = FALSE,
            significant = TRUE,
            exclude_iea = FALSE,
            measure_underrepresentation = FALSE,
            evcodes = TRUE,
            user_threshold = 0.05,
            correction_method = "g_SCS",
            domain_scope = "annotated",
            sources = c("GO:BP", "GO:MF", "GO:CC", "REAC")
        )
        
        if (!is.null(gostres$result)) {
            # Format results for EnrichmentMap
            gem_data <- gostres$result[, c("term_id", "term_name", "p_value", "intersection")]
            colnames(gem_data) <- c("GO.ID", "Description", "p.Val", "Genes")
            gem_data$FDR <- gem_data$p.Val
            gem_data$Phenotype <- "+1"
            
            # Add metadata
            gem_data <- gem_data %>%
                mutate(
                    Disease = if_else(grepl("C9ALS", name), "C9ALS", "sALS"),
                    Modification = case_when(
                        grepl("utr_length", name) ~ "length",
                        grepl("utr_short", name) ~ "short",
                        grepl("intron_APA", name) ~ "intron"
                    ),
                    CellType = str_extract(name, paste(cell_types, collapse = "|")),
                    CellTypeColor = case_when(
                        CellType == "Astro" ~ "#4c72b0",
                        CellType == "Micro" ~ "#8c8c8c",
                        CellType == "Oligo" ~ "#da8bc3",
                        CellType == "OPC" ~ "#937860",
                        CellType == "Exc_upper" ~ "#8172b3",
                        CellType == "Exc_int" ~ "#c44e52",
                        CellType == "Exc_deep" ~ "#55a868",
                        CellType == "Inh" ~ "#dd8452"
                    )
                )
            
            # Save results
            write_tsv(gem_data, 
                     file.path(output_dir, paste0(name, "_enrichment.txt")))
            
            return(gem_data)
        }
        return(NULL)
    }
    
    # Run enrichment for each cell type and condition
    enrichment_results <- list()
    for (cell_type in names(apa_data)) {
        cell_results <- lapply(names(apa_data[[cell_type]]), function(name) {
            run_gprofiler_set(apa_data[[cell_type]][[name]], name)
        })
        names(cell_results) <- names(apa_data[[cell_type]])
        enrichment_results[[cell_type]] <- cell_results
    }
    
    # Save complete results
    saveRDS(enrichment_results, 
            file.path(output_dir, "apa_enrichment_results.rds"))
    
    return(enrichment_results)
}

# Generate summary statistics
generate_apa_summary <- function(apa_data, output_dir = "results/APA") {
    summary_stats <- tibble()
    
    for (cell_type in names(apa_data)) {
        for (condition in names(apa_data[[cell_type]])) {
            data <- apa_data[[cell_type]][[condition]]
            
            stats <- tibble(
                CellType = cell_type,
                Condition = condition,
                Total_Events = nrow(data),
                UTR_Length = sum(grepl("utr_length", rownames(data))),
                UTR_Short = sum(grepl("utr_short", rownames(data))),
                Intron_APA = sum(grepl("intron_APA", rownames(data)))
            )
            
            summary_stats <- bind_rows(summary_stats, stats)
        }
    }
    
    write_tsv(summary_stats, 
              file.path(output_dir, "apa_summary_statistics.txt"))
    
    return(summary_stats)
}

# Main workflow
main <- function() {
    # Process APA data
    apa_data <- process_apa_data()
    
    # Run enrichment analysis
    enrichment_results <- run_apa_enrichment(apa_data)
    
    # Generate summary statistics
    summary_stats <- generate_apa_summary(apa_data)
    
    return(list(
        apa_data = apa_data,
        enrichment = enrichment_results,
        summary = summary_stats
    ))
}

# Run analysis
results <- main()