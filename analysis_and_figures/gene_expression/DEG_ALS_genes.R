library(Seurat)
library(scCustomize)
library(dplyr)
library(tidyr)
set.seed(1234)

snRNA <- readRDS(file = "MAST_round1/snRNA_final.RDS")
table(Idents(snRNA))
snRNA_glia <- subset(snRNA, idents = c("Oligo", "OPC", "Astro", "Micro"))
snRNA_neuron <- subset(snRNA, idents = c("Exc_upper", "Exc_int", "Exc_deep", "Inh"))

celltype_colors <- c("Oligo" = '#DA8BC3', "OPC" = '#937860', "Astro" = '#4C72B0',
                     "Micro" = '#8C8C8C', "Exc_upper" = '#8172B3',
                     "Exc_int" = '#C44E52', "Exc_deep" = '#55A868',
                     "Inh" = '#DD8452')

ALSgenes<-c("TARDBP","C9orf72","SOD1","FUS","NEK1","OPTN","CHCHD10","SQSTM1",
            "TBK1","KIF5A","SETX","UBQLN2","MATR3","VAPB","SIGMAR1","ANXA11",
            "TUBA4A","ALS2","GRN","PFN1","CHMP2B","TIA1","ANG","SPAST","FIG4",
            "SPG11","GLE1","CCNF","ATXN2","VCP")

snRNA <- readRDS("D:/Documents/OneDrive - University of Toronto/ALS_atlas_resubmission/MAST_round1/snRNA_final.RDS")
 
rbps <- read.table('D:/Documents/OneDrive - University of Toronto/ALS_atlas_resubmission/ALS_atlas/revision_scripts/tomtom_rbp_hits_2.txt', header = F)
colnames(rbps) <- c('rbps')
table(snRNA$celltypes)


CPA_genes <- c("PCF11", "CLP1", 
               "CPSF7", "CPSF6",  "NUDT21",
               "CPSF2", "CPSF3", "CPSF4", "FIP1L1", "CPSF1", "CPSF4L",
               "SYMPK", "CSTF1", "CSTF2", "CSTF2T", "CSTF3","RBBP6",
               "LEO1",  "PAF1",  "CTR9", "CDC73", "PHF3",
               "PABPC1", "PABPC4","PABPN1",
               "PAPOLA", "PAPOLG",
               "CPEB1", "SF3A1", "THOC5")

# Function to filter Tanz_DEGs based on ALSgenes
filter_als_genes <- function(df, genes) {
  filtered_df <- df %>%
    filter(gene %in% genes)
  
  return(filtered_df)
}

# Example usage with your Tanz_DEGs dataframe
filtered_Tanz_DEGs <- filter_als_genes(Tanz_DEGs, ALSgenes)
filtered_Tanz_DEGs <- filter_als_genes(Tanz_DEGs, CPAgenes)

ALSgenes_dotplot<-Clustered_DotPlot(snRNA,features=ALSgenes,x_lab_rotate=TRUE, k=7, seed=1234)
ALSgenes_dotplot<-Clustered_DotPlot(snRNA,features=CPA_genes,x_lab_rotate=TRUE, k=7, seed=1234)
celltype_colors <- c("#DA8BC3", "#937860", "#4C72B0", "#8C8C8C", "#8172B3", "#C44E52", "#55A868", "#DD8452")
RBPgenes_dotplot<-Clustered_DotPlot(snRNA, features=unique(rbps$rbps),colors_use_idents = celltype_colors,x_lab_rotate=FALSE,k=8,seed=1234)
CPAgenes_dotplot<-Clustered_DotPlot(snRNA, features=CPA_genes,colors_use_idents = celltype_colors,x_lab_rotate=FALSE,k=8,seed=1234)

pdf(file = "RBP_expn.pdf", width = 6, height = 15)
print(RBPgenes_dotplot[[2]])
dev.off()

pdf(file = "CPA_expn.pdf", width = 7.5, height = 6)
print(CPAgenes_dotplot[[2]])
dev.off()

Idents(snRNA_glia) <- 'subclass_DE'
table(Idents(snRNA_glia))

Idents(snRNA_neuron) <- 'subclass_DE'
table(Idents(snRNA_neuron))

ALSgenes_glia <-Clustered_DotPlot(snRNA_glia, features=ALSgenes,x_lab_rotate=TRUE, k=7, seed=1234)

ALSgenes_neuron <-Clustered_DotPlot(snRNA_neuron, features=ALSgenes,x_lab_rotate=TRUE, k=7, seed=1234)

print(ALSgenes_neuron)


# Initialize an empty list to store the results
result_list <- list()

# Recursive function to search through nested lists
search_genes <- function(deg_list, ALSgenes, parent_name = NULL) {
  results <- list()
  for (name in names(deg_list)) {
    current_item <- deg_list[[name]]
    if (is.list(current_item)) {
      # Recursive call if the current item is a list
      results <- c(results, search_genes(current_item, ALSgenes, name))
    } else if (is.data.frame(current_item)) {
      # Check for presence of ALSgenes in the data frame
      for (gene in ALSgenes) {
        if (gene %in% current_item$gene) {
          if (is.null(results[[gene]])) {
            results[[gene]] <- list()
          }
          df_name <- if (is.null(parent_name)) name else paste(parent_name, name, sep = "_")
          results[[gene]] <- c(results[[gene]], df_name)
        }
      }
    }
  }
  return(results)
}

# Get the results
results <- search_genes(nested_list, ALSgenes)



# Function to filter Tanz_DEGs based on ALSgenes and select specific columns
create_filtered_dataframe <- function(df, genes) {
  filtered_df <- df %>%
    filter(gene %in% genes) %>%
    select(gene, comparison, celltype, z_score)
  
  return(filtered_df)
}

# Example usage with your Tanz_DEGs dataframe
filtered_Tanz_DEGs <- create_filtered_dataframe(Tanz_DEGs, CPA_genes)

# Remove any rows containing "L23-IT" or "L4-IT" in the celltype column
filtered_Tanz_DEGs <- filtered_Tanz_DEGs %>%
  filter(!celltype %in% c("L23-IT", "L4-IT","In-SNCG", "In-SST"))

# Create a unique identifier by combining comparison and celltype
filtered_Tanz_DEGs <- filtered_Tanz_DEGs %>%
  mutate(comparison_celltype = paste(comparison, celltype, sep = "_"))

# Create the matrix for the heatmap using the combined comparison_celltype and gene
heatmap_matrix <- filtered_Tanz_DEGs %>%
  select(comparison_celltype, gene, z_score) %>%
  spread(key = comparison_celltype, value = z_score)

# Convert the matrix to appropriate format for pheatmap
heatmap_matrix <- as.data.frame(heatmap_matrix)
rownames(heatmap_matrix) <- heatmap_matrix$gene
heatmap_matrix$gene <- NULL

# Create a separate matrix indicating which values were NA
na_matrix <- ifelse(is.na(heatmap_matrix), "NA", "")

# Replace NA values with 0 for clustering (we'll visualize them separately)
heatmap_matrix[is.na(heatmap_matrix)] <- 0

# Define the color palette with white at zero and blue/red at min/max
min_z <- min(heatmap_matrix)
max_z <- max(heatmap_matrix)
midpoint <- 0
colors <- colorRampPalette(c("blue", "white", "red"))(100)

# Create breaks for the color palette
breaks <- c(seq(min_z, midpoint, length.out = 51), 
            seq(midpoint + (max_z - min_z) / 100, max_z, length.out = 50))

# Create an annotation for columns
annotation_col <- data.frame(
  Comparison = sub("_.*", "", colnames(heatmap_matrix))
)
rownames(annotation_col) <- colnames(heatmap_matrix)
annotation_colors <- list(
  Comparison = c(C9ALS = "#009E72", sALS = "#0072B2", C9ALSFTLD = "#FF0000")
)

# Create the heatmap with custom annotation for NA values
pheatmap(heatmap_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         scale = "none", 
         main = "Heatmap of ALS Genes Z-Scores",
         show_rownames = TRUE, 
         show_colnames = TRUE,
         na_col = "lightgray",
         display_numbers = na_matrix,
         number_color = "black",
         fontsize_number = 8,
         color = colors,
         breaks = breaks,
         annotation_col = annotation_col,
         annotation_colors = annotation_colors)