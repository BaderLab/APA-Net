# Function to write gmt file to disk
write_gmt <- function(gmt_list, file_name) {
  conn <- file(file_name, "w")
  for (gene_set_name in names(gmt_list)) {
    gene_set <- gmt_list[[gene_set_name]]
    line <- paste(gene_set_name, gene_set_name, paste(gene_set, collapse = "\t"), sep = "\t")
    writeLines(line, conn)
  }
  close(conn)
}
