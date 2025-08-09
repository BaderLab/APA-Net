# function to read GMT file from disk
read_gmt <- function(file_path) {
  con <- file(file_path, "r")
  gmt_list <- list()
  while(TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) { break }
    split_line <- strsplit(line, "\t")[[1]]
    gmt_list[[split_line[1]]] <- split_line[-c(1, 2)]
  }
  close(con)
  return(gmt_list)
}