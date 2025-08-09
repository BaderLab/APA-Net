library(openxlsx)

# Set the directory path
directory_path <- "gProfiler_analysis/gProfiler_APA/new"

for (file in file_names) {
data <- read.table(file, sep = "\t", header = TRUE)
sheet_name <- create_sheet_name(file)
addWorksheet(wb, sheet_name)
writeData(wb, sheet = sheet_name, x = data)
}

wb <- createWorkbook()

for (file in file_names) {
data <- tryCatch({
  read.table(file, sep = "\t", header = TRUE, fill = TRUE, stringsAsFactors = FALSE)
  }, error = function(e) {
  message(paste("Error reading file:", file))
  NULL
})

if (!is.null(data)) {
  sheet_name <- create_sheet_name(file)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = data)
}
}

saveWorkbook(wb, "pathway_results.xlsx", overwrite = TRUE)