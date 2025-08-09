library(readxl)
library(purrr)

# Function to read all sheets from a file
read_all_sheets <- function(file_path) {
  sheets <- excel_sheets(file_path)
  sheets_data <- lapply(sheets, function(sheet) {
    read_excel(file_path, sheet = sheet)
  })
  names(sheets_data) <- sheets
  return(sheets_data)
}
