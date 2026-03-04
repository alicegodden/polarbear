# Set your folder path
folder_path <- "~/Desktop/IMMLER_PAPERS/polarbear/newgenome/denovo_assembly_ASMv1lib_repmask/divs"

# Get all divsum files in folder
files <- list.files(folder_path, pattern = "\\.divsum$", full.names = TRUE)

for (f in files) {
  # Read file as lines
  lines <- readLines(f)
  
  # Find the first "Div" line (allow leading spaces)
  start_idx <- grep("^\\s*Div", lines)[1]
  
  if (is.na(start_idx)) {
    message("No 'Div' line found in: ", f)
    next
  }
  
  # Keep from that line onward
  lines <- lines[start_idx:length(lines)]
  
  # Read into data frame
  df <- read.table(text = lines, header = TRUE, sep = "", stringsAsFactors = FALSE)
  
  # Get the original header line from file and normalize spacing to one space
  header_line <- gsub("\\s+", " ", lines[1])
  
  # Convert table rows to strings with single space separation
  data_lines <- apply(df, 1, function(row) paste(row, collapse = " "))
  
  # Combine header and data lines
  output_lines <- c(header_line, data_lines)
  
  # Write output without column names or quotes
  outname <- sub("\\.divsum$", ".csv", f)
  writeLines(output_lines, outname)
  
  message("Processed: ", f)
}
