remove_file_ids<- function(file_path){
  sample_name <- strsplit(basename(file_path), "\\.")[[1]][1]
  cat("Using sample name:", sample_name, "\n\n")
  
  # Build full paths the same way your Shiny app does
  info <- h5ls(file_path)
  info$full_path <- ifelse(info$group == "/",
                           paste0("/", info$name),
                           paste0(info$group, "/", info$name))
  
  # Print all paths containing "sample_name" so you can verify
  cat("Found sample_name datasets:\n")
  matches <- info[grepl("sample_name", info$name), ]
  print(matches$full_path)
  cat("\n")
  
  paths <- c(
    "/assays/protein_read_counts/ra/sample_name",
    "/assays/protein_read_counts/metadata/sample_name",
    "/assays/dna_read_counts/ra/sample_name",
    "/assays/dna_read_counts/metadata/sample_name",
    "/assays/dna_variants/ra/sample_name",
    "/assays/dna_variants/metadata/sample_name",
    "/all_barcodes/dna_read_counts/metadata/sample_name",
    "/all_barcodes/dna_read_counts/ra/sample_name",
    "/all_barcodes/protein_read_counts/ra/sample_name",
    "/all_barcodes/protein_read_counts/metadata/sample_name"
  )
  
  for (path in paths) {
    if (path %in% info$full_path) {
      existing <- h5read(file_path, path)
      h5write(rep(sample_name, length(existing)), file_path, path)
      cat("Updated:", path, "\n")
    } else {
      cat("NOT FOUND:", path, "\n")
    }
  }
  
  cat("\nDone!\n")
}