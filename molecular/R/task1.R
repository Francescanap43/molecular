#' task 1: loading the data using DropletUtils, retaining only the protein coding genes in #' the data
#'
task1_load_and_filter_data <- function(gtf_path, matrix_zip_path, matrix_extract_dir) {
  # setting the working directory
  # extracting and loading the datasets provvided
  gtf_path <-"C:/Users/Utente/Desktop/magistrale/programming/exam/Homo_sapiens.GRCh38.111.gt   f.gz"
  gtf <- read.delim(gzfile(gtf_path), header = FALSE, comment.char = "#")

  unzip("C:/Users/Utente/Desktop/magistrale/programming/exam/filtered_feature_bc_matrix.zip"   ,exdir = "C:/Users/Utente/Desktop/magistrale/programming/exam/filtered_feature_bc_matrix")

  # Load required packages
  if (!requireNamespace("DropletUtils", quietly = TRUE)) {
    BiocManager::install("DropletUtils")
  }
  library(DropletUtils)
  library(stringr)

  # Read GTF file
  gtf <- read.table(gzfile(gtf_path), header = FALSE, sep = "\t",
                    comment.char = "#", stringsAsFactors = FALSE, quote = "")

  # Extract gene_id and gene_biotype from column 9
  extract_characteristics <- function(attr_string, key) {
    pattern <- paste0(key, " \"([^\"]+)\"")
    match <- regmatches(attr_string, regexec(pattern, attr_string))
    sapply(match, function(x) if(length(x) > 1) x[2] else NA)
  }

  gtf$gene_id <- extract_characteristics(gtf$V9, "gene_id")
  gtf$gene_biotype <- extract_characteristics(gtf$V9, "gene_biotype")

  # Identify protein-coding genes
  protein_coding_genes <- unique(gtf$gene_id[gtf$gene_biotype == "protein_coding"])

  # Unzip the matrix if not already extracted
  if (!dir.exists(matrix_extract_dir)) {
    unzip(matrix_zip_path, exdir = dirname(matrix_extract_dir))
  }

  # Read 10x matrix
  sce <- read10xCounts(matrix_extract_dir)

  # Filter matrix to keep only protein-coding genes
  filtered_matrix <- sce[rownames(sce) %in% protein_coding_genes, ]

  # Return both original and filtered matrices for further use
  return(list(
    sce = sce,
    filtered_matrix = filtered_matrix
  ))
}
