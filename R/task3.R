task3_filter_unwanted_genes <- function(counts_matrix, gtf) {
  # Load stringr just in case it's needed
  library(stringr)

  # Helper function to extract values from the 9th GTF column
  extract_characteristics <- function(attr_string, key) {
    pattern <- paste0(key, " \"([^\"]+)\"")
    match <- regmatches(attr_string, regexec(pattern, attr_string))
    sapply(match, function(x) if(length(x) > 1) x[2] else NA)
  }

  # Add gene_name if missing
  if (!"gene_name" %in% colnames(gtf)) {
    gtf$gene_name <- extract_characteristics(gtf$V9, "gene_name")
  }

  # Identify genes to exclude
  ribosomal_proteins <- unique(gtf$gene_id[grepl("^RP[SL]", gtf$gene_name)])
  ribosomal_pseudogenes <- unique(gtf$gene_id[grepl("^RP[SL]", gtf$gene_name) & grepl("pseudogene", gtf$gene_biotype)])
  mitochondrial_genes <- unique(gtf$gene_id[grepl("^MT-", gtf$gene_name)])

  # Combine all genes to remove
  all_genes_to_remove <- unique(c(ribosomal_proteins, ribosomal_pseudogenes, mitochondrial_genes))

  # Create summary table
  summary_table <- data.frame(
    Category = c("Ribosomal Proteins", "Ribosomal Pseudogenes", "Mitochondrial Genes"),
    Genes_Removed = c(length(ribosomal_proteins), length(ribosomal_pseudogenes), length(mitochondrial_genes))
  )

  # Filter out unwanted genes
  valid_genes <- !is.na(rownames(counts_matrix))
  filtered_matrix <- counts_matrix[valid_genes & !(rownames(counts_matrix) %in% all_genes_to_remove), ]

  # Return filtered matrix and summary table
  return(list(
    filtered_matrix = filtered_matrix,
    summary_table = summary_table
  ))
}
