# molecular
#' task 1: loading the data using DropletUtils, retaining only the protein coding genes in #' the data
#'
#' This function extracts the data provvided and load them in the enviroment
#' it has a function to extract from gtf the geneid, gene biotype
#' to search for the protein coding genes. Filtering the matrix
#'
#' @param gtf_path
#' @param matrix_dir to extract the 10X matrix
#' @return a filtered matrix only with the protein coding genes
#' @examples
#' \dontrun{
#' gtf_path = "C:/Users/Utente/Desktop/magistrale/programming/exam/Homo_sapiens.GRCh38.111.g#' tf.gz",
#' matrix_zip_path = "C:/Users/Utente/Desktop/magistrale/programming/exam/filtered_feature_b#' c_matrix.zip",
#' matrix_extract_dir = "C:/Users/Utente/Desktop/magistrale/programming/exam/filtered_feat #' ure_bc_matrix/filtered_feature_bc_matrix"
#')
#'sce <- result$sce
#'filtered_matrix <- result$filtered_matrix
#'dim(sce)
#'dim(filtered_matrix)
#' }
#' @export
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


#' Title: distribution of genes with ≥3 UMIs
#'
#' here there is a function tant finds the genes with a ≥3 UMIs and displays it in a violin #' plot
#'
#' @param sce_object loadeding SingleCellExperiment
#' @return violin plot describing the distribution
#' @examples
#' \dontrun{
#' task2_plot_genes_over_3_UMIs(filtered_matrix)
#' }
#' my_function(arg1, arg2)
#' @export
task2_plot_genes_over_3_UMIs <- function(sce_object) {
   # Load required libraries
  library(SingleCellExperiment)
  library(ggplot2)
  
  # Get counts matrix from the SCE object
  counts_matrix <- counts(sce_object)
  
  # Calculate how many genes have ≥3 UMIs for each cell (i.e., column)
  genes_over_3 <- colSums(counts_matrix >= 3)
  
  # Create a data frame for plotting
  df <- data.frame(GenesOver3UMI = genes_over_3)
  
  # Create the violin plot
  plot <- ggplot(df, aes(x = "", y = GenesOver3UMI)) +
    geom_violin(fill = "blue", alpha = 0.6) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    labs(
      title = "Number of Genes with ≥3 UMIs",
      y = "Genes with ≥3 UMIs",
      x = ""
    ) +
    theme_minimal(base_size = 12)
  
  # Show plot and return summary values
  print(plot)
  return(summary(genes_over_3))
}

#' filtering the gtf genes for ribosomial protein, ribosomial pseudogenes and mitoconrial 
#' genes
#'
#' searching for the ribosomial protein, ribosomial pseudogenes and mitoconrial 
#' genes biotype, removing them
#'
#' @param counts_matrix
#' @param gtf
#' @return a summary table of the eliminated genes
#' @examples
#' \dontrun{
#' filtered_matrix_step3 <- result$filtered_matrix
#' print(result$summary_table)
#' }
#' @export
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


#' running a PCA 
#'
#' This function runs a PCA and plots it
#'
#' @param sce_object 
#' @return the plot of the PCA 
#' @examples
#' \dontrun{
#' result <- task4_run_pca_and_plot(filtered_matrix_step3)
#' sce <- result$sce
#' variance_explained <- result$variance_explained
#' }
#' @export
task4_run_pca_and_plot <- function(sce_object) {
  # Load necessary libraries
  library(scater)
  library(SingleCellExperiment)
  
  # Normalize the data
  sce_object <- logNormCounts(sce_object)
  
  # Run PCA
  sce_object <- runPCA(sce_object)
  
  # Extract variance explained
  var_explained <- attr(reducedDim(sce_object, "PCA"), "percentVar")
  
  # Plot variance explained by the first 20 principal components
  barplot(
    var_explained[1:20],
    names.arg = 1:20,
    xlab = "Principal Component",
    ylab = "Percentage of Variance Explained",
    main = "Variance Explained by First 20 PCs",
    col = "red"
  )
  
  # Return the updated SCE object and the variance explained
  return(list(
    sce = sce_object,
    variance_explained = var_explained
  ))
}


#' running a UMAP based on the PCA results 
#'
#' This function runs a UMAP and plots it
#'
#' @param sce_object 
#' @param n_pcs = 6
#' @return the plot of the UMAP
#' @examples
#' \dontrun{
#' result <- task5_run_umap_and_plot(sce, n_pcs = 6)
#' sce <- result$sce
#' umap_coords <- result$umap_coords
#' }
#' @export
task5_run_umap_and_plot <- function(sce_object, n_pcs = 6) {
 # Load required libraries
  library(scater)
  library(ggplot2)
  
  # Run UMAP on top of PCA
  sce_object <- runUMAP(sce_object, dimred = "PCA", n_dimred = n_pcs)
  
  # Extract UMAP coordinates
  umap_coords <- as.data.frame(reducedDim(sce_object, "UMAP"))
  colnames(umap_coords) <- c("UMAP1", "UMAP2")
  
  # Create UMAP plot
  plot <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2)) +
    geom_point(size = 0.3, alpha = 1) +
    labs(
      title = paste("UMAP Based on First", n_pcs, "Principal Components"),
      x = "UMAP1",
      y = "UMAP2"
    ) +
    theme_minimal(base_size = 12)
  
  print(plot)
  
  # Return updated SCE object and UMAP coordinates
  return(list(
    sce = sce_object,
    umap_coords = umap_coords
  ))
}

#' running a UMAP based on clusters
#'
#' This function runs a UMAP on clusters and plots them
#'
#' @param sce_object 
#' @param k=10
#' @return the plot of the UMAP
#' @examples
#' \dontrun{
#' result <- task6_cluster_cells_and_plot(sce, k = 10)
#' sce <- result$sce
#' clusters <- result$clusters
#' }
#' @export
task6_cluster_cells_and_plot <- function(sce_object, k = 10) {
  # Load required libraries
  library(scran)
  library(scater)
  library(SingleCellExperiment)
  library(igraph)
  library(ggplot2)
  
  # Build shared nearest neighbor graph
  snn_graph <- buildSNNGraph(sce_object, use.dimred = "PCA", k = k)
  
  # Perform clustering using Walktrap
  clusters <- igraph::cluster_walktrap(snn_graph)$membership
  
  # Store clusters in SCE
  colLabels(sce_object) <- factor(clusters)
  
  # Extract UMAP coordinates and add cluster info
  umap_df <- as.data.frame(reducedDim(sce_object, "UMAP"))
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$Cluster <- colLabels(sce_object)
  
  # Plot UMAP colored by cluster
  plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
    geom_point(size = 0.3, alpha = 1) +
    labs(title = paste("UMAP Colored by Clusters (Walktrap, k =", k, ")")) +
    theme_minimal(base_size = 10)
  
  print(plot)
  
  # Return updated SCE object with clustering
  return(list(
    sce = sce_object,
    clusters = colLabels(sce_object)
  ))
}

#' running a UMAP based on predicted clusters
#'
#' This function runs a UMAP on predicted clusters and plots them
#'
#' @param sce_object 
#' @return the plot of the UMAP
#' @examples
#' \dontrun{
#' result <- task7_annotate_cells_with_singleR(sce)
#' sce <- result$sce
#' predictions <- result$predictions
#' }
#' @export
task7_annotate_cells_with_singleR <- function(sce_object) {
 # Load required libraries
  library(org.Hs.eg.db)
  library(SingleCellExperiment)
  library(SingleR)
  library(celldex)
  library(AnnotationDbi)
  library(ggplot2)
  
  # Clean Ensembl IDs: remove version numbers
  ensembl_ids <- rownames(sce_object)
  ensembl_ids <- sub("\\..*", "", ensembl_ids)
  
  # Map to gene symbols
  gene_symbols <- mapIds(
    org.Hs.eg.db,
    keys = ensembl_ids,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  # Update rownames with gene symbols
  rownames(sce_object) <- gene_symbols
  
  # Filter out genes with NA rownames
  sce_object <- sce_object[!is.na(rownames(sce_object)), ]
  
  # Load reference data and run SingleR
  ref <- celldex::HumanPrimaryCellAtlasData()
  pred <- SingleR(test = sce_object, ref = ref, labels = ref$label.main)
  
  # Add predicted labels to metadata
  sce_object$predicted_cell_type <- pred$labels
  
  # UMAP plot colored by predicted cell types
  umap_df <- as.data.frame(reducedDim(sce_object, "UMAP"))
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$CellType <- sce_object$predicted_cell_type
  
  plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = CellType)) +
    geom_point(size = 0.3, alpha = 1) +
    labs(title = "UMAP Colored by Predicted Cell Types") +
    theme_minimal(base_size = 10)
  
  print(plot)
  
  # Return updated SCE and predictions
  return(list(
    sce = sce_object,
    predictions = pred
  ))
}

#' confronting the gene symbols of the gtf with a real tissue data base 
#'
#' This function connects the gene symbols with the ones in the Panglao database
#' plotting all the tissues that found a match between the data sets and a plot 
#' that shows only the top 25 tisuues most found
#'
#' @param gtf_data 
#' @param panglao_db 
#' @return the plots of the tissues
#' @examples
#' \dontrun{
#' result <- task8_infer_tissue_origin(gtf, panglao)
#' }
#' @export
task8_infer_tissue_origin <- function(gtf_data, panglao_db) {
  library(ggplot2)
  
  # Ensure gene symbols are uppercase and trimmed
  markers_genes <- toupper(trimws(gtf_data$gene_name))
  panglao_genes <- toupper(trimws(panglao_db$official.gene.symbol))
  
  # Find intersecting genes
  matched_genes <- intersect(markers_genes, panglao_genes)
  
  # Filter PanglaoDB for matched genes
  matched_table <- panglao_db[toupper(panglao_db$official.gene.symbol) %in% matched_genes, ]
  
  # Count frequency of cell types
  celltype_counts <- as.data.frame(table(matched_table$cell.type))
  colnames(celltype_counts) <- c("CellType", "Count")
  celltype_counts <- celltype_counts[order(-celltype_counts$Count), ]
  
  print(celltype_counts)
  
  # Plot all cell types frequency
  p_all <- ggplot(celltype_counts, aes(x = reorder(CellType, -Count), y = Count)) +
    geom_bar(stat = "identity", fill = "yellow") +
    theme_minimal() +
    labs(title = "Cell Type Frequency (from PanglaoDB markers)",
         x = "Cell Type",
         y = "Number of Marker Genes") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p_all)
  
  # Plot top 25 cell types frequency
  top25 <- head(celltype_counts, 25)
  
  p_top25 <- ggplot(top25, aes(x = reorder(CellType, -Count), y = Count)) +
    geom_bar(stat = "identity", fill = "yellow") +
    theme_minimal() +
    labs(title = "Top 25 Cell Types by Marker Frequency",
         x = "Cell Type",
         y = "Number of Marker Genes") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p_top25)
  
  # Return list for further use
  return(list(
    matched_table = matched_table,
    celltype_counts = celltype_counts,
    top25 = top25
  ))
}
