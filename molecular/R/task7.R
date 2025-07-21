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
