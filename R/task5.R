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
