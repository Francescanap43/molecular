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
