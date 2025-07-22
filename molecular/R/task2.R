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
