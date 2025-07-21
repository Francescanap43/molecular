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
