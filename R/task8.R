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
