#molecular
I am submitting this as my project for the programming exam, this is the read me file

## packages needed

List any required packages:

```r
install.packages(c("ggplot2", "stringr", "SingleCellExperiment", "DropletUtils", "scater","scran","igraph", "org.Hs.eg.db", "celldex", "SingleR", "BiocManager"))
```

#main script
## setting the working directory

## extracting and loading the datasets provvided
gtf_path <-"C:/Users/Utente/Desktop/magistrale/programming/exam/molecular/Homo_sapiens.GRCh38.111.gtf.gz"
gtf <- read.delim(gzfile(gtf_path), header = FALSE, comment.char = "#")

unzip("C:/Users/Utente/Desktop/magistrale/programming/exam/molecular/filtered_feature_bc_matrix.zip",exdir = "C:/Users/Utente/Desktop/magistrale/programming/exam/molecular/filtered_feature_bc_matrix")

## installing DropletUtils to extract the features
if (!requireNamespace("DropletUtils", quietly = TRUE)) {
    BiocManager::install("DropletUtils")
}
library(DropletUtils)

## setting the path to matrix folder to read the features
matrix_dir <-"C:/Users/Utente/Desktop/magistrale/programming/exam/molecular/filtered_feature_bc_matrix/filtered_feature_bc_matrix"
sce <- read10xCounts(matrix_dir)

## checking everything

sce

# 1.Gene Annotation: Identify and retain only the protein-coding genesfrom the dataset, based on the GTF file.

## creating a table for gtf
gtf <-read.table(gzfile("C:/Users/Utente/Desktop/magistrale/programming/exam/molecular/Homo_sapiens.GRCh38.111.gtf.gz"),
header = FALSE, sep = "\t", comment.char = "#", stringsAsFactors =
FALSE, quote = "")

head(gtf)

## creating a function to estract the characteristics of the genes fromthe 9th column
library(stringr)

extract_characteristics <- function(attr_string, key) {
  pattern <- paste0(key, " \"([^\"]+)\"")
  match <- regmatches(attr_string, regexec(pattern, attr_string))
  sapply(match, function(x) if(length(x) > 1) x[2] else NA)
}

## extracting id and biotype in their column
gtf$gene_id <- extract_characteristics(gtf$V9, "gene_id")
gtf$gene_biotype <- extract_characteristics(gtf$V9,"gene_biotype")

## finding the protein coding genes making it searchig in the gene biotype column
protein_coding_genes <- unique(gtf$gene_id[gtf$gene_biotype == "protein_coding"])

## creating a matrix only with the protein coding genes
filtered_matrix <- sce[rownames(sce) %in% protein_coding_genes, ]

## cheking it worked
ls()
class(sce)
dim(sce)
dim(filtered_matrix)

# 2.Gene Expression Summary:For each cell, calculate the number of genes with expression ≥3 UMIs. Display the distribution of these counts using a violin plot.

## calculating the counts to create a sparse matrix
library(SingleCellExperiment)
counts_matrix <- counts(sce)

## selecting only the genes with expression ≥3 UMIs
genes_over_3 <-colSums(counts_matrix >= 3)

## checking it worked
summary(genes_over_3)

## getting the ggplot2 library to create the plot
library(ggplot2)

## converting to data frame for plotting
df <- data.frame(GenesOver3UMI = genes_over_3)

## violin plot
ggplot(df, aes(x = "", y = GenesOver3UMI)) +
geom_violin(fill = "blue", alpha = 0.6) + geom_boxplot(width = 0.1,
outlier.shape = NA) + labs( title = "Number of Genes with ≥3 UMIs", y =
"Genes with ≥3 UMIs", x = "" ) + theme_minimal(12)

# 3.Gene Filtering:Exclude the following gene categories:Ribosomal proteins, Ribosomal pseudogenes, Mitochondrial genes. Provide a summary table listing the number of genes removed from each category.

## creting a gene name column to identify the family of the genes
gtf$gene_name <- extract_characteristics(gtf$V9, "gene_name")

## finding ribosomal proteins: gene names starting with "RPS" or "RPL"
ribosomal_proteins <-unique(gtf$gene_id[grepl("^RP[SL]", gtf$gene_name)])

## finding ribosomal pseudogenes: gene_type contains "pseudogene" and name starts with RPS or RPL
ribosomal_pseudogenes <- unique(gtf$gene_id[grepl("^RP[SL]", gtf$gene_name) & grepl("pseudogene",gtf$gene_biotype) ])

## finding mitochondrial genes: gene names starting with "MT-"
mitochondrial_genes <- unique(gtf$gene_id[grepl("^MT-", gtf$gene_name)])

## creating the summary table with all the genes found
all_genes_to_remove <- unique(c(ribosomal_proteins, ribosomal_pseudogenes, mitochondrial_genes))
summary_table <- data.frame( Category =c("Ribosomal Proteins", "Ribosomal Pseudogenes", "Mitochondrial Genes"),Genes_Removed = c(length(ribosomal_proteins),length(ribosomal_pseudogenes), length(mitochondrial_genes)) )

print(summary_table)

## filtereing them out
valid_genes <- !is.na(rownames(counts_matrix))
filtered_matrix <- counts_matrix[ valid_genes &!(rownames(counts_matrix) %in% all_genes_to_remove),]

## cheking
dim(sce)
dim(filtered_matrix)

# 4.Principal Component Analysis (PCA):Assess the variance explained bythe first 20 principal components and visualize this using a histogram.

install.packages("BiocManager")
BiocManager::install(c("DropletUtils","scater","scran","SingleCellExperiment"))
a
library(DropletUtils)
library(scater)
library(SingleCellExperiment)


## normalising the data
sce <- logNormCounts(sce)

## running PCA
sce <- runPCA(sce)

## calculating the variance explained by each PC
var_explained <-attr(reducedDim(sce, "PCA"), "percentVar")

## ploting histogram of the first 20 PCs
barplot(var_explained[1:20],
names.arg = 1:20, xlab = "Principal Component", ylab = "Percentage of
Variance Explained", main = "Variance Explained by First 20 PCs", col =
"red")

# 5.UMAP Visualization:Generate a UMAP using a subset of principal components you think are sufficiently informative. Justify your choice of components in the report.

## visualising the variance to choose the subset
var_explained

## deciding on the 6th componment since it is the last higher than 1
## running the UMAP
BiocManager::install(c("scater", "scran"))
a
library(scater)
sce <- runUMAP(sce, dimred = "PCA", n_dimred = 6)

library(ggplot2)
umap_coords <- as.data.frame(reducedDim(sce, "UMAP"))
colnames(umap_coords) <- c("UMAP1", "UMAP2")

## plotting the UMAP
ggplot(umap_coords, aes(x = UMAP1, y = UMAP2)) +
geom_point(size = 0.3, alpha = 1) + labs(title = "UMAP Based on First 6
Principal Components") + theme_minimal(12)

# 6.Clustering:Apply a clustering algorithm of your choice and visualize the clusters. Clearly describe:The clustering method used, Parameters and resolution settings, Interpretation of the resulting clusters.Building shared nearest neighbor graph

library(SingleCellExperiment)
library(scran)
library(scater)
library(igraph)
snn_graph <-buildSNNGraph(sce, use.dimred = "PCA", k = 10)

## clustering with the Walktrap algorithm
clusters <- igraph::cluster_walktrap(snn_graph)$membership

## storing clusters in sce
colLabels(sce) <- factor(clusters)

umap_df <- as.data.frame(reducedDim(sce, "UMAP"))
umap_df$Cluster <- colLabels(sce)

## plotting the UMAP based on the clusters found
library(ggplot2)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
geom_point(size = 0.3, alpha = 1) + labs(title = "UMAP Colored by
Clusters (Walktrap, k=10)") + theme_minimal(10)

# 7.Cell Type Annotation:Using a cell annotation tool of your choice annotate the dataset and discuss:Concordance or discordance between clusters and predicted cell types,Evidence of heterogeneity or homogeneity within clusters

## using sindleR
BiocManager::install("org.Hs.eg.db")
a
## selecting all
library(org.Hs.eg.db)
library(SingleCellExperiment)

## loading reference dataset
ensembl_ids <- rownames(sce)

## removing version numbers
ensembl_ids <- sub("\\..*", "", ensembl_ids)

## mapping Ensembl to symbols
gene_symbols <- mapIds( org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first" )

## updating rownames
rownames(sce) <- gene_symbols

## filtering out genes with NA
sce <- sce[!is.na(rownames(sce)), ]

## running SingleR
library(celldex)
library(SingleR)
ref <- celldex::HumanPrimaryCellAtlasData()
pred <- SingleR(test = sce, ref = ref, labels = ref$label.main)
sce$predicted_cell_type <- pred$labels

## visualising the difference
library(ggplot2)
umap_df <-as.data.frame(reducedDim(sce, "UMAP"))
umap_df$CellType <- sce$predicted_cell_type

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = CellType)) +
geom_point(size = 0.3, alpha = 1) + labs(title = "UMAP Colored by
Predicted Cell Types") + theme_minimal(10)

### the 2 plots enlight the fact thata the clusters are many more than just the Human Primary Cell Atlas Data so there are other clusters that were not considereded in the data set I used, in fact there are many other types of cells beyond the immune system ones.

# 8.Tissue Origin Inference:Based on gene markers, annotation, and cluster structure, propose a hypothesis about the tissue of origin for the dataset. Justify your guess with biological reasoning.

## ensuring both columns are clean
markers_genes <- toupper(trimws(gtf$gene_symbol))
panglao_genes <- toupper(trimws(panglao$official.gene.symbol))

## matching by intersection
matched_genes <- intersect(markers_genes, panglao_genes)

## filtering PanglaoDB by matched genes
matched_table <- panglao[toupper(panglao$official.gene.symbol) %in% matched_genes, ]

## counting how many matched each cell type
celltype_counts <- as.data.frame(table(matched_table$cell.type))
celltype_counts <- celltype_counts[order(-celltype_counts$Freq), ]
print(celltype_counts)

## plotting histogram
library(ggplot2)
colnames(celltype_counts) <- c("CellType", "Count")

ggplot(celltype_counts, aes(x = reorder(CellType, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "yellow") +
  theme_minimal() +
  labs(title = "Cell Type Frequency (from PanglaoDB markers)",
       x = "Cell Type",
       y = "Number of Marker Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## they are too much to visualise, let's see the top 25 tissues
library(ggplot2)

## selecting top 25 most frequent cell types
top25 <-head(celltype_counts[order(-celltype_counts$Count), ], 25)

## plotting histogram
ggplot(top25, aes(x = reorder(CellType, -Count), y =
Count)) + geom_bar(stat = "identity", fill = "yellow") +
theme_minimal() + labs(title = "Top 25 Cell Types by Marker Frequency",
x = "Cell Type", y = "Number of Marker Genes") + theme(axis.text.x =
element_text(angle = 45, hjust = 1))

### As I was saying in the previus point, here we can appresciate how the most present cell types are neuronal cells, entero cells, fibroblasts, hepatocytes, there are a higher specific numer of genes for neuronal cells, so it is natural to see a higher numer of gene symbles in a general data set

