# Step 1: Extract gene expression data
expression_data <- lapply(list(E16.fs, E18.fs, P10.fs, P28.fs), function(obj) {
  # Get the gene expression matrix using GetAssayData
  gene_expr <- as.data.frame(GetAssayData(obj, assay = "RNA")$counts)
  # Reset row names to ensure uniqueness
  rownames(gene_expr) <- seq_len(nrow(gene_expr))
  # Add the Seurat object name as a column
  gene_expr$seurat_object_name <- obj@meta.data$seurat_object_name
  return(gene_expr)
})

# Step 2: Merge data frames by common genes
merged_expression_data <- Reduce(function(x, y) merge(x, y, by = "row.names", all = FALSE), expression_data)

# Check the dimensions of the merged data frame
cat("Merged data frame dimensions:", dim(merged_expression_data), "\n")

# Step 3: Calculate average expression for each gene
average_expression <- lapply(expression_data, function(data) {
  apply(data, 2, mean, na.rm = TRUE)
})

# Step 4: Identify genes that continuously decrease
decreasing_genes <- Reduce(intersect, lapply(average_expression[-1], function(expr, prev_expr) {
  decreasing_genes <- names(expr[expr < prev_expr])
}))

# Step 5: Visualize expression trends if desired
# For example, you can create a boxplot for each gene to observe the expression trend
for (gene in decreasing_genes) {
  gene_expression <- lapply(expression_data, function(data) data[, gene, drop = FALSE])
  gene_data <- do.call(rbind, gene_expression)
  gene_data$Age <- factor(rownames(gene_data), levels = c("E16.fs", "E18.fs", "P10.fs", "P28.fs"))
  
  ggplot(gene_data, aes(x = Age, y = !!sym(gene), fill = Age)) +
    geom_boxplot() +
    labs(title = paste("Expression Trend of", gene)) +
    theme_minimal()
}
