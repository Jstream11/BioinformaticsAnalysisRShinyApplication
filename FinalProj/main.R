library(tidyverse)
library(ggplot2)

#counts functions
counts_filter <- function(counts_df, var_percentile, nonZero_threshold){
  #make sure counts_matrix is matrix
  counts_matrix <- data.matrix(counts_df)
  rownames(counts_matrix)= counts_df[,1] #keep gene names
  counts_matrix <- counts_matrix[,-1]
  
  #filter out genes with all nonzero genes
  non_zero_genes <- counts_matrix[rowSums(counts_matrix>0) >0,]
  
  #calculate variance
  gene_variances <- apply(non_zero_genes, 1, var)
  variance_threshold <- quantile(gene_variances, var_percentile) #calculate variance threshold 
  var_filtered_counts <- non_zero_genes[gene_variances >= variance_threshold,]
  
  #filter for nonzero genes -> removing genes with more zeros than threshold
  final_filtered_counts <- var_filtered_counts[rowSums(var_filtered_counts == 0) <= nonZero_threshold,]
  return(final_filtered_counts)
}

scatter_plot_var <- function(filtered_counts_matrix, counts_df){
  # - Plot median count vs. variance.
  # - Use a log scale for the variance
  # - Mark genes passing filters in a darker color and genes filtered out in a lighter color.
  #making counts df matrix and adding back gene names
  counts_matrix <- data.matrix(counts_df)
  rownames(counts_matrix)= counts_df[,1] #keep gene names
  counts_matrix <- counts_matrix[,-1]
  
  #adding column with filtered status
  filtered_genes <- rownames(counts_matrix)[!(rownames(counts_matrix) %in% rownames(filtered_counts_matrix))]
  counts_matrix <- as.data.frame(counts_matrix)
  counts_matrix$filtered_status <- ifelse(rownames(counts_matrix) %in% filtered_genes, "Filtered Out", "Not Filtered")
  
  #calculating variance and medians
  variance <- apply(counts_matrix[, 1:(ncol(counts_matrix) - 1)], 1, var)
  gene_medians <- apply (counts_matrix[, 1:(ncol(counts_matrix) - 1)], 1, median)
  
  #creating plot
  plot <- ggplot(counts_matrix, aes(x = gene_medians, y = log10(variance), color=filtered_status)) +
    geom_jitter() +
    labs(title = "Gene count median vs. Log10(gene variance)", y = 'log10(variance)', x = 'count median') 
  return(plot)
}

scatter_plot_zeros <- function(filtered_counts_matrix, counts_df){
  # - Mark genes passing filters in a darker color and genes filtered out in a lighter color.
  # - Plot median count vs. the number of zeros.
  
  #making counts df matrix and adding back gene names
  counts_matrix <- data.matrix(counts_df)
  rownames(counts_matrix)= counts_df[,1] #keep gene names
  counts_matrix <- counts_matrix[,-1]
  
  #adding column with filtered status
  filtered_genes <- rownames(counts_matrix)[!(rownames(counts_matrix) %in% rownames(filtered_counts_matrix))]
  counts_matrix <- as.data.frame(counts_matrix)
  counts_matrix$filtered_status <- ifelse(rownames(counts_matrix) %in% filtered_genes, "Filtered Out", "Not Filtered")
  
  #calculating number of zeros/gene and medians
  n_zeros <- rowSums(counts_matrix[, 1:(ncol(counts_matrix) - 1)] == 0)
  gene_medians <- apply (counts_matrix[, 1:(ncol(counts_matrix) - 1)], 1, median)
  
  #creating plot
  plot <- ggplot(counts_matrix, aes(x = gene_medians, y = n_zeros, color=filtered_status)) +
    geom_jitter() +
    labs(title = 'Gene count median vs. Number of Zero Measured Samples in Gene', y = 'number of zeros', x = 'count median') 
  return(plot)
}

clustered_heatmap <- function(filtered_counts_matrix){
  library(gplots)
  filtered_counts_matrix <- log2(filtered_counts_matrix + 1)
  heatmap.2(filtered_counts_matrix,
            xlab = "Sample",
            ylab = "Gene"
  )
}

plot_pca <- function(filtered_counts_matrix, x, y) {
  #Allow the user to select which principal components to plot (e.g., PC1 vs PC2) or plot the top N principal components as a beeswarm plot.
  #Include the % variance explained by each component in the plot labels.
  # prcomp performs PCA
  filtered_counts_df <- as.data.frame(filtered_counts_matrix)
  pca <- prcomp(filtered_counts_df)
  
  # Extract the PC1 and PC2 projections
  pc_x <- pca$x[, x]
  pc_y <- pca$x[, y]
  
  # Create a scatter plot of the PC1 and PC2 projections
  plot <- ggplot(filtered_counts_df, aes(x = pc_x, y = pc_y)) +
    geom_point() +
    labs(x = paste("PC", x), y = paste("PC", y), title = "PCA Projections of Filtered Counts Matrix")
  
  return(plot)
}

#differential expression functions
#plot raw pvals
plot_pvals <- function(de_data) {
  de_data <- as.data.frame(de_data)
  plot <- ggplot(de_data, aes(x = pvalue)) +
    geom_histogram(fill = "skyblue", color = "blue", binwidth=0.02) +
    labs(title = 'Histogram of raw p-values obtained from DE analysis (vP0 vs vAd)')
  
  return(plot)
}

#plot log2fc histogram, need to filter based on padj threshold
plot_log2fc <- function(de_data, padj_threshold) {
  de_data <- as_tibble(de_data)
  de_data_filtered <- de_data %>%
    filter(padj<padj_threshold)
  plot <- ggplot(de_data_filtered, aes(x = log2FoldChange)) +
    geom_histogram(fill = "skyblue", color = "blue", binwidth = 0.2) +
    labs(title = 'Histogram of Log2FoldChange obtained from DE analysis (vP0 vs vAd)')
  return(plot)
}

# volcano plot
plot_volcano <- function(de_data, padj_threshold){
  de_data_tibble <- as_tibble(de_data) %>% mutate(
    volc_plot_status = ifelse(
      padj < padj_threshold & log2FoldChange > 0, "UP",
      ifelse(
        padj < padj_threshold & log2FoldChange < 0, "DOWN",
        "NS")
  ))
  results_no_NA <- drop_na(de_data_tibble)
  
  plot <- ggplot(results_no_NA, aes(x = log2FoldChange, y = -log10(padj), color=volc_plot_status)) +
    geom_point() +
    labs(title = 'Volcano Plot of DESeq2 differential expression results (vP0 vs. vAd)', y = '-log10(padj)', x = 'log2FoldChange')
  return(plot)
}

