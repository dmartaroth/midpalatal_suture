# ## ######################################## ## #
#                     FUNCTIONS                  #
# ## ######################################## ## #

# Date: Sun Mar 31 15:34:16 2024 ------------------
# For Visium analysis

# Convenience function to save the current plot
convenient_save_plot <- function(plot, name, number = plot_number, width = 10, height = 10, dir = results_folder, file_format = "png") {
  # Increment plot number
  number <- number + 1
  
  # Save the plot based on file format
  if (tolower(file_format) == "pdf") {
    ggsave(filename = file.path(dir, sprintf("%02d_%s.pdf", number, name)), 
           width = width, height = height, plot)
  } else {
    ggsave(filename = file.path(dir, sprintf("%02d_%s.png", number, name)), 
           width = width, height = height, plot)
  }
  
  # Update plot_number in the global environment
  assign("plot_number", number, envir = .GlobalEnv)
  
  return(NULL)  # Return NULL since we're not returning a plot object
}




# Define a function to generate SpatialFeaturePlots with iterative gene selection
generate_spatial_feature_plots <- function(seurat_obj = WT1R1ps, cluster_number, top_genes_file, top_genes_number, increment, dir) {
  # Read the top_genes data from the CSV file
  top_genes_df <- read.csv(top_genes_file)
  
  # Check if cluster_number exists in the data frame
  if (!cluster_number %in% unique(top_genes_df$cluster)) {
    stop("Invalid cluster number. Please provide a valid cluster number.")
  }
  
  # Filter the top genes data frame for the specified cluster number
  cluster_genes <- top_genes_df[top_genes_df$cluster == cluster_number, "gene"]
  
  # Check if any genes were found for the cluster
  if (length(cluster_genes) == 0) {
    stop("No genes found for the specified cluster number.")
  }
  
  # Extract top gene names based on top_genes_number
  gene_names <- head(cluster_genes, top_genes_number)
  
  # Initialize a list to store plots for each increment
  all_plots <- list()
  
  # Loop through gene names and generate SpatialFeaturePlots for each increment
  for (i in seq(1, length(gene_names), increment)) {
    # Get the gene set for the current increment
    gene_set <- gene_names[i:min(i + increment - 1, length(gene_names))]
    
    # Generate individual plots for the gene set
    plots_list <- list()
    for (gene_name in gene_set) {
      # Generate the SpatialFeaturePlot for the current gene from the Seurat object
      plot <- SpatialFeaturePlot(seurat_obj, features = gene_name, pt.size.factor = 6, stroke = NA) +
        theme(aspect.ratio = 0.4, legend.position = "right", legend.text = element_text(size = 8),
              legend.title = element_text(size = 7, face = "italic"), legend.justification = "center",
              axis.text.y = element_text(face = "italic", size = 4, angle = 90),
              axis.text.x = element_text(face = "italic", size = 4, angle = 90)) +
        ggplot2::scale_fill_gradientn(colors = c("white", "#FFE1FF", "#DDA0DD", "#CD69C9", "#9932CC"))
      
      # Store the plot in the list
      plots_list[[gene_name]] <- plot
    }
    
    # Combine individual plots into a grid for the current increment
    combined_plot <- cowplot::plot_grid(plotlist = plots_list, ncol = 2)
    
    # Append the combined plot to the list of all plots
    all_plots[[paste0("topmmgenes_Combined_SpatialFeaturePlots_", cluster_number, "_", i, "_to_", min(i + increment - 1, length(gene_names)))]] <- combined_plot
  }
  
  # Save all combined plots into a single PDF
  for (plot_name in names(all_plots)) {
    convenient_save_plot(
      all_plots[[plot_name]],
      name = plot_name,
      dir = dir,
      height = 10,
      width = 16,
      file_format = "pdf"
    )
  }
}

# Example usage of the function
# cluster_number <- 3  # specify the cluster number you want to visualize
# top_genes_number <- 40  # specify the number of top genes to consider
# increment <- 4
# top_genes_file <- here::here(output, "WT1R1_top50_markers_pct.csv")  # replace with the actual path to your CSV file
# dir <- here::here(results_folder)  # replace with your actual directory path
# 
# # Call the function with the Seurat object 'WT1R1ps' as the argument
# generate_spatial_feature_plots(seurat_obj = WT1R1ps, cluster_number = cluster_number, 
#                                top_genes_file = top_genes_file, top_genes_number = top_genes_number, 
#                                increment = increment, dir = dir)




# Define a function to generate individual violin plots for each gene with descriptive filenames
generate_individual_violin_plots <- function(seurat_obj_list, gene_list, output_dir, suture, file_format = "pdf") {
  # Loop through each gene in the gene list
  for (gene_name in gene_list) {
    # Initialize data frame to store combined data for all ages
    combined_data <- data.frame()
    
    # Loop through Seurat objects and combine data for the current gene across ages
    for (seurat_obj_name in names(seurat_obj_list)) {
      seurat_obj <- seurat_obj_list[[seurat_obj_name]]
      
      # Fetch data for the current gene from the Seurat object
      data_gene <- FetchData(seurat_obj, vars = c(gene_name, "age"), layer = "data")
      
      # Filter data to include only the specified gene and age information
      data_gene <- data_gene[, c(gene_name, "age")]
      
      # Rename columns for clarity
      colnames(data_gene) <- c("Expression", "Age")
      
      # Add dataset information to the data
      data_gene$Dataset <- seurat_obj_name
      
      # Combine data for the current gene with the combined data frame
      combined_data <- rbind(combined_data, data_gene)
    }
    
    # Generate the single violin plot for the current gene across ages
    plot <- ggplot(combined_data, aes(x = Age, y = Expression, fill = Age)) +
      geom_violin(scale = "width", width = 0.7, alpha = 0.8) +
      geom_jitter(width = 0.2, alpha = 0.6) +
      scale_fill_manual(values = my_colors) +
      labs(x = "Age", y = "Expression", fill = "Age") +
      theme_minimal() +
      theme(legend.position = "top")
    
    # Save the single violin plot with a descriptive filename
    filename <- paste0(gene_name, "_violin_plot_", suture, ".", file_format)
    ggsave(filename = filename, plot = plot, path = output_dir, device = file_format)
  }
}

# # Example usage of the function
# seurat_obj_list <- list(E16_frontal_suture = E16.fs,
#                         E18_frontal_suture = E18.fs,
#                         P10_frontal_suture = P10.fs,
#                         P28_frontal_suture = P28.fs)
# gene_list <- c("Sfrp1", "Tnn", "Sfrp2")  # Example gene list
# output_dir <- mm_figs  # Specify the output directory
# filename_part <- "individual"  # Define part of the filename
# my_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")  # Example colors
# 
# # Call the function with the list of Seurat objects, gene list, output directory, filename part, and other parameters
# generate_individual_violin_plots(seurat_obj_list, gene_list, output_dir, filename_part)


# Load necessary packages
library(ggplot2)
library(ggpubr)
library(ggrepel)

# Define a function to generate individual violin plots for top genes from a file
generate_individual_violin_plots_from_file <- function(seurat_obj_list, top_genes_file, output_dir, suture, cluster_number, top_genes_number, file_format = "pdf") {
  # Read the top_genes data from the CSV file
  top_genes_df <- read.csv(top_genes_file)
  
  # Check if 'suture' is specified as part of the filename
  if (is.null(suture) || suture == "") {
    stop("Please specify part of the filename using the 'suture' argument.")
  }
  
  # Filter the top genes data frame for the specified cluster number
  cluster_genes <- top_genes_df[top_genes_df$cluster == cluster_number, "gene"]
  
  # Check if any genes were found for the cluster
  if (length(cluster_genes) == 0) {
    stop("No genes found for the specified cluster number.")
  }
  
  # Extract top gene names based on top_genes_number
  top_genes <- head(cluster_genes, top_genes_number)
  
  # Loop through each gene in the top genes list
  for (gene_name in top_genes) {
    # Initialize data frame to store combined data for all ages
    combined_data <- data.frame()
    
    # Loop through Seurat objects and combine data for the current gene across ages
    for (seurat_obj_name in names(seurat_obj_list)) {
      seurat_obj <- seurat_obj_list[[seurat_obj_name]]
      
      # Use try-catch to handle cases where the gene is not found in the data
      tryCatch({
        # Fetch data for the current gene from the Seurat object
        data_gene <- FetchData(seurat_obj, vars = c(gene_name, "age"), layer = "data")
        
        # Filter data to include only the specified gene and age information
        data_gene <- data_gene[, c(gene_name, "age")]
        
        # Rename columns for clarity
        colnames(data_gene) <- c("Expression", "Age")
        
        # Add dataset information to the data
        data_gene$Dataset <- seurat_obj_name
        
        # Combine data for the current gene with the combined data frame
        combined_data <- rbind(combined_data, data_gene)
      }, error = function(e) {
        # Print a message if the gene is not found
        message(paste("Gene", gene_name, "not found in", seurat_obj_name))
      })
    }
    
    # Check if combined_data is empty (gene not found in any dataset)
    if (nrow(combined_data) == 0) {
      message(paste("Skipping gene", gene_name, "as it was not found in any dataset."))
      next  # Skip processing for this gene
    }
    
    
    # Generate the single violin plot for the current gene across ages
    plot <- ggplot(combined_data[combined_data$Expression != 0, ], aes(x = Age, y = Expression, fill = Age, color = Age)) +
      geom_violin(alpha = 0.5, trim = TRUE, scale = "count") +  # Set trim to TRUE
      geom_jitter(width = 0.1, alpha = 0.6, size = 0.5, aes(color = Age)) +  # Adjust jitter parameters and color
      geom_boxplot(width = 0.3, alpha = 0.8, fill = NA, color = "black", outlier.shape = NA) +  # Add black box plot without outliers
      scale_fill_manual(values = pastel_palette) +  # Use the pastel color palette for fill
      scale_color_manual(values = pastel_palette) +  # Set outline color to match fill color
      labs(x = "Age", y = "Expression", fill = "Age", color = "Age") +
      theme_minimal() +
      theme(legend.position = "top",
            panel.grid = element_blank())
    
    # Save the violin plot
    filename <- paste0(gene_name, "_violin_plot_", suture, ".",file_format)
    ggsave(filename = filename, plot = plot, path = output_dir, device = file_format)
  }
}

# # Example usage of the function
# seurat_obj_list <- list(E16_fs = E16.fs,
#                         E18_fs = E18.fs,
#                         P10_fs = P10.fs,
#                         P28_fs = P28.fs)
# 
# 
# 
# # Specify the path to the CSV file containing top genes
# top_genes_file <- here::here("midpalatal-sutures/visium/e15/data-output/WT1R1_all_markers_pct.csv")
# 
# # Specify the output directory (replace mm_figs with your desired output directory)
# output_dir <- mm_figs
# 
# # Specify the suture information, cluster number, and number of top genes to consider
# suture <- "frontal_suture"
# cluster_number <- 3
# top_genes_number <- 25
# 
# # Call the function with the list of Seurat objects, top genes file, output directory, suture, cluster number, and top genes number
# generate_individual_violin_plots_from_file(seurat_obj_list, top_genes_file, output_dir, suture, cluster_number, top_genes_number)
