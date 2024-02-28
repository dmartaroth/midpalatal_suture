# ## ######################################## ## #
#                 FUNCTIONS IN DEV               #
# ## ######################################## ## #


# Initialize plot number
plot_number <- 0




# Function to extract plot function name from ggplot object
extract_plot_function <- function(plot) {
  plot_function_name <- deparse(substitute(plot))
  cleaned_name <- gsub("^\\s*|\\s*$", "", plot_function_name)
  cleaned_name <- gsub("_", "", cleaned_name)  # Remove underscores
  return(cleaned_name)
}



# Function to extract variable name from assignment statement
extract_variable_name <- function(expr) {
  as.character(substitute(expr))
}



# Function to generate filename with sequential numbering, variable name, section, and region
generate_filename <- function(plot_function_name, section, region, plot_number, variable_name) {
  plot_number_formatted <- sprintf("%02d", plot_number)
  filename <- paste0(plot_number_formatted, "_", variable_name, "_", plot_function_name, "_", section, "_", region, ".png")
  return(filename)
}




# Function to save plot with sequential numbering and naming
save_plot <- function(plot, section, region, save_directory, plot_number, variable_name) {
  # Extract plot function name from object name
  plot_function_name <- extract_plot_function(plot)
  
  # Increment plot number
  plot_number <- plot_number + 1
  
  # Generate filename
  filename <- generate_filename(plot_function_name, section, region, plot_number, variable_name)
  
  # Save the plot with the constructed filename and directory specified
  ggsave(file.path(save_directory, filename), plot, width = 6, height = 4, dpi = 300)
  
  # Return the incremented plot number
  return(plot_number)
}

# Example usage:
# Run the plot function and store the plot in a variable


title <- paste0(section, " ", region)
