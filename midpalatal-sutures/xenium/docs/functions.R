# ## ######################################## ## #
#                    FUNCTIONS                   #
# ## ######################################## ## #

# Function to generate plot with axis labels
generate_plot <- function(gobject, x_min, x_max, y_min, y_max, section, region) {
  subset <- subsetGiottoLocs(gobject, x_min = x_min, x_max = x_max, y_min = y_min, y_max = y_max)
  p <- ggplot(data = subset@spatial_locs$cell$raw@coordinates, aes(x = sdimx, y = sdimy)) +
    geom_point() +
    labs(title = paste0(section, " ", region), x = "X Axis Label", y = "Y Axis Label") +
    theme_minimal()
  return(p)
}

title <- paste0(section, " ", region)

# Define the output folder as a global variable


log_step <- function(step_title) {
  # Open the log file in append mode
  outs <- here(section_folder, "data-output")
  outs <- file(here::here(outs, "log.txt"), "a")
  
  # Get the current number of lines in the log file
  num_lines <- length(readLines(outs))
  
  # Write the step title with a sequential number
  cat(paste(num_lines + 1, ". ", step_title, "\n"), file = outs, append = TRUE)
  
  # Close the log file
  close(outs)
}
