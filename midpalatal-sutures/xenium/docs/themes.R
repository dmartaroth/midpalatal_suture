# ## ######################################## ## #
#                      THEMES                    #
# ## ######################################## ## #


# Define a custom ggplot2 theme function
custom_spatplot_theme <- function() {
  theme_minimal() +
    theme(
      text = element_text(size = 8),
      plot.background = element_rect(fill = "white", color = NA),  # Remove plot border
      axis.text = element_text(size = 6),
      plot.title = element_text(size = 8, face = "bold"),
      panel.grid = element_blank(),
      axis.text.x = element_text(margin = margin(t = 5)),
      axis.text.y = element_text(margin = margin(r = 5)),
      axis.title = element_text(size = 8, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 0)),  # Reduce bottom margin of x-axis title
      axis.title.y = element_text(margin = margin(r = 0)),  # Reduce right margin of y-axis title
      plot.title.position = "plot"
    )
}

# Define a custom ggplot2 scatterplot function
my_colors <- c("#FFB6C1", "#ADD8E6", "#FFD700", "#98FB98", "#FFA07A")
custom_scatter_theme <- function() {
  theme_minimal() +
    theme(
      text = element_text(size = 8),
      plot.background = element_rect(fill = "white", color = NA),  # Remove plot border
      axis.text = element_text(size = 8),  # Increase size of axis text
      axis.line = element_line(color = "black"),  # Set color of axis lines to black
      plot.title = element_text(size = 8, face = "bold", hjust = 0.5),  # Center the plot title
      panel.grid = element_blank(),
      axis.text.x = element_text(margin = margin(t = 5)),
      axis.text.y = element_text(margin = margin(r = 5)),
      axis.title = element_text(size = 8, face = "bold"),
      axis.title.x = element_text(margin = margin(t = 0)),  # Reduce bottom margin of x-axis title
      axis.title.y = element_text(margin = margin(r = 0)),  # Reduce right margin of y-axis title
      plot.title.position = "plot",
      plot.caption = element_blank(),  # Remove plot caption
      legend.position = "none",  # Remove legend
      legend.title = element_blank()  # Remove legend title
    ) 
}
