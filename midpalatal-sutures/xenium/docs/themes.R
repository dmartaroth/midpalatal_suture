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

# Custom theme for barplot
custom_theme_bar <- function() {
  theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(color = "white", fill = "white"),
      panel.grid.major = element_blank(),           # No major gridlines
      panel.grid.minor = element_blank(),           # No minor gridlines
      legend.position = "none",                     # No legend
      axis.ticks = element_line(color = "black"),
      axis.ticks.x = element_line(color = NA),
      axis.title = element_text(color = "black",size = 8),
      axis.text = element_text(size = 8),           # Size of axis text
      axis.line.y = element_line(colour = NA), # Set y-axis line color
      axis.line.x = element_line(colour = NA),  
      axis.text.y = element_text(colour = "black"), # Set y-axis text color
      axis.ticks.length = unit(0.2, "cm"),          # Shorten tick length
      panel.border = element_blank(),   
      panel.grid.major.y = element_line(colour = "gray", linetype = 3),  # Dashed horizontal gridlines
      panel.grid.major.x = element_blank(),         # No vertical gridlines
      text = element_text(colour = "black", size = 8),                        # Regular text
      plot.title = element_text(face = "bold", size = 8,hjust = 0),     # Bold plot title
      # plot.margin = margin(20, 20, 20, 20),         # Adjust plot margins
      plot.caption = element_text(color = "black")  # Caption color and position
    )
}
feat_colors <- c("firebrick1", "darkorchid4", "greenyellow","dodgerblue", "aquamarine", 
                 "magenta", "lawngreen", "palevioletred2", "gold", "cadetblue3", 
                 "lightcyan", "violetred1", "lawngreen", "darkorchid4",
                 "firebrick1", "blue", "orchid4", "salmon1", "seagreen3",
                 "steelblue4", "royalblue2", "pink4", "cadetblue3", "red2",
                 "chocolate1", "dodgerblue4", "darkolivegreen4", "magenta2",
                 "skyblue2", "seagreen", "palevioletred2", "mediumpurple3",
                 "aquamarine")
