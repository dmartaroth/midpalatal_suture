# ## ######################################## ## #
#                    FUNCTIONS                   #
# ## ######################################## ## #
# Set up Giotto environment -----------------------------------------------
# Set Giotto python path
python_path = NULL
if(is.null(python_path)) {
  installGiottoEnvironment()
}


waitForInput <- function(message) {
  shinyApp(
    ui = fluidPage(
      style = "background-color: #fffcf1; font-size: 200%; width: 4in; height: 2in;",  # Set background color to #fffcf1, increase font size by 2 times, and adjust width and height
      fluidRow(
        column(12, h5(message, style = "font-size: 16px; color: #8a76b2; text-align: center; font-style: italic;")),  # Set message text size to 16px, color to #8a76b2, center the text, and italicize it
        column(12, offset = 4, align = "center", actionButton("doneButton", "DONE", style = "background-color: #f9f6ff; color: #4f4366; font-weight: bold; padding: 16px 50px; border: 4px solid #c6a9ff; border-radius: 8px; cursor: pointer; min-width: 120px; width: auto;"))  # Center the button below the message text
      ),
      tags$style("body { font-size: 22px; } .btn { font-size: 19px; }"),  # Increase font size for all text and button text
      width = "66%",  # Set the width of the page to be reduced by 1/3
      height = "100px"  # Set the height of the page to be smaller
    ),
    server = function(input, output, session) {
      observeEvent(input$doneButton, {
        stopApp()  # Stop the Shiny app when the button is clicked
      })
    }
  )
}



# Function to generate plot with axis labels
generate_plot <- function(gobject, x_min, x_max, y_min, y_max, section, region) {
  subset <- subsetGiottoLocs(gobject, x_min = x_min, x_max = x_max, y_min = y_min, y_max = y_max)
  p <- ggplot(data = subset@spatial_locs$cell$raw@coordinates, aes(x = sdimx, y = sdimy)) +
    geom_point() +
    labs(title = paste0(section, " ", region), x = "X Axis Label", y = "Y Axis Label") +
    theme_minimal()
  return(p)
  title <- paste0(section, " ", region)
}



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



