# Define UI
ui <- fluidPage(
  titlePanel("Spatial Transcriptomics Plot"),
  sidebarLayout(
    sidebarPanel(
      textInput("features", "Enter features (comma-separated):", placeholder = "e.g., gene1, gene2"),
      actionButton("plotBtn", "Generate Plot"),
      br(),
      downloadButton("downloadPlot", "Save Plot as PNG")
    ),
    mainPanel(
      plotOutput("spatialPlotPreview")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  plot <- NULL
  
  observeEvent(input$plotBtn, {
    # Split user input features by comma and trim whitespace
    features <- trimws(strsplit(input$features, ",")[[1]])
    
    # Generate spatial plot using spatFeatPlot2D
    plot <- spatFeatPlot2D(gobject, expression_values = 'normalized', 
                           feats = features,
                           point_shape = 'no_border',
                           gradient_midpoint = 0,
                           cell_color_gradient = c("skyblue", "bisque1", "red3"),
                           show_network = FALSE, point_size = 1.5,
                           cow_n_col = 3,
                           axis_text = 5,
                           axis_title = 5)
    
    # Display the plot in the preview area
    output$spatialPlotPreview <- renderPlot({
      plot
    })
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      # Set filename based on input features
      if (!is.null(input$features)) {
        features <- trimws(strsplit(input$features, ",")[[1]])
        filename <- paste(features, collapse = "_") %>% paste0("_spatial_plot.png")
      } else {
        filename <- "spatial_plot.png"
      }
      
      filename
    },
    content = function(file) {
      # Split user input features by comma and trim whitespace
      features <- trimws(strsplit(input$features, ",")[[1]])
      
      # Generate spatial plot using spatFeatPlot2D
      plot <- spatFeatPlot2D(gobject, expression_values = 'normalized', 
                             feats = features,
                             point_shape = 'no_border',
                             gradient_midpoint = 0,
                             cell_color_gradient = c("skyblue", "bisque1", "red3"),
                             show_network = FALSE, point_size = 3,
                             cow_n_col = 1,
                             axis_text = 5,
                             axis_title = 5)
      
      # Save plot as a PNG file in the specified directory with the specified filename
      png(file, width = 300, height = 175)
      print(plot)
      dev.off()
    }
  )
}

# Run the Shiny app
shinyApp(ui = ui, server = server)

