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


load_xenium_data <- function() {
  # Create Giotto instructions for saving
  instrs <- createGiottoInstructions(save_dir = results_folder,
                                     save_plot = TRUE,
                                     show_plot = FALSE,
                                     return_plot = TRUE)
  
  # Check if features-blank.tsv.gz file exists
  if (!file.exists(file.path(xenium_folder, "cell_feature_matrix", "features-blank.tsv.gz"))) {
    # Read features.tsv.gz and rename "Unassigned Codeword" to "Blank Codeword"
    features <- read_tsv(here(xenium_folder, "cell_feature_matrix", "features.tsv.gz"))
    features[features == "Unassigned Codeword"] <- "Blank Codeword"
    
    # Write the renamed features to features-blank.tsv.gz
    write_tsv(x = features,
              file = here(xenium_folder, "cell_feature_matrix", "features-blank.tsv.gz"))
    
    # Print messages
    cat(bold(magenta("This package looks for Blank Codeword in the place of Unassigned Codeword, so first we have to rename them in the features.tsv.doc.\n")))
    cat(bold(magenta("It looks like you haven't run this code before. Blank codewords have now been successfully renamed and the file has been saved as features-blank.tsv.gz.\n")))
  } else {
    # Print messages
    cat(bold(magenta("This package looks for Blank Codeword in the place of Unassigned Codeword, so first we have to rename them from the features.tsv.doc.\n")))
    cat(bold(magenta("The file 'features-blank.tsv.gz' already exists. Skipping renaming.\n")))
  }
  
  # General files (some are supplemental files)
  settings_path <- paste0(xenium_folder, '/experiment.xenium')
  he_img_path <- paste0(xenium_folder, '/pyramidalhe.ome.tif')
  panel_meta_path <- paste0(xenium_folder, "/xenium_panel.tsv") # (optional)
  
  # Files (SUBCELLULAR): 
  cell_bound_path <- paste0(xenium_folder, '/cell_boundaries.csv.gz')
  nuc_bound_path <- paste0(xenium_folder, '/nucleus_boundaries.csv.gz')
  tx_path <- paste0(xenium_folder, '/transcripts.csv.gz')
  feat_meta_path <- paste0(xenium_folder, '/cell_feature_matrix/features-blank.tsv.gz')
  
  # Files (AGGREGATE):
  expr_mat_path <- paste0(xenium_folder, '/cell_feature_matrix')
  cell_meta_path <- paste0(xenium_folder, '/cells.csv.gz') # contains spatlocs
  
  # Load features metadata
  feature_dt <- data.table::fread(feat_meta_path, header = FALSE)
  colnames(feature_dt) <- c('ensembl_ID', 'feat_name', 'feat_type')
  
  # Find the feature IDs that belong to each feature type
  feature_dt[, table(feat_type)]
  feat_types <- names(feature_dt[, table(feat_type)])
  feat_types_IDs <- lapply(feat_types, function(type)
    feature_dt[feat_type == type, unique(feat_name)])
  names(feat_types_IDs) <- feat_types
  
  # Load transcript-level data
  tx_dt <- data.table::fread(tx_path)
  data.table::setnames(
    x = tx_dt,
    old = c('feature_name', 'x_location', 'y_location'),
    new = c('feat_ID', 'x', 'y')
  )
  cat(
    'Transcripts info available:\n ',
    paste0('"', colnames(tx_dt), '"'),
    '\n',
    'with',
    tx_dt[, .N],
    'unfiltered detections\n'
  )
  
  # Filter by qv (Phred score)
  tx_dt_filtered <- tx_dt[qv >= 20]
  cat('and', tx_dt_filtered[, .N], 'filtered detections\n\n')
  
  # Separate detections by feature type
  tx_dt_types <- lapply(feat_types_IDs, function(types)
    tx_dt_filtered[feat_ID %in% types])
  
  invisible(lapply(seq_along(tx_dt_types), function(x) {
    cat(names(tx_dt_types)[[x]], ' detections: ', tx_dt_types[[x]][, .N], '\n')
  }))
  
  # Return the paths and data
  return(list(
    settings_path = settings_path,
    he_img_path = he_img_path,
    panel_meta_path = panel_meta_path,
    cell_bound_path = cell_bound_path,
    nuc_bound_path = nuc_bound_path,
    tx_path = tx_path,
    feat_meta_path = feat_meta_path,
    expr_mat_path = expr_mat_path,
    cell_meta_path = cell_meta_path,
    feature_dt = feature_dt,
    feat_types_IDs = feat_types_IDs,
    tx_dt = tx_dt,
    tx_dt_filtered = tx_dt_filtered,
    tx_dt_types = tx_dt_types
  ))
}


load_polygon_data <- function() {
  # Load cell polygon data
  cellPoly_dt <- data.table::fread(cell_bound_path)
  data.table::setnames(
    cellPoly_dt,
    old = c('cell_id', 'vertex_x', 'vertex_y'),
    new = c('poly_ID', 'x', 'y')
  )
  
  # Load nucleus polygon data
  nucPoly_dt <- data.table::fread(nuc_bound_path)
  data.table::setnames(
    nucPoly_dt,
    old = c('cell_id', 'vertex_x', 'vertex_y'),
    new = c('poly_ID', 'x', 'y')
  )
  
  # Create polygons from data frames
  gpoly_cells <- createGiottoPolygonsFromDfr(segmdfr = cellPoly_dt,
                                             name = 'cell',
                                             calc_centroids = TRUE)
  gpoly_nucs <- createGiottoPolygonsFromDfr(segmdfr = nucPoly_dt,
                                            name = 'nucleus',
                                            calc_centroids = TRUE)
  
  # Return the created polygons
  return(list(gpoly_cells = gpoly_cells, gpoly_nucs = gpoly_nucs))
}
