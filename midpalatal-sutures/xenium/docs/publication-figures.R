# ## ######################################## ## #
#                    FIGURE PLOTS                #
# ## ######################################## ## #

# Date: Mon Mar 25 14:17:17 2024 ------------------

# Load required libraries ----------------------------------------------------------
source("midpalatal-sutures/xenium/docs/packages.R")
source("midpalatal-sutures/xenium/docs/themes.R")
source("midpalatal-sutures/xenium/docs/functions.R")


# Set up Giotto environment -----------------------------------------------
# Set Giotto python path
python_path = NULL
if(is.null(python_path)) {
  installGiottoEnvironment()
}

dir.create(pubfig_folder <-
             here("pubfigs"),
           recursive = TRUE)

# Figure 1 ----------------------------------------------------------------


# Figure 2 ----------------------------------------------------------------


# Figure 3 ----------------------------------------------------------------

source("midpalatal-sutures/xenium/docs/spatial_correlation.R")
## panel I -----------------------------------------------------------------
# Sfrp and Tnn expression at e15.5
spatInSituPlotPoints(gobject,
                     show_image = FALSE,
                     feats = list('rna' = c("Sfrp2","Tnn")),
                     feats_color_code = c("Sfrp2" = "dodgerblue3","Tnn" = "green3"),
                     point_size = 1.5, 
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     show_legend=TRUE, 
                     polygon_line_size = 0.1,polygon_color = "lightpink",
                     polygon_alpha = 0.1,
                     axis_text=7,axis_title=7,
                     polygon_fill_as_factor = TRUE, 
                     polygon_fill = "cell_types",
                     polygon_fill_code = colorcode, 
                     background_color="white",
                     coord_fix_ratio = 1,
                     plot_last = c("points"),
                     save_param = list(
                       save_name = paste0("Fig3I_","R7_e15c", "_Sfrp_Tnn_spatplot"),
                       save_dir = pubfig_folder) )

## panel J -----------------------------------------------------------------
# Dcn and Eln expression at e15.5
spatInSituPlotPoints(gobject,
                     show_image = FALSE,
                     feats = list('rna' = c("Dcn","Eln")),
                     feats_color_code = c("Dcn" = "orange","Eln" = "violet"),
                     point_size = 1.5, 
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     show_legend=TRUE, 
                     polygon_line_size = 0.1,polygon_color = "lightpink",
                     polygon_alpha = 0.1,
                     axis_text=7,axis_title=7,
                     polygon_fill_as_factor = TRUE, 
                     polygon_fill = "cell_types",
                     polygon_fill_code = colorcode, 
                     background_color="white",
                     coord_fix_ratio = 1,
                     plot_last = c("points"),
                     save_param = list(
                       save_name = paste0("Fig3J_","R7_e15c", "_Dcn_Eln_spatplot"),
                       save_dir = pubfig_folder) )


# Figure 5 ----------------------------------------------------------------


## panel B -----------------------------------------------------------------
# proof of tenocytes in Visium


## panel C -----------------------------------------------------------------

# tenocyte markers in Xenium

# Subset entire palate instead of just left shelf
gobject <- loadGiotto(here("midpalatal-sutures/xenium/region-7/region-7_preprocessing_giotto-object"))

# Extract spatial data
spatial_data <- gobject@spatial_locs$cell$raw@coordinates
x_coordinates <- spatial_data$sdimx
y_coordinates <- spatial_data$sdimy

# Calculate initial parameters
x_min <- min(x_coordinates)
x_max <- max(x_coordinates)
y_min <- min(y_coordinates)
y_max <- max(y_coordinates)

# Generate and display initial plot with axis labels
initial_plot <- generate_plot(gobject, x_min, x_max, y_min, y_max, section, region)
print(initial_plot)

# Pause to examine the plot
# Next lines not working properly, debug
# cat("Examine the plot and then press Enter to continue...") 
# invisible(readline(prompt = ""))

# Loop to adjust parameters and generate new plots
while (TRUE) {
  # Ask for new parameters
  x_min <- as.numeric(readline("Enter new x_min: "))
  x_max <- as.numeric(readline("Enter new x_max: "))
  y_min <- as.numeric(readline("Enter new y_min: "))
  y_max <- as.numeric(readline("Enter new y_max: "))
  
  # Generate and display new plot
  new_plot <- generate_plot(gobject, x_min, x_max, y_min, y_max, section, region)
  print(new_plot)
  
  # Ask if the user wants to continue adjusting parameters
  continue_response <- readline("Adjust parameters again? (yes/no): ")
  if (tolower(continue_response) != "yes") {
    cat("Exiting plot adjustment.\n")
    break
  }
}

# Enter new x_min: 4850
# Enter new x_max: 6600
# Enter new y_min: 2800
# Enter new y_max: 3750

subset <-
  subsetGiottoLocs(
    gobject,
    x_min = x_min,
    x_max = x_max,
    y_min = y_min,
    y_max = y_max
  )

subset = calculateOverlapRaster(subset,
                                spatial_info = 'cell',
                                feat_info = 'rna')

showGiottoSpatialInfo(subset)

gobject <- subset
gobject <- overlapToMatrix(gobject,
                           poly_info = 'cell',
                           feat_info = 'rna',
                           name = 'raw')
showGiottoExpression(gobject)

panel_meta = data.table::fread(paste0(home.path,"/xenium_panel.tsv"))

# Append metadata
gobject <- addFeatMetadata(gobject = gobject,
                           feat_type = 'rna',
                           spat_unit = 'cell',
                           new_metadata = panel_meta,
                           by_column = TRUE,
                           column_feat_ID = 'feat_ID')


### Filter data and add stats -----------------------------------------------
gobject = filterGiotto(gobject = gobject,
                       spat_unit = 'cell',
                       poly_info = 'cell',
                       expression_threshold = 1,
                       feat_det_in_min_cells = 3,
                       min_det_feats_per_cell = 5)

gobject = addStatistics(gobject, expression_values = 'raw')

showGiottoCellMetadata(gobject)
showGiottoFeatMetadata(gobject)


### Normalize ---------------------------------------------------------------
gobject = normalizeGiotto(gobject = gobject,
                          spat_unit = 'cell',
                          scalefactor = 5000,
                          verbose = T)

spatInSituPlotDensity(gobject,
                      feats = c("Lum","Six2","Mkx", "Col3a1"),
                      feat_type = "rna",
                      polygon_alpha = 1,
                      polygon_color = "white",
                      background_color = "white",
                      cow_n_col = 1,
                      save_param = list(
                        save_name = paste0("Fig5C_","R7_e15c", "_Tenocyte_density"),
                        save_dir = pubfig_folder,save_format = "pdf"))

spatInSituPlotPoints(gobject,
                     show_image = FALSE,
                     feats = list('rna' = c("Dmp1","Mkx","Krt14","Chodl","Col2a1")),
                     feats_color_code = c("Mkx" = "green3","Dmp1" = "red2",
                                          "Krt14"="mediumpurple2","Chodl"="gold","Col2a1"="lightpink"),
                     point_size = 1.5,
                    show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     show_legend=TRUE, 
                     polygon_line_size = 0.1,polygon_color = "powderblue",
                     polygon_alpha = 0.05,
                    stroke = 0.1,
                     axis_text=7,axis_title=7,
                     polygon_fill_as_factor = FALSE, 
                     background_color="white",
                    polygon_bg_color = "white",
                     coord_fix_ratio = 1,
                     plot_last = c("points"),
                     save_param = list(
                       save_name = paste0("Fig5D_","R7_e15c", "_bone_Tenocyte_spatplot"),
                       save_dir = pubfig_folder, save_format="pdf") )



# Figure 6 ----------------------------------------------------------------


## panel F -----------------------------------------------------------------

# Midline and PNT enriched GFs and ECM feature plots plus spatplot

spatFeatPlot2D(gobject, expression_values = 'normalized', 
               feats = "Tnn",
               point_shape = 'no_border',
               gradient_midpoint = 0,
               cell_color_gradient = c("skyblue", "bisque1", "red3"),
               show_network = FALSE, point_size = 2,
               cow_n_col = 1,
               axis_text = 10,
               axis_title = 10,
               save_param = list(
                 save_name = paste0("Fig6F_","R7_e15c", "_Tnn_featplot"),
                 save_dir = pubfig_folder, save_format="pdf", base_width = 4, base_height = 3))

spatFeatPlot2D(gobject, expression_values = 'normalized', 
               feats = "Mkx",
               point_shape = 'no_border',
               gradient_midpoint = 0,
               cell_color_gradient = c("skyblue", "bisque1", "red3"),
               show_network = FALSE, point_size = 2,
               cow_n_col = 1,
               axis_text = 10,
               axis_title = 10,
               save_param = list(
                 save_name = paste0("Fig6F_","R7_e15c", "_Mkx_featplot"),
                 save_dir = pubfig_folder, save_format="pdf", base_width = 4, base_height = 3))

spatFeatPlot2D(gobject, expression_values = 'normalized', 
               feats = "Tgfb2",
               point_shape = 'no_border',
               gradient_midpoint = 0,
               cell_color_gradient = c("skyblue", "bisque1", "red3"),
               show_network = FALSE, point_size = 2,
               cow_n_col = 1,
               axis_text = 10,
               axis_title = 10,
               save_param = list(
                 save_name = paste0("Fig6F_","R7_e15c", "_Tgfb2_featplot"),
                 save_dir = pubfig_folder, save_format="pdf", base_width = 4, base_height = 3))

spatFeatPlot2D(gobject, expression_values = 'normalized', 
               feats = "Thbs1",
               point_shape = 'no_border',
               gradient_midpoint = 0,
               cell_color_gradient = c("skyblue", "bisque1", "red3"),
               show_network = FALSE, point_size = 2,
               cow_n_col = 1,
               axis_text = 10,
               axis_title = 10,
               save_param = list(
                 save_name = paste0("Fig6F_","R7_e15c", "_Thbs1_featplot"),
                 save_dir = pubfig_folder, save_format="pdf", base_width = 4, base_height = 3))

spatFeatPlot2D(gobject, expression_values = 'normalized', 
               feats = "Bmp2",
               point_shape = 'no_border',
               gradient_midpoint = 0,
               cell_color_gradient = c("skyblue", "bisque1", "red3"),
               show_network = FALSE, point_size = 2,
               cow_n_col = 1,
               axis_text = 10,
               axis_title = 10,
               save_param = list(
                 save_name = paste0("Fig6F_","R7_e15c", "_Bmp2_featplot"),
                 save_dir = pubfig_folder, save_format="pdf", base_width = 4, base_height = 3))


spatInSituPlotPoints(gobject,
                     show_image = FALSE,
                     feats = list('rna' = c("Mkx","Tnn","Tgfb2","Thbs1","Gsc","Dcn","Bmp2")),
                     feats_color_code = c("Mkx" = "white","Tnn" = "steelblue",
                                          "Tgfb2"="magenta","Thbs1"="gold",
                                          "Gsc"="greenyellow","Dcn"="cyan",
                                          "Bmp2"="red"),
                     point_size = 1.5,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     show_legend=TRUE, 
                     polygon_line_size = 0.4,polygon_color = "dodgerblue4",
                     polygon_alpha = 0.01,
                     stroke = 0.5,
                     axis_text=7,axis_title=7,
                     polygon_fill_as_factor = FALSE, 
                     legend_text = 7,
                     background_color="navyblue",
                     feat_shape_code = list('rna'=10),
                     polygon_bg_color = "white",
                     coord_fix_ratio = 1,
                     plot_last = c("points"),
                     plot_method = "ggplot",
                     save_param = list(
                       save_name = paste0("Fig6F_","R7_e15c", "_PNT-mm_spatplot"),
                       save_dir = pubfig_folder, save_format="pdf") )

