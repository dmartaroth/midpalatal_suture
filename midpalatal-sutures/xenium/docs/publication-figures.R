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
             here(home.path, "pubfigs"),
           recursive = TRUE)

# Figure 1 ----------------------------------------------------------------


# Figure 2 ----------------------------------------------------------------


# Figure 3 ----------------------------------------------------------------


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



