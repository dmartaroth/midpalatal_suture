# ## ###################### ## #
#  SUBSETTING XENIUM REGIONS   #
# ## ###################### ## #

# Run this code after script 01 to subset regions into individual sections
# This code currently has to be run manually one line at a time; need to make
# corrections before whole script can just be sourced.
# Updated by Daniela M. Roth on
# Date: Mon Feb 26 14:48:07 2024 ------------------


# Load required libraries ----------------------------------------------------------
source("midpalatal-sutures/xenium/docs/packages.R")
source("midpalatal-sutures/xenium/docs/themes.R")
source("midpalatal-sutures/xenium/docs/functions.R")

# Definitions and directories ---------------------------------------------

# # Set paths for overall Xenium data
# here() # check where top level is located
# home.path <- here("midpalatal-sutures", "xenium")


# Set up Giotto environment -----------------------------------------------
# Set Giotto python path
python_path = NULL
if(is.null(python_path)) {
  installGiottoEnvironment()
}


# Region 5 ----------------------------------------------------------------
# Created whole-region gobject.RDS in last script

region <- "region-5" # replace region-x with correct region name for each
home.path <- here("midpalatal-sutures","xenium")
prepro.folder <- paste0(region,"_preprocessing_","giotto-object")

## Section P0a -------------------------------------------------------------
# Directories
section <- "P0a"

dir.create(section_folder <-
             here(home.path, region, section),
           recursive = TRUE)

dir.create(results_folder <-
             here(section_folder, "figs"),
           recursive = TRUE)

dir.create(output <- here(section_folder, "data-output"))
# # Function to generate plot with axis labels
# generate_plot <- function(gobject, x_min, x_max, y_min, y_max, section, region) {
#   subset <- subsetGiottoLocs(gobject, x_min = x_min, x_max = x_max, y_min = y_min, y_max = y_max)
#   p <- ggplot(data = subset@spatial_locs$cell$raw@coordinates, aes(x = sdimx, y = sdimy)) +
#     geom_point() +
#     labs(title = paste0(section, " ", region), x = "X Axis Label", y = "Y Axis Label") +
#     theme_minimal()
#   return(p)
# }

# Load gobject
region <- "region-5"
# prepro.folder <- paste0(region,"_preprocessing_","giotto-object")
gobject <- loadGiotto(here(home.path, region, "data-output", prepro.folder))
# 
# # Define section
# section <- "P0a"

# Create directories
# section_folder <- here(home.path, region, section)
# dir.create(section_folder, recursive = TRUE)
# results_folder <- here(section_folder, "figs")
# dir.create(results_folder, recursive = TRUE)
# output <- here(section_folder, "data-output")

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


# # Record parameters here for future use:
# x_min <- 400
# x_max <- 3100
# y_min <- 300
# y_max <- 2900


### Subset section ----------------------------------------------------------
subset <-
  subsetGiottoLocs(
    gobject,
    x_min = x_min,
    x_max = x_max,
    y_min = y_min,
    y_max = y_max
  )

# 
# # Define a custom ggplot2 theme function
# custom_spatplot_theme <- function() {
#   theme_minimal() +
#     theme(
#       text = element_text(size = 8),
#       plot.background = element_rect(fill = "white", color = NA),  # Remove plot border
#       axis.text = element_text(size = 6),
#       plot.title = element_text(size = 8, face = "bold"),
#       panel.grid = element_blank(),
#       axis.text.x = element_text(margin = margin(t = 5)),
#       axis.text.y = element_text(margin = margin(r = 5)),
#       axis.title = element_text(size = 8, face = "bold"),
#       axis.title.x = element_text(margin = margin(t = 0)),  # Reduce bottom margin of x-axis title
#       axis.title.y = element_text(margin = margin(r = 0)),  # Reduce right margin of y-axis title
#       plot.title.position = "plot"
#     )
# }

# Run the spatPlot2D function and store the plot in a variable
filename <- paste0("01_", section, "_", region, "_subset-spatPlot2D.png")
# title <- paste0(section, " ", region)
spatPlot <- spatPlot2D(
  subset,
  spat_unit = 'cell',
  title = title,
  point_shape = 'no_border',
  point_size = 0.5,
  point_alpha = 0.4,
  return_plot = TRUE,
  save_plot = FALSE)

# Apply the custom theme to the plot
spatPlot <- spatPlot + custom_spatplot_theme()

# Now you can save the plot with the custom theme applied
ggsave(file.path(results_folder, filename), spatPlot, width = 6, height = 4, dpi = 300)

### Store subset metadata ---------------------------------------------------

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


### Calculate Highly Variable Features --------------------------------------
gobject <- calculateHVF(gobject = gobject,
                        spat_unit = 'cell',
                        save_param = list(save_name = paste0("02_", section, "_", region, "_HVF"),
                                          save_dir = results_folder),
                        return_plot = TRUE)


cat(fDataDT(gobject)[, sum(hvf == 'yes')], 'hvf found')


### Dimensional reduction ---------------------------------------------------

# # Define a custom ggplot2 scatterplot function
# my_colors <- c("#FFB6C1", "#ADD8E6", "#FFD700", "#98FB98", "#FFA07A")
# custom_scatter_theme <- function() {
#   theme_minimal() +
#     theme(
#       text = element_text(size = 8),
#       plot.background = element_rect(fill = "white", color = NA),  # Remove plot border
#       axis.text = element_text(size = 8),  # Increase size of axis text
#       axis.line = element_line(color = "black"),  # Set color of axis lines to black
#       plot.title = element_text(size = 8, face = "bold", hjust = 0.5),  # Center the plot title
#       panel.grid = element_blank(),
#       axis.text.x = element_text(margin = margin(t = 5)),
#       axis.text.y = element_text(margin = margin(r = 5)),
#       axis.title = element_text(size = 8, face = "bold"),
#       axis.title.x = element_text(margin = margin(t = 0)),  # Reduce bottom margin of x-axis title
#       axis.title.y = element_text(margin = margin(r = 0)),  # Reduce right margin of y-axis title
#       plot.title.position = "plot",
#       plot.caption = element_blank(),  # Remove plot caption
#       legend.position = "none",  # Remove legend
#       legend.title = element_blank()  # Remove legend title
#     ) 
# }


#### PCA ---------------------------------------------------------------------

gobject = runPCA(gobject = gobject,
                     spat_unit = 'cell',
                     expression_values = 'scaled',
                     feats_to_use = NULL,
                     scale_unit = F,
                     center = F)


# Visualize Screeplot and PCA
# Create the plot
p1 <- screePlot(gobject,
                ncp = 20,
                save_plot = TRUE,
                save_param = list(save_name = paste0("03a_", section, "_", region, "_screePlot"),
                                  save_dir = results_folder),
                return_plot = TRUE)

showGiottoDimRed(gobject)

p2 <- plotPCA(gobject,
              spat_unit = 'cell',
              dim_reduction_name = 'pca',
              dim1_to_use = 1,
              dim2_to_use = 2, 
              save_param = list(save_name = paste0("03b_", section, "_", region, "_PCA"),
                                                 save_dir = results_folder),
              return_plot = T,
              save_plot = T)


# library(cowplot)
screeplot_pca <- plot_grid(p1,p2,ncol = 1)

filename <- paste0("03_", section, "_", region, "_combined_screePlot_plotPCA.png")
ggsave(file.path(results_folder, filename), screeplot_pca, width = 4, height = 8, dpi = 300)


#### tSNE and UMAP -----------------------------------------------------------
gobject = runtSNE(gobject,
                      dimensions_to_use = 1:10,
                      spat_unit = 'cell',check_duplicates=FALSE)

gobject = runUMAP(gobject,
                      dimensions_to_use = 1:10,
                      spat_unit = 'cell')

p1 <- plotTSNE(gobject,
               point_size = 0.005,
               save_param = list(
                 save_name = paste0("04a_", section, "_", region, "_tSNE"),
                                 save_dir = results_folder),return_plot=T)

p2 <- plotUMAP(gobject,
               point_size = 0.005,
               save_param = list(
                 save_name = paste0("04b_", section, "_", region, "_UMAP"),
                 save_dir = results_folder),return_plot=T)

tSNE_UMAP <- plot_grid(p1,p2,ncol = 2)
filename <- paste0("04_", section, "_", region, "_combined_tSNE_UMAP.png")
ggsave(file.path(results_folder, filename), tSNE_UMAP, width = 5, height = 3, dpi = 300)



log_step("Normalization")










# Region 6 ----------------------------------------------------------------

## Create region-specific directories and set paths ----------------------------
region <- "region-6" # replace region-x with correct region name for each
dir.create(here(home.path, region))
dir.create(results_folder <- here(home.path,region,"figs"))
dir.create(xenium_folder <-
             here(home.path, region, "raw-data")) # place raw data here

# Create Giotto instructions for saving
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  return_plot = TRUE)

# This package looks for "Blank Codeword" in the place of "Unassigned Codeword"
# Rename in features.tsv.doc before proceeding to next path definitions
# Uncomment next lines to rename Codewords
# features <-
#   read_tsv(here(xenium_folder, "cell_feature_matrix", "features.tsv.gz"))
# View(features)
# features[features == "Unassigned Codeword"] <- "Blank Codeword"
# write_tsv(
#   x = features,
#   file = here(
#     xenium_folder,
#     "cell_feature_matrix",
#     "features-blank.tsv.gz"
#   )
# )

# General files (some are supplemental files)
settings_path = paste0(xenium_folder, '/experiment.xenium')
he_img_path = paste0(xenium_folder, '/pyramidalhe.ome.tif')
panel_meta_path = paste0(xenium_folder, "/xenium_panel.tsv") # (optional)

# Files (SUBCELLULAR): 
cell_bound_path = paste0(xenium_folder, '/cell_boundaries.csv.gz')
nuc_bound_path = paste0(xenium_folder, '/nucleus_boundaries.csv.gz')
tx_path = paste0(xenium_folder, '/transcripts.csv.gz')
feat_meta_path = paste0(xenium_folder, '/cell_feature_matrix/features-blank.tsv.gz')

# Files (AGGREGATE):
expr_mat_path = paste0(xenium_folder, '/cell_feature_matrix')
cell_meta_path = paste0(xenium_folder, '/cells.csv.gz') # contains spatlocs


## Load features metadata --------------------------------------------------
# Make sure cell_feature_matrix folder is unpacked
feature_dt = data.table::fread(feat_meta_path, header = FALSE)
colnames(feature_dt) = c('ensembl_ID', 'feat_name', 'feat_type')

# Find the feature IDs that belong to each feature type
feature_dt[, table(feat_type)]
feat_types = names(feature_dt[, table(feat_type)])
feat_types_IDs = lapply(feat_types, function(type)
  feature_dt[feat_type == type, unique(feat_name)])
names(feat_types_IDs) = feat_types


## Load transcript-level data ----------------------------------------------
tx_dt = data.table::fread(tx_path)
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
tx_dt_filtered = tx_dt[qv >= 20]
cat('and', tx_dt_filtered[, .N], 'filtered detections\n\n')

# Separate detections by feature type
tx_dt_types = lapply(feat_types_IDs, function(types)
  tx_dt_filtered[feat_ID %in% types])

invisible(lapply(seq_along(tx_dt_types), function(x) {
  cat(names(tx_dt_types)[[x]], 'detections: ', tx_dt_types[[x]][, .N], '\n')
}))


## Preview region ----------------------------------------------------------
gpoints_list = lapply(tx_dt_types, function(x)
  createGiottoPoints(x = x))

# Preview QC probe detections
plot(gpoints_list$`Blank Codeword`,
     point_size = 1,
     main = 'Blank Codeword')
plot(gpoints_list$`Negative Control Codeword`,
     point_size = 1,
     main = 'Negative Control Codeword')
plot(gpoints_list$`Negative Control Probe`,
     point_size = 1,
     main = 'Negative Control Probe')

# Preview selected genes
mygenes <- c("Krt14","Col1a1","Tnn","Lum")

filename <-
  paste0(
    results_folder,
    "/01_preview-region_",
    str_c(mygenes, "_", collapse = ""),
    region,
    ".pdf"
  )
pdf(file =filename, width=4, height=4)
preview.spatplot <- plot(gpoints_list$`Gene Expression`,
                         feats = mygenes)+theme(legend.position = "bottom")
dev.off()

tx_dt_types$`Gene Expression`[feat_ID %in% mygenes, table(feat_ID)]


## Load polygon data -------------------------------------------------------
cellPoly_dt = data.table::fread(cell_bound_path)
nucPoly_dt = data.table::fread(nuc_bound_path)

data.table::setnames(
  cellPoly_dt,
  old = c('cell_id', 'vertex_x', 'vertex_y'),
  new = c('poly_ID', 'x', 'y')
)
data.table::setnames(
  nucPoly_dt,
  old = c('cell_id', 'vertex_x', 'vertex_y'),
  new = c('poly_ID', 'x', 'y')
)

gpoly_cells = createGiottoPolygonsFromDfr(segmdfr = cellPoly_dt,
                                          name = 'cell',
                                          calc_centroids = TRUE)
gpoly_nucs = createGiottoPolygonsFromDfr(segmdfr = nucPoly_dt,
                                         name = 'nucleus',
                                         calc_centroids = TRUE)


## Create Giotto Object for entire slide -----------------------------------
gobject = createGiottoObjectSubcellular(
  gpoints = list(
    rna = gpoints_list$`Gene Expression`,
    blank_code = gpoints_list$`Blank Codeword`,
    neg_code = gpoints_list$`Negative Control Codeword`,
    neg_probe = gpoints_list$`Negative Control Probe`
  ),
  gpolygons = list(cell = gpoly_cells,
                   nucleus = gpoly_nucs),
  instructions = instrs
)

saveGiotto(gobject,dir = here(home.path, region),foldername = paste0(region,"_preprocessing_","giotto-object"), overwrite = TRUE)
# File will save as gobject.RDS so descriptive folder naming is essential


# Region 7 ----------------------------------------------------------------

## Create region-specific directories and set paths ----------------------------
region <- "region-7" # replace region-x with correct region name for each
dir.create(here(home.path, region))
dir.create(results_folder <- here(home.path,region,"figs"))
dir.create(xenium_folder <-
             here(home.path, region, "raw-data")) # place raw data here

# Create Giotto instructions for saving
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  return_plot = TRUE)

# This package looks for "Blank Codeword" in the place of "Unassigned Codeword"
# Rename in features.tsv.doc before proceeding to next path definitions
# Uncomment next lines to rename Codewords
# features <-
#   read_tsv(here(xenium_folder, "cell_feature_matrix", "features.tsv.gz"))
# View(features)
# features[features == "Unassigned Codeword"] <- "Blank Codeword"
# write_tsv(
#   x = features,
#   file = here(
#     xenium_folder,
#     "cell_feature_matrix",
#     "features-blank.tsv.gz"
#   )
# )

# General files (some are supplemental files)
settings_path = paste0(xenium_folder, '/experiment.xenium')
he_img_path = paste0(xenium_folder, '/pyramidalhe.ome.tif')
panel_meta_path = paste0(xenium_folder, "/xenium_panel.tsv") # (optional)

# Files (SUBCELLULAR): 
cell_bound_path = paste0(xenium_folder, '/cell_boundaries.csv.gz')
nuc_bound_path = paste0(xenium_folder, '/nucleus_boundaries.csv.gz')
tx_path = paste0(xenium_folder, '/transcripts.csv.gz')
feat_meta_path = paste0(xenium_folder, '/cell_feature_matrix/features-blank.tsv.gz')

# Files (AGGREGATE):
expr_mat_path = paste0(xenium_folder, '/cell_feature_matrix')
cell_meta_path = paste0(xenium_folder, '/cells.csv.gz') # contains spatlocs


## Load features metadata --------------------------------------------------
# Make sure cell_feature_matrix folder is unpacked
feature_dt = data.table::fread(feat_meta_path, header = FALSE)
colnames(feature_dt) = c('ensembl_ID', 'feat_name', 'feat_type')

# Find the feature IDs that belong to each feature type
feature_dt[, table(feat_type)]
feat_types = names(feature_dt[, table(feat_type)])
feat_types_IDs = lapply(feat_types, function(type)
  feature_dt[feat_type == type, unique(feat_name)])
names(feat_types_IDs) = feat_types


## Load transcript-level data ----------------------------------------------
tx_dt = data.table::fread(tx_path)
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
tx_dt_filtered = tx_dt[qv >= 20]
cat('and', tx_dt_filtered[, .N], 'filtered detections\n\n')

# Separate detections by feature type
tx_dt_types = lapply(feat_types_IDs, function(types)
  tx_dt_filtered[feat_ID %in% types])

invisible(lapply(seq_along(tx_dt_types), function(x) {
  cat(names(tx_dt_types)[[x]], 'detections: ', tx_dt_types[[x]][, .N], '\n')
}))


## Preview region ----------------------------------------------------------
gpoints_list = lapply(tx_dt_types, function(x)
  createGiottoPoints(x = x))

# Preview QC probe detections
plot(gpoints_list$`Blank Codeword`,
     point_size = 1,
     main = 'Blank Codeword')
plot(gpoints_list$`Negative Control Codeword`,
     point_size = 1,
     main = 'Negative Control Codeword')
plot(gpoints_list$`Negative Control Probe`,
     point_size = 1,
     main = 'Negative Control Probe')

# Preview selected genes
mygenes <- c("Krt14","Col1a1","Tnn","Lum")

filename <-
  paste0(
    results_folder,
    "/01_preview-region_",
    str_c(mygenes, "_", collapse = ""),
    region,
    ".pdf"
  )
pdf(file =filename, width=4, height=4)
preview.spatplot <- plot(gpoints_list$`Gene Expression`,
                         feats = mygenes)+theme(legend.position = "bottom")
dev.off()

tx_dt_types$`Gene Expression`[feat_ID %in% mygenes, table(feat_ID)]


## Load polygon data -------------------------------------------------------
cellPoly_dt = data.table::fread(cell_bound_path)
nucPoly_dt = data.table::fread(nuc_bound_path)

data.table::setnames(
  cellPoly_dt,
  old = c('cell_id', 'vertex_x', 'vertex_y'),
  new = c('poly_ID', 'x', 'y')
)
data.table::setnames(
  nucPoly_dt,
  old = c('cell_id', 'vertex_x', 'vertex_y'),
  new = c('poly_ID', 'x', 'y')
)

gpoly_cells = createGiottoPolygonsFromDfr(segmdfr = cellPoly_dt,
                                          name = 'cell',
                                          calc_centroids = TRUE)
gpoly_nucs = createGiottoPolygonsFromDfr(segmdfr = nucPoly_dt,
                                         name = 'nucleus',
                                         calc_centroids = TRUE)


## Create Giotto Object for entire slide -----------------------------------
gobject = createGiottoObjectSubcellular(
  gpoints = list(
    rna = gpoints_list$`Gene Expression`,
    blank_code = gpoints_list$`Blank Codeword`,
    neg_code = gpoints_list$`Negative Control Codeword`,
    neg_probe = gpoints_list$`Negative Control Probe`
  ),
  gpolygons = list(cell = gpoly_cells,
                   nucleus = gpoly_nucs),
  instructions = instrs
)

saveGiotto(gobject,dir = here(home.path, region),foldername = paste0(region,"_preprocessing_","giotto-object"), overwrite = TRUE)
# File will save as gobject.RDS so descriptive folder naming is essential
