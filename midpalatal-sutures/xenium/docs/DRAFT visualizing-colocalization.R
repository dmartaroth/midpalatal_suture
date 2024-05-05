# ## ######################################## ## #
#                 STEM CELLS                     #
# ## ######################################## ## #

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


# Region 7 ----------------------------------------------------------------


region <- "region-7" # replace region-x with correct region name for each
home.path <- here("midpalatal-sutures","xenium")
prepro.folder <- paste0(region,"_preprocessing_","giotto-object")

## Section e15c --------------------------------------------------------------------
section <- "e15c"

dir.create(section_folder <-
             here(home.path, region, section),
           recursive = TRUE)

dir.create(results_folder <-
             here(section_folder, "figs"),
           recursive = TRUE)

dir.create(spatial_folder <-
             here(results_folder, "spatial"),
           recursive = TRUE)

dir.create(output <- here(section_folder, "data-output"))


### Load preprocessed Giotto object -----------------------------------------
annot_foldername = paste0(region,"_",section,"_annotated_","giotto-object")

gobject <- loadGiotto(here(output, annot_foldername))


# Browse for stem cell genes
runApp("midpalatal-sutures/xenium/docs/spatplot-shinyapp.R")



# I would like to subset the cells by gene expression
# Best way to do this appears to be by annotating meta.data columns
# Where does Gli1 come up in the dataset?
# grep('Gli1',gobject@feat_ID$rna) # 126
# grep('Prrx1',gobject@feat_ID$rna) # 219
# grep('Ctsk',gobject@feat_ID$rna) # 70
# 
# 
# 
# # Assign Neg annotation to metadata if genes are 0 expression
# if(gobject@feat_ID$rna[126] <=0) {
#   gobject@cell_metadata$cell$rna$Gli1 <-
#     'Neg'
# } else{
#   gobject@cell_metadata$cell$rna$Gli1 <- 'Pos'
# }



grep('Gli1',gobject@feat_metadata$cell$rna@metaDT$feat_ID) # 44
# Assign Neg annotation to metadata if genes are 0 expression
if(gobject@feat_metadata$cell$rna@metaDT$mean_expr[44] <=0) {
  gobject@feat_metadata$cell$rna@metaDT$Gli1 <-
    'Neg'
} else{
  gobject@feat_metadata$cell$rna@metaDT$Gli1 <- 'Pos'
}




# # Subset only positive cells; change ident first
# Ident <- gobject@feat_metadata$cell$rna@metaDT$Gli1
# Gli1 <- subsetGiotto(gobject,cell_ids = "Gli1")
# 
# 



# # subset only the positive cells; change the ident first.
# dat <- SetAllIdent(dat, id = 'cd34')
# cd34 <- SubsetData(object = dat, ident.use = 'Pos')
# 
# 


plotfeats = c("Gli1")
spatInSituPlotPoints(gobject,
                     show_image = FALSE,
                     feats = list(plotfeats),
                     feats_color_code = feat_colors,
                     point_size = 1,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     
                     show_legend = TRUE,
                     polygon_alpha = 1,
                     polygon_color = 'bisque',
                     background_color = "white",
                     axis_text = 8,
                     axis_title = 9,
                     polygon_line_size = 0.01,
                     polygon_fill = "Gli1",
                     polygon_fill_as_factor = TRUE,
                     coord_fix_ratio = TRUE,
                     polygon_fill_code = colorcode,
                     plot_last = "points",
                     save_param = list(
                       save_name = paste0("13_", section, "_", region, "_Ctsk-in-MSC_spatplot"),
                       save_dir = spatial_folder),return_plot=T)

spatInSituPlotPoints(gobject,
                     show_image = FALSE,
                     feats = list(plotfeats),
                     feats_color_code = feat_colors,
                     point_size = 1,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     show_legend = TRUE,
                     polygon_alpha = 0.05,
                     polygon_color = 'bisque',
                     background_color = "white",
                     axis_text = 8,
                     axis_title = 9,
                     polygon_line_size = 0.01,
                     polygon_fill = 'cell_types',
                     polygon_fill_as_factor = TRUE,
                     coord_fix_ratio = TRUE,
                     polygon_fill_code = colorcode,
                     save_param = list(
                       save_name = paste0("15_", section, "_", region, "_Ctsk-in-all_spatplot"),
                       save_dir = spatial_folder),return_plot=T)
