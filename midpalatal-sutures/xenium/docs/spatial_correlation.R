# ## ######################################## ## #
#            SPATIAL METAFEAT ANALYSIS           #
# ## ######################################## ## #

# Updated by: Daniela M. Roth
# Date: Wed Mar 20 10:04:38 2024 ------------------


# needs to be situated properly in scripts
# this is rough analysis to expedite hi-plex experiment targets
# manually moved giotto objects folder from script 02 to data-output folder for
# region/section folder for now until previous script updated for correct 
# save location


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

colorcode = c("#9F79EE","orchid","darkolivegreen3","greenyellow","goldenrod1","dodgerblue","salmon2","red","#006400")

# Visualize clusters
p1 <- spatInSituPlotPoints(gobject,
                           show_image = FALSE,
                           feats = NULL,
                           point_size = 0.005,
                           show_polygon = TRUE,
                           polygon_feat_type = 'cell',
                           polygon_alpha = 1,
                           polygon_color = 'black',
                           polygon_line_size = 0.01,
                           polygon_fill = 'cell_types',
                           polygon_fill_as_factor = TRUE,
                           coord_fix_ratio = TRUE,
                           polygon_fill_code = colorcode,
                           save_param = list(
                             save_name = paste0("01a_", section, "_", region, "_celltypes_poly"),
                             save_dir = spatial_folder),return_plot=F)


p2 <- plotUMAP(gobject = gobject,
               spat_unit = 'cell',
               cell_color = 'cell_types',
               cell_color_code=colorcode,
               show_legend = TRUE,
               point_size = 0.005,
               point_shape = 'no_border',
               save_param = list(
                 save_name = paste0("01b_", section, "_", region, "_UMAP"),
                 base_width = 3,
                 base_height = 3,
                 save_dir = spatial_folder),return_plot=F)

polys_UMAP <- plot_grid(p1,p2,ncol = 2)
filename <- paste0("01_", section, "_", region, "_combined_polys_UMAP.png")
ggsave(file.path(here(spatial_folder, filename)), plot = polys_UMAP, width = 9, height = 3, dpi = 300)


### Delaunay triangulation --------------------------------------------------
plotStatDelaunayNetwork(gobject = gobject, maximum_distance = 15,
                        save_param = list(
                          save_name = paste0("02_", section, "_", region, "_stat_Delaunay"),
                          save_dir = spatial_folder),return_plot=F)

gobject = createSpatialNetwork(gobject = gobject,
                               minimum_k = 6,
                               maximum_distance_delaunay = 12)

spatPlot(gobject = gobject,
         show_network = T,
         point_shape = "no_border",
         network_color = "lightpink",
         spatial_network_name = "Delaunay_network",
         point_size=0.5,
         cell_color_code = colorcode,
         cell_color = "cell_types",
         coord_fix_ratio = 1,
         save_param = list(
           save_name = paste0("03_", section, "_", region, "_Delaunay"),
           save_dir = spatial_folder),return_plot=F)+theme(legend.position = "bottom")


### Rank spatial feats ------------------------------------------------------
rank_spatialfeats = binSpect(gobject = gobject, bin_method = 'rank')

spatFeatPlot2D(gobject, expression_values = 'normalized', 
               feats = rank_spatialfeats[1:18]$feats,
               point_shape = 'border', point_border_stroke = 0.1,
               show_network = F, network_color = 'lightpink', point_size = 2,
               cow_n_col = 3,
               save_param = list(
                 save_name = paste0("04_", section, "_", region, "_rankedfeats"),
                 save_dir = spatial_folder),return_plot=F)

ext_spatial_feats = rank_spatialfeats[1:200]$feats
spat_cor_netw_DT = detectSpatialCorFeats(gobject,
                                         method = "network",
                                         spatial_network_name = "Delaunay_network",
                                         subset_feats = ext_spatial_feats)

spat_cor_netw_DT = clusterSpatialCorFeats(spat_cor_netw_DT,
                                          name = 'spat_netw_clus',
                                          k=6)
heatmSpatialCorFeats(gobject = gobject,
                     spatCorObject = spat_cor_netw_DT,
                     use_clus_name = 'spat_netw_clus',
                     save_param = list(
                       save_name = paste0("05_", section, "_", region, "_heatmSpatialCorFeats"),
                       save_dir = spatial_folder),return_plot=F)


netw_ranks = rankSpatialCorGroups(gobject, spatCorObject = spat_cor_netw_DT, 
                                  use_clus_name = 'spat_netw_clus',
                                  show_plot = TRUE,
                                  save_param = list(
                                    save_name = paste0("06_", section, "_", region, "_rankSpatialCorGroups"),
                                    save_dir = spatial_folder),return_plot=F)

top_netw_spat_cluster = showSpatialCorFeats(spat_cor_netw_DT, 
                                            use_clus_name = 'spat_netw_clus', 
                                            show_top_feats = 3)
cluster_feats = top_netw_spat_cluster$clus; names(cluster_feats) = top_netw_spat_cluster$feat_ID

gobject = createMetafeats(gobject, feat_clusters = cluster_feats, name = 'cluster_metafeat')


### Plot network rank patterns ----------------------------------------------
spatCellPlot(gobject,
             point_border_stroke=0.1,
             spat_enr_names = 'cluster_metafeat',point_shape="no_border",
             cell_annotation_values = netw_ranks$clusters,cell_color_gradient=c("cyan3","cornsilk","red2"),
             point_size = 0.5, cow_n_col = 1,
             save_param = list(
               save_name = paste0("07_", section, "_", region, "_netwrankpatterns"),
               save_dir = spatial_folder),return_plot=F )


spatcorfeats <- showSpatialCorFeats(spat_cor_netw_DT, 
                                    use_clus_name = 'spat_netw_clus')
top_netw_spat_cluster = showSpatialCorFeats(spat_cor_netw_DT, 
                                            use_clus_name = 'spat_netw_clus',
                                            selected_clusters = 1:6, 
                                            show_top_feats = 1)

write.csv(top_netw_spat_cluster, file = here::here(section_folder, "data-output","top_netw_spat_cluster.csv"))

# Using top hits from this file, generate spatial expression of genes with shiny app
runApp("midpalatal-sutures/xenium/docs/spatplot-shinyapp.R")

# Plot feats chosen for Hi-Plex ACD assay
plotfeats <- c("Pax9","Alx4","Chodl","Col12a1","Dcn","Gsc","Sox9","Thbs1","Tnn")
feat_colors <- c("seagreen3","goldenrod1","dodgerblue","red2","darkorchid4","dodgerblue","white","blue","magenta","lightgrey","red3")

spatInSituPlotPoints(gobject,
                     show_image = FALSE,
                     feats = list(plotfeats),
                     feats_color_code = feat_colors,
                     point_size = 1,
                     show_polygon = FALSE,
                     polygon_feat_type = 'cell',
                     show_legend = TRUE,
                     polygon_alpha = 0.5,
                     polygon_color = 'bisque',
                     background_color = "black",
                     axis_text = 8,
                     axis_title = 9,
                     polygon_line_size = 0.01,
                     polygon_fill = 'cell_types',
                     polygon_fill_as_factor = TRUE,
                     coord_fix_ratio = TRUE,
                     polygon_fill_code = colorcode,
                     save_param = list(
                       save_name = paste0("08_", section, "_", region, "_hiplex_9_spatplot"),
                       save_dir = spatial_folder),return_plot=F)

