# ## ######################################## ## #
#                   ANNOTATION                   #
# ## ######################################## ## #

# Date: Wed Mar 20 10:15:19 2024 ------------------

# Rough documentation for e15c, to be updated for main scripts

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

dir.create(output <- here(section_folder, "data-output"))

dir.create(annot_folder <-
             here(results_folder, "annotation"),
           recursive = TRUE)

### Load preprocessed Giotto object -----------------------------------------
gobject <- loadGiotto(here(output, "e15c_giotto_objects"))

my_colors <- c("greenyellow", "slateblue3", "seagreen2","goldenrod1", "salmon2", "mediumpurple1", "seagreen", "lavender", "red3", "lightpink", "lightcyan", "magenta3", "lawngreen", "darkorchid4", "firebrick1", "blue", "orchid4", "salmon1", "seagreen3", "steelblue4", "royalblue2", "pink4", "cadetblue3", "red2", "chocolate1", "dodgerblue4", "darkolivegreen4", "magenta2", "skyblue2", "seagreen", "palevioletred2", "mediumpurple3", "aquamarine")

feat_colors <- c("firebrick1", "darkorchid4", "greenyellow","dodgerblue", "aquamarine", "magenta", "lawngreen", "palevioletred2", "gold", "cadetblue3", "lightcyan", "violetred1", "lawngreen", "darkorchid4", "firebrick1", "blue", "orchid4", "salmon1", "seagreen3", "steelblue4", "royalblue2", "pink4", "cadetblue3", "red2", "chocolate1", "dodgerblue4", "darkolivegreen4", "magenta2", "skyblue2", "seagreen", "palevioletred2", "mediumpurple3", "aquamarine")

colorcode = my_colors
featcolor = my_colors


spatInSituPlotPoints(gobject = gobject,
                     show_image = FALSE,
                     feats = list('rna' ="Alpl"),
                     feats_color_code = feat_colors,
                     point_size = 1.5,show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     show_legend = TRUE,
                     polygon_alpha = 0.1,
                     polygon_color = 'pink2',
                     background_color = "floralwhite",
                     axis_text = 8,
                     axis_title = 9,
                     polygon_line_size = 0.1,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = TRUE,
                     coord_fix_ratio = TRUE,
                     polygon_fill_code = colorcode,
                     save_param = list(
                       save_name = paste0("01_", section, "_", region, "_leidenclus"),
                       save_dir = annot_folder),return_plot=T)


### Cluster analysis --------------------------------------------------------
showClusterHeatmap(gobject = gobject, cluster_column = 'leiden_clus',
                   save_param = list(
                     save_name = paste0("02_", section, "_", 
                                        region, "_clusterheatmap"),
                     save_dir = annot_folder),return_plot=T)

# See cluster relationships in a dendrogram

p1 <- showClusterDendrogram(gobject = gobject,
                            h = 0.5, rotate = T, cluster_column = 'leiden_clus',
                            save_param = list(
                              save_name = paste0("03a_", section, "_", region,
                                                 "_cluster_dendrogram"),
                              save_dir = annot_folder),return_plot=T)

p2 <- plotUMAP(gobject = gobject,
               spat_unit = 'cell',
               cell_color = 'leiden_clus',
               cell_color_code=my_colors,
               show_legend = TRUE,
               point_size = 0.005,
               point_shape = 'no_border',
               save_param = list(
                 save_name = paste0("03b_", section, "_", region, "_UMAP"),
                 save_dir = annot_folder),return_plot=T)

dendrogram_UMAP <- plot_grid(p1,p2,ncol = 2)

filename <- paste0("03_", section, "_", region, "_combined_dendrogram_UMAP.png")
ggsave(file.path(annot_folder, filename), dendrogram_UMAP, width = 9, height = 3, dpi = 300)



# Find markers ------------------------------------------------------------
# Add code to save top markers to csv for each

markers <- findScranMarkers_one_vs_all(gobject = gobject,
                                       pval=0.01,logFC=0.25,
                                       cluster_column = "leiden_clus")

top10 <- markers %>% 
  group_by(cluster) %>% 
  slice_max(n=10,order_by = logFC)

top3 <- markers %>% 
  group_by(cluster) %>% 
  slice_max(n=3,order_by = logFC)

top1 <- markers %>% 
  group_by(cluster) %>% 
  slice_max(n=1,order_by = logFC)

violinPlot(gobject,expression_values = "scaled",
           feats = top1$feats,
           cluster_column = "leiden_clus",color_violin = "cluster",
           cluster_color_code = my_colors,
           save_param = list(
             save_name = paste0("04_", section, "_", region, "_topviolin"),
             base_width = 3,
             base_height = 7,
             save_dir = annot_folder),return_plot=F)

# Visualize top1 feats spatially
spatInSituPlotPoints(gobject,
                     show_image = FALSE,
                     feats = list('rna' =top1$feats),
                     feats_color_code = feat_colors,
                     point_size = 0.5,show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     show_legend = TRUE,
                     polygon_alpha = 0.1,
                     polygon_color = 'pink2',
                     background_color = "floralwhite",
                     axis_text = 8,
                     axis_title = 9,
                     polygon_line_size = 0.1,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = TRUE,
                     coord_fix_ratio = TRUE,
                     polygon_fill_code = colorcode,
                     save_param = list(
                       save_name = paste0("05_", section, "_", region, "_top1_spatplot"),
                       save_dir = annot_folder),return_plot=F)

# Scran markers
scran_markers_subclusters <- findMarkers_one_vs_all(gobject = gobject,
                                                        method = "scran",
                                                        expression_values = "normalized",
                                                        cluster_column = "leiden_clus")


top3 <- scran_markers_subclusters %>% 
  group_by(cluster) %>% 
  slice_max(n=3,order_by = logFC)

p1 <- plotMetaDataHeatmap(gobject = gobject,
                          selected_feats = top3$feats,
                          metadata_cols = c("leiden_clus"),
                          save_param = list(
                            save_name = paste0("06a_", section, "_", region, "_top3_metadataheatmap"),
                            save_dir = annot_folder),return_plot=F)
p2 <- plotUMAP(gobject = gobject,
               spat_unit = 'cell',
               cell_color = 'leiden_clus',
               cell_color_code=my_colors,
               show_legend = TRUE,
               point_size = 0.005,
               point_shape = 'no_border',
               save_param = list(
                 save_name = paste0("06b_", section, "_", region, "_UMAP"),
                 save_dir = annot_folder),return_plot=F)

top3heatmap_UMAP <- plot_grid(p1,p2,ncol = 2)

filename <- paste0("06_", section, "_", region, "_combined_top3_heatmap_UMAP.png")
ggsave(file.path(annot_folder, filename), top3heatmap_UMAP, width = 7, height = 4, dpi = 300)



### Annotate clusters -------------------------------------------------------

# show leiden clustering results
cell_metadata <- pDataDT(gobject)
cell_metadata[['leiden_clus']]


p1 <- spatInSituPlotPoints(gobject,
                     show_image = FALSE,
                     feats = NULL,
                     point_size = 0.005,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_alpha = 1,
                     polygon_color = 'black',
                     polygon_line_size = 0.01,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = TRUE,
                     coord_fix_ratio = TRUE,
                     polygon_fill_code = colorcode,
                     save_param = list(
                       save_name = paste0("08a_", section, "_", region, "_leidenclus_poly"),
                       save_dir = annot_folder),return_plot=F)

p2 <- plotUMAP(gobject = gobject,
               spat_unit = 'cell',
               cell_color = 'leiden_clus',
               cell_color_code=my_colors,
               show_legend = TRUE,
               point_size = 0.005,
               point_shape = 'no_border',
               save_param = list(
                 save_name = paste0("08b_", section, "_", region, "_UMAP"),
                 base_width = 3,
                 base_height = 3,
                 save_dir = annot_folder),return_plot=F)

spatplot_UMAP <- plot_grid(p1,p2,ncol = 2)

filename <- paste0("08_", section, "_", region, "_combined_spatplot_UMAP.png")
ggsave(file.path(annot_folder, filename), spatplot_UMAP, width = 9, height = 3, dpi = 300)

# Create annotation file and script as in parse analysis scripts
# For now,
# cluster - prediction - annotation - color
# 1 - lateral mesenchyme - mes.2 - "greenyellow"
# 2 - osteogenic front - mes.3 - "goldenrod1"
# 3 - epithelium - epith.1 - "#9F79EE"
# 4 - midline mesenchyme - mes.4 - "dodgerblue"
# 5 - vasculature - vasc - "#006400"
# 6 - bone - osteo - "red"
# 7 - mesenchyme above oral epithelium - mes.5 - "salmon2"
# 8 - ciliated nasal epithelium - epith.2 - "orchid"
# 9 - rare mesenchyme - mes.1 - "darkolivegreen3"
colorcode = c("#9F79EE","orchid","darkolivegreen3","greenyellow","goldenrod1","dodgerblue","salmon2","red","#006400")

# create vector with cell type names as names of the vector
cluster_cell_types <- c("mes.2","mes.3","epith.1","mes.4","vasc","osteo","mes.5","epith.2","mes.1")
names(cluster_cell_types) = 1:9

# convert cluster results into annotations and add to cell metadata
gobject <- annotateGiotto(gobject,annotation_vector = cluster_cell_types,cluster_column = "leiden_clus",name = "cell_types")

# inspect new annotation column
pDataDT(gobject)

# visualize annotation results
# annotation name is cell_types as provided in the previous command
spatDimPlot(gobject,
            cell_color = "cell_types",
            dim_show_legend=TRUE,
            cell_color_code=colorcode,
            spat_point_size = 0.5, dim_point_size = 0.5,
            dim_label_size = 2,
            spat_label_size = 2,
            spat_point_alpha = 0.7,
            axis_title = 5,
            dim_point_shape = "no_border",
            spat_point_shape = "no_border",
            save_param = list(
              save_name = paste0("07_", section, "_", region, "_spatdimplot_annotated"),
              save_dir = annot_folder),return_plot=F)

spatInSituPlotPoints(gobject,
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
                       save_name = paste0("09_", section, "_", region, "_celltypes_poly"),
                       save_dir = annot_folder),return_plot=F)


### Save annotated gobject --------------------------------------------------
saveGiotto(gobject,dir = here(output),foldername = paste0(region,"_",section,"_annotated_","giotto-object"), overwrite = TRUE)
# File will save as gobject.RDS so descriptive folder naming is essential
# This is for left palatal shelf - I need to do this for left, right, and whole

