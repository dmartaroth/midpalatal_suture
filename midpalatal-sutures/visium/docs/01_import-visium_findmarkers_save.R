# ## ######################################## ## #
#             IMPORT VISIUM AND SAVE             #
# ## ######################################## ## #

# Run this code to import and process Visium data

# Date: Sun Mar 31 15:23:08 2024 ------------------
# Updated by Daniela M. Roth

# Make directories for overall Visium data
library(here)
dir.create(here(home.path <-
                  here("midpalatal-sutures", "visium")), recursive = TRUE) 
dir.create(src <- here(home.path, "src"), recursive = TRUE) 



# e15.5 -------------------------------------------------------------------
plot_number <- 0
age <- "e15"

source(here::here("midpalatal-sutures/visium/docs/directories.R"))
source(here::here("midpalatal-sutures/visium/docs/packages.R"))
source(here::here("midpalatal-sutures/visium/docs/functions.R"))
source(here::here("midpalatal-sutures/visium/docs/themes.R"))


WT1_dir <- here::here(visium_folder,"JOP_WT_E15_5_WT","outs")
outs<- here::here(output,"WT_1")
WT1 <- Load10X_Spatial(data.dir = WT1_dir, slice="WT1")

p1 <- SpatialFeaturePlot(WT1, features = "nCount_Spatial")
convenient_save_plot(p1, "WT1_nCount_Spatial_prenormalization")


WT1<-NormalizeData(WT1, normalization.method = "RC", scale.factor = 10000)
WT1_coords<-GetTissueCoordinates(WT1)
p1 <- SpatialFeaturePlot(WT1, features = "nCount_Spatial")
convenient_save_plot(p1, "WT1_nCount_Spatial_postnormalization")


WT2_dir <- here::here(visium_folder, "Vx1_E15_WT","outs")
outs<- here::here(output,"WT_2")
WT2 <- Load10X_Spatial(data.dir = WT2_dir, slice="WT2")

p1 <- SpatialFeaturePlot(WT2, features = "nCount_Spatial")
convenient_save_plot(p1, "WT2_nCount_Spatial_prenormalization")
WT2<-NormalizeData(WT2, normalization.method = "RC", scale.factor = 10000)
WT2_coords<-GetTissueCoordinates(WT2)
p1 <- SpatialFeaturePlot(WT2, features = "nCount_Spatial")
convenient_save_plot(p1, "WT2_nCount_Spatial_postnormalization")



# Subset WT1 regions 1 and 2 ----------------------------------------------

WT1R1psPar<-subset(WT1, subset = wt1_imagerow>195 & wt1_imagerow<223&
                     wt1_imagecol>120&wt1_imagecol<190)

WT1R2psPar<-subset(WT1,subset = wt1_imagerow>200 & wt1_imagerow<260&
                     wt1_imagecol>300&wt1_imagecol<340)

p1 <- SpatialPlot(WT1R1psPar, features = "Col1a1",crop = TRUE, interactive=FALSE,pt.size.factor = 6,stroke=NA,max.cutoff = 500)+
  theme(aspect.ratio = 0.4,legend.position="right",legend.text = element_text(size = 8),legend.title = element_text(size=7,face = "italic"),legend.justification = "center",axis.text.y = element_text(face = "italic",size = 4,angle = 90),axis.text.x = element_text(face = "italic",size = 4,angle = 90))+ggplot2::scale_fill_gradientn(colors=c("#FFFFFF", "#436EEE", "#27408B"))

p2 <- SpatialPlot(WT1R2psPar, features = "Col1a1",crop = TRUE, interactive=FALSE,pt.size.factor = 6,stroke=NA,max.cutoff = 500)+
  theme(aspect.ratio = 1.2,legend.position="right",legend.text = element_text(size = 8),legend.title = element_text(size=7,face = "italic"),legend.justification = "center",axis.text.y = element_text(face = "italic",size = 4,angle = 90),axis.text.x = element_text(face = "italic",size = 4,angle = 90))+ggplot2::scale_fill_gradientn(colors=c("#FFFFFF", "#436EEE", "#27408B"))

plot <- plot_grid(p1,p2,ncol = 1)
convenient_save_plot(plot,"Col1a1_WT1_R1_R2", height = 8, width = 4)



# Cluster WT1R1 -----------------------------------------------------------

WT1R1ps <- ScaleData(WT1R1psPar)
WT1R1ps <- FindVariableFeatures(WT1R1ps)
WT1R1ps <- RunPCA(WT1R1ps)
WT1R1ps <- FindNeighbors(WT1R1ps)
WT1R1ps <- FindClusters(WT1R1ps,graph.name = "Spatial_snn",resolution=2.5,algorithm = 1)
table(WT1R1ps[[]]$seurat_clusters)
WT1R1ps <- RunUMAP(WT1R1ps,reduction="pca",dims=1:30)

(p1 <- SpatialDimPlot(WT1R1ps,
                     label=FALSE,
                     pt.size.factor = 6,
                     stroke=NA,
                     cols = visiumcolors)+
  theme(aspect.ratio = 0.4))

(p2 <- DimPlot(WT1R1ps,label=T,cols = visiumcolors,pt.size = 3))

plot <- plot_grid(p1,p2,ncol=1,rel_heights = c(1,2))
convenient_save_plot(plot,"spatdimplot_WT1R1", height = 6, width = 5)

# Find markers
# Define your minimum average log2FC threshold
min_avg_log2FC <- 0.25  # Adjust this value as needed
min_avg_expression <- 0.5  # Adjust this value as needed

# Find markers
all_markers_pct <- FindAllMarkers(WT1R1ps, verbose = TRUE) %>%
  Add_Pct_Diff()

# Filter markers based on minimum average log2FC and minimum average expression
filtered_markers <- all_markers_pct %>%
  group_by(cluster, gene) %>%
  filter(mean(avg_log2FC) >= min_avg_log2FC & mean(pct.1, pct.2) >= min_avg_expression) %>%
  select(gene) %>%
  unlist()

write.csv(all_markers_pct, file = here(output, "WT1R1_all_markers_pct.csv"))

# Extract the top N marker genes per cluster for plotting
top_5 <- Extract_Top_Markers(marker_dataframe = all_markers_pct, num_genes = 5, rank_by = "avg_log2FC")

# top 50 markers
top50_markers_pct <- all_markers_pct %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>% 
  top_n(n=50, wt = avg_log2FC)%>%
  arrange(cluster)

write.csv(top50_markers_pct, file = here(output, "WT1R1_top50_markers_pct.csv"))


(plot <- DotPlot(
  object = WT1R1ps,
  features =   top_5,
  scale.by = "size",
  dot.scale = 7,
  split.by = NULL,
  cluster.idents = FALSE,
) + scale_colour_gradient2(low = "dodgerblue",
                           mid = "floralwhite",
                           high = "red2") +  custom_dotplot_theme() +RotatedAxis())

convenient_save_plot(plot,"dotplot_top5_WT1R1", height = 4, width = 13)

# top gene spatial plots
cluster_number <- 3  # specify the cluster number you want to visualize
top_genes_number <- 40  # specify the number of top genes to consider
increment <- 4
top_genes_file <- here::here(output, "WT1R1_all_markers_pct.csv")  # replace with the actual path to your CSV file
dir <- here::here(results_folder)  # replace with your actual directory path

# Call the function with the Seurat object 'WT1R1ps' as the argument
generate_spatial_feature_plots(seurat_obj = WT1R1ps, cluster_number = cluster_number, 
                               top_genes_file = top_genes_file, top_genes_number = top_genes_number, 
                               increment = increment, dir = dir)



# Save object -------------------------------------------------------------

saveRDS(WT1R1ps, file = paste0(output, "/processed_WT1R1ps.Rds"))
