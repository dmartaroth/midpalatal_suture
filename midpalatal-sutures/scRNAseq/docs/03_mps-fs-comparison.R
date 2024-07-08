# ## ######################################## ## #
#                MPS - FS COMPARISON             #
# ## ######################################## ## #

# Comparison of e15.5 xenium, e15.5 mid-palatal suture scRNA-seq, and e18
## frontal suture scRNA-seq

# Update scripts with here() package for reproducible saving and compatibility
## with other Xenium import scripts - currently redundant


# Load Xenium data --------------------------------------------------------

#load libraries
library(Giotto)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(scran)

# 1. ** SET WORKING DIRECTORY WHERE PROJECT OUPUTS WILL SAVE TO **
results_folder = 'C:/Users/dmrot/Documents/midpalatal_suture/scRNAseq/comparison/data-output'

# 2. set giotto python path
# set python path to your preferred python version path
# set python path to NULL if you want to automatically install (only the 1st time) and use the giotto miniconda environment
python_path = NULL
if(is.null(python_path)) {
  installGiottoEnvironment()
}

# 3. Create Giotto instructions
# Directly saving plots to the working directory without rendering them in the editor saves time.
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  return_plot = TRUE)

# load xenium
xenium_gobj <- loadGiotto("C:/Users/dmrot/Documents/midpalatal_suture/Visium/docs/RMarkdown/annotated_e15c_ps")

e15_xenium_Lshelf = subsetGiottoLocs(xenium_gobj,
                                     x_min = 5150, x_max = 5800,
                                     y_min = 3150, y_max = 3500)


# Define colors -----------------------------------------------------------

xenium_colors <- c("#FFE4E1", "#BF3EFF", "#9F79EE","dodgerblue","greenyellow","orange","darkolivegreen4","navyblue","turquoise" ,"red", "forestgreen", "navyblue", "greenyellow", "orange", "#006400")

leiden_colors <- c("greenyellow","#9F79EE","dodgerblue","orange","darkolivegreen4","#BF3EFF","forestgreen","red","#FFE4E1","turquoise","navyblue")

feat_colors <- c("firebrick1", "darkorchid4", "greenyellow","dodgerblue", "aquamarine", "magenta", "lawngreen", "palevioletred2", "gold", "cadetblue3", "lightcyan", "violetred1", "lawngreen", "darkorchid4", "firebrick1", "blue", "orchid4", "salmon1", "seagreen3", "steelblue4", "royalblue2", "pink4", "cadetblue3", "red2", "chocolate1", "dodgerblue4", "darkolivegreen4", "magenta2", "skyblue2", "seagreen", "palevioletred2", "mediumpurple3", "aquamarine")

e15sc_colors <- c("darkolivegreen2", "orange", "purple", "lightcoral", "skyblue","maroon2","slateblue","gold","dodgerblue3","plum1","darkseagreen3","pink","violetred4","tomato1","orchid3","darkolivegreen3")

colorcode = xenium_colors
featcolor = feat_colors


# Examine Xenium subset ---------------------------------------------------

spatDimPlot(e15_xenium_Lshelf,
            cell_color = "cell_types2",
            dim_show_legend=TRUE,
            cell_color_code=xenium_colors,
            spat_point_size = 0.5, dim_point_size = 0.5,
            dim_label_size = 2,
            spat_label_size = 2,
            spat_point_alpha = 0.7,
            axis_title = 5,
            dim_point_shape = "no_border",
            spat_point_shape = "no_border",
            save_para = list(
              save_name = '2_spatdimplot_annotated_e15_xeniumLshelf'), return_plot = T)


# Identify top markers for e15 Xenium clusters ----------------------------
xenium_scran_markers <-
  findMarkers_one_vs_all(
    gobject = xenium_gobj,
    method = "scran",
    expression_values = "normalized",
    cluster_column = "leiden_clus"
  )


top10 <- xenium_scran_markers %>% 
  group_by(cluster) %>% 
  slice_max(n=10,order_by = logFC)

top3 <- xenium_scran_markers %>% 
  group_by(cluster) %>% 
  slice_max(n=3,order_by = logFC)

plotMetaDataHeatmap(gobject = xenium_gobj,
                    selected_feats = top3$feats,
                    metadata_cols = c("cell_types2"),
                    custom_cluster_order = c("mes.1","mes.2","mes.3","mes.4","mes.5","mes.6","osteo.2","vasc","chondro","ciliated","neu","epith.1"),
                    save_para = list(
                      save_name = '3_e15_xenium_top3_heatmap'), return_plot = T)+coord_flip()

# Save top10 markers
# top10.markers <-
  # write.csv(top10, file = "C:/Users/dmrot/Documents/midpalatal_suture/scRNAseq/comparison/data-output/top10.e15xenium.csv")



# Annotation --------------------------------------------------------------

# check which clusters are which cell_types
p1 <- spatDimPlot(e15_xenium_Lshelf,
                  cell_color = "cell_types2",
                  dim_show_legend=TRUE,
                  cell_color_code=xenium_colors,
                  spat_point_size = 0.5, dim_point_size = 0.5,
                  dim_label_size = 2,
                  spat_label_size = 2,
                  spat_point_alpha = 0.7,
                  axis_title = 5,
                  dim_point_shape = "no_border",
                  spat_point_shape = "no_border",
                  save_para = list(
                    save_name = '2_spatdimplot_annotated_e15_xeniumLshelf'), return_plot = T)

p2 <- spatDimPlot(e15_xenium_Lshelf,
                  cell_color = "leiden_clus",
                  dim_show_legend=TRUE,
                  cell_color_code=leiden_colors,
                  spat_point_size = 0.5, dim_point_size = 0.5,
                  dim_label_size = 2,
                  spat_label_size = 2,
                  spat_point_alpha = 0.7,
                  axis_title = 5,
                  dim_point_shape = "no_border",
                  spat_point_shape = "no_border",
                  save_para = list(
                    save_name = '2c_spatdimplot_annotated_e15_xeniumLshelf'), return_plot = T)

plot_grid(p1,p2,ncol=2)


chondro.top10 <- c("Col11a1","Sox9","Mgp","Col16a1","Igf2","Acan","Col2a1","Postn","Bgn","Col4a1")
ciliated.top10 <- c("Lrrc23","Foxj1","Deup1","Crocc2","Mrc1","Dynlrb2","Capsl","Top2a","4833427G06Rik","S100a6")
epith.1.top10 <- c("Krt14","Krt5","Hsp90ab1","Sox11","Mapk13","Net1","Notch1","Hspa9","Hmgn2","Tfap2a")
mes.1.top10 <- c("Dlk1","Acta1","Mki67","Lmnb1","Col4a1","Top2a","Atf5","Ranbp1","Nasp","Brd2")
mes.2.top10 <- c("Sfrp2","Postn","Dlk1","Eln","Col3a1","Itm2a","Fn1","Fstl1","Igf1","Pdgfra")
mes.3.top10 <- c("Ptn","Col12a1","Col16a1","Col1a1","Tnn","Col3a1","Vim","Igf2","Sfrp2","Ogn")
mes.4.top10 <- c("Ptch1","Frzb","Col3a1","Inhba","Pdgfra","Postn","Slit3","Net1","Fstl1","Gli1")
neu.top10 <- c("Atf5","Ncam1","Hsp90ab1","Tuba1b","Rock1","Tubb3","Cd200","Fgfr3","Mmp15","Calcrl")
osteo.2.top10 <- c("Ibsp","Sparc","Serpinh1","Creb3l1","Alpl","Col1a1","Sp7","Cd63","Hspa5","Tgfb1")
vasc.top10 <- c("Cdh5","Col4a1","Vim","Sparcl1","Sparc","Pecam1","Flt1","Esam","Notch1","Hspg2")



# Load e15.5 mps scRNA-seq ------------------------------------------------


library(Seurat)
library(dittoSeq)
e15mps_sc <-
  readRDS(
    "C:/Users/dmrot/Documents/midpalatal_suture/scRNAseq/E15_mps/data-output/mps_annotated.Rds"
  )

dittoDimPlot(
  e15mps_sc,
  "ident",
  do.label = TRUE,
  labels.repel = TRUE,
  labels.size = 2,
  legend.size = 1,
  size = 0.5,
  labels.highlight = FALSE,
  opacity = 0.9,
  do.ellipse = TRUE,
  color.panel = c(
    "darkolivegreen2",
    "orange",
    "purple",
    "lightcoral",
    "skyblue",
    "maroon2",
    "slateblue",
    "gold",
    "dodgerblue3",
    "plum1",
    "darkseagreen3",
    "pink",
    "violetred4",
    "tomato1",
    "orchid3",
    "darkolivegreen3"
  )
) + labs(title = 'E15.5 mid-palatal suture') + theme(
  title = element_text(size = 5),
  legend.text = element_text(size = 5),
  axis.title = element_text(size = 4)
)


# Visualize top Xenium markers in scRNAseq --------------------------------

# mes.3 cluster
p1 <- FeaturePlot(e15mps_sc,features=c(mes.3.top10),pt.size = 0.5,ncol = 3)

title <- ggdraw()+draw_label("top10 mes.3 xenium cluster e15scRNAseq feature plots",x=0,hjust=0)+theme(plot.margin = margin(0,0,0,7))

plot_grid(title,p1,ncol=1,rel_heights=c(0.1,1))


# osteo.2 cluster
p1 <- FeaturePlot(e15mps_sc,features=c(osteo.2.top10),pt.size = 0.5,ncol = 3)

title <- ggdraw()+draw_label("top10 osteo.2 xenium cluster e15scRNAseq feature plots",x=0,hjust=0)+theme(plot.margin = margin(0,0,0,7))

plot_grid(title,p1,ncol=1,rel_heights=c(0.1,1))

# mes.2 cluster
p1 <- FeaturePlot(e15mps_sc,features=c(mes.2.top10),pt.size = 0.5,ncol = 3)

title <- ggdraw()+draw_label("top10 mes.2 xenium cluster e15scRNAseq feature plots",x=0,hjust=0)+theme(plot.margin = margin(0,0,0,7))

plot_grid(title,p1,ncol=1,rel_heights=c(0.1,1))

# Visualize with dotplot

# Mesenchymal clusters
library(scCustomize)
p1 <- DotPlot(
  object = e15mps_sc, features =   mes.1.top10,scale.by = "size"
) + scale_colour_gradient2(low = "red", mid = "white", high ="blue")+RotatedAxis() +ggtitle("E15.5 sc mid-palatal suture mes.1 xenium")+
  theme(plot.title = element_text(color="blue",size=9,face="bold"),plot.subtitle=element_text(color="red",face="italic"),text = element_text(size = 4))

p2 <- DotPlot(
  object = e15mps_sc, features =   mes.2.top10,scale.by = "size"
) + scale_colour_gradient2(low = "red", mid = "white", high ="blue")+RotatedAxis() +ggtitle("E15.5 sc mid-palatal suture mes.2 xenium")+
  theme(plot.title = element_text(color="blue",size=9,face="bold"),plot.subtitle=element_text(color="red",face="italic"),text = element_text(size = 4))

p3 <- DotPlot(
  object = e15mps_sc, features =   mes.3.top10,scale.by = "size"
) + scale_colour_gradient2(low = "red", mid = "white", high ="blue")+RotatedAxis() +ggtitle("E15.5 sc mid-palatal suture mes.3 xenium")+
  theme(plot.title = element_text(color="blue",size=9,face="bold"),plot.subtitle=element_text(color="red",face="italic"),text = element_text(size = 4))

p4 <- DotPlot(
  object = e15mps_sc, features =   mes.4.top10,scale.by = "size"
) + scale_colour_gradient2(low = "red", mid = "white", high ="blue")+RotatedAxis() +ggtitle("E15.5 sc mid-palatal suture mes.4 xenium")+
  theme(plot.title = element_text(color="blue",size=9,face="bold"),plot.subtitle=element_text(color="red",face="italic"),text = element_text(size = 4))

plot_grid(p1,p2,p3,p4,ncol=2)

# Epithelial, ciliated, vascular, neuronal clusters
p1 <- DotPlot(
  object = e15mps_sc, features =   epith.1.top10,scale.by = "size"
) + scale_colour_gradient2(low = "red", mid = "white", high ="blue")+RotatedAxis() +ggtitle("E15.5 sc mid-palatal suture epith.1 xenium")+
  theme(plot.title = element_text(color="blue",size=9,face="bold"),plot.subtitle=element_text(color="red",face="italic"),text = element_text(size = 4))

p2 <- DotPlot(
  object = e15mps_sc, features =   ciliated.top10,scale.by = "size"
) + scale_colour_gradient2(low = "red", mid = "white", high ="blue")+RotatedAxis() +ggtitle("E15.5 sc mid-palatal suture ciliated xenium")+
  theme(plot.title = element_text(color="blue",size=9,face="bold"),plot.subtitle=element_text(color="red",face="italic"),text = element_text(size = 4))

p3 <- DotPlot(
  object = e15mps_sc, features =   vasc.top10,scale.by = "size"
) + scale_colour_gradient2(low = "red", mid = "white", high ="blue")+RotatedAxis() +ggtitle("E15.5 sc mid-palatal suture vasc xenium")+
  theme(plot.title = element_text(color="blue",size=9,face="bold"),plot.subtitle=element_text(color="red",face="italic"),text = element_text(size = 4))

p4 <- DotPlot(
  object = e15mps_sc, features =   neu.top10,scale.by = "size"
) + scale_colour_gradient2(low = "red", mid = "white", high ="blue")+RotatedAxis() +ggtitle("E15.5 sc mid-palatal suture neu xenium")+
  theme(plot.title = element_text(color="blue",size=9,face="bold"),plot.subtitle=element_text(color="red",face="italic"),text = element_text(size = 4))

plot_grid(p1,p2,p3,p4,ncol=2)

# Osteochondro clusters
p1 <- DotPlot(
  object = e15mps_sc, features =   osteo.2.top10,scale.by = "size"
) + scale_colour_gradient2(low = "red", mid = "white", high ="blue")+RotatedAxis() +ggtitle("E15.5 sc mid-palatal suture osteo.2 xenium")+
  theme(plot.title = element_text(color="blue",size=9,face="bold"),plot.subtitle=element_text(color="red",face="italic"),text = element_text(size = 4))

p2 <- DotPlot(
  object = e15mps_sc, features =   chondro.top10,scale.by = "size"
) + scale_colour_gradient2(low = "red", mid = "white", high ="blue")+RotatedAxis() +ggtitle("E15.5 sc mid-palatal suture chondro xenium")+
  theme(plot.title = element_text(color="blue",size=9,face="bold"),plot.subtitle=element_text(color="red",face="italic"),text = element_text(size = 4))

plot_grid(p1,p2,ncol=2)

# Compare xenium with scRNAseq mps dimplot
p1 <- dimPlot2D(e15_xenium_Lshelf,cell_color = "cell_types2",cell_color_code = xenium_colors,point_shape = "no_border",point_size = 0.5)

p2 <- DimPlot(e15mps_sc,cols = e15sc_colors,reduction="umap")

plot_grid(p1,p2,ncol=2)


# Load e18 frontal suture dataset -----------------------------------------

library(Seurat)
library(dittoSeq)
e18fs_sc <-
  readRDS(
    "C:/Users/dmrot/Documents/midpalatal_suture/scRNAseq/E18_fs/data-output/fs_annotated.Rds"
  )

dittoDimPlot(
  e18fs_sc,
  "ident",
  do.label = TRUE,
  labels.repel = TRUE,
  labels.size = 2,
  legend.size = 1,
  size = 0.5,
  labels.highlight = FALSE,
  opacity = 0.9,
  do.ellipse = TRUE,
  color.panel = c(
    "darkolivegreen2",
    "orange",
    "purple",
    "lightcoral",
    "skyblue",
    "maroon2",
    "slateblue",
    "gold",
    "dodgerblue3",
    "plum1",
    "darkseagreen3",
    "pink",
    "violetred4",
    "tomato1",
    "orchid3",
    "darkolivegreen3"
  )
) + labs(title = 'e18 frontal suture') + theme(
  title = element_text(size = 5),
  legend.text = element_text(size = 5),
  axis.title = element_text(size = 4)
)


# Plot xenium top markers on fs scRNA-seq data ----------------------------

# mes.3 cluster
p1 <- FeaturePlot(e18fs_sc,features=c(mes.3.top10),pt.size = 0.5,ncol = 3)

title <- ggdraw()+draw_label("top10 mes.3 xenium cluster e18fs scRNAseq feature plots",x=0,hjust=0)+theme(plot.margin = margin(0,0,0,7))

plot_grid(title,p1,ncol=1,rel_heights=c(0.1,1))

# osteo.2 cluster
p1 <- FeaturePlot(e18fs_sc,features=c(osteo.2.top10),pt.size = 0.5,ncol = 3)

title <- ggdraw()+draw_label("top10 osteo.2 xenium cluster e18fs scRNAseq feature plots",x=0,hjust=0)+theme(plot.margin = margin(0,0,0,7))

plot_grid(title,p1,ncol=1,rel_heights=c(0.1,1))

# mes.2 cluster
p1 <- FeaturePlot(e18fs_sc,features=c(mes.2.top10),pt.size = 0.5,ncol = 3)

title <- ggdraw()+draw_label("top10 mes.2 xenium cluster e18fs scRNAseq feature plots",x=0,hjust=0)+theme(plot.margin = margin(0,0,0,7))

plot_grid(title,p1,ncol=1,rel_heights=c(0.1,1))

# Visualize with dot plots
library(scCustomize)
p1 <- DotPlot(
  object = e18fs_sc, features =   mes.1.top10,scale.by = "size"
) + scale_colour_gradient2(low = "red", mid = "white", high ="blue")+RotatedAxis() +ggtitle("E18.5 sc frontal suture mes.1 xenium")+
  theme(plot.title = element_text(color="blue",size=9,face="bold"),plot.subtitle=element_text(color="red",face="italic"),text = element_text(size = 4))

p2 <- DotPlot(
  object = e18fs_sc, features =   mes.2.top10,scale.by = "size"
) + scale_colour_gradient2(low = "red", mid = "white", high ="blue")+RotatedAxis() +ggtitle("E18.5 sc frontal suture mes.2 xenium")+
  theme(plot.title = element_text(color="blue",size=9,face="bold"),plot.subtitle=element_text(color="red",face="italic"),text = element_text(size = 4))

p3 <- DotPlot(
  object = e18fs_sc, features =   mes.3.top10,scale.by = "size"
) + scale_colour_gradient2(low = "red", mid = "white", high ="blue")+RotatedAxis() +ggtitle("E18.5 sc frontal suture mes.3 xenium")+
  theme(plot.title = element_text(color="blue",size=9,face="bold"),plot.subtitle=element_text(color="red",face="italic"),text = element_text(size = 4))

p4 <- DotPlot(
  object = e18fs_sc, features =   mes.4.top10,scale.by = "size"
) + scale_colour_gradient2(low = "red", mid = "white", high ="blue")+RotatedAxis() +ggtitle("E18.5 sc frontal suture mes.4 xenium")+
  theme(plot.title = element_text(color="blue",size=9,face="bold"),plot.subtitle=element_text(color="red",face="italic"),text = element_text(size = 4))

plot_grid(p1,p2,p3,p4,ncol=2)


# Plot top genes (dotplot) ------------------------------------------------

#reorder levels on mpssc dataset to somewhat follow xenium
# this didn't work, debug later
#my_levels <- c("mes.1","mes.2","mes.3","musc","vasc","cilia","neu.1","neu.2","neu.3","neu.4","neu.5","epith","eryth")

#factor(Idents(e15mps_sc),levels = my_levels)

#Idents(e15mps_sc) <- factor(Idents(e15mps_sc),levels = my_levels)


markers_to_plot <- top3$feats
markers_to_plot <- unique(markers_to_plot)

#top3 is a lot of dots for this plot. Let's do top2 for these
top2 <- xenium_scran_markers %>% 
  group_by(cluster) %>% 
  slice_max(n=2,order_by = logFC)

markers_to_plot <- top2$feats
markers_to_plot <- unique(markers_to_plot)

p1 <- DotPlot(object = e15mps_sc, features =   markers_to_plot,scale.by = "size" )+
  scale_colour_gradient2(low="dodgerblue",mid="white",high = "red2")+
  theme(legend.position="bottom",
        legend.text = element_text(size = 6),
        legend.title = element_text(size=6),
        legend.justification = "left",
        axis.text = element_text(size = 8),
        axis.text.x = element_text(face = "italic",angle=45,hjust = 1),
        axis.title = element_text(size = 8),
        axis.title.x = element_text(hjust = 0.5))  

p2 <- DotPlot(  object = e18fs_sc, features =   markers_to_plot,scale.by = "size"
)+  scale_colour_gradient2(low="dodgerblue",mid="white",high = "red2")+
  theme(legend.position="bottom",
        legend.text = element_text(size = 6),
        legend.title = element_text(size=6),
        legend.justification = "left",
        axis.text = element_text(size = 8),
        axis.text.x = element_text(face = "italic",angle=45,hjust = 1),
        axis.title = element_text(size = 8),
        axis.title.x = element_text(hjust = 0.5))

plot_grid(p1,p2,ncol=1)

xenium_colors <- c("#FFE4E1", "#BF3EFF", "#9F79EE","dodgerblue","greenyellow","orange","darkolivegreen4","navyblue","lightblue","turquoise" ,"red", "forestgreen", "navyblue", "greenyellow", "orange", "#006400")

dimPlot2D(xenium_gobj,feat_type = 'rna',cell_color = "cell_types2",point_size = 0.1,
          label_size = 2,axis_title = 5,
          point_shape = "no_border",cell_color_code=xenium_colors)

e15sc_colors <- c("#BF3EFF","#9F79EE","red4","#FFE4B5", "#FAFAD2","dodgerblue","greenyellow","orange","darkolivegreen4","navyblue","turquoise","turquoise3","aquamarine","aquamarine3","aquamarine4","#006400")

e18fs_colors <- c("red4","#FFE4B5", "#FAFAD2","bisque3","beige","dodgerblue","greenyellow","orange","darkolivegreen4","gold","lightgrey","navyblue","turquoise","#006400","green4")


# Visualize annotated clusters  ------------------------------------

DimPlot_scCustom(
  e15mps_sc,
  colors_use = e15sc_colors,
  reduction = "umap",
  pt.size = 0.1,
  shuffle = TRUE,
  label.size = 2
) + theme(
  text = element_text(size = 6),
  axis.title = element_text(size = 4),
  axis.text = element_text(size = 8)
)

DimPlot_scCustom(
  e18fs_sc,
  colors_use = e18fs_colors,
  reduction = "umap",
  pt.size = 0.1,
  shuffle = TRUE,
  label.size = 2
) + theme(
  text = element_text(size = 6),
  axis.title = element_text(size = 4),
  axis.text = element_text(size = 8)
)

xenium_colors <-
  c(
    "#FFE4E1",
    "#BF3EFF",
    "#9F79EE",
    "dodgerblue",
    "greenyellow",
    "orange",
    "darkolivegreen4",
    "lightblue",
    "turquoise" ,
    "red",
    "forestgreen",
    "navyblue",
    "greenyellow",
    "orange",
    "#006400"
  )

spatInSituPlotPoints(e15_xenium_Lshelf,
                     show_image = FALSE,
                     feats = NULL,
                     point_size = 0.6,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_alpha = 1,
                     polygon_fill_code = xenium_colors,
                     polygon_color = 'black',
                     polygon_line_size = 0.01,
                     polygon_fill = 'cell_types2',
                     polygon_fill_as_factor = TRUE,
                     coord_fix_ratio = TRUE,
                     save_para = list(
                       save_name = '32_e15Lshelf_spat_cell_annot'), return_plot = T)


# Plot craniosynostosis genes in e15.5 scRNA-seq mps ----------------------

CS.genes <- c("Efnb1","Fgfr1","Fgfr2","Fgfr3","Twist1")

DoHeatmap(e15mps_sc,group.colors = e15sc_colors,features = CS.genes,label = FALSE,angle=0,size = 3,lines.width = 1,group.bar.height = 0.08)+ scale_fill_gradientn(colors=c("steelblue2","white","red", "firebrick3"))+
  theme(legend.position="bottom",
        legend.text = element_text(size = 6),
        legend.title = element_text(size=6),
        legend.justification = "left",
        axis.text.y.left = element_text(face = "italic",size = 7))


DoHeatmap(e18fs_sc,group.colors = e18fs_colors,features = CS.genes,label = FALSE,angle=0,size = 3,lines.width = 1,group.bar.height = 0.08)+ scale_fill_gradientn(colors=c("steelblue2","white","red", "firebrick3"))+
  theme(legend.position="bottom",
        legend.text = element_text(size = 6),
        legend.title = element_text(size=6),
        legend.justification = "left",
        axis.text.y.left = element_text(face = "italic",size = 7))

feat_colors <- c("limegreen", "orchid","black", "dodgerblue2")

spatInSituPlotPoints(e15_xenium_Lshelf,
                     show_image = FALSE,
                     feats = list('rna' =CS.genes),
                     feats_color_code = feat_colors,
                     point_size = 1,show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     show_legend = TRUE,
                     polygon_alpha = 0.1,
                     polygon_color = 'pink2',
                     background_color = "#FFFAFABD",
                     axis_text = 8,
                     axis_title = 9,
                     polygon_line_size = 0.1,
                     polygon_fill = 'cell_types2',
                     polygon_fill_as_factor = TRUE,
                     coord_fix_ratio = TRUE,
                     polygon_fill_code = xenium_colors,
                     save_para = list(
                       save_name = '36b_e15xeniumL_lopoly_craniosynostosis'), return_plot = T)


