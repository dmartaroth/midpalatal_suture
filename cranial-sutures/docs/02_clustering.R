# ## ######################################## ## #
#             CLUSTER (LOWRES)                   #
# ## ######################################## ## #

# Date: Wed Apr 24 18:48:10 2024 ------------------

library(here)

source(here::here("cranial-sutures","docs","packages.R")) # load packages
source(here::here("cranial-sutures","docs","directories.R")) # load file paths/directories
source(here::here("cranial-sutures","docs","functions.R")) # load functions
source(here::here("cranial-sutures","docs","themes.R")) # load themes

# Load data
E16.fs <- readRDS(here::here("cranial-sutures/E16-fs/data-output/integrated_filtered_E16-fs.Rds"))
E18.fs <- readRDS(here::here("cranial-sutures/E18-fs/data-output/integrated_filtered_E18-fs.Rds"))
P10.fs <- readRDS(here::here("cranial-sutures/P10-fs/data-output/integrated_filtered_P10-fs.Rds"))
P28.fs <- readRDS(here::here("cranial-sutures/P28-fs/data-output/filtered_clustered_P28-fs_P28-fs-2-55H2.Rds"))



# Annotate  -----------------------------------------------------

## Assign clustering resolution --------------------------------------------

# E16
suture <- "E16-fs"
obj <- E16.fs
data.output <- here::here(cranial_sutures,suture,"data-output")

(plot <- DimPlot(obj, reduction = "umap", group.by = c( "RNA_snn_res.0.15"))+
    umap_theme() + scale_color_manual(values = pastel_palette)) 



E16.fs <- obj
Idents(E16.fs) <- E16.fs$RNA_snn_res.0.15

# E18
suture <- "E18-fs"

obj <- E18.fs
data.output <- here::here(cranial_sutures,suture,"data-output")


(plot <- DimPlot(obj, reduction = "umap", group.by = c( "RNA_snn_res.0.08"))+
    umap_theme() + scale_color_manual(values = pastel_palette)) 


E18.fs <- obj
Idents(E18.fs) <- E18.fs$RNA_snn_res.0.08

# P10
suture <- "P10-fs"

obj <- P10.fs
data.output <- here::here(cranial_sutures,suture,"data-output")


(plot <- DimPlot(obj, reduction = "umap", group.by = c( "RNA_snn_res.0.08"))+
    umap_theme() + scale_color_manual(values = pastel_palette)) 


P10.fs <- obj
Idents(P10.fs) <- P10.fs$RNA_snn_res.0.08

# P28
suture <- "P28-fs"

(plot <- DimPlot(P28.fs, reduction = "umap", group.by = c( "RNA_snn_res.0.1"))+
    umap_theme() + scale_color_manual(values = pastel_palette)) 

Idents(P28.fs) <- P28.fs$RNA_snn_res.0.1

## Annotation --------------------------------------------------------------

# List of genes for feature plots
progenitors <- c("Axin2","Gli1","Prrx1","Six2")
osteogenic <- c("Crabp1","Runx2","Sp7","Dmp1")
chondrogenic <- c("Col2a1","Acan","Mgp","Sox9")
osteoclasts <- c("Ctsk","Mmp9","Pheta1","Cd44")
vascular <- c("Mcam","Vwf","Pecam1","Pdgfrb")
myeloid_lymphocyte <- c("Pou2f2","Il1rl1","Gata2")
neurons_gli1 <- c("Neurod1","Cplx3","Otx2","Gfra3","Sox10","Foxd3")
erythrocytes <- c("Hba-a1","Hba-a2","Hbb-bs","Gypa","Gybp","Alas2","Klf1","Slc25a37","Slc2a1")
smoothmuscle <- c("Acta2","Tagln","Myh11","Des")


### E16 ---------------------------------------------------------------------

suture <- "E16-fs"
obj <- E16.fs

data.output <- here::here(cranial_sutures,suture,"data-output")
(plot <- DimPlot(obj, reduction = "umap")+
    umap_theme() + scale_color_manual(values = pastel_palette)) 

# Percent Difference in Expression
# Basic FindAllMarkers DE test
all_markers_pct <- FindAllMarkers(obj,verbose = T) %>% 
  Add_Pct_Diff()

all_markers_pct <- all_markers_pct %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  arrange(cluster)

write.csv(all_markers_pct, file = here(data.output, "all_markers_pct.csv"))

# Extract the top N marker genes per cluster for plotting
top_5 <- Extract_Top_Markers(marker_dataframe = all_markers_pct, num_genes = 5, rank_by = "avg_log2FC")


top50_markers_pct <- all_markers_pct %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>% 
  top_n(n=50, wt = avg_log2FC)%>%
  arrange(cluster)

write.csv(top50_markers_pct, file = here(data.output, "top50_markers_pct.csv"))


(plot <- DotPlot(
  object = obj,
  features =   top_5,
  scale.by = "radius",
  dot.scale = 8,
  split.by = NULL,
  cluster.idents = FALSE,
) + scale_colour_gradient2(low = "dodgerblue",
                           mid = "floralwhite",
                           high = "red2") +  custom_dotplot_theme() +RotatedAxis())

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = progenitors,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = osteogenic,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = chondrogenic,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = osteoclasts,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = vascular,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = myeloid_lymphocyte,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = neurons_gli1,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = erythrocytes,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = smoothmuscle,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Create simple annotation files
# Create_Cluster_Annotation_File(file_path = data.output, file_name = "cluster_annotation")

annotation_info <- Pull_Cluster_Annotation(annotation = here(data.output,"cluster_annotation.csv"))

# Rename clusters
obj_annot <- Rename_Clusters(seurat_object = obj, new_idents = annotation_info$new_cluster_idents, meta_col_name = "mes_annotation")

Idents(obj_annot) <- factor(x = Idents(obj_annot), levels = sort(levels(obj_annot)))
obj_annot$mesenchyme <- Idents(obj_annot)

(plot <- DimPlot(obj_annot, reduction = "umap",label = TRUE,repel = TRUE,label.size = 3,label.box = TRUE,cols = pastel_palette)+
    umap_theme()) 

# Save annotated object
saveRDS(obj_annot, file = paste0(data.output,"/annot_",suture,".Rds"))



### E18 ---------------------------------------------------------------------


suture <- "E18-fs"
obj <- E18.fs

data.output <- here::here(cranial_sutures,suture,"data-output")
(plot <- DimPlot(obj, reduction = "umap")+
   umap_theme() + scale_color_manual(values = pastel_palette)) 

# Percent Difference in Expression
# Basic FindAllMarkers DE test
all_markers_pct <- FindAllMarkers(obj,verbose = T) %>% 
  Add_Pct_Diff()

all_markers_pct <- all_markers_pct %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  arrange(cluster)

write.csv(all_markers_pct, file = here(data.output, "all_markers_pct.csv"))

# Extract the top N marker genes per cluster for plotting
top_5 <- Extract_Top_Markers(marker_dataframe = all_markers_pct, num_genes = 5, rank_by = "avg_log2FC")


top50_markers_pct <- all_markers_pct %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>% 
  top_n(n=50, wt = avg_log2FC)%>%
  arrange(cluster)

write.csv(top50_markers_pct, file = here(data.output, "top50_markers_pct.csv"))


(plot <- DotPlot(
  object = obj,
  features =   top_5,
  scale.by = "radius",
  dot.scale = 8,
  split.by = NULL,
  cluster.idents = FALSE,
) + scale_colour_gradient2(low = "dodgerblue",
                           mid = "floralwhite",
                           high = "red2") +  custom_dotplot_theme() +RotatedAxis())

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = progenitors,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = osteogenic,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = chondrogenic,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = osteoclasts,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = vascular,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = myeloid_lymphocyte,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = neurons_gli1,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = erythrocytes,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = smoothmuscle,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Create simple annotation files
# Create_Cluster_Annotation_File(file_path = data.output, file_name = "cluster_annotation")

annotation_info <- Pull_Cluster_Annotation(annotation = here(data.output,"cluster_annotation.csv"))

# Rename clusters
obj_annot <- Rename_Clusters(seurat_object = obj, new_idents = annotation_info$new_cluster_idents, meta_col_name = "mes_annotation")

Idents(obj_annot) <- factor(x = Idents(obj_annot), levels = sort(levels(obj_annot)))
obj_annot$mesenchyme <- Idents(obj_annot)

(plot <- DimPlot(obj_annot, reduction = "umap",label = TRUE,repel = TRUE,label.size = 3,label.box = TRUE,cols = pastel_palette)+
    umap_theme()) 

# Save annotated object
saveRDS(obj_annot, file = paste0(data.output,"/annot_",suture,".Rds"))



### P10 ---------------------------------------------------------------------
suture <- "P10-fs"
obj <- P10.fs

data.output <- here::here(cranial_sutures,suture,"data-output")
(plot <- DimPlot(obj, reduction = "umap")+
    umap_theme() + scale_color_manual(values = pastel_palette)) 

# Percent Difference in Expression
# Basic FindAllMarkers DE test
all_markers_pct <- FindAllMarkers(obj,verbose = T) %>% 
  Add_Pct_Diff()

all_markers_pct <- all_markers_pct %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  arrange(cluster)

write.csv(all_markers_pct, file = here(data.output, "all_markers_pct.csv"))

# Extract the top N marker genes per cluster for plotting
top_5 <- Extract_Top_Markers(marker_dataframe = all_markers_pct, num_genes = 5, rank_by = "avg_log2FC")


top50_markers_pct <- all_markers_pct %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>% 
  top_n(n=50, wt = avg_log2FC)%>%
  arrange(cluster)

write.csv(top50_markers_pct, file = here(data.output, "top50_markers_pct.csv"))


(plot <- DotPlot(
  object = obj,
  features =   top_5,
  scale.by = "radius",
  dot.scale = 8,
  split.by = NULL,
  cluster.idents = FALSE,
) + scale_colour_gradient2(low = "dodgerblue",
                           mid = "floralwhite",
                           high = "red2") +  custom_dotplot_theme() +RotatedAxis())

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = progenitors,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = osteogenic,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = chondrogenic,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = osteoclasts,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = vascular,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = myeloid_lymphocyte,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = neurons_gli1,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = erythrocytes,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = smoothmuscle,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Create simple annotation files
# Create_Cluster_Annotation_File(file_path = data.output, file_name = "cluster_annotation")

annotation_info <- Pull_Cluster_Annotation(annotation = here(data.output,"cluster_annotation.csv"))

# Rename clusters
obj_annot <- Rename_Clusters(seurat_object = obj, new_idents = annotation_info$new_cluster_idents, meta_col_name = "mes_annotation")

Idents(obj_annot) <- factor(x = Idents(obj_annot), levels = sort(levels(obj_annot)))
obj_annot$mesenchyme <- Idents(obj_annot)

(plot <- DimPlot(obj_annot, reduction = "umap",label = TRUE,repel = TRUE,label.size = 3,label.box = TRUE,cols = pastel_palette)+
    umap_theme()) 

# Save annotated object
saveRDS(obj_annot, file = paste0(data.output,"/annot_",suture,".Rds"))


### P28 ---------------------------------------------------------------------
suture <- "P28-fs"
obj <- P28.fs

data.output <- here::here(cranial_sutures,suture,"data-output")
(plot <- DimPlot(obj, reduction = "umap")+
    umap_theme() + scale_color_manual(values = pastel_palette)) 

# Percent Difference in Expression
# Basic FindAllMarkers DE test
all_markers_pct <- FindAllMarkers(obj,verbose = T) %>% 
  Add_Pct_Diff()

all_markers_pct <- all_markers_pct %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  arrange(cluster)

write.csv(all_markers_pct, file = here(data.output, "all_markers_pct.csv"))

# Extract the top N marker genes per cluster for plotting
top_5 <- Extract_Top_Markers(marker_dataframe = all_markers_pct, num_genes = 5, rank_by = "avg_log2FC")


top50_markers_pct <- all_markers_pct %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>% 
  top_n(n=50, wt = avg_log2FC)%>%
  arrange(cluster)

write.csv(top50_markers_pct, file = here(data.output, "top50_markers_pct.csv"))


(plot <- DotPlot(
  object = obj,
  features =   top_5,
  scale.by = "radius",
  dot.scale = 8,
  split.by = NULL,
  cluster.idents = FALSE,
) + scale_colour_gradient2(low = "dodgerblue",
                           mid = "floralwhite",
                           high = "red2") +  custom_dotplot_theme() +RotatedAxis())

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = progenitors,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = osteogenic,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = chondrogenic,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = osteoclasts,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = vascular,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = myeloid_lymphocyte,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = neurons_gli1,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = erythrocytes,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

(plot <-
    FeaturePlot_scCustom(
      seurat_object = obj,
      reduction = "umap",
      na_cutoff = 0,
      features = smoothmuscle,
      colors_use = c(
        "floralwhite",
        "lavenderblush",
        "plum1",
        "orchid",
        "orchid4",
        "darkorchid4")))

# Create simple annotation files
# Create_Cluster_Annotation_File(file_path = data.output, file_name = "cluster_annotation")

annotation_info <- Pull_Cluster_Annotation(annotation = here(data.output,"cluster_annotation.csv"))

# Rename clusters
obj_annot <- Rename_Clusters(seurat_object = obj, new_idents = annotation_info$new_cluster_idents, meta_col_name = "mes_annotation")

Idents(obj_annot) <- factor(x = Idents(obj_annot), levels = sort(levels(obj_annot)))
obj_annot$mesenchyme <- Idents(obj_annot)

(plot <- DimPlot(obj_annot, reduction = "umap",label = TRUE,repel = TRUE,label.size = 3,label.box = TRUE,cols = pastel_palette)+
    umap_theme()) 

# Save annotated object
saveRDS(obj_annot, file = paste0(data.output,"/annot_",suture,".Rds"))



