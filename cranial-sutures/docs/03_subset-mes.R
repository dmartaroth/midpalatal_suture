# ## ######################################## ## #
#                 SUBSET MES CLUSTERS            #
# ## ######################################## ## #

# Date: Thu Apr 25 17:57:18 2024 ------------------

library(here)

source(here::here("cranial-sutures","docs","packages.R")) # load packages
source(here::here("cranial-sutures","docs","directories.R")) # load file paths/directories
source(here::here("cranial-sutures","docs","functions.R")) # load functions
source(here::here("cranial-sutures","docs","themes.R")) # load themes

# Load data
E16.fs <- readRDS(here::here("cranial-sutures/E16-fs/data-output/annot_E16-fs.Rds"))
E18.fs <- readRDS(here::here("cranial-sutures/E18-fs/data-output/annot_E18-fs.Rds"))
P10.fs <- readRDS(here::here("cranial-sutures/P10-fs/data-output/annot_P10-fs.Rds"))
P28.fs <- readRDS(here::here("cranial-sutures/P28-fs/data-output/annot_P28-fs.Rds"))


colors <- c("violet","darkolivegreen2")
(p1 <- DimPlot(E16.fs, reduction = "umap",label = FALSE,repel = TRUE,label.size = 3,label.box = TRUE,cols = colors)+
  umap_theme())
(p2 <- DimPlot(E18.fs, reduction = "umap",label = FALSE,repel = TRUE,label.size = 3,label.box = TRUE,cols = colors)+
    umap_theme())
(p3 <- DimPlot(P10.fs, reduction = "umap",label = FALSE,repel = TRUE,label.size = 3,label.box = TRUE,cols = colors)+
    umap_theme())
(p4 <- DimPlot(P28.fs, reduction = "umap",label = FALSE,repel = TRUE,label.size = 3,label.box = TRUE,cols = colors)+
    umap_theme())

annotated_umaps <- plot_grid(p1,p2,p3,p4,
          ncol = 2)


dir.create(suppl_data <-
             here("pubfigs","supplementary"),
           recursive = TRUE)
filename <- paste0("S5_frontal-sutures_mes-other-annotation_UMAP.png")
ggsave(file.path(suppl_data, filename), annotated_umaps, width = 7, height = 6, dpi = 300)


# Subset by mes cluster identity in mes_annotation column -----------------

E16.fs.mes <- subset(E16.fs, idents= "mes")
E18.fs.mes <- subset(E18.fs, idents = "mes") 
P10.fs.mes <- subset(P10.fs, idents = "mes")
P28.fs.mes <- subset(P28.fs, idents = "mes")


p1 <- DimPlot(E16.fs.mes, reduction = "umap",label = FALSE,repel = TRUE,label.size = 3,label.box = TRUE,cols = pastel_palette)+
    umap_theme()

p2 <- DimPlot(E18.fs.mes, reduction = "umap",label = FALSE,repel = TRUE,label.size = 3,label.box = TRUE,cols = pastel_palette)+
  umap_theme()

p3 <- DimPlot(P10.fs.mes, reduction = "umap",label = FALSE,repel = TRUE,label.size = 3,label.box = TRUE,cols = pastel_palette)+
  umap_theme()

p4 <- DimPlot(P28.fs.mes, reduction = "umap",label = FALSE,repel = TRUE,label.size = 3,label.box = TRUE,cols = pastel_palette)+
  umap_theme()

plot_grid(p1,p2,p3,p4,
          ncol = 2)


# Save subsets ------------------------------------------------------------

saveRDS(E16.fs.mes, file = here::here("cranial-sutures/data-output/mes-subset_E16-fs.Rds"))
saveRDS(E18.fs.mes, file = here::here("cranial-sutures/data-output/mes-subset_E18-fs.Rds"))
saveRDS(P10.fs.mes, file = here::here("cranial-sutures/data-output/mes-subset_P10-fs.Rds"))
saveRDS(P28.fs.mes, file = here::here("cranial-sutures/data-output/mes-subset_P28-fs.Rds"))
