# ## ######################################## ## #
#         IMPORT FRONTAL SUTURE SCRNASEQ         #
# ## ######################################## ## #

# Date: Wed Apr 24 17:01:41 2024 ------------------

library(here)

source(here::here("cranial-sutures","docs","packages.R")) # load packages
source(here::here("cranial-sutures","docs","directories.R")) # load file paths/directories
source(here::here("cranial-sutures","docs","functions.R")) # load functions
source(here::here("cranial-sutures","docs","themes.R")) # load themes

# E16 frontal suture ------------------------------------------------------
# Sample 1
suture <- "E16-fs"
identifier <- "E16-fs-1-P588"

(dir.create(output_dir <- here::here(cranial_sutures,suture,"figs")))
(dir.create(data.output <- here::here(cranial_sutures,suture,"data-output")))
data.dir <- here::here(cranial_sutures,suture,"raw-data",identifier)

# Load 10x files. If not in raw-data folder nested within identifier folder,
# create and move appropriately.

obj <- Read10X(data.dir = data.dir)
obj <- CreateSeuratObject(counts = obj, project = suture,min.cells = 3, min.features = 200)


prepro.plots(obj,output_dir = output_dir)

obj <- add_percent_mito(obj)


filtered_obj <- subset(x = obj,
                        subset = (nFeature_RNA > 200) & # value from seurat vignette
                          (nFeature_RNA < 10000) &
                          (nCount_RNA > 200) &
                          (nCount_RNA < 100000) &
                          (percent.mito < 0.15))

filtered_obj <- genelvlfilt(filtered_obj)

filt.plots(filtered_obj, output_dir = output_dir)

E16a <- filtered_obj

saveRDS(filtered_obj, file = paste0(data.output,"/filtered_",suture,"_",identifier,".Rds"))


# Sample 2
suture <- "E16-fs"
identifier <- "E16-fs-1-P58A"


(dir.create(output_dir <- here::here(cranial_sutures,suture,"figs", identifier)))
data.dir <- here::here(cranial_sutures,suture,"raw-data",identifier)
(dir.create(data.output <- here::here(cranial_sutures,suture,"data-output")))

# Load 10x files. If not in raw-data folder nested within identifier folder,
# create and move appropriately.

obj <- Read10X(data.dir = data.dir)
obj <- CreateSeuratObject(counts = obj, project = suture,min.cells = 3, min.features = 200)


prepro.plots(obj,output_dir = output_dir)

obj <- add_percent_mito(obj)


filtered_obj <- subset(x = obj,
                       subset = (nFeature_RNA > 200) & # value from seurat vignette
                         (nFeature_RNA < 2000) &
                         (nCount_RNA > 200) &
                         (nCount_RNA < 50000) &
                         (percent.mito < 0.4))

filtered_obj <- genelvlfilt(filtered_obj)

filt.plots(filtered_obj, output_dir = output_dir)

E16b <- filtered_obj

saveRDS(filtered_obj, file = paste0(data.output,"/filtered_",suture,"_",identifier,".Rds"))

# Integrate
E16a$replicate <- "1"
E16b$replicate <- "2"
merged_E16 <- merge(E16a, y = E16b, add.cell.ids = c("1","2"), project = suture)
merged_E16 <- JoinLayers(merged_E16)
merged_E16[["RNA"]] <- split(merged_E16[["RNA"]], f = merged_E16$replicate)

# Analysis without integration
merged_E16 <- NormalizeData(merged_E16)
merged_E16 <- FindVariableFeatures(merged_E16)
merged_E16 <- ScaleData(merged_E16)
merged_E16 <- RunPCA(merged_E16)

merged_E16 <- FindNeighbors(merged_E16, dims = 1:30, reduction = "pca")
merged_E16 <- FindClusters(merged_E16, resolution = 1, cluster.name = "unintegrated_clusters")

# View unintegrated UMAP
merged_E16 <- RunUMAP(merged_E16, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
(DimPlot(merged_E16, reduction = "umap.unintegrated", group.by = c("replicate"),alpha = 0.5) +
  umap_theme() + scale_color_manual(values = c("orchid","olivedrab2"))   )

# Integration
merged_E16 <-
  IntegrateLayers(
    object = merged_E16,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    verbose = FALSE)

merged_E16[["RNA"]] <- JoinLayers(merged_E16[["RNA"]]) # re-join layers after integration

merged_E16 <- FindNeighbors(merged_E16, reduction = "integrated.cca", dims = 1:30)

merged_E16 <- FindClusters(merged_E16, verbose=FALSE,
                               resolution= c(0.05, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
)

(plot <- clustree(merged_E16, prefix = "RNA_snn_res.") )

merged_E16 <- RunUMAP(merged_E16, dims = 1:30, reduction = "integrated.cca")
(plot <- DimPlot(merged_E16, reduction = "umap", group.by = c("replicate", "RNA_snn_res.0.15"))+
    umap_theme() + scale_color_manual(values = pastel_palette)) 

(DimPlot(merged_E16, reduction = "umap", split.by = "replicate", group.by = "RNA_snn_res.0.15")+
  umap_theme() + scale_color_manual(values = pastel_palette))

saveRDS(merged_E16, file = paste0(data.output,"/integrated_filtered_",suture,".Rds"))

# E18 frontal suture ------------------------------------------------------

# Replicate 1
suture <- "E18-fs"
identifier <- "E18-fs-1-C3QM"


(dir.create(output_dir <- here::here(cranial_sutures,suture,"figs",identifier)))
(dir.create(data.output <- here::here(cranial_sutures,suture,"data-output")))
data.dir <- here::here(cranial_sutures,suture,"raw-data",identifier)

# Load 10x files. If not in raw-data folder nested within identifier folder,
# create and move appropriately.

obj <- Read10X(data.dir = data.dir)
obj <- CreateSeuratObject(counts = obj, project = suture,min.cells = 3, min.features = 200)


prepro.plots(obj,output_dir = output_dir)

obj <- add_percent_mito(obj)


filtered_obj <- subset(x = obj,
                       subset = (nFeature_RNA > 200) & # value from seurat vignette
                         (nFeature_RNA < 6000) &
                         (nCount_RNA > 200) &
                         (nCount_RNA < 40000) &
                         (percent.mito < 0.15))

filtered_obj <- genelvlfilt(filtered_obj)

filt.plots(filtered_obj, output_dir = output_dir)

E18a <- filtered_obj

saveRDS(filtered_obj, file = paste0(data.output,"/filtered_",suture,"_",identifier,".Rds"))

# Replicate 2
suture <- "E18-fs"
identifier <- "E18-fs-1-8DVP"


(dir.create(output_dir <- here::here(cranial_sutures,suture,"figs",identifier)))
(dir.create(data.output <- here::here(cranial_sutures,suture,"data-output")))
data.dir <- here::here(cranial_sutures,suture,"raw-data",identifier)

# Load 10x files. If not in raw-data folder nested within identifier folder,
# create and move appropriately.

raw_counts <- read.table(here::here(data.dir,"counts.count"))
obj <- CreateSeuratObject(counts = raw_counts,project = suture,min.cells = 3, min.features = 200)


prepro.plots(obj,output_dir = output_dir)

obj <- add_percent_mito(obj)


filtered_obj <- subset(x = obj,
                       subset = (nFeature_RNA > 200) & # value from seurat vignette
                         (nFeature_RNA < 7500) &
                         (nCount_RNA > 200) &
                         (nCount_RNA < 100000) &
                         (percent.mito < 0.2))

filtered_obj <- genelvlfilt(filtered_obj)

filt.plots(filtered_obj, output_dir = output_dir)

E18b <- filtered_obj

saveRDS(filtered_obj, file = paste0(data.output,"/filtered_",suture,"_",identifier,".Rds"))


# Integrate
E18a$replicate <- "1"
E18b$replicate <- "2"
merged_E18 <- merge(E18a, y = E18b, add.cell.ids = c("1","2"), project = suture)
merged_E18 <- JoinLayers(merged_E18)
merged_E18[["RNA"]] <- split(merged_E18[["RNA"]], f = merged_E18$replicate)

# Analysis without integration
merged_E18 <- NormalizeData(merged_E18)
merged_E18 <- FindVariableFeatures(merged_E18)
merged_E18 <- ScaleData(merged_E18)
merged_E18 <- RunPCA(merged_E18)

merged_E18 <- FindNeighbors(merged_E18, dims = 1:30, reduction = "pca")
merged_E18 <- FindClusters(merged_E18, resolution = 1, cluster.name = "unintegrated_clusters")

# View unintegrated UMAP
merged_E18 <- RunUMAP(merged_E18, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
(DimPlot(merged_E18, reduction = "umap.unintegrated", group.by = c("replicate"),alpha = 0.5) +
    umap_theme() + scale_color_manual(values = c("orchid","olivedrab2"))   )

# Integration
merged_E18 <-
  IntegrateLayers(
    object = merged_E18,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    verbose = FALSE)

merged_E18[["RNA"]] <- JoinLayers(merged_E18[["RNA"]]) # re-join layers after integration

merged_E18 <- FindNeighbors(merged_E18, reduction = "integrated.cca", dims = 1:30)

merged_E18 <- FindClusters(merged_E18, verbose=FALSE,
                           resolution= c(0.05, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
)

(plot <- clustree(merged_E18, prefix = "RNA_snn_res.") )

merged_E18 <- RunUMAP(merged_E18, dims = 1:30, reduction = "integrated.cca")

(DimPlot(merged_E18, reduction = "umap", group.by = c("replicate"),alpha = 0.5) +
    umap_theme() + scale_color_manual(values = c("orchid","olivedrab2"))   )

(DimPlot(merged_E18, reduction = "umap", split.by = "replicate", group.by = "RNA_snn_res.0.08")+
    umap_theme() + scale_color_manual(values = pastel_palette))

saveRDS(merged_E18, file = paste0(data.output,"/integrated_filtered_",suture,".Rds"))


# P10 frontal suture ------------------------------------------------------

# Replicate 1
suture <- "P10-fs"
identifier <- "P10-fs-2-55HC"

(dir.create(output_dir <- here::here(cranial_sutures,suture,"figs",identifier)))
(dir.create(data.output <- here::here(cranial_sutures,suture,"data-output")))
data.dir <- here::here(cranial_sutures,suture,"raw-data",identifier)

# Load 10x files. If not in raw-data folder nested within identifier folder,
# create and move appropriately.

obj <- Read10X(data.dir = data.dir)
obj <- CreateSeuratObject(counts = obj, project = suture,min.cells = 3, min.features = 200)


prepro.plots(obj,output_dir = output_dir)

obj <- add_percent_mito(obj)


filtered_obj <- subset(x = obj,
                       subset = (nFeature_RNA > 200) & # value from seurat vignette
                         (nFeature_RNA < 5000) &
                         (nCount_RNA > 200) &
                         (nCount_RNA < 40000) &
                         (percent.mito < 0.1))

filtered_obj <- genelvlfilt(filtered_obj)

filt.plots(filtered_obj, output_dir = output_dir)

P10a <- filtered_obj
saveRDS(filtered_obj, file = paste0(data.output,"/filtered_",suture,"_",identifier,".Rds"))

# Replicate 2
suture <- "P10-fs"
identifier <- "P10-fs-2-55HE"

(dir.create(output_dir <- here::here(cranial_sutures,suture,"figs",identifier)))
(dir.create(data.output <- here::here(cranial_sutures,suture,"data-output")))
data.dir <- here::here(cranial_sutures,suture,"raw-data",identifier)

# Load 10x files. If not in raw-data folder nested within identifier folder,
# create and move appropriately.

obj <- Read10X(data.dir = data.dir)
obj <- CreateSeuratObject(counts = obj, project = suture,min.cells = 3, min.features = 200)


prepro.plots(obj,output_dir = output_dir)

obj <- add_percent_mito(obj)


filtered_obj <- subset(x = obj,
                       subset = (nFeature_RNA > 200) & # value from seurat vignette
                         (nFeature_RNA < 6000) &
                         (nCount_RNA > 200) &
                         (nCount_RNA < 30000) &
                         (percent.mito < 0.08))

filtered_obj <- genelvlfilt(filtered_obj)

filt.plots(filtered_obj, output_dir = output_dir)

P10b <- filtered_obj
saveRDS(filtered_obj, file = paste0(data.output,"/filtered_",suture,"_",identifier,".Rds"))


# Integrate
P10a$replicate <- "1"
P10b$replicate <- "2"
merged_P10 <- merge(P10a, y = P10b, add.cell.ids = c("1","2"), project = suture)
merged_P10 <- JoinLayers(merged_P10)
merged_P10[["RNA"]] <- split(merged_P10[["RNA"]], f = merged_P10$replicate)

# Analysis without integration
merged_P10 <- NormalizeData(merged_P10)
merged_P10 <- FindVariableFeatures(merged_P10)
merged_P10 <- ScaleData(merged_P10)
merged_P10 <- RunPCA(merged_P10)

merged_P10 <- FindNeighbors(merged_P10, dims = 1:30, reduction = "pca")
merged_P10 <- FindClusters(merged_P10, resolution = 1, cluster.name = "unintegrated_clusters")

# View unintegrated UMAP
merged_P10 <- RunUMAP(merged_P10, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
(DimPlot(merged_P10, reduction = "umap.unintegrated", group.by = c("replicate"),alpha = 0.5) +
    umap_theme() + scale_color_manual(values = c("orchid","olivedrab2"))   )

# Integration
merged_P10 <-
  IntegrateLayers(
    object = merged_P10,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    verbose = FALSE)

merged_P10[["RNA"]] <- JoinLayers(merged_P10[["RNA"]]) # re-join layers after integration

merged_P10 <- FindNeighbors(merged_P10, reduction = "integrated.cca", dims = 1:30)

merged_P10 <- FindClusters(merged_P10, verbose=FALSE,
                           resolution= c(0.05, 0.08, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 1)
)

(plot <- clustree(merged_P10, prefix = "RNA_snn_res.") )

merged_P10 <- RunUMAP(merged_P10, dims = 1:30, reduction = "integrated.cca")

(DimPlot(merged_P10, reduction = "umap", group.by = c("replicate"),alpha = 0.5) +
    umap_theme() + scale_color_manual(values = c("orchid","olivedrab2"))   )

(DimPlot(merged_P10, reduction = "umap", split.by = "replicate", group.by = "RNA_snn_res.0.08")+
    umap_theme() + scale_color_manual(values = pastel_palette))

saveRDS(merged_P10, file = paste0(data.output,"/integrated_filtered_",suture,".Rds"))


# P28 frontal suture ------------------------------------------------------

suture <- "P28-fs"
identifier <- "P28-fs-2-55H2"

(dir.create(output_dir <- here::here(cranial_sutures,suture,"figs",identifier)))
(dir.create(data.output <- here::here(cranial_sutures,suture,"data-output")))
data.dir <- here::here(cranial_sutures,suture,"raw-data",identifier)

# Load 10x files. If not in raw-data folder nested within identifier folder,
# create and move appropriately.

obj <- Read10X(data.dir = data.dir)
obj <- CreateSeuratObject(counts = obj, project = suture,min.cells = 3, min.features = 200)


prepro.plots(obj,output_dir = output_dir)

obj <- add_percent_mito(obj)


filtered_obj <- subset(x = obj,
                       subset = (nFeature_RNA > 200) & # value from seurat vignette
                         (nFeature_RNA < 5000) &
                         (nCount_RNA > 200) &
                         (nCount_RNA < 45000) &
                         (percent.mito < 0.2))

filtered_obj <- genelvlfilt(filtered_obj)

filt.plots(filtered_obj, output_dir = output_dir)

saveRDS(filtered_obj, file = paste0(data.output,"/filtered_",suture,"_",identifier,".Rds"))

# Only 1 P28 sample so no integration necessary
obj <- filtered_obj
data.output <- here::here(cranial_sutures,suture,"data-output")

obj <-
  NormalizeData(obj,
                normalization.method = "LogNormalize",
                scale.factor = 10000)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
ElbowPlot(obj)

obj <- FindNeighbors(obj, dims = 1:19, reduction = "pca")
obj <- FindClusters(obj,resolution = c(0.1, 0.3, 0.5, 0.7, 0.8, 1.0))
obj <- RunUMAP(obj, dims = 1:19, reduction = "pca")
(plot <- clustree(obj, prefix = "RNA_snn_res."))
(plot <- DimPlot(obj, reduction = "umap", group.by = c( "RNA_snn_res.0.1"))+
    umap_theme() + scale_color_manual(values = pastel_palette)) 

saveRDS(obj, file = paste0(data.output,"/filtered_clustered_",suture,"_",identifier,".Rds"))
