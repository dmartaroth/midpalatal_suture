# Preprocessing with PopsicleR --------------------------------------------

# Path definitions --------------------------------------------------------

# Load libraries
library(here)
library(Seurat)
library(tidyverse)
#library(popsicleR)
library(dplyr)
library(patchwork)
library(ggplot2)
library(crayon)
library(gprofiler2)
library(clustree)
library(SingleR)


# Make subdirectory for suture
suture <- "P28-fs" # replace quoted text with name of suture being analyzed
dir.create(home.path <- here("cranial-sutures",suture), recursive = TRUE)

# Create directories for data-output, docs, figures, src, and raw-data
dir.create(data.output <- here(home.path,"data-output"), recursive = TRUE)
dir.create(docs <- here(home.path,"docs"), recursive = TRUE)
dir.create(figs <- here(home.path,"figures"), recursive = TRUE)
dir.create(src <- here(home.path,"src"), recursive = TRUE)
dir.create(data.path <- here(home.path,"raw-data"), recursive = TRUE)



#Downloaded scWFP28_1655703_matrix.count, scWFP28_1655703_features.tsv, and scWFP28_1655703_barcodes.tsv from Facebase (ID 2-55GE[scWFP28])
#Files had different names. Opened in txt editor and saved as "barcodes.tsv.gz", "features.tsv.gz", and "matrix.mtx.gz"


# Feb 13 2024 can't use PopsicleR anymore due to Seuratv5; rewrite this code for seurat analysis


#Load count matrices
sample.name=suture
input.data.dir=data.path
outs = here(figs,"PopsicleR-preprocessing")
sample.umi = PrePlots(sample = sample.name, input_data = input.data.dir, percentage = 0.1, gene_filter = 200, cellranger=TRUE,organism = "mouse", out_folder = here(figs,"PopsicleR-preprocessing"))

# Apply parameters for data filtering based on violin plots above
sample.umi=FilterPlots(UMI=sample.umi,G_RNA_low = 200,G_RNA_hi = 5000,U_RNA_low = 200,U_RNA_hi = 45000,percent_mt_hi = 20,percent_ribo_hi = 40,percent_disso_hi = 5,out_folder = outs)

sample.umi <- CalculateDoublets(UMI = sample.umi, method = "scrublet", dbs_thr ='none', dbs_remove = FALSE, out_folder = outs) # calculate doublets
sample.umi <- CalculateDoublets(UMI = sample.umi, method = "scrublet", dbs_thr =0.27, dbs_remove = TRUE, out_folder = outs) # set doublet threshold based on previous plot

# Data normalization and cell cycle scoring
sample.umi <- Normalize(UMI = sample.umi, variable_genes = 2000, out_folder = outs)
sample.umi <- ApplyRegression(UMI = sample.umi, organism = "mouse", variables = "none", explore_PC = FALSE, out_folder = outs)
sample.umi <- ApplyRegression(UMI = sample.umi, organism = "mouse", variables = c("S.Score", "G2M.Score"), explore_PC = TRUE, out_folder = outs)

#input elbow value below for pcs
pcs.use <- 10
pre.cluster <- sample.umi
post.cluster <- CalculateCluster(UMI = pre.cluster, dim_pca = pcs.use, organism = "mouse", marker.list = "none", PCA = TRUE, cluster_res = 0.6, out_folder=outs)

output.data.dir <- here(data.output,"RDS objects")
if (!file.exists(output.data.dir)){dir.create(output.data.dir)}
saveRDS(post.cluster,file = file.path(output.data.dir,paste0(sample.name,"_clus.Rds")))


