# ## ######################################## ## #
#         PREPROCESSING E15 MPS SCRNA-SEQ        #
# ## ######################################## ## #

# This pre-processing was performed using Seurat v4 and the popsicleR package
# PopsicleR is not compatible with Seurat v5; for further analysis beyond
## generated annotated object, functions must be rewritten for compatibility


#Make subdirectories
# Edit June 2024: Should be edited for reproducibility using here() package

dir.create("data-output")
dir.create("docs")
dir.create("figures")
dir.create("src")
dir.create("raw-data")

#make seurat object for E15 wt midpalatal suture
library(Seurat)
library(tidyverse)
library(popsicleR)
library(dplyr)
library(patchwork)
library(ggplot2)

#Files had different names. Opened in txt editor and saved as "barcodes.tsv.gz",
#"features.tsv.gz", and "matrix.mtx.gz"

#Load count matrices
sample.name="mps"
input.data.dir=file.path("raw-data")
sample.umi = PrePlots(sample = "mps", input_data = input.data.dir, percentage = 0.1, gene_filter = 200, cellranger=TRUE,organism = "mouse", out_folder = "data-output/PopsicleR preprocessing")
outs="data-output/PopsicleR preprocessing"

sample.umi=FilterPlots(UMI=sample.umi,G_RNA_low = 200,G_RNA_hi = 8000,U_RNA_low = 200,U_RNA_hi = 30000,percent_mt_hi = 30,percent_ribo_hi = 50,percent_disso_hi = 7,out_folder = outs)

sample.umi <- CalculateDoublets(UMI = sample.umi, method = "scrublet", dbs_thr ='none', dbs_remove = FALSE, out_folder = outs)
sample.umi <- CalculateDoublets(UMI = sample.umi, method = "scrublet", dbs_thr =0.39, dbs_remove = TRUE, out_folder = outs)

sample.umi <- Normalize(UMI = sample.umi, variable_genes = 2000, out_folder = outs)
sample.umi <- ApplyRegression(UMI = sample.umi, organism = "mouse", variables = "none", explore_PC = FALSE, out_folder = outs)
sample.umi <- ApplyRegression(UMI = sample.umi, organism = "mouse", variables = c("S.Score", "G2M.Score"), explore_PC = TRUE, out_folder = outs)

sample.umi <- CalculateCluster(UMI = sample.umi, dim_pca = 10, organism = "mouse", marker.list = "none", PCA = TRUE, cluster_res = 0.4, out_folder=outs)

output.data.dir <- file.path("data-output/RDS objects")
if (!file.exists(output.data.dir)){dir.create(output.data.dir)}
saveRDS(sample.umi,file = file.path(output.data.dir,paste0(sample.name,".Rds")))


