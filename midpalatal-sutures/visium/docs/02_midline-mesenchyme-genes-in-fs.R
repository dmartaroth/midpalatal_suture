# ## ######################################## ## #
#  MIDLINE MESENCHYME GENES IN CRANIAL SUTURES   #
# ## ######################################## ## #

# Date: Mon Apr 01 13:07:13 2024 ------------------
# Updated by: Daniela M. Roth

# Still need to make script reproducible
library(here)
source(here::here("cranial-sutures/docs/directories.R"))
source(here::here("midpalatal-sutures/visium/docs/packages.R"))
source(here::here("midpalatal-sutures/visium/docs/functions.R"))
source(here::here("midpalatal-sutures/visium/docs/themes.R"))

# Create subdirectory "Midline-mesenchyme"

results_folder <- figs
mm_figs <- file.path(figs, "Midline-mesenchyme")
if (!dir.exists(mm_figs)) {
  dir.create(mm_figs)
}


E16.fs <- readRDS(here::here("cranial-sutures/E16-fs/data-output/sc_E16_fs_clus.Rds"))
E18.fs <- readRDS(here::here("cranial-sutures/E18-fs/data-output/sc_E18_fs_annot.Rds"))
P10.fs <- readRDS(here::here("cranial-sutures/P10-fs/data-output/sc_P10_fs_annot.Rds"))
P28.fs <- readRDS(here::here("cranial-sutures/P28-fs/data-output/RDS objects/P28-fs_clus.Rds"))


# Load e15 and P1 mps datasets --------------------------------------------

E15.mps <- readRDS(here::here("midpalatal-sutures/scRNAseq/E15_mps/data-output/sc_E15_mps_annot.Rds"))
P1.mps <- readRDS(here::here("midpalatal-sutures/scRNAseq/P1_mps/data-output/sc_P1_mps_annot.Rds"))



# Compare midline mesenchyme genes over age -------------------------------
E16.fs$age <- "E16"
E18.fs$age <- "E18"
P10.fs$age <- "P10"
P28.fs$age <- "P28"

E15.mps$age <- "E15"
P1.mps$age <- "P1"

# Specify the path to the CSV file containing top genes
top_genes_file <- here::here("midpalatal-sutures/visium/e15/data-output/WT1R1_all_markers_pct.csv")
# Specify the output directory (replace mm_figs with your desired output directory)
output_dir <- mm_figs
cluster_number <- 3
top_genes_number <- 25

# Frontal sutures
seurat_obj_list <- list(E16_fs = E16.fs,
                        E18_fs = E18.fs,
                        P10_fs = P10.fs,
                        P28_fs = P28.fs)





# Specify the suture information, cluster number, and number of top genes to consider
suture <- "frontal_suture"


# Call the function with the list of Seurat objects, top genes file, output directory, suture, cluster number, and top genes number
generate_individual_violin_plots_from_file(seurat_obj_list, top_genes_file, output_dir, suture, cluster_number, top_genes_number)

# palate
# Specify the suture information, cluster number, and number of top genes to consider
suture <- "midpalatal_suture"

seurat_obj_list <- list(E15_mps = E15.mps,
                        P1_mps = P1.mps)

# Call the function with the list of Seurat objects, top genes file, output directory, suture, cluster number, and top genes number
generate_individual_violin_plots_from_file(seurat_obj_list, top_genes_file, output_dir, suture, cluster_number, top_genes_number)



