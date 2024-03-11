# ## ##################### ## #
#  IMPORTING XENIUM REGIONS   #
# ## ##################### ## #


# Run this code to load and process Xenium data
# All xenium regions loaded here
# Updated by Daniela M. Roth on
# Date: Mon Mar 11 05:58:04 2024 ------------------


# Make directories for overall Xenium data
library(here)
dir.create(here(home.path <-
                  here("midpalatal-sutures", "xenium")), recursive = TRUE) 
dir.create(src <- here(home.path, "src"), recursive = TRUE) 



# Region 5 ----------------------------------------------------------------

## Create region-specific directories and set paths ----------------------------
region <- "region-5" # replace region-x with correct region name for each

source(here::here("midpalatal-sutures","xenium","docs","packages.R"))
source(here::here("midpalatal-sutures/xenium/docs/functions.R"))
source(here::here("midpalatal-sutures/xenium/docs/directories.R"))

load_xenium_data()

## Preview region ----------------------------------------------------------
gpoints_list = lapply(tx_dt_types, function(x)
  createGiottoPoints(x = x))

# Preview QC probe detections
plot(gpoints_list$`Blank Codeword`,
     point_size = 1,
     main = 'Blank Codeword')
plot(gpoints_list$`Negative Control Codeword`,
     point_size = 1,
     main = 'Negative Control Codeword')
plot(gpoints_list$`Negative Control Probe`,
     point_size = 1,
     main = 'Negative Control Probe')

# Preview selected genes
mygenes <- c("Krt14","Col1a1","Tnn","Lum")

filename <-
  paste0(
    results_folder,
    "/01_preview-region_",
    str_c(mygenes, "_", collapse = ""),
    region,
    ".pdf"
  )
pdf(file =filename, width=4, height=4)
preview.spatplot <- plot(gpoints_list$`Gene Expression`,
     feats = mygenes)+theme(legend.position = "bottom")
dev.off()

tx_dt_types$`Gene Expression`[feat_ID %in% mygenes, table(feat_ID)]


load_polygon_data()

## Create Giotto Object for entire slide -----------------------------------
gobject = createGiottoObjectSubcellular(
  gpoints = list(
    rna = gpoints_list$`Gene Expression`,
    blank_code = gpoints_list$`Blank Codeword`,
    neg_code = gpoints_list$`Negative Control Codeword`,
    neg_probe = gpoints_list$`Negative Control Probe`
  ),
  gpolygons = list(cell = gpoly_cells,
                   nucleus = gpoly_nucs),
  instructions = instrs
)

saveGiotto(gobject,dir = here(output),foldername = paste0(region,"_preprocessing_","giotto-object"), overwrite = TRUE)
# File will save as gobject.RDS so descriptive folder naming is essential


# Region 6 ----------------------------------------------------------------

## Create region-specific directories and set paths ----------------------------
region <- "region-6" # replace region-x with correct region name for each

source(here::here("midpalatal-sutures","xenium","docs","packages.R"))
source(here::here("midpalatal-sutures/xenium/docs/functions.R"))
source(here::here("midpalatal-sutures/xenium/docs/directories.R"))

load_xenium_data()

## Preview region ----------------------------------------------------------
gpoints_list = lapply(tx_dt_types, function(x)
  createGiottoPoints(x = x))

# Preview QC probe detections
plot(gpoints_list$`Blank Codeword`,
     point_size = 1,
     main = 'Blank Codeword')
plot(gpoints_list$`Negative Control Codeword`,
     point_size = 1,
     main = 'Negative Control Codeword')
plot(gpoints_list$`Negative Control Probe`,
     point_size = 1,
     main = 'Negative Control Probe')

# Preview selected genes
mygenes <- c("Krt14","Col1a1","Tnn","Lum")

filename <-
  paste0(
    results_folder,
    "/01_preview-region_",
    str_c(mygenes, "_", collapse = ""),
    region,
    ".pdf"
  )
pdf(file =filename, width=4, height=4)
preview.spatplot <- plot(gpoints_list$`Gene Expression`,
                         feats = mygenes)+theme(legend.position = "bottom")
dev.off()

tx_dt_types$`Gene Expression`[feat_ID %in% mygenes, table(feat_ID)]


load_polygon_data()

## Create Giotto Object for entire slide -----------------------------------
gobject = createGiottoObjectSubcellular(
  gpoints = list(
    rna = gpoints_list$`Gene Expression`,
    blank_code = gpoints_list$`Blank Codeword`,
    neg_code = gpoints_list$`Negative Control Codeword`,
    neg_probe = gpoints_list$`Negative Control Probe`
  ),
  gpolygons = list(cell = gpoly_cells,
                   nucleus = gpoly_nucs),
  instructions = instrs
)

saveGiotto(gobject,dir = here(output),foldername = paste0(region,"_preprocessing_","giotto-object"), overwrite = TRUE)
# File will save as gobject.RDS so descriptive folder naming is essential


# Region 7 ----------------------------------------------------------------

## Create region-specific directories and set paths ----------------------------
region <- "region-7" # replace region-x with correct region name for each


source(here::here("midpalatal-sutures","xenium","docs","packages.R"))
source(here::here("midpalatal-sutures/xenium/docs/functions.R"))
source(here::here("midpalatal-sutures/xenium/docs/directories.R"))

load_xenium_data()

## Preview region ----------------------------------------------------------
gpoints_list = lapply(tx_dt_types, function(x)
  createGiottoPoints(x = x))

# Preview QC probe detections
plot(gpoints_list$`Blank Codeword`,
     point_size = 1,
     main = 'Blank Codeword')
plot(gpoints_list$`Negative Control Codeword`,
     point_size = 1,
     main = 'Negative Control Codeword')
plot(gpoints_list$`Negative Control Probe`,
     point_size = 1,
     main = 'Negative Control Probe')

# Preview selected genes
mygenes <- c("Krt14","Col1a1","Tnn","Lum")

filename <-
  paste0(
    results_folder,
    "/01_preview-region_",
    str_c(mygenes, "_", collapse = ""),
    region,
    ".pdf"
  )
pdf(file =filename, width=4, height=4)
preview.spatplot <- plot(gpoints_list$`Gene Expression`,
                         feats = mygenes)+theme(legend.position = "bottom")
dev.off()

tx_dt_types$`Gene Expression`[feat_ID %in% mygenes, table(feat_ID)]

load_polygon_data()

## Create Giotto Object for entire slide -----------------------------------
gobject = createGiottoObjectSubcellular(
  gpoints = list(
    rna = gpoints_list$`Gene Expression`,
    blank_code = gpoints_list$`Blank Codeword`,
    neg_code = gpoints_list$`Negative Control Codeword`,
    neg_probe = gpoints_list$`Negative Control Probe`
  ),
  gpolygons = list(cell = gpoly_cells,
                   nucleus = gpoly_nucs),
  instructions = instrs
)

saveGiotto(gobject,dir = here(output),foldername = paste0(region,"_preprocessing_","giotto-object"), overwrite = TRUE)
# File will save as gobject.RDS so descriptive folder naming is essential

