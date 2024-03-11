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

# Create Giotto instructions for saving
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  return_plot = TRUE)


if (!file.exists(file.path(xenium_folder, "cell_feature_matrix", "features-blank.tsv.gz"))) {
  features <-
    read_tsv(here(xenium_folder, "cell_feature_matrix", "features.tsv.gz"))
  View(features)
  features[features == "Unassigned Codeword"] <- "Blank Codeword"
  write_tsv(
    x = features,
    file = here(
      xenium_folder,
      "cell_feature_matrix",
      "features-blank.tsv.gz"
    )
  )
  cat(bold(magenta("This package looks for Blank Codeword in the place of Unassigned Codeword, so first we have to rename them in the features.tsv.doc.\n")))
  cat(bold(magenta("It looks like you haven't run this code before. Blank codewords have now been successfully renamed and the file has been saved as features-blank.tsv.gz.\n")))
} else {
  cat(bold(magenta("This package looks for Blank Codeword in the place of Unassigned Codeword, so first we have to rename them from the features.tsv.doc.\n")))
    cat(bold(magenta("The file 'features-blank.tsv.gz' already exists. Skipping renaming.\n")))
}

# General files (some are supplemental files)
settings_path = paste0(xenium_folder, '/experiment.xenium')
he_img_path = paste0(xenium_folder, '/pyramidalhe.ome.tif')
panel_meta_path = paste0(xenium_folder, "/xenium_panel.tsv") # (optional)

# Files (SUBCELLULAR): 
cell_bound_path = paste0(xenium_folder, '/cell_boundaries.csv.gz')
nuc_bound_path = paste0(xenium_folder, '/nucleus_boundaries.csv.gz')
tx_path = paste0(xenium_folder, '/transcripts.csv.gz')
feat_meta_path = paste0(xenium_folder, '/cell_feature_matrix/features-blank.tsv.gz')

# Files (AGGREGATE):
expr_mat_path = paste0(xenium_folder, '/cell_feature_matrix')
cell_meta_path = paste0(xenium_folder, '/cells.csv.gz') # contains spatlocs


## Load features metadata --------------------------------------------------
# Make sure cell_feature_matrix folder is unpacked
feature_dt = data.table::fread(feat_meta_path, header = FALSE)
colnames(feature_dt) = c('ensembl_ID', 'feat_name', 'feat_type')

# Find the feature IDs that belong to each feature type
feature_dt[, table(feat_type)]
feat_types = names(feature_dt[, table(feat_type)])
feat_types_IDs = lapply(feat_types, function(type)
  feature_dt[feat_type == type, unique(feat_name)])
names(feat_types_IDs) = feat_types


## Load transcript-level data ----------------------------------------------
tx_dt = data.table::fread(tx_path)
data.table::setnames(
  x = tx_dt,
  old = c('feature_name', 'x_location', 'y_location'),
  new = c('feat_ID', 'x', 'y')
)
cat(
  'Transcripts info available:\n ',
  paste0('"', colnames(tx_dt), '"'),
  '\n',
  'with',
  tx_dt[, .N],
  'unfiltered detections\n'
)

# Filter by qv (Phred score)
tx_dt_filtered = tx_dt[qv >= 20]
cat('and', tx_dt_filtered[, .N], 'filtered detections\n\n')

# Separate detections by feature type
tx_dt_types = lapply(feat_types_IDs, function(types)
  tx_dt_filtered[feat_ID %in% types])

invisible(lapply(seq_along(tx_dt_types), function(x) {
  cat(names(tx_dt_types)[[x]], ' detections: ', tx_dt_types[[x]][, .N], '\n')
}))


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


## Load polygon data -------------------------------------------------------
cellPoly_dt = data.table::fread(cell_bound_path)
nucPoly_dt = data.table::fread(nuc_bound_path)

data.table::setnames(
  cellPoly_dt,
  old = c('cell_id', 'vertex_x', 'vertex_y'),
  new = c('poly_ID', 'x', 'y')
)
data.table::setnames(
  nucPoly_dt,
  old = c('cell_id', 'vertex_x', 'vertex_y'),
  new = c('poly_ID', 'x', 'y')
)

gpoly_cells = createGiottoPolygonsFromDfr(segmdfr = cellPoly_dt,
                                          name = 'cell',
                                          calc_centroids = TRUE)
gpoly_nucs = createGiottoPolygonsFromDfr(segmdfr = nucPoly_dt,
                                         name = 'nucleus',
                                         calc_centroids = TRUE)


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

if (length(list.files(xenium_folder)) == 0) {
    message <- paste("Please move the raw-data for", region, "into the newly created raw-data directory:", xenium_folder)
    waitForInput(message = message) # custom function to wait for data transfer to finish
} else {
  cat(bold(magenta("There is already data in the raw-data folder. Proceeding with the script.\n")))
  # Continue with the rest of your script here
}


if (!file.exists(file.path(xenium_folder, "cell_feature_matrix", "features-blank.tsv.gz"))) {
  features <-
    read_tsv(here(xenium_folder, "cell_feature_matrix", "features.tsv.gz"))
  View(features)
  features[features == "Unassigned Codeword"] <- "Blank Codeword"
  write_tsv(
    x = features,
    file = here(
      xenium_folder,
      "cell_feature_matrix",
      "features-blank.tsv.gz"
    )
  )
  cat(bold(magenta("This package looks for Blank Codeword in the place of Unassigned Codeword, so first we have to rename them in the features.tsv.doc.\n")))
  cat(bold(magenta("It looks like you haven't run this code before. Blank codewords have now been successfully renamed and the file has been saved as features-blank.tsv.gz.\n")))
} else {
  cat(bold(magenta("This package looks for Blank Codeword in the place of Unassigned Codeword, so first we have to rename them from the features.tsv.doc.\n")))
  cat(bold(magenta("The file 'features-blank.tsv.gz' already exists. Skipping renaming.\n")))
}

# General files (some are supplemental files)
settings_path = paste0(xenium_folder, '/experiment.xenium')
he_img_path = paste0(xenium_folder, '/pyramidalhe.ome.tif')
panel_meta_path = paste0(xenium_folder, "/xenium_panel.tsv") # (optional)

# Files (SUBCELLULAR): 
cell_bound_path = paste0(xenium_folder, '/cell_boundaries.csv.gz')
nuc_bound_path = paste0(xenium_folder, '/nucleus_boundaries.csv.gz')
tx_path = paste0(xenium_folder, '/transcripts.csv.gz')
feat_meta_path = paste0(xenium_folder, '/cell_feature_matrix/features-blank.tsv.gz')

# Files (AGGREGATE):
expr_mat_path = paste0(xenium_folder, '/cell_feature_matrix')
cell_meta_path = paste0(xenium_folder, '/cells.csv.gz') # contains spatlocs


## Load features metadata --------------------------------------------------
# Make sure cell_feature_matrix folder is unpacked
feature_dt = data.table::fread(feat_meta_path, header = FALSE)
colnames(feature_dt) = c('ensembl_ID', 'feat_name', 'feat_type')

# Find the feature IDs that belong to each feature type
feature_dt[, table(feat_type)]
feat_types = names(feature_dt[, table(feat_type)])
feat_types_IDs = lapply(feat_types, function(type)
  feature_dt[feat_type == type, unique(feat_name)])
names(feat_types_IDs) = feat_types


## Load transcript-level data ----------------------------------------------
tx_dt = data.table::fread(tx_path)
data.table::setnames(
  x = tx_dt,
  old = c('feature_name', 'x_location', 'y_location'),
  new = c('feat_ID', 'x', 'y')
)
cat(
  'Transcripts info available:\n ',
  paste0('"', colnames(tx_dt), '"'),
  '\n',
  'with',
  tx_dt[, .N],
  'unfiltered detections\n'
)

# Filter by qv (Phred score)
tx_dt_filtered = tx_dt[qv >= 20]
cat('and', tx_dt_filtered[, .N], 'filtered detections\n\n')

# Separate detections by feature type
tx_dt_types = lapply(feat_types_IDs, function(types)
  tx_dt_filtered[feat_ID %in% types])

invisible(lapply(seq_along(tx_dt_types), function(x) {
  cat(names(tx_dt_types)[[x]], ' detections: ', tx_dt_types[[x]][, .N], '\n')
}))


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


## Load polygon data -------------------------------------------------------
cellPoly_dt = data.table::fread(cell_bound_path)
nucPoly_dt = data.table::fread(nuc_bound_path)

data.table::setnames(
  cellPoly_dt,
  old = c('cell_id', 'vertex_x', 'vertex_y'),
  new = c('poly_ID', 'x', 'y')
)
data.table::setnames(
  nucPoly_dt,
  old = c('cell_id', 'vertex_x', 'vertex_y'),
  new = c('poly_ID', 'x', 'y')
)

gpoly_cells = createGiottoPolygonsFromDfr(segmdfr = cellPoly_dt,
                                          name = 'cell',
                                          calc_centroids = TRUE)
gpoly_nucs = createGiottoPolygonsFromDfr(segmdfr = nucPoly_dt,
                                         name = 'nucleus',
                                         calc_centroids = TRUE)


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

if (length(list.files(xenium_folder)) == 0) {
  message <- paste("Please move the raw-data for", region, "into the newly created raw-data directory:", xenium_folder)
  waitForInput(message = message) # custom function to wait for data transfer to finish
} else {
  cat(bold(magenta("There is already data in the raw-data folder. Proceeding with the script.\n")))
  # Continue with the rest of your script here
}

if (!file.exists(file.path(xenium_folder, "cell_feature_matrix", "features-blank.tsv.gz"))) {
  features <-
    read_tsv(here(xenium_folder, "cell_feature_matrix", "features.tsv.gz"))
  View(features)
  features[features == "Unassigned Codeword"] <- "Blank Codeword"
  write_tsv(
    x = features,
    file = here(
      xenium_folder,
      "cell_feature_matrix",
      "features-blank.tsv.gz"
    )
  )
  cat(bold(magenta("This package looks for Blank Codeword in the place of Unassigned Codeword, so first we have to rename them in the features.tsv.doc.\n")))
  cat(bold(magenta("It looks like you haven't run this code before. Blank codewords have now been successfully renamed and the file has been saved as features-blank.tsv.gz.\n")))
} else {
  cat(bold(magenta("This package looks for Blank Codeword in the place of Unassigned Codeword, so first we have to rename them from the features.tsv.doc.\n")))
  cat(bold(magenta("The file 'features-blank.tsv.gz' already exists. Skipping renaming.\n")))
}


# General files (some are supplemental files)
settings_path = paste0(xenium_folder, '/experiment.xenium')
he_img_path = paste0(xenium_folder, '/pyramidalhe.ome.tif')
panel_meta_path = paste0(xenium_folder, "/xenium_panel.tsv") # (optional)

# Files (SUBCELLULAR): 
cell_bound_path = paste0(xenium_folder, '/cell_boundaries.csv.gz')
nuc_bound_path = paste0(xenium_folder, '/nucleus_boundaries.csv.gz')
tx_path = paste0(xenium_folder, '/transcripts.csv.gz')
feat_meta_path = paste0(xenium_folder, '/cell_feature_matrix/features-blank.tsv.gz')

# Files (AGGREGATE):
expr_mat_path = paste0(xenium_folder, '/cell_feature_matrix')
cell_meta_path = paste0(xenium_folder, '/cells.csv.gz') # contains spatlocs


## Load features metadata --------------------------------------------------
# Make sure cell_feature_matrix folder is unpacked
feature_dt = data.table::fread(feat_meta_path, header = FALSE)
colnames(feature_dt) = c('ensembl_ID', 'feat_name', 'feat_type')

# Find the feature IDs that belong to each feature type
feature_dt[, table(feat_type)]
feat_types = names(feature_dt[, table(feat_type)])
feat_types_IDs = lapply(feat_types, function(type)
  feature_dt[feat_type == type, unique(feat_name)])
names(feat_types_IDs) = feat_types


## Load transcript-level data ----------------------------------------------
tx_dt = data.table::fread(tx_path)
data.table::setnames(
  x = tx_dt,
  old = c('feature_name', 'x_location', 'y_location'),
  new = c('feat_ID', 'x', 'y')
)
cat(
  'Transcripts info available:\n ',
  paste0('"', colnames(tx_dt), '"'),
  '\n',
  'with',
  tx_dt[, .N],
  'unfiltered detections\n'
)

# Filter by qv (Phred score)
tx_dt_filtered = tx_dt[qv >= 20]
cat(bold(magenta('and', tx_dt_filtered[, .N], 'filtered detections\n\n')))

# Separate detections by feature type
tx_dt_types = lapply(feat_types_IDs, function(types)
  tx_dt_filtered[feat_ID %in% types])

invisible(lapply(seq_along(tx_dt_types), function(x) {
  cat(names(tx_dt_types)[[x]], ' detections: ', tx_dt_types[[x]][, .N], '\n')
}))


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


## Load polygon data -------------------------------------------------------
cellPoly_dt = data.table::fread(cell_bound_path)
nucPoly_dt = data.table::fread(nuc_bound_path)

data.table::setnames(
  cellPoly_dt,
  old = c('cell_id', 'vertex_x', 'vertex_y'),
  new = c('poly_ID', 'x', 'y')
)
data.table::setnames(
  nucPoly_dt,
  old = c('cell_id', 'vertex_x', 'vertex_y'),
  new = c('poly_ID', 'x', 'y')
)

gpoly_cells = createGiottoPolygonsFromDfr(segmdfr = cellPoly_dt,
                                          name = 'cell',
                                          calc_centroids = TRUE)
gpoly_nucs = createGiottoPolygonsFromDfr(segmdfr = nucPoly_dt,
                                         name = 'nucleus',
                                         calc_centroids = TRUE)


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

