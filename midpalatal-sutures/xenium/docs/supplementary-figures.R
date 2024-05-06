# ## ######################################## ## #
#                 SUPPLEMENTARY FIGURES          #
# ## ######################################## ## #

# Date: Mon May 06 12:59:11 2024 ------------------


# Xenium setup ------------------------------------------------------------

# Load required libraries ----------------------------------------------------------
source("midpalatal-sutures/xenium/docs/packages.R")
source("midpalatal-sutures/xenium/docs/themes.R")
source("midpalatal-sutures/xenium/docs/functions.R")


# Set up Giotto environment -----------------------------------------------
# Set Giotto python path
python_path = NULL
if(is.null(python_path)) {
  installGiottoEnvironment()
}


# Figure 2 ----------------------------------------------------------------

# Spatial expression of Fgfr1/2/3, Twist1 in other planes of sectioning at e15.5
here()
xenium_folder <- here::here("midpalatal-sutures/xenium")
dir.create(suppl_data <-
             here(xenium_folder, "pubfigs","supplementary"),
           recursive = TRUE)

## e15a (anterior) ---------------------------------------------------------


e15a <- loadGiotto(here::here(xenium_folder,"region-5/e15a/data-output/e15a_giotto_objects/gobject.RDS"))
my_colors <- c("darkolivegreen2","dodgerblue3","red2","goldenrod1","orange","mediumpurple2","pink","dodgerblue","mediumorchid","mediumpurple","mintcream","blue3")
feat_colors <- c("green3","violet","black","dodgerblue")

genes <-  list('rna'=c("Fgfr1","Fgfr2","Fgfr3","Twist1"))
spatInSituPlotPoints(gobject = e15a,
                     show_image = FALSE,
                     feats = genes,
                     feats_color_code = feat_colors,
                     point_size = 1.5,show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     show_legend = TRUE,
                     polygon_alpha = 0.1,
                     polygon_color = 'pink2',
                     background_color = "floralwhite",
                     axis_text = 8,
                     axis_title = 9,
                     polygon_line_size = 0.1,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = TRUE,
                     coord_fix_ratio = TRUE,
                     polygon_fill_code = colorcode,
                     save_param = list(
                       save_name = paste0("S2_e15a_leidenclus"),
                       save_dir = suppl_data),return_plot=T)



# Figure 4 ----------------------------------------------------------------

# Spatial expression of Mkx, Chodl, Col3a1, Lum, Six2 in other planes of sectioning at e15.5


## e15a (anterior) ---------------------------------------------------------

my_colors <- c("darkolivegreen2","dodgerblue3","red2","goldenrod1","orange","mediumpurple2","pink","dodgerblue","mediumorchid","mediumpurple","mintcream","blue3")
feat_colors <- c("green3","violet","black","turquoise","red")
genes <-  list('rna'=c("Mkx","Chodl","Col3a1","Lum","Six2"))
spatInSituPlotPoints(gobject = e15a,
                     show_image = FALSE,
                     feats = genes,
                     feats_color_code = feat_colors,
                     point_size = 1.5,show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     show_legend = TRUE,
                     polygon_alpha = 0.1,
                     polygon_color = 'pink2',
                     background_color = "floralwhite",
                     axis_text = 8,
                     axis_title = 9,
                     polygon_line_size = 0.1,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = TRUE,
                     coord_fix_ratio = TRUE,
                     polygon_fill_code = colorcode,
                     save_param = list(
                       save_name = paste0("S4_e15a_leidenclus"),
                       save_dir = suppl_data),return_plot=T)

spatInSituPlotDensity(e15a,
                      feats = c("Lum","Six2","Mkx", "Col3a1"),
                      feat_type = "rna",
                      polygon_alpha = 1,
                      polygon_color = "white",
                      background_color = "white",
                      cow_n_col = 2,
                      save_param = list(
                        save_name = paste0("S4_e15a_Tenocyte_density"),
                        save_dir = suppl_data,save_format = "pdf"))


# Figure 5 ----------------------------------------------------------------

# Cluster markers for Visium section WT1R1, generated in
# 01_import-visium_findmarkers_save.R
# read.csv(here::here("midpalatal-sutures/visium/e15/data-output/WT1R1_all_markers_pct.csv"))
# Pdfs of spatial expression of top genes in "midpalatal-sutures/visium/e15/figs"
# Merged pdf added to supplementary data folder


# Annotated UMAPs for mes or other frontal suture datasets generated in 03_subset-mes.R
# Transferred plot_grid to this folder


# Figure 6 ----------------------------------------------------------------
e15c <- loadGiotto(here::here("midpalatal-sutures/xenium/region-7/e15c/data-output/e15c_giotto_objects"))


spatInSituPlotDensity(e15c,
                      feats = c("Bmp2","Bmp3","Bmp5","Bmp6","Dlx3","Smad9"),
                      feat_type = "rna",
                      polygon_alpha = 0.9,
                      polygon_color = "white",
                      background_color = "white",
                      cow_n_col = 2,
                      save_param = list(
                        save_name = paste0("S6_e15c_Bmp_density"),
                        save_dir = suppl_data,save_format = "pdf"))

spatInSituPlotDensity(e15c,
                      feats = c("Fgf2","Fgfr2","Fgfr3"),
                      feat_type = "rna",
                      polygon_alpha = 0.9,
                      polygon_color = "white",
                      background_color = "white",
                      cow_n_col = 3,
                      save_param = list(
                        save_name = paste0("S6_e15c_Fgf_density"),
                        save_dir = suppl_data,save_format = "pdf"))

spatInSituPlotDensity(e15c,
                      feats = c("Gli1","Kif7","Ptch1","Ptch2","Smo","Sufu"),
                      feat_type = "rna",
                      polygon_alpha = 0.9,
                      polygon_color = "white",
                      background_color = "white",
                      cow_n_col = 2,
                      save_param = list(
                        save_name = paste0("S6_e15c_Hh_density"),
                        save_dir = suppl_data,save_format = "pdf"))

spatInSituPlotDensity(e15c,
                      feats = c("Igf1","Igf2","Igfbp4"),
                      feat_type = "rna",
                      polygon_alpha = 0.9,
                      polygon_color = "white",
                      background_color = "white",
                      cow_n_col = 3,
                      save_param = list(
                        save_name = paste0("S6_e15c_Igf_density"),
                        save_dir = suppl_data,save_format = "pdf"))

spatInSituPlotDensity(e15c,
                      feats = c("Mmp2","Mmp9","Mmp15","Mmp23","Mmp14","Timp1","Timp2"),
                      feat_type = "rna",
                      polygon_alpha = 0.9,
                      polygon_color = "white",
                      background_color = "white",
                      cow_n_col = 2,
                      save_param = list(
                        save_name = paste0("S6_e15c_Mmp_density"),
                        save_dir = suppl_data,save_format = "pdf"))

spatInSituPlotDensity(e15c,
                      feats = c("Ltbp2","Smad2","Smad3","Smad6","Tgfb1","Tgfb2","Tgfb3","Tgfbi","Serpine1"),
                      feat_type = "rna",
                      polygon_alpha = 0.9,
                      polygon_color = "white",
                      background_color = "white",
                      cow_n_col = 3,
                      save_param = list(
                        save_name = paste0("S6_e15c_Tgfb_density"),
                        save_dir = suppl_data,save_format = "pdf"))

spatInSituPlotDensity(e15c,
                      feats = c("Axin2","Dkk2","Frzb","Lrp5","Pax9","Rspo1","Sfrp2","Tnik","Wif1","Lef1","Wnt5a"),
                      feat_type = "rna",
                      polygon_alpha = 0.9,
                      polygon_color = "white",
                      background_color = "white",
                      cow_n_col = 3,
                      save_param = list(
                        save_name = paste0("S6_e15c_Wnt_density"),
                        save_dir = suppl_data,save_format = "pdf"))
