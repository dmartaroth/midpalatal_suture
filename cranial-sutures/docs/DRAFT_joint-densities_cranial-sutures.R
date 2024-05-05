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


# Load data
E16.fs <- readRDS(here::here("cranial-sutures/E16-fs/data-output/annot_E16-fs.Rds"))
E18.fs <- readRDS(here::here("cranial-sutures/E18-fs/data-output/annot_E18-fs.Rds"))
P10.fs <- readRDS(here::here("cranial-sutures/P10-fs/data-output/annot_P10-fs.Rds"))
P28.fs <- readRDS(here::here("cranial-sutures/P28-fs/data-output/annot_P28-fs.Rds"))





p1 <- FeaturePlot(E16.fs,features=c("Sfrp2","Tnn"),blend=TRUE,pt.size=1.5,
            label=TRUE,label.size = 3,label.color = "white",alpha = 0.5,cols = c("navyblue","red","green"),order =  TRUE)
p2 <- FeaturePlot(E18.fs,features=c("Sfrp2","Tnn"),blend=TRUE,pt.size=1.5,
                  label=TRUE,label.size = 3,label.color = "white",alpha = 0.5,cols = c("navyblue","red","green"),order =  TRUE)
p3 <- FeaturePlot(P10.fs,features=c("Sfrp2","Tnn"),blend=TRUE,pt.size=1.5,
                  label=TRUE,label.size = 3,label.color = "white",alpha = 0.5,cols = c("navyblue","red","green"),order =  TRUE)
p4<-FeaturePlot(P28.fs,features=c("Sfrp2","Tnn"),blend=TRUE,pt.size=1.5,
                label=TRUE,label.size = 3,label.color = "white",alpha = 0.5,cols = c("navyblue","red","green"),order =  TRUE)

(plots <- plot_grid(p1,p2,p3,p4,ncol = 1))

convenient_save_plot(
  plots,
  "joint-density_Sfrp2-Tnn",
  number = 01,
  width = 11,
  height = 12,
  dir = results_folder,
  file_format = "pdf"
)
