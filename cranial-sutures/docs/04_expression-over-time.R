# ## ######################################## ## #
#                 EXPRESSION OVER TIME           #
# ## ######################################## ## #

# Date: Thu Apr 25 18:47:14 2024 ------------------

library(here)

source(here::here("cranial-sutures","docs","packages.R")) # load packages
source(here::here("cranial-sutures","docs","directories.R")) # load file paths/directories
source(here::here("cranial-sutures","docs","functions.R")) # load functions
source(here::here("cranial-sutures","docs","themes.R")) # load themes

# Load data
E16.fs <- readRDS(here::here("cranial-sutures/data-output/mes-subset_E16-fs.Rds"))
E18.fs <- readRDS(here::here("cranial-sutures/data-output/mes-subset_E18-fs.Rds"))
P10.fs <- readRDS(here::here("cranial-sutures/data-output/mes-subset_P10-fs.Rds"))
P28.fs <- readRDS(here::here("cranial-sutures/data-output/mes-subset_P28-fs.Rds"))

# Switch to 02_midline-mesenchyme-genes-in-fs.R in visium folder