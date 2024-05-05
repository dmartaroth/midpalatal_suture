# ## ######################################## ## #
#                 LOAD PACKAGES                  #
# ## ######################################## ## #

# List of packages to check and load
packages_to_load <- c("dplyr", "Seurat", "cowplot", "ggplot2", "readr", 
                      "tidyverse", "hdf5r", "patchwork", "qs", "scCustomize")

# Loop through each package
for (package in packages_to_load) {
  # Check if the package is installed
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)  # Install the package if not installed
  }
  # Load the package
  library(package, character.only = TRUE)
}
