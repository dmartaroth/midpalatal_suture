# ## ######################################## ## #
#                 LOAD PACKAGES                  #
# ## ######################################## ## #

# List of packages to check, install if not available, and load
packages_to_load <- c("here", "Giotto", "tidyverse", "crayon", "utils", "cowplot")

# Loop through each package
for (package in packages_to_load) {
  # Check if the package is installed
  if (!requireNamespace(package, quietly = TRUE)) {
    # If not installed, install it
    if (package == "devtools") {
      install.packages("devtools")
    } else if (package == "Giotto") {
      devtools::install_github("drieslab/Giotto@suite")
    } else {
      install.packages(package)
    }
  }
  # Load the package
  library(package, character.only = TRUE)
}