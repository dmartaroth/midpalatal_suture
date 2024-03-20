# ## ######################################## ## #
#                     DIRECTORIES                #
# ## ######################################## ## #

dir.create(here(home.path, region))
dir.create(results_folder <- here(home.path,region,"figs"))
dir.create(xenium_folder <-
             here(home.path, region, "raw-data")) # place raw data here
dir.create(output <- here(home.path, region, "data-output"))


if (length(list.files(xenium_folder)) == 0) {
  message <- paste("Please move the raw-data for", region, "into the newly created raw-data directory:", xenium_folder)
  waitForInput(message = message) # custom function to wait for data transfer to finish
} else {
  cat(bold(magenta("There is already data in the raw-data folder. Proceeding with the script.\n")))
  # Continue with the rest of your script here
}
