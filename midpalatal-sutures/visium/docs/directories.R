# ## ######################################## ## #
#                     DIRECTORIES                #
# ## ######################################## ## #
dir.create(home.path <- here::here("midpalatal-sutures"))
dir.create(here(home.path, age))
dir.create(results_folder <- here(home.path,age,"figs"))
dir.create(visium_folder <-
             here(home.path, age, "raw-data")) # place raw data here
dir.create(output <- here(home.path, age, "data-output"))


if (length(list.files(visium_folder)) == 0) {
  message <- paste("Please move the raw-data for", age, "into the newly created raw-data directory:", visium_folder)
  waitForInput(message = message) # custom function to wait for data transfer to finish
} else {
  cat(bold(magenta("There is already data in the raw-data folder. Proceeding with the script.\n")))
  # Continue with the rest of your script here
}
