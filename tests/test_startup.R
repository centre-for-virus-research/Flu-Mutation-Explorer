
# Validates that the Shiny app can be loaded without errors
tryCatch({
  message("Attempting to load app.R...")
  # Source global.R and R/ directory to ensure environment is set up like Shiny
  source("global.R")
  sapply(list.files("R", full.names = TRUE), source)

  # Source the app file. Run this script from the project root.
  app <- source("app.R", chdir = TRUE)$value
  
  if (inherits(app, "shiny.appobj")) {
    message("SUCCESS: App loaded successfully!")
  } else {
    stop("app.R did not return a shiny app object")
  }
}, error = function(e) {
  message("FAILURE: Error loading app")
  message(e$message)
  quit(status = 1)
})
