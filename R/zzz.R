.onLoad <- function(libname, pkgname) {

  options(repos = c(CRAN = "https://cloud.r-project.org/"))  # Set CRAN mirror

  # Set the cache directory for pak
  if (Sys.getenv("R_USER_CACHE_DIR", unset = "") == "") {
    Sys.setenv(R_USER_CACHE_DIR = tempdir())  # Use temporary directory for checks
  }

  # List of required packages
  required_packages <- c("scCustomize", "stringdist", "cowplot", "ggplot2",
                         "reshape2", "scales", "plotly", "lubridate", "Seurat", "svMisc")

  # Function to conditionally install and load 'pak'
  ensure_pak <- function() {
    if (!requireNamespace("pak", quietly = TRUE)) {
      if (interactive()) {
        answer <- readline("The 'pak' package is not installed. Do you want to install it now? (yes/no): ")
        if (tolower(answer) == "yes") {
          install.packages("pak")
        } else {
          message("Skipping 'pak' installation.")
          return(FALSE)
        }
      } else {
        message("Skipping 'pak' installation (non-interactive mode).")
        return(FALSE)
      }
    }
    library(pak, quietly = TRUE)  # Load 'pak' if installed
    return(TRUE)
  }

  # Install missing required packages
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing missing package:", pkg))
      if (ensure_pak()) {
        pak::pkg_install(pkg)
      } else {
        message(paste("Skipping installation of", pkg, "because 'pak' is missing."))
      }
    }
  }
}
