.onLoad <- function(libname, pkgname) {

  options(repos = c(CRAN = "https://cloud.r-project.org/"))  # Set CRAN mirror

  # Set the cache directory for pak (if necessary, though now we're not using pak)
  if (Sys.getenv("R_USER_CACHE_DIR", unset = "") == "") {
    Sys.setenv(R_USER_CACHE_DIR = tempdir())  # Use temporary directory for checks
  }

  # List of required packages
  required_packages <- c("scCustomize", "stringdist", "cowplot", "ggplot2",
                         "reshape2", "scales", "plotly", "lubridate", "Seurat", "svMisc")

  # Function to install a package if it is not already installed
  install_if_missing <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing missing package:", pkg))
      install.packages(pkg)
    }
  }

  # Install missing required packages from CRAN
  for (pkg in required_packages) {
    install_if_missing(pkg)
  }

  # For GitHub packages (e.g., scCustomize), use remotes if necessary
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }

  # Install scCustomize from GitHub if missing
  if (!requireNamespace("scCustomize", quietly = TRUE)) {
    message("Installing scCustomize from GitHub")
    remotes::install_github("user/scCustomize")
  }
}
