.onLoad <- function(libname, pkgname) {

  options(repos = c(CRAN = "https://cloud.r-project.org/"))  # Set CRAN mirror
  # Set the cache directory for pak
  if (Sys.getenv("R_USER_CACHE_DIR", unset = "") == "") {
    Sys.setenv(R_USER_CACHE_DIR = tempdir())  # Use temporary directory for checks
  }
  # List of required packages
  required_packages <- c("scCustomize", "stringdist", "cowplot", "ggplot2", "reshape2",
                         "scales", "plotly", "lubridate", "Seurat","svMisc")

  # Install pak if not already installed
  if (!requireNamespace("pak", quietly = TRUE)) {
    install.packages("pak")
  }

  # Install missing required packages using pak
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing missing package:", pkg))
      pak::pkg_install(pkg)
    }
  }
}
