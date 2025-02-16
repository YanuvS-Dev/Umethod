.onLoad <- function(libname, pkgname) {

  options(repos = c(CRAN = "https://cloud.r-project.org/"))  # Set CRAN mirror

  # List of required packages
  required_packages <- c("scCustomize", "stringdist", "cowplot", "ggplot2",
                         "reshape2", "scales", "plotly", "lubridate", "Seurat", "svMisc")

  # Check if each package is installed, and install missing ones
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing missing package:", pkg))
      install.packages(pkg)  # Install the missing package from CRAN
    } else {
      message(paste(pkg, "is already installed. Skipping installation."))
    }
  }

  # If any of the packages are not available from CRAN (e.g., GitHub versions),
  # you can add custom installation commands for them here.
  if (!requireNamespace("scCustomize", quietly = TRUE)) {
    message("Installing scCustomize from GitHub")
    remotes::install_github("samuel-marsh/scCustomize")
  }

}
