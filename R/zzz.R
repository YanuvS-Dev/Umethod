.onLoad <- function(libname, pkgname) {

  required_packages <- c("stringdist", "cowplot", "ggplot2", "reshape2",
                         "scales", "plotly", "lubridate", "Seurat", "svMisc","Matrix")

  # Install missing packages without prompting (using utils::install.packages)
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      packageStartupMessage(paste("Installing missing package:", pkg))
      utils::install.packages(pkg) # Use utils::install.packages
    }
  }
}
