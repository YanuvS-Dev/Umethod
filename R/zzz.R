.onLoad <- function(libname, pkgname) {

  options(repos = c(CRAN = "https://cloud.r-project.org/"))  # Set CRAN mirror

  # List of required packages
  required_packages <- c("scCustomize", "stringdist", "cowplot", "ggplot2",
                         "reshape2", "scales", "plotly", "lubridate", "Seurat", "svMisc")

  # Ask user if they want to install dependencies
  install_dependencies <- readline("Do you want to install missing dependencies? (yes/no): ")

  # Function to install missing packages
  install_missing_packages <- function(packages) {
    for (pkg in packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        message(paste("Installing missing package:", pkg))
        tryCatch({
          if (pkg == "scCustomize" && !requireNamespace(pkg, quietly = TRUE)) {
            remotes::install_github("samuel-marsh/scCustomize")
          } else {
            install.packages(pkg)
          }
        }, error = function(e) {
          message(paste("Failed to install package:", pkg))
        })
      } else {
        message(paste(pkg, "is already installed."))
      }
    }
  }

  # Install dependencies if user agrees
  if (tolower(install_dependencies) == "yes") {
    missing_packages <- required_packages[sapply(required_packages, function(pkg) !requireNamespace(pkg, quietly = TRUE))]
    if (length(missing_packages) > 0) {
      install_missing_packages(missing_packages)
    } else {
      message("All dependencies are already installed.")
    }
  }

  # Now install Umethod package if not already installed
  if (!requireNamespace("Umethod", quietly = TRUE)) {
    message("Installing Umethod package...")
    devtools::install_github("YanuvS/Umethod")
  } else {
    message("Umethod package is already installed.")
  }
}
