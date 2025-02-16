.onLoad <- function(libname, pkgname) {
  options(repos = c(CRAN = "https://cloud.r-project.org/"))  # Set CRAN mirror

  # Required dependencies
  required_packages <- c("scCustomize", "stringdist", "cowplot", "ggplot2",
                         "reshape2", "scales", "plotly", "lubridate", "Seurat", "svMisc")

  # Ask user if they want to install dependencies (pauses execution)
  repeat {
    install_dependencies <- readline(prompt = "Do you want to install missing dependencies? (yes/no): ")
    if (tolower(install_dependencies) %in% c("yes", "no")) break
    message("Please enter 'yes' or 'no'.")
  }

  # Function to install missing packages
  install_missing_packages <- function(packages) {
    for (pkg in packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        message(paste("Installing missing package:", pkg))
        tryCatch({
          if (pkg == "scCustomize") {
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
  } else {
    message("Skipping dependency installation.")
  }

  # Always install Umethod (force reinstallation)
  message("Installing Umethod package...")
  devtools::install_github("YanuvS/Umethod", force = TRUE)
}
