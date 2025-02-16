.onLoad <- function(libname, pkgname) {

  options(repos = c(CRAN = "https://cloud.r-project.org/"))  # Set CRAN mirror
  # Set the cache directory for pak
  if (Sys.getenv("R_USER_CACHE_DIR", unset = "") == "") {
    Sys.setenv(R_USER_CACHE_DIR = tempdir())  # Use temporary directory for checks
  }

  # List of required packages
  required_packages <- c("scCustomize", "stringdist", "cowplot", "ggplot2", "reshape2",
                         "scales", "plotly", "lubridate", "Seurat", "svMisc")

  # Check if pak is installed
  if (!requireNamespace("pak", quietly = TRUE)) {
    install_pak <- readline(prompt = "Do you want to install the pak package? (yes/no): ")
    if (tolower(install_pak) == "yes") {
      install.packages("pak")
      message("pak package installed.")
    }
  }

  # Check if scCustomize is installed
  if (!requireNamespace("scCustomize", quietly = TRUE)) {
    install_scCustomize <- readline(prompt = "Do you want to install scCustomize? (yes/no): ")
    if (tolower(install_scCustomize) == "yes") {
      install_method <- readline(prompt = "Choose method for scCustomize installation (pak/remotes/skip): ")

      if (tolower(install_method) == "pak") {
        pak::pkg_install("samuel-marsh/scCustomize")
        message("scCustomize installed using pak.")
      } else if (tolower(install_method) == "remotes") {
        remotes::install_github("samuel-marsh/scCustomize")
        message("scCustomize installed using remotes.")
      } else {
        message("Skipping scCustomize installation.")
      }
    }
  }

  # Ask if to install missing dependencies
  install_dependencies <- readline(prompt = "Do you want to install missing dependencies? (yes/no): ")
  if (tolower(install_dependencies) == "yes") {
    for (pkg in required_packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        message(paste("Installing missing package:", pkg))
        pak::pkg_install(pkg)
      }
    }
  } else {
    message("Skipping installation of missing dependencies.")
  }

  # Check if Umethod is installed and install if not
  if (!requireNamespace("Umethod", quietly = TRUE)) {
    message("Umethod package is not installed. Installing Umethod...")
    devtools::install_github("YanuvS/Umethod", force = TRUE)
    message("Umethod package installed successfully.")
  } else {
    message("Umethod package is already installed.")
  }
}
