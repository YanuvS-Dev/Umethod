.onLoad <- function(libname, pkgname) {

  # List of required packages
  required_packages <- c("scCustomize", "stringdist", "cowplot", "ggplot2", "reshape2",
                         "scales", "plotly", "lubridate", "Seurat", "svMisc")

  # Prevent multiple prompts during package installation
  already_installed <- FALSE

  # Check if pak is installed
  if (!requireNamespace("pak", quietly = TRUE)) {
    install_pak <- readline(prompt = "Do you want to install the pak package? (yes/no): ")
    if (tolower(install_pak) == "yes") {
      install.packages("pak")
      packageStartupMessage("pak package installed.")
    }
  }

  # Check if scCustomize is installed
  if (!requireNamespace("scCustomize", quietly = TRUE)) {
    install_scCustomize <- readline(prompt = "Do you want to install scCustomize? (yes/no): ")
    if (tolower(install_scCustomize) == "yes") {
      install_method <- readline(prompt = "Choose method for scCustomize installation (pak/remotes/skip): ")

      if (tolower(install_method) == "pak") {
        pak::pkg_install("samuel-marsh/scCustomize")
        packageStartupMessage("scCustomize installed using pak.")
      } else if (tolower(install_method) == "remotes") {
        remotes::install_github("samuel-marsh/scCustomize")
        packageStartupMessage("scCustomize installed using remotes.")
      } else {
        packageStartupMessage("Skipping scCustomize installation.")
      }
    }
  }

  # Ask if to install missing dependencies
  install_dependencies <- readline(prompt = "Do you want to install missing dependencies? (yes/no): ")
  if (tolower(install_dependencies) == "yes") {
    for (pkg in required_packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        packageStartupMessage(paste("Installing missing package:", pkg))
        pak::pkg_install(pkg)
      }
    }
  } else {
    packageStartupMessage("Skipping installation of missing dependencies.")
  }

}
