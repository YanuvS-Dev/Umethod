.onLoad <- function(libname, pkgname) {

  options(repos = c(CRAN = "https://cloud.r-project.org/"))  # Set CRAN mirror
  # Set the cache directory for pak
  if (Sys.getenv("R_USER_CACHE_DIR", unset = "") == "") {
    Sys.setenv(R_USER_CACHE_DIR = tempdir())  # Use temporary directory for checks
  }

  # List of required packages
  required_packages <- c("scCustomize", "stringdist", "cowplot", "ggplot2", "reshape2",
                         "scales", "plotly", "lubridate", "Seurat", "svMisc")

  # Step 1: Ask if to install `pak` package
  install_pak <- readline(prompt = "Do you want to install the pak package? (yes/no): ")
  if (tolower(install_pak) == "yes" && !requireNamespace("pak", quietly = TRUE)) {
    install.packages("pak")
    message("pak package installed.")
  }

  # Step 2: Ask if to install `scCustomize` package and choose installation method
  install_scCustomize <- readline(prompt = "Do you want to install scCustomize? (yes/no): ")
  if (tolower(install_scCustomize) == "yes") {
    install_method <- readline(prompt = "Choose method for scCustomize installation (pak/remotes/skip): ")

    if (tolower(install_method) == "pak" && !requireNamespace("scCustomize", quietly = TRUE)) {
      pak::pkg_install("samuel-marsh/scCustomize")
      message("scCustomize installed using pak.")
    } else if (tolower(install_method) == "remotes" && !requireNamespace("scCustomize", quietly = TRUE)) {
      remotes::install_github("samuel-marsh/scCustomize")
      message("scCustomize installed using remotes.")
    } else {
      message("Skipping scCustomize installation.")
    }
  }

  # Step 3: Install other missing dependencies using pak
  install_dependencies <- readline(prompt = "Do you want to install missing dependencies? (yes/no): ")
  if (tolower(install_dependencies) == "yes") {
    if (requireNamespace("pak", quietly = TRUE)) {
      for (pkg in required_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
          message(paste("Installing missing package:", pkg))
          pak::pkg_install(pkg)
        }
      }
    } else {
      message("pak is not installed. Skipping dependency installation.")
    }
  }

  # Step 4: Install Umethod package (always)
  message("Installing Umethod package...")
  devtools::install_github("YanuvS/Umethod", force = TRUE)
  message("Umethod package installed successfully.")

}
