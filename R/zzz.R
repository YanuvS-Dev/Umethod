.onLoad <- function(libname, pkgname) {
  required_packages <- c("stringdist", "cowplot", "ggplot2", "reshape2","arrow",
                         "scales", "plotly", "lubridate", "Seurat", "svMisc", "Matrix")

  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

  if (length(missing_packages) > 0) {
    packageStartupMessage("Please install the following missing packages: ", paste(missing_packages, collapse = ", "))
  }
}
