.onLoad <- function(libname, pkgname) {
  # Must-have packages
  required_packages <- c("Seurat", "Matrix", "reshape2", "scales")

  # Optional but recommended for downstream usage (e.g., plotting or reading parquet)
  optional_packages <- c("ggplot2", "cowplot", "arrow")

  # Check required
  missing_required <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_required) > 0) {
    packageStartupMessage("The following packages are REQUIRED by Umethod but not installed: ",
                          paste(missing_required, collapse = ", "))
  }

  # Check optional
  missing_optional <- optional_packages[!sapply(optional_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_optional) > 0) {
    packageStartupMessage("Optional packages missing (used for plotting or file I/O): ",
                          paste(missing_optional, collapse = ", "))
  }
}
