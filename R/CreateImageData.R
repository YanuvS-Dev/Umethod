#' U-Method: CreateImageData Function
#'
#' This function extracts read counts from the Visium HD dataset for each spatial spot and selected markers, integrating them into a Seurat object.
#'
#' @param bc_matrix A 10x filtered gene expression matrix, loadable via Read10X(".../Visium_HD__binned_outputsXum/filtered_feature_bc_matrix/"). X can be 2, 8, or 16.
#' @param poaraq A parquet file containing spatial barcode metadata, loadable via read_parquet(".../Visium_HD__binned_outputX2um/spatial/tissue_positions.parquet").
#' @param markers A vector of marker gene names to visualize on Visium HD data.
#' @param print.table Logical; if TRUE, prints a table showing the number of spots per marker.
#' @param project A string specifying the project name for the Seurat object. Defaults to "Visium & U-Method Project" if NULL.
#' @param Signituretable A table of marker names with cluster column names (genesetlong). Defaults to NULL.
#' @return A Seurat object containing expression data and spatial metadata for the selected markers.
#' @export

CreateImageData <- function(bc_matrix, poaraq, markers, project = NULL, print.table = FALSE, Signituretable = NULL) {

  # Create Seurat object from Visium HD output
  seurat_object <- CreateSeuratObject(
    counts = bc_matrix,
    assay = "Spatial",
    project = ifelse(is.null(project), "Visium & U-Method Project", project)
  )

  # Extract marker expression and format into dataframe
  marker_expr <- bc_matrix[rownames(bc_matrix) %in% markers, , drop = FALSE]
  df <- as.data.frame(t(marker_expr))
  df$cell_id <- rownames(df)

  # Merge spatial metadata with expression data
  cell_metadata <- poaraq %>%
    filter(barcode %in% colnames(seurat_object)) %>%
    rename(cell_id = barcode) %>%
    left_join(df, by = "cell_id")

  # Add metadata to Seurat object
  seurat_object <- AddMetaData(seurat_object, metadata = cell_metadata)

  # Classify cells based on marker expression
  seurat_object$Class <- apply(cell_metadata[, markers, drop = FALSE], 1, function(x) {
    detected_markers <- markers[x > 0]
    if (length(detected_markers) == 0) NA else paste(detected_markers, collapse = "&")
  })

  # Convert Class to factor
  seurat_object$Class <- factor(seurat_object$Class)

  # Compute overlap index (number of expressed markers per spot)
  seurat_object$ind_overlap <- rowSums(cell_metadata[, markers, drop = FALSE] > 0)

  # Print marker classification table if requested
  if (print.table) {
    print(table(seurat_object$Class, useNA = "ifany"))
  }

  return(seurat_object)
}
