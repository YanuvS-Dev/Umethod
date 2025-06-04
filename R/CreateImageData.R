#' U-Method: CreateImageData Function
#'
#' Creates a Seurat object from Visium HD data, optionally integrating marker expression.
#'
#' @param counts_matrix A 10x filtered gene expression matrix, loadable via Read10X().
#' @param poaraq A data.frame from tissue_positions.parquet with spatial metadata.
#' @param markers Optional vector of gene names to visualize. If NULL, no marker metadata is added.
#' @param print.table Logical; if TRUE, prints a table of marker combinations per spot.
#' @param project A string for the Seurat object project name. Default is "Visium & U-Method Project".
#'
#' @return A Seurat object with spatial metadata and optionally marker classifications.
#' @export
CreateImageData <- function(counts_matrix,
                            poaraq,
                            markers = NULL,
                            project = NULL,
                            print.table = FALSE) {

  # Create base Seurat object
  seurat_object <- CreateSeuratObject(
    counts = counts_matrix,
    assay = "Spatial",
    project = ifelse(is.null(project), "Visium & U-Method Project", project)
  )

  # Align metadata using match for ordered merge
  cell_metadata <- poaraq[match(colnames(seurat_object), poaraq$barcode), ]
  colnames(cell_metadata)[colnames(cell_metadata) == "barcode"] <- "cell_id"
  rownames(cell_metadata) <- cell_metadata$cell_id

  # If markers are provided, apply correction logic
  if (!is.null(markers)) {
    # Generalize .number suffix cleanup (e.g., ECSCR.1 -> ECSCR)
    suffix_match <- grep("\\.[0-9]+$", markers)
    if (length(suffix_match) > 0) {
      for (i in suffix_match) {
        base <- sub("\\.[0-9]+$", "", markers[i])
        if (base %in% rownames(counts_matrix)) {
          markers[i] <- base
        } else {
          markers[i] <- NA
        }
      }
    }

    # Final filtering to valid genes in the matrix
    markers <- na.omit(markers)
    markers <- markers[markers %in% rownames(counts_matrix)]

    if (length(markers) > 0) {
      # Extract marker expression and transpose into cell x gene
      marker_expr <- counts_matrix[markers, , drop = FALSE]
      marker_df <- as.data.frame(Matrix::t(marker_expr))
      marker_df$cell_id <- rownames(marker_df)

      # Merge expression into spatial metadata
      cell_metadata <- merge(cell_metadata, marker_df, by = "cell_id", all.x = TRUE)
      rownames(cell_metadata) <- cell_metadata$cell_id

      # Add all metadata
      seurat_object <- AddMetaData(seurat_object, metadata = cell_metadata)

      # Marker-based classification
      seurat_object$Class <- apply(cell_metadata[, markers, drop = FALSE], 1, function(x) {
        detected <- markers[x > 0]
        if (length(detected) == 0) NA else paste(detected, collapse = "&")
      })
      seurat_object$Class <- factor(seurat_object$Class)

      if (print.table) {
        print(table(seurat_object$Class, useNA = "ifany"))
      }
    }
  }

  return(seurat_object)
}
