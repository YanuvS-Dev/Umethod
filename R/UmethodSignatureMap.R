#' UmethodSignatureMap: Visualize Signature-Level Spatial Classification
#'
#' Computes average expression of top marker sets per cluster and identifies the top-scoring signature per pixel.
#'
#' @param seurat_object A Seurat object containing expression and spatial metadata.
#' @param gene_set A matrix of top marker genes per cluster (from \link{FindUniqueMarkers} output, element \code{[[2]]}).
#' @param color_palette Optional color vector; defaults to hue-based palette sized to non-empty clusters.
#' @param statistics Logical; if TRUE, prints overlap statistics.
#'
#' @return A list with:
#'   - highest_score_data: one signature per pixel (max expression)
#'   - signatureLong: long-format melted data for all signatures
#'   - Classlist: names of signatures
#'   - signature_colors: named color vector
#' @export
UmethodSignatureMap <- function(seurat_object,
                                gene_set,
                                color_palette = NULL,
                                statistics = TRUE) {

  # Drop empty columns (where all markers are NA)
  gene_set_clean <- gene_set[, apply(gene_set, 2, function(col) any(!is.na(col))), drop = FALSE]

  # Compute average expression per signature
  signature_df <- as.data.frame(matrix(NA, nrow = nrow(seurat_object@meta.data), ncol = ncol(gene_set_clean)))

  for (i in seq_len(ncol(gene_set_clean))) {
    marker_genes <- gene_set_clean[, i]
    marker_genes <- marker_genes[!is.na(marker_genes)]

    expr_data <- seurat_object@meta.data[, colnames(seurat_object@meta.data) %in% marker_genes, drop = FALSE]

    if (is.null(dim(expr_data))) {
      signature_df[, i] <- expr_data
      if (statistics) {
        cat("Only used", colnames(seurat_object@meta.data)[colnames(seurat_object@meta.data) %in% marker_genes],
            "for", colnames(gene_set_clean)[i], "\n")
      }
    } else {
      signature_df[, i] <- apply(expr_data, 1, mean, na.rm = TRUE)
    }

    if (statistics) {
      cat("Finished average expression calculation of", colnames(gene_set_clean)[i], "\n")
    }
  }

  colnames(signature_df) <- paste("Signature", colnames(gene_set_clean), sep = ".")
  rownames(signature_df) <- rownames(seurat_object@meta.data)

  # Add to metadata
  seurat_object <- AddMetaData(seurat_object, signature_df)

  # Build beforemelt data frame
  beforemelt <- seurat_object@meta.data[, c("pxl_row_in_fullres", "pxl_col_in_fullres", colnames(signature_df))]

  if (statistics) {
    signature_overlap_count <- apply(beforemelt[, grep(names(beforemelt), pattern = "Signature")], 1, function(x) {
      length(which(x > 0))
    })
    if (sum(signature_overlap_count > 0) > 0) {
      ratio <- sum(signature_overlap_count == 1) / sum(signature_overlap_count > 0)
      cat(round(ratio, 3), "Unique spot probability (non-overlapping)\n")
    }
  }

  # Melt into long format
  signatureLong <- reshape2::melt(beforemelt, id.vars = c("pxl_row_in_fullres", "pxl_col_in_fullres"))
  colnames(signatureLong)[3] <- "Class"
  signatureLong$value <- ifelse(signatureLong$value == 0, NA, signatureLong$value)

  # Top signature per pixel (base R)
  valid_rows <- !is.na(signatureLong$value)
  filtered <- signatureLong[valid_rows, ]

  pixel_ids <- paste(filtered$pxl_col_in_fullres, filtered$pxl_row_in_fullres, sep = "_")
  split_list <- split(seq_len(nrow(filtered)), pixel_ids)

  max_rows <- unlist(lapply(split_list, function(idxs) {
    vals <- filtered$value[idxs]
    idxs[which.max(vals)]
  }))

  highest_score_data <- filtered[max_rows, , drop = FALSE]

  # Assign color palette
  if (is.null(color_palette)) {
    color_palette <- scales::hue_pal()(ncol(gene_set_clean))
  }

  signature_colors <- color_palette
  names(signature_colors) <- colnames(signature_df)
  signature_colors <- signature_colors[!is.na(names(signature_colors))]

  highest_score_data$Class <- factor(highest_score_data$Class, levels = names(signature_colors))

  return(list(
    highest_score_data = highest_score_data,
    signatureLong = signatureLong,
    Classlist = colnames(signature_df),
    signature_colors = signature_colors
  ))
}
