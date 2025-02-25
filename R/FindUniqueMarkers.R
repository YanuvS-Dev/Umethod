
#' Find Unique Markers - U Method
#'
#' This function identifies the most unique markers for each cluster in a Seurat object.
#' It supports both numerical and character-based cluster labels. If the clusters are labeled as characters, set the `matchnames` parameter to `FALSE`.
#'
#' @param obj A Seurat object after clustering.
#' @param group_by The cluster labels in the metadata of `obj`.
#' @param p.threshold The p-value threshold for the U method. Defaults to a value of 0.05.
#' @param threshold The expression threshold above which a cell is considered positive for a gene. Default is 0.
#' @param P_in_thresh The minimum probability of expressing a gene inside the cluster for filtering the final results. Default keeps the full gene list.
#' @param P_out_thresh The maximum probability of expressing a gene in any other cluster for filtering the final results. Default keeps the full gene list.
#' @param varfeatures The features (genes) to use in the analysis. By default, the 2000 most variable genes will be used unless specified otherwise. If NULL, the method uses the 2000 most variable genes.
#' @param smallcluster A vector of cluster names to exclude from the analysis. This is useful for omitting small mixed clusters.
#' @param method The p-value adjustment method. Defaults to "BH" (Benjamini-Hochberg). Set to "none" for raw p-values or choose other methods available in the `p.adjust` function.
#'
#' @return A `data.frame` containing the most uniquely expressed genes per cluster that passed the filtering criteria.
#' @export

FindUniqueMarkers <- function(obj,
                              group_by,
                              p.threshold = 1,
                              threshold = 0,
                              P_in.thersh = 0,
                              P_out.thersh = 1,
                              varfeatures = NULL,
                              smallcluster = NULL,
                              method = "BH")
{
  # Determine expression matrix format
  expr_matrix <- if (class(obj[["RNA"]])[1] == "Assay5") {
    obj@assays$RNA@layers$counts
  }else
  {
    obj@assays$RNA@counts
  }

  groups <- obj@meta.data[[group_by]]

  if (is.null(varfeatures)) {
    varfeatureschoose <- rownames(obj)
  }

  # Binarize expression matrix based on threshold
  binary_matrix <- expr_matrix[row.names(obj) %in% varfeatureschoose, ] > threshold

  # Compute percentage of expressing cells per cluster
  if(any("factor" %in% class(groups)))
  {
    clusters <- levels(groups)

  }else
  {
    clusters <- unique(groups)
  }

  percent_stats <- sapply(clusters, function(cluster) {
    Matrix::rowMeans(binary_matrix[, groups == cluster]) * 100
  })

  rownames(percent_stats) <- varfeatureschoose
  colnames(percent_stats) <- clusters

  if(is.null(smallcluster)){
    ind_column <- 1:length(clusters)
  }else
  {
    ind_column <- c(1:length(clusters))[-which(colnames(percent_stats) %in% smallcluster)]
  }

  # Compute U-scores, Pin, and Pout
  Uscorestats <- t(apply(percent_stats[,ind_column], 1, function(x) {
    sorted_x <- sort(x / 100, decreasing = TRUE)
    x / 100 - ifelse(sorted_x[1] == x / 100, sorted_x[2], sorted_x[1])
  }))

  Pout <- t(apply(percent_stats[,ind_column], 1, function(x) {
    sorted_x <- sort(x / 100, decreasing = TRUE)
    ifelse(sorted_x[1] == x / 100, sorted_x[2], sorted_x[1])
  }))

  Pin <- t(apply(percent_stats, 1, function(x) x / 100))

  # Z-score transformation of U-scores
  Uscoresz <- apply(Uscorestats, 2, function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  })

  # Adjusted p-values
  UscoreszPvalue <- apply(Uscoresz, 2, function(x) {
    p.adjust(1 - pnorm(x), method = method)
  })

  # Compile results into a data frame
  l_out <- do.call(rbind, lapply(colnames(Uscorestats), function(i) {
    data.frame(
      Gene = rownames(Uscorestats),
      Cluster = i,
      Uscore = Uscorestats[, i],
      adj.p.value = UscoreszPvalue[, i],
      P_in = Pin[, i],
      P_out = Pout[, i]
    )
  }))

  # Filter based on thresholds
  l_out <- l_out[l_out$adj.p.value < p.threshold & l_out$P_in >= P_in.thersh & l_out$P_out <= P_out.thersh, ]
  l_out <- l_out[order(l_out$Cluster, -l_out$Uscore), ]
  return(l_out)
}
