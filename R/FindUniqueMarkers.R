
#' U-Method: Find Unique Markers
#'
#' Identifies the most distinct markers for each cluster in a Seurat object using the U-method.
#' Supports both numeric and character-based cluster labels. If cluster labels are character-based,
#'
#' @param obj A Seurat object with cluster assignments.
#' @param group_by The metadata column containing cluster labels.
#' @param p.threshold The adjusted p-value threshold for selecting unique markers. Default is 1.
#' @param threshold The expression cutoff above which a gene is considered expressed in a cell. Default is 0.
#' @param P_in_thresh Minimum probability of gene expression within a cluster for inclusion. Default is 0.
#' @param P_out_thresh Maximum probability of gene expression in other clusters for filtering. Default is 1.
#' @param varfeatures A vector of genes to include in the analysis. Defaults to the 2000 most variable genes if NULL.
#' @param omitCluster A vector of cluster names to exclude, useful for omitting small or mixed clusters.
#' @param method P-value adjustment method. Defaults to "BH" (Benjamini-Hochberg). Other options follow `p.adjust` methods.
#'
#' @return A `data.frame` listing uniquely expressed genes per cluster, filtered based on thresholds.
#' @export

FindUniqueMarkers <- function(obj,
                              group_by,
                              p.threshold = 1,
                              threshold = 0,
                              P_in.thersh = 0,
                              P_out.thersh = 1,
                              varfeatures = NULL,
                              omitCluster = NULL,
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

  if(is.null(omitCluster)){
    ind_column <- 1:length(clusters)
  }else
  {
    ind_column <- c(1:length(clusters))[-which(colnames(percent_stats) %in% omitCluster)]
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
