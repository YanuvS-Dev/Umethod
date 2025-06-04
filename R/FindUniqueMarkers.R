#' U-Method: Find Unique Markers
#'
#' Identifies the most distinct markers for each cluster in a Seurat object using the U-method.
#'
#' @param obj A Seurat object with cluster assignments.
#' @param group_by The metadata column containing cluster labels.
#' @param p.threshold The p-value or adjusted p-value threshold. Default is 1.
#' @param threshold The expression cutoff above which a gene is considered expressed in a cell. Default is 0.
#' @param P_in.thersh Minimum probability of gene expression within a cluster for inclusion. Default is 0.
#' @param P_out.thersh Maximum probability of gene expression in other clusters for filtering. Default is 1.
#' @param varfeatures A vector of genes to include in the analysis. Defaults to the 2000 most variable genes if NULL.
#' @param omitCluster A vector of cluster names to exclude, useful for omitting small or mixed clusters.
#' @param method P-value adjustment method. Set to "none" for raw p-values. Default is "none".
#' @param Uscore Minimum U-score threshold for `gene_set`. Default is 0.
#' @param n Number of genes per cluster to include in `gene_set`. Default is 3.
#'
#' @return A list with:
#'   - gene_list_long: long-form marker table (data.frame)
#'   - gene_set: top `n` markers per cluster (matrix)
#' @export
FindUniqueMarkers <- function(obj,
                              group_by,
                              p.threshold = 1,
                              threshold = 0,
                              P_in.thersh = 0,
                              P_out.thersh = 1,
                              varfeatures = NULL,
                              omitCluster = NULL,
                              method = "none",
                              Uscore = 0,
                              n = 3) {

  expr_matrix <- if (class(obj[["RNA"]])[1] == "Assay5") {
    obj@assays$RNA@layers$counts
  } else {
    obj@assays$RNA@counts
  }

  groups <- obj@meta.data[[group_by]]

  if (is.null(varfeatures)) {
    varfeatureschoose <- rownames(obj)
  }

  binary_matrix <- expr_matrix[row.names(obj) %in% varfeatureschoose, ] > threshold

  if (any("factor" %in% class(groups))) {
    clusters <- levels(groups)
  } else {
    clusters <- unique(groups)
  }

  percent_stats <- sapply(clusters, function(cluster) {
    Matrix::rowMeans(binary_matrix[, groups == cluster]) * 100
  })

  rownames(percent_stats) <- varfeatureschoose
  colnames(percent_stats) <- clusters

  if (is.null(omitCluster)) {
    ind_column <- 1:length(clusters)
  } else {
    ind_column <- c(1:length(clusters))[-which(colnames(percent_stats) %in% omitCluster)]
  }

  Uscorestats <- t(apply(percent_stats[, ind_column], 1, function(x) {
    sorted_x <- sort(x / 100, decreasing = TRUE)
    x / 100 - ifelse(sorted_x[1] == x / 100, sorted_x[2], sorted_x[1])
  }))

  Pout <- t(apply(percent_stats[, ind_column], 1, function(x) {
    sorted_x <- sort(x / 100, decreasing = TRUE)
    ifelse(sorted_x[1] == x / 100, sorted_x[2], sorted_x[1])
  }))

  Pin <- t(apply(percent_stats, 1, function(x) x / 100))

  Uscoresz <- apply(Uscorestats, 2, function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  })

  raw_pvals <- apply(Uscoresz, 2, function(x) 1 - pnorm(x))

  if (method == "none") {
    pval_matrix <- raw_pvals
    pval_column <- "p.value"
  } else {
    pval_matrix <- apply(raw_pvals, 2, function(x) p.adjust(x, method = method))
    pval_column <- "adj.p.value"
  }

  l_out <- do.call(rbind, lapply(colnames(Uscorestats), function(i) {
    df <- data.frame(
      Gene = rownames(Uscorestats),
      Cluster = i,
      Uscore = Uscorestats[, i],
      P_in = Pin[, i],
      P_out = Pout[, i],
      stringsAsFactors = FALSE
    )
    df[[pval_column]] <- pval_matrix[, i]
    return(df)
  }))

  # Filter for gene_list_long
  l_out <- l_out[l_out[[pval_column]] < p.threshold &
                   l_out$P_in >= P_in.thersh &
                   l_out$P_out <= P_out.thersh, ]

  # Cluster ordering
  l_out$Cluster <- factor(l_out$Cluster, levels = levels(obj@meta.data[[group_by]]))
  l_out <- l_out[order(l_out$Cluster, -l_out$Uscore), ]
  rownames(l_out) <- NULL

  # Create gene_set (wide format top-n genes per cluster)
  gene_set_table <- l_out[
    l_out$Uscore > Uscore &
      l_out$P_in >= P_in.thersh &
      l_out$P_out <= P_out.thersh &
      l_out[[pval_column]] < p.threshold, ]

  gene_split <- split(gene_set_table, gene_set_table$Cluster)
  gene_set <- do.call(cbind, lapply(gene_split, function(df) {
    head(df$Gene, n)
  }))

  return(list(
    gene_list_long = l_out,
    gene_set = gene_set
  ))
}
