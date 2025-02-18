
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
#' @param jumpFix An integer specifying how many genes to process in each iteration to avoid potential bugs with large gene lists.
#' @param method The p-value adjustment method. Defaults to "BH" (Benjamini-Hochberg). Set to "none" for raw p-values or choose other methods available in the `p.adjust` function.
#' @param progresstext Boolean indicating whether to display a progress bar during the analysis. Defaults to `TRUE`.
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
                              method = "BH",
                              progresstext = T,
                              jumpFix = 200) # to fix the assay problem
{
  if(is.null(varfeatures))
  {
    # Define variables
    nvarfeatures <- length(row.names(obj))
    jumpind <- seq(0, nvarfeatures, jumpFix)

    # Create an empty data frame with genes as rows and clusters as columns
    clusters <- unique(obj@meta.data[[group_by]])
    percent_stats <- data.frame(matrix(0, nrow = nvarfeatures, ncol = length(clusters)))
    rownames(percent_stats) <- row.names(obj)
    colnames(percent_stats) <- clusters

    # Calculate initial percentages
    initial_percent <- Percent_Expressing(seurat_object = obj, features = row.names(obj)[1:jumpFix], threshold = threshold, group_by = group_by)
    percent_stats[rownames(initial_percent), colnames(initial_percent)] <- initial_percent

    # Loop with progress bar
    for (i in 1:(length(jumpind) - 2)) {
      if(progresstext){progress(i, max.value = length(jumpind) - 2)}

      features <- row.names(obj)[(jumpind[i] + 1):jumpind[i + 1]]
      temp_percent <- Percent_Expressing(seurat_object = obj, features = features, threshold = threshold, group_by = group_by)

      percent_stats[rownames(temp_percent), colnames(temp_percent)] <- temp_percent
    }

    # Finish progress bar
    if(progresstext){progress(length(jumpind) - 2, max.value = length(jumpind) - 2)}
  }else
  {
    # Define variables
    nvarfeatures <- length(varfeatures)
    jumpind <- seq(0, nvarfeatures, jumpFix)

    # Create an empty data frame with genes as rows and clusters as columns
    clusters <- unique(obj@meta.data[[group_by]])
    percent_stats <- data.frame(matrix(0, nrow = nvarfeatures, ncol = length(clusters)))
    rownames(percent_stats) <- varfeatures
    colnames(percent_stats) <- clusters

    # Calculate initial percentages
    initial_percent <- Percent_Expressing(seurat_object = obj, features = varfeatures[1:jumpFix], threshold = threshold, group_by = group_by)
    percent_stats[rownames(initial_percent), colnames(initial_percent)] <- initial_percent

    # Loop with progress bar
    for (i in 1:(length(jumpind) - 2)) {
      progress(i, max.value = length(jumpind) - 2)

      features <- varfeatures[(jumpind[i] + 1):jumpind[i + 1]]
      temp_percent <- Percent_Expressing(seurat_object = obj, features = features, threshold = threshold, group_by = group_by)

      percent_stats[rownames(temp_percent), colnames(temp_percent)] <- temp_percent
    }

    # Finish progress bar
    progress(length(jumpind) - 2, max.value = length(jumpind) - 2)
  }


  # Calculating u-scores, Pin and Pout
  Uscorestats <- t(apply(percent_stats[,which(!(colnames(percent_stats) %in% smallcluster))],1,function(x){ifelse(sort(x/100, decreasing = TRUE)[1] - x/100 == 0,x/100 - sort(x/100, decreasing = TRUE)[2],x/100 - sort(x/100, decreasing = TRUE)[1])})) # (MI-max(MItag))
  Pout <- t(apply(percent_stats[,which(!(colnames(percent_stats) %in% smallcluster))],1,function(x){ifelse(sort(x/100, decreasing = TRUE)[1] - x/100 == 0,sort(x/100, decreasing = TRUE)[2],sort(x/100, decreasing = TRUE)[1])})) # (MI-max(MItag))
  Pin<- t(apply(percent_stats[,which(!(colnames(percent_stats) %in% smallcluster))],1,function(x){x/100})) # (MI-max(MItag))

  Uscoresz <- apply(Uscorestats,2,function(x){(x-mean(x,na.rm = T))/sd(x,na.rm = T)})

  order_ind <- sapply(as.list(colnames(Uscoresz)),function(x){names(sort(Uscoresz[,x],decreasing = T))})

  # adjusted p-value calculation
  UscoreszPvalue <- apply(Uscoresz,2,function(x){p.adjust(1 - pnorm(x),method = method)}) # NEW

  ###
  l <- list()

  for(i in colnames(Uscorestats))
  {
    l[[i]] <- data.frame(Gene = names(sort(Uscorestats[,i],decreasing = T)),
                         Cluster = i,Uscore =sort(Uscorestats[,i],decreasing = T),
                         adj.p.value = UscoreszPvalue[names(sort(Uscorestats[,i],decreasing = T)),i],
                         P_in = Pin[names(sort(Uscorestats[,i],decreasing = T)),i],
                         P_out = Pout[names(sort(Uscorestats[,i],decreasing = T)),i])
  }

  l_out <- l[[1]]
  for(i in colnames(Uscorestats)[-1])
  {
    l_out <- rbind(l_out,l[[i]])
  }

  return(l_out[which( l_out$adj.p.value < p.threshold & l_out$P_in >= P_in.thersh & l_out$P_out <= P_out.thersh ),])

}
