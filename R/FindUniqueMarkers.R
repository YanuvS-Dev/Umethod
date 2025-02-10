
#' Find Unique Markers - U method
#'
#' This function identifies the most unique markers per cluster in a Seurat object.
#' If the cluster are labeled as charachters, and not numbers, matchnames parameter shoud be F.
#'
#' @param obj Seurat object after any clustering
#' @param group_by The cluster labels in metafile of obj
#' @param p.threshold p.value threshold for the Umethod
#' @param threshold what is the expression threshold that a cell is counted as positive for this gene, default is 0.
#' @param P_in.thersh probability of expressing inside the cluster threshold
#' @param varfeatures Which features to enter the analysis. For correct analysis = row.names(seurat.obj),if NULL this method will run on the 2000 most variable genes
#' @param smallcluster a vector of cluster names so they will be omited from the analysis, if there are small mixed clusters it is best if you ommit them here
#' @return A data.frame object with the most uniquely expressed genes per cluster,that passed filtering parameters.
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
      progress(i, max.value = length(jumpind) - 2)

      features <- row.names(obj)[(jumpind[i] + 1):jumpind[i + 1]]
      temp_percent <- Percent_Expressing(seurat_object = obj, features = features, threshold = threshold, group_by = group_by)

      percent_stats[rownames(temp_percent), colnames(temp_percent)] <- temp_percent
    }

    # Finish progress bar
    progress(length(jumpind) - 2, max.value = length(jumpind) - 2)
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
