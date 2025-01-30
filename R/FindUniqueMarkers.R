
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
#' @return A data.frame object with the most uniquely expressed genes per cluster,that passed filtering parameters.
#' @export

FindUniqueMarkers <- function(obj,
                              group_by,
                              p.threshold = 0.2,
                              threshold = 0,
                              P_in.thersh = 0,
                              P_out.thersh = 1,
                              varfeatures = NULL,
                              reduction = "UMAP_on_harmony",
                              onlyOnetime =T,
                              smallcluster = NULL,
                              method = "BH",
                              matchnames = T, ## Correct names - can cause errors if group_by is not factor
                              matchnames_try2=F, #If the first doesnt fit
                              jumpFix = 200) # to fix the assay problem
{
  if(is.null(varfeatures))
  {
    nvarfeatures <- length(obj@assays[[DefaultAssay(object = obj)]]@var.features)
    percent_stats <- Percent_Expressing(seurat_object = obj, features = obj@assays[[DefaultAssay(object = obj)]]@var.features[1:jumpFix], threshold = threshold,group_by = group_by)
    jumpind <- seq(0,nvarfeatures,jumpFix)[-1]
    for(i in 1:(length(seq(0,nvarfeatures,jumpFix)[-1])-1))
    {
      percent_stats <- rbind(percent_stats,Percent_Expressing(seurat_object = obj, features = obj@assays[[DefaultAssay(object = obj)]]@var.features[(1+jumpind[i]):jumpind[i+1]], threshold = 0,group_by = group_by))
    }
  }else
  {
    nvarfeatures <- length(varfeatures)
    percent_stats <- Percent_Expressing(seurat_object = obj, features = varfeatures[1:jumpFix], threshold = threshold,group_by = group_by)
    jumpind <- seq(0,nvarfeatures,jumpFix)[-1]
    for(i in 1:(length(seq(0,nvarfeatures,jumpFix)[-1])-1))
    {
      percent_stats <- rbind(percent_stats,Percent_Expressing(seurat_object = obj, features = varfeatures[(1+jumpind[i]):jumpind[i+1]], threshold = threshold,group_by = group_by))
    }
  }


  ##correct names
  if(matchnames){
    if(matchnames_try2 & all(substr(colnames(percent_stats),1,1) == "X")){
      colnames(percent_stats) <- str_replace(colnames(percent_stats),pattern = "X","")

    }else
    {
      percent_stats <- percent_stats[,amatch(levels(obj@meta.data[,group_by]),colnames(percent_stats) , maxDist = Inf)]
      colnames(percent_stats) <- levels(obj@meta.data[,group_by])

    }
  }

  # Calculating u-scores, deviance and probability ratios
  MinDevstats <- t(apply(percent_stats[,which(!(colnames(percent_stats) %in% smallcluster))],1,function(x){ifelse(sort(x/100, decreasing = TRUE)[1] - x/100 == 0,x/100 - sort(x/100, decreasing = TRUE)[2],x/100 - sort(x/100, decreasing = TRUE)[1])})) # (MI-max(MItag))
  Ustatsnew <- t(apply(percent_stats[,which(!(colnames(percent_stats) %in% smallcluster))],1,function(x){x/100*(ifelse(sort(x/100, decreasing = TRUE)[1] - x/100 == 0,x/100 - sort(x/100, decreasing = TRUE)[2],x/100 - sort(x/100, decreasing = TRUE)[1]))})) # (MI-max(MItag))
  Pout <- t(apply(percent_stats[,which(!(colnames(percent_stats) %in% smallcluster))],1,function(x){ifelse(sort(x/100, decreasing = TRUE)[1] - x/100 == 0,sort(x/100, decreasing = TRUE)[2],sort(x/100, decreasing = TRUE)[1])})) # (MI-max(MItag))
  Pin<- t(apply(percent_stats[,which(!(colnames(percent_stats) %in% smallcluster))],1,function(x){x/100})) # (MI-max(MItag))

  Uznew <- apply(Ustatsnew,2,function(x){(x-mean(x))/sd(x)})

  order_ind <- sapply(as.list(colnames(Uznew)),function(x){names(sort(Uznew[,x],decreasing = T))})

  # adjusted p-value calculation
  UznewPvalue <- apply(Uznew,2,function(x){p.adjust(1 - pnorm(x),method = method)}) # NEW
  UznewPvalue <- round(UznewPvalue[order_ind[,2],],10) #  # when I did the analysis on gene, now its on cluster

  ###
  l <- list()

  for(i in colnames(Ustatsnew))
  {
    l[[i]] <- data.frame(Gene = names(sort(Ustatsnew[,i],decreasing = T)),
                         Cluster = i,Uscore =sort(Ustatsnew[,i],decreasing = T),
                         Dev = MinDevstats[names(sort(Ustatsnew[,i],decreasing = T)),i],
                         adj.p.value = UznewPvalue[names(sort(Ustatsnew[,i],decreasing = T)),i],
                         P_in = Pin[names(sort(Ustatsnew[,i],decreasing = T)),i],
                         P_out = Pout[names(sort(Ustatsnew[,i],decreasing = T)),i])
  }

  l_out <- l[[1]]
  for(i in colnames(Ustatsnew)[-1])
  {
    l_out <- rbind(l_out,l[[i]])
  }
  if(onlyOnetime)
  {
    onlyonce <- names(which(table(l_out$Gene,l_out$adj.p.value < p.threshold)[,2] == 1))
    l_out[which(l_out$Gene %in% onlyonce),]
  }
  return(l_out[which( l_out$adj.p.value < p.threshold & l_out$P_in >= P_in.thersh & l_out$P_out <= P_out.thersh ),])

}
