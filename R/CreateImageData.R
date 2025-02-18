
#' U-Method CreateImageData function
#'
#' This function pulles the read count from the visium dataset for each spot for selected markers
#'
#' @param bc_matrix 10x filtered matrix, can be loaded by Read10X(".../Visium_HD__binned_outputsXum/filtered_feature_bc_matrix/"), X in (2,8,16).
#' @param poaraq poaraq file can be loaded by read_parquet(".../Visium_HD__binned_outputX2um/spatial/tissue_positions.parquet")
#' @param markers marker list to plot on visium HD data
#' @param print.table should the number of spots per marker be plotted
#' @param project name of the project passed to output Seurat object
#' @return Seurat object with all data on markers needed for downstream analysis
#' @export

CreateImageData <- function(bc_matrix,poaraq,markers,project = NULL,print.table = F)
{

  # Creating the seurat object from the Visium output
  seurat_object <- CreateSeuratObject(
    counts = bc_matrix,
    assay = "Spatial",
    project = ifelse(is.null(project),"Visium & U-Method Project",project)
  )

  # Crating a dataframe of markers expression by cells, to merge with seurat_object metadata
  df <- data.frame(t(bc_matrix[which(row.names(bc_matrix) %in% markers),]))
  df[,dim(df)[2]+1] <- row.names(df)
  names(df)[dim(df)[2]] <- "cell_id"

  # Add cell metadata to the Seurat object
  # Ensure cell IDs in metadata match with column names of the expression matrix
  cell_metadata <- poaraq[poaraq$barcode %in% colnames(seurat_object),] %>%
    dplyr::rename(cell_id = barcode)  # Ensure column names match

  cell_metadata <- merge(cell_metadata,df,by = "cell_id")

  # Add metadata to the Seurat object
  seurat_object <- AddMetaData(seurat_object, metadata = as.data.frame(cell_metadata))


  # Classifing to one of the markers
  seurat_object$Class <- sapply((apply(seurat_object@meta.data[,markers],1,function(x){paste(markers[x != 0],sep = "&")})),function(y){ifelse(length(y) < 1,NA,y)})

  if(print.table){
    print(table(seurat_object$Class))
  }

  seurat_object$Class <- as.factor(seurat_object$Class)

  seurat_object$ind_overlap <- apply(cell_metadata[,markers],1,function(x){length(which(x > 0))})

  return(seurat_object)
}
