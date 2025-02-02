# Umethod

**Umethod** is an R package designed for identifying unique markers in
single-cell data sets, Use the FindUniqueMarkers function on a Seurat
object after clustering to get the most unique markers for your
clusters. CreateImageData can be used with u markers (or any else) for
each cluster for downstream analysis, visualizing the markers on Visium
HD spatial data.

## 🚀 Installation

To install Umethod from GitHub:

#### Install devtools if you haven’t already

install.packages(“devtools”)

#### Install Umethod from GitHub

devtools::install\_github(“YanuvS/Umethod”)

#### Load the package

library(Umethod)

## FindUniqueMarkers algorithm

![U method description; scoring each gene for each cluster, then testing
for significant
markers](C:/Migration/R%20projects/Umethod/UmethodImage.png) \### 📈
Example Usage

This example uses published single-cell RNA-seq data from **\[Dataset
Name\]**  
**Reference:** \[Lee, Hae-Ock, et al. Nature genetics 52.6 (2020):
594-603.\]. *“Lineage-dependent gene expression programs influence the
immune landscape of colorectal cancer.”*.

    # Load Umethod
    library(Umethod)
    library("scCustomize")
    library(cowplot)
    library(ggplot2)

    # Load the published data set (replace with the actual data loading code)
    seurat_Full <- readRDS("C:\\Migration\\R projects\\Umethod\\ColonSinglecellDataLeesUmethod.rds")

    # Apply Umethod functions, if there are small/mixed clusters, their name should be added to smallcluster variable to omit them.
    # If the clusters are characters and not numbers, add matchnames = F
    genes_list <- FindUniqueMarkers(obj = seurat_Full,group_by = "Celltype",p.threshold = 0.2,varfeatures = row.names(seurat_Full),method = "none",smallcluster = c("CAFelse","SmallElse"),matchnames = F)

    # Choose thresholds 
    u.threshold <- 0.1
    MinDEVthresh <- 0.2
    P_out <- 0.3

    # Pulling the top 5 markers and the name of the top U marker for each cluster
    genesetshort <- unlist(sapply(split(genes_list[genes_list$adj.p.value < u.threshold & genes_list$Dev > MinDEVthresh & genes_list$P_out < P_out,],genes_list[genes_list$adj.p.value < u.threshold & genes_list$Dev > MinDEVthresh& genes_list$P_out < P_out,]$Cluster),function(x){x[[1]][1]}))
    genesetlong <- unique(unlist(sapply(split(genes_list[genes_list$adj.p.value < u.threshold & genes_list$Dev > MinDEVthresh& genes_list$P_out < P_out,],genes_list[genes_list$adj.p.value < u.threshold & genes_list$Dev > MinDEVthresh& genes_list$P_out < P_out,]$Cluster),function(x){x[[1]][1:5]})))

    # Ordering the clusters that had any umarkers from genesetlong for dotplot
    indclusters <- rep(NA,dim(genesetlong)[2])
    for(i in 1:dim(genesetlong)[2]){indclusters[i] <- which(colnames(genesetlong) == c("CAF","General.Fibro","Adamdec1.Fibro","Normal.Muscle","Sox6..Stroma","CAP.else","Endothelial","Macrofague","T.cells","B.cells","Plasma","Epithelial","Cancer")[i])}

    genesetlong <- genesetlong[,indclusters]
    genesetshort <- genesetshort[indclusters]

# Plotting the results of the top U markers for each cluster

    plot_grid(plot_grid(DimPlot(object = seurat_Full, reduction = "UMAP_on_harmony",pt.size = 0.5, group.by = "Celltype"),
                        DotPlot(seurat_Full,features = c(genesetlong),group.by = "Celltype",scale = F) + theme(axis.text.x = element_text(angle = 90,vjust = -0.0001)),ncol = 1),
              FeaturePlot(object = seurat_Full, features = genesetshort, cols = c("gray", "blue"),reduction = "UMAP_on_harmony",ncol = 3,order = T),ncol = 2)

![](README_files/figure-markdown_strict/unnamed-chunk-2-1.png)
