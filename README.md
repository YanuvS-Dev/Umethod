# Umethod

**Umethod** is an R package designed for identifying unique markers in
single-cell data sets, Use the FindUniqueMarkers function on a Seurat
object after clustering to get the most unique markers for your
clusters. CreateImageData function can be used with U markers (or any
other markers) for each cluster for downstream analysis, visualizing the
markers on Visium HD spatial data.

## 🚀 Installation

To install Umethod from GitHub:

#### Install devtools if you haven’t already

install.packages(“devtools”)

#### Install Umethod from GitHub

devtools::install\_github(“YanuvS/Umethod”)

#### Load the package

library(Umethod)

## FindUniqueMarkers algorithm

<figure>
<img src="images/UmethodImage.png" style="width:30.0%"
alt="process: Scoring each gene for each cluster, then testing for significant markers" />
<figcaption aria-hidden="true">process: Scoring each gene for each
cluster, then testing for significant markers</figcaption>
</figure>

### 📈 Example Usage

**Reference:** \[Lee, Hae-Ock, et al. Nature genetics (2020)\].
*“Lineage-dependent gene expression programs influence the immune
landscape of colorectal cancer.”*.

    # Load Umethod
    library(Umethod)
    library("scCustomize")
    library(cowplot)
    library(ggplot2)
    library(svMisc)

    # Load the published data set (replace with the actual data loading code)
    seurat_Full <- readRDS("C:\\Migration\\R projects\\Umethod\\ColonSinglecellDataLeesUmethod.rds")

    # Apply Umethod functions, if there are small/mixed clusters, their name should be added to smallcluster variable to omit them.
    # The progress bar prints weird massages, in this rmd file I suppress it
    genes_list <- FindUniqueMarkers(
        obj = seurat_Full,
        group_by = "Celltype",
        method = "none",
        smallcluster = c("CAFelse", "SmallElse"),
        progresstext = F)

    # gene_list is the marker list ordered by score and cluster
    head(genes_list)

    ##            Gene Cluster    Uscore  adj.p.value      P_in     P_out
    ## S100P     S100P  Cancer 0.5768754 4.122014e-07 0.8216011 0.2447257
    ## LCN2       LCN2  Cancer 0.4925141 1.169598e-05 0.7696112 0.2770971
    ## ASCL2     ASCL2  Cancer 0.4617432 3.525535e-05 0.5739015 0.1121583
    ## CEACAM6 CEACAM6  Cancer 0.4576173 4.068412e-05 0.7196342 0.2620170
    ## ANXA3     ANXA3  Cancer 0.4400294 7.398427e-05 0.6360709 0.1960415
    ## SCD         SCD  Cancer 0.4278370 1.106696e-04 0.5943754 0.1665384

    # Choose thresholds 
    Uscore <- 0.25
    p_in <- 0.4
    p_out <- 0.3

    # Pulling the top 5 markers and the name of the top U marker for each cluster

    genesetshort <- unlist(sapply(split(genes_list[genes_list$Uscore > Uscore & genes_list$P_in > p_in & genes_list$P_out < p_out,],genes_list[genes_list$Uscore > Uscore & genes_list$P_in > p_in& genes_list$P_out < p_out,]$Cluster),function(x){x[[1]][1]}))

    genesetlong <- unique(unlist(sapply(split(genes_list[genes_list$Uscore > Uscore & genes_list$P_in > p_in& genes_list$P_out < p_out,],genes_list[genes_list$Uscore > Uscore & genes_list$P_in > p_in& genes_list$P_out < p_out,]$Cluster),function(x){x[[1]][1:5]})))

    genesetlong

    ##      Adamdec1.Fibro B.cells     CAF       Cancer    CAP.else Endothelial Epithelial General.Fibro Macrofague Normal.Muscle
    ## [1,] "ADAMDEC1"     "MS4A1"     "COL11A1" "S100P"   "KCNJ8"  "ECSCR.1"   "CA2"      "OGN"         "TYROBP"   "RERGL"      
    ## [2,] "HAPLN1"       "BANK1"     "FAP"     "LCN2"    "HIGD1B" "CLDN5"     "VSIG2"    "PCOLCE2"     "AIF1"     "PLN"        
    ## [3,] "CCL13"        "TNFRSF13C" "PODNL1"  "ASCL2"   "ENPEP"  "PLVAP"     "GUCA2A"   "PI16"        "FCER1G"   "NTRK2"      
    ## [4,] "SFTA1P"       "VPREB3"    "COL10A1" "CEACAM6" "GJC1"   "VWF"       "ADH1C"    "CILP"        "LST1"     "C2orf40"    
    ## [5,] "CCL8"         NA          "TMEM158" "ANXA3"   "EDNRA"  "PCAT19"    "MT1H"     "C1QTNF3"     "FCGR2A"   "ACTG2"      
    ##      Plasma     Sox6..Stroma T.cells
    ## [1,] "MZB1"     "NSG1"       "CD3D" 
    ## [2,] "DERL3"    "ENHO"       "CD3E" 
    ## [3,] "TNFRSF17" "BMP5"       "CD7"  
    ## [4,] "CD27"     "SOX6"       "CD2"  
    ## [5,] "FAM46C"   "EDNRB"      "TRBC1"

    #Order the cluster that had at least one marker that passed threshold, as you want them to be plotted.
    clusterorder<- c("CAF","General.Fibro","Adamdec1.Fibro","Normal.Muscle","Sox6..Stroma","CAP.else","Endothelial","Macrofague","T.cells","B.cells","Plasma","Epithelial","Cancer")
    # Ordering the clusters that had any umarkers from genesetlong for dotplot
    indclusters <- rep(NA,dim(genesetlong)[2])
    for(i in 1:dim(genesetlong)[2]){indclusters[i] <- which(colnames(genesetlong) == clusterorder[i])}

    genesetlong <- genesetlong[,indclusters]
    genesetshort <- genesetshort[indclusters]

# Plotting the results of the top U markers for each cluster

<img src="README_files/figure-markdown_strict/unnamed-chunk-10-1.png" width="100%" />
