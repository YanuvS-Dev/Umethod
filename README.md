# U-method: Identify Unique Markers in Single-Cell and Spatial Data

**U-method** is an R package for identifying unique markers in
single-cell datasets and visualizing them using high-resolution Visium
HD spatial data. The method is designed for fast, interpretable marker
detection and downstream spatial analysis.

## Example Usage

## Reference Dataset

This tutorial uses a reanalyzed dataset from:

-   **Lee, Hae-Ock, et al. Nature Genetics (2020)**  
    *“Lineage-dependent gene expression programs influence the immune
    landscape of colorectal cancer.”*  
    ArrayExpress:
    [E-MTAB-8410](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8410/)

## Installation

To install **U-method** from GitHub:

    # Install devtools if you haven't already
    install.packages("devtools")

    # Install Umethod package from GitHub
    devtools::install_github("YanuvS-Dev/Umethod")

    # Load the package
    library(Umethod)

## Timing the U-method

We include timing benchmarks to demonstrate just how fast the U-method
runs — from identifying robust markers to generating spatial and
single-cell plots. This pipeline is built to scale. **Timing Note:** All
performance timings reported in this tutorial were measured on a **local
laptop running RStudio**, not on a high-performance server or compute
cluster. This reflects realistic, reproducible desktop usage and
emphasizes the efficiency of the U-method pipeline.

## FindUniqueMarkers Algorithm

The `FindUniqueMarkers` function identifies the most **unique markers**
for each cluster in a Seurat object.

<img src="images/UmethodImage.png" width="1567" />

## 1. Load Example Data

    # Load the U-method package
    devtools::load_all(".")
    library(cowplot)
    library(ggplot2)
    library(svMisc)
    library(Seurat)
    library(reshape2)
    library(arrow)

    # Load the published dataset
    rds_url <- "https://github.com/YanuvS-Dev/Umethod/raw/master/inst/extdata/ColonSinglecellDataLeesUmethodSubsampled10.rds"
    seurat_Full <- readRDS(url(rds_url, "rb"))

## 2. Identify Unique Markers

    UmethodResults <- FindUniqueMarkers(
      obj = seurat_Full,
      group_by = "Celltype",
      method = "BH",
      omitCluster = c("CAFelse", "SmallElse")
    )

## 3. Plot Top U-Markers (UMAP, DotPlot, FeaturePlot)

<img src="README_files/figure-markdown_strict/unnamed-chunk-7-1.png" width="100%" />

    cat("Time to run U-method on single-cell and generate UMAP plots: ", round(difftime(Sys.time(), start_time, units = "secs"), 2), "seconds
    ")

    ## Time to run U-method on single-cell and generate UMAP plots:  62.95 seconds

## Visualizing Markers on Visium HD

### Load Visium HD Data

To visualize marker expression spatially, the U-method integrates Seurat
metadata with high-resolution Visium HD spatial transcriptomic data. In
this tutorial, we use two real colorectal cancer samples:

-   **NAT5 (normal adjacent tissue)** — a healthy reference region.
-   **CRC5 (tumor tissue)** — a malignant region from the same cancer
    type.

We use 8µm-binned versions of the data here due to GitHub storage
constraints. However, the U-method is designed to work directly on
full-resolution (2µm) Visium HD data, achieving &gt;95% unique spot
assignment.

To run the spatial pipeline, two essential files are needed for each
sample:

1.  The filtered gene expression matrix (`filtered_feature_bc_matrix/`)
2.  The spatial barcode metadata file (`tissue_positions.parquet`)

Example data used here was originally published by:

**Oliveira, Michelli F., et al. (2024)**  
*Characterization of immune cell populations in the tumor
microenvironment of colorectal cancer using high-definition spatial
profiling.*  
bioRxiv. DOI: 2024-06

## 4. Load Visium HD Data — Normal Sample (NAT5, 8µm)

    counts_matrix <- Read10X("C:/myGithub/Uemethod_Bigfiles/VisiumHDcolon/NAT5/8um/filtered_feature_bc_matrix/")
    poaraq <- read_parquet("C:/myGithub/Uemethod_Bigfiles/VisiumHDcolon/NAT5/8um/spatial/tissue_positions.parquet")

    seurat_object <- CreateImageData(
      counts_matrix = counts_matrix,
      poaraq = poaraq,
      markers = c(UmethodResults$gene_set)
    )

## 5. Compute Signature Scores for Normal Sample

    datainput_control <- UmethodSignatureMap(
      seurat_object = seurat_object,
      gene_set = UmethodResults$gene_set
    )

    ## Finished average expression calculation of CAF 
    ## Finished average expression calculation of General Fibro 
    ## Finished average expression calculation of Adamdec1 Fibro 
    ## Finished average expression calculation of Normal Muscle 
    ## Finished average expression calculation of Sox6+ Stroma 
    ## Finished average expression calculation of CAP else 
    ## Finished average expression calculation of Endothelial 
    ## Finished average expression calculation of Macrofague 
    ## Finished average expression calculation of T-cells 
    ## Finished average expression calculation of B cells 
    ## Finished average expression calculation of Plasma 
    ## Finished average expression calculation of Epithelial 
    ## Finished average expression calculation of Cancer 
    ## 0.669 Unique spot probability (non-overlapping)

## 6. Plot Spatial Signatures — Normal Sample

    g <- list()
    for (i in datainput_control$Classlist) {
      index <- which(datainput_control$Classlist == i)
      g[[index]] <- ggplot(
        datainput_control$signatureLong[
          datainput_control$signatureLong$Class == i &
            !is.na(datainput_control$signatureLong$value),
        ],
        aes(
          x = pxl_col_in_fullres,
          y = pxl_row_in_fullres,
          color = Class,
          alpha = value / max(value)
        )
      ) +
        geom_point(size = 1) +
        theme_void() +
        scale_y_reverse() +
        scale_color_manual(values = datainput_control$signature_colors[index]) +
        theme(
          plot.background = element_rect(fill = "black"),
          legend.position = "right",
          legend.text = element_text(color = "white", face = "bold", size = 20, angle = 90)
        ) +
        scale_alpha(guide = "none")
    }

    # Reorder plot list: put Epithelial and Cancer first (manually found positions)
    plot_grid(g[[12]], g[[13]], g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], g[[7]], g[[8]], g[[9]], g[[10]],g[[11]], ncol = 3)

<img src="README_files/figure-markdown_strict/unnamed-chunk-11-1.png" width="100%" />

## 7. Load Visium HD Data — Tumor Sample (CRC5, 8µm)

    counts_matrix <- Read10X("C:/myGithub/Uemethod_Bigfiles/VisiumHDcolon/CRC5/8um/filtered_feature_bc_matrix/")
    poaraq <- read_parquet("C:/myGithub/Uemethod_Bigfiles/VisiumHDcolon/CRC5/8um/spatial/tissue_positions.parquet")

    seurat_object <- CreateImageData(
      counts_matrix = counts_matrix,
      poaraq = poaraq,
      markers = c(UmethodResults$gene_set)
    )

## 8. Compute Signature Scores for Tumor Sample

    datainput_control <- UmethodSignatureMap(
      seurat_object = seurat_object,
      gene_set = UmethodResults$gene_set
    )

    ## Finished average expression calculation of CAF 
    ## Finished average expression calculation of General Fibro 
    ## Finished average expression calculation of Adamdec1 Fibro 
    ## Finished average expression calculation of Normal Muscle 
    ## Finished average expression calculation of Sox6+ Stroma 
    ## Finished average expression calculation of CAP else 
    ## Finished average expression calculation of Endothelial 
    ## Finished average expression calculation of Macrofague 
    ## Finished average expression calculation of T-cells 
    ## Finished average expression calculation of B cells 
    ## Finished average expression calculation of Plasma 
    ## Finished average expression calculation of Epithelial 
    ## Finished average expression calculation of Cancer 
    ## 0.747 Unique spot probability (non-overlapping)

## 9. Plot Spatial Signatures — Tumor Sample

    g <- list()
    for (i in datainput_control$Classlist) {
      index <- which(datainput_control$Classlist == i)
      g[[index]] <- ggplot(
        datainput_control$signatureLong[
          datainput_control$signatureLong$Class == i &
            !is.na(datainput_control$signatureLong$value),
        ],
        aes(
          x = pxl_col_in_fullres,
          y = pxl_row_in_fullres,
          color = Class,
          alpha = value / max(value)
        )
      ) +
        geom_point(size = 1) +
        theme_void() +
        scale_y_reverse() +
        scale_color_manual(values = datainput_control$signature_colors[index]) +
        theme(
          plot.background = element_rect(fill = "black"),
          legend.position = "right",
          legend.text = element_text(color = "white", face = "bold", size = 20, angle = 90)
        ) +
        scale_alpha(guide = "none")
    }

    # Reorder plot list: put Epithelial and Cancer first (manually found positions)
    plot_grid(g[[12]], g[[13]], g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], g[[7]], g[[8]], g[[9]], g[[10]],g[[11]], ncol = 3)

<img src="README_files/figure-markdown_strict/unnamed-chunk-14-1.png" width="100%" />

    cat("Total time to run U-method and render both Visium HD panels: ", round(difftime(Sys.time(), start_time, units = "secs"), 2), "seconds
    ")

    ## Total time to run U-method and render both Visium HD panels:  1122.2 seconds

## Notes on Signature Expression

-   **Normal Samples:** CAF signature is not expressed.
-   **Tumor Samples:** Cancer signature expression becomes distinct from
    Epithelial only in tumor tissue.
-   The U-method enables robust class detection with minimal
    preprocessing power and time.

## Citation

Stein Y. *The U-method: Leveraging expression probability for robust
biological marker detection.* Department of Biomolecular Sciences,
Weizmann Institute of Science. (Unpublished yet)
