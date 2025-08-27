---
title: "U-method: Installation and Tutorial"
author: "Yaniv Stein"
date: "2025-06-04"
output: html_document
---

# U-method: Identify Unique Markers in Single-Cell and Spatial Data

**U-method** is an R package for identifying unique markers in single-cell datasets and visualizing them using high-resolution Visium HD spatial data. The method is designed for fast, interpretable marker detection and downstream spatial analysis.

## Example Usage

## Reference Dataset

This tutorial uses a reanalyzed dataset from:

- **Lee, Hae-Ock, et al. Nature Genetics (2020)**  
  *"Lineage-dependent gene expression programs influence the immune landscape of colorectal cancer."*  
  ArrayExpress: [E-MTAB-8410](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8410/)

## Installation

To install **U-method** from GitHub:


``` r
# Install devtools if you haven't already
install.packages("devtools")

# Install Umethod package from GitHub
devtools::install_github("YanuvS-Dev/Umethod")

# Load the package
library(Umethod)
```


## Timing the U-method

We include benchmark results to highlight the speed and efficiency of the U-method — 
from detecting robust markers to generating spatial and single-cell classifications. 

Thanks to its lightweight and scalable implementation, the U-method is ideally suited 
for integration into machine learning pipelines and large-scale data workflows.



## FindUniqueMarkers Algorithm

The `FindUniqueMarkers` function identifies the most **unique markers** for each cluster in a Seurat object.

<div class="figure">
<img src="images/UmethodImage.png" alt="plot of chunk unnamed-chunk-2" width="2089" />
<p class="caption">plot of chunk unnamed-chunk-2</p>
</div>

## 1. Load Example Data


``` r
# Load the U-method package
library(Umethod)
library(Seurat)
library(reshape2)
library(scales)
library(Matrix)

# Optional - plotting results
library(cowplot)
library(ggplot2)
library(arrow)

# Load the published dataset
rds_url <- "https://github.com/YanuvS-Dev/Umethod/raw/master/inst/extdata/ColonSinglecellDataLeesUmethodSubsampled10.rds"
seurat_Full <- readRDS(url(rds_url, "rb"))
```



## 2. Identify Unique Markers


``` r
UmethodResults <- FindUniqueMarkers(
  obj = seurat_Full,
  group_by = "Celltype",
  method = "BH",
  omitCluster = c("CAF else", "Small Else")
)
```

## 3. Plot Top U-Markers (UMAP, DotPlot, FeaturePlot)


``` r
plot_grid(
  plot_grid(
    DimPlot(object = seurat_Full, reduction = "UMAP_on_harmony", pt.size = 0.5, group.by = "Celltype"),
    DotPlot(seurat_Full, features = c(UmethodResults$gene_set), group.by = "Celltype", scale = FALSE) +
      theme(axis.text.x = element_text(angle = 90, vjust = -0.0001)),
    ncol = 1
  ),
  FeaturePlot(
    object = seurat_Full,
    features = UmethodResults$gene_set[1, ],
    cols = c("gray", "blue"),
    reduction = "UMAP_on_harmony",
    ncol = 3,
    order = TRUE
  ),
  ncol = 2
)
```

<div class="figure">
<img src="figure/unnamed-chunk-6-1.png" alt="plot of chunk unnamed-chunk-6" width="100%" />
<p class="caption">plot of chunk unnamed-chunk-6</p>
</div>


``` r
cat("Time to load the data, apply the U-method, and generate UMAP plots: ", round(difftime(Sys.time(), start_time, units = "secs"), 2), "seconds
")
```

```
## Time to load the data, apply the U-method, and generate UMAP plots:  23.15 seconds
```


## Visualizing Markers on Visium HD

### Load Visium HD Data

To visualize marker expression spatially, the U-method integrates Seurat metadata with high-resolution Visium HD spatial transcriptomic data. In this tutorial, we use two real colorectal cancer samples:

- **NAT5 (normal adjacent tissue)** — a healthy reference region.
- **CRC5 (tumor tissue)** — a malignant region from the same patient.

We use 8µm-binned versions of the data here due to GitHub storage constraints. However, the U-method is designed to work directly on full-resolution (2µm) Visium HD data, achieving >95% unique spot assignment.

To run the spatial pipeline, two essential files are needed for each sample:

1. The filtered gene expression matrix (`filtered_feature_bc_matrix/`)
2. The spatial barcode metadata file (`tissue_positions.parquet`)

Example data used here was originally published by:

**Oliveira, Michelli F., et al. (2024)**  
*Characterization of immune cell populations in the tumor microenvironment of colorectal cancer using high-definition spatial profiling.*  
bioRxiv. DOI: 2024-06


## 4. Load Visium HD Data — Normal Sample (NAT5, 8µm)


``` r
counts_matrix <- Read10X("C:/myGithub/Uemethod_Bigfiles/VisiumHDcolon/NAT5/8um/filtered_feature_bc_matrix/")
poaraq <- read_parquet("C:/myGithub/Uemethod_Bigfiles/VisiumHDcolon/NAT5/8um/spatial/tissue_positions.parquet")

seurat_object <- CreateImageData(
  counts_matrix = counts_matrix,
  poaraq = poaraq,
  markers = c(UmethodResults$gene_set)
)
```

## 5. Compute Signature Scores for Normal Sample


``` r
datainput_control <- UmethodSignatureMap(
  seurat_object = seurat_object,
  gene_set = UmethodResults$gene_set
)
```

```
## Finished average expression calculation of CAF 
## Finished average expression calculation of General Fibro 
## Finished average expression calculation of Adamdec1 Fibro 
## Finished average expression calculation of Normal Muscle 
## Finished average expression calculation of Sox6+ Fibro 
## Finished average expression calculation of CAP else 
## Finished average expression calculation of Endothelial 
## Finished average expression calculation of Macrophage 
## Finished average expression calculation of T-cells 
## Finished average expression calculation of B-cells 
## Finished average expression calculation of Plasma 
## Finished average expression calculation of Epithelial 
## Finished average expression calculation of Cancer 
## 0.669 Unique spot probability (non-overlapping)
```

## 6. Plot Spatial Signatures — Normal Sample


``` r
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
```

<div class="figure">
<img src="figure/unnamed-chunk-10-1.png" alt="plot of chunk unnamed-chunk-10" width="100%" />
<p class="caption">plot of chunk unnamed-chunk-10</p>
</div>

## 7. Load Visium HD Data — Tumor Sample (CRC5, 8µm)


``` r
counts_matrix <- Read10X("C:/myGithub/Uemethod_Bigfiles/VisiumHDcolon/CRC5/8um/filtered_feature_bc_matrix/")
poaraq <- read_parquet("C:/myGithub/Uemethod_Bigfiles/VisiumHDcolon/CRC5/8um/spatial/tissue_positions.parquet")

seurat_object <- CreateImageData(
  counts_matrix = counts_matrix,
  poaraq = poaraq,
  markers = c(UmethodResults$gene_set)
)
```

## 8. Compute Signature Scores for Tumor Sample


``` r
datainput_crc <- UmethodSignatureMap(
  seurat_object = seurat_object,
  gene_set = UmethodResults$gene_set
)
```

```
## Finished average expression calculation of CAF 
## Finished average expression calculation of General Fibro 
## Finished average expression calculation of Adamdec1 Fibro 
## Finished average expression calculation of Normal Muscle 
## Finished average expression calculation of Sox6+ Fibro 
## Finished average expression calculation of CAP else 
## Finished average expression calculation of Endothelial 
## Finished average expression calculation of Macrophage 
## Finished average expression calculation of T-cells 
## Finished average expression calculation of B-cells 
## Finished average expression calculation of Plasma 
## Finished average expression calculation of Epithelial 
## Finished average expression calculation of Cancer 
## 0.747 Unique spot probability (non-overlapping)
```

## 9. Plot Spatial Signatures — Tumor Sample


``` r
g <- list()
for (i in datainput_crc$Classlist) {
  index <- which(datainput_crc$Classlist == i)
  g[[index]] <- ggplot(
    datainput_crc$signatureLong[
      datainput_crc$signatureLong$Class == i &
        !is.na(datainput_crc$signatureLong$value),
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
    scale_color_manual(values = datainput_crc$signature_colors[index]) +
    theme(
      plot.background = element_rect(fill = "black"),
      legend.position = "right",
      legend.text = element_text(color = "white", face = "bold", size = 20, angle = 90)
    ) +
    scale_alpha(guide = "none")
}

# Reorder plot list: put Epithelial and Cancer first (manually found positions)
plot_grid(g[[12]], g[[13]], g[[1]], g[[2]], g[[3]], g[[4]], g[[5]], g[[7]], g[[8]], g[[9]], g[[10]],g[[11]], ncol = 3)
```

<div class="figure">
<img src="figure/unnamed-chunk-13-1.png" alt="plot of chunk unnamed-chunk-13" width="100%" />
<p class="caption">plot of chunk unnamed-chunk-13</p>
</div>


``` r
cat("Total time to run U-method and render both Visium HD panels: ", round(difftime(Sys.time(), start_time, units = "secs"), 2), "seconds
")
```

```
## Total time to run U-method and render both Visium HD panels:  477.22 seconds
```

## Notes on Signature Expression

- **Normal Samples:** CAF signatures are typically absent.
- **Tumor Samples:** Cancer signatures become distinct from Epithelial only in tumor tissue.
- The U-method enables robust class detection with minimal preprocessing power and time.
- The U-method enables robust cluster-specific marker detection with minimal preprocessing, low computational cost, and fast runtime.


## Guidelines for Use

- Apply FindUniqueMarkers() on a clustered Seurat object to detect uniquely expressed genes (U-markers).

- Start by running the U-method with default parameters (no filters). This provides a first look at all clusters.
[Clusters with no obvious U-markers may represent mixed cell populations or potential misclustering. The U-method can serve as a simple goodness-of-fit check for the clustering]

- For standard marker detection, rerun with a U-score threshold of 0.2 (as used in the manuscript).

- If needed, increase the threshold or apply additional filters to tighten marker selection, and omit ambiguous clusters from the background before reintegration.

- visualization of the results: Dot plots are generated with raw (unscaled) expression.

- For spatial data:

  - Build a Seurat object with CreateImageData().

  - Project signatures and assign clusters with UmethodSignatureMap().



## Citation

Stein Y. *The U-method: Leveraging expression probability for robust biological marker detection.* Department of Biomolecular Sciences, Weizmann Institute of Science. (Unpublished yet)
