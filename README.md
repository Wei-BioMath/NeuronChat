
# NeuronChat

<!-- badges: start -->
<!-- badges: end -->

The goal of NeuronChat is to infer, visualize and analyze neural-specific cell-cell communication from single cell transcriptomics 

## System Requirements

### Hardware requirements

NeuronChat requires a  standard computer and enough memory is recommended to handle large single-cell transcriptomic datasets. 

### Software requirements

#### *** Operating System Requirements

The package has been tested on the following systems: 

``` r
macOS Big Sur, Version 11.5.1  
Windows 10 Pro, version 1909  
```

#### *** R Dependencies (tested and recommended)

``` r
R >= 4.1.0  
dplyr >= 1.0.9
data.table >= 1.14.2  
CellChat >= 1.1.3  
Seurat >= 4.1.0  
SeuratObject >= 4.1.0  
NMF >= 0.23.0  
igraph >= 1.3.4  
ggplot2 >= 3.3.6  
ComplexHeatmap >= 2.8.0  
circlize >= 0.4.14      
ggalluvial >= 0.12.3  
```

## Installation

You can install the development version of NeuronChat like so:

``` r
devtools::install_github("Wei-BioMath/NeuronChat")
```
The time required for installation depends on pre-installed dependencies and the typical time is within 10~20 mins. It's recommended to install the dependencies indicated above before installing NeuronChat. 

## Tutorial 

Please check the vignettes directory of the repo (or click [this link](https://htmlpreview.github.io/?https://github.com/Wei-BioMath/NeuronChat/blob/main/vignettes/NeuronChat-Tutorial.html
)). 

