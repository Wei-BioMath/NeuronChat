
# NeuronChat

<!-- badges: start -->
<!-- badges: end -->

The goal of NeuronChat is to infer, visualize and analyze neural-specific cell-cell communication from single cell transcriptomics or spatially resolved transcriptomics data. 

## Installation

You can install the development version of NeuronChat like so:

``` r
devtools::install_github("Wei-BioMath/NeuronChat")
```
The time required for installation depends on pre-installed dependencies and the typical time is within 10 mins. Dependencies of NeuronChat package can be automatically installed when installing NeuronChat pacakge (When encountering any issues, please see the `Troubleshooting` section below).

## Usage & Tutorial 
### Basic usage of NeuronChat
``` r
library(NeuronChat)
# creat NeuronChat object 
x <- createNeuronChat(normalized_count_mtx,DB='mouse',group.by = cell_group_vector) # use DB='human' for human data
# calculation of communication networks  
x <- run_NeuronChat(x,M=100)
# aggregating networks over interaction pairs
net_aggregated_x <- net_aggregation(x@net,method = 'weight')
# visualization
netVisual_circle_neuron(net_aggregated_x)
```
### Full tutorials
A set of tutorials of NeuronChat are as follows:
- 
- [Update the NeuronChat interaction database]()

Please check the vignettes directory of the repo (or click [this link](https://htmlpreview.github.io/?https://github.com/Wei-BioMath/NeuronChat/blob/main/vignettes/NeuronChat-Tutorial.html
)). 

## System Requirements

### Hardware requirements

NeuronChat requires a  standard computer and enough memory is recommended to handle large single-cell transcriptomic datasets. 

### Software requirements

####  ** Operating System Requirements

The package has been tested on the following systems: 

``` r
macOS Big Sur, Version 11.5.1  
Windows 10 Pro, version 1909  
```

#### ** R Dependencies (tested and recommended)

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

##  Troubleshooting 
Here are some instructions in case users might encounter issues during installation: 

- Install [CellChat](https://github.com/sqjin/CellChat) using `devtools::install_github("sqjin/CellChat")`. See also the `Installation of other dependencies
` in the [README](https://github.com/sqjin/CellChat) of CellChat for any other issues (including `NMF` package installation). 
- Install [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap) using `devtools::install_github("jokergoo/ComplexHeatmap")`. 
- Install [circlize](https://github.com/jokergoo/circlize) using `devtools::install_github("jokergoo/circlize")`. 
- Please see the instructions for installation of `Seurat` or `SeuratObject` via the link [Seurat](https://satijalab.org/seurat/articles/install.html) or [SeuratObject](https://github.com/mojaveazure/seurat-object). 
 
