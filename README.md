# Contents
- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation guide](#Installation-guide)
- [Usage](#Usage)
- [Data Availability](#Data-Availability)
- [License](#license)
- [Issues](#issues)
- [Citation](#citation)
---
## Overview
This repository contains the code used to generate the data and the figures presented in the manuscript <em>“Neural Progenitors as a Novel Pathogenic Mechanism in Microcephaly.”</em>

## Repo Contents
- `R`
  
`cortical_analysis.R` and `rostal_analysis.R` are the main scripts used to analyze the scRNA-seq data. 
`read_human_embryonic_atlas.R` was used to process the data from Braun et al. (2023), *Science* (https://cellxgene.cziscience.com/collections/4d8fed08-2d6d-4692-b5ea-464f1d072077), and to generate the scRNA reference matrix for CIBERSORT bulk deconvolution. This script make use in functions that are in the file `seurat_functions_V2.R`
`label_transfer.R` was used to annotate the cortical scRNA-seq dataset using the Braun et al. dataset as a reference.
## System Requirements

### Hardware requirements
This analysis was performed on a standard workstation.  
For processing large scRNA-seq datasets, at least 16 GB RAM is recommended.

### Software requirements
- R (version 4.4.1)
- Operating system: Linux (Red Hat Enterprise Linux 9.4)

### R package dependencies
The following R packages were used:

- Seurat (5.3.0)
- SeuratObject (5.0.2)
- dplyr (1.1.4)
- ggplot2 (3.5.1)
- patchwork (1.2.0)
- cowplot (1.1.3)
- rio (1.2.3)
- gridExtra (2.3)
- scCustomize (3.0.1)
- SingleR (2.6.0)
- SingleCellExperiment (1.26.0)
- zellkonverter (to convert H5AD to seurat)
## Installation guide
All analyses were performed in R using RStudio.

To run the code, install R (version 4.4.1 or later) and the required R packages listed above. The scripts can then be executed directly within RStudio.
## Usage
The scripts `cortical_analysis.R` and `rostal_analysis.R` accept Cell Ranger count matrices in Matrix Market format (mtx, e.g., from the filtered_feature_bc_matrix directory). An example matrix (W3_E3) is provided in this repository.

To run the scripts, edit lines 47–50 in `cortical_analysis.R` and lines 30–47 in `rostal_analysis.R` to match your local file structure.

The output is an .rds file with seurat object.

To run the script `read_human_embryonic_atlas.R` you need to download the file b40faec8-21d7-4d6c-aa69-4146503d3c64.h5ad from cellxgene (https://cellxgene.cziscience.com/collections/4d8fed08-2d6d-4692-b5ea-464f1d072077)


## Data Availability

All raw and processed data of the scRNAseq as well as the rostal bulk RNAseq used in this study are publicly available from the Gene Expression Omnibus (GEO):GSE229988

Download the data and place it in your local home directory.
## License
This project is licensed under the MIT License – see the LICENSE file for details.
## issues
## citation
