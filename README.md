# Contents
- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [License](#license)
- [Issues](#issues)
- [Citation](#citation)
---
## Overview
This repository contains the code used to generate the data and the figures presented in the manuscript <em>“Neural Progenitors as a Novel Pathogenic Mechanism in Microcephaly.”</em>

## Repo Contents
- `R`
  
`cortical_analysis.R` and `cortical_analysis.R` are the main scripts used to analyze the scRNA-seq data. 
`read_human_embryonic_atlas.R` was used to process the data from Braun et al. (2023), *Science* (https://cellxgene.cziscience.com/collections/4d8fed08-2d6d-4692-b5ea-464f1d072077), and to generate the scRNA reference matrix for CIBERSORT bulk deconvolution.
`label_transfer.R` was used to annotate the cortical scRNA-seq dataset using the Braun et al. dataset as a reference.



## system-requirements
## License
This project is licensed under the MIT License – see the LICENSE file for details.
## issues
## citation
