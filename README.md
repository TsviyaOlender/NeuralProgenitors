# Contents
- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [License](#license)
- [Issues](#issues)
- [Citation](#citation)
---
## Overview
This repository contains the code used to generate the data presented in the manuscript <em>“Neural Progenitors as a Novel Pathogenic Mechanism in Microcephaly.”</em>
The following scripts are included:

## Repo Contents

The following scripts:

- `cortical_analysis.R`
- `rostal_analysis.R`
- `label_transfer.R`
were used to analyze the scRNA-seq data using Seurat.

The script `read_human_embryonic_atlas.R` was used to process the data from Braun et al. (2023), *Science* (downloaded from https://cellxgene.cziscience.com/collections/4d8fed08-2d6d-4692-b5ea-464f1d072077), and to generate the scRNA reference matrix for CIBERSORT bulk deconvolution.
The script <code>label_transfer.R</code> was used to annotate the corticalscRNA-seq dataset using the Braun et al. dataset as a reference.

## system-requirements
## License
This project is licensed under the MIT License – see the LICENSE file for details.
## issues
## citation
