# ARC-ICN-Aphasia

Analysis pipeline for: **ICN Engagement is associated with Aphasia Severity in Chronic Stroke**

## Overview
MATLAB scripts for analysing intrinsic connectivity network (ICN) engagement during picture naming in the Aphasia Recovery Cohort (ARC).

## Requirements
- MATLAB R2025a (or compatible)
- SPM12 (Wellcome Centre for Human Neuroimaging)
- ICN_Atlas toolbox (version r20180306; Kozák et al., 2017)

## Data
The analysis uses the ARC dataset (OpenNeuro accession ds004884, version 1.0.2). Raw BIDS images, first-level SPM outputs, and demographic tables are not redistributed here. Some scripts may expect a file named `participants.txt`, depending on how the raw ARC metadata is organised locally. If it is missing, create it from the original ARC `participants.tsv` by selecting the relevant columns and saving as tab-delimited text. No demographic data are redistributed in this repository.

The Aphasia Recovery Cohort is publicly available on OpenNeuro:  
https://openneuro.org/datasets/ds004884/versions/1.0.2

doi:10.18112/openneuro.ds004884.v1.0.2

### Atlases and tools
- BrainMap20
- SMITH10
- JHU white matter atlas
- Neurosynth v7 meta-analytic maps
- ICN_Atlas toolbox, Z-threshold = 3

## Scripts

| Script | Description |
|--------|-------------|
| ARC_01_DataInventory.mlx | Data inventory and initial setup |
| ARC_02b_xATL_Diagnostics.m | Atlas diagnostics |
| ARC_02c_ConFile_Inventory.m | Contrast file inventory |
| ARC_03a_Individual_ICN_Extraction.m | ICN metric extraction per participant |
| ARC_03b_v3_Compile_Results.m | Compile results and apply IRi QC |
| ARC_03c_v3_Statistical_Analysis.m | Statistical analysis (disease, subtype, severity, contrast specificity) |
| ARC_04_v3_Lesion_Network_Analysis.m | Lesion-network analysis (presence, disconnection, mediation) |
| ARC_05_v3_VolumeStratified_Analysis.mlx | Volume stratification with meta-analysis |
| ARC_06_v3_Smith10_Robustness.mlx | Cross-atlas validation (Smith10) |
| config_template.m | Configuration template |

### Configuration
`config_local.m` is user-created. Copy `config_template.m` to `config_local.m` and edit it with your local paths.

## Setup
1. Download ARC dataset from OpenNeuro
2. Copy `config_template.m` to `config_local.m`
3. Edit `config_local.m` with your local paths
4. Run scripts in numerical order
5. Consult data release on OpenNeuro for participants.txt

## What is not included
Raw imaging data, SPM first-level outputs, and demographic tables are not redistributed. To regenerate the full pipeline, download the ARC dataset from OpenNeuro and follow the setup steps in the README. Derived summary tables and pipeline scripts are provided for reproducibility.

## Citation
If you use this repository, please cite the ARC dataset:

> Gibson M, Newman-Norlund R, Bonilha L, Fridriksson J, Hickok G, Hillis AE, Den Ouden DB, Rorden C. The Aphasia Recovery Cohort, an open-source chronic stroke repository. *Scientific Data*. 2024;11(1):981.

Manuscript citation: [Manuscript in preparation]

## Contact
jkissane22@gmail.com

## Acknowledgements
ICN_Atlas toolbox: Kozák LR et al. (2017) NeuroImage 163:319–341

## License
MIT License
