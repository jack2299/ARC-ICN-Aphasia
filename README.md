# ARC-ICN-Aphasia

Analysis pipeline for: **ICN Engagement Predicts Aphasia Severity in Chronic Stroke**

## Overview
MATLAB scripts for analysing intrinsic connectivity network (ICN) engagement during covert picture naming in the Aphasia Recovery Cohort (ARC).

## Requirements
- MATLAB R2025a (or compatible)
- SPM12 (Wellcome Centre for Human Neuroimaging)
- ICN_Atlas toolbox (version r20180306; Kozák et al., 2017)

## Data
The Aphasia Recovery Cohort is publicly available on OpenNeuro:  
https://openneuro.org/datasets/ds004884

DOI: 10.18112/openneuro.ds004884.v1.0.2

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
| config_local.m | Local path configuration (template) |
| config_template.m | Configuration template |

## Setup
1. Download ARC dataset from OpenNeuro
2. Copy `config_template.m` to `config_local.m`
3. Edit `config_local.m` with your local paths
4. Run scripts in numerical order

## Citation
[Manuscript in preparation]

## Contact
jkissane22@gmail.com

## Acknowledgements
ICN_Atlas toolbox: Kozák LR et al. (2017) NeuroImage 163:319–341

## License
MIT License
