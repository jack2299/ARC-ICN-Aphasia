# QC Corrective Scripts

This branch contains the quality-control corrective scripts that were run after discovering that 14 rest‑only aphasia participants had spurious ICN engagement values (non‑NaN zeros) in the original merged table.  
These scripts re‑run all contaminated analyses on the **clean aphasia sample** (rest‑only excluded, valid C2 IRi) and produce the corrected results used in the final manuscript.

## Scripts

| Script | Purpose |
|--------|---------|
| `motion_posthoc_extraction_rigorous.m` | Extracts mean rotation (degrees) and mean translation (mm) from realignment parameters for the clean aphasia sample. Saves `Motion_Extracted.mat`. |
| `motion_posthoc_rigorous.m` | Runs partial Spearman correlations between motion metrics and ICN engagement (IRi) for all 18 ICNs, with FDR correction. Generates `MOTION_REPORT.txt` and CSV outputs. |
| `MasterRerun_levels1_to_4_cleansample.m` | Re‑runs Levels 1‑4 (disease effect, subtype, severity, contrast specificity) on the clean sample. Saves `CleanSample_Results_Rigorous.mat`. |
| `MasterRerun_remainingLevels_cleansample.m` | Comprehensive re‑run of all remaining engagement‑based analyses: Level 8 (disconnection), Level 8b (specificity), Level 10b (bootstrap mediation), volume‑stratified engagement meta‑analysis, cross‑atlas validation (Smith10 S06), IRi threshold sensitivity. Produces all `QC_Corrected_*.csv` files and a comparison report. |

## Required data files

These scripts do **not** include the ARC dataset. You must obtain the following files from the [Aphasia Recovery Cohort](https://openneuro.org/datasets/ds004884) (OpenNeuro) or from the processed pipeline outputs:

- `ARC_03b_v3_Master_Wide.mat` – post‑IRi‑QC engagement master table  
- `ARC_04_v3_AllResults.mat` – resample log and merged lesion‑engagement table  
- `ARC_06_v3_AllResults.mat` – Smith10 atlas engagement data (for cross‑atlas validation)  
- `task_fmri_first_session.csv` – run mapping (from the ARC dataset)  
- `ds004884‑rp_files/` – folder containing `rp_*.txt` realignment parameter files  

Place these files in the same folder as the scripts before running.

## How to run

1. Run `motion_posthoc_extraction_rigorous.m` to generate `Motion_Extracted.mat`.  
2. Run `motion_posthoc_rigorous.m` to perform the motion sensitivity analysis.  
3. Run `MasterRerun_levels1_to_4_cleansample.m` to produce `CleanSample_Results_Rigorous.mat`.  
4. Run `MasterRerun_remainingLevels_cleansample.m` to generate all corrected CSV outputs and the comparison report.

All scripts are self‑contained and will save outputs to the current folder. They use the same methods, covariates, and FDR correction as the original pipeline.

## Output files

After running all scripts, the folder will contain:

- `Motion_Extracted.mat`  
- `Motion_Results_Rotation.csv`, `Motion_Results_Translation.csv`, `MOTION_REPORT.txt`  
- `CleanSample_Results_Rigorous.mat`  
- `QC_Corrected_Level8_IRi.csv`, `QC_Corrected_Level8_MANi.csv`, `QC_Corrected_Level8_logVari.csv`  
- `QC_Corrected_Level8b.csv`  
- `QC_Corrected_Level10b.csv`  
- `QC_Corrected_VolumeStrat_Engagement_IRi.csv`, `…_MANi.csv`, `…_logVari.csv`  
- `QC_Corrected_Smith10_S06.csv`  
- `QC_Corrected_IRi_Sensitivity.csv`  
- `QC_Corrected_Comparison_Report.txt`

These are the corrected results used in the manuscript.

## Note

The original (uncorrected) pipeline scripts remain on the `main` branch. This branch exists to provide a transparent, reproducible record of the quality‑control step and to allow others to regenerate the corrected numbers.
