%% ARC_03e_v3_UniqueVariance_Regression.m
% =========================================================================
% Hierarchical regression: unique variance of BM17 engagement beyond
% lesion volume and clinical covariates.
%
% Purpose:
%   Computes the incremental variance (Delta-R2) in WAB-R explained by
%   adding BM17 C2 IRi to a baseline model containing lesion volume, age,
%   sex, and days post-stroke.
%
%   This produces the Delta-R2 = 5.8% figure reported in the Abstract,
%   Results, and Discussion.
%
% Note on N:
%   The primary severity correlation (Partial Spearman) uses N=139.
%   This hierarchical model uses N=138, because one participant lacked a
%   lesion mask.
%
% Requires:
%   - config_local.m (copy of config_template.m)
%   - ARC_03b_v3_Master_Wide.mat (from Script 3b)
%   - ARC_04_v3_AllResults.mat  (lesion volume)
%
% Output:
%   - BM17_UniqueVariance_Results.txt
% =========================================================================

clear; clc;

% Load configuration
config = config_local();
rootDir = config.rootDir;

% Paths
wideFile   = fullfile(rootDir, 'ARC_03_Output', 'ARC_03b_v3_Output', 'ARC_03b_v3_Master_Wide.mat');
lesionFile = fullfile(rootDir, 'ARC_04_v3_AllResults.mat');
outputDir  = fullfile(rootDir, 'ARC_03_Output', 'ARC_03e_v3_Output');

if ~exist(outputDir, 'dir'), mkdir(outputDir); end

% Load data
load(wideFile, 'masterWide');
load(lesionFile, 'mergedTable', 'resampleLog');

% restOnlyList: anonymised ARC IDs of participants whose first-session
% functional acquisition was resting-state only. Hardcoded for portability;
% see branch README for details.
restOnlyList = {'M2088','M2097','M2100','M2101','M2113','M2114',...
    'M2117','M2118','M2122','M2126','M2129','M2131','M2135',...
    'M2140','M2141','M2142','M2143','M2144','M2145','M2146',...
    'M2149','M2150','M2151','M2152','M2153','M2155','M2156',...
    'M2158','M2159','M2160','M2162','M2164','M2165','M2169',...
    'M2184','M2254'};

isRestOnly = ismember(masterWide.PatientID, restOnlyList);
isAphasia  = strcmp(masterWide.GroupRole, 'Aphasia');
hasValidC2 = ~isnan(masterWide.C2_ICN17_IRi);
if ~ismember('SexNumeric', masterWide.Properties.VariableNames)
    masterWide.SexNumeric = double(strcmp(masterWide.Sex, 'M'));
end

clean = masterWide(isAphasia & hasValidC2 & ~isRestOnly, :);

% Attach lesion volume
[~, idxM, idxC] = intersect(mergedTable.PatientID, clean.PatientID, 'stable');
cleanLV = clean(idxC, :);
cleanLV.LesionVolume = mergedTable.LesionVolume_mm3(idxM);

% Drop missing
valid = ~isnan(cleanLV.LesionVolume) & ~isnan(cleanLV.WAB_AQ) & ...
    ~isnan(cleanLV.Age_At_Stroke) & ~isnan(cleanLV.SexNumeric) & ...
    ~isnan(cleanLV.Days_Post_Stroke) & ~isnan(cleanLV.C2_ICN17_IRi);

y   = cleanLV.WAB_AQ(valid);
lv  = cleanLV.LesionVolume(valid);
age = cleanLV.Age_At_Stroke(valid);
sex = cleanLV.SexNumeric(valid);
dps = cleanLV.Days_Post_Stroke(valid);
iri = cleanLV.C2_ICN17_IRi(valid);

fprintf('N for recomputation: %d\n', numel(y));

% Step 1: baseline model
X1 = [lv, age, sex, dps];
mdl1 = fitlm(X1, y);
R2_1 = mdl1.Rsquared.Ordinary;

% Step 2: add BM17 engagement
X2 = [lv, age, sex, dps, iri];
mdl2 = fitlm(X2, y);
R2_2 = mdl2.Rsquared.Ordinary;

deltaR2 = R2_2 - R2_1;

% Report
fprintf('R2 baseline (lesion + covariates): %.4f\n', R2_1);
fprintf('R2 with BM17 IRi:                  %.4f\n', R2_2);
fprintf('Delta R2 (unique variance):        %.4f\n', deltaR2);

% Save
outFile = fullfile(outputDir, 'BM17_UniqueVariance_Results.txt');
fid = fopen(outFile, 'w');
fprintf(fid, 'BM17 UNIQUE VARIANCE BEYOND LESION VOLUME AND COVARIATES\n');
fprintf(fid, '=======================================================\n');
fprintf(fid, 'N: %d\n', numel(y));
fprintf(fid, 'Baseline model: WAB-AQ ~ lesion volume + age + sex + days post-stroke\n');
fprintf(fid, '  R2 = %.4f\n', R2_1);
fprintf(fid, 'Full model: baseline + BM17 C2 IRi\n');
fprintf(fid, '  R2 = %.4f\n', R2_2);
fprintf(fid, 'Delta R2 (unique variance of BM17): %.4f\n', deltaR2);
fclose(fid);
fprintf('\nSaved: %s\n', outFile);
