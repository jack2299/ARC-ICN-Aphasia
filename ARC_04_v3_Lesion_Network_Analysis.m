%% ARC_04_v3_Lesion_Network_Analysis.m 
% =========================================================================
% SCRIPT 4 v3: Lesion-Network Analysis
% =========================================================================
% Purpose:
%   Test mechanistic hypothesis: Lesion location determines aphasia presence;
%   residual network engagement determines severity.
%
% Analyses:
%   Level 7:  Lesion → Aphasia presence (logistic regression per ICN)
%   Level 7b: **NEW** Aphasia Occurrence Prediction (3 models comparing AUC)
%   Level 8:  Lesion → Engagement correlation (disconnection hypothesis)
%   Level 8b: Specificity test (within-network vs between-network)
%   Level 9:  Lesion → Severity (regression)
%   Level 10: Mediation (Lesion → Engagement → Severity)
%   Level 11: White matter disconnection (using JHU atlas)
%
% V3 CHANGES FROM V2:
%   1. Input from ARC_03b_v3_Output (IRi QC-cleaned data)
%   2. Output filenames updated to v3
%   3. **NEW Level 7b**: Aphasia Occurrence Prediction
%      - Model 1: Volume Only
%      - Model 2: Volume + Total ICN Damage
%      - Model 3: Volume + ICN-Specific Damage (LASSO)
%      - Compares AUC to test mass effect vs modularity
%   4. All other analyses unchanged (same logic as v2)
%
% METRICS:
%   - IRi (PRIMARY) - ICNiRelativeInvolvement
%   - MANi (SECONDARY) - NormalisedMeanICNiActivation
%   - logVari (EXPLORATORY) - log-transformed variance of Z-scores
%
% Inputs:
%   - Lesion masks: NIfTI folder (wsub-M*_lesion.nii)
%   - ICN atlas: atlas_BrainMap20.atl
%   - ICN engagement: ARC_03b_v3_Output/ARC_03b_v3_Master_Wide.mat
%   - JHU WM atlas: atlas_JHUWM2mm.atl
%
% Outputs:
%   - ARC_04_v3_Lesion_ICN_Overlap.mat/csv
%   - ARC_04_v3_Level7_LesionVsPresence.csv
%   - ARC_04_v3_Level7b_OccurrencePrediction.csv (NEW)
%   - ARC_04_v3_Level8_DamageEngagement_[Metric].csv
%   - ARC_04_v3_Level9_LesionSeverity.csv
%   - ARC_04_v3_Level10_Mediation_[Metric].csv
%   - ARC_04_v3_Level11_Disconnection.csv
%   - ARC_04_v3_Summary.txt
% =========================================================================

%% SETUP
clear; clc;

% Load configuration
config = config_local();

fprintf('=================================================================\n');
fprintf('ARC_04_v3_Lesion_Network_Analysis.m\n');
fprintf('Lesion-Network Mechanistic Analysis (V3 - IRi QC + Occurrence)\n');
fprintf('=================================================================\n');
fprintf('Started: %s\n\n', datestr(now));

fprintf('V3 CHANGES FROM V2:\n');
fprintf('  - Input from ARC_03b_v3_Output (IRi QC-cleaned)\n');
fprintf('  - NEW Level 7b: Aphasia Occurrence Prediction (3 models)\n');
fprintf('  - All output filenames updated to v3\n');
fprintf('  - Metrics analyzed:\n');
fprintf('      IRi (PRIMARY) - Relative Spatial Involvement\n');
fprintf('      MANi (SECONDARY) - Normalized Mean Activation\n');
fprintf('      logVari (EXPLORATORY) - Log Activation Variance\n\n');

%% DEFINE PATHS
rootDir = config.rootDir;
lesionDir = fullfile(rootDir, 'NIfTI');
atlasDir = config.icnAtlasPath;

% Input files - V3 ALIGNED
icnAtlasFile = fullfile(atlasDir, 'atlas_BrainMap20.atl');
jhuAtlasFile = fullfile(atlasDir, 'atlas_JHUWM2mm.atl');
masterWideFile = fullfile(rootDir, 'ARC_03_Output', 'ARC_03b_v3_Output', 'ARC_03b_v3_Master_Wide.mat');

% Output directory - V3
outputDir = fullfile(rootDir, 'ARC_04_v3_Output');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('Created output directory: %s\n', outputDir);
end

% Add SPM to path
addpath(config.spmPath);

%% DEFINE METRICS (V3 - ALL THREE)
metrics = {'IRi', 'MANi', 'logVari'};
metricLabels = {'Relative Spatial Involvement (PRIMARY)', ...
                'Normalized Mean Activation (SECONDARY)', ...
                'Log Activation Variance (EXPLORATORY)'};
primaryMetric = 'IRi';
nMetrics = length(metrics);

fprintf('Metrics to analyze:\n');
for m = 1:nMetrics
    if strcmp(metrics{m}, primaryMetric)
        fprintf('  %s - %s ***\n', metrics{m}, metricLabels{m});
    else
        fprintf('  %s - %s\n', metrics{m}, metricLabels{m});
    end
end
fprintf('\n');

%% BOOTSTRAP PARAMETERS (for Level 10 mediation)
nBootstrap = 5000;
bootAlpha = 0.05;
bootCI_method = 'BCa';  % Bias-corrected accelerated

fprintf('Bootstrap iterations for mediation: %d\n', nBootstrap);
fprintf('Bootstrap CI method: %s\n\n', bootCI_method);

%% LOAD MASTER DATA FIRST (to get trueICN_idx)
fprintf('-----------------------------------------------------------------\n');
fprintf('Loading ICN engagement data from Script 3b v3 (IRi QC-cleaned)\n');
fprintf('-----------------------------------------------------------------\n');

% Verify input file exists
if ~exist(masterWideFile, 'file')
    error('Master wide file not found: %s\nRun Script 3b_v3 first.', masterWideFile);
end

% Load ALL relevant variables including trueICN_idx
masterData = load(masterWideFile);
masterWide = masterData.masterWide;
trueICN_labels = masterData.trueICN_labels;
nICN = masterData.nICN;

% Check if trueICN_idx was saved
if isfield(masterData, 'trueICN_idx')
    trueICN_idx = masterData.trueICN_idx;
    fprintf('Loaded trueICN_idx from Script 3b v3: %s\n', mat2str(trueICN_idx));
else
    % Will need to identify from atlas labels
    trueICN_idx = [];
    fprintf('trueICN_idx not found in Script 3b v3 output - will identify from atlas\n');
end

% Check if IRi_threshold was saved (for documentation)
if isfield(masterData, 'IRi_threshold')
    IRi_threshold = masterData.IRi_threshold;
    fprintf('IRi QC threshold used: %.1f\n', IRi_threshold);
else
    IRi_threshold = 0.8;  % Default
    fprintf('IRi QC threshold (assumed): %.1f\n', IRi_threshold);
end

fprintf('Loaded masterWide: %d patients\n', height(masterWide));
fprintf('Number of ICNs: %d\n', nICN);
fprintf('ICN labels:\n');
for n = 1:min(5, nICN)
    fprintf('  ICN%02d: %s\n', n, trueICN_labels{n});
end
if nICN > 5
    fprintf('  ... and %d more\n', nICN - 5);
end

%% LOG-TRANSFORM VARI (same as Script 3c - columns not saved by 3b)
% =========================================================================
fprintf('\nLog-transforming Vari columns (replicating Script 3c transformation)...\n');

% First verify Vari columns exist
sampleVari = sprintf('C2_ICN%02d_Vari', 1);
if ~ismember(sampleVari, masterWide.Properties.VariableNames)
    error('Vari columns not found in masterWide. Check Script 3b output.');
end

% Add log-transformed Vari columns
for n = 1:nICN
    % C1
    colC1 = sprintf('C1_ICN%02d_Vari', n);
    colC1_log = sprintf('C1_ICN%02d_logVari', n);
    if ismember(colC1, masterWide.Properties.VariableNames)
        masterWide.(colC1_log) = log(masterWide.(colC1) + 1);
    end
    
    % C2
    colC2 = sprintf('C2_ICN%02d_Vari', n);
    colC2_log = sprintf('C2_ICN%02d_logVari', n);
    if ismember(colC2, masterWide.Properties.VariableNames)
        masterWide.(colC2_log) = log(masterWide.(colC2) + 1);
    end
end

fprintf('Added logVari columns to masterWide (%d ICNs × 2 contrasts)\n', nICN);

% Verify all metric columns exist
fprintf('\nVerifying metric columns exist...\n');
for m = 1:nMetrics
    metric = metrics{m};
    sampleCol = sprintf('C2_ICN01_%s', metric);
    if ~ismember(sampleCol, masterWide.Properties.VariableNames)
        error('Expected column %s not found in masterWide. Check metric name alignment.', sampleCol);
    else
        fprintf('  %s columns: ✓\n', metric);
    end
end
fprintf('\n');

%% LOAD ICN ATLAS
fprintf('-----------------------------------------------------------------\n');
fprintf('Loading ICN Atlas (BrainMap20)\n');
fprintf('-----------------------------------------------------------------\n');

atl = load(icnAtlasFile, '-mat');
icnDataFull = atl.atl_mapdatamatrix;  % 91 x 109 x 91 x 20
icnLabelsFull = atl.atl_maplabels;

fprintf('Atlas dimensions: %s\n', mat2str(size(icnDataFull)));
fprintf('Total atlas components: %d\n', length(icnLabelsFull));

% Display all atlas labels
fprintf('Atlas labels:\n');
for n = 1:length(icnLabelsFull)
    fprintf('  %2d: %s\n', n, icnLabelsFull{n});
end
fprintf('\n');

% Identify true ICN indices if not loaded from Script 3b
if isempty(trueICN_idx)
    fprintf('Identifying true ICNs by label pattern...\n');
    
    artifactIdx = [];
    for n = 1:length(icnLabelsFull)
        label = lower(icnLabelsFull{n});
        if contains(label, 'artifact') || contains(label, 'outside') || ...
           contains(label, 'noise') || contains(label, 'edge') || ...
           contains(label, 'artefact')
            artifactIdx(end+1) = n;
            fprintf('  Artifact at index %d: %s\n', n, icnLabelsFull{n});
        end
    end
    
    trueICN_idx = setdiff(1:size(icnDataFull, 4), artifactIdx);
    fprintf('True ICN indices identified: %s\n', mat2str(trueICN_idx));
end

% Verify we have the expected number of ICNs
if length(trueICN_idx) ~= nICN
    warning('Mismatch! trueICN_idx has %d elements but nICN = %d', length(trueICN_idx), nICN);
    fprintf('Using first %d indices from trueICN_idx\n', nICN);
    trueICN_idx = trueICN_idx(1:nICN);
end

% Extract only true ICN data
icnData = icnDataFull(:,:,:,trueICN_idx);
icnLabels = icnLabelsFull(trueICN_idx);

fprintf('Using %d true ICNs (indices: %s)\n\n', nICN, mat2str(trueICN_idx));

% Verify alignment with Script 3b labels
fprintf('Verifying label alignment with Script 3b v3:\n');
for n = 1:min(3, nICN)
    atlasLabel = icnLabels{n};
    script3bLabel = trueICN_labels{n};
    fprintf('  ICN%02d: Atlas="%s" | 3b="%s"\n', n, atlasLabel(1:min(30,end)), script3bLabel(1:min(30,end)));
end
fprintf('\n');

% Atlas spatial info
atlasDim = [91, 109, 91];
atlasVoxelSize = 2;  % mm

%% LOAD JHU WHITE MATTER ATLAS (for Level 11)
fprintf('-----------------------------------------------------------------\n');
fprintf('Loading JHU White Matter Atlas\n');
fprintf('-----------------------------------------------------------------\n');

if exist(jhuAtlasFile, 'file')
    jhu = load(jhuAtlasFile, '-mat');
    jhuDataRaw = jhu.atl_mapdatamatrix;
    jhuLabels = jhu.atl_maplabels;
    
    fprintf('JHU raw data dimensions: %s\n', mat2str(size(jhuDataRaw)));
    fprintf('JHU tract labels: %d\n', length(jhuLabels));
    
    % Check if it's a labeled volume (3D) or probability maps (4D)
    if ndims(jhuDataRaw) == 3
        fprintf('JHU atlas is a LABELED VOLUME (integer labels per voxel)\n');
        
        % Find unique non-zero labels
        uniqueLabels = unique(jhuDataRaw(:));
        uniqueLabels = uniqueLabels(uniqueLabels > 0);
        nTracts = length(uniqueLabels);
        fprintf('Unique tract labels found: %d\n', nTracts);
        
        % Convert to binary masks (4D)
        fprintf('Converting to binary tract masks...\n');
        jhuData = false(size(jhuDataRaw, 1), size(jhuDataRaw, 2), size(jhuDataRaw, 3), nTracts);
        jhuTractLabels = cell(nTracts, 1);
        
        for t = 1:nTracts
            jhuData(:,:,:,t) = (jhuDataRaw == uniqueLabels(t));
            if uniqueLabels(t) <= length(jhuLabels)
                jhuTractLabels{t} = jhuLabels{uniqueLabels(t)};
            else
                jhuTractLabels{t} = sprintf('Tract_%d', uniqueLabels(t));
            end
        end
        
        jhuLabels = jhuTractLabels;
        fprintf('Converted to %d binary tract masks\n', nTracts);
        hasJHU = true;
        
    elseif ndims(jhuDataRaw) == 4
        fprintf('JHU atlas is already 4D probability maps\n');
        jhuData = jhuDataRaw;
        nTracts = size(jhuData, 4);
        hasJHU = true;
    else
        fprintf('Unexpected JHU atlas format (ndims = %d)\n', ndims(jhuDataRaw));
        hasJHU = false;
    end
    
    if hasJHU
        fprintf('JHU atlas ready: %d tracts\n\n', nTracts);
    end
else
    fprintf('JHU atlas not found - Level 11 will be skipped\n\n');
    hasJHU = false;
    nTracts = 0;
end

%% INVENTORY LESION FILES
fprintf('-----------------------------------------------------------------\n');
fprintf('Inventorying lesion files\n');
fprintf('-----------------------------------------------------------------\n');

lesionFiles = dir(fullfile(lesionDir, 'wsub-M*_lesion.nii'));
fprintf('Lesion files found: %d\n', length(lesionFiles));

% Show sample filenames for verification
fprintf('Sample filenames:\n');
for i = 1:min(5, length(lesionFiles))
    fprintf('  %s\n', lesionFiles(i).name);
end
fprintf('\n');

% Initialize table properly
nLesionFiles = length(lesionFiles);
lesionInventory = table('Size', [nLesionFiles, 5], ...
    'VariableTypes', {'cell', 'cell', 'cell', 'cell', 'cell'}, ...
    'VariableNames', {'PatientID', 'SessionInfo', 'FileName', 'FilePath', 'ParseStatus'});

% Parse each filename (handle both patterns)
nParsed = 0;
nUnparsed = 0;

for i = 1:nLesionFiles
    fname = lesionFiles(i).name;
    
    % Pattern 1: wsub-M####_ses-####x####_lesion.nii (with x separator)
    tokens = regexp(fname, 'wsub-(M\d+)_ses-(\d+x\d+)_lesion\.nii', 'tokens');
    
    if isempty(tokens)
        % Pattern 2: wsub-M####_ses-####_lesion.nii (no x separator)
        tokens = regexp(fname, 'wsub-(M\d+)_ses-(\d+)_lesion\.nii', 'tokens');
    end
    
    if ~isempty(tokens)
        lesionInventory.PatientID{i} = tokens{1}{1};
        lesionInventory.SessionInfo{i} = tokens{1}{2};
        lesionInventory.FileName{i} = fname;
        lesionInventory.FilePath{i} = fullfile(lesionFiles(i).folder, fname);
        lesionInventory.ParseStatus{i} = 'OK';
        nParsed = nParsed + 1;
    else
        lesionInventory.PatientID{i} = '';
        lesionInventory.SessionInfo{i} = '';
        lesionInventory.FileName{i} = fname;
        lesionInventory.FilePath{i} = fullfile(lesionFiles(i).folder, fname);
        lesionInventory.ParseStatus{i} = 'UNPARSED';
        nUnparsed = nUnparsed + 1;
        fprintf('  Warning: Could not parse: %s\n', fname);
    end
end

fprintf('Parsed successfully: %d\n', nParsed);
fprintf('Unparsed: %d\n', nUnparsed);

% Remove unparsed rows
lesionInventory = lesionInventory(strcmp(lesionInventory.ParseStatus, 'OK'), :);
fprintf('Valid lesion entries: %d\n', height(lesionInventory));

% Get unique patients
uniqueLesionPatients = unique(lesionInventory.PatientID);
fprintf('Unique patients with lesions: %d\n\n', length(uniqueLesionPatients));

%% MATCH PATIENTS BETWEEN DATASETS
fprintf('-----------------------------------------------------------------\n');
fprintf('Matching patients between lesion and engagement data\n');
fprintf('-----------------------------------------------------------------\n');

analysisPatients = masterWide.PatientID;

fprintf('Sample patient IDs from masterWide: %s\n', strjoin(analysisPatients(1:min(5,end)), ', '));
fprintf('Sample patient IDs from lesions: %s\n', strjoin(uniqueLesionPatients(1:min(5,end)), ', '));

% Find intersection
matchedPatients = intersect(analysisPatients, uniqueLesionPatients);
fprintf('\nMatched patients: %d\n', length(matchedPatients));

if isempty(matchedPatients)
    error('NO PATIENTS MATCHED! Check patient ID format between datasets.');
end

% Show some matched IDs
fprintf('Sample matched IDs: %s\n\n', strjoin(matchedPatients(1:min(5,end)), ', '));

%% CREATE RESAMPLED LESION DIRECTORY
fprintf('-----------------------------------------------------------------\n');
fprintf('Resampling lesion masks to 2mm resolution\n');
fprintf('-----------------------------------------------------------------\n');

resampledDir = fullfile(outputDir, 'Lesions_2mm');
if ~exist(resampledDir, 'dir')
    mkdir(resampledDir);
end

% Target dimensions (2mm MNI)
targetDim = atlasDim;  % [91, 109, 91]

% Get transformation matrix from atlas if available
if isfield(atl, 'atl_Ta')
    targetMat = atl.atl_Ta;
    fprintf('Using transformation matrix from atlas\n');
else
    targetMat = [-2 0 0 92; 0 2 0 -128; 0 0 2 -74; 0 0 0 1];  % Standard 2mm MNI
    fprintf('Using standard 2mm MNI transformation matrix\n');
end

% Process each matched patient
nMatched = length(matchedPatients);
resampleLog = table('Size', [nMatched, 6], ...
    'VariableTypes', {'cell', 'cell', 'cell', 'cell', 'double', 'double'}, ...
    'VariableNames', {'PatientID', 'OriginalFile', 'ResampledFile', 'Status', 'OriginalVoxels', 'ResampledVoxels'});

fprintf('Resampling %d lesion masks...\n', nMatched);
progressStep = max(1, round(nMatched / 10));

for p = 1:nMatched
    pid = matchedPatients{p};
    
    % Progress
    if mod(p, progressStep) == 0
        fprintf('  %d/%d (%.0f%%)\n', p, nMatched, 100*p/nMatched);
    end
    
    % Find lesion file for this patient
    lesionIdx = find(strcmp(lesionInventory.PatientID, pid), 1);
    
    resampleLog.PatientID{p} = pid;
    
    if isempty(lesionIdx)
        resampleLog.Status{p} = 'No lesion file';
        resampleLog.OriginalFile{p} = '';
        resampleLog.ResampledFile{p} = '';
        continue;
    end
    
    lesionPath = lesionInventory.FilePath{lesionIdx};
    resampledPath = fullfile(resampledDir, sprintf('%s_lesion_2mm.nii', pid));
    
    resampleLog.OriginalFile{p} = lesionPath;
    resampleLog.ResampledFile{p} = resampledPath;
    
    % Skip if already resampled
    if exist(resampledPath, 'file')
        resampleLog.Status{p} = 'Already exists';
        
        % Load to get voxel counts
        try
            V = spm_vol(resampledPath);
            vol = spm_read_vols(V);
            resampleLog.ResampledVoxels(p) = sum(vol(:) > 0.5);
        catch
            resampleLog.ResampledVoxels(p) = NaN;
        end
        continue;
    end
    
    try
        % Load original lesion mask
        V = spm_vol(lesionPath);
        lesionOrig = spm_read_vols(V);
        
        resampleLog.OriginalVoxels(p) = sum(lesionOrig(:) > 0);
        
        % Resample to 2mm (91x109x91)
        if exist('imresize3', 'file')
            lesionResampled = imresize3(double(lesionOrig), targetDim, 'nearest');
        else
            % Manual resampling using interp3
            [X, Y, Z] = ndgrid(linspace(1, size(lesionOrig,1), targetDim(1)), ...
                               linspace(1, size(lesionOrig,2), targetDim(2)), ...
                               linspace(1, size(lesionOrig,3), targetDim(3)));
            lesionResampled = interp3(double(lesionOrig), Y, X, Z, 'nearest', 0);
        end
        
        % Binarize
        lesionResampled = double(lesionResampled > 0.5);
        
        resampleLog.ResampledVoxels(p) = sum(lesionResampled(:) > 0);
        
        % Create output header (2mm MNI space)
        Vout = V;
        Vout.fname = resampledPath;
        Vout.dim = targetDim;
        Vout.mat = targetMat;
        Vout.dt = [spm_type('float32') 0];
        
        % Write resampled volume
        spm_write_vol(Vout, lesionResampled);
        
        resampleLog.Status{p} = 'Success';
        
    catch ME
        resampleLog.Status{p} = sprintf('Error: %s', ME.message);
        resampleLog.OriginalVoxels(p) = NaN;
        resampleLog.ResampledVoxels(p) = NaN;
    end
end

fprintf('  %d/%d (100%%)\n', nMatched, nMatched);

% Summary
nSuccess = sum(strcmp(resampleLog.Status, 'Success') | strcmp(resampleLog.Status, 'Already exists'));
nFailed = sum(contains(resampleLog.Status, 'Error'));
fprintf('Successfully resampled: %d/%d\n', nSuccess, nMatched);
if nFailed > 0
    fprintf('Failed: %d\n', nFailed);
end
fprintf('\n');

%% COMPUTE LESION-ICN OVERLAP
fprintf('-----------------------------------------------------------------\n');
fprintf('Computing lesion-ICN overlap\n');
fprintf('-----------------------------------------------------------------\n');

% Initialize overlap matrix: patients × ICNs
overlapMatrix = NaN(nMatched, nICN);
lesionVolumeVec = NaN(nMatched, 1);
patientIDs = matchedPatients;

% Threshold ICN maps (Z > 3 typically used)
icnThreshold = 3;
icnMasks = icnData > icnThreshold;  % Binary masks for each ICN

% Calculate total voxels per ICN (for reference)
icnVoxelCounts = zeros(nICN, 1);
for n = 1:nICN
    icnVoxelCounts(n) = sum(sum(sum(icnMasks(:,:,:,n))));
end
fprintf('ICN voxel counts (Z>3):\n');
for n = 1:min(5, nICN)
    fprintf('  ICN%02d: %d voxels\n', n, icnVoxelCounts(n));
end
fprintf('  ... (mean: %.0f voxels)\n\n', mean(icnVoxelCounts));

fprintf('Computing overlap for %d patients × %d ICNs...\n', nMatched, nICN);

for p = 1:nMatched
    pid = patientIDs{p};
    resampledPath = fullfile(resampledDir, sprintf('%s_lesion_2mm.nii', pid));
    
    if ~exist(resampledPath, 'file')
        continue;
    end
    
    try
        % Load resampled lesion
        V = spm_vol(resampledPath);
        lesion = spm_read_vols(V) > 0.5;
        
        % Total lesion volume (in mm³)
        lesionVolumeVec(p) = sum(lesion(:)) * (atlasVoxelSize^3);
        
        % Check dimensions match
        if ~isequal(size(lesion), atlasDim)
            warning('Dimension mismatch for patient %s: lesion=%s, atlas=%s', ...
                pid, mat2str(size(lesion)), mat2str(atlasDim));
            continue;
        end
        
        % Compute overlap with each ICN
        for n = 1:nICN
            icnMask = icnMasks(:,:,:,n);
            
            % Overlap = proportion of ICN damaged by lesion
            icnVoxels = icnVoxelCounts(n);
            overlapVoxels = sum(lesion(:) & icnMask(:));
            
            if icnVoxels > 0
                overlapMatrix(p, n) = overlapVoxels / icnVoxels;
            else
                overlapMatrix(p, n) = 0;
            end
        end
        
    catch ME
        warning('Error processing patient %s: %s', pid, ME.message);
    end
    
    % Progress
    if mod(p, progressStep) == 0
        fprintf('  %d/%d (%.0f%%)\n', p, nMatched, 100*p/nMatched);
    end
end

fprintf('  %d/%d (100%%)\n', nMatched, nMatched);

% Summary statistics
validLesions = ~isnan(lesionVolumeVec);
fprintf('\nLesion-ICN overlap summary:\n');
fprintf('  Patients with valid data: %d\n', sum(validLesions));
fprintf('  Mean lesion volume: %.0f mm³ (SD: %.0f)\n', ...
    mean(lesionVolumeVec(validLesions)), std(lesionVolumeVec(validLesions)));
fprintf('  Mean ICN overlap: %.4f (SD: %.4f)\n', ...
    mean(overlapMatrix(:), 'omitnan'), std(overlapMatrix(:), 'omitnan'));
fprintf('  Max ICN overlap: %.4f\n', max(overlapMatrix(:)));
fprintf('\n');

%% CREATE LESION-ICN OVERLAP TABLE WITH ALL METRICS
fprintf('-----------------------------------------------------------------\n');
fprintf('Creating merged lesion overlap table (all metrics)\n');
fprintf('-----------------------------------------------------------------\n');

% Build table by merging lesion data with clinical/engagement data
% Pre-allocate with proper size
nValid = sum(validLesions);
mergedTable = table('Size', [nValid, 0]);

% Find valid patients and merge
mergedIdx = 0;
for p = 1:nMatched
    pid = patientIDs{p};
    masterIdx = find(strcmp(masterWide.PatientID, pid), 1);
    
    if ~isempty(masterIdx) && ~isnan(lesionVolumeVec(p))
        mergedIdx = mergedIdx + 1;
        
        % Patient ID and lesion data
        mergedTable.PatientID{mergedIdx} = pid;
        mergedTable.LesionVolume_mm3(mergedIdx) = lesionVolumeVec(p);
        
        % ICN damage
        for n = 1:nICN
            mergedTable.(sprintf('ICN%02d_Damage', n))(mergedIdx) = overlapMatrix(p, n);
        end
        
        % Clinical data from master table
        mergedTable.AphasiaType{mergedIdx} = masterWide.AphasiaType{masterIdx};
        mergedTable.GroupRole{mergedIdx} = masterWide.GroupRole{masterIdx};
        mergedTable.WAB_AQ(mergedIdx) = masterWide.WAB_AQ(masterIdx);
        mergedTable.Age_At_Stroke(mergedIdx) = masterWide.Age_At_Stroke(masterIdx);
        mergedTable.Sex{mergedIdx} = masterWide.Sex{masterIdx};
        mergedTable.Days_Post_Stroke(mergedIdx) = masterWide.Days_Post_Stroke(masterIdx);
        
        % ICN engagement - ALL THREE METRICS
        for m = 1:nMetrics
            metric = metrics{m};
            for n = 1:nICN
                engCol = sprintf('C2_ICN%02d_%s', n, metric);
                outCol = sprintf('ICN%02d_Eng_%s', n, metric);
                if ismember(engCol, masterWide.Properties.VariableNames)
                    mergedTable.(outCol)(mergedIdx) = masterWide.(engCol)(masterIdx);
                else
                    mergedTable.(outCol)(mergedIdx) = NaN;
                end
            end
        end
    end
end

fprintf('Merged table: %d patients with complete lesion + engagement data\n', height(mergedTable));
fprintf('Engagement columns per ICN: %d metrics × %d ICNs = %d columns\n\n', nMetrics, nICN, nMetrics*nICN);

% Convert Sex to numeric for regression
mergedTable.SexNumeric = double(strcmp(mergedTable.Sex, 'M'));

% Define group indices
isAphasia = strcmp(mergedTable.GroupRole, 'Aphasia');
isControl = strcmp(mergedTable.GroupRole, 'Stroke Control (No Aphasia)');

fprintf('Group breakdown in merged data:\n');
fprintf('  Aphasia patients: %d\n', sum(isAphasia));
fprintf('  Stroke controls: %d\n\n', sum(isControl));

%% HELPER FUNCTION: FDR CORRECTION
    function [h, crit_p, adj_p] = fdr_bh(pvals, q)
        if nargin < 2, q = 0.05; end
        
        pvals = pvals(:);
        m = length(pvals);
        
        % Handle all NaN case
        if all(isnan(pvals))
            h = false(size(pvals));
            crit_p = 0;
            adj_p = NaN(size(pvals));
            return;
        end
        
        [sorted_p, sort_idx] = sort(pvals);
        
        thresh = (1:m)' / m * q;
        below = sorted_p <= thresh & ~isnan(sorted_p);
        
        if any(below)
            crit_p = max(sorted_p(below));
        else
            crit_p = 0;
        end
        
        adj_p = NaN(m, 1);
        validP = ~isnan(sorted_p);
        adj_p(sort_idx(validP)) = sorted_p(validP) .* m ./ find(validP);
        
        for i = m-1:-1:1
            if ~isnan(adj_p(sort_idx(i))) && ~isnan(adj_p(sort_idx(min(i+1,m))))
                adj_p(sort_idx(i)) = min(adj_p(sort_idx(i)), adj_p(sort_idx(i+1)));
            end
        end
        adj_p = min(adj_p, 1);
        
        h = pvals <= crit_p & ~isnan(pvals);
    end

%% HELPER FUNCTION: COMPUTE AUC
    function auc = computeAUC(labels, scores)
        % Compute AUC using trapezoidal rule (no toolbox required)
        % labels: binary (0/1)
        % scores: continuous predictions
        
        % Handle NaN
        validIdx = ~isnan(labels) & ~isnan(scores);
        labels = labels(validIdx);
        scores = scores(validIdx);
        
        if isempty(labels) || length(unique(labels)) < 2
            auc = NaN;
            return;
        end
        
        % Sort by scores descending
        [~, sortIdx] = sort(scores, 'descend');
        labels = labels(sortIdx);
        
        % Count positives and negatives
        nPos = sum(labels == 1);
        nNeg = sum(labels == 0);
        
        if nPos == 0 || nNeg == 0
            auc = NaN;
            return;
        end
        
        % Compute AUC via Mann-Whitney U statistic
        tpr = cumsum(labels == 1) / nPos;
        fpr = cumsum(labels == 0) / nNeg;
        
        % Trapezoidal integration
        auc = trapz([0; fpr; 1], [0; tpr; 1]);
    end

%% =========================================================================
%% LEVEL 7: LESION → APHASIA PRESENCE (LOGISTIC REGRESSION PER ICN)
%% =========================================================================
fprintf('=================================================================\n');
fprintf('LEVEL 7: Lesion Damage → Aphasia Presence (Per ICN)\n');
fprintf('=================================================================\n');

% Outcome: Aphasia (1) vs Stroke Control (0)
Y_L7 = double(isAphasia);

% Covariates
age = mergedTable.Age_At_Stroke;
sex = mergedTable.SexNumeric;
days = mergedTable.Days_Post_Stroke;
lesionVol = mergedTable.LesionVolume_mm3;

% Check for sufficient variation in outcome
fprintf('Outcome distribution: Aphasia=%d, Control=%d\n', sum(Y_L7==1), sum(Y_L7==0));

if sum(Y_L7==0) < 5 || sum(Y_L7==1) < 5
    warning('Insufficient variation in outcome for logistic regression');
end

% Storage for Level 7 results
Level7_Results = table('Size', [nICN, 9], ...
    'VariableTypes', {'double', 'cell', 'double', 'double', 'double', 'double', 'double', 'double', 'logical'}, ...
    'VariableNames', {'ICN', 'ICN_Label', 'Beta', 'SE', 'Z', 'p_uncorrected', 'OddsRatio', 'N', 'Significant_FDR'});

for n = 1:nICN
    damageCol = sprintf('ICN%02d_Damage', n);
    damage = mergedTable.(damageCol);
    
    Level7_Results.ICN(n) = n;
    Level7_Results.ICN_Label{n} = trueICN_labels{n};
    
    % Valid observations
    validIdx = ~isnan(damage) & ~isnan(age) & ~isnan(sex) & ~isnan(days) & ~isnan(lesionVol);
    
    % Check for variation in predictor
    if sum(validIdx) < 20 || var(damage(validIdx)) < 1e-10
        Level7_Results.Beta(n) = NaN;
        Level7_Results.SE(n) = NaN;
        Level7_Results.Z(n) = NaN;
        Level7_Results.p_uncorrected(n) = NaN;
        Level7_Results.OddsRatio(n) = NaN;
        Level7_Results.N(n) = sum(validIdx);
        continue;
    end
    
    % Predictors: ICN damage + covariates (including lesion volume)
    X = [damage(validIdx), lesionVol(validIdx), age(validIdx), sex(validIdx), days(validIdx)];
    Y = Y_L7(validIdx);
    
    try
        % Fit logistic regression
        [B, ~, stats] = glmfit(X, Y, 'binomial', 'link', 'logit');
        
        % ICN damage coefficient (index 2, after intercept)
        beta_icn = B(2);
        se_icn = stats.se(2);
        z_icn = beta_icn / se_icn;
        p_icn = 2 * (1 - normcdf(abs(z_icn)));
        or_icn = exp(beta_icn);
        
        Level7_Results.Beta(n) = beta_icn;
        Level7_Results.SE(n) = se_icn;
        Level7_Results.Z(n) = z_icn;
        Level7_Results.p_uncorrected(n) = p_icn;
        Level7_Results.OddsRatio(n) = or_icn;
        Level7_Results.N(n) = sum(validIdx);
        
    catch ME
        fprintf('  ICN%02d: Logistic regression failed - %s\n', n, ME.message);
        Level7_Results.Beta(n) = NaN;
        Level7_Results.SE(n) = NaN;
        Level7_Results.Z(n) = NaN;
        Level7_Results.p_uncorrected(n) = NaN;
        Level7_Results.OddsRatio(n) = NaN;
        Level7_Results.N(n) = sum(validIdx);
    end
end

% FDR correction
[Level7_Results.Significant_FDR, ~, Level7_Results.p_FDR] = fdr_bh(Level7_Results.p_uncorrected, 0.05);

% Sort by p-value
Level7_Results = sortrows(Level7_Results, 'p_uncorrected');

% Display results
fprintf('\nLevel 7 Results: ICN Damage → Aphasia Presence\n');
fprintf('Controlling for: Lesion Volume, Age, Sex, Days Post-Stroke\n\n');

sigL7 = Level7_Results(Level7_Results.Significant_FDR == true, :);
fprintf('Significant predictors (FDR q<0.05): %d / %d\n\n', height(sigL7), nICN);

if height(sigL7) > 0
    fprintf('%-6s %-35s %8s %8s %8s\n', 'ICN', 'Label', 'OR', 'Z', 'p_FDR');
    fprintf('%s\n', repmat('-', 1, 70));
    for i = 1:height(sigL7)
        fprintf('ICN%02d  %-35s %8.2f %8.2f %8.4f\n', ...
            sigL7.ICN(i), sigL7.ICN_Label{i}(1:min(35,end)), ...
            sigL7.OddsRatio(i), sigL7.Z(i), sigL7.p_FDR(i));
    end
else
    fprintf('No significant predictors after FDR correction.\n');
end
fprintf('\n');

%% =========================================================================
%% LEVEL 7b: APHASIA OCCURRENCE PREDICTION (NEW IN V3)
%% =========================================================================
fprintf('=================================================================\n');
fprintf('LEVEL 7b: Aphasia Occurrence Prediction (NEW IN V3)\n');
fprintf('         Comparing Mass Effect vs Modularity\n');
fprintf('=================================================================\n');
fprintf('Question: What predicts developing aphasia after left hemisphere stroke?\n');
fprintf('Models:\n');
fprintf('  1. Volume Only: Lesion volume alone\n');
fprintf('  2. Volume + Total ICN Damage: Adds sum of all ICN overlaps\n');
fprintf('  3. Volume + ICN-Specific: LASSO across individual ICNs\n\n');

% Outcome
Y_7b = double(isAphasia);
nTotal = length(Y_7b);

% Prepare data - ensure no NaN in covariates
validIdx_7b = ~isnan(mergedTable.LesionVolume_mm3) & ...
              ~isnan(mergedTable.Age_At_Stroke) & ...
              ~isnan(mergedTable.SexNumeric) & ...
              ~isnan(mergedTable.Days_Post_Stroke);

% Also check all ICN damage columns are valid
for n = 1:nICN
    damageCol = sprintf('ICN%02d_Damage', n);
    validIdx_7b = validIdx_7b & ~isnan(mergedTable.(damageCol));
end

fprintf('Valid observations for occurrence prediction: %d/%d\n', sum(validIdx_7b), nTotal);
fprintf('  Aphasia: %d, Control: %d\n\n', sum(Y_7b(validIdx_7b)==1), sum(Y_7b(validIdx_7b)==0));

% Extract valid data
Y_valid = Y_7b(validIdx_7b);
lesionVol_valid = mergedTable.LesionVolume_mm3(validIdx_7b);
age_valid = mergedTable.Age_At_Stroke(validIdx_7b);
sex_valid = mergedTable.SexNumeric(validIdx_7b);
days_valid = mergedTable.Days_Post_Stroke(validIdx_7b);

% Create ICN damage matrix
icnDamage_valid = NaN(sum(validIdx_7b), nICN);
for n = 1:nICN
    damageCol = sprintf('ICN%02d_Damage', n);
    icnDamage_valid(:, n) = mergedTable.(damageCol)(validIdx_7b);
end

% Total ICN damage (sum of all overlaps)
totalICNDamage = sum(icnDamage_valid, 2);

% Standardize continuous predictors for model stability
lesionVol_z = (lesionVol_valid - mean(lesionVol_valid)) / std(lesionVol_valid);
totalICN_z = (totalICNDamage - mean(totalICNDamage)) / std(totalICNDamage);
age_z = (age_valid - mean(age_valid)) / std(age_valid);
days_z = (days_valid - mean(days_valid)) / std(days_valid);

% Standardize ICN columns
icnDamage_z = (icnDamage_valid - mean(icnDamage_valid, 1)) ./ std(icnDamage_valid, 0, 1);
icnDamage_z(isnan(icnDamage_z)) = 0;  % Handle zero-variance columns

% Initialize results storage
Level7b_Results = table('Size', [3, 10], ...
    'VariableTypes', {'cell', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'cell'}, ...
    'VariableNames', {'Model', 'N', 'N_Aphasia', 'N_Control', 'N_Predictors', 'AUC', 'AUC_SE', ...
                      'LogLikelihood', 'AIC', 'Interpretation'});

% -------------------------------------------------------------------------
% MODEL 1: Volume Only
% -------------------------------------------------------------------------
fprintf('--- Model 1: Volume Only ---\n');
try
    X_m1 = [lesionVol_z, age_z, sex_valid, days_z];
    [B_m1, ~, stats_m1] = glmfit(X_m1, Y_valid, 'binomial', 'link', 'logit');
    
    % Predicted probabilities
    eta_m1 = [ones(length(Y_valid), 1), X_m1] * B_m1;
    prob_m1 = 1 ./ (1 + exp(-eta_m1));
    
    % AUC
    auc_m1 = computeAUC(Y_valid, prob_m1);
    
    % Log-likelihood and AIC
    ll_m1 = sum(Y_valid .* log(prob_m1 + 1e-10) + (1 - Y_valid) .* log(1 - prob_m1 + 1e-10));
    k_m1 = length(B_m1);
    aic_m1 = -2 * ll_m1 + 2 * k_m1;
    
    % AUC SE (DeLong approximation)
    nPos = sum(Y_valid == 1);
    nNeg = sum(Y_valid == 0);
    auc_se_m1 = sqrt((auc_m1 * (1 - auc_m1) + (nPos - 1) * (auc_m1 / (2 - auc_m1) - auc_m1^2) + ...
                      (nNeg - 1) * (2 * auc_m1^2 / (1 + auc_m1) - auc_m1^2)) / (nPos * nNeg));
    
    fprintf('  AUC = %.3f (SE = %.3f)\n', auc_m1, auc_se_m1);
    fprintf('  AIC = %.1f\n', aic_m1);
    fprintf('  Volume beta = %.3f, p = %.4f\n', B_m1(2), stats_m1.p(2));
    
    Level7b_Results.Model{1} = 'Volume Only';
    Level7b_Results.N(1) = length(Y_valid);
    Level7b_Results.N_Aphasia(1) = sum(Y_valid == 1);
    Level7b_Results.N_Control(1) = sum(Y_valid == 0);
    Level7b_Results.N_Predictors(1) = k_m1 - 1;  % Exclude intercept
    Level7b_Results.AUC(1) = auc_m1;
    Level7b_Results.AUC_SE(1) = auc_se_m1;
    Level7b_Results.LogLikelihood(1) = ll_m1;
    Level7b_Results.AIC(1) = aic_m1;
    Level7b_Results.Interpretation{1} = 'Baseline: mass effect only';
    
catch ME
    fprintf('  ERROR: %s\n', ME.message);
    Level7b_Results.Model{1} = 'Volume Only';
    Level7b_Results.AUC(1) = NaN;
end
fprintf('\n');

% -------------------------------------------------------------------------
% MODEL 2: Volume + Total ICN Damage
% -------------------------------------------------------------------------
fprintf('--- Model 2: Volume + Total ICN Damage ---\n');
try
    X_m2 = [lesionVol_z, totalICN_z, age_z, sex_valid, days_z];
    [B_m2, ~, stats_m2] = glmfit(X_m2, Y_valid, 'binomial', 'link', 'logit');
    
    % Predicted probabilities
    eta_m2 = [ones(length(Y_valid), 1), X_m2] * B_m2;
    prob_m2 = 1 ./ (1 + exp(-eta_m2));
    
    % AUC
    auc_m2 = computeAUC(Y_valid, prob_m2);
    
    % Log-likelihood and AIC
    ll_m2 = sum(Y_valid .* log(prob_m2 + 1e-10) + (1 - Y_valid) .* log(1 - prob_m2 + 1e-10));
    k_m2 = length(B_m2);
    aic_m2 = -2 * ll_m2 + 2 * k_m2;
    
    % AUC SE
    auc_se_m2 = sqrt((auc_m2 * (1 - auc_m2) + (nPos - 1) * (auc_m2 / (2 - auc_m2) - auc_m2^2) + ...
                      (nNeg - 1) * (2 * auc_m2^2 / (1 + auc_m2) - auc_m2^2)) / (nPos * nNeg));
    
    fprintf('  AUC = %.3f (SE = %.3f)\n', auc_m2, auc_se_m2);
    fprintf('  AIC = %.1f (Δ from Model 1: %.1f)\n', aic_m2, aic_m2 - aic_m1);
    fprintf('  Volume beta = %.3f, p = %.4f\n', B_m2(2), stats_m2.p(2));
    fprintf('  Total ICN Damage beta = %.3f, p = %.4f\n', B_m2(3), stats_m2.p(3));
    
    Level7b_Results.Model{2} = 'Volume + Total ICN Damage';
    Level7b_Results.N(2) = length(Y_valid);
    Level7b_Results.N_Aphasia(2) = sum(Y_valid == 1);
    Level7b_Results.N_Control(2) = sum(Y_valid == 0);
    Level7b_Results.N_Predictors(2) = k_m2 - 1;
    Level7b_Results.AUC(2) = auc_m2;
    Level7b_Results.AUC_SE(2) = auc_se_m2;
    Level7b_Results.LogLikelihood(2) = ll_m2;
    Level7b_Results.AIC(2) = aic_m2;
    
    % Interpretation
    auc_diff_m2 = auc_m2 - auc_m1;
    if auc_diff_m2 > 0.02 && stats_m2.p(3) < 0.05
        Level7b_Results.Interpretation{2} = 'ICN coverage adds beyond volume';
    elseif auc_diff_m2 < 0.01
        Level7b_Results.Interpretation{2} = 'Minimal improvement over volume alone';
    else
        Level7b_Results.Interpretation{2} = 'Modest improvement';
    end
    
catch ME
    fprintf('  ERROR: %s\n', ME.message);
    Level7b_Results.Model{2} = 'Volume + Total ICN Damage';
    Level7b_Results.AUC(2) = NaN;
end
fprintf('\n');

% -------------------------------------------------------------------------
% MODEL 3: Volume + ICN-Specific Damage (LASSO or Ridge if LASSO unavailable)
% -------------------------------------------------------------------------
fprintf('--- Model 3: Volume + ICN-Specific Damage ---\n');

% Check if lassoglm is available
hasLASSO = exist('lassoglm', 'file') == 2;

if hasLASSO
    fprintf('Using LASSO regularization (lassoglm)\n');
    try
        % Combine predictors: [volume, ICN1, ICN2, ..., ICN18, age, sex, days]
        X_m3 = [lesionVol_z, icnDamage_z, age_z, sex_valid, days_z];
        
        % Run LASSO with cross-validation
        [B_lasso, FitInfo] = lassoglm(X_m3, Y_valid, 'binomial', ...
            'CV', 5, 'Alpha', 1, 'NumLambda', 50);
        
        % Use lambda with minimum deviance (or 1SE rule)
        idxLambdaMinDeviance = FitInfo.IndexMinDeviance;
        B_m3_selected = B_lasso(:, idxLambdaMinDeviance);
        intercept_m3 = FitInfo.Intercept(idxLambdaMinDeviance);
        
        % Count non-zero ICN coefficients (indices 2 to nICN+1)
        icnCoefs = B_m3_selected(2:nICN+1);
        nSelectedICN = sum(abs(icnCoefs) > 1e-6);
        
        % Predicted probabilities
        eta_m3 = intercept_m3 + X_m3 * B_m3_selected;
        prob_m3 = 1 ./ (1 + exp(-eta_m3));
        
        % AUC
        auc_m3 = computeAUC(Y_valid, prob_m3);
        
        % Log-likelihood and AIC (approximate - use effective df)
        ll_m3 = sum(Y_valid .* log(prob_m3 + 1e-10) + (1 - Y_valid) .* log(1 - prob_m3 + 1e-10));
        k_m3 = sum(abs(B_m3_selected) > 1e-6) + 1;  % Non-zero coefs + intercept
        aic_m3 = -2 * ll_m3 + 2 * k_m3;
        
        % AUC SE
        auc_se_m3 = sqrt((auc_m3 * (1 - auc_m3) + (nPos - 1) * (auc_m3 / (2 - auc_m3) - auc_m3^2) + ...
                          (nNeg - 1) * (2 * auc_m3^2 / (1 + auc_m3) - auc_m3^2)) / (nPos * nNeg));
        
        fprintf('  AUC = %.3f (SE = %.3f)\n', auc_m3, auc_se_m3);
        fprintf('  AIC = %.1f (Δ from Model 1: %.1f)\n', aic_m3, aic_m3 - aic_m1);
        fprintf('  Non-zero ICN predictors: %d / %d\n', nSelectedICN, nICN);
        
        % List selected ICNs
        if nSelectedICN > 0
            fprintf('  Selected ICNs:\n');
            for n = 1:nICN
                if abs(icnCoefs(n)) > 1e-6
                    fprintf('    ICN%02d: beta = %.3f (%s)\n', n, icnCoefs(n), trueICN_labels{n}(1:min(30,end)));
                end
            end
        end
        
        Level7b_Results.Model{3} = 'Volume + ICN-Specific (LASSO)';
        Level7b_Results.N(3) = length(Y_valid);
        Level7b_Results.N_Aphasia(3) = sum(Y_valid == 1);
        Level7b_Results.N_Control(3) = sum(Y_valid == 0);
        Level7b_Results.N_Predictors(3) = k_m3 - 1;
        Level7b_Results.AUC(3) = auc_m3;
        Level7b_Results.AUC_SE(3) = auc_se_m3;
        Level7b_Results.LogLikelihood(3) = ll_m3;
        Level7b_Results.AIC(3) = aic_m3;
        
        % Store LASSO details for later
        Level7b_LASSO = struct();
        Level7b_LASSO.Coefficients = B_m3_selected;
        Level7b_LASSO.ICN_Coefficients = icnCoefs;
        Level7b_LASSO.SelectedICNs = find(abs(icnCoefs) > 1e-6);
        Level7b_LASSO.Lambda = FitInfo.Lambda(idxLambdaMinDeviance);
        
    catch ME
        fprintf('  LASSO ERROR: %s\n', ME.message);
        fprintf('  Falling back to simple logistic regression with all ICNs\n');
        hasLASSO = false;
    end
end

if ~hasLASSO
    % Fallback: Simple logistic with all ICNs (no regularization)
    fprintf('Using unregularized logistic regression (LASSO unavailable)\n');
    try
        X_m3 = [lesionVol_z, icnDamage_z, age_z, sex_valid, days_z];
        [B_m3, ~, stats_m3] = glmfit(X_m3, Y_valid, 'binomial', 'link', 'logit');
        
        eta_m3 = [ones(length(Y_valid), 1), X_m3] * B_m3;
        prob_m3 = 1 ./ (1 + exp(-eta_m3));
        
        auc_m3 = computeAUC(Y_valid, prob_m3);
        
        ll_m3 = sum(Y_valid .* log(prob_m3 + 1e-10) + (1 - Y_valid) .* log(1 - prob_m3 + 1e-10));
        k_m3 = length(B_m3);
        aic_m3 = -2 * ll_m3 + 2 * k_m3;
        
        auc_se_m3 = sqrt((auc_m3 * (1 - auc_m3) + (nPos - 1) * (auc_m3 / (2 - auc_m3) - auc_m3^2) + ...
                          (nNeg - 1) * (2 * auc_m3^2 / (1 + auc_m3) - auc_m3^2)) / (nPos * nNeg));
        
        fprintf('  AUC = %.3f (SE = %.3f)\n', auc_m3, auc_se_m3);
        fprintf('  AIC = %.1f (Δ from Model 1: %.1f)\n', aic_m3, aic_m3 - aic_m1);
        fprintf('  WARNING: Unregularized model may be overfit with %d predictors\n', k_m3 - 1);
        
        Level7b_Results.Model{3} = 'Volume + ICN-Specific (Unregularized)';
        Level7b_Results.N(3) = length(Y_valid);
        Level7b_Results.N_Aphasia(3) = sum(Y_valid == 1);
        Level7b_Results.N_Control(3) = sum(Y_valid == 0);
        Level7b_Results.N_Predictors(3) = k_m3 - 1;
        Level7b_Results.AUC(3) = auc_m3;
        Level7b_Results.AUC_SE(3) = auc_se_m3;
        Level7b_Results.LogLikelihood(3) = ll_m3;
        Level7b_Results.AIC(3) = aic_m3;
        
        Level7b_LASSO = struct();
        Level7b_LASSO.Coefficients = B_m3;
        Level7b_LASSO.ICN_Coefficients = B_m3(3:nICN+2);  % After intercept and volume
        Level7b_LASSO.SelectedICNs = 1:nICN;  % All included
        Level7b_LASSO.Lambda = NaN;
        
    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        Level7b_Results.Model{3} = 'Volume + ICN-Specific (FAILED)';
        Level7b_Results.AUC(3) = NaN;
        Level7b_LASSO = struct();
    end
end
fprintf('\n');

% -------------------------------------------------------------------------
% Level 7b Summary and Interpretation
% -------------------------------------------------------------------------
fprintf('--- Level 7b Summary: Model Comparison ---\n');
fprintf('%-35s %8s %8s %10s\n', 'Model', 'AUC', 'SE', 'AIC');
fprintf('%s\n', repmat('-', 1, 65));
for i = 1:height(Level7b_Results)
    if ~isnan(Level7b_Results.AUC(i))
        fprintf('%-35s %8.3f %8.3f %10.1f\n', ...
            Level7b_Results.Model{i}, Level7b_Results.AUC(i), ...
            Level7b_Results.AUC_SE(i), Level7b_Results.AIC(i));
    end
end
fprintf('\n');

% Statistical comparison: AUC difference significance (approximate z-test)
if ~isnan(Level7b_Results.AUC(1)) && ~isnan(Level7b_Results.AUC(2))
    auc_diff_12 = Level7b_Results.AUC(2) - Level7b_Results.AUC(1);
    se_diff_12 = sqrt(Level7b_Results.AUC_SE(1)^2 + Level7b_Results.AUC_SE(2)^2);
    z_12 = auc_diff_12 / se_diff_12;
    p_12 = 2 * (1 - normcdf(abs(z_12)));
    fprintf('Model 2 vs Model 1: ΔAUC = %.3f, z = %.2f, p = %.4f\n', auc_diff_12, z_12, p_12);
end

if ~isnan(Level7b_Results.AUC(1)) && ~isnan(Level7b_Results.AUC(3))
    auc_diff_13 = Level7b_Results.AUC(3) - Level7b_Results.AUC(1);
    se_diff_13 = sqrt(Level7b_Results.AUC_SE(1)^2 + Level7b_Results.AUC_SE(3)^2);
    z_13 = auc_diff_13 / se_diff_13;
    p_13 = 2 * (1 - normcdf(abs(z_13)));
    fprintf('Model 3 vs Model 1: ΔAUC = %.3f, z = %.2f, p = %.4f\n', auc_diff_13, z_13, p_13);
end

% Overall interpretation
fprintf('\n--- Level 7b Interpretation ---\n');
if ~isnan(Level7b_Results.AUC(1)) && ~isnan(Level7b_Results.AUC(3))
    auc_improvement = Level7b_Results.AUC(3) - Level7b_Results.AUC(1);
    
    if auc_improvement < 0.02
        fprintf('CONCLUSION: MASS EFFECT dominates.\n');
        fprintf('  ICN-specific damage adds minimal prediction beyond total lesion volume.\n');
        fprintf('  Aphasia occurrence is primarily determined by lesion SIZE, not LOCATION.\n');
        Level7b_Results.Interpretation{3} = 'MASS EFFECT: Volume dominates occurrence prediction';
    elseif auc_improvement >= 0.02 && auc_improvement < 0.05
        fprintf('CONCLUSION: MODEST MODULARITY effect.\n');
        fprintf('  ICN-specific damage adds some prediction, but volume remains primary.\n');
        Level7b_Results.Interpretation{3} = 'MODEST MODULARITY: Some ICN-specific contribution';
    else
        fprintf('CONCLUSION: MODULARITY effect detected.\n');
        fprintf('  ICN-specific damage substantially improves occurrence prediction.\n');
        fprintf('  Lesion LOCATION matters for aphasia development.\n');
        Level7b_Results.Interpretation{3} = 'MODULARITY: ICN location improves prediction';
    end
else
    fprintf('CONCLUSION: Unable to compare models (missing AUC values).\n');
end
fprintf('\n');

%% =========================================================================
%% LEVEL 8: LESION → ENGAGEMENT CORRELATION (DISCONNECTION) - ALL METRICS
%% =========================================================================
fprintf('=================================================================\n');
fprintf('LEVEL 8: Lesion Damage → Network Engagement (Disconnection)\n');
fprintf('         Analyzing ALL THREE METRICS\n');
fprintf('=================================================================\n');

% Use aphasia patients only (stroke controls may have minimal lesions)
aphasiaTable = mergedTable(isAphasia, :);
fprintf('Aphasia patients for Level 8: %d\n\n', height(aphasiaTable));

% Storage for Level 8 results - one table per metric
Level8_Results = struct();

for m = 1:nMetrics
    metric = metrics{m};
    
    fprintf('--- Metric: %s ---\n', metric);
    
    % Initialize results table for this metric
    L8_metric = table('Size', [nICN, 6], ...
        'VariableTypes', {'double', 'cell', 'double', 'double', 'double', 'logical'}, ...
        'VariableNames', {'ICN', 'ICN_Label', 'Spearman_rho', 'p_uncorrected', 'N', 'Significant_FDR'});
    
    for n = 1:nICN
        damageCol = sprintf('ICN%02d_Damage', n);
        engageCol = sprintf('ICN%02d_Eng_%s', n, metric);
        
        L8_metric.ICN(n) = n;
        L8_metric.ICN_Label{n} = trueICN_labels{n};
        
        if ~ismember(damageCol, aphasiaTable.Properties.VariableNames) || ...
           ~ismember(engageCol, aphasiaTable.Properties.VariableNames)
            L8_metric.Spearman_rho(n) = NaN;
            L8_metric.p_uncorrected(n) = NaN;
            L8_metric.N(n) = 0;
            continue;
        end
        
        damage = aphasiaTable.(damageCol);
        engagement = aphasiaTable.(engageCol);
        
        % Valid observations
        validIdx = ~isnan(damage) & ~isnan(engagement);
        L8_metric.N(n) = sum(validIdx);
        
        if sum(validIdx) < 10 || var(damage(validIdx)) < 1e-10
            L8_metric.Spearman_rho(n) = NaN;
            L8_metric.p_uncorrected(n) = NaN;
            continue;
        end
        
        % Spearman correlation: Damage → Engagement
        [rho, pval] = corr(damage(validIdx), engagement(validIdx), 'Type', 'Spearman');
        
        L8_metric.Spearman_rho(n) = rho;
        L8_metric.p_uncorrected(n) = pval;
    end
    
    % FDR correction
    [L8_metric.Significant_FDR, ~, L8_metric.p_FDR] = fdr_bh(L8_metric.p_uncorrected, 0.05);
    
    % Sort by rho
    L8_metric = sortrows(L8_metric, 'Spearman_rho');
    
    % Store in struct
    Level8_Results.(metric) = L8_metric;
    
    % Display summary
    sigL8_m = L8_metric(L8_metric.Significant_FDR == true, :);
    fprintf('Significant (FDR q<0.05): %d / %d | Mean rho: %.3f\n', ...
        height(sigL8_m), nICN, mean(L8_metric.Spearman_rho, 'omitnan'));
    
    if height(sigL8_m) > 0
        for i = 1:height(sigL8_m)
            fprintf('  ICN%02d: rho=%.3f, p_FDR=%.4f\n', ...
                sigL8_m.ICN(i), sigL8_m.Spearman_rho(i), sigL8_m.p_FDR(i));
        end
    end
    fprintf('\n');
end

% Display comparison across metrics
fprintf('--- METRIC COMPARISON (Level 8) ---\n');
fprintf('%-10s %12s %12s %12s\n', 'Metric', 'Mean rho', 'N sig', 'Interpretation');
fprintf('%s\n', repmat('-', 1, 55));
for m = 1:nMetrics
    metric = metrics{m};
    L8_m = Level8_Results.(metric);
    meanRho = mean(L8_m.Spearman_rho, 'omitnan');
    nSig = sum(L8_m.Significant_FDR);
    
    if strcmp(metric, primaryMetric)
        interp = 'PRIMARY';
    else
        interp = '';
    end
    
    fprintf('%-10s %12.3f %12d %12s\n', metric, meanRho, nSig, interp);
end
fprintf('\n');

%% LEVEL 8b: SPECIFICITY TEST (PRIMARY METRIC ONLY)
fprintf('-----------------------------------------------------------------\n');
fprintf('Level 8b: Specificity Test (Within vs Between Network)\n');
fprintf('         Using PRIMARY metric: %s\n', primaryMetric);
fprintf('-----------------------------------------------------------------\n');

Level8b_Results = table('Size', [nICN, 5], ...
    'VariableTypes', {'double', 'cell', 'double', 'double', 'double'}, ...
    'VariableNames', {'ICN', 'ICN_Label', 'Rho_Within', 'Rho_Between_Mean', 'Specificity'});

for n = 1:nICN
    damageCol = sprintf('ICN%02d_Damage', n);
    
    Level8b_Results.ICN(n) = n;
    Level8b_Results.ICN_Label{n} = trueICN_labels{n};
    
    if ~ismember(damageCol, aphasiaTable.Properties.VariableNames)
        Level8b_Results.Rho_Within(n) = NaN;
        Level8b_Results.Rho_Between_Mean(n) = NaN;
        Level8b_Results.Specificity(n) = NaN;
        continue;
    end
    
    damage = aphasiaTable.(damageCol);
    
    % Within-network correlation (using primary metric)
    engageCol_within = sprintf('ICN%02d_Eng_%s', n, primaryMetric);
    if ismember(engageCol_within, aphasiaTable.Properties.VariableNames)
        engage_within = aphasiaTable.(engageCol_within);
        validIdx = ~isnan(damage) & ~isnan(engage_within);
        
        if sum(validIdx) >= 10 && var(damage(validIdx)) > 1e-10
            [rho_within, ~] = corr(damage(validIdx), engage_within(validIdx), 'Type', 'Spearman');
        else
            rho_within = NaN;
        end
    else
        rho_within = NaN;
    end
    
    Level8b_Results.Rho_Within(n) = rho_within;
    
    % Between-network correlations
    rho_between = NaN(nICN-1, 1);
    betweenIdx = 0;
    for other = 1:nICN
        if other == n
            continue;
        end
        betweenIdx = betweenIdx + 1;
        engageCol_between = sprintf('ICN%02d_Eng_%s', other, primaryMetric);
        
        if ismember(engageCol_between, aphasiaTable.Properties.VariableNames)
            engage_between = aphasiaTable.(engageCol_between);
            validIdx_b = ~isnan(damage) & ~isnan(engage_between);
            
            if sum(validIdx_b) >= 10 && var(damage(validIdx_b)) > 1e-10
                rho_between(betweenIdx) = corr(damage(validIdx_b), engage_between(validIdx_b), 'Type', 'Spearman');
            end
        end
    end
    
    mean_rho_between = mean(rho_between, 'omitnan');
    Level8b_Results.Rho_Between_Mean(n) = mean_rho_between;
    
    % Specificity: within should be MORE NEGATIVE than between
    Level8b_Results.Specificity(n) = rho_within - mean_rho_between;
end

% Sort by specificity
Level8b_Results = sortrows(Level8b_Results, 'Specificity');

fprintf('Specificity = Within_rho - Mean_Between_rho\n');
fprintf('Negative specificity = damage more strongly reduces same-network engagement\n\n');

fprintf('%-6s %-30s %10s %10s %10s\n', 'ICN', 'Label', 'Within', 'Between', 'Specific');
fprintf('%s\n', repmat('-', 1, 70));
for i = 1:height(Level8b_Results)
    row = Level8b_Results(i, :);
    if ~isnan(row.Rho_Within)
        fprintf('ICN%02d  %-30s %10.3f %10.3f %10.3f\n', ...
            row.ICN, row.ICN_Label{1}(1:min(30,end)), ...
            row.Rho_Within, row.Rho_Between_Mean, row.Specificity);
    end
end
fprintf('\n');

%% =========================================================================
%% LEVEL 9: LESION → SEVERITY (REGRESSION)
%% =========================================================================
fprintf('=================================================================\n');
fprintf('LEVEL 9: Lesion Damage → Severity (WAB-AQ)\n');
fprintf('=================================================================\n');

Y_L9 = aphasiaTable.WAB_AQ;

% Covariates
age_aph = aphasiaTable.Age_At_Stroke;
sex_aph = aphasiaTable.SexNumeric;
days_aph = aphasiaTable.Days_Post_Stroke;
lesionVol_aph = aphasiaTable.LesionVolume_mm3;

% Storage for Level 9 results
Level9_Results = table('Size', [nICN, 9], ...
    'VariableTypes', {'double', 'cell', 'double', 'double', 'double', 'double', 'double', 'double', 'logical'}, ...
    'VariableNames', {'ICN', 'ICN_Label', 'Beta', 'SE', 't_stat', 'p_uncorrected', 'R2', 'N', 'Significant_FDR'});

for n = 1:nICN
    damageCol = sprintf('ICN%02d_Damage', n);
    
    Level9_Results.ICN(n) = n;
    Level9_Results.ICN_Label{n} = trueICN_labels{n};
    
    if ~ismember(damageCol, aphasiaTable.Properties.VariableNames)
        Level9_Results.Beta(n) = NaN;
        Level9_Results.SE(n) = NaN;
        Level9_Results.t_stat(n) = NaN;
        Level9_Results.p_uncorrected(n) = NaN;
        Level9_Results.R2(n) = NaN;
        Level9_Results.N(n) = 0;
        continue;
    end
    
    damage = aphasiaTable.(damageCol);
    
    % Valid observations
    validIdx = ~isnan(damage) & ~isnan(Y_L9) & ~isnan(age_aph) & ...
               ~isnan(sex_aph) & ~isnan(days_aph) & ~isnan(lesionVol_aph);
    
    Level9_Results.N(n) = sum(validIdx);
    
    if sum(validIdx) < 20 || var(damage(validIdx)) < 1e-10
        Level9_Results.Beta(n) = NaN;
        Level9_Results.SE(n) = NaN;
        Level9_Results.t_stat(n) = NaN;
        Level9_Results.p_uncorrected(n) = NaN;
        Level9_Results.R2(n) = NaN;
        continue;
    end
    
    % Build regression table
    tbl = table(Y_L9(validIdx), damage(validIdx), lesionVol_aph(validIdx), ...
                age_aph(validIdx), sex_aph(validIdx), days_aph(validIdx), ...
                'VariableNames', {'WAB', 'Damage', 'LesionVol', 'Age', 'Sex', 'Days'});
    
    try
        mdl = fitlm(tbl, 'WAB ~ Damage + LesionVol + Age + Sex + Days');
        
        damageRow = find(strcmp(mdl.CoefficientNames, 'Damage'));
        
        Level9_Results.Beta(n) = mdl.Coefficients.Estimate(damageRow);
        Level9_Results.SE(n) = mdl.Coefficients.SE(damageRow);
        Level9_Results.t_stat(n) = mdl.Coefficients.tStat(damageRow);
        Level9_Results.p_uncorrected(n) = mdl.Coefficients.pValue(damageRow);
        Level9_Results.R2(n) = mdl.Rsquared.Ordinary;
        
    catch ME
        fprintf('  ICN%02d: Regression failed - %s\n', n, ME.message);
        Level9_Results.Beta(n) = NaN;
        Level9_Results.SE(n) = NaN;
        Level9_Results.t_stat(n) = NaN;
        Level9_Results.p_uncorrected(n) = NaN;
        Level9_Results.R2(n) = NaN;
    end
end

% FDR correction
[Level9_Results.Significant_FDR, ~, Level9_Results.p_FDR] = fdr_bh(Level9_Results.p_uncorrected, 0.05);

% Sort by p-value
Level9_Results = sortrows(Level9_Results, 'p_uncorrected');

% Display results
fprintf('Level 9 Results: ICN Damage → WAB-AQ\n');
fprintf('Controlling for: Lesion Volume, Age, Sex, Days Post-Stroke\n\n');

sigL9 = Level9_Results(Level9_Results.Significant_FDR == true, :);
fprintf('Significant predictors (FDR q<0.05): %d / %d\n\n', height(sigL9), nICN);

fprintf('%-6s %-30s %10s %8s %8s %8s\n', 'ICN', 'Label', 'Beta', 't', 'p_unc', 'p_FDR');
fprintf('%s\n', repmat('-', 1, 80));
for i = 1:min(10, height(Level9_Results))
    row = Level9_Results(i, :);
    if ~isnan(row.Beta)
        if row.Significant_FDR
            sigMark = '*';
        else
            sigMark = '';
        end
        fprintf('ICN%02d  %-30s %10.2f %8.2f %8.4f %8.4f %s\n', ...
            row.ICN, row.ICN_Label{1}(1:min(30,end)), ...
            row.Beta, row.t_stat, row.p_uncorrected, row.p_FDR, sigMark);
    end
end
fprintf('\n');

%% =========================================================================
%% LEVEL 10: MEDIATION ANALYSIS - ALL METRICS
%% =========================================================================
fprintf('=================================================================\n');
fprintf('LEVEL 10: Mediation (Lesion → Engagement → Severity)\n');
fprintf('         Analyzing ALL THREE METRICS\n');
fprintf('=================================================================\n');

fprintf('Testing if engagement mediates the lesion → severity relationship\n\n');

% Storage for Level 10 results - one table per metric
Level10_Results = struct();

for m = 1:nMetrics
    metric = metrics{m};
    
    fprintf('--- Metric: %s ---\n', metric);
    
    % Initialize results table for this metric
    L10_metric = table('Size', [nICN, 10], ...
        'VariableTypes', {'double', 'cell', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
        'VariableNames', {'ICN', 'ICN_Label', 'Path_a', 'Path_b', 'Path_c', 'Path_c_prime', 'Indirect_ab', 'Mediation_Pct', 'N', 'Abs_Indirect'});
    
    for n = 1:nICN
        damageCol = sprintf('ICN%02d_Damage', n);
        engageCol = sprintf('ICN%02d_Eng_%s', n, metric);
        
        L10_metric.ICN(n) = n;
        L10_metric.ICN_Label{n} = trueICN_labels{n};
        
        if ~ismember(damageCol, aphasiaTable.Properties.VariableNames) || ...
           ~ismember(engageCol, aphasiaTable.Properties.VariableNames)
            L10_metric.Path_a(n) = NaN;
            L10_metric.Path_b(n) = NaN;
            L10_metric.Path_c(n) = NaN;
            L10_metric.Path_c_prime(n) = NaN;
            L10_metric.Indirect_ab(n) = NaN;
            L10_metric.Mediation_Pct(n) = NaN;
            L10_metric.N(n) = 0;
            L10_metric.Abs_Indirect(n) = NaN;
            continue;
        end
        
        X = aphasiaTable.(damageCol);        % Lesion damage (predictor)
        M = aphasiaTable.(engageCol);        % Engagement (mediator)
        Y = aphasiaTable.WAB_AQ;             % Severity (outcome)
        
        % Valid observations
        validIdx = ~isnan(X) & ~isnan(M) & ~isnan(Y);
        L10_metric.N(n) = sum(validIdx);
        
        if sum(validIdx) < 30 || var(X(validIdx)) < 1e-10
            L10_metric.Path_a(n) = NaN;
            L10_metric.Path_b(n) = NaN;
            L10_metric.Path_c(n) = NaN;
            L10_metric.Path_c_prime(n) = NaN;
            L10_metric.Indirect_ab(n) = NaN;
            L10_metric.Mediation_Pct(n) = NaN;
            L10_metric.Abs_Indirect(n) = NaN;
            continue;
        end
        
        X_v = X(validIdx);
        M_v = M(validIdx);
        Y_v = Y(validIdx);
        
        % Standardize for comparability
        X_z = (X_v - mean(X_v)) / std(X_v);
        M_z = (M_v - mean(M_v)) / std(M_v);
        Y_z = (Y_v - mean(Y_v)) / std(Y_v);
        
        try
            % Path a: X → M
            mdl_a = fitlm(X_z, M_z);
            a = mdl_a.Coefficients.Estimate(2);
            
            % Path c: X → Y (total effect)
            mdl_c = fitlm(X_z, Y_z);
            c = mdl_c.Coefficients.Estimate(2);
            
            % Paths b and c': X + M → Y
            mdl_bc = fitlm([X_z, M_z], Y_z);
            c_prime = mdl_bc.Coefficients.Estimate(2);  % Direct effect
            b = mdl_bc.Coefficients.Estimate(3);        % M → Y controlling for X
            
            % Indirect effect
            ab = a * b;
            
            % Mediation percentage (capped at reasonable values)
            if abs(c) > 0.001
                mediation_pct = (ab / c) * 100;
                % Cap at ±200% to handle suppression effects
                mediation_pct = max(min(mediation_pct, 200), -200);
            else
                mediation_pct = NaN;
            end
            
            L10_metric.Path_a(n) = a;
            L10_metric.Path_b(n) = b;
            L10_metric.Path_c(n) = c;
            L10_metric.Path_c_prime(n) = c_prime;
            L10_metric.Indirect_ab(n) = ab;
            L10_metric.Mediation_Pct(n) = mediation_pct;
            L10_metric.Abs_Indirect(n) = abs(ab);
            
        catch ME
            L10_metric.Path_a(n) = NaN;
            L10_metric.Path_b(n) = NaN;
            L10_metric.Path_c(n) = NaN;
            L10_metric.Path_c_prime(n) = NaN;
            L10_metric.Indirect_ab(n) = NaN;
            L10_metric.Mediation_Pct(n) = NaN;
            L10_metric.Abs_Indirect(n) = NaN;
        end
    end
    
    % Sort by absolute indirect effect
    L10_metric = sortrows(L10_metric, 'Abs_Indirect', 'descend');
    
    % Store in struct
    Level10_Results.(metric) = L10_metric;
    
    % Display top results
    fprintf('Top 3 mediation effects (by |ab|):\n');
    for i = 1:min(3, height(L10_metric))
        row = L10_metric(i, :);
        if ~isnan(row.Indirect_ab)
            fprintf('  ICN%02d: ab=%.3f, mediation=%.1f%%\n', ...
                row.ICN, row.Indirect_ab, row.Mediation_Pct);
        end
    end
    fprintf('\n');
end

% Display comparison across metrics
fprintf('--- METRIC COMPARISON (Level 10) ---\n');
fprintf('%-10s %15s %15s\n', 'Metric', 'Mean |ab|', 'Max Mediation%');
fprintf('%s\n', repmat('-', 1, 45));
for m = 1:nMetrics
    metric = metrics{m};
    L10_m = Level10_Results.(metric);
    meanAb = mean(L10_m.Abs_Indirect, 'omitnan');
    maxMed = max(abs(L10_m.Mediation_Pct));
    fprintf('%-10s %15.4f %15.1f\n', metric, meanAb, maxMed);
end

%% LEVEL 10b: BOOTSTRAP MEDIATION FOR KEY ICNs (NEW IN V3)
fprintf('=================================================================\n');
fprintf('LEVEL 10b: Bootstrap Mediation with BCa CIs (NEW IN V3)\n');
fprintf('=================================================================\n');
fprintf('Purpose: Formal significance test for indirect effects\n');
fprintf('Method: %d bootstrap iterations, %s confidence intervals\n', nBootstrap, bootCI_method);
fprintf('Focus: ICN17 (motor network - key severity predictor from Script 4)\n\n');

% Key ICNs to test (from Script 4 findings)
keyICNs_mediation = [17];  % ICN17 showed strongest severity prediction
% Can expand: keyICNs_mediation = [17, 9, 10];

% Storage for bootstrap results
Level10b_Bootstrap = table('Size', [length(keyICNs_mediation) * nMetrics, 14], ...
    'VariableTypes', {'double', 'cell', 'cell', 'double', 'double', 'double', ...
                      'double', 'double', 'double', 'double', 'double', ...
                      'logical', 'double', 'double'}, ...
    'VariableNames', {'ICN', 'ICN_Label', 'Metric', 'N', 'Path_a', 'Path_b', ...
                      'Indirect_ab', 'Indirect_SE', 'CI_Lower', 'CI_Upper', ...
                      'Mediation_Pct', 'Significant', 'Boot_p', 'N_Boot_Valid'});

bootIdx = 0;

for n = keyICNs_mediation
    for m = 1:nMetrics
        metric = metrics{m};
        bootIdx = bootIdx + 1;
        
        fprintf('Bootstrap: ICN%02d - %s\n', n, metric);
        
        damageCol = sprintf('ICN%02d_Damage', n);
        engageCol = sprintf('ICN%02d_Eng_%s', n, metric);
        
        Level10b_Bootstrap.ICN(bootIdx) = n;
        Level10b_Bootstrap.ICN_Label{bootIdx} = trueICN_labels{n};
        Level10b_Bootstrap.Metric{bootIdx} = metric;
        
        if ~ismember(damageCol, aphasiaTable.Properties.VariableNames) || ...
           ~ismember(engageCol, aphasiaTable.Properties.VariableNames)
            fprintf('  Columns not found - skipping\n');
            continue;
        end
        
        X = aphasiaTable.(damageCol);        % Lesion damage
        M = aphasiaTable.(engageCol);        % Engagement
        Y = aphasiaTable.WAB_AQ;             % Severity
        
        % Valid observations
        validIdx = ~isnan(X) & ~isnan(M) & ~isnan(Y);
        nValid = sum(validIdx);
        Level10b_Bootstrap.N(bootIdx) = nValid;
        
        if nValid < 30
            fprintf('  N=%d too small - skipping\n', nValid);
            continue;
        end
        
        X_v = X(validIdx);
        M_v = M(validIdx);
        Y_v = Y(validIdx);
        
        % Standardize
        X_z = (X_v - mean(X_v)) / std(X_v);
        M_z = (M_v - mean(M_v)) / std(M_v);
        Y_z = (Y_v - mean(Y_v)) / std(Y_v);
        
        % Point estimates (same as Level 10)
        mdl_a = fitlm(X_z, M_z);
        a_hat = mdl_a.Coefficients.Estimate(2);
        
        mdl_c = fitlm(X_z, Y_z);
        c_hat = mdl_c.Coefficients.Estimate(2);
        
        mdl_bc = fitlm([X_z, M_z], Y_z);
        b_hat = mdl_bc.Coefficients.Estimate(3);
        
        ab_hat = a_hat * b_hat;
        
        Level10b_Bootstrap.Path_a(bootIdx) = a_hat;
        Level10b_Bootstrap.Path_b(bootIdx) = b_hat;
        Level10b_Bootstrap.Indirect_ab(bootIdx) = ab_hat;
        
        % Bootstrap
        ab_boot = NaN(nBootstrap, 1);
        nBootValid = 0;
        
        rng(42);  % Reproducibility
        
        for b = 1:nBootstrap
            % Resample with replacement
            bootSample = randsample(nValid, nValid, true);
            
            X_b = X_z(bootSample);
            M_b = M_z(bootSample);
            Y_b = Y_z(bootSample);
            
            % Check variance in bootstrap sample
            if var(X_b) < 1e-10 || var(M_b) < 1e-10
                continue;
            end
            
            try
                % Path a
                mdl_a_b = fitlm(X_b, M_b);
                a_b = mdl_a_b.Coefficients.Estimate(2);
                
                % Path b (controlling for X)
                mdl_bc_b = fitlm([X_b, M_b], Y_b);
                b_b = mdl_bc_b.Coefficients.Estimate(3);
                
                ab_boot(b) = a_b * b_b;
                nBootValid = nBootValid + 1;
            catch
                % Skip failed iterations
                continue;
            end
        end
        
        Level10b_Bootstrap.N_Boot_Valid(bootIdx) = nBootValid;
        
        if nBootValid < nBootstrap * 0.9
            fprintf('  Warning: Only %d/%d valid bootstrap iterations\n', nBootValid, nBootstrap);
        end
        
        % Remove NaN
        ab_boot_valid = ab_boot(~isnan(ab_boot));
        
        if length(ab_boot_valid) < 100
            fprintf('  Insufficient valid bootstrap samples - skipping\n');
            continue;
        end
        
        % Bootstrap SE
        boot_se = std(ab_boot_valid);
        Level10b_Bootstrap.Indirect_SE(bootIdx) = boot_se;
        
        % BCa Confidence Intervals
        % Bias correction factor
        z0 = norminv(mean(ab_boot_valid < ab_hat));
        
        % Acceleration factor (jackknife)
        jackknife_ab = NaN(nValid, 1);
        for j = 1:nValid
            jk_idx = [1:j-1, j+1:nValid];
            X_jk = X_z(jk_idx);
            M_jk = M_z(jk_idx);
            Y_jk = Y_z(jk_idx);
            
            if var(X_jk) < 1e-10 || var(M_jk) < 1e-10
                continue;
            end
            
            try
                mdl_a_jk = fitlm(X_jk, M_jk);
                a_jk = mdl_a_jk.Coefficients.Estimate(2);
                mdl_bc_jk = fitlm([X_jk, M_jk], Y_jk);
                b_jk = mdl_bc_jk.Coefficients.Estimate(3);
                jackknife_ab(j) = a_jk * b_jk;
            catch
                continue;
            end
        end
        
        jackknife_valid = jackknife_ab(~isnan(jackknife_ab));
        if length(jackknife_valid) > 10
            theta_dot = mean(jackknife_valid);
            num = sum((theta_dot - jackknife_valid).^3);
            denom = 6 * (sum((theta_dot - jackknife_valid).^2))^1.5;
            if denom > 0
                acc = num / denom;
            else
                acc = 0;
            end
        else
            acc = 0;  % Fall back to percentile method
        end
        
        % BCa adjusted percentiles
        alpha_lo = bootAlpha / 2;
        alpha_hi = 1 - bootAlpha / 2;
        
        if ~isnan(z0) && ~isinf(z0)
            z_lo = norminv(alpha_lo);
            z_hi = norminv(alpha_hi);
            
            p_lo = normcdf(z0 + (z0 + z_lo) / (1 - acc * (z0 + z_lo)));
            p_hi = normcdf(z0 + (z0 + z_hi) / (1 - acc * (z0 + z_hi)));
            
            % Bound percentiles
            p_lo = max(0.001, min(0.999, p_lo));
            p_hi = max(0.001, min(0.999, p_hi));
        else
            % Fall back to percentile
            p_lo = alpha_lo;
            p_hi = alpha_hi;
        end
        
        ci_lower = quantile(ab_boot_valid, p_lo);
        ci_upper = quantile(ab_boot_valid, p_hi);
        
        Level10b_Bootstrap.CI_Lower(bootIdx) = ci_lower;
        Level10b_Bootstrap.CI_Upper(bootIdx) = ci_upper;
        
        % Significance: CI excludes zero
        isSignificant = (ci_lower > 0) || (ci_upper < 0);
        Level10b_Bootstrap.Significant(bootIdx) = isSignificant;
        
        % Bootstrap p-value (proportion of samples with opposite sign or zero)
        if ab_hat >= 0
            boot_p = 2 * mean(ab_boot_valid <= 0);
        else
            boot_p = 2 * mean(ab_boot_valid >= 0);
        end
        boot_p = min(boot_p, 1);
        Level10b_Bootstrap.Boot_p(bootIdx) = boot_p;
        
        % Mediation percentage (capped)
        if abs(c_hat) > 0.001
            med_pct = (ab_hat / c_hat) * 100;
            med_pct = max(min(med_pct, 200), -200);
        else
            med_pct = NaN;
        end
        Level10b_Bootstrap.Mediation_Pct(bootIdx) = med_pct;
        
        % Display
        sigStr = '';
        if isSignificant
            sigStr = '*** SIGNIFICANT';
        end
        fprintf('  ab = %.4f, 95%% BCa CI [%.4f, %.4f], p = %.4f %s\n', ...
            ab_hat, ci_lower, ci_upper, boot_p, sigStr);
    end
end

fprintf('\n--- Level 10b Summary ---\n');
fprintf('%-6s %-10s %10s %20s %10s %12s\n', 'ICN', 'Metric', 'ab', '95% BCa CI', 'p', 'Significant');
fprintf('%s\n', repmat('-', 1, 75));

for i = 1:height(Level10b_Bootstrap)
    row = Level10b_Bootstrap(i, :);
    if ~isnan(row.Indirect_ab)
        sigStr = '';
        if row.Significant
            sigStr = '***';
        end
        fprintf('ICN%02d  %-10s %10.4f [%8.4f, %8.4f] %10.4f %12s\n', ...
            row.ICN, row.Metric{1}, row.Indirect_ab, ...
            row.CI_Lower, row.CI_Upper, row.Boot_p, sigStr);
    end
end

% Interpretation
nSigBoot = sum(Level10b_Bootstrap.Significant);
fprintf('\nBootstrap mediation results: %d/%d significant indirect effects\n', ...
    nSigBoot, height(Level10b_Bootstrap));

if nSigBoot == 0
    fprintf('INTERPRETATION: No significant mediation effects.\n');
    fprintf('  Engagement does NOT mediate the lesion → severity relationship.\n');
    fprintf('  Consistent with mass effect interpretation from Levels 7b and 11.\n');
else
    fprintf('INTERPRETATION: %d significant mediation effects found.\n', nSigBoot);
    fprintf('  Some engagement variables mediate the lesion-severity link.\n');
end
fprintf('\n');

%% =========================================================================
%% LEVEL 11: WHITE MATTER DISCONNECTION (JHU ATLAS)
%% =========================================================================
fprintf('=================================================================\n');
fprintf('LEVEL 11: White Matter Disconnection Analysis\n');
fprintf('=================================================================\n');
fprintf('Using PARTIAL correlation controlling for total lesion volume\n\n');

if hasJHU
    fprintf('Using JHU White Matter Atlas\n');
    fprintf('JHU tracts: %d\n\n', nTracts);
   
    % Compute lesion-tract overlap for each aphasia patient
    tractOverlapMatrix = NaN(height(aphasiaTable), nTracts);
    
    fprintf('Computing tract damage for %d aphasia patients...\n', height(aphasiaTable));
    progressStep_L11 = max(1, round(height(aphasiaTable) / 10));
    
    for p = 1:height(aphasiaTable)
        pid = aphasiaTable.PatientID{p};
        resampledPath = fullfile(resampledDir, sprintf('%s_lesion_2mm.nii', pid));
        
        if ~exist(resampledPath, 'file')
            continue;
        end
        
        try
            V = spm_vol(resampledPath);
            lesion = spm_read_vols(V) > 0.5;
            
            % Check dimension match
            if ~isequal(size(lesion), size(jhuData(:,:,:,1)))
                if ~isequal(size(lesion), atlasDim)
                    continue;
                end
            end
            
            for t = 1:nTracts
                tractMask = jhuData(:,:,:,t) > 0;
                
                tractVoxels = sum(tractMask(:));
                overlapVoxels = sum(lesion(:) & tractMask(:));
                
                if tractVoxels > 0
                    tractOverlapMatrix(p, t) = overlapVoxels / tractVoxels;
                else
                    tractOverlapMatrix(p, t) = 0;
                end
            end
        catch ME
            % Skip silently
            continue;
        end
        
        % Progress
        if mod(p, progressStep_L11) == 0
            fprintf('  %d/%d (%.0f%%)\n', p, height(aphasiaTable), 100*p/height(aphasiaTable));
        end
    end
    
    fprintf('  %d/%d (100%%)\n', height(aphasiaTable), height(aphasiaTable));
    
    % Get lesion volume for aphasia patients (covariate for partial correlation)
    lesionVol_aph = aphasiaTable.LesionVolume_mm3;
    
    % Initialize table with BOTH bivariate and partial rho columns
    Level11_Results = table('Size', [nTracts, 9], ...
        'VariableTypes', {'double', 'cell', 'double', 'double', 'double', 'double', 'double', 'logical', 'double'}, ...
        'VariableNames', {'Tract', 'Tract_Label', 'Spearman_rho', 'Partial_rho', 'p_uncorrected', 'N', 'Mean_Damage', 'Significant_FDR', 'p_FDR'});
    
    validTractCount = 0;
    
    for t = 1:nTracts
        tractDamage = tractOverlapMatrix(:, t);
        wab = aphasiaTable.WAB_AQ;
        
        % Also require valid lesion volume for partial correlation
        validIdx = ~isnan(tractDamage) & ~isnan(wab) & ~isnan(lesionVol_aph) & tractDamage > 0;
        
        Level11_Results.Tract(t) = t;
        if t <= length(jhuLabels)
            Level11_Results.Tract_Label{t} = jhuLabels{t};
        else
            Level11_Results.Tract_Label{t} = sprintf('Tract_%d', t);
        end
        
        Level11_Results.N(t) = sum(validIdx);
        Level11_Results.Mean_Damage(t) = mean(tractDamage, 'omitnan');
        
        if sum(validIdx) < 20 || var(tractDamage(validIdx)) < 1e-10
            Level11_Results.Spearman_rho(t) = NaN;
            Level11_Results.Partial_rho(t) = NaN;
            Level11_Results.p_uncorrected(t) = NaN;
            continue;
        end
        
        % Bivariate correlation (kept for comparison)
        [rho_bivar, ~] = corr(tractDamage(validIdx), wab(validIdx), 'Type', 'Spearman');
        
        % Partial correlation controlling for total lesion volume
        try
            [rho_partial, pval_partial] = partialcorr(tractDamage(validIdx), wab(validIdx), ...
                lesionVol_aph(validIdx), 'Type', 'Spearman', 'Rows', 'complete');
        catch ME
            % Fallback if partialcorr fails (e.g., collinearity)
            warning('Tract %d: partialcorr failed - %s', t, ME.message);
            rho_partial = NaN;
            pval_partial = NaN;
        end
        
        Level11_Results.Spearman_rho(t) = rho_bivar;      % Bivariate (for comparison)
        Level11_Results.Partial_rho(t) = rho_partial;     % Lesion-volume-controlled
        Level11_Results.p_uncorrected(t) = pval_partial;  % Use partial p-value
        validTractCount = validTractCount + 1;
    end
    
    fprintf('\nTracts with sufficient data for analysis: %d / %d\n', validTractCount, nTracts);
    
    % Remove rows with NaN correlations for FDR
    validRows = ~isnan(Level11_Results.Partial_rho);
    
    if sum(validRows) > 0
        % FDR correction on valid rows only (using partial correlation p-values)
        pvals_valid = Level11_Results.p_uncorrected(validRows);
        [sig_valid, ~, pfdr_valid] = fdr_bh(pvals_valid, 0.05);
        
        Level11_Results.Significant_FDR(validRows) = sig_valid;
        Level11_Results.p_FDR(validRows) = pfdr_valid;
        
        % Sort by partial rho (most negative first - strongest severity predictors)
        Level11_Results = sortrows(Level11_Results, 'Partial_rho');
        
        fprintf('\nWhite matter tract damage → Severity correlations:\n');
        fprintf('Controlling for: Total lesion volume (partial correlation)\n\n');
        
        sigL11 = Level11_Results(Level11_Results.Significant_FDR == true, :);
        fprintf('Significant tracts (FDR q<0.05): %d / %d\n\n', height(sigL11), sum(validRows));
        
        % Display header with both rho values
        fprintf('%-40s %10s %10s %8s %8s\n', 'Tract', 'rho_bivar', 'rho_partial', 'p_unc', 'p_FDR');
        fprintf('%s\n', repmat('-', 1, 85));
        for i = 1:min(15, height(Level11_Results))
            row = Level11_Results(i, :);
            if ~isnan(row.Partial_rho)
                if row.Significant_FDR
                    sigMark = '*';
                else
                    sigMark = '';
                end
                labelStr = row.Tract_Label{1};
                if length(labelStr) > 40
                    labelStr = [labelStr(1:37) '...'];
                end
                fprintf('%-40s %10.3f %10.3f %8.4f %8.4f %s\n', ...
                    labelStr, row.Spearman_rho, row.Partial_rho, row.p_uncorrected, row.p_FDR, sigMark);
            end
        end
        
        % Report on attenuation from bivariate to partial
        fprintf('\n--- Attenuation Analysis ---\n');
        validBoth = ~isnan(Level11_Results.Spearman_rho) & ~isnan(Level11_Results.Partial_rho);
        if sum(validBoth) > 0
            mean_bivar = mean(abs(Level11_Results.Spearman_rho(validBoth)));
            mean_partial = mean(abs(Level11_Results.Partial_rho(validBoth)));
            attenuation_pct = (1 - mean_partial/mean_bivar) * 100;
            fprintf('Mean |rho| bivariate: %.3f\n', mean_bivar);
            fprintf('Mean |rho| partial (lesion vol controlled): %.3f\n', mean_partial);
            fprintf('Attenuation: %.1f%% (effect explained by total lesion volume)\n', attenuation_pct);
        end
    else
        fprintf('No tracts had sufficient data for correlation analysis.\n');
    end
else
    Level11_Results = table();
    fprintf('Level 11 skipped - JHU atlas not available\n');
end
fprintf('\n');

%% =========================================================================
%% SAVE ALL RESULTS
%% =========================================================================
fprintf('=================================================================\n');
fprintf('Saving results (v3 filenames)\n');
fprintf('=================================================================\n');

% Create lesionOverlapTable (subset of mergedTable with just lesion data)
damageColNames = arrayfun(@(n) sprintf('ICN%02d_Damage', n), 1:nICN, 'UniformOutput', false);
lesionOverlapTable = mergedTable(:, ['PatientID', 'LesionVolume_mm3', damageColNames]);

% Save lesion overlap data
save(fullfile(outputDir, 'ARC_04_v3_Lesion_ICN_Overlap.mat'), ...
    'lesionOverlapTable', 'mergedTable', 'overlapMatrix', 'patientIDs', 'lesionVolumeVec', ...
    'trueICN_idx', 'trueICN_labels', 'nICN', 'metrics', 'primaryMetric', 'IRi_threshold');
fprintf('Saved: ARC_04_v3_Lesion_ICN_Overlap.mat\n');

writetable(mergedTable, fullfile(outputDir, 'ARC_04_v3_Lesion_ICN_Overlap.csv'));
fprintf('Saved: ARC_04_v3_Lesion_ICN_Overlap.csv\n');

% Save resample log
writetable(resampleLog, fullfile(outputDir, 'ARC_04_v3_Resample_Log.csv'));
fprintf('Saved: ARC_04_v3_Resample_Log.csv\n');

% Save Level 7
writetable(Level7_Results, fullfile(outputDir, 'ARC_04_v3_Level7_LesionVsPresence.csv'));
fprintf('Saved: ARC_04_v3_Level7_LesionVsPresence.csv\n');

% Save Level 7b (NEW in v3)
writetable(Level7b_Results, fullfile(outputDir, 'ARC_04_v3_Level7b_OccurrencePrediction.csv'));
fprintf('Saved: ARC_04_v3_Level7b_OccurrencePrediction.csv\n');

% Save Level 7b LASSO details if available
if exist('Level7b_LASSO', 'var') && ~isempty(fieldnames(Level7b_LASSO))
    save(fullfile(outputDir, 'ARC_04_v3_Level7b_LASSO_Details.mat'), 'Level7b_LASSO');
    fprintf('Saved: ARC_04_v3_Level7b_LASSO_Details.mat\n');
end

% Save Level 8 - one file per metric
for m = 1:nMetrics
    metric = metrics{m};
    L8_file = fullfile(outputDir, sprintf('ARC_04_v3_Level8_DamageEngagement_%s.csv', metric));
    writetable(Level8_Results.(metric), L8_file);
    fprintf('Saved: ARC_04_v3_Level8_DamageEngagement_%s.csv\n', metric);
end

% Save Level 8b
writetable(Level8b_Results, fullfile(outputDir, 'ARC_04_v3_Level8b_Specificity.csv'));
fprintf('Saved: ARC_04_v3_Level8b_Specificity.csv\n');

% Save Level 9
writetable(Level9_Results, fullfile(outputDir, 'ARC_04_v3_Level9_LesionSeverity.csv'));
fprintf('Saved: ARC_04_v3_Level9_LesionSeverity.csv\n');

% Save Level 10 - one file per metric
for m = 1:nMetrics
    metric = metrics{m};
    L10_file = fullfile(outputDir, sprintf('ARC_04_v3_Level10_Mediation_%s.csv', metric));
    writetable(Level10_Results.(metric), L10_file);
    fprintf('Saved: ARC_04_v3_Level10_Mediation_%s.csv\n', metric);
end

% Save Level 10b Bootstrap (NEW in v3)
writetable(Level10b_Bootstrap, fullfile(outputDir, 'ARC_04_v3_Level10b_Bootstrap_Mediation.csv'));
fprintf('Saved: ARC_04_v3_Level10b_Bootstrap_Mediation.csv\n');

% Save Level 11 (includes Partial_rho)
if exist('Level11_Results', 'var') && height(Level11_Results) > 0
    writetable(Level11_Results, fullfile(outputDir, 'ARC_04_v3_Level11_Disconnection.csv'));
    fprintf('Saved: ARC_04_v3_Level11_Disconnection.csv\n');
end

% Save all results to .mat
save(fullfile(outputDir, 'ARC_04_v3_AllResults.mat'), ...
    'Level7_Results', 'Level7b_Results', 'Level8_Results', 'Level8b_Results', ...
    'Level9_Results', 'Level10_Results', 'Level10b_Bootstrap', 'Level11_Results', ...
    'mergedTable', 'trueICN_labels', 'nICN', 'trueICN_idx', ...
    'resampleLog', 'icnVoxelCounts', 'metrics', 'primaryMetric', 'nMetrics', ...
    'IRi_threshold', 'nBootstrap');
fprintf('Saved: ARC_04_v3_AllResults.mat\n\n');

%% GENERATE SUMMARY REPORT
fprintf('-----------------------------------------------------------------\n');
fprintf('Generating summary report\n');
fprintf('-----------------------------------------------------------------\n');

reportFile = fullfile(outputDir, 'ARC_04_v3_Summary.txt');
fid = fopen(reportFile, 'w');

fprintf(fid, '=================================================================\n');
fprintf(fid, 'ARC_04_v3 LESION-NETWORK ANALYSIS SUMMARY\n');
fprintf(fid, '=================================================================\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));

fprintf(fid, 'V3 CHANGES FROM V2:\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Input: ARC_03b_v3_Output (IRi QC-cleaned, threshold=%.1f)\n', IRi_threshold);
fprintf(fid, 'NEW: Level 7b - Aphasia Occurrence Prediction (3 model comparison)\n');
fprintf(fid, 'All output filenames updated to v3\n\n');

fprintf(fid, 'Metrics analyzed:\n');
for m = 1:nMetrics
    if strcmp(metrics{m}, primaryMetric)
        fprintf(fid, '  %s - %s ***PRIMARY***\n', metrics{m}, metricLabels{m});
    else
        fprintf(fid, '  %s - %s\n', metrics{m}, metricLabels{m});
    end
end
fprintf(fid, '\n');

fprintf(fid, 'DATA SUMMARY\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Lesion files found: %d\n', nLesionFiles);
fprintf(fid, 'Successfully parsed: %d\n', nParsed);
fprintf(fid, 'Matched with engagement data: %d\n', length(matchedPatients));
fprintf(fid, 'Successfully resampled: %d\n', nSuccess);
fprintf(fid, 'Final merged table: %d patients\n', height(mergedTable));
fprintf(fid, '  Aphasia patients: %d\n', sum(isAphasia));
fprintf(fid, '  Stroke controls: %d\n', sum(isControl));
fprintf(fid, 'ICNs analyzed: %d\n', nICN);
fprintf(fid, 'ICN indices used: %s\n\n', mat2str(trueICN_idx));

fprintf(fid, 'LESION CHARACTERISTICS\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Mean lesion volume: %.0f mm³ (SD: %.0f)\n', ...
    mean(mergedTable.LesionVolume_mm3, 'omitnan'), std(mergedTable.LesionVolume_mm3, 'omitnan'));
fprintf(fid, 'Lesion volume range: %.0f - %.0f mm³\n', ...
    min(mergedTable.LesionVolume_mm3), max(mergedTable.LesionVolume_mm3));
fprintf(fid, 'Mean ICN overlap: %.4f (range: %.4f - %.4f)\n\n', ...
    mean(overlapMatrix(:), 'omitnan'), min(overlapMatrix(:)), max(overlapMatrix(:)));

fprintf(fid, 'LEVEL 7: LESION → APHASIA PRESENCE (PER ICN)\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Analysis: Logistic regression (aphasia vs stroke control)\n');
fprintf(fid, 'Covariates: Lesion volume, Age, Sex, Days post-stroke\n');
fprintf(fid, 'Significant ICNs (FDR q<0.05): %d / %d\n', sum(Level7_Results.Significant_FDR), nICN);
if sum(Level7_Results.Significant_FDR) > 0
    sigL7_rep = Level7_Results(Level7_Results.Significant_FDR, :);
    for i = 1:height(sigL7_rep)
        fprintf(fid, '  ICN%02d: OR=%.2f, p_FDR=%.4f\n', sigL7_rep.ICN(i), sigL7_rep.OddsRatio(i), sigL7_rep.p_FDR(i));
    end
else
    fprintf(fid, '  No significant predictors\n');
end
fprintf(fid, '\n');

fprintf(fid, 'LEVEL 7b: APHASIA OCCURRENCE PREDICTION (NEW IN V3)\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Question: What predicts developing aphasia after left hemisphere stroke?\n');
fprintf(fid, 'Models compared:\n');
fprintf(fid, '  1. Volume Only: Lesion volume + covariates\n');
fprintf(fid, '  2. Volume + Total ICN Damage: Adds sum of all ICN overlaps\n');
fprintf(fid, '  3. Volume + ICN-Specific: LASSO across individual ICNs\n\n');

fprintf(fid, 'Model Comparison:\n');
fprintf(fid, '%-35s %8s %10s\n', 'Model', 'AUC', 'AIC');
fprintf(fid, '%s\n', repmat('-', 1, 55));
for i = 1:height(Level7b_Results)
    if ~isnan(Level7b_Results.AUC(i))
        fprintf(fid, '%-35s %8.3f %10.1f\n', ...
            Level7b_Results.Model{i}, Level7b_Results.AUC(i), Level7b_Results.AIC(i));
    end
end
fprintf(fid, '\n');

fprintf(fid, 'Interpretation: %s\n\n', Level7b_Results.Interpretation{3});

fprintf(fid, 'LEVEL 8: LESION → ENGAGEMENT (DISCONNECTION)\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Analysis: Spearman correlation (damage vs engagement, same ICN)\n');
fprintf(fid, 'Sample: Aphasia patients only (N=%d)\n', height(aphasiaTable));
fprintf(fid, 'Hypothesis: Negative correlation (damage reduces engagement)\n\n');

for m = 1:nMetrics
    metric = metrics{m};
    L8_m = Level8_Results.(metric);
    nSig = sum(L8_m.Significant_FDR);
    meanRho = mean(L8_m.Spearman_rho, 'omitnan');
    
    if strcmp(metric, primaryMetric)
        fprintf(fid, '%s (PRIMARY):\n', metric);
    else
        fprintf(fid, '%s:\n', metric);
    end
    fprintf(fid, '  Significant (FDR q<0.05): %d / %d\n', nSig, nICN);
    fprintf(fid, '  Mean within-network rho: %.3f\n', meanRho);
    
    if nSig > 0
        sigL8_rep = L8_m(L8_m.Significant_FDR, :);
        for i = 1:height(sigL8_rep)
            fprintf(fid, '    ICN%02d: rho=%.3f, p_FDR=%.4f\n', ...
                sigL8_rep.ICN(i), sigL8_rep.Spearman_rho(i), sigL8_rep.p_FDR(i));
        end
    end
    fprintf(fid, '\n');
end

fprintf(fid, 'LEVEL 8b: SPECIFICITY TEST\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Question: Is damage-engagement correlation specific to same network?\n');
fprintf(fid, 'Metric used: %s (PRIMARY)\n', primaryMetric);
fprintf(fid, 'Mean within-network rho: %.3f\n', mean(Level8b_Results.Rho_Within, 'omitnan'));
fprintf(fid, 'Mean between-network rho: %.3f\n', mean(Level8b_Results.Rho_Between_Mean, 'omitnan'));
fprintf(fid, 'Mean specificity: %.3f\n\n', mean(Level8b_Results.Specificity, 'omitnan'));

fprintf(fid, 'LEVEL 9: LESION → SEVERITY\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Analysis: Linear regression (ICN damage → WAB-AQ)\n');
fprintf(fid, 'Covariates: Lesion volume, Age, Sex, Days post-stroke\n');
fprintf(fid, 'Significant ICNs (FDR q<0.05): %d / %d\n', sum(Level9_Results.Significant_FDR), nICN);
if sum(Level9_Results.Significant_FDR) > 0
    sigL9_rep = Level9_Results(Level9_Results.Significant_FDR, :);
    for i = 1:height(sigL9_rep)
        fprintf(fid, '  ICN%02d: beta=%.2f, p_FDR=%.4f\n', sigL9_rep.ICN(i), sigL9_rep.Beta(i), sigL9_rep.p_FDR(i));
    end
end
fprintf(fid, '\n');

fprintf(fid, 'LEVEL 10: MEDIATION ANALYSIS\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Model: Lesion → Engagement → Severity\n');
fprintf(fid, 'Question: Does engagement mediate lesion-severity relationship?\n\n');

for m = 1:nMetrics
    metric = metrics{m};
    L10_m = Level10_Results.(metric);
    
    if strcmp(metric, primaryMetric)
        fprintf(fid, '%s (PRIMARY):\n', metric);
    else
        fprintf(fid, '%s:\n', metric);
    end
    
    fprintf(fid, '  Top 3 mediation effects (by |ab|):\n');
    for i = 1:min(3, height(L10_m))
        row = L10_m(i, :);
        if ~isnan(row.Indirect_ab)
            fprintf(fid, '    ICN%02d: ab=%.3f, mediation=%.1f%%\n', ...
                row.ICN, row.Indirect_ab, row.Mediation_Pct);
        end
    end
    fprintf(fid, '\n');
end

fprintf(fid, 'LEVEL 10b: BOOTSTRAP MEDIATION (NEW IN V3)\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Method: %d bootstrap iterations, BCa confidence intervals\n', nBootstrap);
fprintf(fid, 'Focus: ICN17 (motor network)\n\n');

fprintf(fid, 'Results:\n');
for i = 1:height(Level10b_Bootstrap)
    row = Level10b_Bootstrap(i, :);
    if ~isnan(row.Indirect_ab)
        sigStr = '';
        if row.Significant
            sigStr = ' ***SIGNIFICANT';
        end
        fprintf(fid, '  ICN%02d %s: ab=%.4f, 95%% CI [%.4f, %.4f], p=%.4f%s\n', ...
            row.ICN, row.Metric{1}, row.Indirect_ab, ...
            row.CI_Lower, row.CI_Upper, row.Boot_p, sigStr);
    end
end
fprintf(fid, '\nSignificant indirect effects: %d/%d\n', nSigBoot, height(Level10b_Bootstrap));
if nSigBoot == 0
    fprintf(fid, 'Interpretation: No mediation - engagement does not explain lesion-severity link.\n');
end
fprintf(fid, '\n');

if hasJHU && height(Level11_Results) > 0
    fprintf(fid, 'LEVEL 11: WHITE MATTER DISCONNECTION\n');
    fprintf(fid, '-----------------------------------------------------------------\n');
    fprintf(fid, 'Analysis: Tract damage → WAB-AQ PARTIAL correlation\n');
    fprintf(fid, 'Controlling for: Total lesion volume\n');
    fprintf(fid, 'Tracts analyzed: %d\n', sum(~isnan(Level11_Results.Partial_rho)));
    fprintf(fid, 'Significant tracts (FDR q<0.05): %d\n\n', sum(Level11_Results.Significant_FDR));
    
    if sum(Level11_Results.Significant_FDR) > 0
        sigL11_rep = Level11_Results(Level11_Results.Significant_FDR, :);
        fprintf(fid, 'Significant tracts:\n');
        for i = 1:height(sigL11_rep)
            fprintf(fid, '  %s:\n', sigL11_rep.Tract_Label{i}(1:min(40,end)));
            fprintf(fid, '    rho_bivariate = %.3f\n', sigL11_rep.Spearman_rho(i));
            fprintf(fid, '    rho_partial   = %.3f (controlling for lesion volume)\n', sigL11_rep.Partial_rho(i));
            fprintf(fid, '    p_FDR = %.4f\n', sigL11_rep.p_FDR(i));
        end
    else
        fprintf(fid, 'No significant tracts after lesion volume control.\n');
    end
    
    % Report attenuation
    validBoth = ~isnan(Level11_Results.Spearman_rho) & ~isnan(Level11_Results.Partial_rho);
    if sum(validBoth) > 0
        mean_bivar = mean(abs(Level11_Results.Spearman_rho(validBoth)));
        mean_partial = mean(abs(Level11_Results.Partial_rho(validBoth)));
        attenuation_pct = (1 - mean_partial/mean_bivar) * 100;
        fprintf(fid, '\nAttenuation analysis:\n');
        fprintf(fid, '  Mean |rho| bivariate: %.3f\n', mean_bivar);
        fprintf(fid, '  Mean |rho| partial: %.3f\n', mean_partial);
        fprintf(fid, '  Attenuation: %.1f%% (variance explained by total lesion volume)\n', attenuation_pct);
    end
    fprintf(fid, '\n');
end

fprintf(fid, 'KEY FINDINGS SUMMARY\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, '1. Lesion → Aphasia Presence (L7): ');
if sum(Level7_Results.Significant_FDR) > 0
    fprintf(fid, '%d ICNs predict aphasia diagnosis\n', sum(Level7_Results.Significant_FDR));
else
    fprintf(fid, 'No specific ICN damage predicts aphasia (after lesion volume control)\n');
end

fprintf(fid, '\n2. Aphasia Occurrence (L7b NEW): %s\n', Level7b_Results.Interpretation{3});
fprintf(fid, '   Best model AUC: %.3f\n', max(Level7b_Results.AUC));

fprintf(fid, '\n3. Disconnection Effect (L8):\n');
for m = 1:nMetrics
    metric = metrics{m};
    L8_m = Level8_Results.(metric);
    nSig = sum(L8_m.Significant_FDR);
    meanRho = mean(L8_m.Spearman_rho, 'omitnan');
    fprintf(fid, '   %s: %d significant, mean rho=%.3f\n', metric, nSig, meanRho);
end

fprintf(fid, '\n4. Lesion → Severity (L9): ');
if sum(Level9_Results.Significant_FDR) > 0
    fprintf(fid, '%d ICNs'' damage predicts severity\n', sum(Level9_Results.Significant_FDR));
else
    fprintf(fid, 'No specific ICN damage predicts severity (beyond total lesion volume)\n');
end

fprintf(fid, '\n5. Mediation (L10):\n');
for m = 1:nMetrics
    metric = metrics{m};
    L10_m = Level10_Results.(metric);
    topRow = L10_m(1, :);
    if ~isnan(topRow.Mediation_Pct)
        fprintf(fid, '   %s: Top=ICN%02d (%.1f%% mediation)\n', metric, topRow.ICN, topRow.Mediation_Pct);
    else
        fprintf(fid, '   %s: Limited mediation effects\n', metric);
    end
end

fprintf(fid, '\n5b. Bootstrap Mediation (L10b NEW):\n');
fprintf(fid, '   Significant indirect effects: %d/%d\n', nSigBoot, height(Level10b_Bootstrap));

fprintf(fid, '\n6. White Matter (L11): ');
if hasJHU && height(Level11_Results) > 0
    nSigL11 = sum(Level11_Results.Significant_FDR);
    if nSigL11 > 0
        fprintf(fid, '%d tracts predict severity after lesion volume control\n', nSigL11);
    else
        fprintf(fid, 'No tracts significant after controlling for lesion volume\n');
    end
else
    fprintf(fid, 'Not analyzed (JHU atlas unavailable)\n');
end

fprintf(fid, '\nOUTPUT FILES\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'ARC_04_v3_Lesion_ICN_Overlap.mat/csv - Merged lesion + engagement data\n');
fprintf(fid, 'ARC_04_v3_Resample_Log.csv - Lesion resampling status\n');
fprintf(fid, 'ARC_04_v3_Level7_LesionVsPresence.csv - Logistic regression results\n');
fprintf(fid, 'ARC_04_v3_Level7b_OccurrencePrediction.csv - NEW: 3-model comparison\n');
fprintf(fid, 'ARC_04_v3_Level7b_LASSO_Details.mat - NEW: LASSO coefficients\n');
for m = 1:nMetrics
    fprintf(fid, 'ARC_04_v3_Level8_DamageEngagement_%s.csv - Disconnection\n', metrics{m});
end
fprintf(fid, 'ARC_04_v3_Level8b_Specificity.csv - Within vs between network\n');
fprintf(fid, 'ARC_04_v3_Level9_LesionSeverity.csv - Damage → severity regression\n');
for m = 1:nMetrics
fprintf(fid, 'ARC_04_v3_Level10_Mediation_%s.csv - Mediation analysis\n', metrics{m});
end
if hasJHU
    fprintf(fid, 'ARC_04_v3_Level11_Disconnection.csv - White matter analysis\n');
end
fprintf(fid, 'ARC_04_v3_AllResults.mat - All results in one file\n');
fprintf(fid, 'ARC_04_v3_Summary.txt - This report\n');

fprintf(fid, '\n=================================================================\n');
fprintf(fid, 'END OF REPORT\n');
fprintf(fid, '=================================================================\n');

fclose(fid);
fprintf('Saved: ARC_04_v3_Summary.txt\n\n');

%% =========================================================================
%% COMPLETION
%% =========================================================================
fprintf('=================================================================\n');
fprintf('ARC_04_v3_Lesion_Network_Analysis.m - COMPLETE\n');
fprintf('=================================================================\n');
fprintf('End time: %s\n', datestr(now));
fprintf('\nV3 CHANGES FROM V2:\n');
fprintf('  - Input: ARC_03b_v3_Output (IRi QC threshold=%.1f)\n', IRi_threshold);
fprintf('  - NEW: Level 7b Aphasia Occurrence Prediction (3 models)\n');
fprintf('  - All output filenames: ARC_04_v3_*\n');
fprintf('\nMETRICS ANALYZED:\n');
for m = 1:nMetrics
    if strcmp(metrics{m}, primaryMetric)
        fprintf('  %s - %s ***PRIMARY***\n', metrics{m}, metricLabels{m});
    else
        fprintf('  %s - %s\n', metrics{m}, metricLabels{m});
    end
end
fprintf('\nHYPOTHESIS SUMMARY:\n');
fprintf('  Level 7 (Lesion → Presence): %d/%d significant\n', sum(Level7_Results.Significant_FDR), nICN);
fprintf('  Level 7b (Occurrence Prediction): %s\n', Level7b_Results.Interpretation{3});
fprintf('  Level 8 (Disconnection):\n');
for m = 1:nMetrics
    metric = metrics{m};
    L8_m = Level8_Results.(metric);
    fprintf('    %s: %d/%d significant\n', metric, sum(L8_m.Significant_FDR), nICN);
end
fprintf('  Level 9 (Lesion → Severity): %d/%d significant\n', sum(Level9_Results.Significant_FDR), nICN);
if hasJHU && height(Level11_Results) > 0
    fprintf('  Level 11 (White Matter): %d/%d significant\n', sum(Level11_Results.Significant_FDR), sum(~isnan(Level11_Results.Partial_rho)));
end
fprintf('\nOutputs saved to: %s\n', outputDir);
fprintf('\nOutput files:\n');
fprintf('  - ARC_04_v3_Lesion_ICN_Overlap.mat/csv\n');
fprintf('  - ARC_04_v3_Resample_Log.csv\n');
fprintf('  - ARC_04_v3_Level7_LesionVsPresence.csv\n');
fprintf('  - ARC_04_v3_Level7b_OccurrencePrediction.csv (NEW)\n');
fprintf('  - ARC_04_v3_Level7b_LASSO_Details.mat (NEW)\n');
for m = 1:nMetrics
    fprintf('  - ARC_04_v3_Level8_DamageEngagement_%s.csv\n', metrics{m});
end
fprintf('  - ARC_04_v3_Level8b_Specificity.csv\n');
fprintf('  - ARC_04_v3_Level9_LesionSeverity.csv\n');
for m = 1:nMetrics
    fprintf('  - ARC_04_v3_Level10_Mediation_%s.csv\n', metrics{m});
end
if hasJHU
    fprintf('  - ARC_04_v3_Level11_Disconnection.csv\n');
end
fprintf('  - ARC_04_v3_AllResults.mat\n');
fprintf('  - ARC_04_v3_Summary.txt\n');
fprintf('\nNext steps:\n');
fprintf('  1. Review ARC_04_v3_Summary.txt\n');
fprintf('  2. Run ARC_05_v3 for volume-stratified analysis\n');
fprintf('  3. Run ARC_06_v3 for Smith10 robustness check\n');
fprintf('  4. Run ARC_07 for final integrated report\n');
fprintf('=================================================================\n');
