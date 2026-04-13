%% ARC_03b_v3_Compile_Results.m
% =========================================================================
% SCRIPT 3b v3: Compile ICN Results (IRi QC Added)
% =========================================================================
% Purpose:
%   - Parse all individual xATL.mat files from Script 3a
%   - Extract CORRECT metrics (IRi, MANi, Vari) for 18 ICNs
%   - V3: Add IRi Quality Control - exclude subjects with IRi sum < 0.8
%   - Merge with clinical data (WAB-AQ, age, sex, aphasia type)
%   - Create master analysis tables (long and wide formats)
%
% V3 CHANGES FROM V2:
%   - Added IRi Quality Control section (excludes subjects with IRi sum < 0.8)
%   - Generates IRi QC report (ARC_03b_v3_IRi_QC_Report.csv)
%   - All output filenames updated to v3
%   - Documented IRi threshold as principled criterion (>20% outside ICNs)
%
% METRIC DEFINITIONS:
%   - IRi: ICNiRelativeInvolvement (% of whole-brain activation in ICN)
%   - MANi: NormalisedMeanICNiActivation (0-1 normalized)
%   - Vari: var(ZActiveVoxels) (spatial variance of Z-scores)
%
% Input:
%   - Individual xATL files from ARC_03_Output/Individual_ICN/
%   - ARC_inventory.mat from Script 1 (clinical data)
%   - ARC_02c_ConFile_Inventory.mat from Script 2c (file metadata)
%
% Output:
%   - ARC_03b_v3_IRi_QC_Report.csv (QC exclusion report)
%   - ARC_03b_v3_Master_Long.mat/csv (one row per patient-contrast)
%   - ARC_03b_v3_Master_Wide.mat/csv (one row per patient)
%   - ARC_03b_v3_Network_Labels.csv (ICN names reference)
%   - ARC_03b_v3_Summary.txt
% =========================================================================

%% SETUP
% =========================================================================
clear; clc;

% Load configuration
config = config_local();

fprintf('=================================================================\n');
fprintf('ARC_03b_v3_Compile_Results.m (IRi QC ADDED)\n');
fprintf('Compile ICN Results & Merge Clinical Data\n');
fprintf('=================================================================\n');
fprintf('Start time: %s\n\n', datestr(now));

fprintf('*** V3 CHANGES ***\n');
fprintf('IRi QC: Subjects with IRi sum < 0.8 will be EXCLUDED\n');
fprintf('Rationale: >20%% of activation outside defined ICNs indicates\n');
fprintf('           poor atlas coverage or data quality issues.\n');
fprintf('This is a PRINCIPLED criterion, not data-driven.\n\n');

%% V3 CHANGE: Define IRi QC threshold
IRi_threshold = 0.8;  % Principled: >20% outside ICNs is problematic

%% DEFINE PATHS
% =========================================================================
icnOutputDir = fullfile(config.rootDir, 'ARC_03_Output', 'Individual_ICN');
inventoryFile = fullfile(config.rootDir, 'Analysis_Output', 'ARC_inventory.mat');
conInventoryFile = fullfile(config.rootDir, 'ARC_02_Output', 'ARC_02c_ConFile_Inventory.mat');

%% V3 CHANGE: Output directory updated to v3
outputDir = fullfile(config.rootDir, 'ARC_03_Output', 'ARC_03b_v3_Output');

% Create output directory if needed
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('Created output directory: %s\n\n', outputDir);
else
    fprintf('Output directory: %s\n\n', outputDir);
end

%% GROUP LABELING REMINDER
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('GROUP LABELING REMINDER\n');
fprintf('-----------------------------------------------------------------\n');
fprintf('The "Unspecified" folder contains patients with wab_type = "None"\n');
fprintf('These are STROKE PATIENTS WITHOUT DIAGNOSED APHASIA (N=32)\n');
fprintf('They are NOT healthy controls - they had left hemisphere strokes\n');
fprintf('but did not develop aphasia (WAB-AQ typically >93)\n\n');

%% VERIFY INPUTS
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Verifying inputs\n');
fprintf('-----------------------------------------------------------------\n');

% Check xATL directory
if ~exist(icnOutputDir, 'dir')
    error('ICN output directory not found: %s\nRun Script 3a first.', icnOutputDir);
end
xatlFiles = dir(fullfile(icnOutputDir, '*.xATL.mat'));
fprintf('xATL files found: %d\n', length(xatlFiles));

% Check inventory files
if ~exist(inventoryFile, 'file')
    error('Inventory file not found: %s\nRun Script 1 first.', inventoryFile);
end
fprintf('Clinical inventory: ✓\n');

if ~exist(conInventoryFile, 'file')
    error('Con inventory file not found: %s\nRun Script 2c first.', conInventoryFile);
end
fprintf('Con file inventory: ✓\n\n');

%% LOAD CLINICAL DATA
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Loading clinical data\n');
fprintf('-----------------------------------------------------------------\n');

% Load Script 1 inventory (has clinical data per patient)
load(inventoryFile, 'inventory');
fprintf('Loaded inventory: %d rows\n', height(inventory));

% Create patient-level clinical lookup (first row per patient)
[uniquePatients, firstIdx] = unique(inventory.PatientID, 'stable');
clinicalLookup = inventory(firstIdx, :);
fprintf('Unique patients with clinical data: %d\n', height(clinicalLookup));

% Display clinical variables available
clinicalVars = {'PatientID', 'WAB_AQ', 'WAB_Type', 'Age_At_Stroke', 'Sex', 'Days_Post_Stroke'};
fprintf('Clinical variables: %s\n\n', strjoin(clinicalVars, ', '));

% Load Con inventory (has aphasia type from folder structure)
load(conInventoryFile, 'conInventory');
fprintf('Loaded con inventory: %d rows\n\n', height(conInventory));

%% INSPECT XATL STRUCTURE (FIRST FILE)
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Inspecting xATL structure\n');
fprintf('-----------------------------------------------------------------\n');

% Load first file to understand structure
sampleFile = fullfile(xatlFiles(1).folder, xatlFiles(1).name);
sampleData = load(sampleFile);

% Detect structure (handle different possible formats)
if isfield(sampleData, 'xATL')
    xATL = sampleData.xATL;
    structureType = 'direct';
elseif isfield(sampleData, 'ICN_Result')
    if isfield(sampleData.ICN_Result, 'xATL')
        xATL = sampleData.ICN_Result.xATL;
        if isfield(xATL, 'xATL')
            xATL = xATL.xATL;  % Double-nested
        end
    end
    structureType = 'wrapped';
else
    fields = fieldnames(sampleData);
    xATL = sampleData.(fields{1});
    structureType = 'unknown';
end

fprintf('Structure type detected: %s\n', structureType);

% Get network information
if isfield(xATL, 'Atl') && isfield(xATL, 'Labels')
    nNetworksTotal = length(xATL.Atl);
    networkLabels = xATL.Labels;
    fprintf('Total networks in atlas: %d\n', nNetworksTotal);
    
    % Display all labels to identify artifacts
    fprintf('\nNetwork labels:\n');
    for n = 1:nNetworksTotal
        fprintf('  %2d. %s\n', n, networkLabels{n});
    end
else
    error('Unexpected xATL structure. Check file: %s', sampleFile);
end

%% IDENTIFY TRUE ICNs (EXCLUDE ARTIFACTS)
% =========================================================================
fprintf('\n-----------------------------------------------------------------\n');
fprintf('Identifying true ICNs (excluding artifacts)\n');
fprintf('-----------------------------------------------------------------\n');

% Identify artifact/outside components by label
artifactIdx = [];
for n = 1:nNetworksTotal
    label = lower(networkLabels{n});
    if contains(label, 'artifact') || contains(label, 'outside') || ...
       contains(label, 'noise') || contains(label, 'edge')
        artifactIdx(end+1) = n; %#ok<SAGROW>
    end
end

% True ICN indices
trueICN_idx = setdiff(1:nNetworksTotal, artifactIdx);
nICN = length(trueICN_idx);
trueICN_labels = networkLabels(trueICN_idx);

fprintf('Artifact components excluded: %d (indices: %s)\n', length(artifactIdx), mat2str(artifactIdx));
fprintf('True ICNs for analysis: %d\n\n', nICN);

fprintf('ICNs to analyze:\n');
for n = 1:nICN
    fprintf('  ICN%02d: %s\n', n, trueICN_labels{n});
end
fprintf('\n');

%% VERIFY TARGET FIELDS EXIST
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Verifying target metric fields exist\n');
fprintf('-----------------------------------------------------------------\n');

sampleAtl = xATL.Atl(trueICN_idx(1));

% Check IRi
if isfield(sampleAtl, 'ICNiRelativeInvolvement')
    fprintf('✓ IRi field (ICNiRelativeInvolvement): Found\n');
else
    error('IRi field (ICNiRelativeInvolvement) not found in xATL structure');
end

% Check MANi
if isfield(sampleAtl, 'NormalisedMeanICNiActivation')
    fprintf('✓ MANi field (NormalisedMeanICNiActivation): Found\n');
else
    error('MANi field (NormalisedMeanICNiActivation) not found in xATL structure');
end

% Check ZActiveVoxels for Vari computation
if isfield(sampleAtl, 'ZActiveVoxels')
    fprintf('✓ Vari source (ZActiveVoxels): Found\n');
else
    error('ZActiveVoxels not found - cannot compute Vari');
end

fprintf('\nAll required fields verified. Proceeding with extraction.\n\n');

%% EXTRACT METRICS FROM ALL FILES
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Extracting CORRECTED metrics from %d files\n', length(xatlFiles));
fprintf('-----------------------------------------------------------------\n');

% Initialize storage
nFiles = length(xatlFiles);

% Metadata columns
PatientID = cell(nFiles, 1);
SessionID = cell(nFiles, 1);
ContrastNum = cell(nFiles, 1);
ContrastLabel = cell(nFiles, 1);

% CORRECTED Metrics: 3 metrics × nICN networks
IRi_matrix = NaN(nFiles, nICN);     % Relative Spatial Involvement (CORRECTED)
MANi_matrix = NaN(nFiles, nICN);    % Normalized Mean Activation (CORRECTED)
Vari_matrix = NaN(nFiles, nICN);    % Variance of Z-scores (NEW)

% Additional tracking
ExtractionSuccess = false(nFiles, 1);
ErrorMessage = cell(nFiles, 1);
NumActiveVoxels_matrix = NaN(nFiles, nICN);  % For QC

% Progress tracking
progressInterval = max(1, round(nFiles / 10));

for f = 1:nFiles
    % Progress update
    if mod(f, progressInterval) == 0
        fprintf('  Processing: %d/%d (%.0f%%)\n', f, nFiles, (f/nFiles)*100);
    end
    
    % Parse filename: [PatientID]_ses-[Session]_con_[Contrast].xATL.mat
    fileName = xatlFiles(f).name;
    tokens = regexp(fileName, '(\w+)_ses-(\d+)_con_(\d+)\.xATL\.mat', 'tokens');
    
    if ~isempty(tokens)
        PatientID{f} = tokens{1}{1};
        SessionID{f} = tokens{1}{2};
        ContrastNum{f} = tokens{1}{3};
        
        % Assign contrast label
        if strcmp(ContrastNum{f}, '0001')
            ContrastLabel{f} = 'TaskRest';
        elseif strcmp(ContrastNum{f}, '0002')
            ContrastLabel{f} = 'NamingAbstract';
        else
            ContrastLabel{f} = 'Other';
        end
    else
        PatientID{f} = 'PARSE_ERROR';
        SessionID{f} = 'PARSE_ERROR';
        ContrastNum{f} = 'PARSE_ERROR';
        ContrastLabel{f} = 'PARSE_ERROR';
        ErrorMessage{f} = 'Filename parsing failed';
        continue;
    end
    
    % Load xATL file
    try
        filePath = fullfile(xatlFiles(f).folder, fileName);
        data = load(filePath);
        
        % Extract xATL (handle structure variations)
        if isfield(data, 'xATL')
            xATL = data.xATL;
        elseif isfield(data, 'ICN_Result')
            xATL = data.ICN_Result.xATL;
            if isfield(xATL, 'xATL')
                xATL = xATL.xATL;
            end
        else
            fields = fieldnames(data);
            xATL = data.(fields{1});
        end
        
        % Extract CORRECTED metrics for each true ICN
        for n = 1:nICN
            atlIdx = trueICN_idx(n);  % Index in original atlas
            
            if atlIdx <= length(xATL.Atl)
                atlData = xATL.Atl(atlIdx);
                
                % IRi (Relative Spatial Involvement) - CORRECTED
                if isfield(atlData, 'ICNiRelativeInvolvement')
                    IRi_matrix(f, n) = atlData.ICNiRelativeInvolvement;
                end
                
                % MANi (Normalized Mean ICN Activation) - CORRECTED
                if isfield(atlData, 'NormalisedMeanICNiActivation')
                    MANi_matrix(f, n) = atlData.NormalisedMeanICNiActivation;
                end
                
                % Vari (Variance of Z-scores) - NEW, computed
                if isfield(atlData, 'ZActiveVoxels')
                    zVals = atlData.ZActiveVoxels;
                    if length(zVals) > 1
                        Vari_matrix(f, n) = var(zVals);
                    elseif length(zVals) == 1
                        Vari_matrix(f, n) = 0;  % Single voxel = no variance
                    else
                        Vari_matrix(f, n) = NaN;  % No active voxels
                    end
                end
                
                % Store voxel count for QC
                if isfield(atlData, 'NumActiveVoxels')
                    NumActiveVoxels_matrix(f, n) = atlData.NumActiveVoxels;
                end
            end
        end
        
        ExtractionSuccess(f) = true;
        ErrorMessage{f} = '';
        
    catch ME
        ExtractionSuccess(f) = false;
        ErrorMessage{f} = ME.message;
    end
end

fprintf('  Processing: %d/%d (100%%)\n\n', nFiles, nFiles);

% Report extraction success
nSuccess = sum(ExtractionSuccess);
nFailed = sum(~ExtractionSuccess);
fprintf('Extraction results:\n');
fprintf('  Successful: %d\n', nSuccess);
fprintf('  Failed: %d\n\n', nFailed);

if nFailed > 0
    fprintf('Failed files:\n');
    failedIdx = find(~ExtractionSuccess);
    for i = 1:min(10, length(failedIdx))
        idx = failedIdx(i);
        fprintf('  %s: %s\n', xatlFiles(idx).name, ErrorMessage{idx});
    end
    if length(failedIdx) > 10
        fprintf('  ... and %d more\n', length(failedIdx) - 10);
    end
    fprintf('\n');
end

%% V3 CHANGE: IRi QUALITY CONTROL
% =========================================================================
fprintf('=================================================================\n');
fprintf('V3: IRi QUALITY CONTROL\n');
fprintf('=================================================================\n');
fprintf('Criterion: IRi sum < %.1f indicates >%.0f%% of activation outside\n', ...
    IRi_threshold, (1 - IRi_threshold) * 100);
fprintf('           defined ICNs - flagged for EXCLUSION.\n');
fprintf('Rationale: Poor atlas coverage or data quality issues.\n');
fprintf('Status: PRINCIPLED criterion (not data-driven).\n\n');

% Calculate IRi sums
IRi_sums = sum(IRi_matrix, 2, 'omitnan');

% Flag low IRi subjects
lowIRi_flag = IRi_sums < IRi_threshold;

% Also flag extraction failures
lowIRi_flag(~ExtractionSuccess) = true;

fprintf('IRi Distribution (before exclusion):\n');
fprintf('  Mean: %.3f\n', mean(IRi_sums(ExtractionSuccess), 'omitnan'));
fprintf('  Median: %.3f\n', median(IRi_sums(ExtractionSuccess), 'omitnan'));
fprintf('  Range: %.3f - %.3f\n', min(IRi_sums(ExtractionSuccess)), max(IRi_sums(ExtractionSuccess)));
fprintf('  SD: %.3f\n\n', std(IRi_sums(ExtractionSuccess), 'omitnan'));

fprintf('IRi QC Results:\n');
fprintf('  Threshold: IRi sum < %.1f\n', IRi_threshold);
fprintf('  Observations flagged: %d/%d (%.1f%%)\n', ...
    sum(lowIRi_flag), nFiles, 100*sum(lowIRi_flag)/nFiles);

% Create QC report table
IRi_QC = table(PatientID, SessionID, ContrastNum, ContrastLabel, ...
    IRi_sums, ExtractionSuccess, lowIRi_flag, ErrorMessage, ...
    'VariableNames', {'PatientID', 'SessionID', 'ContrastNum', 'ContrastLabel', ...
                      'IRi_Sum', 'ExtractionSuccess', 'Excluded', 'ErrorMessage'});

% Sort by IRi sum to show worst cases first
IRi_QC = sortrows(IRi_QC, 'IRi_Sum');

% Save QC report
qcReportFile = fullfile(outputDir, 'ARC_03b_v3_IRi_QC_Report.csv');
writetable(IRi_QC, qcReportFile);
fprintf('\nSaved QC Report: %s\n', qcReportFile);

% List excluded subjects
excludedQC = IRi_QC(IRi_QC.Excluded, :);
if height(excludedQC) > 0
    fprintf('\nExcluded observations (%d):\n', height(excludedQC));
    fprintf('%-15s %-8s %-8s %-10s %-10s\n', 'PatientID', 'Session', 'Contrast', 'IRi_Sum', 'Reason');
    fprintf('%s\n', repmat('-', 1, 60));
    
    for i = 1:min(30, height(excludedQC))
        row = excludedQC(i, :);
        if row.ExtractionSuccess
            reason = 'Low IRi';
        else
            reason = 'Extract fail';
        end
        fprintf('%-15s %-8s %-8s %-10.3f %-10s\n', ...
            row.PatientID{1}, row.SessionID{1}, row.ContrastNum{1}, ...
            row.IRi_Sum, reason);
    end
    if height(excludedQC) > 30
        fprintf('... and %d more\n', height(excludedQC) - 30);
    end
    
    % Unique excluded patients
    excludedPatients = unique(excludedQC.PatientID);
    fprintf('\nUnique patients with ANY excluded observation: %d\n', length(excludedPatients));
end

% Apply exclusion to matrices
fprintf('\n--- APPLYING EXCLUSIONS ---\n');
fprintf('Before: %d observations\n', nFiles);

IRi_matrix(lowIRi_flag, :) = NaN;
MANi_matrix(lowIRi_flag, :) = NaN;
Vari_matrix(lowIRi_flag, :) = NaN;
NumActiveVoxels_matrix(lowIRi_flag, :) = NaN;

% Update ExtractionSuccess to reflect QC exclusions
ExtractionSuccess(lowIRi_flag) = false;

nAfterQC = sum(~lowIRi_flag);
fprintf('After QC: %d observations (%.1f%% retained)\n\n', ...
    nAfterQC, 100*nAfterQC/nFiles);

%% METRIC QC: VERIFY IRi SUMS AND RANGES (POST-EXCLUSION)
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Metric QC: Verifying extraction (POST-EXCLUSION)\n');
fprintf('-----------------------------------------------------------------\n');

% IRi should sum to approximately 1 across ICNs (it's relative involvement)
IRi_sums_clean = sum(IRi_matrix, 2, 'omitnan');
validIdx_QC = ~isnan(IRi_sums_clean);

fprintf('IRi row sums (after QC, should be ≥%.1f):\n', IRi_threshold);
fprintf('  Mean: %.3f\n', mean(IRi_sums_clean(validIdx_QC), 'omitnan'));
fprintf('  Range: %.3f - %.3f\n', min(IRi_sums_clean(validIdx_QC)), max(IRi_sums_clean(validIdx_QC)));
fprintf('  Note: <1 means some activation outside defined ICNs\n\n');

% MANi should be 0-1
fprintf('MANi range (should be 0-1):\n');
fprintf('  Min: %.4f\n', min(MANi_matrix(:), [], 'omitnan'));
fprintf('  Max: %.4f\n', max(MANi_matrix(:), [], 'omitnan'));
fprintf('\n');

% Vari is unbounded but check for outliers
fprintf('Vari distribution:\n');
fprintf('  Mean: %.3f\n', mean(Vari_matrix(:), 'omitnan'));
fprintf('  Median: %.3f\n', median(Vari_matrix(:), 'omitnan'));
fprintf('  Range: %.3f - %.3f\n', min(Vari_matrix(:), [], 'omitnan'), max(Vari_matrix(:), [], 'omitnan'));
fprintf('  Note: Unbounded, higher = more heterogeneous activation\n\n');

%% MERGE WITH CLINICAL DATA
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Merging with clinical data\n');
fprintf('-----------------------------------------------------------------\n');

% Initialize clinical columns
WAB_AQ = NaN(nFiles, 1);
WAB_Type = cell(nFiles, 1);
Age_At_Stroke = NaN(nFiles, 1);
Sex = cell(nFiles, 1);
Days_Post_Stroke = NaN(nFiles, 1);
AphasiaType_Folder = cell(nFiles, 1);
WAB_Type_Source = cell(nFiles, 1);

% Merge clinical data by PatientID
matchCount = 0;
for f = 1:nFiles
    pid = PatientID{f};
    
    % Find in clinical lookup
    clinIdx = find(strcmp(clinicalLookup.PatientID, pid), 1);
    
    if ~isempty(clinIdx)
        WAB_AQ(f) = clinicalLookup.WAB_AQ(clinIdx);
        WAB_Type{f} = clinicalLookup.WAB_Type{clinIdx};
        Age_At_Stroke(f) = clinicalLookup.Age_At_Stroke(clinIdx);
        Sex{f} = clinicalLookup.Sex{clinIdx};
        Days_Post_Stroke(f) = clinicalLookup.Days_Post_Stroke(clinIdx);
        matchCount = matchCount + 1;
    else
        WAB_Type{f} = '';
        Sex{f} = '';
    end
    
    % Get aphasia type from folder structure (con inventory)
    conIdx = find(strcmp(conInventory.PatientID, pid) & ...
                  strcmp(conInventory.ContrastNum, ContrastNum{f}), 1);
    if ~isempty(conIdx)
        AphasiaType_Folder{f} = conInventory.AphasiaType{conIdx};
    else
        AphasiaType_Folder{f} = '';
    end
end

fprintf('Clinical data matched: %d/%d files (%.1f%%)\n', matchCount, nFiles, (matchCount/nFiles)*100);

% Determine primary aphasia type column and track source
AphasiaType = cell(nFiles, 1);
for f = 1:nFiles
    if ~isempty(AphasiaType_Folder{f})
        AphasiaType{f} = AphasiaType_Folder{f};
        WAB_Type_Source{f} = 'folder';
    elseif ~isempty(WAB_Type{f})
        AphasiaType{f} = WAB_Type{f};
        WAB_Type_Source{f} = 'clinical';
    else
        AphasiaType{f} = '';
        WAB_Type_Source{f} = 'missing';
    end
end

fprintf('Aphasia type assignment: %d from folder, %d from clinical, %d missing\n\n', ...
        sum(strcmp(WAB_Type_Source, 'folder')), ...
        sum(strcmp(WAB_Type_Source, 'clinical')), ...
        sum(strcmp(WAB_Type_Source, 'missing')));

%% ASSIGN GROUP ROLE (APHASIA VS STROKE CONTROL)
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Assigning GroupRole (Aphasia vs Stroke Control)\n');
fprintf('-----------------------------------------------------------------\n');

GroupRole = cell(nFiles, 1);
for f = 1:nFiles
    if strcmp(AphasiaType{f}, 'Unspecified') || strcmp(WAB_Type{f}, 'None')
        GroupRole{f} = 'Stroke Control (No Aphasia)';
    elseif ~isempty(AphasiaType{f})
        GroupRole{f} = 'Aphasia';
    else
        GroupRole{f} = 'Unknown';
    end
end

nAphasia = sum(strcmp(GroupRole, 'Aphasia'));
nStrokeControl = sum(strcmp(GroupRole, 'Stroke Control (No Aphasia)'));
nUnknown = sum(strcmp(GroupRole, 'Unknown'));

fprintf('GroupRole assignment:\n');
fprintf('  Aphasia: %d observations\n', nAphasia);
fprintf('  Stroke Control (No Aphasia): %d observations\n', nStrokeControl);
fprintf('  Unknown: %d observations\n\n', nUnknown);

%% CREATE LONG-FORMAT TABLE
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Creating long-format master table\n');
fprintf('-----------------------------------------------------------------\n');

%% V3 CHANGE: Add QC exclusion flag to table
masterLong = table(PatientID, SessionID, ContrastNum, ContrastLabel, ...
                   AphasiaType, GroupRole, WAB_AQ, WAB_Type, WAB_Type_Source, ...
                   Age_At_Stroke, Sex, Days_Post_Stroke, ...
                   ExtractionSuccess, lowIRi_flag, IRi_sums, ...
                   'VariableNames', {'PatientID', 'SessionID', 'ContrastNum', 'ContrastLabel', ...
                                     'AphasiaType', 'GroupRole', 'WAB_AQ', 'WAB_Type', 'WAB_Type_Source', ...
                                     'Age_At_Stroke', 'Sex', 'Days_Post_Stroke', ...
                                     'ExtractionSuccess', 'IRi_QC_Excluded', 'IRi_Sum'});

% Add CORRECTED ICN metrics with clear column names
for n = 1:nICN
    masterLong.(sprintf('ICN%02d_IRi', n)) = IRi_matrix(:, n);
    masterLong.(sprintf('ICN%02d_MANi', n)) = MANi_matrix(:, n);
    masterLong.(sprintf('ICN%02d_Vari', n)) = Vari_matrix(:, n);
end

fprintf('Long-format table created: %d rows × %d columns\n', height(masterLong), width(masterLong));

%% V3 CHANGE: Filter to successful extractions AND passing QC
masterLong = masterLong(~masterLong.IRi_QC_Excluded, :);
fprintf('After filtering (extraction success + IRi QC): %d rows\n\n', height(masterLong));

%% CREATE WIDE-FORMAT TABLE (ONE ROW PER PATIENT)
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Creating wide-format master table\n');
fprintf('-----------------------------------------------------------------\n');

% Get unique patients from CLEANED long table
uniquePts = unique(masterLong.PatientID);
nPatients = length(uniquePts);

% Initialize wide table
masterWide = table('Size', [nPatients, 0]);

% Add patient-level columns
masterWide.PatientID = uniquePts;
masterWide.SessionID = cell(nPatients, 1);
masterWide.AphasiaType = cell(nPatients, 1);
masterWide.GroupRole = cell(nPatients, 1);
masterWide.WAB_AQ = NaN(nPatients, 1);
masterWide.WAB_Type = cell(nPatients, 1);
masterWide.WAB_Type_Source = cell(nPatients, 1);
masterWide.Age_At_Stroke = NaN(nPatients, 1);
masterWide.Sex = cell(nPatients, 1);
masterWide.Days_Post_Stroke = NaN(nPatients, 1);

% Initialize ICN columns for both contrasts with CORRECTED metric names
for n = 1:nICN
    % Con_0001 (Task > Rest)
    masterWide.(sprintf('C1_ICN%02d_IRi', n)) = NaN(nPatients, 1);
    masterWide.(sprintf('C1_ICN%02d_MANi', n)) = NaN(nPatients, 1);
    masterWide.(sprintf('C1_ICN%02d_Vari', n)) = NaN(nPatients, 1);
    
    % Con_0002 (Naming > Abstract)
    masterWide.(sprintf('C2_ICN%02d_IRi', n)) = NaN(nPatients, 1);
    masterWide.(sprintf('C2_ICN%02d_MANi', n)) = NaN(nPatients, 1);
    masterWide.(sprintf('C2_ICN%02d_Vari', n)) = NaN(nPatients, 1);
end

% Flag columns
masterWide.HasCon0001 = false(nPatients, 1);
masterWide.HasCon0002 = false(nPatients, 1);

% Fill wide table
for p = 1:nPatients
    pid = uniquePts{p};
    patientRows = strcmp(masterLong.PatientID, pid);
    patientData = masterLong(patientRows, :);
    
    % Get first row for clinical data
    masterWide.SessionID{p} = patientData.SessionID{1};
    masterWide.AphasiaType{p} = patientData.AphasiaType{1};
    masterWide.GroupRole{p} = patientData.GroupRole{1};
    masterWide.WAB_AQ(p) = patientData.WAB_AQ(1);
    masterWide.WAB_Type{p} = patientData.WAB_Type{1};
    masterWide.WAB_Type_Source{p} = patientData.WAB_Type_Source{1};
    masterWide.Age_At_Stroke(p) = patientData.Age_At_Stroke(1);
    masterWide.Sex{p} = patientData.Sex{1};
    masterWide.Days_Post_Stroke(p) = patientData.Days_Post_Stroke(1);
    
    % Con_0001 data
    con1Idx = strcmp(patientData.ContrastNum, '0001');
    if any(con1Idx)
        masterWide.HasCon0001(p) = true;
        con1Data = patientData(con1Idx, :);
        for n = 1:nICN
            masterWide.(sprintf('C1_ICN%02d_IRi', n))(p) = con1Data.(sprintf('ICN%02d_IRi', n));
            masterWide.(sprintf('C1_ICN%02d_MANi', n))(p) = con1Data.(sprintf('ICN%02d_MANi', n));
            masterWide.(sprintf('C1_ICN%02d_Vari', n))(p) = con1Data.(sprintf('ICN%02d_Vari', n));
        end
    end
    
    % Con_0002 data
    con2Idx = strcmp(patientData.ContrastNum, '0002');
    if any(con2Idx)
        masterWide.HasCon0002(p) = true;
        con2Data = patientData(con2Idx, :);
        for n = 1:nICN
            masterWide.(sprintf('C2_ICN%02d_IRi', n))(p) = con2Data.(sprintf('ICN%02d_IRi', n));
            masterWide.(sprintf('C2_ICN%02d_MANi', n))(p) = con2Data.(sprintf('ICN%02d_MANi', n));
            masterWide.(sprintf('C2_ICN%02d_Vari', n))(p) = con2Data.(sprintf('ICN%02d_Vari', n));
        end
    end
end

fprintf('Wide-format table created: %d patients × %d columns\n', height(masterWide), width(masterWide));
fprintf('Patients with Con_0001: %d\n', sum(masterWide.HasCon0001));
fprintf('Patients with Con_0002: %d\n', sum(masterWide.HasCon0002));
fprintf('Patients with both contrasts: %d\n\n', sum(masterWide.HasCon0001 & masterWide.HasCon0002));

%% CREATE NETWORK LABELS REFERENCE TABLE
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Creating network labels reference\n');
fprintf('-----------------------------------------------------------------\n');

networkRef = table((1:nICN)', trueICN_idx', trueICN_labels', ...
                   'VariableNames', {'ICN_Number', 'Atlas_Index', 'Label'});

fprintf('Network reference table: %d ICNs\n\n', height(networkRef));
disp(networkRef);
fprintf('\n');

%% SUMMARY STATISTICS
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Summary statistics (POST-IRi QC)\n');
fprintf('-----------------------------------------------------------------\n');

% Group summary
fprintf('=== GROUP ROLE SUMMARY ===\n');
nAphasiaPatients = sum(strcmp(masterWide.GroupRole, 'Aphasia'));
nStrokeControlPatients = sum(strcmp(masterWide.GroupRole, 'Stroke Control (No Aphasia)'));
fprintf('Aphasia patients: N = %d\n', nAphasiaPatients);
fprintf('Stroke Controls (No Aphasia): N = %d\n', nStrokeControlPatients);
fprintf('Total: N = %d\n\n', nPatients);

% By aphasia type
fprintf('=== APHASIA TYPE BREAKDOWN ===\n');
aphasiaTypes = unique(masterWide.AphasiaType(~cellfun(@isempty, masterWide.AphasiaType)));
for a = 1:length(aphasiaTypes)
    aType = aphasiaTypes{a};
    n = sum(strcmp(masterWide.AphasiaType, aType));
    meanWAB = mean(masterWide.WAB_AQ(strcmp(masterWide.AphasiaType, aType)), 'omitnan');
    typeIdx = find(strcmp(masterWide.AphasiaType, aType), 1);
    groupRole = masterWide.GroupRole{typeIdx};
    fprintf('  %-15s: N=%3d, Mean WAB-AQ=%5.1f  [%s]\n', aType, n, meanWAB, groupRole);
end
fprintf('\n');

% WAB-AQ by GroupRole
fprintf('=== WAB-AQ BY GROUP ROLE ===\n');
aphasiaWAB = masterWide.WAB_AQ(strcmp(masterWide.GroupRole, 'Aphasia'));
controlWAB = masterWide.WAB_AQ(strcmp(masterWide.GroupRole, 'Stroke Control (No Aphasia)'));

fprintf('Aphasia (N=%d):\n', sum(~isnan(aphasiaWAB)));
fprintf('  Mean (SD): %.1f (%.1f)\n', mean(aphasiaWAB, 'omitnan'), std(aphasiaWAB, 'omitnan'));
fprintf('  Range: %.1f - %.1f\n', min(aphasiaWAB), max(aphasiaWAB));

fprintf('Stroke Controls (N=%d):\n', sum(~isnan(controlWAB)));
fprintf('  Mean (SD): %.1f (%.1f)\n', mean(controlWAB, 'omitnan'), std(controlWAB, 'omitnan'));
fprintf('  Range: %.1f - %.1f\n\n', min(controlWAB), max(controlWAB));

% Metric summaries
fprintf('=== METRIC SUMMARIES (C2 - Naming>Abstract) ===\n');

% Get C2 IRi matrix for aphasia patients
aphasiaIdx = strcmp(masterWide.GroupRole, 'Aphasia');
fprintf('IRi (Relative Spatial Involvement):\n');
C2_IRi = NaN(nAphasiaPatients, nICN);
for n = 1:nICN
    C2_IRi(:, n) = masterWide.(sprintf('C2_ICN%02d_IRi', n))(aphasiaIdx);
end
fprintf('  Mean across ICNs: %.4f\n', mean(C2_IRi(:), 'omitnan'));
fprintf('  SD across ICNs: %.4f\n', std(C2_IRi(:), 'omitnan'));

fprintf('MANi (Normalized Mean Activation):\n');
C2_MANi = NaN(nAphasiaPatients, nICN);
for n = 1:nICN
    C2_MANi(:, n) = masterWide.(sprintf('C2_ICN%02d_MANi', n))(aphasiaIdx);
end
fprintf('  Mean across ICNs: %.4f\n', mean(C2_MANi(:), 'omitnan'));
fprintf('  SD across ICNs: %.4f\n', std(C2_MANi(:), 'omitnan'));

fprintf('Vari (Activation Variance):\n');
C2_Vari = NaN(nAphasiaPatients, nICN);
for n = 1:nICN
    C2_Vari(:, n) = masterWide.(sprintf('C2_ICN%02d_Vari', n))(aphasiaIdx);
end
fprintf('  Mean across ICNs: %.3f\n', mean(C2_Vari(:), 'omitnan'));
fprintf('  SD across ICNs: %.3f\n', std(C2_Vari(:), 'omitnan'));
fprintf('\n');

%% SAVE OUTPUTS
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Saving outputs (V3 filenames)\n');
fprintf('-----------------------------------------------------------------\n');

%% V3 CHANGE: All filenames updated to v3

% Save long format
longMatFile = fullfile(outputDir, 'ARC_03b_v3_Master_Long.mat');
save(longMatFile, 'masterLong', 'trueICN_labels', 'trueICN_idx', 'nICN', ...
     'IRi_threshold', 'IRi_QC');
fprintf('Saved: %s\n', longMatFile);

longCsvFile = fullfile(outputDir, 'ARC_03b_v3_Master_Long.csv');
writetable(masterLong, longCsvFile);
fprintf('Saved: %s\n', longCsvFile);

% Save wide format
wideMatFile = fullfile(outputDir, 'ARC_03b_v3_Master_Wide.mat');
save(wideMatFile, 'masterWide', 'trueICN_labels', 'trueICN_idx', 'nICN', ...
     'IRi_threshold');
fprintf('Saved: %s\n', wideMatFile);

wideCsvFile = fullfile(outputDir, 'ARC_03b_v3_Master_Wide.csv');
writetable(masterWide, wideCsvFile);
fprintf('Saved: %s\n', wideCsvFile);

% Save network reference
networkRefFile = fullfile(outputDir, 'ARC_03b_v3_Network_Labels.csv');
writetable(networkRef, networkRefFile);
fprintf('Saved: %s\n', networkRefFile);

% Save raw matrices (convenient for analysis) - V3: Include QC info
matricesFile = fullfile(outputDir, 'ARC_03b_v3_ICN_Matrices.mat');
save(matricesFile, 'IRi_matrix', 'MANi_matrix', 'Vari_matrix', 'NumActiveVoxels_matrix', ...
     'PatientID', 'ContrastNum', 'AphasiaType', 'GroupRole', 'WAB_AQ', ...
     'trueICN_labels', 'trueICN_idx', 'nICN', ...
     'IRi_threshold', 'lowIRi_flag', 'IRi_sums', 'IRi_QC');
fprintf('Saved: %s\n', matricesFile);

%% GENERATE SUMMARY REPORT
% =========================================================================
reportFile = fullfile(outputDir, 'ARC_03b_v3_Summary.txt');
fid = fopen(reportFile, 'w');

fprintf(fid, '=================================================================\n');
fprintf(fid, 'ARC_03b_v3 COMPILE RESULTS SUMMARY (IRi QC ADDED)\n');
fprintf(fid, '=================================================================\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));

fprintf(fid, 'V3 CHANGE: IRi QUALITY CONTROL\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Criterion: IRi sum < %.1f (>%.0f%% activation outside ICNs)\n', ...
        IRi_threshold, (1 - IRi_threshold) * 100);
fprintf(fid, 'Status: PRINCIPLED criterion (not data-driven)\n');
fprintf(fid, 'Rationale: Poor atlas coverage or data quality issues\n\n');

fprintf(fid, 'QC RESULTS\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Total observations: %d\n', nFiles);
fprintf(fid, 'Excluded (low IRi or extraction failure): %d (%.1f%%)\n', ...
        sum(lowIRi_flag), 100*sum(lowIRi_flag)/nFiles);
fprintf(fid, 'Retained: %d (%.1f%%)\n\n', nAfterQC, 100*nAfterQC/nFiles);

fprintf(fid, 'See ARC_03b_v3_IRi_QC_Report.csv for full exclusion details.\n\n');

fprintf(fid, 'CORRECTED METRICS\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'IRi:  ICNiRelativeInvolvement (was PropActiveVoxels)\n');
fprintf(fid, '      Measures: Proportion of whole-brain activation in this ICN\n');
fprintf(fid, '      Range: 0-1, sums to ~1 across ICNs\n');
fprintf(fid, '      PRIMARY METRIC - sensitive to extensive activations\n\n');
fprintf(fid, 'MANi: NormalisedMeanICNiActivation (was MeanActiveVoxels)\n');
fprintf(fid, '      Measures: Normalized mean Z-score within ICN\n');
fprintf(fid, '      Range: 0-1\n');
fprintf(fid, '      SECONDARY - more variable, sensitive to focal activations\n\n');
fprintf(fid, 'Vari: var(ZActiveVoxels) (was RANi - incorrect)\n');
fprintf(fid, '      Measures: Variance of Z-scores within ICN\n');
fprintf(fid, '      Range: 0 to unbounded\n');
fprintf(fid, '      EXPLORATORY - activation heterogeneity\n\n');

fprintf(fid, 'DATA OVERVIEW (POST-QC)\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'xATL files processed: %d\n', nFiles);
fprintf(fid, 'After IRi QC: %d observations\n', height(masterLong));
fprintf(fid, 'Unique patients: %d\n', nPatients);
fprintf(fid, 'ICNs extracted: %d\n', nICN);
fprintf(fid, 'Metrics per ICN: 3 (IRi, MANi, Vari)\n\n');

fprintf(fid, 'METRIC QC (POST-EXCLUSION)\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'IRi row sums: Mean=%.3f, Range=%.3f-%.3f\n', ...
        mean(IRi_sums_clean(validIdx_QC), 'omitnan'), ...
        min(IRi_sums_clean(validIdx_QC)), max(IRi_sums_clean(validIdx_QC)));
fprintf(fid, 'MANi range: %.4f - %.4f\n', ...
        min(MANi_matrix(:), [], 'omitnan'), max(MANi_matrix(:), [], 'omitnan'));
fprintf(fid, 'Vari: Mean=%.3f, Median=%.3f, Range=%.3f-%.3f\n\n', ...
        mean(Vari_matrix(:), 'omitnan'), median(Vari_matrix(:), 'omitnan'), ...
        min(Vari_matrix(:), [], 'omitnan'), max(Vari_matrix(:), [], 'omitnan'));

fprintf(fid, 'GROUP ROLE SUMMARY (POST-QC)\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Aphasia patients: N = %d\n', nAphasiaPatients);
fprintf(fid, 'Stroke Controls (No Aphasia): N = %d\n', nStrokeControlPatients);
fprintf(fid, 'Total: N = %d\n\n', nPatients);

fprintf(fid, 'OUTPUT TABLES\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Long format: %d rows × %d columns\n', height(masterLong), width(masterLong));
fprintf(fid, 'Wide format: %d rows × %d columns\n\n', height(masterWide), width(masterWide));

fprintf(fid, 'PATIENTS BY APHASIA TYPE (POST-QC)\n');
fprintf(fid, '-----------------------------------------------------------------\n');
for a = 1:length(aphasiaTypes)
    aType = aphasiaTypes{a};
    n = sum(strcmp(masterWide.AphasiaType, aType));
    meanWAB = mean(masterWide.WAB_AQ(strcmp(masterWide.AphasiaType, aType)), 'omitnan');
    typeIdx = find(strcmp(masterWide.AphasiaType, aType), 1);
    groupRole = masterWide.GroupRole{typeIdx};
    fprintf(fid, '%-15s: N=%3d, Mean WAB-AQ=%5.1f  [%s]\n', aType, n, meanWAB, groupRole);
end

fprintf(fid, '\nWAB-AQ BY GROUP ROLE (POST-QC)\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Aphasia (N=%d): Mean=%.1f, SD=%.1f, Range=%.1f-%.1f\n', ...
        sum(~isnan(aphasiaWAB)), mean(aphasiaWAB, 'omitnan'), std(aphasiaWAB, 'omitnan'), ...
        min(aphasiaWAB), max(aphasiaWAB));
fprintf(fid, 'Stroke Controls (N=%d): Mean=%.1f, SD=%.1f, Range=%.1f-%.1f\n', ...
        sum(~isnan(controlWAB)), mean(controlWAB, 'omitnan'), std(controlWAB, 'omitnan'), ...
        min(controlWAB), max(controlWAB));

fprintf(fid, '\nNETWORK LABELS\n');
fprintf(fid, '-----------------------------------------------------------------\n');
for n = 1:nICN
    fprintf(fid, 'ICN%02d: %s\n', n, trueICN_labels{n});
end

fprintf(fid, '\nCLINICAL DATA COMPLETENESS (POST-QC)\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'WAB-AQ available: %d/%d (%.1f%%)\n', ...
        sum(~isnan(masterWide.WAB_AQ)), nPatients, (sum(~isnan(masterWide.WAB_AQ))/nPatients)*100);
fprintf(fid, 'Age available: %d/%d (%.1f%%)\n', ...
        sum(~isnan(masterWide.Age_At_Stroke)), nPatients, (sum(~isnan(masterWide.Age_At_Stroke))/nPatients)*100);
fprintf(fid, 'Days_Post_Stroke available: %d/%d (%.1f%%)\n', ...
        sum(~isnan(masterWide.Days_Post_Stroke)), nPatients, (sum(~isnan(masterWide.Days_Post_Stroke))/nPatients)*100);

fprintf(fid, '\nOUTPUT FILES\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'ARC_03b_v3_IRi_QC_Report.csv - QC exclusion details\n');
fprintf(fid, 'ARC_03b_v3_Master_Long.mat/csv - Long format (%d rows)\n', height(masterLong));
fprintf(fid, 'ARC_03b_v3_Master_Wide.mat/csv - Wide format (%d patients)\n', height(masterWide));
fprintf(fid, 'ARC_03b_v3_Network_Labels.csv - ICN reference\n');
fprintf(fid, 'ARC_03b_v3_ICN_Matrices.mat - Raw matrices with QC info\n');
fprintf(fid, 'ARC_03b_v3_Summary.txt - This report\n');

fprintf(fid, '\n=================================================================\n');
fprintf(fid, 'END OF REPORT\n');
fprintf(fid, '=================================================================\n');

fclose(fid);
fprintf('Saved: %s\n\n', reportFile);

%% COMPLETION
% =========================================================================
fprintf('=================================================================\n');
fprintf('ARC_03b_v3_Compile_Results.m - COMPLETE\n');
fprintf('=================================================================\n');
fprintf('End time: %s\n', datestr(now));
fprintf('\nV3 IRi QC APPLIED:\n');
fprintf('  Threshold: IRi sum < %.1f\n', IRi_threshold);
fprintf('  Excluded: %d/%d observations (%.1f%%)\n', ...
        sum(lowIRi_flag), nFiles, 100*sum(lowIRi_flag)/nFiles);
fprintf('  Retained: %d observations\n', nAfterQC);
fprintf('\nCORRECTED METRICS:\n');
fprintf('  IRi:  ICNiRelativeInvolvement (PRIMARY)\n');
fprintf('  MANi: NormalisedMeanICNiActivation (SECONDARY)\n');
fprintf('  Vari: var(ZActiveVoxels) (EXPLORATORY)\n');
fprintf('\nOutputs saved to: %s\n', outputDir);
fprintf('\nOutput files:\n');
fprintf('  - QC Report: ARC_03b_v3_IRi_QC_Report.csv\n');
fprintf('  - Long format: ARC_03b_v3_Master_Long.mat/csv (%d rows)\n', height(masterLong));
fprintf('  - Wide format: ARC_03b_v3_Master_Wide.mat/csv (%d patients)\n', height(masterWide));
fprintf('  - Network labels: ARC_03b_v3_Network_Labels.csv\n');
fprintf('  - ICN matrices: ARC_03b_v3_ICN_Matrices.mat\n');
fprintf('  - Summary: ARC_03b_v3_Summary.txt\n');
fprintf('\nGroup breakdown (POST-QC):\n');
fprintf('  - Aphasia: N = %d\n', nAphasiaPatients);
fprintf('  - Stroke Controls (No Aphasia): N = %d\n', nStrokeControlPatients);
fprintf('\nNext: Run ARC_03c_v3_Statistical_Analysis.m with cleaned data\n');
fprintf('=================================================================\n');