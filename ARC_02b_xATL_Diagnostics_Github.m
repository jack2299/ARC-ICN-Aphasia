%% ARC_02b_xATL_Diagnostics.m
% =========================================================================
% SCRIPT 2b: xATL Structure Diagnostic
% =========================================================================
% Purpose:
%   - Load and inspect xATL.mat files from ICN_Atlas extraction
%   - Document actual field structure
%   - Identify where network engagement metrics are stored
%   - Create reference for Script 3 parsing
%
% Output:
%   - Console output showing complete xATL structure
%   - ARC_02b_xATL_Structure.txt (reference document)
% =========================================================================

%% SETUP
clear; clc;

% Load configuration
config = config_local();

fprintf('=================================================================\n');
fprintf('ARC_02b_xATL_Diagnostics.m - xATL Structure Inspection\n');
fprintf('=================================================================\n');
fprintf('Start time: %s\n\n', datestr(now));

%% DEFINE PATHS
icnOutputDir = fullfile(config.rootDir, 'ARC_02_Output');
outputDir = icnOutputDir;

%% FIND xATL FILES
fprintf('-----------------------------------------------------------------\n');
fprintf('STEP 1: Locating xATL files\n');
fprintf('-----------------------------------------------------------------\n');

xatlFiles = dir(fullfile(icnOutputDir, 'ICN_*.mat'));
fprintf('Found %d xATL files in ARC_02_Output\n\n', length(xatlFiles));

if isempty(xatlFiles)
    % Try alternative location (original Conduction folders)
    altDir = fullfile(config.rootDir, 'Conduction aphasia_Group_Analysis_Con_001');
    xatlFiles = dir(fullfile(altDir, '*.xATL.mat'));
    if ~isempty(xatlFiles)
        fprintf('Found xATL files in alternative location:\n  %s\n\n', altDir);
        icnOutputDir = altDir;
    else
        error('No xATL files found. Run Script 2 first.');
    end
end

%% LOAD FIRST FILE FOR INSPECTION
fprintf('-----------------------------------------------------------------\n');
fprintf('STEP 2: Loading first xATL file for inspection\n');
fprintf('-----------------------------------------------------------------\n');

testFile = fullfile(icnOutputDir, xatlFiles(1).name);
fprintf('Loading: %s\n\n', xatlFiles(1).name);

data = load(testFile);

%% RECURSIVE STRUCTURE EXPLORER
fprintf('-----------------------------------------------------------------\n');
fprintf('STEP 3: Complete structure exploration\n');
fprintf('-----------------------------------------------------------------\n');

% Open output file
reportFile = fullfile(outputDir, 'ARC_02b_xATL_Structure.txt');
fid = fopen(reportFile, 'w');

fprintf(fid, '=================================================================\n');
fprintf(fid, 'xATL STRUCTURE REFERENCE DOCUMENT\n');
fprintf(fid, '=================================================================\n');
fprintf(fid, 'Generated: %s\n', datestr(now));
fprintf(fid, 'Source file: %s\n\n', xatlFiles(1).name);

% Explore structure recursively
exploreStruct(data, 'data', 0, fid);

fprintf(fid, '\n=================================================================\n');
fprintf(fid, 'END OF STRUCTURE DOCUMENT\n');
fprintf(fid, '=================================================================\n');

fclose(fid);
fprintf('\nStructure saved to: %s\n\n', reportFile);

%% STEP 4: TARGETED INSPECTION - LIKELY NETWORK DATA LOCATIONS
fprintf('-----------------------------------------------------------------\n');
fprintf('STEP 4: Targeted inspection of likely network data fields\n');
fprintf('-----------------------------------------------------------------\n');

% Check what's at top level
topFields = fieldnames(data);
fprintf('Top-level fields: %s\n\n', strjoin(topFields, ', '));

% Navigate to xATL if nested
if isfield(data, 'xATL')
    xATL = data.xATL;
    fprintf('Navigating into data.xATL...\n\n');
elseif isfield(data, 'xatl')
    xATL = data.xatl;
    fprintf('Navigating into data.xatl...\n\n');
else
    xATL = data;
    fprintf('Using top-level structure as xATL...\n\n');
end

% List xATL fields
if isstruct(xATL)
    xatlFields = fieldnames(xATL);
    fprintf('xATL fields:\n');
    for i = 1:length(xatlFields)
        val = xATL.(xatlFields{i});
        fprintf('  .%s: [%s] %s\n', xatlFields{i}, class(val), mat2str(size(val)));
    end
    fprintf('\n');
end

%% STEP 5: INSPECT KEY FIELDS
fprintf('-----------------------------------------------------------------\n');
fprintf('STEP 5: Detailed inspection of key fields\n');
fprintf('-----------------------------------------------------------------\n');

% Check Labels (network names)
if isfield(xATL, 'Labels')
    fprintf('LABELS (Network Names):\n');
    labels = xATL.Labels;
    if iscell(labels)
        for i = 1:min(length(labels), 25)
            if iscell(labels{i})
                fprintf('  %d: %s\n', i, labels{i}{1});
            else
                fprintf('  %d: %s\n', i, char(labels{i}));
            end
        end
        if length(labels) > 25
            fprintf('  ... and %d more\n', length(labels)-25);
        end
    end
    fprintf('\n');
end

% Check Clu (Cluster data)
if isfield(xATL, 'Clu')
    fprintf('CLU (Cluster Data):\n');
    Clu = xATL.Clu;
    if isstruct(Clu)
        cluFields = fieldnames(Clu);
        fprintf('  Clu fields: %s\n', strjoin(cluFields, ', '));
        for i = 1:length(cluFields)
            val = Clu.(cluFields{i});
            fprintf('    .%s: [%s] %s\n', cluFields{i}, class(val), mat2str(size(val)));
            if isnumeric(val) && numel(val) <= 20
                fprintf('      Values: %s\n', mat2str(val));
            elseif isnumeric(val) && numel(val) > 20
                fprintf('      First 10: %s...\n', mat2str(val(1:min(10,numel(val)))));
            end
        end
    elseif iscell(Clu)
        fprintf('  Clu is cell array: %s\n', mat2str(size(Clu)));
        for c = 1:min(length(Clu), 3)
            fprintf('  Clu{%d}: [%s] %s\n', c, class(Clu{c}), mat2str(size(Clu{c})));
        end
    end
    fprintf('\n');
end

% Check Atl (Atlas data)
if isfield(xATL, 'Atl')
    fprintf('ATL (Atlas Data):\n');
    Atl = xATL.Atl;
    if isstruct(Atl)
        atlFields = fieldnames(Atl);
        fprintf('  Atl fields: %s\n', strjoin(atlFields, ', '));
        for i = 1:length(atlFields)
            val = Atl.(atlFields{i});
            fprintf('    .%s: [%s] %s\n', atlFields{i}, class(val), mat2str(size(val)));
            if isnumeric(val) && numel(val) <= 20
                fprintf('      Values: %s\n', mat2str(val));
            elseif isnumeric(val) && numel(val) > 20
                fprintf('      First 10: %s...\n', mat2str(val(1:min(10,numel(val)))));
            end
        end
    end
    fprintf('\n');
end

% Check for direct network engagement fields
directFields = {'nClu', 'nVox', 'meanZ', 'maxZ', 'Vol', 'Mean', 'Max', 'Peak', ...
                'voxels', 'volume', 'activation', 'engagement'};
fprintf('CHECKING DIRECT ENGAGEMENT FIELDS:\n');
for i = 1:length(directFields)
    if isfield(xATL, directFields{i})
        val = xATL.(directFields{i});
        fprintf('  .%s: [%s] %s\n', directFields{i}, class(val), mat2str(size(val)));
        if isnumeric(val) && numel(val) <= 30
            fprintf('    Values: %s\n', mat2str(val));
        end
    end
end
fprintf('\n');

%% STEP 6: CROSS-CHECK WITH ORIGINAL CONDUCTION xATL
fprintf('-----------------------------------------------------------------\n');
fprintf('STEP 6: Compare with original Conduction xATL (if different)\n');
fprintf('-----------------------------------------------------------------\n');

origXatlPath = fullfile(config.rootDir, 'Conduction aphasia_Group_Analysis_Con_002', 'spmT_0001.nii.xATL.mat');
if exist(origXatlPath, 'file')
    fprintf('Loading original Conduction xATL for comparison...\n');
    origData = load(origXatlPath);
    
    origFields = fieldnames(origData);
    fprintf('Original xATL top-level fields: %s\n\n', strjoin(origFields, ', '));
    
    if isfield(origData, 'xATL')
        origXATL = origData.xATL;
        fprintf('Original xATL subfields:\n');
        origSubFields = fieldnames(origXATL);
        for i = 1:length(origSubFields)
            val = origXATL.(origSubFields{i});
            fprintf('  .%s: [%s] %s\n', origSubFields{i}, class(val), mat2str(size(val)));
        end
    end
else
    fprintf('Original Conduction xATL not found at expected path.\n');
end

%% STEP 7: EXTRACT SAMPLE DATA FOR VERIFICATION
fprintf('\n-----------------------------------------------------------------\n');
fprintf('STEP 7: Attempting to extract sample network engagement data\n');
fprintf('-----------------------------------------------------------------\n');

extracted = false;

if isfield(xATL, 'Atl') && isstruct(xATL.Atl)
    Atl = xATL.Atl;
    if isfield(Atl, 'nVox') || isfield(Atl, 'Vol')
        fprintf('PATTERN 1: Network data in xATL.Atl\n');
        if isfield(Atl, 'nVox')
            fprintf('  Voxel counts (nVox): %s\n', mat2str(Atl.nVox));
        end
        if isfield(Atl, 'Vol')
            fprintf('  Volume (Vol): %s\n', mat2str(Atl.Vol));
        end
        if isfield(Atl, 'Mean')
            fprintf('  Mean activation (Mean): %s\n', mat2str(Atl.Mean));
        end
        if isfield(Atl, 'Max')
            fprintf('  Max activation (Max): %s\n', mat2str(Atl.Max));
        end
        extracted = true;
    end
end

if isfield(xATL, 'Clu') && isstruct(xATL.Clu) && ~extracted
    Clu = xATL.Clu;
    fprintf('PATTERN 2: Cluster data in xATL.Clu\n');
    cluFields = fieldnames(Clu);
    fprintf('  Available fields: %s\n', strjoin(cluFields, ', '));
    
    if isfield(Clu, 'Atlas')
        fprintf('  Atlas assignments: %s\n', mat2str(size(Clu.Atlas)));
    end
    if isfield(Clu, 'nVox')
        fprintf('  Voxels per cluster: %s\n', mat2str(Clu.nVox));
    end
    extracted = true;
end

if ~extracted
    fprintf('PATTERN 3: Checking direct arrays at xATL level...\n');
    if isfield(xATL, 'nClu')
        fprintf('  nClu (number of clusters): %s\n', mat2str(xATL.nClu));
    end
end

%% STEP 8: SUMMARY AND RECOMMENDATIONS
fprintf('\n-----------------------------------------------------------------\n');
fprintf('STEP 8: Summary and parsing recommendations\n');
fprintf('-----------------------------------------------------------------\n');

fprintf('Based on inspection, recommended parsing approach:\n\n');

if isfield(xATL, 'Atl') && isstruct(xATL.Atl) && isfield(xATL.Atl, 'nVox')
    fprintf('STRUCTURE TYPE: Atlas-level aggregation\n');
    fprintf('PARSING CODE:\n');
    fprintf('  networkNames = xATL.Labels;\n');
    fprintf('  voxelCounts = xATL.Atl.nVox;\n');
    fprintf('  meanActivation = xATL.Atl.Mean;\n');
    fprintf('  peakActivation = xATL.Atl.Max;\n');
elseif isfield(xATL, 'Clu') && isstruct(xATL.Clu)
    fprintf('STRUCTURE TYPE: Cluster-level data (needs aggregation by atlas)\n');
    fprintf('PARSING CODE: Aggregate Clu data by atlas assignment\n');
else
    fprintf('STRUCTURE TYPE: Unknown - manual inspection required\n');
    fprintf('Check ARC_02b_xATL_Structure.txt for full field listing\n');
end

%% COMPLETION
fprintf('\n=================================================================\n');
fprintf('ARC_02b_xATL_Diagnostics.m - COMPLETE\n');
fprintf('=================================================================\n');
fprintf('End time: %s\n', datestr(now));
fprintf('Reference document: %s\n', reportFile);
fprintf('=================================================================\n');

%% HELPER FUNCTION: Recursive Structure Explorer
function exploreStruct(s, name, depth, fid)
    indent = repmat('  ', 1, depth);
    
    if isstruct(s)
        fields = fieldnames(s);
        for i = 1:length(fields)
            f = fields{i};
            val = s.(f);
            
            if depth == 0
                fprintf('%s.%s: [%s] %s\n', name, f, class(val), mat2str(size(val)));
                fprintf(fid, '%s.%s: [%s] %s\n', name, f, class(val), mat2str(size(val)));
            else
                fprintf('%s.%s: [%s] %s\n', indent, f, class(val), mat2str(size(val)));
                fprintf(fid, '%s.%s: [%s] %s\n', indent, f, class(val), mat2str(size(val)));
            end
            
            if isstruct(val) && depth < 4
                exploreStruct(val, [name '.' f], depth + 1, fid);
            elseif iscell(val) && ~isempty(val) && isstruct(val{1}) && depth < 3
                fprintf('%s  {1}: [%s] %s\n', indent, class(val{1}), mat2str(size(val{1})));
                fprintf(fid, '%s  {1}: [%s] %s\n', indent, class(val{1}), mat2str(size(val{1})));
                exploreStruct(val{1}, [name '.' f '{1}'], depth + 2, fid);
            elseif isnumeric(val) && numel(val) <= 10
                fprintf('%s    = %s\n', indent, mat2str(val));
                fprintf(fid, '%s    = %s\n', indent, mat2str(val));
            elseif iscell(val) && numel(val) <= 5 && all(cellfun(@ischar, val) | cellfun(@isstring, val))
                fprintf('%s    = {%s}\n', indent, strjoin(string(val), ', '));
                fprintf(fid, '%s    = {%s}\n', indent, strjoin(string(val), ', '));
            end
        end
    end
end