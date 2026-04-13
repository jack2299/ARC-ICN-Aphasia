%% ARC_03a_Individual_ICN_Extraction.m
% =========================================================================
% SCRIPT 3a: Individual Patient ICN Extraction
% =========================================================================
% Purpose:
%   - Run ICN_Atlas on each individual patient's contrast images
%   - Process both Con_0001 (Task>Rest) and Con_0002 (Naming>Abstract)
%   - First session only for statistical independence
%   - Save individual xATL.mat files to centralized output folder
%
% Input:
%   - ARC_02c_ConFile_Inventory.mat (from Script 2c)
%
% Output:
%   - Individual xATL files: ARC_03_Output/Individual_ICN/[PatientID]_[Session]_[Contrast].xATL.mat
%   - Extraction log: ARC_03a_Extraction_Log.csv
%   - Summary report: ARC_03a_Summary.txt
%
% Dependencies:
%   - SPM12 (on path)
%   - ICN_Atlas toolbox (on path)
%   - ARC_02c outputs
% =========================================================================

%% SETUP
clear; clc;

% Load configuration
config = config_local();

fprintf('=================================================================\n');
fprintf('ARC_03a_Individual_ICN_Extraction.m\n');
fprintf('Individual Patient ICN Network Extraction\n');
fprintf('=================================================================\n');
fprintf('Start time: %s\n\n', datestr(now));
startTime = tic;

%% DEFINE PATHS
% =========================================================================
% Get paths from config (user must add these to config_local.m)
if ~isfield(config, 'icnAtlasPath') || isempty(config.icnAtlasPath)
    error('config_local.m must define config.icnAtlasPath (path to ICN_atlas_public folder)');
end
if ~isfield(config, 'spmPath') || isempty(config.spmPath)
    error('config_local.m must define config.spmPath (path to SPM12 folder)');
end

% Add toolboxes to path
addpath(config.icnAtlasPath);
addpath(config.spmPath);

% Input from Script 2c
inventoryFile = fullfile(config.rootDir, 'ARC_02_Output', 'ARC_02c_ConFile_Inventory.mat');

% Output directory
outputDir = fullfile(config.rootDir, 'ARC_03_allSessionsCheckOutput');
icnOutputDir = fullfile(outputDir, 'Individual_ICN');

% Create output directories
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('Created: %s\n', outputDir);
end
if ~exist(icnOutputDir, 'dir')
    mkdir(icnOutputDir);
    fprintf('Created: %s\n', icnOutputDir);
end
fprintf('\n');

%% VERIFY DEPENDENCIES
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Verifying dependencies\n');
fprintf('-----------------------------------------------------------------\n');

if ~exist('ICN_atlas', 'file')
    error('ICN_Atlas toolbox not found. Check path: %s', config.icnAtlasPath);
else
    fprintf('ICN_Atlas toolbox: ✓\n');
end

if ~exist('spm', 'file')
    error('SPM12 not found. Check path: %s', config.spmPath);
else
    fprintf('SPM12: ✓\n');
end

if ~exist(inventoryFile, 'file')
    error('Inventory file not found: %s\nRun Script 2c first.', inventoryFile);
else
    fprintf('Inventory file: ✓\n');
end
fprintf('\n');

%% LOAD AND FILTER INVENTORY
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Loading and filtering inventory\n');
fprintf('-----------------------------------------------------------------\n');

load(inventoryFile, 'conInventory');
fprintf('Loaded inventory: %d total files\n', height(conInventory));

firstSessionIdx = conInventory.IsFirstSession;
analysisFiles = conInventory(firstSessionIdx, :);
fprintf('First session files: %d\n', height(analysisFiles));

contrastIdx = strcmp(analysisFiles.ContrastNum, '0001') | ...
              strcmp(analysisFiles.ContrastNum, '0002');
analysisFiles = analysisFiles(contrastIdx, :);
fprintf('After contrast filter (0001 + 0002): %d\n', height(analysisFiles));

nCon001 = sum(strcmp(analysisFiles.ContrastNum, '0001'));
nCon002 = sum(strcmp(analysisFiles.ContrastNum, '0002'));
fprintf('  Con_0001 (Task>Rest): %d files\n', nCon001);
fprintf('  Con_0002 (Naming>Abstract): %d files\n', nCon002);

uniquePatients = unique(analysisFiles.PatientID);
fprintf('Unique patients: %d\n\n', length(uniquePatients));

fprintf('Files by aphasia type:\n');
aphasiaTypes = unique(analysisFiles.AphasiaType);
for a = 1:length(aphasiaTypes)
    aType = aphasiaTypes{a};
    n001 = sum(strcmp(analysisFiles.AphasiaType, aType) & strcmp(analysisFiles.ContrastNum, '0001'));
    n002 = sum(strcmp(analysisFiles.AphasiaType, aType) & strcmp(analysisFiles.ContrastNum, '0002'));
    fprintf('  %-15s: Con_0001=%d, Con_0002=%d\n', aType, n001, n002);
end
fprintf('\n');

%% ICN_ATLAS SETTINGS
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('ICN_Atlas extraction settings\n');
fprintf('-----------------------------------------------------------------\n');

atlasName = 'BM20';           % BrainMap20 - 18 functional ICNs
atlasThreshold = 3;           % Z-threshold for atlas
mapThreshold = 3;             % Z-threshold for activation maps
processingMode = 'act';       % Activations only
normMode = 'zero';            % Zero-based normalization
withinSessionNorm = 0;        % No within-session normalization
clusterBased = 1;             % Cluster-based processing
outputMode = 'xa2f';          % Save xATL to file

fprintf('Atlas: %s (BrainMap20 - 18 ICNs)\n', atlasName);
fprintf('Atlas threshold: Z = %d\n', atlasThreshold);
fprintf('Map threshold: Z = %d\n', mapThreshold);
fprintf('Processing mode: %s\n', processingMode);
fprintf('Normalization: %s\n', normMode);
fprintf('\n');

%% CHECK FOR ALREADY-PROCESSED FILES
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Checking for already-processed files\n');
fprintf('-----------------------------------------------------------------\n');

analysisFiles.OutputFileName = cell(height(analysisFiles), 1);
analysisFiles.OutputFilePath = cell(height(analysisFiles), 1);
analysisFiles.AlreadyProcessed = false(height(analysisFiles), 1);

for i = 1:height(analysisFiles)
    pid = analysisFiles.PatientID{i};
    ses = analysisFiles.SessionID{i};
    con = analysisFiles.ContrastNum{i};
    
    outputFileName = sprintf('%s_ses-%s_con_%s.xATL.mat', pid, ses, con);
    outputFilePath = fullfile(icnOutputDir, outputFileName);
    
    analysisFiles.OutputFileName{i} = outputFileName;
    analysisFiles.OutputFilePath{i} = outputFilePath;
    analysisFiles.AlreadyProcessed(i) = exist(outputFilePath, 'file') == 2;
end

nAlreadyDone = sum(analysisFiles.AlreadyProcessed);
nToProcess = sum(~analysisFiles.AlreadyProcessed);

fprintf('Already processed: %d\n', nAlreadyDone);
fprintf('To process: %d\n\n', nToProcess);

if nAlreadyDone > 0 && nToProcess > 0
    fprintf('Skipping %d already-processed files.\n', nAlreadyDone);
    fprintf('To reprocess all, delete files in: %s\n\n', icnOutputDir);
end

%% EXTRACTION LOOP
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Starting ICN_Atlas extraction\n');
fprintf('-----------------------------------------------------------------\n');

if nToProcess == 0
    fprintf('All files already processed. Nothing to do.\n\n');
else
    fprintf('Processing %d files...\n', nToProcess);
    fprintf('Estimated time: %.0f - %.0f minutes\n\n', nToProcess*3/60, nToProcess*8/60);
end

extractionLog = table('Size', [height(analysisFiles), 7], ...
    'VariableTypes', {'string', 'string', 'string', 'string', 'logical', 'string', 'double'}, ...
    'VariableNames', {'PatientID', 'SessionID', 'ContrastNum', 'AphasiaType', ...
                      'Success', 'ErrorMessage', 'ProcessingTime'});

successCount = 0;
failCount = 0;
skipCount = 0;
progressInterval = max(1, round(nToProcess / 20));
lastProgress = 0;

for i = 1:height(analysisFiles)
    extractionLog.PatientID(i) = analysisFiles.PatientID{i};
    extractionLog.SessionID(i) = analysisFiles.SessionID{i};
    extractionLog.ContrastNum(i) = analysisFiles.ContrastNum{i};
    extractionLog.AphasiaType(i) = analysisFiles.AphasiaType{i};
    
    if analysisFiles.AlreadyProcessed(i)
        extractionLog.Success(i) = true;
        extractionLog.ErrorMessage(i) = "Skipped - already processed";
        extractionLog.ProcessingTime(i) = 0;
        skipCount = skipCount + 1;
        continue;
    end
    
    inputFile = analysisFiles.FilePath{i};
    outputFile = analysisFiles.OutputFilePath{i};
    
    currentProgress = successCount + failCount;
    if currentProgress - lastProgress >= progressInterval || currentProgress == 1
        elapsed = toc(startTime);
        if currentProgress > 0
            estimatedTotal = elapsed / currentProgress * nToProcess;
            remaining = estimatedTotal - elapsed;
            fprintf('[%s] Processing %d/%d (%.0f%%) - Est. remaining: %.1f min\n', ...
                    datestr(now, 'HH:MM:SS'), currentProgress, nToProcess, ...
                    (currentProgress/nToProcess)*100, remaining/60);
        end
        lastProgress = currentProgress;
    end
    
    fileTimer = tic;
    
    try
        if ~exist(inputFile, 'file')
            error('Input file not found: %s', inputFile);
        end
        
        ICN_atlas(atlasName, atlasThreshold, processingMode, mapThreshold, ...
                  normMode, withinSessionNorm, clusterBased, outputMode, {inputFile});
        
        sourceXATL = [inputFile, '.xATL.mat'];
        pause(0.5);
        
        if exist(sourceXATL, 'file')
            copyfile(sourceXATL, outputFile);
            extractionLog.Success(i) = true;
            extractionLog.ErrorMessage(i) = "";
            successCount = successCount + 1;
        else
            error('xATL file not created by ICN_Atlas');
        end
        
    catch ME
        extractionLog.Success(i) = false;
        extractionLog.ErrorMessage(i) = string(ME.message);
        failCount = failCount + 1;
        
        fprintf('  ✗ ERROR [%s ses-%s con_%s]: %s\n', ...
                analysisFiles.PatientID{i}, analysisFiles.SessionID{i}, ...
                analysisFiles.ContrastNum{i}, ME.message);
    end
    
    extractionLog.ProcessingTime(i) = toc(fileTimer);
end

%% EXTRACTION SUMMARY
% =========================================================================
fprintf('\n-----------------------------------------------------------------\n');
fprintf('Extraction complete\n');
fprintf('-----------------------------------------------------------------\n');

totalTime = toc(startTime);

fprintf('Results:\n');
fprintf('  Successful: %d\n', successCount);
fprintf('  Skipped (already done): %d\n', skipCount);
fprintf('  Failed: %d\n', failCount);
fprintf('  Total time: %.1f minutes\n', totalTime/60);

if successCount > 0
    avgTime = mean(extractionLog.ProcessingTime(extractionLog.Success & extractionLog.ProcessingTime > 0));
    fprintf('  Average time per file: %.1f seconds\n', avgTime);
end
fprintf('\n');

if failCount > 0
    fprintf('Failed extractions:\n');
    failedIdx = find(~extractionLog.Success & extractionLog.ErrorMessage ~= "Skipped - already processed");
    for f = 1:min(10, length(failedIdx))
        idx = failedIdx(f);
        fprintf('  %s ses-%s con_%s: %s\n', ...
                extractionLog.PatientID(idx), extractionLog.SessionID(idx), ...
                extractionLog.ContrastNum(idx), extractionLog.ErrorMessage(idx));
    end
    if length(failedIdx) > 10
        fprintf('  ... and %d more failures\n', length(failedIdx)-10);
    end
    fprintf('\n');
end

%% VERIFY OUTPUT FILES
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Verifying output files\n');
fprintf('-----------------------------------------------------------------\n');

outputFiles = dir(fullfile(icnOutputDir, '*.xATL.mat'));
fprintf('Total xATL files in output folder: %d\n', length(outputFiles));

con001Files = dir(fullfile(icnOutputDir, '*_con_0001.xATL.mat'));
con002Files = dir(fullfile(icnOutputDir, '*_con_0002.xATL.mat'));
fprintf('  Con_0001 files: %d\n', length(con001Files));
fprintf('  Con_0002 files: %d\n', length(con002Files));

outputPatients = regexp({outputFiles.name}, '^(\w+)_ses-', 'tokens', 'once');
outputPatients = [outputPatients{:}];
uniqueOutputPatients = unique(outputPatients);
fprintf('Unique patients with extractions: %d\n\n', length(uniqueOutputPatients));

%% SAVE OUTPUTS
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Saving outputs\n');
fprintf('-----------------------------------------------------------------\n');

% Save extraction log
logFile = fullfile(outputDir, 'ARC_03a_Extraction_Log.csv');
writetable(extractionLog, logFile);
fprintf('Saved: %s\n', logFile);

% Save summary report
reportFile = fullfile(outputDir, 'ARC_03a_Summary.txt');
fid = fopen(reportFile, 'w');

fprintf(fid, '=================================================================\n');
fprintf(fid, 'ARC_03a INDIVIDUAL ICN EXTRACTION SUMMARY\n');
fprintf(fid, '=================================================================\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));

fprintf(fid, 'INPUT\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Files requested: %d\n', height(analysisFiles));
fprintf(fid, '  Con_0001: %d\n', nCon001);
fprintf(fid, '  Con_0002: %d\n', nCon002);
fprintf(fid, 'Unique patients: %d\n\n', length(uniquePatients));

fprintf(fid, 'RESULTS\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Successful: %d\n', successCount);
fprintf(fid, 'Skipped (already done): %d\n', skipCount);
fprintf(fid, 'Failed: %d\n', failCount);
fprintf(fid, 'Total time: %.1f minutes\n', totalTime/60);

if successCount > 0
    fprintf(fid, 'Average time per file: %.1f seconds\n', avgTime);
end

fprintf(fid, '\nOUTPUT FILES\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'xATL files in: %s\n', icnOutputDir);
fprintf(fid, '  Con_0001: %d\n', length(con001Files));
fprintf(fid, '  Con_0002: %d\n', length(con002Files));
fprintf(fid, 'Unique patients: %d\n', length(uniqueOutputPatients));
fprintf(fid, '\nExtraction log: %s\n', logFile);

fprintf(fid, '\n=================================================================\n');
fprintf(fid, 'END OF REPORT\n');
fprintf(fid, '=================================================================\n');

fclose(fid);
fprintf('Saved: %s\n', reportFile);

%% COMPLETION
% =========================================================================
fprintf('\n=================================================================\n');
fprintf('ARC_03a_Individual_ICN_Extraction.m - COMPLETE\n');
fprintf('=================================================================\n');
fprintf('End time: %s\n', datestr(now));
fprintf('xATL files saved to: %s\n', icnOutputDir);
fprintf('Next: Run ARC_03b_Compile_Results.m\n');
fprintf('=================================================================\n');
