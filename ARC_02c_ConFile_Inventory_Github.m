%% ARC_02c_ConFile_Inventory.m
% =========================================================================
% SCRIPT 2c: Individual Contrast File Inventory
% =========================================================================
% Purpose:
%   - Locate all individual con_*.nii files across second-level folders
%   - Extract patient ID, session ID, contrast number
%   - Count files per aphasia type and contrast
%   - Identify files available for inferential statistics
%
% Output:
%   - ARC_02c_ConFile_Inventory.mat
%   - ARC_02c_ConFile_Inventory.csv
%   - ARC_02c_Summary.txt
% =========================================================================

%% SETUP
clear; clc;

% Load configuration
config = config_local();

fprintf('=================================================================\n');
fprintf('ARC_02c_ConFile_Inventory.m - Individual Contrast File Inventory\n');
fprintf('=================================================================\n');
fprintf('Start time: %s\n\n', datestr(now));

%% DEFINE PATHS
outputDir = fullfile(config.rootDir, 'ARC_02_Output');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

%% STEP 1: DEFINE SEARCH LOCATIONS
fprintf('-----------------------------------------------------------------\n');
fprintf('STEP 1: Identifying folders to search\n');
fprintf('-----------------------------------------------------------------\n');

% All second-level analysis folders
searchFolders = {
    'second_level_analysis_Anomic'
    'second_level_analysis_Broca'
    'second_level_analysis_Conduction'
    'second_level_analysis_Global'
    'second_level_analysis_Unspecified'
    'second_level_analysis_Wernicke'
    'Conduction aphasia_Group_Analysis_Con_001'
    'Conduction aphasia_Group_Analysis_Con_002'
    'Conduction aphasia_Group_Analysis_Con_003'
};

% Contrast subfolders
contrastFolders = {'contrast_0001', 'contrast_0002', 'contrast_0003'};

fprintf('Folders to search:\n');
for i = 1:length(searchFolders)
    folderPath = fullfile(config.rootDir, searchFolders{i});
    if exist(folderPath, 'dir')
        fprintf('  ✓ %s\n', searchFolders{i});
    else
        fprintf('  ✗ %s (NOT FOUND)\n', searchFolders{i});
    end
end
fprintf('\n');

%% STEP 2: SEARCH FOR INDIVIDUAL CON FILES
fprintf('-----------------------------------------------------------------\n');
fprintf('STEP 2: Searching for individual con files\n');
fprintf('-----------------------------------------------------------------\n');

% Initialize storage
allFiles = {};
allPatientID = {};
allSessionID = {};
allContrastNum = {};
allAphasiaType = {};
allFolderSource = {};
allFilePath = {};

fileCount = 0;

for f = 1:length(searchFolders)
    folderName = searchFolders{f};
    folderPath = fullfile(config.rootDir, folderName);
    
    if ~exist(folderPath, 'dir')
        continue;
    end
    
    % Determine aphasia type from folder name
    if contains(folderName, 'Anomic')
        aphasiaType = 'Anomic';
    elseif contains(folderName, 'Broca')
        aphasiaType = 'Broca';
    elseif contains(folderName, 'Conduction')
        aphasiaType = 'Conduction';
    elseif contains(folderName, 'Global')
        aphasiaType = 'Global';
    elseif contains(folderName, 'Unspecified')
        aphasiaType = 'Unspecified';
    elseif contains(folderName, 'Wernicke')
        aphasiaType = 'Wernicke';
    else
        aphasiaType = 'Unknown';
    end
    
    % Search in contrast subfolders
    for c = 1:length(contrastFolders)
        contrastPath = fullfile(folderPath, contrastFolders{c});
        
        if ~exist(contrastPath, 'dir')
            contrastPath = folderPath;
        end
        
        % Find individual con files (pattern: sub-M*_ses-*_con_*.nii)
        conFiles = dir(fullfile(contrastPath, 'sub-M*_ses-*_con_*.nii'));
        
        for i = 1:length(conFiles)
            fileCount = fileCount + 1;
            
            % Parse filename
            tokens = regexp(conFiles(i).name, 'sub-(\w+)_ses-(\d+)_con_(\d+)\.nii', 'tokens');
            
            if ~isempty(tokens)
                patientID = tokens{1}{1};
                sessionID = tokens{1}{2};
                contrastNum = tokens{1}{3};
            else
                patientID = 'PARSE_ERROR';
                sessionID = 'PARSE_ERROR';
                contrastNum = 'PARSE_ERROR';
            end
            
            % Store
            allFiles{end+1} = conFiles(i).name;
            allPatientID{end+1} = patientID;
            allSessionID{end+1} = sessionID;
            allContrastNum{end+1} = contrastNum;
            allAphasiaType{end+1} = aphasiaType;
            allFolderSource{end+1} = folderName;
            allFilePath{end+1} = fullfile(contrastPath, conFiles(i).name);
        end
    end
end

fprintf('Total individual con files found: %d\n\n', fileCount);

%% STEP 3: CREATE INVENTORY TABLE
fprintf('-----------------------------------------------------------------\n');
fprintf('STEP 3: Creating inventory table\n');
fprintf('-----------------------------------------------------------------\n');

conInventory = table(allPatientID', allSessionID', allContrastNum', ...
                     allAphasiaType', allFolderSource', allFiles', allFilePath', ...
                     'VariableNames', {'PatientID', 'SessionID', 'ContrastNum', ...
                                       'AphasiaType', 'FolderSource', 'FileName', 'FilePath'});

conInventory = sortrows(conInventory, {'AphasiaType', 'PatientID', 'ContrastNum'});

fprintf('Inventory table created: %d rows\n\n', height(conInventory));

%% STEP 4: SUMMARY BY GROUP AND CONTRAST
fprintf('-----------------------------------------------------------------\n');
fprintf('STEP 4: Summary by aphasia type and contrast\n');
fprintf('-----------------------------------------------------------------\n');

aphasiaTypes = unique(conInventory.AphasiaType);
contrastNums = unique(conInventory.ContrastNum);

fprintf('%-15s', 'Aphasia Type');
for c = 1:length(contrastNums)
    fprintf('Con_%s\t', contrastNums{c});
end
fprintf('Total\tUnique_Pts\n');
fprintf('%s\n', repmat('-', 1, 70));

for a = 1:length(aphasiaTypes)
    aType = aphasiaTypes{a};
    fprintf('%-15s', aType);
    
    rowTotal = 0;
    for c = 1:length(contrastNums)
        cNum = contrastNums{c};
        count = sum(strcmp(conInventory.AphasiaType, aType) & ...
                    strcmp(conInventory.ContrastNum, cNum));
        fprintf('%d\t', count);
        rowTotal = rowTotal + count;
    end
    
    groupIdx = strcmp(conInventory.AphasiaType, aType);
    uniquePts = length(unique(conInventory.PatientID(groupIdx)));
    
    fprintf('%d\t%d\n', rowTotal, uniquePts);
end

fprintf('%s\n', repmat('-', 1, 70));

fprintf('%-15s', 'TOTAL');
grandTotal = 0;
for c = 1:length(contrastNums)
    cNum = contrastNums{c};
    count = sum(strcmp(conInventory.ContrastNum, cNum));
    fprintf('%d\t', count);
    grandTotal = grandTotal + count;
end
totalUniquePts = length(unique(conInventory.PatientID));
fprintf('%d\t%d\n\n', grandTotal, totalUniquePts);

%% STEP 5: IDENTIFY CONTRAST 2 (NAMING > ABSTRACT) FOR PRIMARY ANALYSIS
fprintf('-----------------------------------------------------------------\n');
fprintf('STEP 5: Primary contrast (Con_0002 = Naming > Abstract) inventory\n');
fprintf('-----------------------------------------------------------------\n');

primaryContrast = conInventory(strcmp(conInventory.ContrastNum, '0002'), :);

fprintf('Con_0002 (Naming > Abstract) files: %d\n\n', height(primaryContrast));

fprintf('By aphasia type:\n');
for a = 1:length(aphasiaTypes)
    aType = aphasiaTypes{a};
    count = sum(strcmp(primaryContrast.AphasiaType, aType));
    uniquePts = length(unique(primaryContrast.PatientID(strcmp(primaryContrast.AphasiaType, aType))));
    fprintf('  %-15s: %d files, %d unique patients\n', aType, count, uniquePts);
end
fprintf('\n');

%% STEP 6: CHECK FOR ADDITIONAL FILES
fprintf('-----------------------------------------------------------------\n');
fprintf('STEP 6: Checking for additional analysis files\n');
fprintf('-----------------------------------------------------------------\n');

existingICN = dir(fullfile(config.rootDir, 'second_level_analysis_*', '*', 'ICN_Atlas*.xlsx'));
fprintf('Existing ICN_Atlas Excel files: %d\n', length(existingICN));
for i = 1:min(5, length(existingICN))
    fprintf('  %s\n', fullfile(existingICN(i).folder, existingICN(i).name));
end
if length(existingICN) > 5
    fprintf('  ... and %d more\n', length(existingICN)-5);
end
fprintf('\n');

betaFiles = dir(fullfile(config.rootDir, 'second_level_analysis_*', '*', 'beta_*.nii'));
fprintf('Beta files found: %d\n', length(betaFiles));

spmFiles = dir(fullfile(config.rootDir, 'second_level_analysis_*', '*', 'SPM.mat'));
fprintf('SPM.mat files in contrast folders: %d\n\n', length(spmFiles));

%% STEP 7: FLAG FIRST SESSION PER PATIENT
fprintf('-----------------------------------------------------------------\n');
fprintf('STEP 7: Identifying first session per patient\n');
fprintf('-----------------------------------------------------------------\n');

uniquePatients = unique(conInventory.PatientID);
firstSessionFlag = false(height(conInventory), 1);

for p = 1:length(uniquePatients)
    pid = uniquePatients{p};
    patientIdx = find(strcmp(conInventory.PatientID, pid));
    
    if ~isempty(patientIdx)
        sessions = str2double(conInventory.SessionID(patientIdx));
        [~, minIdx] = min(sessions);
        
        firstSession = conInventory.SessionID{patientIdx(minIdx)};
        for idx = patientIdx'
            if strcmp(conInventory.SessionID{idx}, firstSession)
                firstSessionFlag(idx) = true;
            end
        end
    end
end

conInventory.IsFirstSession = firstSessionFlag;

fprintf('Files from first session: %d\n', sum(firstSessionFlag));
fprintf('Files from later sessions: %d\n\n', sum(~firstSessionFlag));

firstSessionPrimary = conInventory(conInventory.IsFirstSession & ...
                                   strcmp(conInventory.ContrastNum, '0002'), :);
fprintf('First-session Con_0002 files by group:\n');
for a = 1:length(aphasiaTypes)
    aType = aphasiaTypes{a};
    count = sum(strcmp(firstSessionPrimary.AphasiaType, aType));
    fprintf('  %-15s: %d patients\n', aType, count);
end
fprintf('  %-15s: %d patients\n\n', 'TOTAL', height(firstSessionPrimary));

%% STEP 8: SAVE OUTPUTS
fprintf('-----------------------------------------------------------------\n');
fprintf('STEP 8: Saving outputs\n');
fprintf('-----------------------------------------------------------------\n');

matFile = fullfile(outputDir, 'ARC_02c_ConFile_Inventory.mat');
save(matFile, 'conInventory', 'primaryContrast', 'firstSessionPrimary');
fprintf('Saved: %s\n', matFile);

csvFile = fullfile(outputDir, 'ARC_02c_ConFile_Inventory.csv');
writetable(conInventory, csvFile);
fprintf('Saved: %s\n', csvFile);

primaryCsvFile = fullfile(outputDir, 'ARC_02c_FirstSession_Con0002.csv');
writetable(firstSessionPrimary, primaryCsvFile);
fprintf('Saved: %s\n', primaryCsvFile);

%% STEP 9: GENERATE SUMMARY REPORT
reportFile = fullfile(outputDir, 'ARC_02c_Summary.txt');
fid = fopen(reportFile, 'w');

fprintf(fid, '=================================================================\n');
fprintf(fid, 'INDIVIDUAL CONTRAST FILE INVENTORY SUMMARY\n');
fprintf(fid, '=================================================================\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));

fprintf(fid, 'TOTAL FILES\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Total con files found: %d\n', height(conInventory));
fprintf(fid, 'Unique patients: %d\n', totalUniquePts);
fprintf(fid, 'Con_0002 (Naming>Abstract) files: %d\n\n', height(primaryContrast));

fprintf(fid, 'FILES BY APHASIA TYPE (Con_0002 only)\n');
fprintf(fid, '-----------------------------------------------------------------\n');
for a = 1:length(aphasiaTypes)
    aType = aphasiaTypes{a};
    count = sum(strcmp(primaryContrast.AphasiaType, aType));
    fprintf(fid, '%-15s: %d\n', aType, count);
end

fprintf(fid, '\nFIRST SESSION ONLY (Con_0002) - FOR INFERENTIAL STATS\n');
fprintf(fid, '-----------------------------------------------------------------\n');
for a = 1:length(aphasiaTypes)
    aType = aphasiaTypes{a};
    count = sum(strcmp(firstSessionPrimary.AphasiaType, aType));
    fprintf(fid, '%-15s: %d\n', aType, count);
end
fprintf(fid, '%-15s: %d\n', 'TOTAL', height(firstSessionPrimary));

fprintf(fid, '\n=================================================================\n');
fprintf(fid, 'END OF REPORT\n');
fprintf(fid, '=================================================================\n');

fclose(fid);
fprintf('Saved: %s\n\n', reportFile);

%% COMPLETION
fprintf('=================================================================\n');
fprintf('ARC_02c_ConFile_Inventory.m - COMPLETE\n');
fprintf('=================================================================\n');
fprintf('End time: %s\n', datestr(now));
fprintf('\nKey outputs:\n');
fprintf('  - Total con files: %d\n', height(conInventory));
fprintf('  - First-session Con_0002 (for analysis): %d patients\n', height(firstSessionPrimary));
fprintf('\nReady for Script 3: Individual-level ICN extraction\n');
fprintf('=================================================================\n');