%% MOTION EXTRACTION – CLEAN SAMPLE (PORTABLE VERSION)
% Extracts mean rotation (degrees) and mean translation (mm) from
% realignment parameters for the clean aphasia sample.
%
% REQUIRED DATA FILES (not included in the repository – obtain from ARC):
%   1. ARC_03b_v3_Master_Wide.mat   – post‑IRi‑QC engagement master table
%   2. task_fmri_first_session.csv  – run mapping (from the ARC dataset)
%   3. ds004884‑rp_files/          – folder containing rp_*.txt files
%      (realignment parameters, one file per functional run)
%
% The script saves Motion_Extracted.mat with the motion metrics and the
% filtered cleanAphasia table.
%
% =========================================================================
clear; clc;

%% ---- 0. USER‑CONFIGURABLE PATHS ----
masterWideFile = 'ARC_03b_v3_Master_Wide.mat';   % adjust if needed
runMappingCSV  = 'task_fmri_first_session.csv';   % from ARC
rpFolder       = 'ds004884-rp_files';            % folder with rp_*.txt
outputFile     = 'Motion_Extracted.mat';

%% ---- 1. LOAD DATA AND DEFINE CLEAN SAMPLE ----
load(masterWideFile, 'masterWide');

% Create SexNumeric if not present
if ~ismember('SexNumeric', masterWide.Properties.VariableNames)
    SexNumeric = NaN(height(masterWide), 1);
    SexNumeric(strcmp(masterWide.Sex, 'F')) = 0;
    SexNumeric(strcmp(masterWide.Sex, 'M')) = 1;
    masterWide.SexNumeric = SexNumeric;
end

% Rest‑only aphasia participants (IDs from ARC dataset; anonymized)
restOnlyList = {'M2088','M2097','M2100','M2101','M2113','M2114',...
    'M2117','M2118','M2122','M2126','M2129','M2131','M2135',...
    'M2140','M2141','M2142','M2143','M2144','M2145','M2146',...
    'M2149','M2150','M2151','M2152','M2153','M2155','M2156',...
    'M2158','M2159','M2160','M2162','M2164','M2165','M2169',...
    'M2184','M2254'};

isRestOnly = ismember(masterWide.PatientID, restOnlyList);
isAphasia  = strcmp(masterWide.GroupRole, 'Aphasia');
hasValidC2 = ~isnan(masterWide.C2_ICN17_IRi);
cleanAphasia = masterWide(isAphasia & hasValidC2 & ~isRestOnly, :);
fprintf('Clean aphasia sample: N = %d\n', height(cleanAphasia));

%% ---- 2. PARSE RUN NUMBERS FROM CSV ----
fid = fopen(runMappingCSV, 'r');
lines = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
fclose(fid);
lines = lines{1};

taskParticipants = {};
taskRuns = [];
for i = 2:length(lines)
    line = lines{i};
    subToken = regexp(line, 'sub-(\w+)', 'tokens');
    if ~isempty(subToken)
        taskParticipants{end+1} = subToken{1}{1};
    end
    runToken = regexp(line, 'run-(\d+)', 'tokens');
    if ~isempty(runToken)
        taskRuns(end+1, 1) = str2double(runToken{1}{1});
    end
end

runLookup = containers.Map();
for i = 1:length(taskParticipants)
    if ~isempty(taskParticipants{i}) && ~isnan(taskRuns(i))
        runLookup(taskParticipants{i}) = taskRuns(i);
    end
end

%% ---- 3. FIND RP FILES AND COMPUTE MOTION METRICS ----
rpFiles = dir(fullfile(rpFolder, 'rp_*.txt'));

motion_rot_valid   = NaN(height(cleanAphasia), 1);
motion_trans_valid = NaN(height(cleanAphasia), 1);

fprintf('Matching participants to motion files...\n');
for i = 1:height(cleanAphasia)
    pid = cleanAphasia.PatientID{i};
    
    if ~isKey(runLookup, pid)
        fprintf('  Warning: No run number for %s\n', pid);
        continue;
    end
    targetRun = runLookup(pid);
    
    pattern = sprintf('sub-%s.*run-%d_bold', pid, targetRun);
    matchingFiles = rpFiles(~cellfun(@isempty, regexp({rpFiles.name}, pattern)));
    
    if isempty(matchingFiles)
        fprintf('  Warning: No rp file for %s (run %d)\n', pid, targetRun);
        continue;
    end
    
    rpData = load(fullfile(rpFolder, matchingFiles(1).name));
    
    % Mean rotation (degrees) – cols 4‑6, radians -> degrees
    mean_rot_rad = mean(abs(rpData(:, 4:6)), 'all');
    mean_rot_deg = mean_rot_rad * (180 / pi);
    
    % Mean translation (mm) – cols 1‑3
    mean_trans = mean(abs(rpData(:, 1:3)), 'all');
    
    motion_rot_valid(i)   = mean_rot_deg;
    motion_trans_valid(i) = mean_trans;
end

validMotion = ~isnan(motion_rot_valid) & ~isnan(motion_trans_valid);
fprintf('\nSuccessfully matched %d / %d participants\n', sum(validMotion), height(cleanAphasia));

motion_rot_valid   = motion_rot_valid(validMotion);
motion_trans_valid = motion_trans_valid(validMotion);
cleanAphasia       = cleanAphasia(validMotion, :);
fprintf('Final clean sample with motion data: N = %d\n', height(cleanAphasia));

%% ---- 4. SAVE ----
save(outputFile, 'motion_rot_valid', 'motion_trans_valid', 'cleanAphasia');
fprintf('Saved: %s\n', outputFile);
fprintf('Now run the motion post‑hoc script.\n');