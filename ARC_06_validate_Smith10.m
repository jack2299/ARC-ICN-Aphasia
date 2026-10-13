%% ARC_06_validate_Smith10.m
% =========================================================================
% Reproducible validation of the Smith10 sensorimotor engagement finding.
%
% This script does NOT re-run the per-subject xATL extraction; those
% inputs are not redistributed. It reads the preserved master tables
% and reproduces the manuscript's reported S06 correlation from
% documented inputs.
%
% Expected output:
%   S06 sensorimotor: rho = +0.318, p = 0.0002, N = 139
%
% Data requirements
% -----------------
% Two master tables are needed and are not committed to this repository
% (participant-level derived data):
%
%   ARC_06_v3_Smith10_Master.mat   (per-subject Smith10 IRi metrics)
%   ARC_03b_v3_Master_Wide.mat     (BM20 primary sample + engagement)
%
% Both are available from the project's Zenodo archive:
%   DOI: 10.5281/zenodo.XXXXXXX
%
% Place them in either of these locations:
%   - alongside this script
%   - in a subfolder named "data" next to this script
%
% Alternatively, create a config_local.m in the same folder as this
% script containing:
%   masterFileOverride = 'full/path/to/ARC_06_v3_Smith10_Master.mat';
%   bm20FileOverride   = 'full/path/to/ARC_03b_v3_Master_Wide.mat';
% config_local.m is git-ignored and never committed.
%
% Author: Jack Kissane
% =========================================================================

clear; clc;

%% ------------------------------------------------------------------------
% PATHS
% -------------------------------------------------------------------------
scriptDir = fileparts(mfilename('fullpath'));

% Optional local override (git-ignored)
if exist('config_local.m', 'file')
    run('config_local.m');
end

% Smith10 master: search order = local override, sibling, data/ subfolder
if exist('masterFileOverride', 'var') && exist(masterFileOverride, 'file') == 2
    masterFile = masterFileOverride;
else
    candidates = {
        fullfile(scriptDir, 'ARC_06_v3_Smith10_Master.mat')
        fullfile(scriptDir, 'data', 'ARC_06_v3_Smith10_Master.mat')
    };
    masterFile = '';
    for i = 1:numel(candidates)
        if exist(candidates{i}, 'file') == 2
            masterFile = candidates{i};
            break;
        end
    end
    assert(~isempty(masterFile), ...
        ['ARC_06_v3_Smith10_Master.mat not found. Looked in:\n  %s\n' ...
         'Download from Zenodo (DOI: 10.5281/zenodo.XXXXXXX) and place ' ...
         'alongside this script or in data/.'], ...
        strjoin(candidates, sprintf('\n  ')));
end

% BM20 master: same pattern
if exist('bm20FileOverride', 'var') && exist(bm20FileOverride, 'file') == 2
    bm20File = bm20FileOverride;
else
    candidatesBM20 = {
        fullfile(scriptDir, 'ARC_03b_v3_Master_Wide.mat')
        fullfile(scriptDir, 'data', 'ARC_03b_v3_Master_Wide.mat')
    };
    bm20File = '';
    for i = 1:numel(candidatesBM20)
        if exist(candidatesBM20{i}, 'file') == 2
            bm20File = candidatesBM20{i};
            break;
        end
    end
    assert(~isempty(bm20File), ...
        ['ARC_03b_v3_Master_Wide.mat not found. Looked in:\n  %s\n' ...
         'Download from Zenodo (DOI: 10.5281/zenodo.XXXXXXX).'], ...
        strjoin(candidatesBM20, sprintf('\n  ')));
end

fprintf('Smith10 master: %s\n', masterFile);
fprintf('BM20 master:    %s\n\n', bm20File);

%% ------------------------------------------------------------------------
% LOAD
% -------------------------------------------------------------------------
load(masterFile, 'Smith10_Wide');
load(bm20File,   'masterWide');
bm20 = masterWide;

%% ------------------------------------------------------------------------
% SAMPLE MASK
% -------------------------------------------------------------------------
% Rest-only exclusion: 36 participants whose first session contained
% resting-state data only. See project documentation for derivation.
restOnly = {'M2088','M2097','M2100','M2101','M2113','M2114','M2117','M2118','M2122','M2126','M2129','M2131','M2135',...
    'M2140','M2141','M2142','M2143','M2144','M2145','M2146','M2149','M2150','M2151','M2152','M2153','M2155','M2156',...
    'M2158','M2159','M2160','M2162','M2164','M2165','M2169','M2184','M2254'};
assert(numel(restOnly) == 36, 'Rest-only list should contain 36 IDs, has %d', numel(restOnly));

isAph   = strcmp(Smith10_Wide.GroupRole, 'Aphasia');
isClean = ~ismember(Smith10_Wide.PatientID, restOnly);
validC2 = ~isnan(bm20.C2_ICN17_IRi);

D = Smith10_Wide(isAph & isClean, :);
D = D(ismember(D.PatientID, bm20.PatientID(validC2)), :);

fprintf('Sample: aphasia ∩ clean ∩ valid-C2   N = %d   (expected 139)\n\n', height(D));

%% ------------------------------------------------------------------------
% PARTIAL SPEARMAN
% -------------------------------------------------------------------------
% Partial Spearman correlation of network IRi with WAB-AQ, controlling
% for age, sex, days post-stroke. Covariate adjustment uses MATLAB's
% partialcorr with 'Type','Spearman' (rank-based).
% -------------------------------------------------------------------------
fprintf('%-5s %-22s %8s %10s %6s\n', 'ICN', 'Label', 'rho', 'p', 'N');
fprintf('%s\n', repmat('-', 1, 60));

labels = {'Medial_Visual','Occipital_Visual','Lateral_Visual','Default_Mode','Cerebellum',...
    'Sensorimotor','Auditory','Executive_Control','Right_Frontoparietal','Left_Frontoparietal'};

results = table();
for n = 1:10
    col = sprintf('S%02d_IRi', n);
    x   = D.(col);  w = D.WAB_AQ;
    age = D.Age_At_Stroke;  sex = D.SexNumeric;  dys = D.Days_Post_Stroke;
    v = ~isnan(x) & ~isnan(w) & ~isnan(age) & ~isnan(sex) & ~isnan(dys);
    [r, p] = partialcorr(x(v), w(v), [age(v), sex(v), dys(v)], ...
        'Type', 'Spearman', 'Rows', 'complete');

    fprintf('S%02d   %-22s %+0.3f   %0.4f   %d\n', n, labels{n}, r, p, sum(v));

    row = table(n, labels(n), r, p, sum(v), ...
        'VariableNames', {'ICN','Label','Partial_rho','p_uncorrected','N'});
    results = [results; row]; %#ok<AGROW>
end

%% ------------------------------------------------------------------------
% SANITY CHECK vs manuscript values
% -------------------------------------------------------------------------
target = results(strcmp(results.Label, 'Sensorimotor'), :);

assert(abs(target.Partial_rho - 0.318) < 0.001, ...
    'Sensorimotor rho = %.4f does not match manuscript value 0.318', target.Partial_rho);
assert(target.N == 139, ...
    'Sensorimotor N = %d, expected 139', target.N);

fprintf('\nVALIDATION PASSED: S06 sensorimotor rho = %+0.3f, N = %d\n', ...
    target.Partial_rho, target.N);

%% ------------------------------------------------------------------------
% SAVE OUTPUT
% -------------------------------------------------------------------------
outFile = fullfile(scriptDir, 'ARC_06_validate_Smith10_output.csv');
writetable(results, outFile);
fprintf('Saved: %s\n', outFile);
