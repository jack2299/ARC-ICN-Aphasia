%% ARC_06_validate_Smith10.m
% =========================================================================
% Reproducible validation of the Smith10 sensorimotor engagement finding.
%
% This script does NOT re-run the per-subject xATL extraction (those
% inputs are not redistributed). It reads the preserved master table
% produced by ARC_06_v3_Smith10_Robustness.m and reproduces the
% manuscript's reported S06 correlation from documented inputs.
%
% Expected output:
%   S06 sensorimotor: rho = +0.318, p = 0.0002, N = 139
%
% Author: Jack Kissane
% =========================================================================

clear; clc;

rootDir = fileparts(mfilename('fullpath'));
% Adjust if the data folder sits elsewhere:
masterFile = fullfile(rootDir, 'data', 'ARC_06_v3_Smith10_Master.mat');
bm20File   = fullfile(rootDir, 'data', 'ARC_03b_v3_Master_Wide.mat');

assert(exist(masterFile, 'file') == 2, 'Smith10 master not found: %s', masterFile);
assert(exist(bm20File,   'file') == 2, 'BM20 master not found: %s',   bm20File);

load(masterFile, 'Smith10_Wide');
load(bm20File,   'masterWide');
bm20 = masterWide;

% Rest-only exclusion (36 IDs; see Handoff Section 3.2)
restOnly = {'M2088','M2097','M2100','M2101','M2113','M2114','M2117','M2118','M2122','M2126','M2129','M2131','M2135',...
    'M2140','M2141','M2142','M2143','M2144','M2145','M2146','M2149','M2150','M2151','M2152','M2153','M2155','M2156',...
    'M2158','M2159','M2160','M2162','M2164','M2165','M2169','M2184','M2254'};

isAph   = strcmp(Smith10_Wide.GroupRole, 'Aphasia');
isClean = ~ismember(Smith10_Wide.PatientID, restOnly);
validC2 = ~isnan(bm20.C2_ICN17_IRi);

D = Smith10_Wide(isAph & isClean, :);
D = D(ismember(D.PatientID, bm20.PatientID(validC2)), :);

fprintf('Sample: aphasia ∩ clean ∩ valid-C2  N = %d  (expected 139)\n\n', height(D));

% Partial Spearman, matching ARC_06_v3 exactly
fprintf('%-5s %10s %8s %10s %8s\n', 'ICN', 'Label', 'rho', 'p', 'N');
fprintf('%s\n', repmat('-', 1, 70));

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

    fprintf('S%02d   %-18s %+0.3f   %0.4f   %d\n', n, labels{n}, r, p, sum(v));

    row = table(n, labels(n), r, p, sum(v), ...
        'VariableNames', {'ICN','Label','Partial_rho','p_uncorrected','N'});
    results = [results; row]; %#ok<AGROW>
end

% Sanity check against manuscript values
target = results(strcmp(results.Label, 'Sensorimotor'), :);
assert(abs(target.Partial_rho - 0.318) < 0.005, ...
    'Sensorimotor rho = %.3f does not match manuscript value 0.318', target.Partial_rho);
assert(target.N == 139, 'Sensorimotor N = %d, expected 139', target.N);

fprintf('\nVALIDATION PASSED: S06 sensorimotor rho = %+0.3f, N = %d\n', ...
    target.Partial_rho, target.N);

writetable(results, fullfile(rootDir, 'data', 'ARC_06_validate_Smith10_output.csv'));