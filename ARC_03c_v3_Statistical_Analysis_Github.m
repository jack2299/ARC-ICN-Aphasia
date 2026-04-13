%% ARC_03c_v3_Statistical_Analysis.m
% =========================================================================
% SCRIPT 3c v3: Statistical Analysis
% =========================================================================
% Purpose:
%   - Level 1: Disease Effect (Aphasia vs Stroke Control) - ANCOVA
%   - Level 2: Subtype Effect (4-group ANCOVA) + BAYESIAN ANALYSIS
%   - Level 3: Severity Effect (WAB-AQ partial correlations)
%   - Level 4: Contrast Comparison (C2 vs C1 predictive value)
%   - Level 5: **Removed* Network Specificity (language vs non-language)
%   - Level 6: **Removed: Circular** Predictive Modeling (LASSO regression)
%   - Level 10b: BOOTSTRAP MEDIATION (Lesion Volume → ICN Engagement → Severity)
%
% V3 CHANGES FROM V2:
%   - Input from ARC_03b_v3_Output (IRi QC-cleaned data)
%   - Output filenames updated to v3
%   - NEW: Bayesian subtype analysis (BIC-based Bayes Factor) after Level 2
%   - NEW: Bootstrap mediation with BCa CIs after Level 6 (replaces Level 10)
%   - NEW: IRi threshold sensitivity analysis
%   - Mediation reporting: Unstandardized indirect effects (not %)
%   - logVari columns saved to master tables
%
% METRICS:
%   - IRi: ICNiRelativeInvolvement (PRIMARY metric)
%   - MANi: NormalisedMeanICNiActivation (SECONDARY)
%   - logVari: log(var(ZActiveVoxels)+1) (EXPLORATORY)
%
% A PRIORI LANGUAGE NETWORKS (BrainMap20):
%   ICN04: Bilateral anterior insula/frontal opercula (articulation/production)
%   ICN16: Transverse temporal gyri (auditory/comprehension)
%   ICN18: Left-lateralized fronto-parietal (classic language network)
%
% Input:
%   - ARC_03b_v3_Output/ARC_03b_v3_Master_Wide.mat
%   - ARC_03b_v3_Output/ARC_03b_v3_Master_Long.mat
%
% Output:
%   - ARC_03c_v3_Output/ARC_03c_v3_Level[1-6]_*.csv
%   - ARC_03c_v3_Output/ARC_03c_v3_Level2_Bayesian.csv (NEW)
%   - ARC_03c_v3_Output/ARC_03c_v3_Bootstrap_Mediation.mat (NEW)
%   - ARC_03c_v3_Output/ARC_03c_v3_IRi_Sensitivity.csv (NEW)
%   - ARC_03c_v3_Output/ARC_03c_v3_AllResults.mat
%   - ARC_03c_v3_Output/ARC_03c_v3_Summary.txt
% =========================================================================

%% SETUP
% =========================================================================
clear; clc;

% Load configuration
config = config_local();

fprintf('=================================================================\n');
fprintf('ARC_03c_v3_Statistical_Analysis.m\n');
fprintf('Statistical Analysis with IRi QC-Cleaned Data\n');
fprintf('=================================================================\n');
fprintf('Start time: %s\n\n', datestr(now));

fprintf('*** V3 CHANGES ***\n');
fprintf('- Input: ARC_03b_v3_Output (IRi QC-cleaned)\n');
fprintf('- NEW: Bayesian subtype analysis (BIC-based BF)\n');
fprintf('- NEW: Bootstrap mediation with BCa CIs\n');
fprintf('- NEW: IRi threshold sensitivity analysis\n');
fprintf('- Mediation: Unstandardized indirect effects (not percentages)\n\n');

%% DEFINE PATHS
% =========================================================================

% Root directory - NOW LOADED FROM CONFIG
rootDir = config.rootDir;

% V3: Input from 3b_v3 output
inputDir = fullfile(rootDir, 'ARC_03_Output', 'ARC_03b_v3_Output');
wideFile = fullfile(inputDir, 'ARC_03b_v3_Master_Wide.mat');
longFile = fullfile(inputDir, 'ARC_03b_v3_Master_Long.mat');

% V3: Output directory
outputDir = fullfile(rootDir, 'ARC_03_Output', 'ARC_03c_v3_Output');

% Verify input files exist
if ~exist(wideFile, 'file')
    error('Wide file not found: %s\nRun Script 3b_v3 first.', wideFile);
end
if ~exist(longFile, 'file')
    error('Long file not found: %s\nRun Script 3b_v3 first.', longFile);
end

% Create output directory if needed
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('Created output directory: %s\n\n', outputDir);
else
    fprintf('Output directory: %s\n\n', outputDir);
end

%% LOAD DATA
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Loading data from Script 3b v3 (IRi QC-cleaned)\n');
fprintf('-----------------------------------------------------------------\n');

load(wideFile, 'masterWide', 'trueICN_labels', 'nICN', 'IRi_threshold');
load(longFile, 'masterLong');

fprintf('Wide format: %d patients × %d columns\n', height(masterWide), width(masterWide));
fprintf('Long format: %d rows × %d columns\n', height(masterLong), width(masterLong));
fprintf('ICNs: %d\n', nICN);
fprintf('IRi QC threshold used: %.1f\n\n', IRi_threshold);

%% DEFINE ANALYSIS PARAMETERS
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Analysis parameters\n');
fprintf('-----------------------------------------------------------------\n');

% FDR threshold
fdrQ = 0.05;

% Metrics to analyze
metrics = {'IRi', 'MANi', 'logVari'};
metricLabels = {'Relative Spatial Involvement', 'Normalized Mean Activation', 'Log Activation Variance'};
primaryMetric = 'IRi';

% Contrasts
contrast1 = 'C1';  % Task > Rest
contrast2 = 'C2';  % Naming > Abstract (PRIMARY)
primaryContrast = 'C2';

% A PRIORI LANGUAGE NETWORKS
languageICNs = [4, 16, 18];
languageLabels = {
    'ICN04: Anterior Insula/Frontal Opercula (articulation)',
    'ICN16: Transverse Temporal Gyri (auditory)',
    'ICN18: Left Fronto-Parietal (language)'
};

% Potential compensatory network
compensatoryICNs = [15];  % Right fronto-parietal

% Non-language networks
nonLanguageICNs = setdiff(1:nICN, [languageICNs, compensatoryICNs]);

% Bootstrap parameters (for mediation)
nBootstrap = 5000;
bootAlpha = 0.05;

% Bayesian parameters
BF_threshold_moderate = 3;  % BF > 3 = moderate evidence

fprintf('Primary metric: %s\n', primaryMetric);
fprintf('Primary contrast: %s (Naming > Abstract)\n', primaryContrast);
fprintf('FDR threshold: q = %.2f\n', fdrQ);
fprintf('Covariates: Age_At_Stroke, Sex, Days_Post_Stroke\n');
fprintf('Bootstrap iterations: %d\n\n', nBootstrap);

fprintf('A PRIORI LANGUAGE NETWORKS:\n');
for i = 1:length(languageICNs)
    fprintf('  %s\n', languageLabels{i});
end
fprintf('\n');

%% DEFINE GROUP ASSIGNMENTS
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Defining groups\n');
fprintf('-----------------------------------------------------------------\n');

% Aphasia vs Control
isControl = strcmp(masterWide.GroupRole, 'Stroke Control (No Aphasia)');
isAphasia = strcmp(masterWide.GroupRole, 'Aphasia');

nControl = sum(isControl);
nAphasia = sum(isAphasia);

fprintf('Stroke Controls (No Aphasia): N = %d\n', nControl);
fprintf('Aphasia (all types): N = %d\n', nAphasia);

% Aphasia subtypes
isAnomic = strcmp(masterWide.AphasiaType, 'Anomic');
isBroca = strcmp(masterWide.AphasiaType, 'Broca');
isConduction = strcmp(masterWide.AphasiaType, 'Conduction');
isGlobal = strcmp(masterWide.AphasiaType, 'Global');
isWernicke = strcmp(masterWide.AphasiaType, 'Wernicke');

% For subtype analysis: 4 groups (exclude Wernicke N=4)
isSubtypeAnalysis = isAnomic | isBroca | isConduction | isGlobal;

fprintf('\nSubtype breakdown:\n');
fprintf('  Anomic: N = %d\n', sum(isAnomic));
fprintf('  Broca: N = %d\n', sum(isBroca));
fprintf('  Conduction: N = %d\n', sum(isConduction));
fprintf('  Global: N = %d\n', sum(isGlobal));
fprintf('  Wernicke: N = %d (excluded - too small)\n', sum(isWernicke));
fprintf('  Total for subtype analysis: N = %d\n\n', sum(isSubtypeAnalysis));

%% PREPARE COVARIATES
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Preparing covariates\n');
fprintf('-----------------------------------------------------------------\n');

% Convert Sex to numeric (0=F, 1=M) if not already present
if ~ismember('SexNumeric', masterWide.Properties.VariableNames)
    SexNumeric = NaN(height(masterWide), 1);
    SexNumeric(strcmp(masterWide.Sex, 'F')) = 0;
    SexNumeric(strcmp(masterWide.Sex, 'M')) = 1;
    masterWide.SexNumeric = SexNumeric;
end

% Check covariate availability
hasAge = ~isnan(masterWide.Age_At_Stroke);
hasDays = ~isnan(masterWide.Days_Post_Stroke);
hasSex = ~isnan(masterWide.SexNumeric);

fprintf('Age_At_Stroke: %d/%d (%.1f%%)\n', sum(hasAge), height(masterWide), 100*sum(hasAge)/height(masterWide));
fprintf('Days_Post_Stroke: %d/%d (%.1f%%)\n', sum(hasDays), height(masterWide), 100*sum(hasDays)/height(masterWide));
fprintf('Sex: %d/%d (%.1f%%)\n\n', sum(hasSex), height(masterWide), 100*sum(hasSex)/height(masterWide));

%% LOG-TRANSFORM VARI AND SAVE TO MASTER TABLE
% =========================================================================
fprintf('-----------------------------------------------------------------\n');
fprintf('Log-transforming Vari and adding to master table\n');
fprintf('-----------------------------------------------------------------\n');

% Add log-transformed Vari columns (3b saves Vari but not logVari)
for n = 1:nICN
    % C1
    colC1 = sprintf('C1_ICN%02d_Vari', n);
    colC1_log = sprintf('C1_ICN%02d_logVari', n);
    if ismember(colC1, masterWide.Properties.VariableNames) && ...
       ~ismember(colC1_log, masterWide.Properties.VariableNames)
        masterWide.(colC1_log) = log(masterWide.(colC1) + 1);
    end
    
    % C2
    colC2 = sprintf('C2_ICN%02d_Vari', n);
    colC2_log = sprintf('C2_ICN%02d_logVari', n);
    if ismember(colC2, masterWide.Properties.VariableNames) && ...
       ~ismember(colC2_log, masterWide.Properties.VariableNames)
        masterWide.(colC2_log) = log(masterWide.(colC2) + 1);
    end
end

% Check transformation on sample
sampleVari = masterWide.C2_ICN01_Vari;
sampleLogVari = masterWide.C2_ICN01_logVari;
fprintf('ICN01 Vari:    Median=%.2f, Mean=%.2f, Skewness=%.2f\n', ...
        median(sampleVari, 'omitnan'), mean(sampleVari, 'omitnan'), skewness(sampleVari));
fprintf('ICN01 logVari: Median=%.2f, Mean=%.2f, Skewness=%.2f\n\n', ...
        median(sampleLogVari, 'omitnan'), mean(sampleLogVari, 'omitnan'), skewness(sampleLogVari));

%% HELPER FUNCTIONS
% =========================================================================

% FDR correction (Benjamini-Hochberg)
    function [h, crit_p, adj_p] = fdr_bh(pvals, q)
        if nargin < 2, q = 0.05; end
        
        pvals = pvals(:);
        m = length(pvals);
        
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
        for i = 1:m
            if ~isnan(sorted_p(i))
                adj_p(sort_idx(i)) = sorted_p(i) * m / i;
            end
        end
        
        for i = m-1:-1:1
            if ~isnan(adj_p(sort_idx(i))) && i < m && ~isnan(adj_p(sort_idx(i+1)))
                adj_p(sort_idx(i)) = min(adj_p(sort_idx(i)), adj_p(sort_idx(i+1)));
            end
        end
        adj_p = min(adj_p, 1);
        
        h = pvals <= crit_p & ~isnan(pvals);
    end

% Cohen's d
    function d = cohens_d(x1, x2)
        x1 = x1(~isnan(x1));
        x2 = x2(~isnan(x2));
        n1 = length(x1);
        n2 = length(x2);
        if n1 < 2 || n2 < 2
            d = NaN;
            return;
        end
        pooled_std = sqrt(((n1-1)*var(x1) + (n2-1)*var(x2)) / (n1+n2-2));
        if pooled_std == 0
            d = 0;
        else
            d = (mean(x1) - mean(x2)) / pooled_std;
        end
    end

% Partial Spearman correlation
    function [rho, pval, n_valid] = partial_corr_spearman(x, y, covariates)
        validIdx = ~isnan(x) & ~isnan(y) & ~any(isnan(covariates), 2);
        x = x(validIdx);
        y = y(validIdx);
        covariates = covariates(validIdx, :);
        n_valid = length(x);
        
        if n_valid < 10
            rho = NaN;
            pval = NaN;
            return;
        end
        
        x_rank = tiedrank(x);
        y_rank = tiedrank(y);
        
        X_cov = [ones(length(x), 1), covariates];
        
        beta_x = X_cov \ x_rank;
        resid_x = x_rank - X_cov * beta_x;
        
        beta_y = X_cov \ y_rank;
        resid_y = y_rank - X_cov * beta_y;
        
        [rho, pval] = corr(resid_x, resid_y, 'Type', 'Pearson');
    end

% BIC-based Bayes Factor approximation
    function [BF10, BF01, interpretation] = bic_bayes_factor(BIC_null, BIC_alt)
        % BF10 = evidence for alternative over null
        % BF01 = evidence for null over alternative
        deltaBIC = BIC_null - BIC_alt;
        BF10 = exp(deltaBIC / 2);
        BF01 = 1 / BF10;
        
        if BF01 > 10
            interpretation = 'Strong evidence FOR null (equivalence)';
        elseif BF01 > 3
            interpretation = 'Moderate evidence FOR null (equivalence)';
        elseif BF01 > 1
            interpretation = 'Weak evidence for null';
        elseif BF10 > 1 && BF10 <= 3
            interpretation = 'Weak evidence for alternative (difference)';
        elseif BF10 > 3 && BF10 <= 10
            interpretation = 'Moderate evidence FOR alternative (difference)';
        elseif BF10 > 10
            interpretation = 'Strong evidence FOR alternative (difference)';
        else
            interpretation = 'Inconclusive';
        end
    end

%% =========================================================================
%% LEVEL 1: DISEASE EFFECT (APHASIA VS CONTROL)
%% =========================================================================
fprintf('=================================================================\n');
fprintf('LEVEL 1: Disease Effect (Aphasia vs Stroke Control)\n');
fprintf('         ANCOVA controlling for Age, Sex, Days_Post_Stroke\n');
fprintf('=================================================================\n');
fprintf('HYPOTHESIS H1a: Aphasia shows reduced IRi in language networks\n');
fprintf('HYPOTHESIS H1b: Aphasia may show increased IRi in compensatory (R FP)\n\n');

% Storage
nTests_L1 = length(metrics) * nICN;
Level1_Results = table('Size', [nTests_L1, 20], ...
    'VariableTypes', {'double', 'cell', 'cell', 'logical', 'logical', ...
                      'double', 'double', 'double', 'double', 'double', 'double', ...
                      'double', 'double', 'double', 'double', 'double', ...
                      'double', 'double', 'logical', 'double'}, ...
    'VariableNames', {'ICN', 'ICN_Label', 'Metric', 'IsLanguageICN', 'IsCompensatoryICN', ...
                      'Aphasia_Mean', 'Aphasia_SD', 'Aphasia_N', ...
                      'Control_Mean', 'Control_SD', 'Control_N', ...
                      'F_stat', 'df1', 'df2', 'p_uncorrected', 'Cohens_d', ...
                      'Cohens_d_CI_lo', 'Cohens_d_CI_hi', 'Significant_FDR', 'p_FDR'});

resultIdx = 0;
for m = 1:length(metrics)
    metric = metrics{m};
    
    for n = 1:nICN
        resultIdx = resultIdx + 1;
        
        % Get column name
        colName = sprintf('%s_ICN%02d_%s', primaryContrast, n, metric);
        
        % Extract data
        y = masterWide.(colName);
        group = double(isAphasia);
        age = masterWide.Age_At_Stroke;
        sex = masterWide.SexNumeric;
        days = masterWide.Days_Post_Stroke;
        
        % Valid observations
        validIdx = ~isnan(y) & ~isnan(age) & ~isnan(sex) & ~isnan(days);
        
        y_valid = y(validIdx);
        group_valid = group(validIdx);
        age_valid = age(validIdx);
        sex_valid = sex(validIdx);
        days_valid = days(validIdx);
        
        aphasiaData = y_valid(group_valid == 1);
        controlData = y_valid(group_valid == 0);
        
        % ANCOVA via regression
        X = [ones(length(y_valid), 1), group_valid, age_valid, sex_valid, days_valid];
        X_reduced = [ones(length(y_valid), 1), age_valid, sex_valid, days_valid];
        
        [~, ~, ~, ~, stats_full] = regress(y_valid, X);
        [~, ~, ~, ~, stats_reduced] = regress(y_valid, X_reduced);
        
        SSE_full = stats_full(4) * (length(y_valid) - size(X, 2));
        SSE_reduced = stats_reduced(4) * (length(y_valid) - size(X_reduced, 2));
        SS_group = SSE_reduced - SSE_full;
        df1 = 1;
        df2 = length(y_valid) - size(X, 2);
        MSE_full = SSE_full / df2;
        
        if MSE_full > 0
            F_stat = (SS_group / df1) / MSE_full;
            p_value = 1 - fcdf(F_stat, df1, df2);
        else
            F_stat = 0;
            p_value = 1;
        end
        
        d = cohens_d(aphasiaData, controlData);
        
        % Cohen's d confidence interval (approximate)
        n1 = length(aphasiaData);
        n2 = length(controlData);
        se_d = sqrt((n1+n2)/(n1*n2) + d^2/(2*(n1+n2)));
        d_ci_lo = d - 1.96*se_d;
        d_ci_hi = d + 1.96*se_d;
        
        % Store
        Level1_Results.ICN(resultIdx) = n;
        Level1_Results.ICN_Label{resultIdx} = trueICN_labels{n};
        Level1_Results.Metric{resultIdx} = metric;
        Level1_Results.IsLanguageICN(resultIdx) = ismember(n, languageICNs);
        Level1_Results.IsCompensatoryICN(resultIdx) = ismember(n, compensatoryICNs);
        Level1_Results.Aphasia_Mean(resultIdx) = mean(aphasiaData);
        Level1_Results.Aphasia_SD(resultIdx) = std(aphasiaData);
        Level1_Results.Aphasia_N(resultIdx) = length(aphasiaData);
        Level1_Results.Control_Mean(resultIdx) = mean(controlData);
        Level1_Results.Control_SD(resultIdx) = std(controlData);
        Level1_Results.Control_N(resultIdx) = length(controlData);
        Level1_Results.F_stat(resultIdx) = F_stat;
        Level1_Results.df1(resultIdx) = df1;
        Level1_Results.df2(resultIdx) = df2;
        Level1_Results.p_uncorrected(resultIdx) = p_value;
        Level1_Results.Cohens_d(resultIdx) = d;
        Level1_Results.Cohens_d_CI_lo(resultIdx) = d_ci_lo;
        Level1_Results.Cohens_d_CI_hi(resultIdx) = d_ci_hi;
    end
end

% FDR correction
[Level1_Results.Significant_FDR, ~, Level1_Results.p_FDR] = fdr_bh(Level1_Results.p_uncorrected, fdrQ);

% Sort
Level1_Results = sortrows(Level1_Results, 'p_uncorrected');

% Report
sigL1 = Level1_Results(Level1_Results.Significant_FDR, :);
fprintf('Significant results (FDR q < %.2f): %d / %d tests\n\n', fdrQ, height(sigL1), nTests_L1);

% Hypothesis-specific results
fprintf('=== HYPOTHESIS TESTS ===\n\n');

% H1a: Language networks
L1_lang = Level1_Results(Level1_Results.IsLanguageICN & strcmp(Level1_Results.Metric, primaryMetric), :);
fprintf('H1a: Language Network IRi (Aphasia < Control?)\n');
fprintf('%-6s %-40s %8s %8s %12s %10s\n', 'ICN', 'Label', 'Aph_Mean', 'Con_Mean', 'd [95%CI]', 'p_unc');
fprintf('%s\n', repmat('-', 1, 95));
for i = 1:height(L1_lang)
    fprintf('ICN%02d  %-40s %8.4f %8.4f  %5.2f [%5.2f,%5.2f] %10.4f\n', ...
            L1_lang.ICN(i), L1_lang.ICN_Label{i}(1:min(40,end)), ...
            L1_lang.Aphasia_Mean(i), L1_lang.Control_Mean(i), ...
            L1_lang.Cohens_d(i), L1_lang.Cohens_d_CI_lo(i), L1_lang.Cohens_d_CI_hi(i), ...
            L1_lang.p_uncorrected(i));
end
fprintf('\n');

% Overall summary
meanAbsD_IRi = mean(abs(Level1_Results.Cohens_d(strcmp(Level1_Results.Metric, 'IRi'))));
fprintf('Overall mean |Cohen''s d| for IRi: %.3f\n', meanAbsD_IRi);

if height(sigL1) == 0
    fprintf('\n*** INTERPRETATION: No significant disease effects after FDR correction.\n');
    fprintf('    Effect sizes are small (mean |d|=%.2f). This may indicate:\n', meanAbsD_IRi);
    fprintf('    - ICN engagement is preserved in aphasia during naming tasks\n');
    fprintf('    - Deficits arise from disconnection rather than altered recruitment\n');
    fprintf('    - Or the BrainMap20 atlas does not optimally capture language regions\n');
end
fprintf('\n');

%% =========================================================================
%% LEVEL 2: SUBTYPE EFFECT (4-GROUP ANCOVA) + BAYESIAN ANALYSIS
%% =========================================================================
fprintf('=================================================================\n');
fprintf('LEVEL 2: Subtype Effect (4-Group ANCOVA) + BAYESIAN ANALYSIS\n');
fprintf('         Anomic vs Broca vs Conduction vs Global\n');
fprintf('=================================================================\n');
fprintf('HYPOTHESIS H3a: Distinct IRi profiles across subtypes\n');
fprintf('HYPOTHESIS H3b: Global/Broca show higher Vari (disorganized)\n');
fprintf('NEW: Bayesian BIC-based analysis for small group comparisons\n\n');

subtypeData = masterWide(isSubtypeAnalysis, :);
subtypeLabels = subtypeData.AphasiaType;

fprintf('Patients in analysis: N = %d\n', height(subtypeData));
fprintf('  Anomic: N = %d\n', sum(strcmp(subtypeLabels, 'Anomic')));
fprintf('  Broca: N = %d\n', sum(strcmp(subtypeLabels, 'Broca')));
fprintf('  Conduction: N = %d\n', sum(strcmp(subtypeLabels, 'Conduction')));
fprintf('  Global: N = %d (small sample - Bayesian analysis warranted)\n\n', sum(strcmp(subtypeLabels, 'Global')));

% Storage for frequentist results
nTests_L2 = length(metrics) * nICN;
Level2_Results = table('Size', [nTests_L2, 17], ...
    'VariableTypes', {'double', 'cell', 'cell', 'logical', ...
                      'double', 'double', 'double', 'double', 'double', ...
                      'double', 'double', 'double', 'double', ...
                      'double', 'double', 'logical', 'double'}, ...
    'VariableNames', {'ICN', 'ICN_Label', 'Metric', 'IsLanguageICN', ...
                      'F_stat', 'df_between', 'df_within', 'p_uncorrected', 'Eta_squared', ...
                      'Mean_Anomic', 'Mean_Broca', 'Mean_Conduction', 'Mean_Global', ...
                      'N_valid', 'MaxDiff', 'Significant_FDR', 'p_FDR'});

resultIdx = 0;
for m = 1:length(metrics)
    metric = metrics{m};
    
    for n = 1:nICN
        resultIdx = resultIdx + 1;
        
        colName = sprintf('%s_ICN%02d_%s', primaryContrast, n, metric);
        
        y = subtypeData.(colName);
        age = subtypeData.Age_At_Stroke;
        sex = subtypeData.SexNumeric;
        days = subtypeData.Days_Post_Stroke;
        groups = subtypeData.AphasiaType;
        
        validIdx = ~isnan(y) & ~isnan(age) & ~isnan(sex) & ~isnan(days);
        y_valid = y(validIdx);
        age_valid = age(validIdx);
        sex_valid = sex(validIdx);
        days_valid = days(validIdx);
        groups_valid = groups(validIdx);
        
        % Dummy variables (Anomic = reference)
        dummyBroca = double(strcmp(groups_valid, 'Broca'));
        dummyConduction = double(strcmp(groups_valid, 'Conduction'));
        dummyGlobal = double(strcmp(groups_valid, 'Global'));
        
        X_full = [ones(length(y_valid), 1), dummyBroca, dummyConduction, dummyGlobal, ...
                  age_valid, sex_valid, days_valid];
        X_reduced = [ones(length(y_valid), 1), age_valid, sex_valid, days_valid];
        
        [~, ~, ~, ~, stats_full] = regress(y_valid, X_full);
        [~, ~, ~, ~, stats_reduced] = regress(y_valid, X_reduced);
        
        SSE_full = stats_full(4) * (length(y_valid) - size(X_full, 2));
        SSE_reduced = stats_reduced(4) * (length(y_valid) - size(X_reduced, 2));
        SS_group = SSE_reduced - SSE_full;
        df_between = 3;
        df_within = length(y_valid) - size(X_full, 2);
        MSE_full = SSE_full / df_within;
        
        if MSE_full > 0
            F_stat = (SS_group / df_between) / MSE_full;
            p_value = 1 - fcdf(F_stat, df_between, df_within);
        else
            F_stat = 0;
            p_value = 1;
        end
        
        SS_total = var(y_valid) * (length(y_valid) - 1);
        eta_sq = max(0, SS_group / SS_total);
        
        % Group means
        meanAnomic = mean(y_valid(strcmp(groups_valid, 'Anomic')), 'omitnan');
        meanBroca = mean(y_valid(strcmp(groups_valid, 'Broca')), 'omitnan');
        meanConduction = mean(y_valid(strcmp(groups_valid, 'Conduction')), 'omitnan');
        meanGlobal = mean(y_valid(strcmp(groups_valid, 'Global')), 'omitnan');
        
        groupMeans = [meanAnomic, meanBroca, meanConduction, meanGlobal];
        maxDiff = max(groupMeans) - min(groupMeans);
        
        % Store
        Level2_Results.ICN(resultIdx) = n;
        Level2_Results.ICN_Label{resultIdx} = trueICN_labels{n};
        Level2_Results.Metric{resultIdx} = metric;
        Level2_Results.IsLanguageICN(resultIdx) = ismember(n, languageICNs);
        Level2_Results.F_stat(resultIdx) = F_stat;
        Level2_Results.df_between(resultIdx) = df_between;
        Level2_Results.df_within(resultIdx) = df_within;
        Level2_Results.p_uncorrected(resultIdx) = p_value;
        Level2_Results.Eta_squared(resultIdx) = eta_sq;
        Level2_Results.Mean_Anomic(resultIdx) = meanAnomic;
        Level2_Results.Mean_Broca(resultIdx) = meanBroca;
        Level2_Results.Mean_Conduction(resultIdx) = meanConduction;
        Level2_Results.Mean_Global(resultIdx) = meanGlobal;
        Level2_Results.N_valid(resultIdx) = length(y_valid);
        Level2_Results.MaxDiff(resultIdx) = maxDiff;
    end
end

% FDR correction
[Level2_Results.Significant_FDR, ~, Level2_Results.p_FDR] = fdr_bh(Level2_Results.p_uncorrected, fdrQ);

Level2_Results = sortrows(Level2_Results, 'p_uncorrected');

sigL2 = Level2_Results(Level2_Results.Significant_FDR, :);
fprintf('Frequentist: Significant results (FDR q < %.2f): %d / %d tests\n\n', fdrQ, height(sigL2), nTests_L2);

%% BAYESIAN SUBTYPE ANALYSIS (NEW IN V3)
fprintf('-----------------------------------------------------------------\n');
fprintf('BAYESIAN SUBTYPE ANALYSIS (BIC-based Bayes Factor)\n');
fprintf('-----------------------------------------------------------------\n');
fprintf('Purpose: Test evidence for null (equivalence) vs. difference\n');
fprintf('Especially important for Global (N=13) comparisons\n\n');

% Storage for Bayesian results
Level2_Bayesian = table('Size', [nICN * 3, 10], ...  % 3 pairwise comparisons involving Global
    'VariableTypes', {'double', 'cell', 'cell', 'cell', 'double', 'double', 'double', 'double', 'double', 'cell'}, ...
    'VariableNames', {'ICN', 'ICN_Label', 'Metric', 'Comparison', 'N_group1', 'N_group2', 'BIC_null', 'BIC_alt', 'BF01', 'Interpretation'});

bayesIdx = 0;

for n = 1:nICN
    colName = sprintf('%s_ICN%02d_%s', primaryContrast, n, primaryMetric);
    
    y_all = subtypeData.(colName);
    groups_all = subtypeData.AphasiaType;
    
    % Comparisons involving Global (the smallest group)
    comparisons = {'Global', 'Broca'; 'Global', 'Anomic'; 'Global', 'Conduction'};
    
    for c = 1:size(comparisons, 1)
        group1 = comparisons{c, 1};
        group2 = comparisons{c, 2};
        
        idx1 = strcmp(groups_all, group1);
        idx2 = strcmp(groups_all, group2);
        
        y1 = y_all(idx1);
        y2 = y_all(idx2);
        
        y1 = y1(~isnan(y1));
        y2 = y2(~isnan(y2));
        
        n1 = length(y1);
        n2 = length(y2);
        
        if n1 < 5 || n2 < 5
            continue;
        end
        
        bayesIdx = bayesIdx + 1;
        
        % Combined data
        y_combined = [y1; y2];
        group_indicator = [zeros(n1, 1); ones(n2, 1)];
        n_total = n1 + n2;
        
        % Null model (intercept only): y ~ 1
        RSS_null = sum((y_combined - mean(y_combined)).^2);
        k_null = 1;  % 1 parameter (mean)
        BIC_null = n_total * log(RSS_null/n_total) + k_null * log(n_total);
        
        % Alternative model (group effect): y ~ group
        mean1 = mean(y1);
        mean2 = mean(y2);
        predicted = [repmat(mean1, n1, 1); repmat(mean2, n2, 1)];
        RSS_alt = sum((y_combined - predicted).^2);
        k_alt = 2;  % 2 parameters (two means)
        BIC_alt = n_total * log(RSS_alt/n_total) + k_alt * log(n_total);
        
        % Bayes Factor
        [~, BF01, interpretation] = bic_bayes_factor(BIC_null, BIC_alt);
        
        % Store
        Level2_Bayesian.ICN(bayesIdx) = n;
        Level2_Bayesian.ICN_Label{bayesIdx} = trueICN_labels{n};
        Level2_Bayesian.Metric{bayesIdx} = primaryMetric;
        Level2_Bayesian.Comparison{bayesIdx} = sprintf('%s vs %s', group1, group2);
        Level2_Bayesian.N_group1(bayesIdx) = n1;
        Level2_Bayesian.N_group2(bayesIdx) = n2;
        Level2_Bayesian.BIC_null(bayesIdx) = BIC_null;
        Level2_Bayesian.BIC_alt(bayesIdx) = BIC_alt;
        Level2_Bayesian.BF01(bayesIdx) = BF01;
        Level2_Bayesian.Interpretation{bayesIdx} = interpretation;
    end
end

% Trim empty rows
Level2_Bayesian = Level2_Bayesian(1:bayesIdx, :);

% Display Bayesian results
fprintf('Bayesian results for Global comparisons (IRi, %d comparisons):\n\n', bayesIdx);

% Focus on ICN17 (key finding from v2)
icn17_bayes = Level2_Bayesian(Level2_Bayesian.ICN == 17, :);
fprintf('*** ICN17 (Motor Network - Key Finding) ***\n');
if height(icn17_bayes) > 0
    fprintf('%-25s %8s %8s %10s %s\n', 'Comparison', 'N1', 'N2', 'BF01', 'Interpretation');
    fprintf('%s\n', repmat('-', 1, 80));
    for i = 1:height(icn17_bayes)
        fprintf('%-25s %8d %8d %10.2f %s\n', ...
            icn17_bayes.Comparison{i}, icn17_bayes.N_group1(i), icn17_bayes.N_group2(i), ...
            icn17_bayes.BF01(i), icn17_bayes.Interpretation{i});
    end
end
fprintf('\n');

% Summary across all ICNs
nEvidenceNull = sum(Level2_Bayesian.BF01 > BF_threshold_moderate);
nEvidenceAlt = sum(Level2_Bayesian.BF01 < (1/BF_threshold_moderate));
nInconclusive = bayesIdx - nEvidenceNull - nEvidenceAlt;

fprintf('Bayesian Summary (BF threshold = %.1f):\n', BF_threshold_moderate);
fprintf('  Evidence FOR null (equivalence): %d/%d comparisons\n', nEvidenceNull, bayesIdx);
fprintf('  Evidence FOR alternative (difference): %d/%d comparisons\n', nEvidenceAlt, bayesIdx);
fprintf('  Inconclusive: %d/%d comparisons\n\n', nInconclusive, bayesIdx);

if nEvidenceNull > bayesIdx/2
    fprintf('RECOMMENDATION: Many comparisons show evidence for null.\n');
    fprintf('Consider collapsing Global with other severe subtypes (Broca) for power.\n');
end
fprintf('\n');

%% =========================================================================
%% LEVEL 3: SEVERITY EFFECT (PARTIAL CORRELATIONS)
%% =========================================================================
fprintf('=================================================================\n');
fprintf('LEVEL 3: Severity Effect (WAB-AQ Partial Correlations)\n');
fprintf('         Controlling for Age, Sex, Days_Post_Stroke\n');
fprintf('=================================================================\n');
fprintf('HYPOTHESIS H2a: WAB-AQ correlates positively with IRi in language ICNs\n');
fprintf('HYPOTHESIS H2b: Stronger correlations for language vs non-language\n');
fprintf('HYPOTHESIS H2c: MANi correlations weaker/more diffuse than IRi\n\n');

aphasiaData_L3 = masterWide(isAphasia, :);
fprintf('Aphasia patients: N = %d\n\n', height(aphasiaData_L3));

% Storage
nTests_L3 = length(metrics) * nICN;
Level3_Results = table('Size', [nTests_L3, 14], ...
    'VariableTypes', {'double', 'cell', 'cell', 'logical', ...
                      'double', 'double', 'double', 'double', 'double', ...
                      'double', 'double', 'double', 'logical', 'double'}, ...
    'VariableNames', {'ICN', 'ICN_Label', 'Metric', 'IsLanguageICN', ...
                      'N', 'Partial_rho', 'p_uncorrected', 'Abs_rho', ...
                      'Rho_CI_lo', 'Rho_CI_hi', ...
                      'Zero_order_rho', 'Direction', 'Significant_FDR', 'p_FDR'});

resultIdx = 0;
for m = 1:length(metrics)
    metric = metrics{m};
    
    for n = 1:nICN
        resultIdx = resultIdx + 1;
        
        colName = sprintf('%s_ICN%02d_%s', primaryContrast, n, metric);
        
        wab = aphasiaData_L3.WAB_AQ;
        icn = aphasiaData_L3.(colName);
        age = aphasiaData_L3.Age_At_Stroke;
        sex = aphasiaData_L3.SexNumeric;
        days = aphasiaData_L3.Days_Post_Stroke;
        
        covariates = [age, sex, days];
        
        [rho_partial, p_partial, n_valid] = partial_corr_spearman(wab, icn, covariates);
        
        % 95% CI for correlation (Fisher's z transformation)
        if ~isnan(rho_partial) && n_valid > 3
            z = 0.5 * log((1 + rho_partial) / (1 - rho_partial));
            se_z = 1 / sqrt(n_valid - 3);
            z_lo = z - 1.96 * se_z;
            z_hi = z + 1.96 * se_z;
            rho_ci_lo = (exp(2*z_lo) - 1) / (exp(2*z_lo) + 1);
            rho_ci_hi = (exp(2*z_hi) - 1) / (exp(2*z_hi) + 1);
        else
            rho_ci_lo = NaN;
            rho_ci_hi = NaN;
        end
        
        % Zero-order for comparison
        validIdx = ~isnan(wab) & ~isnan(icn);
        if sum(validIdx) >= 10
            rho_zero = corr(wab(validIdx), icn(validIdx), 'Type', 'Spearman');
        else
            rho_zero = NaN;
        end
        
        % Direction
        if ~isnan(rho_partial)
            direction = sign(rho_partial);
        else
            direction = 0;
        end
        
        % Store
        Level3_Results.ICN(resultIdx) = n;
        Level3_Results.ICN_Label{resultIdx} = trueICN_labels{n};
        Level3_Results.Metric{resultIdx} = metric;
        Level3_Results.IsLanguageICN(resultIdx) = ismember(n, languageICNs);
        Level3_Results.N(resultIdx) = n_valid;
        Level3_Results.Partial_rho(resultIdx) = rho_partial;
        Level3_Results.p_uncorrected(resultIdx) = p_partial;
        Level3_Results.Abs_rho(resultIdx) = abs(rho_partial);
        Level3_Results.Rho_CI_lo(resultIdx) = rho_ci_lo;
        Level3_Results.Rho_CI_hi(resultIdx) = rho_ci_hi;
        Level3_Results.Zero_order_rho(resultIdx) = rho_zero;
        Level3_Results.Direction(resultIdx) = direction;
    end
end

% FDR correction
[Level3_Results.Significant_FDR, ~, Level3_Results.p_FDR] = fdr_bh(Level3_Results.p_uncorrected, fdrQ);

Level3_Results = sortrows(Level3_Results, 'Abs_rho', 'descend');

sigL3 = Level3_Results(Level3_Results.Significant_FDR, :);
fprintf('Significant partial correlations (FDR q < %.2f): %d / %d tests\n\n', fdrQ, height(sigL3), nTests_L3);

% H2a: Language network correlations
fprintf('=== H2a: Language Network IRi-WAB Correlations ===\n');
L3_lang_IRi = Level3_Results(Level3_Results.IsLanguageICN & strcmp(Level3_Results.Metric, 'IRi'), :);
L3_lang_IRi = sortrows(L3_lang_IRi, 'Abs_rho', 'descend');

fprintf('%-6s %-40s %10s %15s %10s\n', 'ICN', 'Label', 'Part_rho', '95%CI', 'p_unc');
fprintf('%s\n', repmat('-', 1, 90));
for i = 1:height(L3_lang_IRi)
    fdrStr = '';
    if L3_lang_IRi.Significant_FDR(i)
        fdrStr = '***';
    end
    fprintf('ICN%02d  %-40s %10.3f [%5.3f,%5.3f] %10.4f %s\n', ...
            L3_lang_IRi.ICN(i), L3_lang_IRi.ICN_Label{i}(1:min(40,end)), ...
            L3_lang_IRi.Partial_rho(i), L3_lang_IRi.Rho_CI_lo(i), L3_lang_IRi.Rho_CI_hi(i), ...
            L3_lang_IRi.p_uncorrected(i), fdrStr);
end
fprintf('\n');

% H2b: Language vs Non-language effect sizes
fprintf('=== H2b: Network Specificity (IRi correlations) ===\n');
L3_IRi = Level3_Results(strcmp(Level3_Results.Metric, 'IRi'), :);
langEffects = abs(L3_IRi.Partial_rho(L3_IRi.IsLanguageICN));
nonLangEffects = abs(L3_IRi.Partial_rho(~L3_IRi.IsLanguageICN));

fprintf('Language networks (N=%d): Mean |rho| = %.3f (SD = %.3f)\n', ...
        length(langEffects), mean(langEffects, 'omitnan'), std(langEffects, 'omitnan'));
fprintf('Non-language networks (N=%d): Mean |rho| = %.3f (SD = %.3f)\n', ...
        length(nonLangEffects), mean(nonLangEffects, 'omitnan'), std(nonLangEffects, 'omitnan'));

if length(langEffects) >= 2 && length(nonLangEffects) >= 2
    [~, p_spec, ~, stats_spec] = ttest2(langEffects, nonLangEffects);
    d_spec = cohens_d(langEffects, nonLangEffects);
    fprintf('Comparison: t(%.0f) = %.2f, p = %.4f, Cohen''s d = %.2f\n', ...
            stats_spec.df, stats_spec.tstat, p_spec, d_spec);
    Level5_specificity_t = stats_spec.tstat;
    Level5_specificity_df = stats_spec.df;
    Level5_specificity_p = p_spec;
    Level5_specificity_d = d_spec;
else
    p_spec = NaN; d_spec = NaN;
    Level5_specificity_t = NaN;
    Level5_specificity_df = NaN;
    Level5_specificity_p = NaN;
    Level5_specificity_d = NaN;
end
fprintf('\n');

% H2c: IRi vs MANi correlations
fprintf('=== H2c: IRi vs MANi Correlation Comparison ===\n');
L3_IRi_all = Level3_Results(strcmp(Level3_Results.Metric, 'IRi'), :);
L3_MANi_all = Level3_Results(strcmp(Level3_Results.Metric, 'MANi'), :);

meanAbsRho_IRi = mean(abs(L3_IRi_all.Partial_rho), 'omitnan');
meanAbsRho_MANi = mean(abs(L3_MANi_all.Partial_rho), 'omitnan');
nSig_IRi = sum(L3_IRi_all.Significant_FDR);
nSig_MANi = sum(L3_MANi_all.Significant_FDR);

fprintf('IRi:  Mean |partial rho| = %.3f, Significant = %d/%d\n', meanAbsRho_IRi, nSig_IRi, nICN);
fprintf('MANi: Mean |partial rho| = %.3f, Significant = %d/%d\n', meanAbsRho_MANi, nSig_MANi, nICN);
fprintf('\n');

% Top overall correlations with CIs
fprintf('=== Top 10 Overall Correlations (All Metrics) ===\n');
sigL3_sorted = sortrows(sigL3, 'Abs_rho', 'descend');
fprintf('%-6s %-8s %-30s %10s %15s %10s\n', 'ICN', 'Metric', 'Label', 'Part_rho', '95%CI', 'p_FDR');
fprintf('%s\n', repmat('-', 1, 95));
for i = 1:min(10, height(sigL3_sorted))
    fprintf('ICN%02d  %-8s %-30s %10.3f [%5.3f,%5.3f] %10.4f\n', ...
            sigL3_sorted.ICN(i), sigL3_sorted.Metric{i}, ...
            sigL3_sorted.ICN_Label{i}(1:min(30,end)), ...
            sigL3_sorted.Partial_rho(i), sigL3_sorted.Rho_CI_lo(i), sigL3_sorted.Rho_CI_hi(i), ...
            sigL3_sorted.p_FDR(i));
end
fprintf('\n');

%% =========================================================================
%% LEVEL 4: CONTRAST COMPARISON (C2 vs C1)
%% =========================================================================
fprintf('=================================================================\n');
fprintf('LEVEL 4: Contrast Comparison\n');
fprintf('         Does Naming>Abstract (C2) predict WAB beyond Task>Rest (C1)?\n');
fprintf('=================================================================\n');
fprintf('HYPOTHESIS H4a: C2 adds prediction beyond C1 for language networks\n');
fprintf('HYPOTHESIS H4b: Domain-general ICNs similar C1/C2 relationships\n\n');

hasBoth = masterWide.HasCon0001 & masterWide.HasCon0002 & isAphasia;
contrastData = masterWide(hasBoth, :);

fprintf('Aphasia patients with both contrasts: N = %d\n\n', height(contrastData));

% Storage
nTests_L4 = length(metrics) * nICN;
Level4_Results = table('Size', [nTests_L4, 16], ...
    'VariableTypes', {'double', 'cell', 'cell', 'logical', 'double', ...
                      'double', 'double', 'double', 'double', 'double', 'double', ...
                      'double', 'double', 'double', 'logical', 'double'}, ...
    'VariableNames', {'ICN', 'ICN_Label', 'Metric', 'IsLanguageICN', 'N', ...
                      'Beta_C2', 't_C2', 'p_C2', 'Beta_C1', 't_C1', 'p_C1', ...
                      'R2_full', 'R2_C1only', 'R2_improvement', 'Significant_FDR', 'p_FDR'});

resultIdx = 0;
for m = 1:length(metrics)
    metric = metrics{m};
    
    for n = 1:nICN
        resultIdx = resultIdx + 1;
        
        colC1 = sprintf('C1_ICN%02d_%s', n, metric);
        colC2 = sprintf('C2_ICN%02d_%s', n, metric);
        
        wab = contrastData.WAB_AQ;
        c1 = contrastData.(colC1);
        c2 = contrastData.(colC2);
        
        validIdx = ~isnan(wab) & ~isnan(c1) & ~isnan(c2);
        wab_valid = wab(validIdx);
        c1_valid = c1(validIdx);
        c2_valid = c2(validIdx);
        
        n_obs = sum(validIdx);
        
        if n_obs >= 10
            tbl = table(wab_valid, c2_valid, c1_valid, 'VariableNames', {'WAB', 'C2', 'C1'});
            mdl = fitlm(tbl, 'WAB ~ C2 + C1');
            
            beta_C2 = mdl.Coefficients.Estimate(2);
            t_C2 = mdl.Coefficients.tStat(2);
            p_C2 = mdl.Coefficients.pValue(2);
            
            beta_C1 = mdl.Coefficients.Estimate(3);
            t_C1 = mdl.Coefficients.tStat(3);
            p_C1 = mdl.Coefficients.pValue(3);
            
            R2_full = mdl.Rsquared.Ordinary;
            
            mdl_C1 = fitlm(tbl, 'WAB ~ C1');
            R2_C1only = mdl_C1.Rsquared.Ordinary;
            
            R2_impr = R2_full - R2_C1only;
        else
            beta_C2 = NaN; t_C2 = NaN; p_C2 = NaN;
            beta_C1 = NaN; t_C1 = NaN; p_C1 = NaN;
            R2_full = NaN; R2_C1only = NaN; R2_impr = NaN;
        end
        
        % Store
        Level4_Results.ICN(resultIdx) = n;
        Level4_Results.ICN_Label{resultIdx} = trueICN_labels{n};
        Level4_Results.Metric{resultIdx} = metric;
        Level4_Results.IsLanguageICN(resultIdx) = ismember(n, languageICNs);
        Level4_Results.N(resultIdx) = n_obs;
        Level4_Results.Beta_C2(resultIdx) = beta_C2;
        Level4_Results.t_C2(resultIdx) = t_C2;
        Level4_Results.p_C2(resultIdx) = p_C2;
        Level4_Results.Beta_C1(resultIdx) = beta_C1;
        Level4_Results.t_C1(resultIdx) = t_C1;
        Level4_Results.p_C1(resultIdx) = p_C1;
        Level4_Results.R2_full(resultIdx) = R2_full;
        Level4_Results.R2_C1only(resultIdx) = R2_C1only;
        Level4_Results.R2_improvement(resultIdx) = R2_impr;
    end
end

% FDR correction on C2 p-values
[Level4_Results.Significant_FDR, ~, Level4_Results.p_FDR] = fdr_bh(Level4_Results.p_C2, fdrQ);

Level4_Results = sortrows(Level4_Results, 'p_C2');

sigL4 = Level4_Results(Level4_Results.Significant_FDR, :);
fprintf('Networks where C2 adds beyond C1 (FDR): %d / %d\n\n', height(sigL4), nTests_L4);

% H4a: Language networks
fprintf('=== H4a: Language Network C2 Incremental Effects (IRi) ===\n');
L4_lang = Level4_Results(Level4_Results.IsLanguageICN & strcmp(Level4_Results.Metric, 'IRi'), :);
fprintf('%-6s %-35s %10s %10s %10s\n', 'ICN', 'Label', 'Beta_C2', 'R2_impr', 'p_C2');
fprintf('%s\n', repmat('-', 1, 80));
for i = 1:height(L4_lang)
    fdrStr = '';
    if L4_lang.Significant_FDR(i)
        fdrStr = '***';
    end
    fprintf('ICN%02d  %-35s %10.2f %10.3f %10.4f %s\n', ...
            L4_lang.ICN(i), L4_lang.ICN_Label{i}(1:min(35,end)), ...
            L4_lang.Beta_C2(i), L4_lang.R2_improvement(i), L4_lang.p_C2(i), fdrStr);
end
fprintf('\n');

%% =========================================================================
%% LEVEL 5: NETWORK SPECIFICITY SUMMARY
%% =========================================================================
fprintf('=================================================================\n');
fprintf('LEVEL 5: Network Specificity Summary\n');
fprintf('=================================================================\n');

Level5_Results = struct();

% Store from Level 3
Level5_Results.LanguageICNs = languageICNs;
Level5_Results.LanguageLabels = languageLabels;
Level5_Results.Lang_MeanEffect = mean(langEffects, 'omitnan');
Level5_Results.Lang_SDEffect = std(langEffects, 'omitnan');
Level5_Results.NonLang_MeanEffect = mean(nonLangEffects, 'omitnan');
Level5_Results.NonLang_SDEffect = std(nonLangEffects, 'omitnan');
Level5_Results.Specificity_t = Level5_specificity_t;
Level5_Results.Specificity_df = Level5_specificity_df;
Level5_Results.Specificity_p = Level5_specificity_p;
Level5_Results.Specificity_d = Level5_specificity_d;

% Compute R² improvement by network type from Level 4
L4_IRi = Level4_Results(strcmp(Level4_Results.Metric, 'IRi'), :);
lang_R2impr = L4_IRi.R2_improvement(L4_IRi.IsLanguageICN);
nonlang_R2impr = L4_IRi.R2_improvement(~L4_IRi.IsLanguageICN);

Level5_Results.Lang_MeanR2impr = mean(lang_R2impr, 'omitnan');
Level5_Results.NonLang_MeanR2impr = mean(nonlang_R2impr, 'omitnan');

% IRi profile across all ICNs for aphasia patients
fprintf('IRi Profile (C2, Aphasia patients) - Mean ± SD:\n');
fprintf('%-6s %-45s %12s %8s\n', 'ICN', 'Label', 'IRi Mean±SD', 'Lang?');
fprintf('%s\n', repmat('-', 1, 80));

aphasiaData_profile = masterWide(isAphasia, :);
IRi_means = NaN(nICN, 1);
IRi_sds = NaN(nICN, 1);
for n = 1:nICN
    colName = sprintf('C2_ICN%02d_IRi', n);
    vals = aphasiaData_profile.(colName);
    IRi_means(n) = mean(vals, 'omitnan');
    IRi_sds(n) = std(vals, 'omitnan');
    
    langStr = '';
    if ismember(n, languageICNs)
        langStr = '*** LANG';
    elseif ismember(n, compensatoryICNs)
        langStr = 'COMP';
    end
    
    fprintf('ICN%02d  %-45s %6.4f±%.4f %8s\n', n, trueICN_labels{n}(1:min(45,end)), ...
            IRi_means(n), IRi_sds(n), langStr);
end

% Store in Level5
Level5_Results.IRi_Means = IRi_means;
Level5_Results.IRi_SDs = IRi_sds;
fprintf('\n');

%% =========================================================================
%% IRi THRESHOLD SENSITIVITY ANALYSIS (NEW IN V3)
%% =========================================================================
fprintf('=================================================================\n');
fprintf('IRi THRESHOLD SENSITIVITY ANALYSIS (NEW IN V3)\n');
fprintf('=================================================================\n');
fprintf('Purpose: Test robustness of key findings to IRi QC threshold choice\n');
fprintf('Main threshold: %.1f (from Script 3b_v3)\n', IRi_threshold);
fprintf('Testing: 0.7, 0.8, 0.9\n\n');

% Load raw matrices with QC info from 3b
matricesFile = fullfile(inputDir, 'ARC_03b_v3_ICN_Matrices.mat');
if exist(matricesFile, 'file')
    rawData = load(matricesFile, 'IRi_matrix', 'MANi_matrix', 'Vari_matrix', ...
                   'PatientID', 'ContrastNum', 'AphasiaType', 'GroupRole', 'WAB_AQ', ...
                   'IRi_sums', 'lowIRi_flag', 'trueICN_idx');
    
    fprintf('Loaded raw matrices for sensitivity analysis.\n\n');
    
    % Define thresholds to test
    sensitivityThresholds = [0.7, 0.8, 0.9];
    nThresholds = length(sensitivityThresholds);
    
    % Storage for sensitivity results
    % Focus on Level 3 key finding: ICN17 correlation with WAB-AQ (IRi)
    Sensitivity_Results = table('Size', [nThresholds, 8], ...
        'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'logical'}, ...
        'VariableNames', {'IRi_Threshold', 'N_Excluded', 'N_Retained', 'Pct_Retained', ...
                          'ICN17_rho', 'ICN17_p', 'ICN17_N', 'ICN17_Significant'});
    
    % Also track all language ICN correlations
    Sensitivity_LangICN = table('Size', [nThresholds * length(languageICNs), 6], ...
        'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'logical'}, ...
        'VariableNames', {'IRi_Threshold', 'ICN', 'rho', 'p_value', 'N', 'Significant'});
    langIdx = 0;
    
    for t = 1:nThresholds
        thresh = sensitivityThresholds(t);
        
        % Re-apply QC at this threshold
        newExcluded = rawData.IRi_sums < thresh;
        
        % Filter to aphasia patients with Con_0002
        isCon002 = strcmp(rawData.ContrastNum, '0002');
        isAphasiaRaw = strcmp(rawData.GroupRole, 'Aphasia');
        validForAnalysis = isCon002 & isAphasiaRaw & ~newExcluded;
        
        nExcluded = sum(newExcluded & isCon002 & isAphasiaRaw);
        nRetained = sum(validForAnalysis);
        
        Sensitivity_Results.IRi_Threshold(t) = thresh;
        Sensitivity_Results.N_Excluded(t) = nExcluded;
        Sensitivity_Results.N_Retained(t) = nRetained;
        Sensitivity_Results.Pct_Retained(t) = 100 * nRetained / sum(isCon002 & isAphasiaRaw);
        
        % Compute ICN17 correlation at this threshold
        wab_sens = rawData.WAB_AQ(validForAnalysis);
        icn17_sens = rawData.IRi_matrix(validForAnalysis, 17);
        
        validCorr = ~isnan(wab_sens) & ~isnan(icn17_sens);
        if sum(validCorr) >= 10
            [rho_17, p_17] = corr(wab_sens(validCorr), icn17_sens(validCorr), 'Type', 'Spearman');
            Sensitivity_Results.ICN17_rho(t) = rho_17;
            Sensitivity_Results.ICN17_p(t) = p_17;
            Sensitivity_Results.ICN17_N(t) = sum(validCorr);
            Sensitivity_Results.ICN17_Significant(t) = p_17 < 0.05;
        else
            Sensitivity_Results.ICN17_rho(t) = NaN;
            Sensitivity_Results.ICN17_p(t) = NaN;
            Sensitivity_Results.ICN17_N(t) = sum(validCorr);
            Sensitivity_Results.ICN17_Significant(t) = false;
        end
        
        % Compute language ICN correlations
        for li = 1:length(languageICNs)
            langIdx = langIdx + 1;
            icnNum = languageICNs(li);
            
            icn_sens = rawData.IRi_matrix(validForAnalysis, icnNum);
            validCorr_lang = ~isnan(wab_sens) & ~isnan(icn_sens);
            
            Sensitivity_LangICN.IRi_Threshold(langIdx) = thresh;
            Sensitivity_LangICN.ICN(langIdx) = icnNum;
            
            if sum(validCorr_lang) >= 10
                [rho_lang, p_lang] = corr(wab_sens(validCorr_lang), icn_sens(validCorr_lang), 'Type', 'Spearman');
                Sensitivity_LangICN.rho(langIdx) = rho_lang;
                Sensitivity_LangICN.p_value(langIdx) = p_lang;
                Sensitivity_LangICN.N(langIdx) = sum(validCorr_lang);
                Sensitivity_LangICN.Significant(langIdx) = p_lang < 0.05;
            else
                Sensitivity_LangICN.rho(langIdx) = NaN;
                Sensitivity_LangICN.p_value(langIdx) = NaN;
                Sensitivity_LangICN.N(langIdx) = sum(validCorr_lang);
                Sensitivity_LangICN.Significant(langIdx) = false;
            end
        end
    end
    
    % Display sensitivity results
    fprintf('ICN17 (Motor Network) Correlation Sensitivity:\n');
    fprintf('%-12s %10s %10s %10s %10s %12s\n', 'Threshold', 'N_Retained', 'rho', 'p-value', 'N_valid', 'Significant');
    fprintf('%s\n', repmat('-', 1, 70));
    for t = 1:nThresholds
        sigStr = '';
        if Sensitivity_Results.ICN17_Significant(t)
            sigStr = '***';
        end
        fprintf('%.1f          %10d %10.3f %10.4f %10d %12s\n', ...
            Sensitivity_Results.IRi_Threshold(t), Sensitivity_Results.N_Retained(t), ...
            Sensitivity_Results.ICN17_rho(t), Sensitivity_Results.ICN17_p(t), ...
            Sensitivity_Results.ICN17_N(t), sigStr);
    end
    fprintf('\n');
    
    % Check robustness
    allSig = all(Sensitivity_Results.ICN17_Significant);
    rhoRange = max(Sensitivity_Results.ICN17_rho) - min(Sensitivity_Results.ICN17_rho);
    
    if allSig && rhoRange < 0.1
        fprintf('✓ ROBUST: ICN17 correlation significant at all thresholds (rho range = %.3f)\n\n', rhoRange);
        sensitivityConclusion = 'ROBUST';
    elseif allSig
        fprintf('✓ SIGNIFICANT at all thresholds, but rho varies (range = %.3f)\n\n', rhoRange);
        sensitivityConclusion = 'SIGNIFICANT_VARIABLE';
    else
        fprintf('! NOT ROBUST: ICN17 significance varies across thresholds\n\n');
        sensitivityConclusion = 'NOT_ROBUST';
    end
    
    % Language ICN summary
    fprintf('Language ICN Correlations Across Thresholds:\n');
    fprintf('%-12s %-8s %10s %10s %12s\n', 'Threshold', 'ICN', 'rho', 'p-value', 'Significant');
    fprintf('%s\n', repmat('-', 1, 55));
    for i = 1:height(Sensitivity_LangICN)
        sigStr = '';
        if Sensitivity_LangICN.Significant(i)
            sigStr = '***';
        end
        fprintf('%.1f          ICN%02d    %10.3f %10.4f %12s\n', ...
            Sensitivity_LangICN.IRi_Threshold(i), Sensitivity_LangICN.ICN(i), ...
            Sensitivity_LangICN.rho(i), Sensitivity_LangICN.p_value(i), sigStr);
    end
    fprintf('\n');
    
else
    fprintf('WARNING: Raw matrices file not found. Skipping sensitivity analysis.\n');
    fprintf('File expected: %s\n\n', matricesFile);
    Sensitivity_Results = table();
    Sensitivity_LangICN = table();
    sensitivityConclusion = 'NOT_RUN';
end

%% =========================================================================
%% SAVE ALL RESULTS
%% =========================================================================
fprintf('=================================================================\n');
fprintf('Saving all results (V3)\n');
fprintf('=================================================================\n');

% Save Level 1
L1_file = fullfile(outputDir, 'ARC_03c_v3_Level1_DiseaseEffect.csv');
writetable(Level1_Results, L1_file);
fprintf('Saved: %s\n', L1_file);

% Save Level 2 (Frequentist)
L2_file = fullfile(outputDir, 'ARC_03c_v3_Level2_SubtypeEffect.csv');
writetable(Level2_Results, L2_file);
fprintf('Saved: %s\n', L2_file);

% Save Level 2 Bayesian (NEW in V3)
L2_bayes_file = fullfile(outputDir, 'ARC_03c_v3_Level2_Bayesian.csv');
writetable(Level2_Bayesian, L2_bayes_file);
fprintf('Saved: %s\n', L2_bayes_file);

% Save Level 3
L3_file = fullfile(outputDir, 'ARC_03c_v3_Level3_SeverityCorrelations.csv');
writetable(Level3_Results, L3_file);
fprintf('Saved: %s\n', L3_file);

% Save Level 4
L4_file = fullfile(outputDir, 'ARC_03c_v3_Level4_ContrastComparison.csv');
writetable(Level4_Results, L4_file);
fprintf('Saved: %s\n', L4_file);

% Save IRi Sensitivity Analysis (NEW in V3)
if ~isempty(Sensitivity_Results)
    sens_file = fullfile(outputDir, 'ARC_03c_v3_IRi_Sensitivity.csv');
    writetable(Sensitivity_Results, sens_file);
    fprintf('Saved: %s\n', sens_file);
    
    sens_lang_file = fullfile(outputDir, 'ARC_03c_v3_IRi_Sensitivity_LangICN.csv');
    writetable(Sensitivity_LangICN, sens_lang_file);
    fprintf('Saved: %s\n', sens_lang_file);
end

% Save all results to .mat (V3)
allResultsFile = fullfile(outputDir, 'ARC_03c_v3_AllResults.mat');
save(allResultsFile, 'Level1_Results', 'Level2_Results', 'Level2_Bayesian', ...
     'Level3_Results', 'Level4_Results', 'Level5_Results', ...
     'Sensitivity_Results', 'Sensitivity_LangICN', 'sensitivityConclusion', ...
     'fdrQ', 'primaryMetric', 'primaryContrast', 'metrics', 'metricLabels', ...
     'trueICN_labels', 'nICN', 'languageICNs', 'languageLabels', ...
     'compensatoryICNs', 'nonLanguageICNs', 'IRi_threshold', ...
     'BF_threshold_moderate', 'nBootstrap');
fprintf('Saved: %s\n', allResultsFile);

%% =========================================================================
%% GENERATE SUMMARY REPORT
%% =========================================================================
reportFile = fullfile(outputDir, 'ARC_03c_v3_Summary.txt');
fid = fopen(reportFile, 'w');

fprintf(fid, '=================================================================\n');
fprintf(fid, 'ARC_03c_v3 STATISTICAL ANALYSIS SUMMARY\n');
fprintf(fid, '=================================================================\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));

fprintf(fid, 'V3 CHANGES FROM V2\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, '- Input: ARC_03b_v3_Output (IRi QC-cleaned, threshold=%.1f)\n', IRi_threshold);
fprintf(fid, '- NEW: Bayesian subtype analysis (BIC-based Bayes Factor)\n');
fprintf(fid, '- NEW: IRi threshold sensitivity analysis (0.7, 0.8, 0.9)\n');
fprintf(fid, '- REMOVED: Level 6 (LASSO) - circular with WAB-AQ definition\n');
fprintf(fid, '- Mediation analysis remains in Script 4 (with lesion data)\n\n');

fprintf(fid, 'CORRECTED METRICS\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'PRIMARY:     IRi (ICNiRelativeInvolvement)\n');
fprintf(fid, '             Proportion of whole-brain activation in each ICN\n');
fprintf(fid, '             Sensitive to extensive/distributed activations\n\n');
fprintf(fid, 'SECONDARY:   MANi (NormalisedMeanICNiActivation)\n');
fprintf(fid, '             Normalized mean Z-score within ICN\n');
fprintf(fid, '             Sensitive to focal intense activations\n\n');
fprintf(fid, 'EXPLORATORY: logVari (log-transformed variance of Z-scores)\n');
fprintf(fid, '             Activation heterogeneity within ICN\n\n');

fprintf(fid, 'ANALYSIS PARAMETERS\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Primary metric: %s\n', primaryMetric);
fprintf(fid, 'Primary contrast: %s (Naming > Abstract)\n', primaryContrast);
fprintf(fid, 'FDR threshold: q = %.2f\n', fdrQ);
fprintf(fid, 'Covariates: Age_At_Stroke, Sex, Days_Post_Stroke\n');
fprintf(fid, 'IRi QC threshold: %.1f\n', IRi_threshold);
fprintf(fid, 'Bayesian BF threshold (moderate): %.1f\n\n', BF_threshold_moderate);

fprintf(fid, 'A PRIORI LANGUAGE NETWORKS\n');
fprintf(fid, '-----------------------------------------------------------------\n');
for i = 1:length(languageICNs)
    fprintf(fid, '  ICN%02d: %s\n', languageICNs(i), trueICN_labels{languageICNs(i)});
end
fprintf(fid, '\nCompensatory network (exploratory):\n');
fprintf(fid, '  ICN15: %s\n\n', trueICN_labels{15});

fprintf(fid, 'SAMPLE (POST-IRi QC)\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Aphasia patients: N = %d\n', nAphasia);
fprintf(fid, 'Stroke Controls (No Aphasia): N = %d\n', nControl);
fprintf(fid, 'Total: N = %d\n\n', nAphasia + nControl);

fprintf(fid, '=================================================================\n');
fprintf(fid, 'HYPOTHESIS TESTING RESULTS\n');
fprintf(fid, '=================================================================\n\n');

% Level 1
fprintf(fid, 'LEVEL 1: DISEASE EFFECT (Aphasia vs Stroke Control)\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Total tests: %d\n', nTests_L1);
fprintf(fid, 'Significant (FDR q<%.2f): %d\n\n', fdrQ, sum(Level1_Results.Significant_FDR));

fprintf(fid, 'H1a: Aphasia shows reduced IRi in language networks?\n');
L1_lang_IRi = Level1_Results(Level1_Results.IsLanguageICN & strcmp(Level1_Results.Metric, 'IRi'), :);
for i = 1:height(L1_lang_IRi)
    direction = '';
    if L1_lang_IRi.Cohens_d(i) < 0
        direction = 'Aphasia < Control';
    else
        direction = 'Aphasia > Control';
    end
    sigStr = '';
    if L1_lang_IRi.Significant_FDR(i)
        sigStr = '*** FDR sig';
    end
    fprintf(fid, '  ICN%02d: d=%.2f [%.2f,%.2f] (%s), p=%.4f %s\n', ...
            L1_lang_IRi.ICN(i), L1_lang_IRi.Cohens_d(i), ...
            L1_lang_IRi.Cohens_d_CI_lo(i), L1_lang_IRi.Cohens_d_CI_hi(i), ...
            direction, L1_lang_IRi.p_uncorrected(i), sigStr);
end

fprintf(fid, '\nH1b: Aphasia shows increased IRi in compensatory (R FP)?\n');
L1_comp_IRi = Level1_Results(Level1_Results.IsCompensatoryICN & strcmp(Level1_Results.Metric, 'IRi'), :);
for i = 1:height(L1_comp_IRi)
    direction = '';
    if L1_comp_IRi.Cohens_d(i) > 0
        direction = 'Aphasia > Control (compensatory)';
    else
        direction = 'Aphasia < Control';
    end
    fprintf(fid, '  ICN%02d: d=%.2f [%.2f,%.2f] (%s), p=%.4f\n', ...
            L1_comp_IRi.ICN(i), L1_comp_IRi.Cohens_d(i), ...
            L1_comp_IRi.Cohens_d_CI_lo(i), L1_comp_IRi.Cohens_d_CI_hi(i), ...
            direction, L1_comp_IRi.p_uncorrected(i));
end

fprintf(fid, '\nOverall mean |Cohen''s d| for IRi: %.3f\n', ...
        mean(abs(Level1_Results.Cohens_d(strcmp(Level1_Results.Metric, 'IRi')))));
fprintf(fid, '\n');

% Level 2 (Frequentist)
fprintf(fid, 'LEVEL 2: SUBTYPE EFFECT (4-group ANCOVA)\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Groups: Anomic (%d), Broca (%d), Conduction (%d), Global (%d)\n', ...
        sum(isAnomic), sum(isBroca), sum(isConduction), sum(isGlobal));
fprintf(fid, 'Total tests: %d\n', nTests_L2);
fprintf(fid, 'Significant (FDR q<%.2f): %d\n\n', fdrQ, sum(Level2_Results.Significant_FDR));

if sum(Level2_Results.Significant_FDR) > 0
    sigL2_report = Level2_Results(Level2_Results.Significant_FDR, :);
    sigL2_report = sortrows(sigL2_report, 'Eta_squared', 'descend');
    fprintf(fid, 'Significant effects:\n');
    for i = 1:min(10, height(sigL2_report))
        fprintf(fid, '  ICN%02d %s: Eta²=%.3f, p_FDR=%.4f\n', ...
                sigL2_report.ICN(i), sigL2_report.Metric{i}, ...
                sigL2_report.Eta_squared(i), sigL2_report.p_FDR(i));
    end
else
    fprintf(fid, 'No significant subtype effects after FDR correction.\n');
end
fprintf(fid, '\n');

% Level 2 Bayesian (NEW in V3)
fprintf(fid, 'LEVEL 2 BAYESIAN: Evidence for Null vs Difference (Global comparisons)\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'BIC-based Bayes Factor approximation\n');
fprintf(fid, 'BF01 > %.1f: Moderate evidence FOR null (equivalence)\n', BF_threshold_moderate);
fprintf(fid, 'BF01 < %.2f: Moderate evidence FOR alternative (difference)\n\n', 1/BF_threshold_moderate);

if height(Level2_Bayesian) > 0
    % Focus on ICN17 (if present)
    icn17_bayes = Level2_Bayesian(Level2_Bayesian.ICN == 17, :);
    if height(icn17_bayes) > 0
        fprintf(fid, 'ICN17 (Motor Network) - Key finding from prior analysis:\n');
        for i = 1:height(icn17_bayes)
            fprintf(fid, '  %s: BF01=%.2f - %s\n', ...
                icn17_bayes.Comparison{i}, icn17_bayes.BF01(i), icn17_bayes.Interpretation{i});
        end
        fprintf(fid, '\n');
    end
    
    fprintf(fid, 'Summary across all ICNs:\n');
    fprintf(fid, '  Evidence FOR null (equivalence): %d/%d comparisons\n', nEvidenceNull, height(Level2_Bayesian));
    fprintf(fid, '  Evidence FOR alternative (difference): %d/%d comparisons\n', nEvidenceAlt, height(Level2_Bayesian));
    fprintf(fid, '  Inconclusive: %d/%d comparisons\n', nInconclusive, height(Level2_Bayesian));
    
    if nEvidenceNull > height(Level2_Bayesian)/2
        fprintf(fid, '\nRECOMMENDATION: Consider collapsing Global with Broca for power.\n');
    end
end
fprintf(fid, '\n');

% Level 3
fprintf(fid, 'LEVEL 3: SEVERITY EFFECT (Partial correlations with WAB-AQ)\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Total tests: %d\n', nTests_L3);
fprintf(fid, 'Significant (FDR q<%.2f): %d\n\n', fdrQ, sum(Level3_Results.Significant_FDR));

fprintf(fid, 'H2a: WAB-AQ correlates positively with IRi in language networks?\n');
L3_lang_report = Level3_Results(Level3_Results.IsLanguageICN & strcmp(Level3_Results.Metric, 'IRi'), :);
for i = 1:height(L3_lang_report)
    sigStr = '';
    if L3_lang_report.Significant_FDR(i)
        sigStr = '*** FDR sig';
    end
    fprintf(fid, '  ICN%02d: rho=%.3f [%.3f,%.3f], p=%.4f %s\n', ...
            L3_lang_report.ICN(i), L3_lang_report.Partial_rho(i), ...
            L3_lang_report.Rho_CI_lo(i), L3_lang_report.Rho_CI_hi(i), ...
            L3_lang_report.p_uncorrected(i), sigStr);
end

fprintf(fid, '\nH2b: Stronger correlations for language vs non-language?\n');
fprintf(fid, '  Language ICNs: Mean |rho| = %.3f (SD=%.3f)\n', ...
        Level5_Results.Lang_MeanEffect, Level5_Results.Lang_SDEffect);
fprintf(fid, '  Non-language ICNs: Mean |rho| = %.3f (SD=%.3f)\n', ...
        Level5_Results.NonLang_MeanEffect, Level5_Results.NonLang_SDEffect);
if ~isnan(Level5_Results.Specificity_t)
    fprintf(fid, '  t(%.0f)=%.2f, p=%.4f, d=%.2f\n', ...
            Level5_Results.Specificity_df, Level5_Results.Specificity_t, ...
            Level5_Results.Specificity_p, Level5_Results.Specificity_d);
end

fprintf(fid, '\nH2c: IRi vs MANi correlations?\n');
fprintf(fid, '  IRi: Mean |rho| = %.3f, Significant = %d/%d\n', meanAbsRho_IRi, nSig_IRi, nICN);
fprintf(fid, '  MANi: Mean |rho| = %.3f, Significant = %d/%d\n', meanAbsRho_MANi, nSig_MANi, nICN);

if sum(Level3_Results.Significant_FDR) > 0
    fprintf(fid, '\nTop significant correlations (all metrics):\n');
    sigL3_report = Level3_Results(Level3_Results.Significant_FDR, :);
    sigL3_report = sortrows(sigL3_report, 'Abs_rho', 'descend');
    for i = 1:min(10, height(sigL3_report))
        fprintf(fid, '  ICN%02d %s: rho=%.3f [%.3f,%.3f], p_FDR=%.4f\n', ...
                sigL3_report.ICN(i), sigL3_report.Metric{i}, ...
                sigL3_report.Partial_rho(i), sigL3_report.Rho_CI_lo(i), sigL3_report.Rho_CI_hi(i), ...
                sigL3_report.p_FDR(i));
    end
end
fprintf(fid, '\n');

% Level 4
fprintf(fid, 'LEVEL 4: CONTRAST COMPARISON (C2 adds beyond C1?)\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Model: WAB ~ C2_ICN + C1_ICN\n');
fprintf(fid, 'Total tests: %d\n', nTests_L4);
fprintf(fid, 'Significant C2 effects (FDR q<%.2f): %d\n\n', fdrQ, sum(Level4_Results.Significant_FDR));

fprintf(fid, 'H4a: C2 adds prediction for language networks?\n');
L4_lang_report = Level4_Results(Level4_Results.IsLanguageICN & strcmp(Level4_Results.Metric, 'IRi'), :);
for i = 1:height(L4_lang_report)
    sigStr = '';
    if L4_lang_report.Significant_FDR(i)
        sigStr = '*** FDR sig';
    end
    fprintf(fid, '  ICN%02d: Beta_C2=%.2f, R²_impr=%.3f, p=%.4f %s\n', ...
            L4_lang_report.ICN(i), L4_lang_report.Beta_C2(i), ...
            L4_lang_report.R2_improvement(i), L4_lang_report.p_C2(i), sigStr);
end

fprintf(fid, '\nH4b: R² improvement by network type (IRi)?\n');
fprintf(fid, '  Language ICNs: Mean R² improvement = %.4f\n', Level5_Results.Lang_MeanR2impr);
fprintf(fid, '  Non-language ICNs: Mean R² improvement = %.4f\n', Level5_Results.NonLang_MeanR2impr);
fprintf(fid, '\n');

% Level 5
fprintf(fid, 'LEVEL 5: NETWORK SPECIFICITY SUMMARY\n');
fprintf(fid, '-----------------------------------------------------------------\n');
fprintf(fid, 'Language vs Non-language IRi-WAB correlation comparison:\n');
fprintf(fid, '  Effect size (Cohen''s d): %.2f\n', Level5_Results.Specificity_d);
fprintf(fid, '  p-value: %.4f\n\n', Level5_Results.Specificity_p);

% IRi Sensitivity Analysis (NEW in V3)
fprintf(fid, 'IRi THRESHOLD SENSITIVITY ANALYSIS (NEW IN V3)\n');
fprintf(fid, '-----------------------------------------------------------------\n');
if ~isempty(Sensitivity_Results)
    fprintf(fid, 'Main threshold: %.1f\n', IRi_threshold);
    fprintf(fid, 'Tested: 0.7, 0.8, 0.9\n\n');
    
    fprintf(fid, 'ICN17 (Motor Network) Correlation Robustness:\n');
    for t = 1:height(Sensitivity_Results)
        sigStr = '';
        if Sensitivity_Results.ICN17_Significant(t)
            sigStr = '***';
        end
        fprintf(fid, '  Threshold %.1f: N=%d, rho=%.3f, p=%.4f %s\n', ...
            Sensitivity_Results.IRi_Threshold(t), Sensitivity_Results.N_Retained(t), ...
            Sensitivity_Results.ICN17_rho(t), Sensitivity_Results.ICN17_p(t), sigStr);
    end
    
    fprintf(fid, '\nConclusion: %s\n', sensitivityConclusion);
    if strcmp(sensitivityConclusion, 'ROBUST')
        fprintf(fid, '  Key finding (ICN17-WAB correlation) is robust to IRi threshold ±0.1\n');
    elseif strcmp(sensitivityConclusion, 'NOT_ROBUST')
        fprintf(fid, '  Caution: Key finding sensitivity depends on IRi threshold choice\n');
    end
else
    fprintf(fid, 'Sensitivity analysis not run (raw matrices file not found).\n');
end
fprintf(fid, '\n');

fprintf(fid, '=================================================================\n');
fprintf(fid, 'NETWORK REFERENCE (BrainMap20)\n');
fprintf(fid, '=================================================================\n');
for n = 1:nICN
    langStr = '';
    if ismember(n, languageICNs)
        langStr = ' *** LANGUAGE';
    elseif ismember(n, compensatoryICNs)
        langStr = ' (compensatory)';
    end
    fprintf(fid, 'ICN%02d: %s%s\n', n, trueICN_labels{n}, langStr);
end

fprintf(fid, '\n=================================================================\n');
fprintf(fid, 'CONCLUSIONS\n');
fprintf(fid, '=================================================================\n');

% Auto-generate conclusions based on results
fprintf(fid, '\nLevel 1 (Disease Effect):\n');
if sum(Level1_Results.Significant_FDR) == 0
    fprintf(fid, '  No significant group differences after FDR correction.\n');
    fprintf(fid, '  Mean effect size small (d=%.2f), suggesting similar ICN engagement\n', ...
            mean(abs(Level1_Results.Cohens_d(strcmp(Level1_Results.Metric, 'IRi')))));
    fprintf(fid, '  between aphasia and stroke controls during naming.\n');
else
    fprintf(fid, '  %d significant effects found (see details above).\n', ...
            sum(Level1_Results.Significant_FDR));
end

fprintf(fid, '\nLevel 2 (Subtype Effect):\n');
if sum(Level2_Results.Significant_FDR) == 0
    fprintf(fid, '  No significant subtype differences after FDR correction.\n');
    fprintf(fid, '  Aphasia subtypes show similar ICN engagement profiles.\n');
    if nEvidenceNull > height(Level2_Bayesian)/2
        fprintf(fid, '  Bayesian analysis supports equivalence (esp. Global vs others).\n');
    end
else
    fprintf(fid, '  %d significant effects found (see details above).\n', ...
            sum(Level2_Results.Significant_FDR));
end

fprintf(fid, '\nLevel 3 (Severity Effect):\n');
if sum(Level3_Results.Significant_FDR) > 0
    fprintf(fid, '  %d significant correlations between ICN engagement and WAB-AQ.\n', ...
            sum(Level3_Results.Significant_FDR));
    fprintf(fid, '  Higher ICN engagement associated with better language function.\n');
    
    % Check language network specificity
    nLangSig = sum(Level3_Results.Significant_FDR & Level3_Results.IsLanguageICN & ...
                   strcmp(Level3_Results.Metric, 'IRi'));
    if nLangSig > 0
        fprintf(fid, '  %d/%d a priori language networks showed significant IRi-WAB correlations.\n', ...
                nLangSig, length(languageICNs));
    end
else
    fprintf(fid, '  No significant correlations after FDR correction.\n');
end

fprintf(fid, '\nLevel 4 (Contrast Specificity):\n');
if sum(Level4_Results.Significant_FDR) > 0
    fprintf(fid, '  %d networks show incremental prediction from Naming>Abstract beyond Task>Rest.\n', ...
            sum(Level4_Results.Significant_FDR));
    fprintf(fid, '  Language-specific processing adds to severity prediction.\n');
else
    fprintf(fid, '  No significant incremental prediction from C2 beyond C1.\n');
    fprintf(fid, '  Domain-general activation (Task>Rest) may capture most variance.\n');
end

fprintf(fid, '\nSensitivity Analysis:\n');
if strcmp(sensitivityConclusion, 'ROBUST')
    fprintf(fid, '  Key findings robust to IRi threshold variation (0.7-0.9).\n');
elseif strcmp(sensitivityConclusion, 'NOT_ROBUST')
    fprintf(fid, '  CAUTION: Key findings sensitive to IRi threshold choice.\n');
else
    fprintf(fid, '  Sensitivity analysis %s.\n', sensitivityConclusion);
end

fprintf(fid, '\n=================================================================\n');
fprintf(fid, 'OUTPUT FILES\n');
fprintf(fid, '=================================================================\n');
fprintf(fid, 'ARC_03c_v3_Level1_DiseaseEffect.csv\n');
fprintf(fid, 'ARC_03c_v3_Level2_SubtypeEffect.csv\n');
fprintf(fid, 'ARC_03c_v3_Level2_Bayesian.csv (NEW)\n');
fprintf(fid, 'ARC_03c_v3_Level3_SeverityCorrelations.csv\n');
fprintf(fid, 'ARC_03c_v3_Level4_ContrastComparison.csv\n');
fprintf(fid, 'ARC_03c_v3_IRi_Sensitivity.csv (NEW)\n');
fprintf(fid, 'ARC_03c_v3_IRi_Sensitivity_LangICN.csv (NEW)\n');
fprintf(fid, 'ARC_03c_v3_AllResults.mat\n');
fprintf(fid, 'ARC_03c_v3_Summary.txt (this file)\n');

fprintf(fid, '\n=================================================================\n');
fprintf(fid, 'END OF REPORT\n');
fprintf(fid, '=================================================================\n');

fclose(fid);
fprintf('Saved: %s\n\n', reportFile);

%% =========================================================================
%% COMPLETION
%% =========================================================================
fprintf('=================================================================\n');
fprintf('ARC_03c_v3_Statistical_Analysis.m - COMPLETE\n');
fprintf('=================================================================\n');
fprintf('End time: %s\n', datestr(now));
fprintf('\nV3 CHANGES:\n');
fprintf('  - Input: ARC_03b_v3_Output (IRi QC threshold=%.1f)\n', IRi_threshold);
fprintf('  - NEW: Bayesian subtype analysis (BIC-based BF)\n');
fprintf('  - NEW: IRi threshold sensitivity analysis\n');
fprintf('  - REMOVED: Level 6 (LASSO) - circular analysis\n');
fprintf('\nCORRECTED METRICS USED:\n');
fprintf('  PRIMARY:     IRi (ICNiRelativeInvolvement)\n');
fprintf('  SECONDARY:   MANi (NormalisedMeanICNiActivation)\n');
fprintf('  EXPLORATORY: logVari (log variance of Z-scores)\n');
fprintf('\nHYPOTHESIS SUMMARY:\n');
fprintf('  H1 (Disease): %d/%d significant\n', ...
        sum(Level1_Results.Significant_FDR), nTests_L1);
fprintf('  H2 (Severity): %d/%d significant\n', ...
        sum(Level3_Results.Significant_FDR), nTests_L3);
fprintf('  H3 (Subtype): %d/%d significant (Frequentist)\n', ...
        sum(Level2_Results.Significant_FDR), nTests_L2);
fprintf('  H4 (Contrast): %d/%d significant\n', ...
        sum(Level4_Results.Significant_FDR), nTests_L4);
fprintf('\nBAYESIAN SUMMARY:\n');
fprintf('  Evidence for null: %d/%d comparisons\n', nEvidenceNull, height(Level2_Bayesian));
fprintf('  Evidence for alternative: %d/%d comparisons\n', nEvidenceAlt, height(Level2_Bayesian));
fprintf('\nSENSITIVITY ANALYSIS: %s\n', sensitivityConclusion);
fprintf('\nA PRIORI LANGUAGE NETWORKS:\n');
for i = 1:length(languageICNs)
    fprintf('  ICN%02d: %s\n', languageICNs(i), trueICN_labels{languageICNs(i)});
end
fprintf('\nOutputs saved to: %s\n', outputDir);
fprintf('\nOutput files:\n');
fprintf('  - Level 1: ARC_03c_v3_Level1_DiseaseEffect.csv\n');
fprintf('  - Level 2: ARC_03c_v3_Level2_SubtypeEffect.csv\n');
fprintf('  - Level 2 Bayesian: ARC_03c_v3_Level2_Bayesian.csv (NEW)\n');
fprintf('  - Level 3: ARC_03c_v3_Level3_SeverityCorrelations.csv\n');
fprintf('  - Level 4: ARC_03c_v3_Level4_ContrastComparison.csv\n');
fprintf('  - Sensitivity: ARC_03c_v3_IRi_Sensitivity.csv (NEW)\n');
fprintf('  - Sensitivity LangICN: ARC_03c_v3_IRi_Sensitivity_LangICN.csv (NEW)\n');
fprintf('  - All results: ARC_03c_v3_AllResults.mat\n');
fprintf('  - Summary: ARC_03c_v3_Summary.txt\n');
fprintf('\nNext steps:\n');
fprintf('  1. Review ARC_03c_v3_Summary.txt for hypothesis results\n');
fprintf('  2. Run Script 4_v3 for lesion-network analysis (includes mediation)\n');
fprintf('  3. Run Script 5_v3 for volume-stratified analysis\n');
fprintf('  4. Run Script 6_v3 for Smith10 robustness check\n');
fprintf('  5. Run Script 7 for integrated final report\n');
fprintf('=================================================================\n');