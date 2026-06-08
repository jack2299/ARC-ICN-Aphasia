%% CLEAN SAMPLE RE-RUN – LEVELS 1–4 (PORTABLE VERSION)
% Re‑runs Levels 1, 2, 3, 4 on the clean aphasia sample (rest‑only excluded).
% Uses identical methods, covariates, and FDR correction as the original pipeline.
% Computes logVari = log(Vari + 1) for all ICNs and contrasts.
%
% REQUIRED INPUT FILE (not included in repository – obtain from ARC):
%   ARC_03b_v3_Master_Wide.mat  – post‑IRi‑QC engagement master table
%                                  (must contain: PatientID, GroupRole, AphasiaType,
%                                   WAB_AQ, Age_At_Stroke, Sex, Days_Post_Stroke,
%                                   C1_ICN*_*, C2_ICN*_*, trueICN_labels)
%
% OUTPUT:
%   CleanSample_Results_Rigorous.mat  – Level1_Clean through Level4_Clean tables
%   CleanSample_Analysis_Report.txt   – human‑readable summary
% =========================================================================
clear; clc;

fprintf('=================================================================\n');
fprintf('CLEAN SAMPLE ANALYSIS (Excluding Rest‑Only Participants)\n');
fprintf('=================================================================\n');
fprintf('Start time: %s\n\n', datestr(now));

%% ---- 0. USER‑CONFIGURABLE PATHS ----
masterWideFile = 'ARC_03b_v3_Master_Wide.mat';
outputMatFile  = 'CleanSample_Results_Rigorous.mat';
outputReport   = 'CleanSample_Analysis_Report.txt';

%% ---- 1. LOAD DATA ----
if ~exist(masterWideFile, 'file')
    error('Master file not found: %s', masterWideFile);
end
load(masterWideFile, 'masterWide', 'trueICN_labels');

% Verify critical variables
requiredCols = {'PatientID', 'GroupRole', 'AphasiaType', 'WAB_AQ', ...
                'Age_At_Stroke', 'Sex', 'Days_Post_Stroke', 'C2_ICN17_IRi'};
missingCols = setdiff(requiredCols, masterWide.Properties.VariableNames);
if ~isempty(missingCols)
    error('Missing columns in masterWide: %s', strjoin(missingCols, ', '));
end

%% ---- 2. CREATE DERIVED VARIABLES (matching original pipeline) ----

% SexNumeric
if ~ismember('SexNumeric', masterWide.Properties.VariableNames)
    masterWide.SexNumeric = double(strcmp(masterWide.Sex, 'M'));
end

% logVari = log(Vari + 1) for C1 and C2
fprintf('Computing logVari = log(Vari + 1) for all ICNs...\n');
for n = 1:18
    for contrast = {'C1', 'C2'}
        c = contrast{1};
        variCol    = sprintf('%s_ICN%02d_Vari', c, n);
        logVariCol = sprintf('%s_ICN%02d_logVari', c, n);
        if ismember(variCol, masterWide.Properties.VariableNames) && ...
           ~ismember(logVariCol, masterWide.Properties.VariableNames)
            masterWide.(logVariCol) = log(masterWide.(variCol) + 1);
        end
    end
end
fprintf('logVari columns added.\n\n');

%% ---- 3. DEFINE CLEAN SAMPLE ----

% Rest‑only participants (anonymized ARC IDs – safe for public repository)
restOnlyList = {'M2088','M2097','M2100','M2101','M2113','M2114',...
    'M2117','M2118','M2122','M2126','M2129','M2131','M2135',...
    'M2140','M2141','M2142','M2143','M2144','M2145','M2146',...
    'M2149','M2150','M2151','M2152','M2153','M2155','M2156',...
    'M2158','M2159','M2160','M2162','M2164','M2165','M2169',...
    'M2184','M2254'};

isRestOnly = ismember(masterWide.PatientID, restOnlyList);
isAphasia  = strcmp(masterWide.GroupRole, 'Aphasia');
isControl  = strcmp(masterWide.GroupRole, 'Stroke Control (No Aphasia)');
hasValidC2 = ~isnan(masterWide.C2_ICN17_IRi);

cleanAphasia = masterWide(isAphasia  & hasValidC2 & ~isRestOnly, :);
cleanControl = masterWide(isControl  & hasValidC2 & ~isRestOnly, :);

fprintf('CLEAN SAMPLE:\n');
fprintf('  Aphasia:  N = %d (original: 153)\n', height(cleanAphasia));
fprintf('  Controls: N = %d (original:  21)\n', height(cleanControl));
fprintf('  Excluded rest‑only aphasia:  14\n');
fprintf('  Excluded rest‑only controls:  7\n\n');

%% ---- 4. PARAMETERS ----
fdrQ           = 0.05;
metrics        = {'IRi', 'MANi', 'logVari'};
primaryMetric  = 'IRi';
nICN           = 18;
nMetrics       = length(metrics);
nTests         = nICN * nMetrics;   % 54 tests per level
languageICNs   = [4, 16, 18];
compensatoryICNs = [15];

%% ---- 5. HELPER FUNCTIONS (same as original) ----

function [h, crit_p, adj_p] = fdr_bh(pvals, q)
    if nargin < 2, q = 0.05; end
    pvals = pvals(:);
    m = length(pvals);
    if all(isnan(pvals))
        h = false(size(pvals)); crit_p = 0; adj_p = NaN(size(pvals)); return;
    end
    [sorted_p, sort_idx] = sort(pvals);
    thresh = (1:m)' / m * q;
    below = sorted_p <= thresh & ~isnan(sorted_p);
    if any(below), crit_p = max(sorted_p(below)); else, crit_p = 0; end
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

function d = cohens_d(x1, x2)
    x1 = x1(~isnan(x1)); x2 = x2(~isnan(x2));
    n1 = length(x1); n2 = length(x2);
    if n1 < 2 || n2 < 2, d = NaN; return; end
    pooled_std = sqrt(((n1-1)*var(x1) + (n2-1)*var(x2)) / (n1+n2-2));
    if pooled_std == 0, d = 0; else, d = (mean(x1) - mean(x2)) / pooled_std; end
end

function [rho, pval, n_valid] = partial_corr_spearman(x, y, covariates)
    validIdx = ~isnan(x) & ~isnan(y) & ~any(isnan(covariates), 2);
    x = x(validIdx); y = y(validIdx); covariates = covariates(validIdx, :);
    n_valid = length(x);
    if n_valid < 10, rho = NaN; pval = NaN; return; end
    x_rank = tiedrank(x); y_rank = tiedrank(y);
    X_cov = [ones(n_valid,1), covariates];
    resid_x = x_rank - X_cov * (X_cov \ x_rank);
    resid_y = y_rank - X_cov * (X_cov \ y_rank);
    [rho, pval] = corr(resid_x, resid_y, 'Type', 'Pearson');
end

function str = ternary(cond, trueStr, falseStr)
    if cond, str = trueStr; else, str = falseStr; end
end

%% ========================================================================
%% LEVEL 1: DISEASE EFFECT (Aphasia vs Control)
%% ========================================================================
fprintf('=================================================================\n');
fprintf('LEVEL 1: Disease Effect (Aphasia vs Stroke Control)\n');
fprintf('         ANCOVA controlling for Age, Sex, Days_Post_Stroke\n');
fprintf('=================================================================\n');

Level1_Clean = table('Size', [nTests, 20], ...
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
for m = 1:nMetrics
    metric = metrics{m};
    for n = 1:nICN
        resultIdx = resultIdx + 1;
        colName = sprintf('C2_ICN%02d_%s', n, metric);
        
        y_aph = cleanAphasia.(colName);
        y_con = cleanControl.(colName);
        age_aph = cleanAphasia.Age_At_Stroke; age_con = cleanControl.Age_At_Stroke;
        sex_aph = cleanAphasia.SexNumeric;    sex_con = cleanControl.SexNumeric;
        days_aph = cleanAphasia.Days_Post_Stroke; days_con = cleanControl.Days_Post_Stroke;
        
        y_valid = [y_aph; y_con];
        group_valid = [ones(height(cleanAphasia),1); zeros(height(cleanControl),1)];
        age_valid = [age_aph; age_con];
        sex_valid = [sex_aph; sex_con];
        days_valid = [days_aph; days_con];
        
        validIdx = ~isnan(y_valid) & ~isnan(age_valid) & ~isnan(sex_valid) & ~isnan(days_valid);
        y_valid = y_valid(validIdx); group_valid = group_valid(validIdx);
        age_valid = age_valid(validIdx); sex_valid = sex_valid(validIdx); days_valid = days_valid(validIdx);
        
        aphasiaData = y_valid(group_valid == 1);
        controlData = y_valid(group_valid == 0);
        
        X = [ones(length(y_valid),1), group_valid, age_valid, sex_valid, days_valid];
        X_reduced = [ones(length(y_valid),1), age_valid, sex_valid, days_valid];
        
        [~,~,~,~,stats_full]    = regress(y_valid, X);
        [~,~,~,~,stats_reduced] = regress(y_valid, X_reduced);
        
        SSE_full = stats_full(4) * (length(y_valid) - size(X,2));
        SSE_reduced = stats_reduced(4) * (length(y_valid) - size(X_reduced,2));
        SS_group = SSE_reduced - SSE_full;
        df1 = 1; df2 = length(y_valid) - size(X,2);
        MSE_full = SSE_full / df2;
        
        if MSE_full > 0
            F_stat = (SS_group/df1) / MSE_full;
            p_value = 1 - fcdf(F_stat, df1, df2);
        else
            F_stat = 0; p_value = 1;
        end
        
        d = cohens_d(aphasiaData, controlData);
        n1 = length(aphasiaData); n2 = length(controlData);
        se_d = sqrt((n1+n2)/(n1*n2) + d^2/(2*(n1+n2)));
        d_ci_lo = d - 1.96*se_d; d_ci_hi = d + 1.96*se_d;
        
        Level1_Clean.ICN(resultIdx) = n;
        Level1_Clean.ICN_Label{resultIdx} = trueICN_labels{n};
        Level1_Clean.Metric{resultIdx} = metric;
        Level1_Clean.IsLanguageICN(resultIdx) = ismember(n, languageICNs);
        Level1_Clean.IsCompensatoryICN(resultIdx) = ismember(n, compensatoryICNs);
        Level1_Clean.Aphasia_Mean(resultIdx) = mean(aphasiaData);
        Level1_Clean.Aphasia_SD(resultIdx) = std(aphasiaData);
        Level1_Clean.Aphasia_N(resultIdx) = length(aphasiaData);
        Level1_Clean.Control_Mean(resultIdx) = mean(controlData);
        Level1_Clean.Control_SD(resultIdx) = std(controlData);
        Level1_Clean.Control_N(resultIdx) = length(controlData);
        Level1_Clean.F_stat(resultIdx) = F_stat;
        Level1_Clean.df1(resultIdx) = df1; Level1_Clean.df2(resultIdx) = df2;
        Level1_Clean.p_uncorrected(resultIdx) = p_value;
        Level1_Clean.Cohens_d(resultIdx) = d;
        Level1_Clean.Cohens_d_CI_lo(resultIdx) = d_ci_lo;
        Level1_Clean.Cohens_d_CI_hi(resultIdx) = d_ci_hi;
    end
end

[Level1_Clean.Significant_FDR, ~, Level1_Clean.p_FDR] = fdr_bh(Level1_Clean.p_uncorrected, fdrQ);
Level1_Clean = sortrows(Level1_Clean, 'p_uncorrected');

fprintf('Significant after FDR: %d / %d tests\n', sum(Level1_Clean.Significant_FDR), nTests);
icn17_l1 = Level1_Clean(Level1_Clean.ICN == 17 & strcmp(Level1_Clean.Metric, 'IRi'), :);
if ~isempty(icn17_l1)
    fprintf('ICN17 IRi: d = %.2f [%.2f, %.2f], p_FDR = %.4f\n\n', ...
        icn17_l1.Cohens_d, icn17_l1.Cohens_d_CI_lo, icn17_l1.Cohens_d_CI_hi, icn17_l1.p_FDR);
end

%% ========================================================================
%% LEVEL 2: SUBTYPE EFFECT (4‑group ANCOVA)
%% ========================================================================
fprintf('=================================================================\n');
fprintf('LEVEL 2: Subtype Effect (Anomic, Broca, Conduction, Global)\n');
fprintf('         ANCOVA controlling for Age, Sex, Days_Post_Stroke\n');
fprintf('=================================================================\n');

cleanSubtype = masterWide(isAphasia & hasValidC2 & ~isRestOnly & ...
    ~strcmp(masterWide.AphasiaType, 'Wernicke'), :);
fprintf('Clean subtype sample: N = %d (Wernicke excluded)\n\n', height(cleanSubtype));

Level2_Clean = table('Size', [nTests, 17], ...
    'VariableTypes', {'double', 'cell', 'cell', 'logical', ...
                      'double', 'double', 'double', 'double', 'double', ...
                      'double', 'double', 'double', 'double', ...
                      'double', 'double', 'logical', 'double'}, ...
    'VariableNames', {'ICN', 'ICN_Label', 'Metric', 'IsLanguageICN', ...
                      'F_stat', 'df_between', 'df_within', 'p_uncorrected', 'Eta_squared', ...
                      'Mean_Anomic', 'Mean_Broca', 'Mean_Conduction', 'Mean_Global', ...
                      'N_valid', 'MaxDiff', 'Significant_FDR', 'p_FDR'});

resultIdx = 0;
for m = 1:nMetrics
    metric = metrics{m};
    for n = 1:nICN
        resultIdx = resultIdx + 1;
        colName = sprintf('C2_ICN%02d_%s', n, metric);
        
        y = cleanSubtype.(colName);
        age = cleanSubtype.Age_At_Stroke; sex = cleanSubtype.SexNumeric;
        days = cleanSubtype.Days_Post_Stroke; groups = cleanSubtype.AphasiaType;
        
        validIdx = ~isnan(y) & ~isnan(age) & ~isnan(sex) & ~isnan(days);
        y_valid = y(validIdx); age_valid = age(validIdx);
        sex_valid = sex(validIdx); days_valid = days(validIdx);
        groups_valid = groups(validIdx);
        
        dummyBroca = double(strcmp(groups_valid, 'Broca'));
        dummyConduction = double(strcmp(groups_valid, 'Conduction'));
        dummyGlobal = double(strcmp(groups_valid, 'Global'));
        
        X_full = [ones(length(y_valid),1), dummyBroca, dummyConduction, dummyGlobal, ...
                  age_valid, sex_valid, days_valid];
        X_reduced = [ones(length(y_valid),1), age_valid, sex_valid, days_valid];
        
        [~,~,~,~,stats_full]    = regress(y_valid, X_full);
        [~,~,~,~,stats_reduced] = regress(y_valid, X_reduced);
        
        SSE_full = stats_full(4) * (length(y_valid) - size(X_full,2));
        SSE_reduced = stats_reduced(4) * (length(y_valid) - size(X_reduced,2));
        SS_group = SSE_reduced - SSE_full;
        df_between = 3; df_within = length(y_valid) - size(X_full,2);
        MSE_full = SSE_full / df_within;
        
        if MSE_full > 0
            F_stat = (SS_group/df_between) / MSE_full;
            p_value = 1 - fcdf(F_stat, df_between, df_within);
        else
            F_stat = 0; p_value = 1;
        end
        
        SS_total = var(y_valid) * (length(y_valid) - 1);
        eta_sq = max(0, SS_group / SS_total);
        
        meanAnomic = mean(y_valid(strcmp(groups_valid,'Anomic')), 'omitnan');
        meanBroca  = mean(y_valid(strcmp(groups_valid,'Broca')), 'omitnan');
        meanConduction = mean(y_valid(strcmp(groups_valid,'Conduction')), 'omitnan');
        meanGlobal = mean(y_valid(strcmp(groups_valid,'Global')), 'omitnan');
        maxDiff = max([meanAnomic, meanBroca, meanConduction, meanGlobal]) - ...
                  min([meanAnomic, meanBroca, meanConduction, meanGlobal]);
        
        Level2_Clean.ICN(resultIdx) = n;
        Level2_Clean.ICN_Label{resultIdx} = trueICN_labels{n};
        Level2_Clean.Metric{resultIdx} = metric;
        Level2_Clean.IsLanguageICN(resultIdx) = ismember(n, languageICNs);
        Level2_Clean.F_stat(resultIdx) = F_stat;
        Level2_Clean.df_between(resultIdx) = df_between;
        Level2_Clean.df_within(resultIdx) = df_within;
        Level2_Clean.p_uncorrected(resultIdx) = p_value;
        Level2_Clean.Eta_squared(resultIdx) = eta_sq;
        Level2_Clean.Mean_Anomic(resultIdx) = meanAnomic;
        Level2_Clean.Mean_Broca(resultIdx) = meanBroca;
        Level2_Clean.Mean_Conduction(resultIdx) = meanConduction;
        Level2_Clean.Mean_Global(resultIdx) = meanGlobal;
        Level2_Clean.N_valid(resultIdx) = length(y_valid);
        Level2_Clean.MaxDiff(resultIdx) = maxDiff;
    end
end

[Level2_Clean.Significant_FDR, ~, Level2_Clean.p_FDR] = fdr_bh(Level2_Clean.p_uncorrected, fdrQ);
Level2_Clean = sortrows(Level2_Clean, 'p_uncorrected');
fprintf('Significant after FDR: %d / %d tests\n\n', sum(Level2_Clean.Significant_FDR), nTests);

%% ========================================================================
%% LEVEL 3: SEVERITY EFFECT (Partial Spearman correlations)
%% ========================================================================
fprintf('=================================================================\n');
fprintf('LEVEL 3: Severity Effect (Partial Spearman vs WAB‑R)\n');
fprintf('         Controlling for Age, Sex, Days_Post_Stroke\n');
fprintf('=================================================================\n');

Level3_Clean = table('Size', [nTests, 14], ...
    'VariableTypes', {'double', 'cell', 'cell', 'logical', ...
                      'double', 'double', 'double', 'double', 'double', ...
                      'double', 'double', 'double', 'logical', 'double'}, ...
    'VariableNames', {'ICN', 'ICN_Label', 'Metric', 'IsLanguageICN', ...
                      'N', 'Partial_rho', 'p_uncorrected', 'Abs_rho', ...
                      'Rho_CI_lo', 'Rho_CI_hi', 'Zero_order_rho', ...
                      'Direction', 'Significant_FDR', 'p_FDR'});

resultIdx = 0;
for m = 1:nMetrics
    metric = metrics{m};
    for n = 1:nICN
        resultIdx = resultIdx + 1;
        colName = sprintf('C2_ICN%02d_%s', n, metric);
        
        wab = cleanAphasia.WAB_AQ; icn = cleanAphasia.(colName);
        age = cleanAphasia.Age_At_Stroke; sex = cleanAphasia.SexNumeric;
        days = cleanAphasia.Days_Post_Stroke;
        covariates = [age, sex, days];
        
        [rho_partial, p_partial, n_valid] = partial_corr_spearman(wab, icn, covariates);
        
        if ~isnan(rho_partial) && n_valid > 3
            z = 0.5 * log((1+rho_partial)/(1-rho_partial));
            se_z = 1/sqrt(n_valid-3);
            rho_ci_lo = tanh(z - 1.96*se_z);
            rho_ci_hi = tanh(z + 1.96*se_z);
        else
            rho_ci_lo = NaN; rho_ci_hi = NaN;
        end
        
        validIdx = ~isnan(wab) & ~isnan(icn);
        if sum(validIdx) >= 10
            rho_zero = corr(wab(validIdx), icn(validIdx), 'Type', 'Spearman');
        else
            rho_zero = NaN;
        end
        
        Level3_Clean.ICN(resultIdx) = n;
        Level3_Clean.ICN_Label{resultIdx} = trueICN_labels{n};
        Level3_Clean.Metric{resultIdx} = metric;
        Level3_Clean.IsLanguageICN(resultIdx) = ismember(n, languageICNs);
        Level3_Clean.N(resultIdx) = n_valid;
        Level3_Clean.Partial_rho(resultIdx) = rho_partial;
        Level3_Clean.p_uncorrected(resultIdx) = p_partial;
        Level3_Clean.Abs_rho(resultIdx) = abs(rho_partial);
        Level3_Clean.Rho_CI_lo(resultIdx) = rho_ci_lo;
        Level3_Clean.Rho_CI_hi(resultIdx) = rho_ci_hi;
        Level3_Clean.Zero_order_rho(resultIdx) = rho_zero;
        Level3_Clean.Direction(resultIdx) = sign(rho_partial);
    end
end

[Level3_Clean.Significant_FDR, ~, Level3_Clean.p_FDR] = fdr_bh(Level3_Clean.p_uncorrected, fdrQ);
Level3_Clean = sortrows(Level3_Clean, 'Abs_rho', 'descend');

fprintf('Significant after FDR: %d / %d tests\n', sum(Level3_Clean.Significant_FDR), nTests);
icn17_l3 = Level3_Clean(Level3_Clean.ICN == 17 & strcmp(Level3_Clean.Metric, 'IRi'), :);
if ~isempty(icn17_l3)
    fprintf('ICN17 IRi (Primary): ρ = %.3f [%.3f, %.3f], p_FDR = %.4f\n\n', ...
        icn17_l3.Partial_rho, icn17_l3.Rho_CI_lo, icn17_l3.Rho_CI_hi, icn17_l3.p_FDR);
end

%% ========================================================================
%% LEVEL 4: CONTRAST SPECIFICITY (C2 vs C1)
%% ========================================================================
fprintf('=================================================================\n');
fprintf('LEVEL 4: Contrast Specificity (Does C2 add beyond C1?)\n');
fprintf('         Hierarchical regression: WAB ~ C2 + C1\n');
fprintf('=================================================================\n');

hasBoth = ~isnan(cleanAphasia.C1_ICN01_IRi) & ~isnan(cleanAphasia.C2_ICN01_IRi);
contrastData = cleanAphasia(hasBoth, :);
fprintf('Clean sample with both contrasts: N = %d\n\n', height(contrastData));

Level4_Clean = table('Size', [nTests, 16], ...
    'VariableTypes', {'double', 'cell', 'cell', 'logical', 'double', ...
                      'double', 'double', 'double', 'double', 'double', 'double', ...
                      'double', 'double', 'double', 'logical', 'double'}, ...
    'VariableNames', {'ICN', 'ICN_Label', 'Metric', 'IsLanguageICN', 'N', ...
                      'Beta_C2', 't_C2', 'p_C2', 'Beta_C1', 't_C1', 'p_C1', ...
                      'R2_full', 'R2_C1only', 'R2_improvement', 'Significant_FDR', 'p_FDR'});

resultIdx = 0;
for m = 1:nMetrics
    metric = metrics{m};
    for n = 1:nICN
        resultIdx = resultIdx + 1;
        colC1 = sprintf('C1_ICN%02d_%s', n, metric);
        colC2 = sprintf('C2_ICN%02d_%s', n, metric);
        
        wab = contrastData.WAB_AQ; c1 = contrastData.(colC1); c2 = contrastData.(colC2);
        validIdx = ~isnan(wab) & ~isnan(c1) & ~isnan(c2);
        wab_valid = wab(validIdx); c1_valid = c1(validIdx); c2_valid = c2(validIdx);
        n_obs = sum(validIdx);
        
        if n_obs >= 10
            tbl = table(wab_valid, c2_valid, c1_valid, 'VariableNames', {'WAB','C2','C1'});
            mdl = fitlm(tbl, 'WAB ~ C2 + C1');
            beta_C2 = mdl.Coefficients.Estimate(2); t_C2 = mdl.Coefficients.tStat(2);
            p_C2 = mdl.Coefficients.pValue(2);
            beta_C1 = mdl.Coefficients.Estimate(3); t_C1 = mdl.Coefficients.tStat(3);
            p_C1 = mdl.Coefficients.pValue(3);
            R2_full = mdl.Rsquared.Ordinary;
            mdl_C1 = fitlm(tbl, 'WAB ~ C1');
            R2_impr = R2_full - mdl_C1.Rsquared.Ordinary;
        else
            beta_C2 = NaN; t_C2 = NaN; p_C2 = NaN;
            beta_C1 = NaN; t_C1 = NaN; p_C1 = NaN;
            R2_full = NaN; R2_impr = NaN;
        end
        
        Level4_Clean.ICN(resultIdx) = n;
        Level4_Clean.ICN_Label{resultIdx} = trueICN_labels{n};
        Level4_Clean.Metric{resultIdx} = metric;
        Level4_Clean.IsLanguageICN(resultIdx) = ismember(n, languageICNs);
        Level4_Clean.N(resultIdx) = n_obs;
        Level4_Clean.Beta_C2(resultIdx) = beta_C2;
        Level4_Clean.t_C2(resultIdx) = t_C2; Level4_Clean.p_C2(resultIdx) = p_C2;
        Level4_Clean.Beta_C1(resultIdx) = beta_C1;
        Level4_Clean.t_C1(resultIdx) = t_C1; Level4_Clean.p_C1(resultIdx) = p_C1;
        Level4_Clean.R2_full(resultIdx) = R2_full;
        Level4_Clean.R2_C1only(resultIdx) = NaN;  % not stored separately
        Level4_Clean.R2_improvement(resultIdx) = R2_impr;
    end
end

[Level4_Clean.Significant_FDR, ~, Level4_Clean.p_FDR] = fdr_bh(Level4_Clean.p_C2, fdrQ);
Level4_Clean = sortrows(Level4_Clean, 'p_C2');

fprintf('Networks where C2 adds beyond C1 (FDR): %d / %d\n', sum(Level4_Clean.Significant_FDR), nTests);
icn17_l4 = Level4_Clean(Level4_Clean.ICN == 17 & strcmp(Level4_Clean.Metric, 'IRi'), :);
if ~isempty(icn17_l4)
    fprintf('ICN17 IRi: ΔR² = %.3f (%.1f%%), p_FDR = %.4f\n\n', ...
        icn17_l4.R2_improvement, icn17_l4.R2_improvement*100, icn17_l4.p_FDR);
end

%% ========================================================================
%% GENERATE TEXT REPORT
%% ========================================================================
fprintf('Generating text report...\n');
fid = fopen(outputReport, 'w');

fprintf(fid, '========================================\n');
fprintf(fid, 'CLEAN SAMPLE ANALYSIS REPORT\n');
fprintf(fid, 'Generated: %s\n', datestr(now));
fprintf(fid, '========================================\n\n');

fprintf(fid, 'SAMPLE\n');
fprintf(fid, '----------------------------------------\n');
fprintf(fid, 'Original aphasia with valid C2:     %d\n', 153);
fprintf(fid, 'Excluded (rest‑only aphasia):        %d\n', 14);
fprintf(fid, 'Clean aphasia sample:                %d\n', height(cleanAphasia));
fprintf(fid, 'Original controls with valid C2:     %d\n', 21);
fprintf(fid, 'Excluded (rest‑only controls):       %d\n', 7);
fprintf(fid, 'Clean control sample:                %d\n\n', height(cleanControl));

fprintf(fid, 'LEVEL 1: DISEASE EFFECT\n');
fprintf(fid, 'Tests: %d | FDR sig: %d\n\n', nTests, sum(Level1_Clean.Significant_FDR));

fprintf(fid, 'LEVEL 2: SUBTYPE EFFECT\n');
fprintf(fid, 'Clean subtype N: %d | FDR sig: %d\n\n', height(cleanSubtype), sum(Level2_Clean.Significant_FDR));

fprintf(fid, 'LEVEL 3: SEVERITY EFFECT (Primary)\n');
fprintf(fid, 'ICN17 IRi: ρ = %.3f [%.3f, %.3f], p_FDR = %.4f, N = %d\n\n', ...
    icn17_l3.Partial_rho, icn17_l3.Rho_CI_lo, icn17_l3.Rho_CI_hi, icn17_l3.p_FDR, icn17_l3.N);

fprintf(fid, 'LEVEL 4: CONTRAST SPECIFICITY\n');
fprintf(fid, 'ICN17 IRi: ΔR² = %.1f%%, p_FDR = %.4f\n\n', icn17_l4.R2_improvement*100, icn17_l4.p_FDR);

fprintf(fid, '========================================\n');
fprintf(fid, 'END OF REPORT\n');
fprintf(fid, '========================================\n');
fclose(fid);
fprintf('Report saved to: %s\n', outputReport);

%% ========================================================================
%% SAVE RESULTS
%% ========================================================================
fprintf('\nSaving results...\n');
save(outputMatFile, 'Level1_Clean', 'Level2_Clean', 'Level3_Clean', 'Level4_Clean');
fprintf('Saved: %s\n', outputMatFile);
fprintf('=================================================================\n');
fprintf('COMPLETE\n');
fprintf('=================================================================\n');
