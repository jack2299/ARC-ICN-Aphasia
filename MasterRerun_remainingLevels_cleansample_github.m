%% QC CORRECTIVE SCRIPT – ALL ENGAGEMENT‑BASED RE‑ANALYSES ON CLEAN SAMPLE
% Re‑runs every analysis that was contaminated by the inclusion of rest‑only
% participants with spurious engagement values. Uses the CLEAN aphasia sample
% (rest‑only excluded, valid C2 IRi).
%
% Included analyses:
%   - Level 8  (disconnection: damage‑engagement Spearman) -> QC_Corrected_Level8_[metric].csv
%   - Level 8b (specificity, IRi only)                     -> QC_Corrected_Level8b.csv
%   - Level 10b (bootstrap mediation, ICN17, all metrics)  -> QC_Corrected_Level10b.csv
%   - Volume‑stratified engagement meta‑analysis (all ICNs, 3 metrics)
%                                                           -> QC_Corrected_VolumeStrat_Engagement_[metric].csv
%   - Cross‑atlas validation (Smith10 S06)                 -> QC_Corrected_Smith10_S06.csv
%   - IRi threshold sensitivity (ICN17)                    -> QC_Corrected_IRi_Sensitivity.csv
%   - Comparison report                                    -> QC_Corrected_Comparison_Report.txt
%
% REQUIRED INPUT FILES (not included in repository – obtain from ARC):
%   ARC_03b_v3_Master_Wide.mat   – post‑IRi‑QC engagement master table
%   ARC_04_v3_AllResults.mat     – resampleLog & mergedTable (lesion overlap)
%   ARC_06_v3_AllResults.mat     – Smith10_Wide (cross‑atlas validation)
%
% All outputs are saved in the current folder.
% =========================================================================
clearvars; clc;

fprintf('=================================================================\n');
fprintf('QC CORRECTIVE SCRIPT – ALL ENGAGEMENT‑BASED RE‑ANALYSES (CLEAN)\n');
fprintf('=================================================================\n');
fprintf('Started: %s\n\n', datestr(now));

%% ---- 0. USER‑CONFIGURABLE PATHS ----
% Input files (must exist)
masterWideFile   = 'ARC_03b_v3_Master_Wide.mat';
lesionFile       = 'ARC_04_v3_AllResults.mat';
smith10File      = 'ARC_06_v3_AllResults.mat';

% Output file prefixes (CSV and TXT will be written to current folder)
outPrefix = 'QC_Corrected';

% Rest‑only participant list (anonymized ARC IDs – safe for public repository)
restOnlyFull = {'M2088','M2097','M2100','M2101','M2113','M2114','M2117','M2118','M2122','M2126','M2129','M2131','M2135',...
    'M2140','M2141','M2142','M2143','M2144','M2145','M2146','M2149','M2150','M2151','M2152','M2153','M2155','M2156',...
    'M2158','M2159','M2160','M2162','M2164','M2165','M2169','M2184','M2254','M2307'};

%% ---- 1. LOAD DATA & VERIFY ----
if ~exist(masterWideFile, 'file')
    error('Master file not found: %s', masterWideFile);
end
load(masterWideFile, 'masterWide');

if ~exist(lesionFile, 'file')
    error('Lesion results file not found: %s', lesionFile);
end
load(lesionFile, 'resampleLog', 'mergedTable');

% Cross‑atlas data is optional – script will skip if missing
hasSmith10 = exist(smith10File, 'file');

%% ---- 2. DEFINE CLEAN SAMPLES ----
isAph = strcmp(masterWide.GroupRole, 'Aphasia');
isRestOnly = ismember(masterWide.PatientID, restOnlyFull);
hasValidC2 = ~isnan(masterWide.C2_ICN17_IRi);

% Ensure SexNumeric exists (needed for regressions)
if ~ismember('SexNumeric', masterWide.Properties.VariableNames)
    masterWide.SexNumeric = double(strcmp(masterWide.Sex, 'M'));
end

% logVari creation (log(Vari + 1))
fprintf('Preparing logVari columns...\n');
for n = 1:18
    c2_vari = sprintf('C2_ICN%02d_Vari', n);
    c2_logvari = sprintf('C2_ICN%02d_logVari', n);
    if ismember(c2_vari, masterWide.Properties.VariableNames)
        masterWide.(c2_logvari) = log(masterWide.(c2_vari) + 1);
    end
end
fprintf('Done.\n');

% Clean aphasia samples
cleanAph_full = masterWide(isAph & ~isRestOnly & hasValidC2, :);
lesionIDs = unique(resampleLog.PatientID(~cellfun('isempty', resampleLog.PatientID)));
hasLesion = ismember(cleanAph_full.PatientID, lesionIDs);
cleanAph_lesion = cleanAph_full(hasLesion, :);
fprintf('Clean aphasia (full, N=139) and with lesion mask (N=%d)\n\n', height(cleanAph_lesion));

% Match cleanAph_lesion to mergedTable for damage columns
[~, idxM, idxC] = intersect(mergedTable.PatientID, cleanAph_lesion.PatientID, 'stable');
damage_mat = mergedTable{idxM, 3:20};   % ICN01_Damage ... ICN18_Damage

%% ---- 3. ANALYSIS PARAMETERS ----
metrics = {'IRi', 'MANi', 'logVari'};
nICN = 18;
keyICN = 17;
nBoot = 5000;
alpha = 0.05;
fdrQ = 0.05;

% Quartile boundaries for volume stratification (from manuscript)
Q1_max = 91556; Q2_max = 173496; Q3_max = 280746;

%% ---- 4. LEVEL 8 – DISCONNECTION (Spearman damage‑engagement) ----
fprintf('=================================================================\n');
fprintf('LEVEL 8: DISCONNECTION (CLEAN)\n');
fprintf('=================================================================\n');

for m = 1:length(metrics)
    metric = metrics{m};
    L8 = table((1:nICN)', 'VariableNames', {'ICN'});
    L8.Spearman_rho = NaN(nICN,1);
    L8.p_uncorrected = NaN(nICN,1);
    L8.N = NaN(nICN,1);
    
    for n = 1:nICN
        engCol = sprintf('C2_ICN%02d_%s', n, metric);
        eng = cleanAph_lesion.(engCol)(idxC);
        damage = damage_mat(:, n);
        valid = ~isnan(damage) & ~isnan(eng);
        if sum(valid) >= 10
            [rho, p] = corr(damage(valid), eng(valid), 'Type', 'Spearman');
            L8.Spearman_rho(n) = rho;
            L8.p_uncorrected(n) = p;
            L8.N(n) = sum(valid);
        end
    end
    [h, ~, p_fdr] = fdr_bh(L8.p_uncorrected, fdrQ);
    L8.Significant_FDR = h;
    L8.p_FDR = p_fdr;
    writetable(L8, sprintf('%s_Level8_%s.csv', outPrefix, metric));
    fprintf('Level 8 %s saved. Mean rho: %.3f\n', metric, mean(L8.Spearman_rho, 'omitnan'));
end

%% ---- 5. LEVEL 8b – SPECIFICITY (IRi only) ----
fprintf('\n=================================================================\n');
fprintf('LEVEL 8b: SPECIFICITY (IRi, CLEAN)\n');
fprintf('=================================================================\n');

L8b = table((1:nICN)', 'VariableNames', {'ICN'});
L8b.Rho_Within = NaN(nICN,1);
L8b.Rho_Between_Mean = NaN(nICN,1);
L8b.Specificity = NaN(nICN,1);

for n = 1:nICN
    eng_within = cleanAph_lesion.(sprintf('C2_ICN%02d_IRi', n))(idxC);
    damage = damage_mat(:, n);
    valid = ~isnan(damage) & ~isnan(eng_within);
    if sum(valid) >= 10
        rho_within = corr(damage(valid), eng_within(valid), 'Type', 'Spearman');
    else
        rho_within = NaN;
    end
    
    rho_between = [];
    for other = setdiff(1:nICN, n)
        eng_other = cleanAph_lesion.(sprintf('C2_ICN%02d_IRi', other))(idxC);
        valid_b = ~isnan(damage) & ~isnan(eng_other);
        if sum(valid_b) >= 10
            rho_between(end+1) = corr(damage(valid_b), eng_other(valid_b), 'Type', 'Spearman');
        end
    end
    L8b.Rho_Within(n) = rho_within;
    L8b.Rho_Between_Mean(n) = mean(rho_between, 'omitnan');
    L8b.Specificity(n) = rho_within - L8b.Rho_Between_Mean(n);
end
writetable(L8b, sprintf('%s_Level8b.csv', outPrefix));
fprintf('Level 8b saved. Mean within: %.3f, mean between: %.3f, mean specificity: %.3f\n\n', ...
    mean(L8b.Rho_Within, 'omitnan'), mean(L8b.Rho_Between_Mean, 'omitnan'), mean(L8b.Specificity, 'omitnan'));

%% ---- 6. LEVEL 10b – BOOTSTRAP MEDIATION (ICN17) ----
fprintf('=================================================================\n');
fprintf('LEVEL 10b: Bootstrap Mediation (ICN17)\n');
fprintf('=================================================================\n');

damage_clean = damage_mat(:, keyICN);
Level10b_Corrected = table();

for m = 1:length(metrics)
    metric = metrics{m};
    engCol = sprintf('C2_ICN%02d_%s', keyICN, metric);
    eng_vals = cleanAph_lesion.(engCol)(idxC);
    
    valid = ~isnan(damage_clean) & ~isnan(eng_vals) & ~isnan(cleanAph_lesion.WAB_AQ(idxC));
    X = damage_clean(valid); M = eng_vals(valid); Y = cleanAph_lesion.WAB_AQ(idxC(valid));
    nVal = sum(valid);
    if nVal < 30, warning('Metric %s: N=%d too small', metric, nVal); continue; end
    
    Xz = (X - mean(X)) / std(X); Mz = (M - mean(M)) / std(M); Yz = (Y - mean(Y)) / std(Y);
    a_hat = Xz \ Mz; c_hat = Xz \ Yz; b_tmp = [Xz, Mz] \ Yz; b_hat = b_tmp(2); ab_hat = a_hat * b_hat;
    
    rng(42);
    ab_boot = NaN(nBoot, 1);
    for b = 1:nBoot
        idx_b = randsample(nVal, nVal, true);
        Xb = Xz(idx_b); Mb = Mz(idx_b); Yb = Yz(idx_b);
        if var(Xb)<1e-10 || var(Mb)<1e-10, continue; end
        a_b = Xb \ Mb; b_b = [Xb, Mb] \ Yb; b_b = b_b(2);
        ab_boot(b) = a_b * b_b;
    end
    ab_boot = ab_boot(~isnan(ab_boot));
    
    % BCa confidence intervals
    z0 = norminv(mean(ab_boot < ab_hat));
    jack_ab = NaN(nVal, 1);
    for j = 1:nVal
        idx_j = setdiff(1:nVal, j);
        Xj = Xz(idx_j); Mj = Mz(idx_j); Yj = Yz(idx_j);
        if var(Xj)<1e-10 || var(Mj)<1e-10, continue; end
        a_j = Xj \ Mj; b_j = [Xj, Mj] \ Yj; b_j = b_j(2);
        jack_ab(j) = a_j * b_j;
    end
    jack_valid = jack_ab(~isnan(jack_ab));
    if length(jack_valid) > 10
        theta_dot = mean(jack_valid);
        num = sum((theta_dot - jack_valid).^3);
        denom = 6 * (sum((theta_dot - jack_valid).^2))^1.5;
        acc = num / denom;
    else, acc = 0; end
    
    if ~isnan(z0) && ~isinf(z0)
        z_lo = norminv(alpha/2); z_hi = norminv(1-alpha/2);
        p_lo = normcdf(z0 + (z0 + z_lo) / (1 - acc*(z0 + z_lo)));
        p_hi = normcdf(z0 + (z0 + z_hi) / (1 - acc*(z0 + z_hi)));
        ci_lo = quantile(ab_boot, p_lo); ci_hi = quantile(ab_boot, p_hi);
    else
        ci_lo = quantile(ab_boot, alpha/2); ci_hi = quantile(ab_boot, 1-alpha/2);
    end
    
    sig = (ci_lo > 0) || (ci_hi < 0);
    boot_p = min(1, 2 * (mean(ab_boot <= 0) * (ab_hat >= 0) + mean(ab_boot >= 0) * (ab_hat < 0)));
    if abs(c_hat) > 0.001
        med_pct = (ab_hat / c_hat) * 100; med_pct = max(-200, min(200, med_pct));
    else, med_pct = NaN; end
    
    Level10b_Corrected = [Level10b_Corrected; table(keyICN, {metric}, nVal, a_hat, b_hat, ab_hat, ...
        ci_lo, ci_hi, boot_p, med_pct, sig, ...
        'VariableNames', {'ICN','Metric','N','Path_a','Path_b','Indirect_ab','CI_lower','CI_upper','Boot_p','Mediation_pct','Significant'})];
    fprintf('  %s: ab=%.4f, CI [%.4f, %.4f], p=%.4f, sig=%d\n', metric, ab_hat, ci_lo, ci_hi, boot_p, sig);
end
writetable(Level10b_Corrected, sprintf('%s_Level10b.csv', outPrefix));

%% ---- 7. VOLUME‑STRATIFIED ENGAGEMENT META‑ANALYSIS ----
fprintf('\n=================================================================\n');
fprintf('VOLUME-STRATIFIED ENGAGEMENT META-ANALYSIS (CLEAN)\n');
fprintf('=================================================================\n');

vol = mergedTable.LesionVolume_mm3(idxM);
q = zeros(size(vol));
q(vol <= Q1_max) = 1;
q(vol > Q1_max & vol <= Q2_max) = 2;
q(vol > Q2_max & vol <= Q3_max) = 3;
q(vol > Q3_max) = 4;
fprintf('Clean subjects per quartile:\n');
for k = 1:4, fprintf('Q%d: %d\n', k, sum(q == k)); end

Meta_Corrected = cell(length(metrics), 1);
for m = 1:length(metrics)
    metric = metrics{m};
    MetaTab = table((1:nICN)', 'VariableNames', {'ICN'});
    MetaTab.Meta_Beta = NaN(nICN,1); MetaTab.Meta_SE = NaN(nICN,1);
    MetaTab.Meta_Z = NaN(nICN,1); MetaTab.Meta_p = NaN(nICN,1);
    MetaTab.I2 = NaN(nICN,1); MetaTab.Significant_FDR = false(nICN,1);
    
    for n = 1:nICN
        engCol = sprintf('C2_ICN%02d_%s', n, metric);
        eng = cleanAph_lesion.(engCol)(idxC);
        wab = cleanAph_lesion.WAB_AQ(idxC);
        age = cleanAph_lesion.Age_At_Stroke(idxC);
        sex = cleanAph_lesion.SexNumeric(idxC);
        days = cleanAph_lesion.Days_Post_Stroke(idxC);
        
        betas = NaN(4,1); ses = NaN(4,1);
        for k = 1:4
            idx_q = (q == k) & ~isnan(eng) & ~isnan(wab) & ~isnan(age) & ~isnan(sex) & ~isnan(days);
            if sum(idx_q) < 15, continue; end
            X = [eng(idx_q), age(idx_q), sex(idx_q), days(idx_q)];
            mdl = fitlm(X, wab(idx_q));
            betas(k) = mdl.Coefficients.Estimate(2); ses(k) = mdl.Coefficients.SE(2);
        end
        
        valid_q = ~isnan(betas) & ses > 0;
        if sum(valid_q) < 2, continue; end
        
        betas_v = betas(valid_q); ses_v = ses(valid_q);
        weights = 1 ./ (ses_v.^2);
        meta_beta = sum(weights .* betas_v) / sum(weights);
        meta_se = sqrt(1 / sum(weights));
        meta_z = meta_beta / meta_se;
        meta_p = 2 * (1 - normcdf(abs(meta_z)));
        Qstat = sum(weights .* (betas_v - meta_beta).^2);
        df_Q = sum(valid_q) - 1;
        I2 = max(0, (Qstat - df_Q) / Qstat * 100);
        
        MetaTab.Meta_Beta(n) = meta_beta; MetaTab.Meta_SE(n) = meta_se;
        MetaTab.Meta_Z(n) = meta_z; MetaTab.Meta_p(n) = meta_p;
        MetaTab.I2(n) = I2;
    end
    [h, ~, p_fdr] = fdr_bh(MetaTab.Meta_p, fdrQ);
    MetaTab.Significant_FDR = h; MetaTab.p_FDR = p_fdr;
    Meta_Corrected{m} = MetaTab;
    writetable(MetaTab, sprintf('%s_VolumeStrat_Engagement_%s.csv', outPrefix, metric));
    icn17 = MetaTab(MetaTab.ICN == 17, :);
    fprintf('  %s ICN17: β=%.2f, Z=%.2f, p=%.4f, I²=%.1f%%, FDR sig=%d\n', ...
        metric, icn17.Meta_Beta, icn17.Meta_Z, icn17.Meta_p, icn17.I2, icn17.Significant_FDR);
end

%% ---- 8. CROSS‑ATLAS VALIDATION (Smith10 S06) ----
fprintf('\n=================================================================\n');
fprintf('CROSS-ATLAS VALIDATION (Smith10 S06)\n');
fprintf('=================================================================\n');
if hasSmith10
    load(smith10File, 'Smith10_Wide');
    [~, idxS, idxCA] = intersect(Smith10_Wide.PatientID, cleanAph_full.PatientID, 'stable');
    x = Smith10_Wide.S06_IRi(idxS);
    wab_ca = cleanAph_full.WAB_AQ(idxCA);
    age_ca = cleanAph_full.Age_At_Stroke(idxCA); sex_ca = cleanAph_full.SexNumeric(idxCA);
    days_ca = cleanAph_full.Days_Post_Stroke(idxCA);
    valid = ~isnan(x) & ~isnan(wab_ca) & ~isnan(age_ca) & ~isnan(sex_ca) & ~isnan(days_ca);
    n_val_ca = sum(valid);
    if n_val_ca > 3
        [rho_ca, p_ca] = partialcorr(x(valid), wab_ca(valid), [age_ca(valid), sex_ca(valid), days_ca(valid)], 'Type', 'Spearman');
        z = atanh(rho_ca); se_z = 1/sqrt(n_val_ca-3);
        CI_lo = tanh(z - 1.96*se_z); CI_hi = tanh(z + 1.96*se_z);
    else
        rho_ca = NaN; p_ca = NaN; CI_lo = NaN; CI_hi = NaN;
    end
    crossAtlas = table(rho_ca, CI_lo, CI_hi, p_ca, n_val_ca, 'VariableNames', {'rho','CI_lo','CI_hi','p_uncorrected','N'});
    writetable(crossAtlas, sprintf('%s_Smith10_S06.csv', outPrefix));
    fprintf('S06: ρ = %.3f [%.3f, %.3f], p = %.4f, N = %d\n', rho_ca, CI_lo, CI_hi, p_ca, n_val_ca);
else
    fprintf('Smith10 file not found – cross‑atlas section skipped.\n');
    crossAtlas = table();
end

%% ---- 9. IRi THRESHOLD SENSITIVITY ----
fprintf('\n=================================================================\n');
fprintf('IRi THRESHOLD SENSITIVITY (ICN17)\n');
fprintf('=================================================================\n');
IRi_cols = cell(18,1);
for n = 1:18, IRi_cols{n} = sprintf('C2_ICN%02d_IRi', n); end
IRi_matrix = cleanAph_full{:, IRi_cols};
IRi_sum = nansum(IRi_matrix, 2);
thresholds = [0.7, 0.8, 0.9];
icn17_iri = cleanAph_full.C2_ICN17_IRi;
age_s = cleanAph_full.Age_At_Stroke; sex_s = cleanAph_full.SexNumeric;
days_s = cleanAph_full.Days_Post_Stroke; wab_s = cleanAph_full.WAB_AQ;
Sensitivity = table();
for i = 1:length(thresholds)
    th = thresholds(i);
    incl = IRi_sum >= th & ~isnan(icn17_iri) & ~isnan(wab_s) & ~isnan(age_s) & ~isnan(sex_s) & ~isnan(days_s);
    n_incl = sum(incl);
    if n_incl > 3
        [rho_s, p_s] = partialcorr(icn17_iri(incl), wab_s(incl), [age_s(incl), sex_s(incl), days_s(incl)], 'Type', 'Spearman');
        z_s = atanh(rho_s); se_s = 1/sqrt(n_incl-3);
        CI_lo_s = tanh(z_s - 1.96*se_s); CI_hi_s = tanh(z_s + 1.96*se_s);
    else
        rho_s = NaN; p_s = NaN; CI_lo_s = NaN; CI_hi_s = NaN;
    end
    Sensitivity = [Sensitivity; table(th, n_incl, rho_s, p_s, CI_lo_s, CI_hi_s, ...
        'VariableNames', {'Threshold','N','rho','p_uncorrected','CI_lo','CI_hi'})];
end
disp(Sensitivity);
writetable(Sensitivity, sprintf('%s_IRi_Sensitivity.csv', outPrefix));

%% ---- 10. COMPARISON REPORT ----
fprintf('\n=================================================================\n');
fprintf('GENERATING COMPARISON REPORT\n');
fprintf('=================================================================\n');

reportFile = sprintf('%s_Comparison_Report.txt', outPrefix);
fid = fopen(reportFile, 'w');
fprintf(fid, 'QC CORRECTIVE ANALYSIS COMPARISON REPORT\n');
fprintf(fid, 'Generated: %s\n\n', datestr(now));

fprintf(fid, 'CLEAN APHASIA SAMPLES:\n');
fprintf(fid, '  Full (no lesion mask required): N = %d\n', height(cleanAph_full));
fprintf(fid, '  With lesion mask:               N = %d\n\n', height(cleanAph_lesion));

fprintf(fid, 'LEVEL 8 (Disconnection, IRi):\n');
fprintf(fid, '  Clean N = %d, Mean rho = %.3f\n\n', L8.N(1), mean(L8.Spearman_rho, 'omitnan'));

fprintf(fid, 'LEVEL 10b MEDIATION (ICN17 IRi):\n');
fprintf(fid, '  Clean: ab = %.4f, 95%% CI [%.4f, %.4f], p = %.4f (N=%d)\n', ...
    Level10b_Corrected.Indirect_ab(1), Level10b_Corrected.CI_lower(1), ...
    Level10b_Corrected.CI_upper(1), Level10b_Corrected.Boot_p(1), Level10b_Corrected.N(1));
fprintf(fid, '  (Previous: N=152, ab=-0.011, CI [-0.060, 0.033], p=0.61)\n\n');

fprintf(fid, 'VOLUME-STRATIFIED ENGAGEMENT (ICN17 IRi):\n');
newMeta = Meta_Corrected{1}(Meta_Corrected{1}.ICN == 17, :);
fprintf(fid, '  Clean: β = %.2f, Z = %.2f, p = %.4f, I² = %.1f%%\n', ...
    newMeta.Meta_Beta, newMeta.Meta_Z, newMeta.Meta_p, newMeta.I2);
fprintf(fid, '  (Previous: β = 72.3, Z = 4.68, p < 0.0001, I² = 41.1%%)\n\n');

if hasSmith10
    fprintf(fid, 'CROSS-ATLAS (Smith10 S06):\n');
    fprintf(fid, '  Clean: ρ = %.3f [%.3f, %.3f], p = %.4f, N = %d\n', ...
        crossAtlas.rho, crossAtlas.CI_lo, crossAtlas.CI_hi, crossAtlas.p_uncorrected, crossAtlas.N);
    fprintf(fid, '  (Previous: ρ = 0.408, p < 0.0001, N = 153)\n\n');
end

fprintf(fid, 'IRi THRESHOLD SENSITIVITY:\n');
old_rhos = [0.418, 0.418, 0.365];
old_Ns   = [153, 153, 92];
for i = 1:height(Sensitivity)
    fprintf(fid, '  %.1f: ρ = %.3f, p = %.4f, N = %d', ...
        Sensitivity.Threshold(i), Sensitivity.rho(i), Sensitivity.p_uncorrected(i), Sensitivity.N(i));
    fprintf(fid, '   (Prev: %.1f: ρ = %.3f, N = %d)\n', ...
        Sensitivity.Threshold(i), old_rhos(i), old_Ns(i));
end
fprintf(fid, '\nAll corrected outputs saved.\n');
fclose(fid);
fprintf('Report saved to: %s\n\n', reportFile);
fprintf('All corrections complete.\n');

%% ---- LOCAL FUNCTIONS ----
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