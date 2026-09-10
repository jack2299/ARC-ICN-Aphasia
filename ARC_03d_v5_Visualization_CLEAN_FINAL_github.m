%% ARC_03d_v5_Visualization_CLEAN_FINAL_v2.m
% =========================================================================
% ARC ICN Analysis Pipeline - Script 3d v5: CLEAN SAMPLE FIGURES (FINAL v2)
% =========================================================================
% Fixes:
%   - Figure 1 x-axis title all caps BRAINMAP20
%   - Figure 2 subscript for IRi in x-axis
%   - Figure 4 y‑axis title moved left and all caps
%   - Figure 5 y‑axis title all caps
%   - Figure 6 right panels restored, label clipping fixed via xlim
% =========================================================================
profile on
%% Setup
clear; clc; close all;

fprintf('=================================================================\n');
fprintf('ARC_03d_v5_Visualization_CLEAN_FINAL_v2.m\n');
fprintf('Publication Figure Generation - CLEAN SAMPLE (N=139 aphasia)\n');
fprintf('=================================================================\n');
fprintf('Started: %s\n\n', datestr(now));

%% Define Paths
rootDir = '/MATLAB Drive';
outputDir = fullfile(rootDir, 'ARC_03_Output', 'Figures_v5_CLEAN_FINAL');

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
    fprintf('Created output directory: %s\n', outputDir);
end

%% ========================================================================
%% LOAD AND VALIDATE CLEAN SAMPLE RESULTS
%% ========================================================================
fprintf('Loading clean sample results...\n');

cleanResultsFile = 'CleanSample_Results_Rigorous.mat';
if ~exist(cleanResultsFile, 'file')
    error('Clean sample results not found: %s', cleanResultsFile);
end

load(cleanResultsFile, 'Level1_Clean', 'Level2_Clean', 'Level3_Clean', 'Level4_Clean');
fprintf('  - Loaded Level3_Clean (%d rows)\n', height(Level3_Clean));
fprintf('  - Loaded Level4_Clean (%d rows)\n', height(Level4_Clean));
fprintf('  - Loaded Level2_Clean (%d rows)\n\n', height(Level2_Clean));

%% ========================================================================
%% LOAD ORIGINAL DATA FILES
%% ========================================================================
fprintf('Loading original data files...\n');

% Master wide file
wideFile = fullfile(rootDir, 'ARC_03b_v3_Master_Wide.mat');
if ~exist(wideFile, 'file')
    error('Master wide file not found: %s', wideFile);
end
load(wideFile, 'masterWide', 'trueICN_labels', 'nICN', 'IRi_threshold');
% Add SexNumeric if missing (needed for volume-stratification figure)
if ~ismember('SexNumeric', masterWide.Properties.VariableNames)
    masterWide.SexNumeric = double(strcmp(masterWide.Sex, 'M'));
end
fprintf('  - Loaded masterWide: %d patients\n', height(masterWide));

% ---- LOAD LESION VOLUME DATA for volume stratification (from mergedTable) ----
lesionFile = 'ARC_04_v3_AllResults.mat';
if ~exist(lesionFile, 'file')
    error('Lesion results file not found: %s', lesionFile);
end
load(lesionFile, 'mergedTable', 'resampleLog');
fprintf('  - Loaded mergedTable (%d patients) and resampleLog\n', height(mergedTable));

%% ========================================================================
%% DEFINE CLEAN SAMPLE
%% ========================================================================
fprintf('\nDefining clean sample...\n');

restOnlyList = {'M2088','M2097','M2100','M2101','M2113','M2114',...
    'M2117','M2118','M2122','M2126','M2129','M2131','M2135',...
    'M2140','M2141','M2142','M2143','M2144','M2145','M2146',...
    'M2149','M2150','M2151','M2152','M2153','M2155','M2156',...
    'M2158','M2159','M2160','M2162','M2164','M2165','M2169',...
    'M2184','M2254'};

isRestOnly = ismember(masterWide.PatientID, restOnlyList);
isAphasia = strcmp(masterWide.GroupRole, 'Aphasia');
isControl = strcmp(masterWide.GroupRole, 'Stroke Control (No Aphasia)');
hasValidC2 = ~isnan(masterWide.C2_ICN17_IRi);

cleanAphasia = masterWide(isAphasia & hasValidC2 & ~isRestOnly, :);
cleanControl = masterWide(isControl & hasValidC2 & ~isRestOnly, :);

fprintf('  - Clean aphasia sample: N = %d\n', height(cleanAphasia));
fprintf('  - Clean control sample: N = %d\n', height(cleanControl));

% Logical indices for clean sample (full masterWide length)
isAphasia_clean = isAphasia & hasValidC2 & ~isRestOnly;
isControl_clean = isControl & hasValidC2 & ~isRestOnly;

% Subtype indices for clean sample
isAnomic_clean = isAphasia_clean & strcmp(masterWide.AphasiaType, 'Anomic');
isBroca_clean = isAphasia_clean & strcmp(masterWide.AphasiaType, 'Broca');
isConduction_clean = isAphasia_clean & strcmp(masterWide.AphasiaType, 'Conduction');
isGlobal_clean = isAphasia_clean & strcmp(masterWide.AphasiaType, 'Global');

fprintf('Clean subtypes: Anomic=%d, Broca=%d, Conduction=%d, Global=%d\n', ...
    sum(isAnomic_clean), sum(isBroca_clean), sum(isConduction_clean), sum(isGlobal_clean));

% ICN matrix for clean sample
X_C2_clean = NaN(height(cleanAphasia), nICN);
for n = 1:nICN
    colC2 = sprintf('C2_ICN%02d_IRi', n);
    if ismember(colC2, masterWide.Properties.VariableNames)
        X_C2_clean(:, n) = masterWide{isAphasia_clean, colC2};
    end
end

%% ========================================================================
%% DEFINE WONG COLORBLIND-SAFE PALETTE
%% ========================================================================
colors.blue = [0, 114, 178] / 255;
colors.vermillion = [213, 94, 0] / 255;
colors.green = [0, 158, 115] / 255;
colors.yellow = [240, 228, 66] / 255;
colors.skyblue = [86, 180, 233] / 255;
colors.orange = [230, 159, 0] / 255;
colors.purple = [204, 121, 167] / 255;
colors.nonsig = [0.6, 0.6, 0.6];

subtypes = {'Anomic', 'Broca', 'Conduction', 'Global'};
subtypeColors = [colors.blue; colors.vermillion; colors.green; colors.purple];

% A priori language networks
languageICNs = [4, 16, 18];

fprintf('Setting up Wong colorblind-safe palette...\n\n');

%% Helper Functions
function saveFigAll(fig, filepath)
    savefig(fig, [filepath '.fig']);
    print(fig, [filepath '.png'], '-dpng', '-r600');
    print(fig, [filepath '.tif'], '-dtiff', '-r600');
    print(fig, [filepath '.svg'], '-dsvg');
    fprintf('    Saved: %s (.fig, .png, .tif, .svg)\n', filepath);
end

%% ========================================================================
%% FIGURE S1 (Script Figure 1): APHASIA VS STROKE CONTROL
%% ========================================================================
fprintf('Creating Figure S1: Aphasia vs Stroke Control (Clean Sample)...\n');

figS1 = figure('Position', [100 100 1400 600], 'Color', 'w');

N_aphasia = height(cleanAphasia);
N_control = height(cleanControl);

hold on;
jitterWidth = 0.15;

for n = 1:nICN
    aphasiaVals = X_C2_clean(:, n);
    aphasiaVals = aphasiaVals(~isnan(aphasiaVals));
    nAph = length(aphasiaVals);
    xAphasia = n - 0.2 + jitterWidth * (rand(nAph, 1) - 0.5);
    scatter(xAphasia, aphasiaVals, 15, colors.vermillion, 'filled', ...
        'MarkerFaceAlpha', 0.3, 'MarkerEdgeColor', 'none');
    
    colC2 = sprintf('C2_ICN%02d_IRi', n);
    controlVals = masterWide{isControl_clean, colC2};
    controlVals = controlVals(~isnan(controlVals));
    nCtrl = length(controlVals);
    xControl = n + 0.2 + jitterWidth * (rand(nCtrl, 1) - 0.5);
    scatter(xControl, controlVals, 15, colors.nonsig, 'filled', ...
        'MarkerFaceAlpha', 0.3, 'MarkerEdgeColor', 'none');
    
    meanAph = mean(aphasiaVals, 'omitnan');
    meanCtrl = mean(controlVals, 'omitnan');
    plot(n - 0.2, meanAph, 'o', 'MarkerSize', 8, 'MarkerFaceColor', colors.vermillion, ...
        'MarkerEdgeColor', 'w', 'LineWidth', 1);
    plot(n + 0.2, meanCtrl, 'o', 'MarkerSize', 8, 'MarkerFaceColor', colors.nonsig, ...
        'MarkerEdgeColor', 'w', 'LineWidth', 1);
end

xlabel('BRAINMAP20 Network', 'FontSize', 10, 'FontWeight', 'bold');   % ALL CAPS
ylabel('Network Engagement (IR_i)', 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'tex'); % subscript 'i'

% Replace x-tick labels with BMxx
set(gca, 'XTick', 1:nICN, 'XTickLabel', ...
    arrayfun(@(x) sprintf('BM%02d', x), 1:nICN, 'UniformOutput', false), ...
    'FontSize', 10);

h1 = scatter(NaN, NaN, 50, colors.vermillion, 'filled');
h2 = scatter(NaN, NaN, 50, colors.nonsig, 'filled');
legend([h1 h2], {sprintf('Aphasia (N=%d)', N_aphasia), sprintf('Stroke Control (N=%d)', N_control)}, ...
    'Location', 'northeast', 'FontSize', 10, 'Box', 'off');

xlim([0.5 nICN+0.5]);
grid on; box off;

saveFigAll(figS1, fullfile(outputDir, 'FigS1_AphasiaVsControl_CLEAN'));


%% ========================================================================
%% FIGURE 4 (Script Figure 2): TOP SEVERITY PREDICTOR (BM17) – TWO-PANEL
%% ========================================================================
fprintf('Creating Figure 2: Top Severity Predictor (BM17) – Two-Panel...\n');

fig4 = figure('Position', [100 100 1400 600], 'Color', 'w');

% ---- Extract clean data ----
topData_clean = masterWide{isAphasia_clean, 'C2_ICN17_IRi'};
wab_r_clean = masterWide.WAB_AQ(isAphasia_clean);
age_clean = masterWide.Age_At_Stroke(isAphasia_clean);
sex_clean = masterWide.SexNumeric(isAphasia_clean);
days_clean = masterWide.Days_Post_Stroke(isAphasia_clean);

validIdx = ~isnan(topData_clean) & ~isnan(wab_r_clean) & ~isnan(age_clean) & ~isnan(sex_clean) & ~isnan(days_clean);
topData_valid = topData_clean(validIdx);
wab_valid = wab_r_clean(validIdx);
age_valid = age_clean(validIdx);
sex_valid = sex_clean(validIdx);
days_valid = days_clean(validIdx);
subtype_valid = masterWide.AphasiaType(isAphasia_clean);
subtype_valid = subtype_valid(validIdx);

% ---- Define five subtypes for severity figure ----
subtypeList5 = {'Anomic','Broca','Conduction','Global','Wernicke'};
subtypeColors5 = [colors.blue; colors.vermillion; colors.green; colors.purple; colors.yellow];

% ========= PANEL A: Raw scatter, no line =========
subplot(1,2,1);
hold on;
for s = 1:length(subtypeList5)
    subMask = strcmp(subtype_valid, subtypeList5{s});
    scatter(topData_valid(subMask), wab_valid(subMask), 50, subtypeColors5(s,:), ...
        'filled', 'MarkerFaceAlpha', 0.75);
end
xlabel('BM17 IR_i (Relative Spatial Involvement)', 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'tex');
ylabel('WAB-R Score (Aphasia Severity)', 'FontSize', 10, 'FontWeight', 'bold');
legend(subtypeList5, 'Location', 'southeast', 'FontSize', 9);
grid on; box off;
set(gca, 'FontSize', 10);
text(0.02, 0.98, 'A', 'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% ========= PANEL B: Partial residual plot =========
subplot(1,2,2);
hold on;

% Rank transform both variables
rank_IRi = tiedrank(topData_valid);
rank_WAB = tiedrank(wab_valid);

% Regress ranks on covariates and get residuals
X_cov = [age_valid, sex_valid, days_valid];
mdl_IRi = fitlm(X_cov, rank_IRi);
mdl_WAB = fitlm(X_cov, rank_WAB);
resid_IRi = mdl_IRi.Residuals.Raw;
resid_WAB = mdl_WAB.Residuals.Raw;

% Scatter residuals by subtype
for s = 1:length(subtypeList5)
    subMask = strcmp(subtype_valid, subtypeList5{s});
    scatter(resid_IRi(subMask), resid_WAB(subMask), 50, subtypeColors5(s,:), ...
        'filled', 'MarkerFaceAlpha', 0.75);
end

% Add least-squares line on residuals
coeffs = polyfit(resid_IRi, resid_WAB, 1);
xFit = linspace(min(resid_IRi), max(resid_IRi), 100);
yFit = polyval(coeffs, xFit);
plot(xFit, yFit, 'k-', 'LineWidth', 3);

xlabel('Rank-based partial residuals: BM17 IR_i', 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'tex');
ylabel('Rank-based partial residuals: WAB-R', 'FontSize', 10, 'FontWeight', 'bold');

grid on; box off;
set(gca, 'FontSize', 10);
text(0.02, 0.98, 'B', 'Units', 'normalized', 'FontSize', 16, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

saveFigAll(fig4, fullfile(outputDir, 'Fig2_TopPredictor_CLEAN'));

%% ========================================================================
%% FIGURE S2 (Script Figure 3): SUBTYPE PROFILE DIFFERENCES
%% ========================================================================
fprintf('Creating Figure S2: Subtype Differences (Clean Sample)...\n');

sigIdx = Level2_Clean.Significant_FDR & strcmp(Level2_Clean.Metric, 'IRi');
sigLevel2ICNs_clean = Level2_Clean.ICN(sigIdx);

subtypes4 = {'Anomic','Broca','Conduction','Global'};
minNForErrorBar = 5;
errorLineWidth = 2.5;
errorCapSize = 12;

if ~isempty(sigLevel2ICNs_clean)
    nSigICN = length(sigLevel2ICNs_clean);
    if nSigICN == 1
        figS2 = figure('Position', [100 100 600 500], 'Color', 'w');
    else
        figS2 = figure('Position', [100 100 1200 500], 'Color', 'w');
    end

    for i = 1:nSigICN
        if nSigICN == 1
            subplot(1, 1, 1);
        else
            subplot(1, nSigICN, i);
        end

        icnNum = sigLevel2ICNs_clean(i);
        colName = sprintf('C2_ICN%02d_IRi', icnNum);
        if ~ismember(colName, masterWide.Properties.VariableNames)
            warning('Column %s not found. Skipping ICN%d.', colName, icnNum);
            continue;
        end

        icnData = masterWide{isAphasia_clean, colName};

        isAnomic_rel = strcmp(masterWide.AphasiaType(isAphasia_clean), 'Anomic');
        isBroca_rel   = strcmp(masterWide.AphasiaType(isAphasia_clean), 'Broca');
        isConduction_rel = strcmp(masterWide.AphasiaType(isAphasia_clean), 'Conduction');
        isGlobal_rel   = strcmp(masterWide.AphasiaType(isAphasia_clean), 'Global');
        groupIdx_list = {isAnomic_rel, isBroca_rel, isConduction_rel, isGlobal_rel};

        hold on;
        jitterWidth = 0.25;
        groupMeans = NaN(1,4);
        groupSEMs = NaN(1,4);
        groupNs = NaN(1,4);

        % Scatter individual participants
        for g = 1:4
            vals = icnData(groupIdx_list{g});
            vals = vals(~isnan(vals));
            nVals = length(vals);
            groupNs(g) = nVals;

            if nVals > 0
                xJitter = g + jitterWidth * (rand(nVals, 1) - 0.5);
                scatter(xJitter, vals, 20, subtypeColors(g, :), 'filled', ...
                    'MarkerFaceAlpha', 0.3, 'MarkerEdgeColor', 'none');
                groupMeans(g) = mean(vals);
                groupSEMs(g) = std(vals) / sqrt(nVals);
            end
        end

        % Error bars only for groups with n >= threshold
        for g = 1:4
            if ~isnan(groupMeans(g)) && ~isnan(groupSEMs(g)) && groupNs(g) >= minNForErrorBar
                errorbar(g, groupMeans(g), groupSEMs(g), 'k', ...
                    'LineStyle', 'none', ...
                    'LineWidth', errorLineWidth, ...
                    'CapSize', errorCapSize);
            end
        end

        % Mean markers for larger groups; line for small Global group
        for g = 1:4
            if ~isnan(groupMeans(g))
                if groupNs(g) >= minNForErrorBar
                    plot(g, groupMeans(g), 'o', 'MarkerSize', 10, ...
                        'MarkerFaceColor', subtypeColors(g, :), ...
                        'MarkerEdgeColor', 'w', 'LineWidth', 2);
                else
                    plot([g-0.15, g+0.15], [groupMeans(g), groupMeans(g)], '-', ...
                        'Color', subtypeColors(g, :), 'LineWidth', 2.5);
                end
            end
        end

        % X-axis labels: subtype names only
        set(gca, 'XTick', 1:4, 'XTickLabel', subtypes4, ...
            'FontSize', 10, 'TickLabelInterpreter', 'none');
        ylabel('Network Engagement (IR_i)', 'FontSize', 9, 'FontWeight', 'bold', 'Interpreter', 'tex');
        xlabel('Aphasia Subtype', 'FontSize', 10, 'FontWeight', 'bold');
        xlim([0.5 4.5]);
        box off; grid on;
    end
else
    figS2 = figure('Position', [100 100 800 500], 'Color', 'w');
    annotation('textbox', [0.3 0.4 0.4 0.2], 'String', ...
        'No significant subtype effects in clean sample', ...
        'HorizontalAlignment', 'center', 'FontSize', 14, 'EdgeColor', 'none');
end

saveFigAll(figS2, fullfile(outputDir, 'FigS2_SubtypeDifferences_CLEAN'));

%% ========================================================================
%% FIGURE 6 (Script Figure 4): FOREST PLOT
%% ========================================================================
fprintf('Creating Figure 3: Forest Plot (Clean Sample)...\n');

fig6 = figure('Position', [100 100 1200 800], 'Color', 'w');

L3_primary = Level3_Clean(strcmp(Level3_Clean.Metric, 'IRi'), :);
if height(L3_primary) > 0
    [effects_sorted, sortIdx] = sort(L3_primary.Partial_rho, 'descend');
    icns_sorted = L3_primary.ICN(sortIdx); sigFlag_sorted = L3_primary.Significant_FDR(sortIdx);
    if ismember('Rho_CI_lo', L3_primary.Properties.VariableNames)
        ci_low = L3_primary.Rho_CI_lo(sortIdx); ci_high = L3_primary.Rho_CI_hi(sortIdx);
    else
        n_subj = L3_primary.N(1); n_cov = 3;
        ci_low = NaN(length(effects_sorted),1); ci_high = NaN(length(effects_sorted),1);
        for i = 1:length(effects_sorted)
            r = effects_sorted(i);
            if abs(r) < 0.999
                z = 0.5*log((1+r)/(1-r)); se = 1/sqrt(n_subj-3-n_cov);
                ci_low(i)=tanh(z-1.96*se); ci_high(i)=tanh(z+1.96*se);
            end
        end
    end
end

nEffects = length(effects_sorted);
effects_plot = flip(effects_sorted); icns_plot = flip(icns_sorted); sigFlag_plot = flip(sigFlag_sorted);
ci_low_plot = flip(ci_low); ci_high_plot = flip(ci_high); yPos = 1:nEffects;

hold on; plot([0 0],[0 nEffects+1],'k--','LineWidth',1);

% Grey band (zero-effect region) – keep handle for legend
hBand = fill([-0.1 0.1 0.1 -0.1],[0 0 nEffects+1 nEffects+1],[0.95 0.95 0.95],'EdgeColor','none');

for i=1:nEffects
    if ismember(icns_plot(i),languageICNs)
        plot(effects_plot(i),yPos(i),'s','MarkerSize',22,'MarkerEdgeColor',colors.skyblue,'LineWidth',2,'MarkerFaceColor','none');
    end
end
for i=1:nEffects
    if sigFlag_plot(i), markerColor=colors.vermillion; else markerColor=colors.nonsig; end
    if ~isnan(ci_low_plot(i)) && ~isnan(ci_high_plot(i))
        plot([ci_low_plot(i) ci_high_plot(i)],[yPos(i) yPos(i)],'-','Color',markerColor,'LineWidth',2);
    end
end
for i=1:nEffects
    if sigFlag_plot(i), markerColor=colors.vermillion; markerSize=14; else markerColor=colors.nonsig; markerSize=10; end
    plot(effects_plot(i),yPos(i),'o','MarkerSize',markerSize,'MarkerFaceColor',markerColor,'MarkerEdgeColor',markerColor,'LineWidth',1);
end

% ---------- normal axes, no artificial padding ----------
ax = gca;
xlim([-0.2 0.55]);
ylim([0.5 nEffects+0.5]);

% Remove automatic y tick labels so we can place them manually and avoid clipping
set(ax, 'YTick', yPos, 'YTickLabel', []);

% Manually place y labels, locked to each row, offset left
xLimits = get(ax, 'XLim');
xLabelPos = xLimits(1) - 0.005 * diff(xLimits);

yLabels = cell(nEffects,1);
for i = 1:nEffects
    yLabels{i} = sprintf('BM%02d', icns_plot(i));
end

for i = 1:nEffects
    text(xLabelPos, yPos(i), yLabels{i}, ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 14, ...
        'FontWeight', 'normal', ...
        'Color', 'k', ...
        'Clipping', 'off');
end

% Move y-axis title further left
ylabel('BRAINMAP20 Network', 'FontSize', 14, 'FontWeight', 'bold');
ylh = get(ax, 'YLabel');
ylh.Position(1) = xLabelPos - 0.10 * diff(xLimits);

grid on; box off;          % do this BEFORE touching the x ruler

% x-axis ticks: use built-in labels with format forced
xticks(ax, [-0.2 -0.1 0 0.1 0.2 0.3 0.4 0.5]);
ax.XAxis.Exponent = 0;                 % no shared exponent
ax.XAxis.TickLabelFormat = '%.1f';     % forces 0.2, never 2

xlabel(ax, 'Partial Correlation (\rho) with WAB-R', 'FontSize', 14, 'FontWeight', 'bold');

% legend
h1 = plot(NaN,NaN,'o','MarkerSize',14,'MarkerFaceColor',colors.vermillion,'MarkerEdgeColor',colors.vermillion);
h2 = plot(NaN,NaN,'o','MarkerSize',10,'MarkerFaceColor',colors.nonsig,'MarkerEdgeColor',colors.nonsig);
h3 = plot(NaN,NaN,'s','MarkerSize',16,'MarkerEdgeColor',colors.skyblue,'LineWidth',2,'MarkerFaceColor','none');
legend([h1 h2 h3 hBand], {'Significant (FDR)','Non-significant','Language ICN','Zero-effect band'}, ...
    'Location','southeast','FontSize',14,'Box','on');

saveFigAll(fig6, fullfile(outputDir, 'Fig3_ForestPlot_CLEAN'));


%% ========================================================================
%% FIGURE 7 (Script Figure 5): C2 INCREMENTAL PREDICTION (Lollipop)
%% ========================================================================
fprintf('Creating Figure 4: C2 Incremental Prediction (Clean Sample)...\n');

fig7 = figure('Position',[100 100 1000 600],'Color','w');
L4_primary=Level4_Clean(strcmp(Level4_Clean.Metric,'IRi'),:);
if height(L4_primary)>0
    [~,sortIdx]=sort(L4_primary.R2_improvement,'descend'); L4_sorted=L4_primary(sortIdx,:);
    hold on; plot([0 0],[0.5 nICN+0.5],'k--','LineWidth',1);
    for n=1:height(L4_sorted)
        yPos=nICN-n+1;
        if L4_sorted.Significant_FDR(n), markerColor=colors.orange; markerSize=12;
        elseif ismember(L4_sorted.ICN(n),languageICNs), markerColor=colors.skyblue; markerSize=10;
        else markerColor=colors.nonsig; markerSize=8; end
        plot([0 L4_sorted.R2_improvement(n)],[yPos yPos],'-','Color',markerColor,'LineWidth',2);
        plot(L4_sorted.R2_improvement(n),yPos,'o','MarkerSize',markerSize,'MarkerFaceColor',markerColor,'MarkerEdgeColor',markerColor);
    end
    yLabels=cell(nICN,1); for n=1:nICN, yLabels{nICN-n+1}=sprintf('BM%02d',L4_sorted.ICN(n)); end
    set(gca,'YTick',1:nICN,'YTickLabel',yLabels,'FontSize',10);
    xlabel('R^2 Improvement','FontSize',10,'FontWeight','bold');
    ylabel('BRAINMAP20 Network','FontSize',10,'FontWeight','bold'); % all caps, no (sorted)
    ylim([0.5 nICN+0.5]); xlim([-0.02 max(L4_sorted.R2_improvement)*1.1]); grid on; box off;
    % Shift y-label left to avoid blocking labels
    ax7 = gca; ylh7 = get(ax7,'YLabel'); xl = get(ax7,'XLim');
    ylh7.Position(1) = xl(1) - 0.08*diff(xl);   % manual offset
    h1=plot(NaN,NaN,'o','MarkerSize',12,'MarkerFaceColor',colors.orange,'MarkerEdgeColor',colors.orange);
    h2=plot(NaN,NaN,'o','MarkerSize',10,'MarkerFaceColor',colors.skyblue,'MarkerEdgeColor',colors.skyblue);
    h3=plot(NaN,NaN,'o','MarkerSize',8,'MarkerFaceColor',colors.nonsig,'MarkerEdgeColor',colors.nonsig);
    legend([h1 h2 h3],{'Significant (FDR)','Language ICN','Non-significant'},'Location','southeast','FontSize',10,'Box','on');
end
saveFigAll(fig7,fullfile(outputDir,'Fig4_ContrastComparison_CLEAN'));

%% ========================================================================
%% FIGURE: VOLUME STRATIFICATION (Script Figure 6)
%% ========================================================================
%% ========================================================================
%% FIGURE 5: VOLUME STRATIFICATION (CORRECTED)
%% ========================================================================
fprintf('Creating Figure 5: Volume Stratification (CORRECTED)...\n');

lesionIDs = unique(resampleLog.PatientID(~cellfun('isempty',resampleLog.PatientID)));
cleanAphIDs = cleanAphasia.PatientID;
hasLesion = ismember(cleanAphIDs, lesionIDs);
cleanAphLesion = cleanAphasia(hasLesion, :);
fprintf('  Clean aphasia with lesion mask: N = %d\n', height(cleanAphLesion));

[~, idxM, idxC] = intersect(mergedTable.PatientID, cleanAphLesion.PatientID, 'stable');
vol = mergedTable.LesionVolume_mm3(idxM);

Q1_max = 91556;
Q2_max = 173496;
Q3_max = 280746;
q = zeros(size(vol));
q(vol <= Q1_max) = 1;
q(vol > Q1_max & vol <= Q2_max) = 2;
q(vol > Q2_max & vol <= Q3_max) = 3;
q(vol > Q3_max) = 4;

eng = cleanAphLesion.C2_ICN17_IRi(idxC);
wab = cleanAphLesion.WAB_AQ(idxC);
age = cleanAphLesion.Age_At_Stroke(idxC);
sex = cleanAphLesion.SexNumeric(idxC);
days = cleanAphLesion.Days_Post_Stroke(idxC);

betas = NaN(4,1);
ses = NaN(4,1);
Ns = NaN(4,1);
for k = 1:4
    idx_q = (q == k) & ~isnan(eng) & ~isnan(wab) & ~isnan(age) & ~isnan(sex) & ~isnan(days);
    Ns(k) = sum(idx_q);
    if sum(idx_q) >= 15
        X = [eng(idx_q), age(idx_q), sex(idx_q), days(idx_q)];
        mdl = fitlm(X, wab(idx_q));
        betas(k) = mdl.Coefficients.Estimate(2);
        ses(k) = mdl.Coefficients.SE(2);
    end
end

% Quartile confidence intervals
CI_lo = betas - 1.96 * ses;
CI_hi = betas + 1.96 * ses;

% Fixed-effects inverse-variance meta-analysis across quartiles
valid = ~isnan(betas) & ~isnan(ses);
w = 1 ./ ses(valid).^2;
meta_Beta = sum(w .* betas(valid)) / sum(w);
meta_SE = sqrt(1 / sum(w));
meta_Z = meta_Beta / meta_SE;
meta_p_unc = 2 * (1 - normcdf(abs(meta_Z)));
Q = sum(w .* (betas(valid) - meta_Beta).^2);
df = sum(valid) - 1;
meta_I2 = max(0, (Q - df) / Q) * 100;
meta_CI_lo = meta_Beta - 1.96 * meta_SE;
meta_CI_hi = meta_Beta + 1.96 * meta_SE;

figVol = figure('Position', [100 100 1100 600], 'Color', 'w');
yPos_Q = [5, 4, 3, 2];
yPos_Meta = 1;

% Main forest panel
ax1 = axes('Position', [0.10 0.12 0.48 0.80]);
hold on;
plot([0 0], [0.5 5.5], 'k--', 'LineWidth', 1);

for i = 1:4
    if ~isnan(betas(i))
        plot([CI_lo(i) CI_hi(i)], [yPos_Q(i) yPos_Q(i)], '-', ...
            'Color', colors.blue, 'LineWidth', 2);
        plot(betas(i), yPos_Q(i), 'o', 'MarkerSize', 12, ...
            'MarkerFaceColor', colors.blue, 'MarkerEdgeColor', colors.blue);
    end
end

plot([meta_CI_lo meta_CI_hi], [yPos_Meta yPos_Meta], '-', ...
    'Color', colors.vermillion, 'LineWidth', 2.5);
plot(meta_Beta, yPos_Meta, 'd', 'MarkerSize', 14, ...
    'MarkerFaceColor', colors.vermillion, 'MarkerEdgeColor', colors.vermillion);
plot([meta_Beta meta_Beta], [0.5 5.5], ':', ...
    'Color', colors.vermillion, 'LineWidth', 1.5);

yLabels = {'Meta-\newlineanalytic', 'Q4', 'Q3', 'Q2', 'Q1'};
set(gca, 'YTick', 1:5, 'YTickLabel', yLabels, 'FontSize', 10, 'TickLabelInterpreter', 'tex');
xlabel('\beta Coefficient (Engagement \rightarrow Severity)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Lesion Volume Stratum', 'FontSize', 10, 'FontWeight', 'bold');
xlim([-50 260]);
ylim([0.5 5.5]);
grid on;
box off;

% Right panels for N and beta values
y_norm_positions = [0.85, 0.68, 0.50, 0.32, 0.10];

ax2 = axes('Position', [0.62 0.12 0.08 0.80]);
axis off;
text(0.5, 0.97, 'N', 'FontSize', 10, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
for i = 1:4
    text(0.5, y_norm_positions(i), sprintf('%d', Ns(i)), ...
        'FontSize', 10, 'HorizontalAlignment', 'center');
end
text(0.5, y_norm_positions(5), sprintf('%d', 138), 'FontSize', 10, ...
    'HorizontalAlignment', 'center', 'Color', colors.vermillion, 'FontWeight', 'bold');
rectangle('Position', [0, 0, 1, 1], 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 1);

ax3 = axes('Position', [0.72 0.12 0.26 0.80]);
axis off;
text(0.5, 0.97, '\beta [95% CI]', 'FontSize', 10, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
for i = 1:4
    if ~isnan(betas(i))
        txt = sprintf('%.1f [%.1f, %.1f]', betas(i), CI_lo(i), CI_hi(i));
        text(0.5, y_norm_positions(i), txt, 'FontSize', 9, 'HorizontalAlignment', 'center');
    end
end
text(0.5, y_norm_positions(5) + 0.05, sprintf('%.1f [%.1f, %.1f]', meta_Beta, meta_CI_lo, meta_CI_hi), ...
    'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', colors.vermillion);

% Update meta-analytic annotation: p_unc in panel, p_FDR noted in caption
text(0.5, y_norm_positions(5) - 0.03, ...
    sprintf('(Z=%.2f, p_{unc}=%.4f, I²=%.1f%%)', meta_Z, meta_p_unc, meta_I2), ...
    'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
    'Color', colors.vermillion);
rectangle('Position', [0, 0, 1, 1], 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 1);


saveFigAll(figVol, fullfile(outputDir, 'Fig5_VolumeStratification_CLEAN'));

%% ========================================================================
%% SUMMARY
%% ========================================================================
fprintf('\n=================================================================\n');
fprintf('VISUALIZATION COMPLETE - CLEAN SAMPLE (FINAL v2)\n');
fprintf('=================================================================\n');
fprintf('Output directory: %s\n\n', outputDir);
figFiles=dir(fullfile(outputDir,'*.png')); fprintf('Generated figures (%d total):\n',length(figFiles));
for f=1:length(figFiles), fprintf('  - %s\n',figFiles(f).name); end
fprintf('\nFIXES APPLIED:\n');
fprintf('  - Figure 1 x-axis all caps BRAINMAP20\n');
fprintf('  - Figure 2 subscript IR_i\n');
fprintf('  - Figure 4 ylabel moved left & all caps BRAINMAP20\n');
fprintf('  - Figure 5 ylabel all caps BRAINMAP20\n');
fprintf('  - Figure 6 right panels restored, xlim extended to 260\n');
close all;
fprintf('\n=================================================================\n');
fprintf('Completed: %s\n', datestr(now));
fprintf('=================================================================\n');