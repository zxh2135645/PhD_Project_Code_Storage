clear all;
close all;

addpath('../function/');
A_clin_conc = readtable('../../For_Keyur/2023/Regression_ROC_Plots_CLIN_CONC.xlsx');
A_clin_velo = readtable('../../For_Keyur/2023/Regression_ROC_Plots_CLIN_VELO.xlsx');
A_preclin_conc = readtable('../../For_Keyur/2023/Regression_ROC_Plots_PRECLIN_CONC.xlsx');
A_preclin_velo = readtable('../../For_Keyur/2023/Regression_ROC_Plots_PRECLIN_VELO.xlsx');

%% BlandAltman (Should skip)
addpath('../function/BlandAltman/');
territories = {'Invivo'};
nterritories = length(territories);

% Patient states during measurement
states = {'Volume (%)'};
nstates = length(states);

BIOM_1 = A_clin_conc.BIOM;
IMH_1 = A_clin_conc.IMH;

data1 = cat(3,BIOM_1);
data2 = cat(3, IMH_1);

% BA plot paramters
tit = 'Invivo-Exvivo Comparison'; % figure title
gnames = {territories, states}; % names of groups in data {dimension 1 and 2}
label = {'Invivo','Ground Truth','%'}; % Names of data sets
corrinfo = {'n','SSE','r2','eq'}; % stats to display of correlation scatter plot
BAinfo = {'RPC(%)','ks'}; % stats to display on Bland-ALtman plot
limits = 'auto'; % how to set the axes limits
if 1 % colors for the data sets may be set as:
	colors = 'br';      % character codes
else
	colors = [0 0 1;... % or RGB triplets
		      1 0 0];
end

% Generate figure with symbols
[cr, fig, statsStruct] = BlandAltman(data1, data2, label, tit, gnames, 'corrInfo', corrinfo, 'baInfo', BAinfo, 'axesLimits', limits, 'colors', colors, 'showFitCI', 'on');
%% Linear regression
BIOM_1 = A_clin_conc.BIOM;
IMH_1 = A_clin_conc.IMH;
mdl = fitlm(BIOM_1, IMH_1);
ci = coefCI(mdl);
% figure();
% plot(mdl)
%% LG-1
figure('Position', [0 100 400 400]);
itercept = mdl.Coefficients.Estimate(1);
slope = mdl.Coefficients.Estimate(2);
x0 = [0 1200];
y0 = slope * x0 + itercept;
yL = ci(2,1) * x0 + ci(1,1);
yU = ci(2,2) * x0 + ci(1,2);
hData = plot(BIOM_1, IMH_1,  'o');
hModel = line(x0, y0);
hCI(1) = line(x0, yL);
hCI(2) = line(x0, yU);
patch([x0 fliplr(x0)], [yL fliplr(yU)], 'r', 'FaceColor', [1 0 0], 'FaceAlpha', 0.1)
set(hData, 'LineStyle', 'none', 'Marker', '.')
set(hModel, 'LineStyle', '-', 'Color', 'r')
set(hCI(1), 'LineStyle', '-.', 'Color', [0 .5 0])
set(hCI(2), 'LineStyle', '-.', 'Color', [0 .5 0])

set(hData, 'Marker', 'o', 'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', [.75 .75 1])
set(hModel, 'LineWidth', 3)
set(hCI(1), 'LineWidth', 2)
set(hCI(2), 'LineWidth', 2)

hTitle = title('CLIN CONC');
hXLabel = xlabel('BIOM');
hYLabel = ylabel('IMH Vol (%LV)');

% Add text
hText = text(800, 10, ...
     sprintf('{\\itY = %0.1gX + %0.1g}', slope, itercept));
hText1 = text(800, 9, ...
     sprintf('{\\itR^2 = %0.4g}', mdl.Rsquared.Ordinary));
hText2 = text(800, 8, ...
     sprintf('{\\itp-value < 10^{-3}}', mdl.Rsquared.Ordinary));
% Add legend
hLegend = legend([hData, hModel, hCI(1)], ...
    'Observed', 'Regression line', '95% CI', ...
    'Location', 'NorthWest');

% Adjust font
set(gca, 'FontName', 'Helvetica')
set([hTitle, hXLabel, hYLabel, hText, hText1, hText2], 'FontName', 'AvantGarde')
% set([hTitle, hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hLegend, gca], 'FontSize', 12)
set([hXLabel, hYLabel, hText, hText1, hText2], 'FontSize', 12)
% set([hXLabel, hYLabel], 'FontSize', 10)
set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')

% Adjust axes properties
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', 0:5:20, ...
    'LineWidth', 1.5);
ylim([0 20]);
%% LG-2
BIOM_2 = A_preclin_conc.BIOMPRE;
IMH_2 = A_preclin_conc.IMH;
mdl = fitlm(BIOM_2, IMH_2);
ci = coefCI(mdl);

figure('Position', [0 100 400 400]);
itercept = mdl.Coefficients.Estimate(1);
slope = mdl.Coefficients.Estimate(2);
x0 = [0 1200];
y0 = slope * x0 + itercept;
yL = ci(2,1) * x0 + ci(1,1);
yU = ci(2,2) * x0 + ci(1,2);
hData = plot(BIOM_2, IMH_2,  'o');
hModel = line(x0, y0);
hCI(1) = line(x0, yL);
hCI(2) = line(x0, yU);
patch([x0 fliplr(x0)], [yL fliplr(yU)], 'r', 'FaceColor', [1 0 0], 'FaceAlpha', 0.1)
set(hData, 'LineStyle', 'none', 'Marker', '.')
set(hModel, 'LineStyle', '-', 'Color', 'r')
set(hCI(1), 'LineStyle', '-.', 'Color', [0 .5 0])
set(hCI(2), 'LineStyle', '-.', 'Color', [0 .5 0])

set(hData, 'Marker', 'o', 'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', [.75 .75 1])
set(hModel, 'LineWidth', 3)
set(hCI(1), 'LineWidth', 2)
set(hCI(2), 'LineWidth', 2)

hTitle = title('PRECLIN CONC');
hXLabel = xlabel('BIOM');
hYLabel = ylabel('IMH Vol (%LV)');

% Add text
hText = text(320, 10, ...
     sprintf('{\\itY = %0.1gX + %0.1g}', slope, itercept));
hText1 = text(320, 9, ...
     sprintf('{\\itR^2 = %0.4g}', mdl.Rsquared.Ordinary));
hText2 = text(320, 8, ...
     sprintf('{\\itp-value < 10^{-3}}', mdl.Rsquared.Ordinary));
% Add legend
hLegend = legend([hData, hModel, hCI(1)], ...
    'Observed', 'Regression line', '95% CI', ...
    'Location', 'NorthWest');

% Adjust font
set(gca, 'FontName', 'Helvetica')
set([hTitle, hXLabel, hYLabel, hText, hText1, hText2], 'FontName', 'AvantGarde')
% set([hTitle, hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hLegend, gca], 'FontSize', 12)
set([hXLabel, hYLabel, hText, hText1, hText2], 'FontSize', 12)
% set([hXLabel, hYLabel], 'FontSize', 10)
set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')

% Adjust axes properties
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', 0:5:20, ...
    'LineWidth', 1.5);
ylim([0 25]); xlim([0 500])

%% LG-1 vs VELO
BIOM_1_v = A_clin_velo.BIOM;
IMH_1_v = A_clin_velo.IMH;
mdl = fitlm(BIOM_1_v, IMH_1_v);
ci = coefCI(mdl);

figure('Position', [0 100 400 400]);
itercept = mdl.Coefficients.Estimate(1);
slope = mdl.Coefficients.Estimate(2);
x0 = [0 1200];
y0 = slope * x0 + itercept;
yL = ci(2,1) * x0 + ci(1,1);
yU = ci(2,2) * x0 + ci(1,2);
hData = plot(BIOM_1_v, IMH_1_v,  'o');
hModel = line(x0, y0);
hCI(1) = line(x0, yL);
hCI(2) = line(x0, yU);
patch([x0 fliplr(x0)], [yL fliplr(yU)], 'r', 'FaceColor', [1 0 0], 'FaceAlpha', 0.1)
set(hData, 'LineStyle', 'none', 'Marker', '.')
set(hModel, 'LineStyle', '-', 'Color', 'r')
set(hCI(1), 'LineStyle', '-.', 'Color', [0 .5 0])
set(hCI(2), 'LineStyle', '-.', 'Color', [0 .5 0])

set(hData, 'Marker', 'o', 'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', [.75 .75 1])
set(hModel, 'LineWidth', 3)
set(hCI(1), 'LineWidth', 2)
set(hCI(2), 'LineWidth', 2)

hTitle = title('CLIN VELO');
hXLabel = xlabel('BIOM');
hYLabel = ylabel('IMH Vol (%LV)');

% Add text
hText = text(320, 25, ...
     sprintf('{\\itY = %0.1gX + %0.1g}', slope, itercept));
hText1 = text(320, 22, ...
     sprintf('{\\itR^2 = %0.4g}', mdl.Rsquared.Ordinary));
hText2 = text(320, 19, ...
     sprintf('{\\itp-value = %0.2g}', mdl.Coefficients.pValue(2)));
% Add legend
hLegend = legend([hData, hModel, hCI(1)], ...
    'Observed', 'Regression line', '95% CI', ...
    'Location', 'NorthWest');

% Adjust font
set(gca, 'FontName', 'Helvetica')
set([hTitle, hXLabel, hYLabel, hText, hText1, hText2], 'FontName', 'AvantGarde')
% set([hTitle, hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hLegend, gca], 'FontSize', 12)
set([hXLabel, hYLabel, hText, hText1, hText2], 'FontSize', 12)
% set([hXLabel, hYLabel], 'FontSize', 10)
set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')

% Adjust axes properties
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', -10:5:50, ...
    'LineWidth', 1.5);
ylim([-10 50]); %xlim([0 800])

%% LG-2 vs VELO
BIOM_2_v = A_preclin_velo.BIOMPRE;
IMH_2_v = A_preclin_velo.IMH;
mdl = fitlm(BIOM_2_v, IMH_2_v);
ci = coefCI(mdl);

figure('Position', [0 100 400 400]);
itercept = mdl.Coefficients.Estimate(1);
slope = mdl.Coefficients.Estimate(2);
x0 = [0 1200];
y0 = slope * x0 + itercept;
yL = ci(2,1) * x0 + ci(1,1);
yU = ci(2,2) * x0 + ci(1,2);
hData = plot(BIOM_2_v, IMH_2_v,  'o');
hModel = line(x0, y0);
hCI(1) = line(x0, yL);
hCI(2) = line(x0, yU);
patch([x0 fliplr(x0)], [yL fliplr(yU)], 'r', 'FaceColor', [1 0 0], 'FaceAlpha', 0.1)
set(hData, 'LineStyle', 'none', 'Marker', '.')
set(hModel, 'LineStyle', '-', 'Color', 'r')
set(hCI(1), 'LineStyle', '-.', 'Color', [0 .5 0])
set(hCI(2), 'LineStyle', '-.', 'Color', [0 .5 0])

set(hData, 'Marker', 'o', 'MarkerSize', 10, ...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', [.75 .75 1])
set(hModel, 'LineWidth', 3)
set(hCI(1), 'LineWidth', 2)
set(hCI(2), 'LineWidth', 2)

hTitle = title('PRECLIN VELO');
hXLabel = xlabel('BIOM');
hYLabel = ylabel('IMH Vol (%LV)');

% Add text
hText = text(90, 10, ...
     sprintf('{\\itY = %0.1gX - %0.1g}', slope, abs(itercept)));
hText1 = text(90, 9, ...
     sprintf('{\\itR^2 = %0.4g}', mdl.Rsquared.Ordinary));
hText2 = text(90, 8, ...
     sprintf('{\\itp-value < 10^{-3}}', mdl.Coefficients.pValue(2)));
% Add legend
hLegend = legend([hData, hModel, hCI(1)], ...
    'Observed', 'Regression line', '95% CI', ...
    'Location', 'NorthWest');

% Adjust font
set(gca, 'FontName', 'Helvetica')
set([hTitle, hXLabel, hYLabel, hText, hText1, hText2], 'FontName', 'AvantGarde')
% set([hTitle, hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hLegend, gca], 'FontSize', 12)
set([hXLabel, hYLabel, hText, hText1, hText2], 'FontSize', 12)
% set([hXLabel, hYLabel], 'FontSize', 10)
set(hTitle, 'FontSize', 12, 'FontWeight' , 'bold')

% Adjust axes properties
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', 0:5:25, ...
    'LineWidth', 1.5);
ylim([0 25]); %xlim([0 800])

%% ROC - 1
figure('Position', [100 0 400 400]);
BIOM_1 = A_clin_conc.BIOM;
IMH_1 = A_clin_conc.IMH;
IMH_1_label1 = IMH_1 > 3;
IMH_1_label2 = IMH_1 > 4;

[X,Y,T,AUC,OPTROCPT] = perfcurve(IMH_1_label1, BIOM_1, 1);
scatter(X, Y, 'filled');
hold on;
plot(X,Y, 'LineWidth', 2);
xlabel('FPR');
ylabel('TPR');
text(0.5,0.5,num2str(round(AUC, 2)),'FontSize',24);
% title('ROC for Classification')
set(gca, 'FontSize', 16);
line([0 1], [0 1], 'Color', [0.5 0.5 0.5], 'LineStyle', '--');

figure('Position', [100 0 400 400]);
[X,Y,T,AUC,OPTROCPT] = perfcurve(IMH_1_label2, BIOM_1, 1);
scatter(X, Y, 'filled');
hold on;
plot(X,Y, 'LineWidth', 2);
xlabel('FPR');
ylabel('TPR');
text(0.5,0.5,num2str(round(AUC, 2)),'FontSize',24);
% title('ROC for Classification')
set(gca, 'FontSize', 16);
line([0 1], [0 1], 'Color', [0.5 0.5 0.5], 'LineStyle', '--');

%% ROC - 1 - PRECON
figure('Position', [100 0 400 400]);
BIOM_2 = A_preclin_conc.BIOMPRE;
IMH_2 = A_preclin_conc.IMH;
IMH_2_label1 = IMH_2 > 3;
IMH_2_label2 = IMH_2 > 4;

[X,Y,T,AUC,OPTROCPT] = perfcurve(IMH_2_label1, BIOM_2, 1);
scatter(X, Y, 'filled');
hold on;
plot(X,Y, 'LineWidth', 2);
xlabel('FPR');
ylabel('TPR');
text(0.5,0.5,num2str(round(AUC, 2)),'FontSize',24);
% title('ROC for Classification')
set(gca, 'FontSize', 16);
line([0 1], [0 1], 'Color', [0.5 0.5 0.5], 'LineStyle', '--');

figure('Position', [100 0 400 400]);
[X,Y,T,AUC,OPTROCPT] = perfcurve(IMH_2_label2, BIOM_2, 1);
scatter(X, Y, 'filled');
hold on;
plot(X,Y, 'LineWidth', 2);
xlabel('FPR');
ylabel('TPR');
text(0.5,0.5,num2str(round(AUC, 2)),'FontSize',24);
% title('ROC for Classification')
set(gca, 'FontSize', 16);
line([0 1], [0 1], 'Color', [0.5 0.5 0.5], 'LineStyle', '--');


%% ROC - 1 - VELO
figure('Position', [100 0 400 400]);
BIOM_1_v = A_clin_velo.BIOM;
IMH_1_v = A_clin_velo.IMH;
IMH_1_v_label1 = IMH_1_v > 3;
IMH_1_v_label2 = IMH_1_v > 4;

[X,Y,T,AUC,OPTROCPT] = perfcurve(IMH_1_v_label1, BIOM_1_v, 1);
scatter(X, Y, 'filled');
hold on;
plot(X,Y, 'LineWidth', 2);
xlabel('FPR');
ylabel('TPR');
text(0.5,0.5,num2str(round(AUC, 2)),'FontSize',24);
% title('ROC for Classification')
set(gca, 'FontSize', 16);
line([0 1], [0 1], 'Color', [0.5 0.5 0.5], 'LineStyle', '--');

figure('Position', [100 0 400 400]);
[X,Y,T,AUC,OPTROCPT] = perfcurve(IMH_1_v_label2, BIOM_1_v, 1);
scatter(X, Y, 'filled');
hold on;
plot(X,Y, 'LineWidth', 2);
xlabel('FPR');
ylabel('TPR');
text(0.5,0.5,num2str(round(AUC, 2)),'FontSize',24);
% title('ROC for Classification')
set(gca, 'FontSize', 16);
line([0 1], [0 1], 'Color', [0.5 0.5 0.5], 'LineStyle', '--');

%% ROC - 2 - VELO
figure('Position', [100 0 400 400]);
BIOM_2_v = A_preclin_velo.BIOMPRE;
IMH_2_v = A_preclin_velo.IMH;
IMH_2_v_label1 = IMH_2_v > 3;
IMH_2_v_label2 = IMH_2_v > 4;

[X,Y,T,AUC,OPTROCPT] = perfcurve(IMH_2_v_label1, BIOM_2_v, 1);
scatter(X, Y, 'filled');
hold on;
plot(X,Y, 'LineWidth', 2);
xlabel('FPR');
ylabel('TPR');
text(0.5,0.5,num2str(round(AUC, 2)),'FontSize',24);
% title('ROC for Classification')
set(gca, 'FontSize', 16);
line([0 1], [0 1], 'Color', [0.5 0.5 0.5], 'LineStyle', '--');

figure('Position', [100 0 400 400]);
[X,Y,T,AUC,OPTROCPT] = perfcurve(IMH_2_v_label2, BIOM_2_v, 1);
scatter(X, Y, 'filled');
hold on;
plot(X,Y, 'LineWidth', 2);
xlabel('FPR');
ylabel('TPR');
text(0.5,0.5,num2str(round(AUC, 2)),'FontSize',24);
% title('ROC for Classification')
set(gca, 'FontSize', 16);
line([0 1], [0 1], 'Color', [0.5 0.5 0.5], 'LineStyle', '--');
FPR = OPTROCPT(1)
TPR = OPTROCPT(2)
% NPV and PPV