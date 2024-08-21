clear all;
close all;

%% 
x_sofia = [3.7	3.5	3.3	2.5	2.1	1.7	1.2	0.9	0.6	0.3	0.3	0.3	0.1	0.2	0.3 nan];
y_sofia = [nan, 6, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,	0, nan, nan,	0];

x_lisbon = [12.5	14.5	10.4	9.6	8.6	7.1	7.9	5.4	5	4	3.6	2.7	2.3	2	2.2 nan];
y_lisbon = [nan, 19, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,	2.4, nan, nan,	1.8];

x_paris = [5.6	4.9	4.2	4.5	3.7	2.8	2.2	2.8	2.6	2.1	2.5	1.7	1.6	1.6	1.5 nan];
y_paris = [nan, 9.1, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,	0.6, nan, nan,	0.8];

x_jesse = [15.2	14.5	14.7	13.3	13.5	13.1	12.8	12.6	12.3	12	11.4	11.5	12	11.9	10.4 nan];
y_jesse = [nan, 15.5, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,	12.3, nan, nan,	10.3];

% Slice 3
x_george = [22.6	22.2	21.5	21.4	20	21	20.9	20.8	20.7	20.5	20.4	18.6	19.8	18.5	18 nan];
y_george = [nan, 21.9, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,	15.2, nan, nan,	13.5];

% Slice 2
x_george = [20	18.5	16.8	16	16.1	15.7	15.4	15.6	14.5	14.3	12.3	13	13.1	12.6	11.9 nan];
y_george = [nan, 21.9, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,	15.2, nan, nan,	13.5];

x_chili = [6.6	4.3	4.2	3.8	4	3.8	3.7	3.8	3.9	3.7	3.7	3.6	3.2	3.2	2.8 nan];
y_chili = [nan, 9.7, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,	4.6, nan, nan,	3];

x_nutmeg = [12.2	11.2	11.1	10	8.4	8.3	6.5	5.4	5.6	4.6	4.8	3.9	3.5	3.4	3.2 nan];
y_nutmeg = [nan, 10.6, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,	6.1, nan, nan,	3.6];

x_ginger = [20.2	19.1	16	14.4	14.8	13.2	11.8	11	10.8	8.1	10.2	5.9	6.3	4.8	5.5 nan];
y_ginger = [nan, 10.4, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,	6.2, nan, nan,	5.3];

x_dave = [15.9	12.1	12	10.8	10.1	9.2	9.1	7.8	6.6	6.5	6.1	5.8	5.6	4.8	4.5 nan];
y_dave = [nan, 12, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,	4.5, nan, nan,	2.9];

x_carlos = [39.9	38.5	37.6	37.4	36.9	34.8	35	34.7	34	33.3	34.2	33.7	33.6	33.4	32.5 nan];
y_carlos = [nan, 36.5, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,	32.9, nan, nan,	30.8];

x_paprika = [9.2	8.9	8.5	7.4	6.6	5.7	4.9	4.1	4.6	4.1	3.8	3.2	3.3	3	3 nan];
y_paprika = [nan, 9, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,	4.7, nan, nan,	2.1];

x_cinnamon = [3.9	3.4	2.2	2.4	2.2	2	2.5	1.9	1.5	1.8	1.8	1.7	1.8	1.4	1.8 nan];
y_cinnamon = [nan, 3.6, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,	1.4, nan, nan,	0.5];


%% 
figure();
subplot(3,4,1);
plot(x_sofia, 'LineWidth', 1.5);
hold on;
plot(y_sofia, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('Sofia');

subplot(3,4,2);
plot(x_lisbon, 'LineWidth', 1.5);
hold on;
plot(y_lisbon, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('Lisbon');

subplot(3,4,3);
plot(x_paris, 'LineWidth', 1.5);
hold on;
plot(y_paris, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('Paris');

subplot(3,4,4);
plot(x_jesse, 'LineWidth', 1.5);
hold on;
plot(y_jesse, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('Jesse');

subplot(3,4,5);
plot(x_george, 'LineWidth', 1.5);
hold on;
plot(y_george, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('George');

subplot(3,4,6);
plot(x_chili, 'LineWidth', 1.5);
hold on;
plot(y_chili, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('Chili');

subplot(3,4,7);
plot(x_nutmeg, 'LineWidth', 1.5);
hold on;
plot(y_nutmeg, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('Nutmeg');

subplot(3,4,8);
plot(x_ginger, 'LineWidth', 1.5);
hold on;
plot(y_ginger, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('Ginger');

subplot(3,4,9);
plot(x_dave, 'LineWidth', 1.5);
hold on;
plot(y_dave, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('Dave');

subplot(3,4,10);
plot(x_carlos, 'LineWidth', 1.5);
hold on;
plot(y_carlos, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('Carlos');

subplot(3,4,11);
plot(x_paprika, 'LineWidth', 1.5);
hold on;
plot(y_paprika, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('Paprika');

subplot(3,4,12);
plot(x_cinnamon, 'LineWidth', 1.5);
hold on;
plot(y_cinnamon, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('Cinnamon');

%%

figure();
subplot(3,4,1);
plot(medfilt1(x_sofia,3), 'LineWidth', 1.5);
hold on;
plot(y_sofia, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('Sofia');

subplot(3,4,2);
plot(medfilt1(x_lisbon,3), 'LineWidth', 1.5);
hold on;
plot(y_lisbon, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('Lisbon');

subplot(3,4,3);
plot(medfilt1(x_paris,3), 'LineWidth', 1.5);
hold on;
plot(y_paris, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('Paris');

subplot(3,4,4);
plot(medfilt1(x_jesse,3), 'LineWidth', 1.5);
hold on;
plot(y_jesse, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('Jesse');

subplot(3,4,5);
plot(medfilt1(x_george,3), 'LineWidth', 1.5);
hold on;
plot(y_george, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('George');

subplot(3,4,6);
plot(medfilt1(x_chili,3), 'LineWidth', 1.5);
hold on;
plot(y_chili, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('Chili');

subplot(3,4,7);
plot(medfilt1(x_nutmeg,3), 'LineWidth', 1.5);
hold on;
plot(y_nutmeg, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('Nutmeg');

subplot(3,4,8);
plot(medfilt1(x_ginger,3), 'LineWidth', 1.5);
hold on;
plot(y_ginger, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('Ginger');

subplot(3,4,9);
plot(medfilt1(x_dave,3), 'LineWidth', 1.5);
hold on;
plot(y_dave, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('Dave');

subplot(3,4,10);
plot(medfilt1(x_carlos,3), 'LineWidth', 1.5);
hold on;
plot(y_carlos, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('Carlos');

subplot(3,4,11);
plot(medfilt1(x_paprika,3), 'LineWidth', 1.5);
hold on;
plot(y_paprika, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('Paprika');

subplot(3,4,12);
plot(medfilt1(x_cinnamon,3), 'LineWidth', 1.5);
hold on;
plot(y_cinnamon, 'o', 'MarkerSize', 10, 'LineWidth',2);
title('Cinnamon');

%%
n = 2;
min2_lrt = [x_sofia(n), x_lisbon(n), x_paris(n), x_jesse(n), x_george(n), x_chili(n), x_nutmeg(n), x_ginger(n), x_dave(n), x_carlos(n), x_paprika(n), x_cinnamon(n)];
min2_ege = [y_sofia(n), y_lisbon(n), y_paris(n), y_jesse(n), y_george(n), y_chili(n), y_nutmeg(n), y_ginger(n), y_dave(n), y_carlos(n), y_paprika(n), y_cinnamon(n)];

n = 13;
min13_lrt = [x_sofia(n), x_lisbon(n), x_paris(n), x_jesse(n), x_george(n), x_chili(n), x_nutmeg(n), x_ginger(n), x_dave(n), x_carlos(n), x_paprika(n), x_cinnamon(n)];
min13_ege = [y_sofia(n), y_lisbon(n), y_paris(n), y_jesse(n), y_george(n), y_chili(n), y_nutmeg(n), y_ginger(n), y_dave(n), y_carlos(n), y_paprika(n), y_cinnamon(n)];

color_cell1 = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell2 = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};

%% MVO
addpath('../function/BlandAltman/');
% Heart muscle territories per patient
territories = {''};
nterritories = length(territories);

% Patient states during measurement
states = {'WK1'};
nstates = length(states);

% Data preparation

% Baseline data with noise
data2_acute = cat(3, min2_lrt(:));
data1_acute = cat(3, min2_ege(:));

color_palette = {[101,101,101]/255, [51,51,51]/255};
color_palette_marker = [[206, 225, 252]/255; [138, 152, 227]/255];

% BA plot paramters
tit = 'EGE MVO Size Comparison'; % figure title
gnames = {territories, states}; % names of groups in data {dimension 1 and 2}
label = {'LRT-MVO','CMR-MVO','%'}; % Names of data sets
corrinfo = {'r2', 'P', 'eq'}; % stats to display of correlation scatter plot
BAinfo = {'RPC(%)','ks'}; % stats to display on Bland-ALtman plot
limits = 'auto'; % how to set the axes limits
if 0 % colors for the data sets may be set as:
	colors = 'br';      % character codes
else
	colors = [color_palette{1};... % or RGB triplets
		      color_palette{2}];
end

% axesLimits = [40 90 40 90];

% BA plot paramters
[cr_acute, fig_acute, statsStruct_acute] = BlandAltman_Customized(data1_acute, data2_acute,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors, 'showFitCI',' on', 'MarkerSize', 48, 'MarkerFaceColors', color_palette_marker);

%% MVO
addpath('../function/BlandAltman/');
% Heart muscle territories per patient
territories = {''};
nterritories = length(territories);

% Patient states during measurement
states = {'WK1'};
nstates = length(states);

% Data preparation

% Baseline data with noise
data2_acute = cat(3, min13_lrt(:));
data1_acute = cat(3, min13_ege(:));

color_palette = {[51,51,51]/255, [101,101,101]/255};
color_palette_marker = [[138, 152, 227]/255; [206, 225, 252]/255];

% BA plot paramters
tit = 'EGE MVO Size Comparison'; % figure title
gnames = {territories, states}; % names of groups in data {dimension 1 and 2}
label = {'LRT-MVO','CMR-MVO','%'}; % Names of data sets
corrinfo = {'r2', 'P', 'eq'}; % stats to display of correlation scatter plot
BAinfo = {'RPC(%)','ks'}; % stats to display on Bland-ALtman plot
limits = 'auto'; % how to set the axes limits
if 0 % colors for the data sets may be set as:
	colors = 'br';      % character codes
else
	colors = [color_palette{1};... % or RGB triplets
		      color_palette{2}];
end
% axesLimits = [40 90 40 90];

% BA plot paramters
[cr_acute, fig_acute, statsStruct_acute] = BlandAltman_Customized(data1_acute, data2_acute,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors, 'showFitCI',' on', 'MarkerSize', 48, 'MarkerFaceColors', color_palette_marker);


%% For Publication (Dave)
color_cell1 = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell2 = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};
x_array = 1:18;
x_dave_lge = [15.9,12.1,12,10.8,10.1,9.2,9.1,7.8,6.6,6.5,6.1,5.8,5.6,4.8,4.5,nan,nan,nan];
y_dave_lge = [nan, 12, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, 4.5, nan, nan,	nan, nan, 2.9];

figure('Position', [0 100 400 600]);

%plot(x_dave_lge, '-', 'LineWidth', 2);
x = x_array(1:15).';
y = x_dave_lge(1:15).';
g = fittype('a-b*exp(-c*x)');
f0 = fit(x,y,g,'StartPoint',[[ones(size(x)), -exp(-x)]\y; 1]);

xx1 = linspace(1, 15, 15);
plot(xx1,f0(xx1),'-k', 'LineWidth', 3);
hold on;
xx2 = linspace(15, 20, 6);
plot(xx2,f0(xx2),'--k', 'LineWidth', 3);


ph1 = scatter(x_array, x_dave_lge, 192, 'o', 'LineWidth', 3, 'MarkerEdgeColor', color_cell1{4}, 'MarkerFaceColor', color_cell1{3});
ph2 = scatter(x_array, y_dave_lge, 192, '^',   'LineWidth', 3, 'MarkerEdgeColor', color_cell2{4}, 'MarkerFaceColor', color_cell2{2});
set(ph1, 'MarkerFaceAlpha', 0.5); 
set(ph2, 'MarkerFaceAlpha', 0.5); 

title('Dave');
grid on;
xlim([0 20]); ylim([0 18]);
set(gca,'box','off');

%% For Publication (Chili)
x_chili_lge = [6.6,4.3,4.2,3.8,4,3.8,3.7,3.8,3.9,3.7,3.7,3.6,3.2,3.2,2.8,nan,nan,nan];
y_chili_lge = [nan, 9.7, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan,	4.6, nan, nan, nan, nan, 3];

figure('Position', [0 100 400 600]);

%plot(x_dave_lge, '-', 'LineWidth', 2);
x = x_array(1:15).';
y = x_chili_lge(1:15).';
g = fittype('a-b*exp(-c*x)');
f0 = fit(x,y,g,'StartPoint',[[ones(size(x)), -exp(-x)]\y; 1]);

xx1 = linspace(1, 15, 15);
plot(xx1,f0(xx1),'-k', 'LineWidth', 3);
hold on;
xx2 = linspace(15, 20, 6);
plot(xx2,f0(xx2),'--k', 'LineWidth', 3);


ph1 = scatter(x_array, x_chili_lge, 192, 'o', 'LineWidth', 3, 'MarkerEdgeColor', color_cell1{4}, 'MarkerFaceColor', color_cell1{3});
ph2 = scatter(x_array, y_chili_lge, 192, '^',   'LineWidth', 3, 'MarkerEdgeColor', color_cell2{4}, 'MarkerFaceColor', color_cell2{2});
set(ph1, 'MarkerFaceAlpha', 0.5); 
set(ph2, 'MarkerFaceAlpha', 0.5); 

title('Chili');
grid on;
xlim([0 20]); ylim([0 18]);
set(gca,'box','off');
