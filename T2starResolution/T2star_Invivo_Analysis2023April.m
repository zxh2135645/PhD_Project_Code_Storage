clear all;
close all;

% 20P10_Exvivo7, 20P11_Exvivo6, 18P90, 18P93, 20P40, 18P92, 18P94_Exvivo3, 18P95, 17P73, 20P48
width_array = [2.40,      0.91,  1.21,  0.53,  0.56,  0.59,          0.92,  0.59,  1.10,  0.69];
SNR_array =   [9.48,	  5.20,	 5.63,	4.69,  7.05,  6.31,	         8.52,	7.80,  9.54,  11.04];

bias_array = [0.160714286, 0.404761905, 0.094147583, -1, -1, -1, -0.090909091, 3.252873563, -0.568788501, 0.916167665];

invivo_array = [2.60, 1.18, 4.3, 0, 0, 0, 0.7, 3.7, 2.1, 3.2];
gt_array = [2.24	0.84	3.93	1.23	0.56	0.9	0.77	0.87	4.87	1.67];
figure();
p1 = plot3(width_array, bias_array, SNR_array, 'o');

%%
sz = 144;
% figure('Position', [0 100 400 400]);
% p1 = scatter(width_array, SNR_array, sz, bias_array,'filled');

bias_array_v2 = bias_array;
bias_array_v2(bias_array>1) = 1;

% bias_array but to penalize false negative more
bias_array_v3 = bias_array;
bias_array_v3(bias_array == -1) = -5;


figure('Position', [0 100 400 400]);
p1 = scatter(width_array, bias_array, sz, SNR_array, 'filled');
caxis([5 15]);
hold on;
yline([0 0]);

figure('Position', [0 100 400 400]);
p1 = scatter(width_array, abs(bias_array), sz, SNR_array, 'filled');
caxis([5 15]);

%% Inverse proportional fitting
sz = 192;
y = abs(bias_array);
x = width_array;
yfit = @(b,x) b(1)./(x + b(2)) + b(3);  % Objective Function
CF = @(b) sum((y-yfit(b,x)).^2);        % Cost Function
b0 = [1, 0, 0];                      % Initial Parameter Estimates
[B, fv] = fminsearch(CF, b0); 

X = 0.35:0.1:2.5;
figure('Position', [0 100 400 400]);
p1 = scatter(width_array, abs(bias_array), sz, SNR_array, 'filled');
caxis([5 15]);
hold on;
plot(X,yfit(B,X), '--', 'LineWidth', 3);
set(gca, 'YDir','reverse');
ylim([-0.2 3.5]);
Rsq2 = 1 - sum((y - yfit(B,x)).^2)/sum((y - mean(y)).^2) % R^2 0.3290

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 2)
%% Version 2 
% yfit = @(b,x) b(1)./(x);  % Objective Function
% CF = @(b) sum((y-yfit(b,x)).^2);        % Cost Function
% b0 = [1];                      % Initial Parameter Estimates
% [B, fv] = fminsearch(CF, b0); 
% 
% X = 0.3:0.1:2.5;
% figure('Position', [0 100 400 400]);
% p1 = scatter(width_array, abs(bias_array), sz, SNR_array, 'filled');
% caxis([5 15]);
% hold on;
% plot(X,yfit(B,X));
sz = 192;
y = abs(bias_array);
x = width_array;
yfit = @(b,x) b(1)./(x + b(2)) + b(3);  % Objective Function
CF = @(b) sum((y-yfit(b,x)).^2);        % Cost Function
b0 = [1, 0, 0];                      % Initial Parameter Estimates
[B, fv] = fminsearch(CF, b0); 

X = 0.5:0.1:2.5;
figure('Position', [0 100 400 400]);
p1 = scatter(width_array, abs(bias_array), sz, SNR_array, 'filled');
caxis([5 15]);
hold on;
plot(X,yfit(B,X), '--', 'LineWidth', 3);

Rsq2 = 1 - sum((y - yfit(B,x)).^2)/sum((y - mean(y)).^2) % R^2 0.3290


%% Version 3
sz = 192;

y = abs(bias_array);
x = width_array;
yfit = @(b,x) b(1)./(x + b(2)) + b(3);  % Objective Function
CF = @(b) sum((y-yfit(b,x)).^2);        % Cost Function
b0 = [1, 0, 0];                      % Initial Parameter Estimates
[B, fv] = fminsearch(CF, b0); 



X = 0.4:0.1:2.5;
figure('Position', [0 100 400 400]);
p1 = scatter(width_array, abs(bias_array), sz, 'filled');
caxis([5 15]);
hold on;
plot(X,yfit(B,X), '--k', 'LineWidth', 3);
set(gca, 'YDir','reverse', 'XAxisLocation','top');
ylim([-0.2 4]);
Rsq2 = 1 - sum((y - yfit(B,x)).^2)/sum((y - mean(y)).^2) % R^2 0.3290

set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 2)
%% Clustering into two groups
[idx,C] = kmeans(abs(bias_array'), 3)
[idx,C] = kmeans(abs(cat(1, bias_array, width_array)'), 3)
% Cluster 1
c1 = [1 2 3 7 9];
c2 = [4 5 6 8 10];
mean(abs(bias_array(c1)))
mean(abs(bias_array(c2)))
std(abs(bias_array(c1)))
std(abs(bias_array(c2)))

mean(abs(width_array(c1)))
mean(abs(width_array(c2)))
std(abs(width_array(c1)))
std(abs(width_array(c2)))
%%  20P40 varying resolution
res_mat =    [0.9       1.2     1.6     1.9
              0.9       1.2     1.6     1.9
              0.9       1.2     1.6     1.9
              0.9       1.2     1.6     1.9];


res_mat_voxel = [0.9; 1.2; 1.6; 1.9] * [2 4 6 8];

bias_mat = [-1              -1              -1              -1    
            -1              -0.107142857    -0.642857143    -1
            -1              0.607142857     -0.464285714    -1
            0.071428571     -0.464285714    -1              -1];

SNR_mat = [0.45503212  0.718887262  1.646226415  3.035294118
           0.719665272  1.669172932  3.265957447  3.038834951
           0.928977273  3.467391304  6.419354839  3.358974359
           2.040650407  3.886792453  3.578947368  6.90625];


figure();
% p1 = plot3(res_mat, bias_mat, SNR_mat, 'o');
p1 = plot(res_mat_voxel, bias_mat);
hold on;
yline(0);

figure();
% p1 = plot3(res_mat, bias_mat, SNR_mat, 'o');
p1 = plot(res_mat_voxel, abs(bias_mat));
hold on;
yline(0);

figure();
p1 = plot(res_mat(1,:), abs(bias_mat(1,:)));
hold on;
p2 = plot(res_mat(2,:), abs(bias_mat(2,:)));
p3 = plot(res_mat(3,:), abs(bias_mat(3,:)));
p4 = plot(res_mat(4,:), abs(bias_mat(4,:)));


figure();
imagesc(abs(bias_mat)); colormap(brewermap([],'*RdYlBu'));

sz = 144;
figure();
p1 = scatter(res_mat(1,:), abs(bias_mat(1,:)), sz, SNR_mat(1,:), 'filled');
hold on;
p2 = scatter(res_mat(2,:), abs(bias_mat(2,:)), sz, SNR_mat(2,:), 'filled');
p3 = scatter(res_mat(3,:), abs(bias_mat(3,:)), sz, SNR_mat(3,:), 'filled');
p4 = scatter(res_mat(4,:), abs(bias_mat(4,:)), sz, SNR_mat(4,:), 'filled');
colormap(brewermap([],'*RdYlBu'));
caxis([2 6])

%% Bland-Altman (width vs bias)
addpath('../function/BlandAltman/');
territories = {'Invivo'};
nterritories = length(territories);

% Patient states during measurement
states = {'Volume (%)'};
nstates = length(states);

data1 = cat(3,invivo_array');
data2 = cat(3, gt_array');

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
