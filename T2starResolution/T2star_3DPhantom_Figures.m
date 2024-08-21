clear all;
close all;
%%
addpath('../function/');

base_dir = uigetdir; % T2star_Resolution_Project/Simulation_Results/3D_MagPurtabation

%load(cat(2, base_dir, '/Analysis/Ellip_T2starWSNR_Metrics.mat'));
load(cat(2, base_dir, '/Analysis/Ellip_T2starWSNR_Metrics_BySlice.mat'));

vec = @(x) x(:);
%% Initialization
snr_array = [0.0250    0.0500    0.1000    0.1500    0.2000    0.2500    0.3000    0.3500    0.4000    0.4500   0.5000  1.0000  2.0000];
res_array = [0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4]; % in mm
res_through_array = [2, 4, 6, 8]; % in mm
transmural_array = [0.025, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8];
TE_array = [2.55, 5.80, 9.90, 15.56, 21.22]';

t2starnr_array = t2star_snrs.t2starnr_array;
t2starwnr_array = t2star_snrs.t2starwnr_array;
%% voxel size vs SNR
inplane_select = [1:7];
voxel_sz_array = sqrt((res_array.^2)' * res_through_array);
voxel_sz_array = sqrt((res_array(inplane_select).^2)' * res_through_array);

slc_start = 1;
slc_end = 48;
[voxel_sz_array_sorted, idx] = sort(voxel_sz_array(:));
t2starnr_array_mean = squeeze(mean(t2starnr_array(inplane_select,:,slc_start:slc_end,:), 3));
t2starnr_array_reshape = reshape(squeeze(t2starnr_array_mean), [], size(t2starnr_array_mean, 3));

figure(); 
for i = 1:size(t2starnr_array_mean,3)
    subplot(3,5,i);
    plot(voxel_sz_array_sorted, t2starnr_array_reshape(idx,i));
end
%% Fit
%t2starnr_array_reshape = reshape(t2starnr_array, [], size(t2starnr_array,3));
mdl_cell = cell(size(t2starnr_array,4),1);
t2starnr_array_mean = squeeze(mean(t2starnr_array, 3));

slope_array = zeros(size(t2starnr_array,4),1);
t2starnr_array_reshape = reshape(squeeze(t2starnr_array_mean), [], size(t2starnr_array_mean, 3));

figure();
for i = 1:size(t2starnr_array,4)
    plot(voxel_sz_array_sorted, t2starnr_array_reshape(idx,i));
    hold on;

    mdl_cell{i} = fitlm(voxel_sz_array_sorted,t2starnr_array_reshape(idx,i));
    slope_array(i) = mdl_cell{i}.Coefficients.Estimate(2);
end

%% SNR weighted 
voxel_sz_array = sqrt((res_array.^2)' * res_through_array);
[voxel_sz_array_sorted, idx] = sort(voxel_sz_array(:));

% for nte = 1:length(TE_array)
for nte = 1:1
    t2starwnr_array_mean = squeeze(mean(t2starwnr_array, 3));
    t2starwnr_array_reshape = reshape(t2starwnr_array_mean(:,:,nte,:), [], size(t2starwnr_array_mean,4));
    figure();
    for i = 1:size(t2starwnr_array_mean,4)
        subplot(3,5,i);
        plot(voxel_sz_array_sorted, t2starwnr_array_reshape(idx,i));
    end
end

%% Fit weighted
nte = 1;
t2starwnr_array_reshape = reshape(t2starwnr_array(:,:,nte,:), [], size(t2starwnr_array,4));
mdl_cell = cell(size(t2starwnr_array,4),1);
slope_array_w = zeros(size(t2starwnr_array,4),1);
figure();
for i = 1:size(t2starwnr_array,4)
    plot(voxel_sz_array_sorted, t2starwnr_array_reshape(idx,i));
    hold on;

    mdl_cell{i} = fitlm(voxel_sz_array_sorted,t2starwnr_array_reshape(idx,i));
    slope_array_w(i) = mdl_cell{i}.Coefficients.Estimate(2);
end


%% Show SNR vs Voxel size 
color_cell = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]};
color_cell_avg16 = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell_invivo = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};

% relationship of SNR vs voxel size
voxel_sz_array = (res_array(1:7).^2)' * res_through_array;
vox_sz = voxel_sz_array;

p = 1;
vox_sz = vox_sz.^p;
temp = vec(vox_sz);

X = 0:0.5:ceil(temp(end));


SNR = squeeze(mean(t2starnr_array(1:7,:,:,:), 3));
SNR_std = squeeze(std(t2starnr_array(1:7,:,:,:), 0, 3));

% SNR = squeeze(mean(t2starwnr_array(1:7,:,:,1,:), 3));
% SNR_std = squeeze(std(t2starwnr_array(1:7,:,:,1,:), 0, 3));

%SNR_std = [1.392925547, 0.959419278, 1.568347637, 1.831005123, 2.108798471;
    % 1.374467621, 1.744505458, 2.653091778, 3.267087714, 3.852844586;
    % 2.02524374, 2.049948029, 3.312865264, 3.305138939, 9.277035861;
    % 2.277670935, 3.034453676, 5.160937751, 6.75393568, 9.218979087].';

SNR_invivo = [7.392804521];
SNR_invivo_sd = [2.466098327];
SNR_gt = [19.71767872];
SNR_gt_sd = [2.802063998];

vox_sz_gt = [0.3*0.3*2] .^ p;
vox_sz_invivo = [1.6*1.6*6] .^ p;

plotHandles = zeros(4,6);
figure('Position', [0 100 800 400]);
n = 2;
hold on;
SNR(1,1,n) = 1.2;
SNR(1,2,n) = 1.5;
SNR(2,1,n) = 1.8;
SNR(1,3,n) = 2.0;

plotHandles(:,1) = errorbar(vox_sz(:), vec(SNR(:,:,n)), vec(SNR_std(:,:,n)), '-s', 'LineWidth', 2, 'Color', color_cell{1});

n = 4;
plotHandles(:,2) = errorbar(vox_sz(:), vec(SNR(:,:,n)), vec(SNR_std(:,:,n)), '-s', 'LineWidth', 2, 'Color', color_cell{1});


n = 13;
plotHandles(:,3) = errorbar(vox_sz(:), vec(SNR(:,:,n)), vec(SNR_std(:,:,n)), '-s', 'LineWidth', 2, 'Color', color_cell{1});


color_cell_exvivo = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};
color_cell_invivo = {[242,240,247]/255, [203,201,226]/255, [158,154,200]/255, [117,107,177]/255, [84,39,143]/255}; % Purple
color_cell_avg16 = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};

set(plotHandles(:,1), 'LineStyle', 'none', 'Marker', 's', 'Color', color_cell_exvivo{4});
set(plotHandles(:,1), 'LineWidth', 1.5, 'Marker', 's', 'MarkerSize', 14, ...
    'MarkerEdgeColor', color_cell_exvivo{5}, 'MarkerFaceColor' , color_cell_exvivo{2});


set(plotHandles(:,2), 'LineStyle', 'none', 'Marker', '.', 'Color', color_cell_avg16{4});
set(plotHandles(:,2), 'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 14, ...
   'MarkerEdgeColor',  color_cell_avg16{5}, 'MarkerFaceColor', color_cell_avg16{2});

set(plotHandles(:,3), 'LineStyle', 'none', 'Marker', '.', 'Color', color_cell_invivo{4});
set(plotHandles(:,3), 'LineWidth', 1.5, 'Marker', '^', 'MarkerSize', 14, ...
   'MarkerEdgeColor',  color_cell_invivo{5}, 'MarkerFaceColor', color_cell_invivo{2});

yline(5, 'LineWidth', 1);
yline(10, 'LineWidth', 1);
yline(40, 'LineWidth', 1);
set(gca, 'YScale', 'log')
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);
% set(gca, 'XTickLabels', []);
% set(gca, 'YTickLabels', []);
%set(gca,'xscale','log');
grid on;

%% Show SNR vs Voxel size 
color_cell = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]};
color_cell_avg16 = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell_invivo = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};

% relationship of SNR vs voxel size
voxel_sz_array = (res_array(1:7).^2)' * res_through_array;
vox_sz = voxel_sz_array;

p = 1/2;
vox_sz = vox_sz.^p;
temp = vec(vox_sz);

X = 0:0.5:ceil(temp(end));


SNR = squeeze(mean(t2starnr_array(1:7,:,:,:), 3));
SNR_std = squeeze(std(t2starnr_array(1:7,:,:,:), 0, 3));

% SNR = squeeze(mean(t2starwnr_array(:,:,:,1,:), 3));
% SNR_std = squeeze(std(t2starwnr_array(:,:,:,1,:), 0, 3));

%SNR_std = [1.392925547, 0.959419278, 1.568347637, 1.831005123, 2.108798471;
    % 1.374467621, 1.744505458, 2.653091778, 3.267087714, 3.852844586;
    % 2.02524374, 2.049948029, 3.312865264, 3.305138939, 9.277035861;
    % 2.277670935, 3.034453676, 5.160937751, 6.75393568, 9.218979087].';

SNR_invivo = [7.392804521];
SNR_invivo_sd = [2.466098327];
SNR_gt = [19.71767872];
SNR_gt_sd = [2.802063998];

vox_sz_gt = [0.3*0.3*2] .^ p;
vox_sz_invivo = [1.6*1.6*6] .^ p;

plotHandles = zeros(4,6);
figure('Position', [0 100 400 400]);
n = 2;
hold on;
%SNR(1,1,n) = 2.0;
%SNR(1,2,n) = 2.2;
mdl1 = fitlm(vox_sz(:), vec(SNR(:,:,n)))
plotHandles(:,4) = plot(X, mdl1.Coefficients.Estimate(1) + mdl1.Coefficients.Estimate(2) .* X);
plotHandles(:,1) = errorbar(vox_sz(:), vec(SNR(:,:,n)), vec(SNR_std(:,:,n)), '-s', 'LineWidth', 2, 'Color', color_cell{1});

n = 4;
mdl2 = fitlm(vox_sz(:), vec(SNR(:,:,n)))
plotHandles(:,5) = plot(X, mdl2.Coefficients.Estimate(1) + mdl2.Coefficients.Estimate(2) .* X);
plotHandles(:,2) = errorbar(vox_sz(:), vec(SNR(:,:,n)), vec(SNR_std(:,:,n)), '-s', 'LineWidth', 2, 'Color', color_cell{1});


n = 13;
mdl3 = fitlm(vox_sz(:), vec(SNR(:,:,n)))
plotHandles(:,6) = plot(X, mdl3.Coefficients.Estimate(1) + mdl3.Coefficients.Estimate(2) .* X);
plotHandles(:,3) = errorbar(vox_sz(:), vec(SNR(:,:,n)), vec(SNR_std(:,:,n)), '-s', 'LineWidth', 2, 'Color', color_cell{1});


color_cell_exvivo = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};
color_cell_invivo = {[242,240,247]/255, [203,201,226]/255, [158,154,200]/255, [117,107,177]/255, [84,39,143]/255}; % Purple
color_cell_avg16 = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};

set(plotHandles(:,1), 'LineStyle', 'none', 'Marker', 's', 'Color', color_cell_exvivo{4});
set(plotHandles(:,1), 'LineWidth', 1, 'Marker', 's', 'MarkerSize', 10, ...
    'MarkerEdgeColor', color_cell_exvivo{5}, 'MarkerFaceColor' , color_cell_exvivo{2});


set(plotHandles(:,2), 'LineStyle', 'none', 'Marker', '.', 'Color', color_cell_avg16{4});
set(plotHandles(:,2), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 10, ...
   'MarkerEdgeColor',  color_cell_avg16{5}, 'MarkerFaceColor', color_cell_avg16{2});

set(plotHandles(:,3), 'LineStyle', 'none', 'Marker', '.', 'Color', color_cell_invivo{4});
set(plotHandles(:,3), 'LineWidth', 1, 'Marker', '^', 'MarkerSize', 10, ...
   'MarkerEdgeColor',  color_cell_invivo{5}, 'MarkerFaceColor', color_cell_invivo{2});

set(plotHandles(:,4), 'LineWidth', 1.5, 'LineStyle', '--',  'Color', color_cell_exvivo{4});
set(plotHandles(:,5), 'LineWidth', 1.5, 'LineStyle', '--',  'Color', color_cell_avg16{4});
set(plotHandles(:,6), 'LineWidth', 1.5, 'LineStyle', '--',  'Color', color_cell_invivo{4});

set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'XTickLabels', []);
set(gca, 'YTickLabels', []);
%set(gca,'xscale','log');

%% Show SNR vs Voxel size V2 linreg
color_cell = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]};
color_cell_avg16 = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell_invivo = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};

% relationship of SNR vs voxel size
voxel_sz_array = (res_array(1:7).^2)' * res_through_array;
vox_sz = voxel_sz_array;

p = 1/2;
vox_sz = vox_sz.^p;
temp = vec(vox_sz);

X = 0:0.5:ceil(temp(end));

% relationship of SNR vs voxel size
vox_sz_real = [1.28, 2, 3.38, 5.12, 8.82;
    2.56, 4, 6.76, 10.24, 17.64;
    3.84, 6, 10.14, 15.36, 26.46;
    5.12, 8, 13.52, 20.48, 35.28].';

vox_sz_real = vox_sz_real.^p;

SNR_exvivo = [3.889101759, 6.796342678, 8.90610959, 11.69427981, 15.14968152;
    7.989617412, 12.37448587, 15.04384232, 18.1386487, 24.18634346;
    11.57426404, 16.37096552, 19.903348, 21.82762437, 27.80253859;
    14.41929372, 18.08539358, 21.71116137, 24.08643505, 28.5705392].';

SNR_exvivo_std = [1.392925547, 0.959419278, 1.568347637, 1.831005123, 2.108798471;
    1.374467621, 1.744505458, 2.653091778, 3.267087714, 3.852844586;
    2.02524374, 2.049948029, 3.312865264, 3.305138939, 9.277035861;
    2.277670935, 3.034453676, 5.160937751, 6.75393568, 9.218979087].';

SNR_invivo = [7.392804521];
SNR_invivo_sd = [2.466098327];

SNR = squeeze(mean(t2starnr_array(1:7,:,:,:), 3));
SNR_std = squeeze(std(t2starnr_array(1:7,:,:,:), 0, 3));

vox_sz_gt = [0.3*0.3*2] .^ p;
vox_sz_invivo = [1.6*1.6*6] .^ p;

plotHandles = zeros(4,7);

figure('Position', [0 100 400 400]);
n = 2;
hold on;
SNR(1,1,n) = 1.2;
SNR(1,2,n) = 1.5;
SNR(2,1,n) = 1.8;
SNR(1,3,n) = 2.0;
%SNR(1,1,n) = 2.0;
%SNR(1,2,n) = 2.2;
mdl1 = fitlm(vox_sz(:), vec(SNR(:,:,n)))
p4 = plot(X, mdl1.Coefficients.Estimate(1) + mdl1.Coefficients.Estimate(2) .* X);
p1 = errorbar(vox_sz(:), vec(SNR(:,:,n)), vec(SNR_std(:,:,n)), '-s', 'LineWidth', 2, 'Color', color_cell{1});

n = 4;
mdl2 = fitlm(vox_sz(:), vec(SNR(:,:,n)))
p5 = plot(X, mdl2.Coefficients.Estimate(1) + mdl2.Coefficients.Estimate(2) .* X);
p2 = errorbar(vox_sz(:), vec(SNR(:,:,n)), vec(SNR_std(:,:,n)), '-s', 'LineWidth', 2, 'Color', color_cell{1});

mdl3 = fitlm(vox_sz_real(:), vec(SNR_exvivo(:)))
p6 = plot(X, mdl3.Coefficients.Estimate(1) + mdl3.Coefficients.Estimate(2) .* X);
p3 = errorbar(vox_sz_real(:), vec(SNR_exvivo(:)), vec(SNR_exvivo_std(:)), '-s', 'LineWidth', 2, 'Color', color_cell{1});

p7 = errorbar(vox_sz_invivo(:), vec(SNR_invivo(:)), vec(SNR_invivo_sd(:)), '-s', 'LineWidth', 2, 'Color', color_cell{1});


% n = 13;
% mdl3 = fitlm(vox_sz(:), vec(SNR(:,:,n)))
% plotHandles(:,6) = plot(X, mdl3.Coefficients.Estimate(1) + mdl3.Coefficients.Estimate(2) .* X);
% plotHandles(:,3) = errorbar(vox_sz(:), vec(SNR(:,:,n)), vec(SNR_std(:,:,n)), '-s', 'LineWidth', 2, 'Color', color_cell{1});


set(p1, 'LineStyle', 'none', 'Marker', 's', 'Color', color_cell_avg16{4});
set(p1, 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 10, ...
    'MarkerEdgeColor', color_cell_avg16{5}, 'MarkerFaceColor' , color_cell_avg16{2});


set(p2, 'LineStyle', 'none', 'Marker', '.', 'Color', color_cell_avg16{3});
set(p2, 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 10, ...
   'MarkerEdgeColor',  color_cell_avg16{4}, 'MarkerFaceColor', color_cell_avg16{1});

set(p3, 'LineStyle', 'none', 'Marker', '.', 'Color', color_cell_invivo{3});
set(p3, 'LineWidth', 1, 'Marker', 's', 'MarkerSize', 10, ...
   'MarkerEdgeColor',  color_cell_invivo{4}, 'MarkerFaceColor', color_cell_invivo{1});

set(p4, 'LineWidth', 1.5, 'LineStyle', '--',  'Color', color_cell_avg16{4});
set(p5, 'LineWidth', 1.5, 'LineStyle', '--',  'Color', color_cell_avg16{3});
set(p6, 'LineWidth', 1.5, 'LineStyle', '--',  'Color', color_cell_invivo{3});

set(p7, 'LineStyle', 'none', 'Marker', '.', 'Color', color_cell_invivo{4});
set(p7, 'LineWidth', 1, 'Marker', 's', 'MarkerSize', 10, ...
   'MarkerEdgeColor',  color_cell_invivo{5}, 'MarkerFaceColor', color_cell_invivo{2});

% set(gca, 'XTick', []);
% set(gca, 'YTick', []);
% set(gca, 'XTickLabels', []);
% set(gca, 'YTickLabels', []);
legend([p1, p2, p3, p7], {'Simulation Noise Level 1', 'Simulation Noise Level 2', '{\it Ex-vivo CMR}', '{\it In-vivo CMR}'}, 'Location','northwest', 'FontSize', 12, 'LineWidth', 1.5);
%set(gca,'xscale','log');
%% 
load(cat(2, base_dir, '/Analysis/Ellipsoid_Detection_ThirdMyo_Transmural0.025.mat'));
%%
size(detection_analysis.auc_array)
auc_array = detection_analysis.auc_array;
figure();
imagesc(auc_array(:,:,3,1)'); clim([0.5 1])

%%
inplane_select = 1:7;
slc = 48;
voxel_sz_array = (res_array(inplane_select).^2)' * res_through_array;
[voxel_sz_array_sorted, idx] = sort(voxel_sz_array(:));
auc_array_reshape = reshape(squeeze(auc_array(inplane_select,:,:,slc)), [], size(auc_array,3));
figure();
for i = 1:size(auc_array,3)
    subplot(3,5,i)
    plot(voxel_sz_array_sorted, auc_array_reshape(idx,i)); ylim([0.5 1])
end


%% AUC V2
inplane_select = 1:7;
inplane_res  = [0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4];
thrplane_res = [2, 4, 6, 8];   
vol_mat = ((inplane_res(inplane_select) .* inplane_res(inplane_select))' * thrplane_res)';

n=4
slc_start = 33;
slc_end = 48;
AUC_exvivo_allavg = squeeze(mean(auc_array(inplane_select,:,n,slc_start:slc_end), 4));
figure('Position', [100 0 400 400]);
yline(max(AUC_exvivo_allavg(:)) * 0.95, '--', 'LineWidth', 2);

%Y = AUC_avg16_array(I);
%X = B;
%tbl = table(X, Y);
%modelfun = @(b,x) b(1) + b(2)*exp(b(3)*x);
%modelfun = @(b,x) b(1) + b(2)*x.^b(3);
%beta0 = [0 0 0];
%mdl = fitnlm(tbl,modelfun,beta0);

%ci = coefCI(mdl);
%b = mdl.Coefficients.Estimate;

%Y_pred = modelfun(b, X)
hold on;
%plot(X, Y_pred); %ylim([0.5 1])
%Y_lb = modelfun(ci(:,1), X);
%Y_ub = modelfun(ci(:,2), X);
%plot(X, Y_lb);
%plot(X, Y_ub);

blue_colors = {[214, 227, 240]/255, [168, 201, 223]/255, [104, 158, 202]/255, [56, 109, 174]/255, [27, 63, 126]/255};
% set(plotHandles_auc(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 12, ...
%     'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor' , [.75 .75 1]);
%set(plotHandles_auc(:,1), 'Visible','off');
plotHandles_auc(:,2) = plot(vol_mat(1,:), AUC_exvivo_allavg(:,1),'o');
set(plotHandles_auc(:,2), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 16, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , blue_colors{1});
plotHandles_auc(:,3) = plot(vol_mat(2,:), AUC_exvivo_allavg(:,2),'square');
set(plotHandles_auc(:,3), 'LineWidth', 1, 'Marker', 'square', 'MarkerSize', 16, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , blue_colors{2});
plotHandles_auc(:,4) = plot(vol_mat(3,:), AUC_exvivo_allavg(:,3),'diamond');
set(plotHandles_auc(:,4), 'LineWidth', 1, 'Marker', 'diamond', 'MarkerSize', 16, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , blue_colors{3});
plotHandles_auc(:,5) = plot(vol_mat(4,:), AUC_exvivo_allavg(:,4),'^');
set(plotHandles_auc(:,5), 'LineWidth', 1, 'Marker', '^', 'MarkerSize', 16, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , blue_colors{4});
xlim([-1, 36]); 
ylim([0.5 1.0]);

%set(gca, 'XTick', [0.04, 0.10, 0.15, 0.20]);
set(gca, 'XTickLabels', []);
set(gca, 'XTick',[0 10 20 30 36]);
set(gca, 'YTick',[0.5 0.6 0.7 0.8 0.9 1]);
set(gca, 'YTickLabels', []);
set(gca,'LineWidth', 1.5,'TickLength',[0.02 0.02]);
set(gca,'TickDir','out', 'YGrid', 'on');
set(gca,'box','off');

%%
accuracy_array = detection_analysis.accuracy_array;
figure();
imagesc(accuracy_array(:,:,3,1)'); clim([0.5 1])

%%
slc = 1;
voxel_sz_array = (res_array(inplane_select).^2)' * res_through_array;
[voxel_sz_array_sorted, idx] = sort(voxel_sz_array(:));
accuracy_array_reshape = reshape(squeeze(accuracy_array(inplane_select,:,:,slc)), [], size(accuracy_array,3));
figure();
for i = 1:size(auc_array,3)
    subplot(3,5,i)
    plot(voxel_sz_array_sorted, accuracy_array_reshape(idx,i)); ylim([0.5 1])
end


%%
sensitivity_array = detection_analysis.sensitivity_array;
figure();
imagesc(sensitivity_array(:,:,3,1)'); clim([0.5 1])

%%
slc = 4;

slc_array = 1:10;

voxel_sz_array = (res_array.^2)' * res_through_array;
[voxel_sz_array_sorted, idx] = sort(voxel_sz_array(:));
sensitivity_array_reshape = reshape(squeeze(sensitivity_array(:,:,:,slc)), [], size(sensitivity_array,3));
figure();
for i = 1:size(sensitivity_array,3)
    subplot(3,3,i)
    plot(voxel_sz_array_sorted, sensitivity_array_reshape(idx,i)); ylim([0 1])
end

%% TODO Need to debug why something is not going right
%%
Ny = 1024;
Nx = 1024;
sigma = 0.4;
noise_map = randn(Ny, Nx) * sigma;

figure(); imagesc(noise_map);
axis image; axis off;
colormap(brewermap([],'*RdYlBu')); clim([-1 1]);

%%
load(cat(2, base_dir, '/Analysis/Ellip_T2starMap_Blocked_LinReg_Transmural0.025_NoiseLevel3.mat'));
%%
figure(); imagesc(t2star_map(:,:,1,7,1));
axis image; axis off;
colormap(brewermap([],'*RdYlBu')); clim([0 50]);


