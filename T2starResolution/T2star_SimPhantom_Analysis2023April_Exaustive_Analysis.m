clear all;
close all;
% Load SimPhantom_2023_analysis_April.mat
addpath('../function/');
%%
base_dir = uigetdir;
fname = 'SimPhantom_2023_analysis_April_2024.mat';
%fname = 'SimPhantom_2023_analysis_April_fwhm.mat';
%fname = 'SimPhantom_2023_analysis_April_fwhm_V2.mat';
f_to_read = cat(2, base_dir, '/', fname);
load(f_to_read);

SNR_mat_04 = SimPhantom_2023_analysis_April.SNR_mat_04;
SNR_mat_06 = SimPhantom_2023_analysis_April.SNR_mat_06;
SNR_mat_08 = SimPhantom_2023_analysis_April.SNR_mat_08;
SNR_mat_10 = SimPhantom_2023_analysis_April.SNR_mat_10;
SNR_mat_15 = SimPhantom_2023_analysis_April.SNR_mat_15;
SNR_mat_20 = SimPhantom_2023_analysis_April.SNR_mat_20;
SNR_mat_25 = SimPhantom_2023_analysis_April.SNR_mat_25;
SNR_mat_30 = SimPhantom_2023_analysis_April.SNR_mat_30;

SD_mat_04 = SimPhantom_2023_analysis_April.SD_mat_04;
SD_mat_06 = SimPhantom_2023_analysis_April.SD_mat_06;
SD_mat_08 = SimPhantom_2023_analysis_April.SD_mat_08;
SD_mat_10 = SimPhantom_2023_analysis_April.SD_mat_10;
SD_mat_15 = SimPhantom_2023_analysis_April.SD_mat_15;
SD_mat_20 = SimPhantom_2023_analysis_April.SD_mat_20;
SD_mat_25 = SimPhantom_2023_analysis_April.SD_mat_25;
SD_mat_30 = SimPhantom_2023_analysis_April.SD_mat_30;

res_array = SimPhantom_2023_analysis_April.res_array;
sigma_array = SimPhantom_2023_analysis_April.sigma_array;
width_max = SimPhantom_2023_analysis_April.width_max;
accu_matt2_04 = SimPhantom_2023_analysis_April.accu_matt2_04;
accu_matt2_06 = SimPhantom_2023_analysis_April.accu_matt2_06;
accu_matt2_08 = SimPhantom_2023_analysis_April.accu_matt2_08;
accu_matt2_10 = SimPhantom_2023_analysis_April.accu_matt2_10;
accu_matt2_15 = SimPhantom_2023_analysis_April.accu_matt2_15;
accu_matt2_20 = SimPhantom_2023_analysis_April.accu_matt2_20;
accu_matt2_25 = SimPhantom_2023_analysis_April.accu_matt2_25;
accu_matt2_30 = SimPhantom_2023_analysis_April.accu_matt2_30;

auc_matt2_04 = SimPhantom_2023_analysis_April.auc_matt2_04;
auc_matt2_06 = SimPhantom_2023_analysis_April.auc_matt2_06;
auc_matt2_08 = SimPhantom_2023_analysis_April.auc_matt2_08;
auc_matt2_10 = SimPhantom_2023_analysis_April.auc_matt2_10;
auc_matt2_15 = SimPhantom_2023_analysis_April.auc_matt2_15;
auc_matt2_20 = SimPhantom_2023_analysis_April.auc_matt2_20;
auc_matt2_25 = SimPhantom_2023_analysis_April.auc_matt2_25;
auc_matt2_30 = SimPhantom_2023_analysis_April.auc_matt2_30;

dice_matt2_04 = SimPhantom_2023_analysis_April.dice_matt2_04;
dice_matt2_06 = SimPhantom_2023_analysis_April.dice_matt2_06;
dice_matt2_08 = SimPhantom_2023_analysis_April.dice_matt2_08;
dice_matt2_10 = SimPhantom_2023_analysis_April.dice_matt2_10;
dice_matt2_15 = SimPhantom_2023_analysis_April.dice_matt2_15;
dice_matt2_20 = SimPhantom_2023_analysis_April.dice_matt2_20;
dice_matt2_25 = SimPhantom_2023_analysis_April.dice_matt2_25;
dice_matt2_30 = SimPhantom_2023_analysis_April.dice_matt2_30;

sens_matt2_04 = SimPhantom_2023_analysis_April.sens_matt2_04;
sens_matt2_06 = SimPhantom_2023_analysis_April.sens_matt2_06;
sens_matt2_08 = SimPhantom_2023_analysis_April.sens_matt2_08;
sens_matt2_10 = SimPhantom_2023_analysis_April.sens_matt2_10;
sens_matt2_15 = SimPhantom_2023_analysis_April.sens_matt2_15;
sens_matt2_20 = SimPhantom_2023_analysis_April.sens_matt2_20;
sens_matt2_25 = SimPhantom_2023_analysis_April.sens_matt2_25;
sens_matt2_30 = SimPhantom_2023_analysis_April.sens_matt2_30;

accu_matt3_04 = SimPhantom_2023_analysis_April.accu_matt3_04;
accu_matt3_06 = SimPhantom_2023_analysis_April.accu_matt3_06;
accu_matt3_08 = SimPhantom_2023_analysis_April.accu_matt3_08;
accu_matt3_10 = SimPhantom_2023_analysis_April.accu_matt3_10;
accu_matt3_15 = SimPhantom_2023_analysis_April.accu_matt3_15;
accu_matt3_20 = SimPhantom_2023_analysis_April.accu_matt3_20;
accu_matt3_25 = SimPhantom_2023_analysis_April.accu_matt3_25;
accu_matt3_30 = SimPhantom_2023_analysis_April.accu_matt3_30;

auc_matt3_04 = SimPhantom_2023_analysis_April.auc_matt3_04;
auc_matt3_06 = SimPhantom_2023_analysis_April.auc_matt3_06;
auc_matt3_08 = SimPhantom_2023_analysis_April.auc_matt3_08;
auc_matt3_10 = SimPhantom_2023_analysis_April.auc_matt3_10;
auc_matt3_15 = SimPhantom_2023_analysis_April.auc_matt3_15;
auc_matt3_20 = SimPhantom_2023_analysis_April.auc_matt3_20;
auc_matt3_25 = SimPhantom_2023_analysis_April.auc_matt3_25;
auc_matt3_30 = SimPhantom_2023_analysis_April.auc_matt3_30;
%%
auc_matt3_avg = (auc_matt3_30+auc_matt3_25+auc_matt3_20+auc_matt3_15+auc_matt3_10+auc_matt3_08+auc_matt3_06+auc_matt3_04)/8;
accu_matt3_avg = (accu_matt3_30+accu_matt3_25+accu_matt3_20+accu_matt3_15+accu_matt3_10+accu_matt3_08+accu_matt3_06+accu_matt3_04)/8;
sens_matt2_avg = (sens_matt2_30+sens_matt2_25+sens_matt2_20+sens_matt2_15+sens_matt2_10+sens_matt2_08+sens_matt2_06+sens_matt2_04)/8;

figure();
imagesc(auc_matt3_avg.');axis image; axis off;
caxis([0.8 1]);
colormap(brewermap([],'*RdYlBu'));
colorbar;

%%
noise_level = 10;
figure('Position', [0 100 900 300]);
p1 = plot(SNR_mat_04(:,noise_level), sens_matt2_04(:,noise_level), 'LineWidth', 3);
hold on;
p2 = plot(SNR_mat_06(:,noise_level), sens_matt2_06(:,noise_level), 'LineWidth', 3);
p3 = plot(SNR_mat_08(:,noise_level), sens_matt2_08(:,noise_level), 'LineWidth', 3);
p4 = plot(SNR_mat_10(:,noise_level), sens_matt2_10(:,noise_level), 'LineWidth', 3);
p5 = plot(SNR_mat_15(:,noise_level), sens_matt2_15(:,noise_level), 'LineWidth', 3);
p6 = plot(SNR_mat_20(:,noise_level), sens_matt2_20(:,noise_level), 'LineWidth', 3);
p7 = plot(SNR_mat_25(:,noise_level), sens_matt2_25(:,noise_level), 'LineWidth', 3);
p8 = plot(SNR_mat_30(:,noise_level), sens_matt2_30(:,noise_level), 'LineWidth', 3);
c1 = p1.Color;c2 = p2.Color;c3 = p3.Color;c4 = p4.Color;c5 = p5.Color;c6 = p6.Color;
c7 = p7.Color; c8 = p8.Color;
ylim([0 1]); xlim([0 25]);
axis off;

noise_level = 5;
figure('Position', [0 100 900 300]);
plot(SNR_mat_04(:,noise_level), sens_matt2_04(:,noise_level),  'LineWidth', 3, 'Color', c1);
hold on;
plot(SNR_mat_06(:,noise_level), sens_matt2_06(:,noise_level),  'LineWidth', 3, 'Color', c2);
plot(SNR_mat_08(:,noise_level), sens_matt2_08(:,noise_level),  'LineWidth', 3, 'Color', c3);
plot(SNR_mat_10(:,noise_level), sens_matt2_10(:,noise_level),  'LineWidth', 3, 'Color', c4);
plot(SNR_mat_15(:,noise_level), sens_matt2_15(:,noise_level),  'LineWidth', 3, 'Color', c5);
plot(SNR_mat_20(:,noise_level), sens_matt2_20(:,noise_level),  'LineWidth', 3, 'Color', c6);
plot(SNR_mat_25(:,noise_level), sens_matt2_25(:,noise_level),  'LineWidth', 3, 'Color', c7);
plot(SNR_mat_30(:,noise_level), sens_matt2_30(:,noise_level),  'LineWidth', 3, 'Color', c8);
ylim([0 1]); xlim([0 25]);
grid on; axis off;

%% thickness vs accuracy
width_mat = repmat(width_max, [length(res_array),1]);
accu_matt3_3d = cat(3,  accu_matt3_04, accu_matt3_06, accu_matt3_08, accu_matt3_10, accu_matt3_15, accu_matt3_20, accu_matt3_25, accu_matt3_30);
accu_matt3_3d = permute(accu_matt3_3d, [1,3,2]);

sens_matt2_3d = cat(3,  sens_matt2_04, sens_matt2_06, sens_matt2_08, sens_matt2_10, sens_matt2_15, sens_matt2_20, sens_matt2_25, sens_matt2_30);
sens_matt2_3d = permute(sens_matt2_3d, [1,3,2]);

SNR_mat_3d = cat(3,  SNR_mat_04, SNR_mat_06, SNR_mat_08, SNR_mat_10, SNR_mat_15, SNR_mat_20, SNR_mat_25, SNR_mat_30);
SNR_mat_3d = permute(SNR_mat_3d, [1,3,2]);

accu_matt2_3d = cat(3,  accu_matt2_04, accu_matt2_06, accu_matt2_08, accu_matt2_10, accu_matt2_15, accu_matt2_20, accu_matt2_25, accu_matt2_30);
accu_matt2_3d = permute(accu_matt2_3d, [1,3,2]);

auc_matt2_3d = cat(3,  auc_matt2_04, auc_matt2_06, auc_matt2_08, auc_matt2_10, auc_matt2_15, auc_matt2_20, auc_matt2_25, auc_matt2_30);
auc_matt2_3d = permute(auc_matt2_3d, [1,3,2]);

auc_matt3_3d = cat(3,  auc_matt3_04, auc_matt3_06, auc_matt3_08, auc_matt3_10, auc_matt3_15, auc_matt3_20, auc_matt3_25, auc_matt3_30);
auc_matt3_3d = permute(auc_matt3_3d, [1,3,2]);

figure();
p1 = scatter(width_mat, sens_matt2_3d(:,:,10));
%% 2024 April
sz = 72;
noise_level = 5;

[max_lvl2 idx_lvl2] = max(auc_matt2_3d(:,:,noise_level));

SNR_correspond_lvl2 = zeros(1, length(idx_lvl2));
for i = 1:length(idx_lvl2)
    SNR_correspond_lvl2(i) = SNR_mat_3d(idx_lvl2(i),i,noise_level);
end


noise_level = 10;
[max_lvl3 idx_lvl3] = max(auc_matt2_3d(:,:,noise_level));
SNR_correspond_lvl3 = zeros(1, length(idx_lvl3));
for i = 1:length(idx_lvl3)
    SNR_correspond_lvl3(i) = SNR_mat_3d(idx_lvl3(i),i,noise_level);
end

figure(); 
p1 = scatter(width_mat(1,:), SNR_correspond_lvl2, sz,'filled');
hold on;
p2 = scatter(width_mat(1,:), SNR_correspond_lvl3, sz,'filled');

%%
sz = 72;
noise_level = 5;

[max_lvl2 idx_lvl2] = max(auc_matt2_3d(:,:,noise_level));


noise_level = 10;
[max_lvl3 idx_lvl3] = max(auc_matt2_3d(:,:,noise_level));


figure(); 
p1 = scatter(width_mat(1,:), max_lvl2, sz,'filled');
hold on;
p2 = scatter(width_mat(1,:), max_lvl3, sz,'filled');


idx_lvl3(3) = 7;
figure(); 
p1 = scatter(width_mat(1,:), res_array(idx_lvl2), sz,'filled');
hold on;
p2 = scatter(width_mat(1,:), res_array(idx_lvl3), sz,'filled');


%%
figure(); 
p1 = scatter(width_mat(1,:), idx_lvl2,'filled');
hold on;
p2 = scatter(width_mat(1,:), idx_lvl3,'filled');
%%
figure();
p1 = plot3(width_mat, sens_matt2_3d(:,:,10), SNR_mat_3d(:,:,10), 'o');
%%
figure();
p1 = plot3(width_mat, SNR_mat_3d(:,:,10), sens_matt2_3d(:,:,10),  'o');
%%
sz = 72;
figure();
p1 = scatter(width_mat(1,:), sens_matt2_3d(1,:,10), sz, SNR_mat_3d(1,:,10),'filled');
hold on;
p2 = scatter(width_mat(2,:), sens_matt2_3d(2,:,10), sz, SNR_mat_3d(2,:,10),'filled');
p3 = scatter(width_mat(3,:), sens_matt2_3d(3,:,10), sz, SNR_mat_3d(3,:,10),'filled');
p4 = scatter(width_mat(4,:), sens_matt2_3d(4,:,10), sz, SNR_mat_3d(4,:,10),'filled');
p5 = scatter(width_mat(5,:), sens_matt2_3d(5,:,10), sz, SNR_mat_3d(5,:,10),'filled');
p6 = scatter(width_mat(6,:), sens_matt2_3d(6,:,10), sz, SNR_mat_3d(6,:,10),'filled');
p7 = scatter(width_mat(7,:), sens_matt2_3d(7,:,10), sz, SNR_mat_3d(7,:,10),'filled');
caxis([5 15]);

sz = 72;
figure();
p1 = scatter(width_mat(1,:), sens_matt2_3d(1,:,5), sz, SNR_mat_3d(1,:,5),'filled');
hold on;
p2 = scatter(width_mat(2,:), sens_matt2_3d(2,:,5), sz, SNR_mat_3d(2,:,5),'filled');
p3 = scatter(width_mat(3,:), sens_matt2_3d(3,:,5), sz, SNR_mat_3d(3,:,5),'filled');
p4 = scatter(width_mat(4,:), sens_matt2_3d(4,:,5), sz, SNR_mat_3d(4,:,5),'filled');
p5 = scatter(width_mat(5,:), sens_matt2_3d(5,:,5), sz, SNR_mat_3d(5,:,5),'filled');
p6 = scatter(width_mat(6,:), sens_matt2_3d(6,:,5), sz, SNR_mat_3d(6,:,5),'filled');
p7 = scatter(width_mat(7,:), sens_matt2_3d(7,:,5), sz, SNR_mat_3d(7,:,5),'filled');
caxis([5 15]);

%% 4 Dimensional plot
%%
noise_level = 1:1:10;
C = zeros([size(sens_matt2_3d,2) size(sens_matt2_3d,3)]);
[X,Y] = meshgrid(noise_level, width_max);
% noise_level = 5;
figure('Position', [0 100 800 400]);
s1 = surf(X, Y, squeeze(SNR_mat_3d(1,:,:)), squeeze(sens_matt2_3d(1,:,:)), 'FaceAlpha',0.5);
hold on;
s2 = surf(X, Y, squeeze(SNR_mat_3d(2,:,:)), squeeze(sens_matt2_3d(2,:,:)), 'FaceAlpha',0.5);
s3 = surf(X, Y, squeeze(SNR_mat_3d(3,:,:)), squeeze(sens_matt2_3d(3,:,:)), 'FaceAlpha',0.5);
s4 = surf(X, Y, squeeze(SNR_mat_3d(4,:,:)), squeeze(sens_matt2_3d(4,:,:)), 'FaceAlpha',0.5);
s5 = surf(X, Y, squeeze(SNR_mat_3d(5,:,:)), squeeze(sens_matt2_3d(5,:,:)), 'FaceAlpha',0.5);
s6 = surf(X, Y, squeeze(SNR_mat_3d(6,:,:)), squeeze(sens_matt2_3d(6,:,:)),  'FaceAlpha',0.5);
s7 = surf(X, Y, squeeze(SNR_mat_3d(7,:,:)), squeeze(sens_matt2_3d(7,:,:)), 'FaceAlpha',0.5);
s1.EdgeColor = 'none'; s2.EdgeColor = 'none'; s3.EdgeColor = 'none'; s4.EdgeColor = 'none';
s5.EdgeColor = 'none'; s6.EdgeColor = 'none'; s7.EdgeColor = 'none';


colormap(brewermap([],'*RdYlBu'));
caxis([0.5 1]);
zlim([0 15]);


width_max_ps = [0.4, 3];
res_array_ps = [0 15];
[X,Y,Z] = meshgrid(noise_level, width_max_ps, res_array_ps);
xslice = [5 10];   
yslice = [];
zslice = [];
sens_matt2_3d_ps = zeros(size(X));
ss = slice(X,Y,Z,sens_matt2_3d_ps,xslice,yslice,zslice)
%ss(1).EdgeColor = 'none';
set(ss(1),'FaceColor','interp',...
    'FaceAlpha',0.2, ...
    'LineWidth', 1.5);
%ss(2).EdgeColor = 'none';
set(ss(2),'FaceColor','interp',...
    'FaceAlpha',0.2,...
    'LineWidth', 1.5);

v = [2 -4 2];
[caz,cel] = view(v);

%%
sz = 72;
noise_level = 10;
figure('Position', [0 100 400 400]);
p1 = scatter(width_mat(1,:), SNR_mat_3d(1,:,noise_level), sz, sens_matt2_3d(1,:,noise_level),'filled');
hold on;
p1 = scatter(width_mat(2,:), SNR_mat_3d(2,:,noise_level), sz, sens_matt2_3d(2,:,noise_level),'filled');
p1 = scatter(width_mat(3,:), SNR_mat_3d(3,:,noise_level), sz, sens_matt2_3d(3,:,noise_level),'filled');
p1 = scatter(width_mat(4,:), SNR_mat_3d(4,:,noise_level), sz, sens_matt2_3d(4,:,noise_level),'filled');
p1 = scatter(width_mat(5,:), SNR_mat_3d(5,:,noise_level), sz, sens_matt2_3d(5,:,noise_level),'filled');
p1 = scatter(width_mat(6,:), SNR_mat_3d(6,:,noise_level), sz, sens_matt2_3d(6,:,noise_level),'filled');
p1 = scatter(width_mat(7,:), SNR_mat_3d(7,:,noise_level), sz, sens_matt2_3d(7,:,noise_level),'filled');
caxis([0.5 1]);

[X,Y] = meshgrid(width_mat(1,:),linspace(1,12,7));
[C,h] = contour(X,Y,sens_matt2_3d(:,:,noise_level), 'LineWidth', 2,'ShowText','on');
clabel(C,h,'FontSize',18,'Color','black');
xlim([0 3]); % ylim([0 24]);
% axis off;
colormap(brewermap([],'*RdYlBu'));


sz = 72;
noise_level = 5;
figure('Position', [0 100 400 400]);
p1 = scatter(width_mat(1,:), SNR_mat_3d(1,:,noise_level), sz, sens_matt2_3d(1,:,noise_level),'filled');
hold on;
p1 = scatter(width_mat(2,:), SNR_mat_3d(2,:,noise_level), sz, sens_matt2_3d(2,:,noise_level),'filled');
p1 = scatter(width_mat(3,:), SNR_mat_3d(3,:,noise_level), sz, sens_matt2_3d(3,:,noise_level),'filled');
p1 = scatter(width_mat(4,:), SNR_mat_3d(4,:,noise_level), sz, sens_matt2_3d(4,:,noise_level),'filled');
p1 = scatter(width_mat(5,:), SNR_mat_3d(5,:,noise_level), sz, sens_matt2_3d(5,:,noise_level),'filled');
p1 = scatter(width_mat(6,:), SNR_mat_3d(6,:,noise_level), sz, sens_matt2_3d(6,:,noise_level),'filled');
p1 = scatter(width_mat(7,:), SNR_mat_3d(7,:,noise_level), sz, sens_matt2_3d(7,:,noise_level),'filled');
caxis([0.5 1]);

[X,Y] = meshgrid(width_mat(1,:),linspace(1,24,7));
[C,h] = contour(X,Y,sens_matt2_3d(:,:,noise_level), 'LineWidth', 2, 'ShowText','on');
xlim([0 3]);
% set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
clabel(C,h,'FontSize',18,'Color','black');
% axis off;
colormap(brewermap([],'*RdYlBu'));

%% Linear regression of SNR vs voxel size
color_cell = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]};
color_cell_avg16 = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell_invivo = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};

sz = 36;
SNR_mat_avg = (SNR_mat_04 + SNR_mat_06 + SNR_mat_08 + SNR_mat_10 + SNR_mat_15 + SNR_mat_20 + ...
    SNR_mat_25 + SNR_mat_30) / 8;
SNR_mat_std = std(cat(3, SNR_mat_04, SNR_mat_06, SNR_mat_08, SNR_mat_10, SNR_mat_15, SNR_mat_20, SNR_mat_25, SNR_mat_30), [], 3);
px_sz = res_array.^2;
%px_sz = res_array;
plotHandles = zeros(4,4);

figure('Position', [0 100 300 400]);
% plotHandles(:,1) = plot(px_sz, SNR_mat_avg(:,5), '.', 'MarkerSize', sz, 'Color', [0, 0.4470, 0.7410]);
plotHandles(:,1) = errorbar(px_sz, SNR_mat_avg(:,5), SNR_mat_std(:,5), '-s', 'LineWidth', 2, 'Color', color_cell{1});
hold on;
%plot(px_sz, SNR_mat_avg(:,5), '-', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2);
plotHandles(:,2) = errorbar(px_sz, SNR_mat_avg(:,10), SNR_mat_std(:,10), '-s', 'LineWidth', 2, 'Color', color_cell{2});
%plot(px_sz, SNR_mat_avg(:,10), '-', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2);
%plotHandles(:,3) = plot(px_sz, SNR_mat_avg(:,5), '-', 'MarkerSize', sz, 'Color', [0, 0.4470, 0.7410]);
%plotHandles(:,4) = plot(px_sz, SNR_mat_avg(:,10), '-', 'MarkerSize', sz, 'Color', [0, 0.4470, 0.7410]);
xlim([0 4.1])
grid on;

set(gca, 'XTickLabels', []);
set(gca, 'YTickLabels', []);

set(plotHandles(:,1), 'LineStyle', 'none', 'Marker', '.', 'Color', color_cell_avg16{4});
set(plotHandles(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 12, ...
    'MarkerEdgeColor',  color_cell_avg16{5}, 'MarkerFaceColor', color_cell_avg16{2});
set(plotHandles(:,2), 'LineStyle', 'none', 'Marker', 's', 'Color', color_cell_invivo{4});
set(plotHandles(:,2), 'LineWidth', 1, 'Marker', 's', 'MarkerSize', 12, ...
    'MarkerEdgeColor', color_cell_invivo{5}, 'MarkerFaceColor' , color_cell_invivo{2});

%set(plotHandles(:,3), 'LineWidth', 1, 'Color', color_cell_avg16{4});
%set(plotHandles(:,4), 'LineWidth', 1, 'Color', color_cell_invivo{4});

%legend({'Medium Noise Level', 'High Noise Level'})
%figure();
%plot(px_sz, SNR_mat_avg(:,10), 'o')

%% Linear regression of SNR vs voxel size (Square)
color_cell = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]};
color_cell_avg16 = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell_invivo = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};

sz = 36;
SNR_mat_avg_sq = (SNR_mat_04.^2 + SNR_mat_06.^2 + SNR_mat_08.^2 + SNR_mat_10.^2 + SNR_mat_15.^2 + SNR_mat_20.^2 + ...
    SNR_mat_25.^2 + SNR_mat_30.^2) / 8;
SNR_mat_std_sq = std(cat(3, SNR_mat_04.^2, SNR_mat_06.^2, SNR_mat_08.^2, SNR_mat_10.^2, SNR_mat_15.^2, SNR_mat_20.^2, SNR_mat_25.^2, SNR_mat_30.^2), [], 3);
px_sz = res_array.^2;
%px_sz = res_array;
plotHandles = zeros(4,4);

figure('Position', [0 100 300 400]);
% plotHandles(:,1) = plot(px_sz, SNR_mat_avg(:,5), '.', 'MarkerSize', sz, 'Color', [0, 0.4470, 0.7410]);
plotHandles(:,1) = errorbar(px_sz, SNR_mat_avg_sq(:,5), SNR_mat_std_sq(:,5), '-s', 'LineWidth', 2, 'Color', color_cell{1});
hold on;
%plot(px_sz, SNR_mat_avg(:,5), '-', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2);
plotHandles(:,2) = errorbar(px_sz, SNR_mat_avg_sq(:,10), SNR_mat_std_sq(:,10), '-s', 'LineWidth', 2, 'Color', color_cell{2});
%plot(px_sz, SNR_mat_avg(:,10), '-', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2);
%plotHandles(:,3) = errorbar(px_sz, SNR_mat_avg_sq(:,1), SNR_mat_std_sq(:,1), '-s', 'LineWidth', 2, 'Color', color_cell{1});
%plotHandles(:,4) = plot(px_sz, SNR_mat_avg(:,10), '-', 'MarkerSize', sz, 'Color', [0, 0.4470, 0.7410]);
xlim([0 4.1])
grid on;

set(gca, 'XTickLabels', []);
set(gca, 'YTickLabels', []);

set(plotHandles(:,1), 'LineStyle', 'none', 'Marker', '.', 'Color', color_cell_avg16{4});
set(plotHandles(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 12, ...
    'MarkerEdgeColor',  color_cell_avg16{5}, 'MarkerFaceColor', color_cell_avg16{2});
set(plotHandles(:,2), 'LineStyle', 'none', 'Marker', 's', 'Color', color_cell_invivo{4});
set(plotHandles(:,2), 'LineWidth', 1, 'Marker', 's', 'MarkerSize', 12, ...
    'MarkerEdgeColor', color_cell_invivo{5}, 'MarkerFaceColor' , color_cell_invivo{2});

%set(plotHandles(:,3), 'LineWidth', 1, 'Color', color_cell_avg16{4});
%set(plotHandles(:,4), 'LineWidth', 1, 'Color', color_cell_invivo{4});

%legend({'Medium Noise Level', 'High Noise Level'})
%figure();
%plot(px_sz, SNR_mat_avg(:,10), 'o')

y = SNR_mat_avg_sq(:,5);
X = [ones(length(px_sz),1) px_sz.'];
b = X\SNR_mat_avg_sq(:,5);
yCalc2 = X*b;
hold on;
plot(px_sz,yCalc2,'--', 'LineWidth', 1.5, 'Color', color_cell{2})
Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2)

y = SNR_mat_avg_sq(:,10);
X = [ones(length(px_sz),1) px_sz.'];
b = X\SNR_mat_avg_sq(:,10);
yCalc2 = X*b;
plot(px_sz,yCalc2,'--', 'LineWidth', 1.5, 'Color', color_cell{1})
Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2)

%% Linear regression of SNR vs pixel size (Square root of pixel size)
color_cell = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]};
color_cell_avg16 = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell_invivo = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};

sz = 36;
SNR_mat_avg_sq = (SNR_mat_04 + SNR_mat_06 + SNR_mat_08 + SNR_mat_10 + SNR_mat_15 + SNR_mat_20 + ...
    SNR_mat_25 + SNR_mat_30) / 8;
SNR_mat_std_sq = std(cat(3, SNR_mat_04, SNR_mat_06, SNR_mat_08, SNR_mat_10, SNR_mat_15, SNR_mat_20, SNR_mat_25, SNR_mat_30), [], 3);
px_sz = log(res_array.^2);
%px_sz = res_array;
plotHandles = zeros(4,4);

figure('Position', [0 100 300 400]);
% plotHandles(:,1) = plot(px_sz, SNR_mat_avg(:,5), '.', 'MarkerSize', sz, 'Color', [0, 0.4470, 0.7410]);
plotHandles(:,1) = errorbar(px_sz, SNR_mat_avg_sq(:,5), SNR_mat_std_sq(:,5), '-s', 'LineWidth', 2, 'Color', color_cell{1});
hold on;
%plot(px_sz, SNR_mat_avg(:,5), '-', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2);
plotHandles(:,2) = errorbar(px_sz, SNR_mat_avg_sq(:,10), SNR_mat_std_sq(:,10), '-s', 'LineWidth', 2, 'Color', color_cell{2});
%plot(px_sz, SNR_mat_avg(:,10), '-', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2);
%plotHandles(:,3) = errorbar(px_sz, SNR_mat_avg_sq(:,1), SNR_mat_std_sq(:,1), '-s', 'LineWidth', 2, 'Color', color_cell{1});
%plotHandles(:,4) = plot(px_sz, SNR_mat_avg(:,10), '-', 'MarkerSize', sz, 'Color', [0, 0.4470, 0.7410]);
xlim([0 2.1])
grid on;

set(gca, 'XTickLabels', []);
set(gca, 'YTickLabels', []);

set(plotHandles(:,1), 'LineStyle', 'none', 'Marker', '.', 'Color', color_cell_avg16{4});
set(plotHandles(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 12, ...
    'MarkerEdgeColor',  color_cell_avg16{5}, 'MarkerFaceColor', color_cell_avg16{2});
set(plotHandles(:,2), 'LineStyle', 'none', 'Marker', 's', 'Color', color_cell_invivo{4});
set(plotHandles(:,2), 'LineWidth', 1, 'Marker', 's', 'MarkerSize', 12, ...
    'MarkerEdgeColor', color_cell_invivo{5}, 'MarkerFaceColor' , color_cell_invivo{2});

%set(plotHandles(:,3), 'LineWidth', 1, 'Color', color_cell_avg16{4});
%set(plotHandles(:,4), 'LineWidth', 1, 'Color', color_cell_invivo{4});

%legend({'Medium Noise Level', 'High Noise Level'})
%figure();
%plot(px_sz, SNR_mat_avg(:,10), 'o')

y = SNR_mat_avg_sq(:,5);
X = [ones(length(px_sz),1) px_sz.'];
b = X\SNR_mat_avg_sq(:,5);
yCalc2 = X*b;
hold on;
plot(px_sz,yCalc2,'--', 'LineWidth', 1.5, 'Color', color_cell{2})
Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2)

y = SNR_mat_avg_sq(:,10);
X = [ones(length(px_sz),1) px_sz.'];
b = X\SNR_mat_avg_sq(:,10);
yCalc2 = X*b;
plot(px_sz,yCalc2,'--', 'LineWidth', 1.5, 'Color', color_cell{1})
Rsq2 = 1 - sum((y - yCalc2).^2)/sum((y - mean(y)).^2)
%% Linear regression of SD vs voxel size
%% Linear regression of SD vs voxel size
color_cell = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]};
color_cell_avg16 = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell_invivo = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};

sz = 36;
SD_mat_avg = (SD_mat_04 + SD_mat_06 + SD_mat_08 + SD_mat_10 + SD_mat_15 + SD_mat_20 + ...
    SD_mat_25 + SD_mat_30) / 8;
SD_mat_std = std(cat(3, SD_mat_04, SD_mat_06, SD_mat_08, SD_mat_10, SD_mat_15, SD_mat_20, SD_mat_25, SD_mat_30), [], 3);
px_sz = res_array.^2;
%px_sz = res_array;
plotHandles = zeros(4,4);

figure('Position', [0 100 300 400]);
% plotHandles(:,1) = plot(px_sz, SNR_mat_avg(:,5), '.', 'MarkerSize', sz, 'Color', [0, 0.4470, 0.7410]);
plotHandles(:,1) = errorbar(px_sz, SD_mat_avg(:,5), SD_mat_std(:,5), '-s', 'LineWidth', 2, 'Color', color_cell{1});
hold on;
%plot(px_sz, SNR_mat_avg(:,5), '-', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2);
plotHandles(:,2) = errorbar(px_sz, SD_mat_avg(:,10), SD_mat_std(:,10), '-s', 'LineWidth', 2, 'Color', color_cell{2});
%plot(px_sz, SNR_mat_avg(:,10), '-', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2);
%plotHandles(:,3) = plot(px_sz, SNR_mat_avg(:,5), '-', 'MarkerSize', sz, 'Color', [0, 0.4470, 0.7410]);
%plotHandles(:,4) = plot(px_sz, SNR_mat_avg(:,10), '-', 'MarkerSize', sz, 'Color', [0, 0.4470, 0.7410]);
xlim([0 4.1])
grid on;

set(gca, 'XTickLabels', []);
set(gca, 'YTickLabels', []);

set(plotHandles(:,1), 'LineStyle', 'none', 'Marker', '.', 'Color', color_cell_avg16{4});
set(plotHandles(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 12, ...
    'MarkerEdgeColor',  color_cell_avg16{5}, 'MarkerFaceColor', color_cell_avg16{2});
set(plotHandles(:,2), 'LineStyle', 'none', 'Marker', 's', 'Color', color_cell_invivo{4});
set(plotHandles(:,2), 'LineWidth', 1, 'Marker', 's', 'MarkerSize', 12, ...
    'MarkerEdgeColor', color_cell_invivo{5}, 'MarkerFaceColor' , color_cell_invivo{2});

%set(plotHandles(:,3), 'LineWidth', 1, 'Color', color_cell_avg16{4});
%set(plotHandles(:,4), 'LineWidth', 1, 'Color', color_cell_invivo{4});

%legend({'Medium Noise Level', 'High Noise Level'})
%figure();
%plot(px_sz, SNR_mat_avg(:,10), 'o')

%% Linear regression of SD vs voxel size 
color_cell = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]};
color_cell_avg16 = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell_invivo = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};

p = -2;
sz = 36;
SD_mat_avg = (SD_mat_04.^p + SD_mat_06.^p + SD_mat_08.^p + SD_mat_10.^p + SD_mat_15.^p + SD_mat_20.^p + ...
    SD_mat_25.^p + SD_mat_30.^p) / 8;
SD_mat_std = std(cat(3, SD_mat_04.^p, SD_mat_06.^p, SD_mat_08.^p, SD_mat_10.^p, SD_mat_15.^p, SD_mat_20.^p, SD_mat_25.^p, SD_mat_30.^p), [], 3);
px_sz = res_array.^2;
%px_sz = res_array;
plotHandles = zeros(4,4);

figure('Position', [0 100 300 400]);
% plotHandles(:,1) = plot(px_sz, SNR_mat_avg(:,5), '.', 'MarkerSize', sz, 'Color', [0, 0.4470, 0.7410]);
plotHandles(:,1) = errorbar(px_sz, SD_mat_avg(:,5), SD_mat_std(:,5), '-s', 'LineWidth', 2, 'Color', color_cell{1});
hold on;
%plot(px_sz, SNR_mat_avg(:,5), '-', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2);
plotHandles(:,2) = errorbar(px_sz, SD_mat_avg(:,10), SD_mat_std(:,10), '-s', 'LineWidth', 2, 'Color', color_cell{2});
%plot(px_sz, SNR_mat_avg(:,10), '-', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2);
%plotHandles(:,3) = plot(px_sz, SNR_mat_avg(:,5), '-', 'MarkerSize', sz, 'Color', [0, 0.4470, 0.7410]);
%plotHandles(:,4) = plot(px_sz, SNR_mat_avg(:,10), '-', 'MarkerSize', sz, 'Color', [0, 0.4470, 0.7410]);
xlim([0 4.1])
grid on;

set(gca, 'XTickLabels', []);
set(gca, 'YTickLabels', []);

set(plotHandles(:,1), 'LineStyle', 'none', 'Marker', '.', 'Color', color_cell_avg16{4});
set(plotHandles(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 12, ...
    'MarkerEdgeColor',  color_cell_avg16{5}, 'MarkerFaceColor', color_cell_avg16{2});
set(plotHandles(:,2), 'LineStyle', 'none', 'Marker', 's', 'Color', color_cell_invivo{4});
set(plotHandles(:,2), 'LineWidth', 1, 'Marker', 's', 'MarkerSize', 12, ...
    'MarkerEdgeColor', color_cell_invivo{5}, 'MarkerFaceColor' , color_cell_invivo{2});


%% Balanced Accuracy (3D surface plot)
noise_level = 1:1:10;
C = zeros([size(sens_matt2_3d,2) size(sens_matt2_3d,3)]);
[X,Y] = meshgrid(noise_level, width_max);
% noise_level = 5;
figure('Position', [0 100 800 400]);
s1 = surf(X, Y, squeeze(SNR_mat_3d(1,:,:)), squeeze(accu_matt3_3d(1,:,:)), 'FaceAlpha',0.5);
hold on;
s2 = surf(X, Y, squeeze(SNR_mat_3d(2,:,:)), squeeze(accu_matt3_3d(2,:,:)), 'FaceAlpha',0.5);
s3 = surf(X, Y, squeeze(SNR_mat_3d(3,:,:)), squeeze(accu_matt3_3d(3,:,:)), 'FaceAlpha',0.5);
s4 = surf(X, Y, squeeze(SNR_mat_3d(4,:,:)), squeeze(accu_matt3_3d(4,:,:)), 'FaceAlpha',0.5);
s5 = surf(X, Y, squeeze(SNR_mat_3d(5,:,:)), squeeze(accu_matt3_3d(5,:,:)), 'FaceAlpha',0.5);
s6 = surf(X, Y, squeeze(SNR_mat_3d(6,:,:)), squeeze(accu_matt3_3d(6,:,:)),  'FaceAlpha',0.5);
s7 = surf(X, Y, squeeze(SNR_mat_3d(7,:,:)), squeeze(accu_matt3_3d(7,:,:)), 'FaceAlpha',0.5);
s1.EdgeColor = 'none'; s2.EdgeColor = 'none'; s3.EdgeColor = 'none'; s4.EdgeColor = 'none';
s5.EdgeColor = 'none'; s6.EdgeColor = 'none'; s7.EdgeColor = 'none';


colormap(brewermap([],'*RdYlBu'));
caxis([0.5 1]);
zlim([0 15]);


width_max_ps = [0.4, 3];
res_array_ps = [0 15];
[X,Y,Z] = meshgrid(noise_level, width_max_ps, res_array_ps);
xslice = [5 10];   
yslice = [];
zslice = [];
sens_matt2_3d_ps = zeros(size(X));
ss = slice(X,Y,Z,sens_matt2_3d_ps,xslice,yslice,zslice)
%ss(1).EdgeColor = 'none';
set(ss(1),'FaceColor','interp',...
    'FaceAlpha',0.2, ...
    'LineWidth', 1.5);
%ss(2).EdgeColor = 'none';
set(ss(2),'FaceColor','interp',...
    'FaceAlpha',0.2,...
    'LineWidth', 1.5);

v = [2 -4 2];
[caz,cel] = view(v);

%% Balanced Accuracy (Dots with contours)
sz = 72;
noise_level = 10;
figure('Position', [0 100 400 400]);
p1 = scatter(width_mat(1,:), SNR_mat_3d(1,:,noise_level), sz, accu_matt3_3d(1,:,noise_level),'filled');
hold on;
p1 = scatter(width_mat(2,:), SNR_mat_3d(2,:,noise_level), sz, accu_matt3_3d(2,:,noise_level),'filled');
p1 = scatter(width_mat(3,:), SNR_mat_3d(3,:,noise_level), sz, accu_matt3_3d(3,:,noise_level),'filled');
p1 = scatter(width_mat(4,:), SNR_mat_3d(4,:,noise_level), sz, accu_matt3_3d(4,:,noise_level),'filled');
p1 = scatter(width_mat(5,:), SNR_mat_3d(5,:,noise_level), sz, accu_matt3_3d(5,:,noise_level),'filled');
p1 = scatter(width_mat(6,:), SNR_mat_3d(6,:,noise_level), sz, accu_matt3_3d(6,:,noise_level),'filled');
p1 = scatter(width_mat(7,:), SNR_mat_3d(7,:,noise_level), sz, accu_matt3_3d(7,:,noise_level),'filled');
caxis([0 1]);

[X,Y] = meshgrid(width_mat(1,:),linspace(1,12,7));
[C,h] = contour(X,Y,accu_matt3_3d(:,:,noise_level), 'LineWidth', 2,'ShowText','on');
clabel(C,h,'FontSize',18,'Color','black');
xlim([0 3]); % ylim([0 24]);
% axis off;
colormap(brewermap([],'*RdYlBu'));


sz = 72;
noise_level = 5;
figure('Position', [0 100 400 400]);
p1 = scatter(width_mat(1,:), SNR_mat_3d(1,:,noise_level), sz, accu_matt3_3d(1,:,noise_level),'filled');
hold on;
p1 = scatter(width_mat(2,:), SNR_mat_3d(2,:,noise_level), sz, accu_matt3_3d(2,:,noise_level),'filled');
p1 = scatter(width_mat(3,:), SNR_mat_3d(3,:,noise_level), sz, accu_matt3_3d(3,:,noise_level),'filled');
p1 = scatter(width_mat(4,:), SNR_mat_3d(4,:,noise_level), sz, accu_matt3_3d(4,:,noise_level),'filled');
p1 = scatter(width_mat(5,:), SNR_mat_3d(5,:,noise_level), sz, accu_matt3_3d(5,:,noise_level),'filled');
p1 = scatter(width_mat(6,:), SNR_mat_3d(6,:,noise_level), sz, accu_matt3_3d(6,:,noise_level),'filled');
p1 = scatter(width_mat(7,:), SNR_mat_3d(7,:,noise_level), sz, accu_matt3_3d(7,:,noise_level),'filled');
caxis([0 1]);

[X,Y] = meshgrid(width_mat(1,:),linspace(1,24,7));
[C,h] = contour(X,Y,accu_matt3_3d(:,:,noise_level), 'LineWidth', 2, 'ShowText','on');
xlim([0 3]);
% set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
clabel(C,h,'FontSize',18,'Color','black');
% axis off;
colormap(brewermap([],'*RdYlBu'));

%% Version 2
sz = 72;
noise_level = 10;
figure('Position', [0 100 400 400]);
p1 = scatter(width_mat(1,:), repmat(res_array(1), [1, size(width_mat, 2)]), sz, sens_matt2_3d(1,:,noise_level),'filled');
hold on;
p1 = scatter(width_mat(2,:), repmat(res_array(2), [1, size(width_mat, 2)]), sz, sens_matt2_3d(2,:,noise_level),'filled');
p1 = scatter(width_mat(3,:), repmat(res_array(3), [1, size(width_mat, 2)]), sz, sens_matt2_3d(3,:,noise_level),'filled');
p1 = scatter(width_mat(4,:), repmat(res_array(4), [1, size(width_mat, 2)]), sz, sens_matt2_3d(4,:,noise_level),'filled');
p1 = scatter(width_mat(5,:), repmat(res_array(5), [1, size(width_mat, 2)]), sz, sens_matt2_3d(5,:,noise_level),'filled');
p1 = scatter(width_mat(6,:), repmat(res_array(6), [1, size(width_mat, 2)]), sz, sens_matt2_3d(6,:,noise_level),'filled');
p1 = scatter(width_mat(7,:), repmat(res_array(7), [1, size(width_mat, 2)]), sz, sens_matt2_3d(7,:,noise_level),'filled');
caxis([0.5 1]);

[X,Y] = meshgrid(width_mat(1,:),res_array);
[C,h] = contour(X,Y,sens_matt2_3d(:,:,noise_level), 'LineWidth', 2,'ShowText','on');
clabel(C,h,'FontSize',18,'Color','black');
xlim([0 3]); ylim([0.3 2]);
% axis off;
colormap(brewermap([],'*RdYlBu'));


sz = 72;
noise_level = 5;
figure('Position', [0 100 400 400]);
p1 = scatter(width_mat(1,:), repmat(res_array(1), [1, size(width_mat, 2)]), sz, sens_matt2_3d(1,:,noise_level),'filled');
hold on;
p1 = scatter(width_mat(2,:), repmat(res_array(2), [1, size(width_mat, 2)]), sz, sens_matt2_3d(2,:,noise_level),'filled');
p1 = scatter(width_mat(3,:), repmat(res_array(3), [1, size(width_mat, 2)]), sz, sens_matt2_3d(3,:,noise_level),'filled');
p1 = scatter(width_mat(4,:), repmat(res_array(4), [1, size(width_mat, 2)]), sz, sens_matt2_3d(4,:,noise_level),'filled');
p1 = scatter(width_mat(5,:), repmat(res_array(5), [1, size(width_mat, 2)]), sz, sens_matt2_3d(5,:,noise_level),'filled');
p1 = scatter(width_mat(6,:), repmat(res_array(6), [1, size(width_mat, 2)]), sz, sens_matt2_3d(6,:,noise_level),'filled');
p1 = scatter(width_mat(7,:), repmat(res_array(7), [1, size(width_mat, 2)]), sz, sens_matt2_3d(7,:,noise_level),'filled');
caxis([0.5 1]);

[X,Y] = meshgrid(width_mat(1,:),res_array);
[C,h] = contour(X,Y,sens_matt2_3d(:,:,noise_level), 'LineWidth', 2, 'ShowText','on');
xlim([0 3]); ylim([0.3 2]);
% set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
clabel(C,h,'FontSize',18,'Color','black');
% axis off;
colormap(brewermap([],'*RdYlBu'));

%% Version 2 (Balanced Accuracy)
sz = 72;
noise_level = 10;
figure('Position', [0 100 400 400]);
p1 = scatter(width_mat(1,:), repmat(res_array(1), [1, size(width_mat, 2)]), sz, accu_matt3_3d(1,:,noise_level),'filled');
hold on;
p1 = scatter(width_mat(2,:), repmat(res_array(2), [1, size(width_mat, 2)]), sz, accu_matt3_3d(2,:,noise_level),'filled');
p1 = scatter(width_mat(3,:), repmat(res_array(3), [1, size(width_mat, 2)]), sz, accu_matt3_3d(3,:,noise_level),'filled');
p1 = scatter(width_mat(4,:), repmat(res_array(4), [1, size(width_mat, 2)]), sz, accu_matt3_3d(4,:,noise_level),'filled');
p1 = scatter(width_mat(5,:), repmat(res_array(5), [1, size(width_mat, 2)]), sz, accu_matt3_3d(5,:,noise_level),'filled');
p1 = scatter(width_mat(6,:), repmat(res_array(6), [1, size(width_mat, 2)]), sz, accu_matt3_3d(6,:,noise_level),'filled');
p1 = scatter(width_mat(7,:), repmat(res_array(7), [1, size(width_mat, 2)]), sz, accu_matt3_3d(7,:,noise_level),'filled');
caxis([0.5 1]);

[X,Y] = meshgrid(width_mat(1,:),res_array);
[C,h] = contour(X,Y,accu_matt3_3d(:,:,noise_level), 'LineWidth', 2,'ShowText','on');
clabel(C,h,'FontSize',18,'Color','black');
xlim([0 3]); ylim([0.3 2]);
% axis off;
colormap(brewermap([],'*RdYlBu'));


sz = 72;
noise_level = 5;
figure('Position', [0 100 400 400]);
p1 = scatter(width_mat(1,:), repmat(res_array(1), [1, size(width_mat, 2)]), sz, accu_matt3_3d(1,:,noise_level),'filled');
hold on;
p1 = scatter(width_mat(2,:), repmat(res_array(2), [1, size(width_mat, 2)]), sz, accu_matt3_3d(2,:,noise_level),'filled');
p1 = scatter(width_mat(3,:), repmat(res_array(3), [1, size(width_mat, 2)]), sz, accu_matt3_3d(3,:,noise_level),'filled');
p1 = scatter(width_mat(4,:), repmat(res_array(4), [1, size(width_mat, 2)]), sz, accu_matt3_3d(4,:,noise_level),'filled');
p1 = scatter(width_mat(5,:), repmat(res_array(5), [1, size(width_mat, 2)]), sz, accu_matt3_3d(5,:,noise_level),'filled');
p1 = scatter(width_mat(6,:), repmat(res_array(6), [1, size(width_mat, 2)]), sz, accu_matt3_3d(6,:,noise_level),'filled');
p1 = scatter(width_mat(7,:), repmat(res_array(7), [1, size(width_mat, 2)]), sz, accu_matt3_3d(7,:,noise_level),'filled');
caxis([0.5 1]);

[X,Y] = meshgrid(width_mat(1,:),res_array);
[C,h] = contour(X,Y,accu_matt3_3d(:,:,noise_level), 'LineWidth', 2, 'ShowText','on');
xlim([0 3]); ylim([0.3 2]);
% set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
clabel(C,h,'FontSize',18,'Color','black');
% axis off;
colormap(brewermap([],'*RdYlBu'));

%% Balanced AUC (3D surface plot)
noise_level = 1:1:10;
C = zeros([size(sens_matt2_3d,2) size(sens_matt2_3d,3)]);
[X,Y] = meshgrid(noise_level, width_max);
% noise_level = 5;
figure('Position', [0 100 800 400]);
s1 = surf(X, Y, squeeze(SNR_mat_3d(1,:,:)), squeeze(auc_matt3_3d(1,:,:)), 'FaceAlpha',0.5);
hold on;
s2 = surf(X, Y, squeeze(SNR_mat_3d(2,:,:)), squeeze(auc_matt3_3d(2,:,:)), 'FaceAlpha',0.5);
s3 = surf(X, Y, squeeze(SNR_mat_3d(3,:,:)), squeeze(auc_matt3_3d(3,:,:)), 'FaceAlpha',0.5);
s4 = surf(X, Y, squeeze(SNR_mat_3d(4,:,:)), squeeze(auc_matt3_3d(4,:,:)), 'FaceAlpha',0.5);
s5 = surf(X, Y, squeeze(SNR_mat_3d(5,:,:)), squeeze(auc_matt3_3d(5,:,:)), 'FaceAlpha',0.5);
s6 = surf(X, Y, squeeze(SNR_mat_3d(6,:,:)), squeeze(auc_matt3_3d(6,:,:)),  'FaceAlpha',0.5);
s7 = surf(X, Y, squeeze(SNR_mat_3d(7,:,:)), squeeze(auc_matt3_3d(7,:,:)), 'FaceAlpha',0.5);
s1.EdgeColor = 'none'; s2.EdgeColor = 'none'; s3.EdgeColor = 'none'; s4.EdgeColor = 'none';
s5.EdgeColor = 'none'; s6.EdgeColor = 'none'; s7.EdgeColor = 'none';


colormap(brewermap([],'*RdYlBu'));
caxis([0.5 1]);
zlim([0 15]);


width_max_ps = [0.4, 3];
res_array_ps = [0 15];
[X,Y,Z] = meshgrid(noise_level, width_max_ps, res_array_ps);
xslice = [5 10];   
yslice = [];
zslice = [];
auc_matt3_3d_ps = zeros(size(X));
ss = slice(X,Y,Z,auc_matt3_3d_ps,xslice,yslice,zslice)
%ss(1).EdgeColor = 'none';
set(ss(1),'FaceColor','interp',...
    'FaceAlpha',0.2, ...
    'LineWidth', 1.5);
%ss(2).EdgeColor = 'none';
set(ss(2),'FaceColor','interp',...
    'FaceAlpha',0.2,...
    'LineWidth', 1.5);

v = [2 -4 2];
[caz,cel] = view(v);
%% Version 2 (AUC)
sz = 72;
noise_level = 10;
figure('Position', [0 100 400 400]);
p1 = scatter(width_mat(1,:), repmat(res_array(1), [1, size(width_mat, 2)]), sz, auc_matt3_3d(1,:,noise_level),'filled');
hold on;
p1 = scatter(width_mat(2,:), repmat(res_array(2), [1, size(width_mat, 2)]), sz, auc_matt3_3d(2,:,noise_level),'filled');
p1 = scatter(width_mat(3,:), repmat(res_array(3), [1, size(width_mat, 2)]), sz, auc_matt3_3d(3,:,noise_level),'filled');
p1 = scatter(width_mat(4,:), repmat(res_array(4), [1, size(width_mat, 2)]), sz, auc_matt3_3d(4,:,noise_level),'filled');
p1 = scatter(width_mat(5,:), repmat(res_array(5), [1, size(width_mat, 2)]), sz, auc_matt3_3d(5,:,noise_level),'filled');
p1 = scatter(width_mat(6,:), repmat(res_array(6), [1, size(width_mat, 2)]), sz, auc_matt3_3d(6,:,noise_level),'filled');
p1 = scatter(width_mat(7,:), repmat(res_array(7), [1, size(width_mat, 2)]), sz, auc_matt3_3d(7,:,noise_level),'filled');
caxis([0.5 1]);

[X,Y] = meshgrid(width_mat(1,:),res_array);
[C,h] = contour(X,Y,auc_matt3_3d(:,:,noise_level), 'LineWidth', 2,'ShowText','on');
clabel(C,h,'FontSize',18,'Color','black');
xlim([0 3]); ylim([0.3 2]);
% axis off;
colormap(brewermap([],'*RdYlBu'));


sz = 72;
noise_level = 5;
figure('Position', [0 100 400 400]);
p1 = scatter(width_mat(1,:), repmat(res_array(1), [1, size(width_mat, 2)]), sz, auc_matt3_3d(1,:,noise_level),'filled');
hold on;
p1 = scatter(width_mat(2,:), repmat(res_array(2), [1, size(width_mat, 2)]), sz, auc_matt3_3d(2,:,noise_level),'filled');
p1 = scatter(width_mat(3,:), repmat(res_array(3), [1, size(width_mat, 2)]), sz, auc_matt3_3d(3,:,noise_level),'filled');
p1 = scatter(width_mat(4,:), repmat(res_array(4), [1, size(width_mat, 2)]), sz, auc_matt3_3d(4,:,noise_level),'filled');
p1 = scatter(width_mat(5,:), repmat(res_array(5), [1, size(width_mat, 2)]), sz, auc_matt3_3d(5,:,noise_level),'filled');
p1 = scatter(width_mat(6,:), repmat(res_array(6), [1, size(width_mat, 2)]), sz, auc_matt3_3d(6,:,noise_level),'filled');
p1 = scatter(width_mat(7,:), repmat(res_array(7), [1, size(width_mat, 2)]), sz, auc_matt3_3d(7,:,noise_level),'filled');
caxis([0.5 1]);

[X,Y] = meshgrid(width_mat(1,:),res_array);
[C,h] = contour(X,Y,auc_matt3_3d(:,:,noise_level), 'LineWidth', 2, 'ShowText','on');
xlim([0 3]); ylim([0.3 2]);
% set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
clabel(C,h,'FontSize',18,'Color','black');
% axis off;
colormap(brewermap([],'*RdYlBu'));

sz = 72;
noise_level = 1;
figure('Position', [0 100 400 400]);
p1 = scatter(width_mat(1,:), repmat(res_array(1), [1, size(width_mat, 2)]), sz, auc_matt3_3d(1,:,noise_level),'filled');
hold on;
p1 = scatter(width_mat(2,:), repmat(res_array(2), [1, size(width_mat, 2)]), sz, auc_matt3_3d(2,:,noise_level),'filled');
p1 = scatter(width_mat(3,:), repmat(res_array(3), [1, size(width_mat, 2)]), sz, auc_matt3_3d(3,:,noise_level),'filled');
p1 = scatter(width_mat(4,:), repmat(res_array(4), [1, size(width_mat, 2)]), sz, auc_matt3_3d(4,:,noise_level),'filled');
p1 = scatter(width_mat(5,:), repmat(res_array(5), [1, size(width_mat, 2)]), sz, auc_matt3_3d(5,:,noise_level),'filled');
p1 = scatter(width_mat(6,:), repmat(res_array(6), [1, size(width_mat, 2)]), sz, auc_matt3_3d(6,:,noise_level),'filled');
p1 = scatter(width_mat(7,:), repmat(res_array(7), [1, size(width_mat, 2)]), sz, auc_matt3_3d(7,:,noise_level),'filled');
caxis([0.5 1]);

[X,Y] = meshgrid(width_mat(1,:),res_array);
[C,h] = contour(X,Y,auc_matt3_3d(:,:,noise_level), [0.88, 0.90, 0.92, 0.94, 0.96, 0.98], 'LineWidth', 2, 'ShowText','on');
xlim([0 3]); ylim([0.3 2]);
% set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
clabel(C,h,'FontSize',18,'Color','black');
% axis off;
colormap(brewermap([],'*RdYlBu'));

%% 2024 AUC for PPT
sz = 72;
noise_level = 10;
figure('Position', [0 100 400 400]);
p1 = scatter(width_mat(1,:), repmat(res_array(1), [1, size(width_mat, 2)]), sz, auc_matt3_3d(1,:,noise_level),'filled');
hold on;
p1 = scatter(width_mat(2,:), repmat(res_array(2), [1, size(width_mat, 2)]), sz, auc_matt3_3d(2,:,noise_level),'filled');
p1 = scatter(width_mat(3,:), repmat(res_array(3), [1, size(width_mat, 2)]), sz, auc_matt3_3d(3,:,noise_level),'filled');
p1 = scatter(width_mat(4,:), repmat(res_array(4), [1, size(width_mat, 2)]), sz, auc_matt3_3d(4,:,noise_level),'filled');
p1 = scatter(width_mat(5,:), repmat(res_array(5), [1, size(width_mat, 2)]), sz, auc_matt3_3d(5,:,noise_level),'filled');
p1 = scatter(width_mat(6,:), repmat(res_array(6), [1, size(width_mat, 2)]), sz, auc_matt3_3d(6,:,noise_level),'filled');
p1 = scatter(width_mat(7,:), repmat(res_array(7), [1, size(width_mat, 2)]), sz, auc_matt3_3d(7,:,noise_level),'filled');
caxis([0.5 1]);


[max_lvl3 idx_lvl3] = max(auc_matt2_3d(:,:,noise_level));
idx_lvl3(3) = 7;

hold on; 
pline = plot(width_mat(1,:), res_array(idx_lvl3), 'k', 'LineWidth', 2);

%[X,Y] = meshgrid(width_mat(1,:),res_array);
%[C,h] = contour(X,Y,auc_matt3_3d(:,:,noise_level), 'LineWidth', 2,'ShowText','on');
%clabel(C,h,'FontSize',18,'Color','black');
xlim([0 3]); ylim([0.3 2]);
% axis off;
colormap(brewermap([],'*RdYlBu'));


sz = 72;
noise_level = 5;
figure('Position', [0 100 400 400]);
p1 = scatter(width_mat(1,:), repmat(res_array(1), [1, size(width_mat, 2)]), sz, auc_matt3_3d(1,:,noise_level),'filled');
hold on;
p1 = scatter(width_mat(2,:), repmat(res_array(2), [1, size(width_mat, 2)]), sz, auc_matt3_3d(2,:,noise_level),'filled');
p1 = scatter(width_mat(3,:), repmat(res_array(3), [1, size(width_mat, 2)]), sz, auc_matt3_3d(3,:,noise_level),'filled');
p1 = scatter(width_mat(4,:), repmat(res_array(4), [1, size(width_mat, 2)]), sz, auc_matt3_3d(4,:,noise_level),'filled');
p1 = scatter(width_mat(5,:), repmat(res_array(5), [1, size(width_mat, 2)]), sz, auc_matt3_3d(5,:,noise_level),'filled');
p1 = scatter(width_mat(6,:), repmat(res_array(6), [1, size(width_mat, 2)]), sz, auc_matt3_3d(6,:,noise_level),'filled');
p1 = scatter(width_mat(7,:), repmat(res_array(7), [1, size(width_mat, 2)]), sz, auc_matt3_3d(7,:,noise_level),'filled');
caxis([0.5 1]);


[max_lvl2 idx_lvl2] = max(auc_matt2_3d(:,:,noise_level));
hold on; 
pline = plot(width_mat(1,:), res_array(idx_lvl2),  'k', 'LineWidth', 2);

%[X,Y] = meshgrid(width_mat(1,:),res_array);
%[C,h] = contour(X,Y,auc_matt3_3d(:,:,noise_level), 'LineWidth', 2, 'ShowText','on');
xlim([0 3]); ylim([0.3 2]);
% set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
%clabel(C,h,'FontSize',18,'Color','black');
% axis off;
colormap(brewermap([],'*RdYlBu'));

sz = 72;
noise_level = 1;
figure('Position', [0 100 400 400]);
p1 = scatter(width_mat(1,:), repmat(res_array(1), [1, size(width_mat, 2)]), sz, auc_matt3_3d(1,:,noise_level),'filled');
hold on;
p1 = scatter(width_mat(2,:), repmat(res_array(2), [1, size(width_mat, 2)]), sz, auc_matt3_3d(2,:,noise_level),'filled');
p1 = scatter(width_mat(3,:), repmat(res_array(3), [1, size(width_mat, 2)]), sz, auc_matt3_3d(3,:,noise_level),'filled');
p1 = scatter(width_mat(4,:), repmat(res_array(4), [1, size(width_mat, 2)]), sz, auc_matt3_3d(4,:,noise_level),'filled');
p1 = scatter(width_mat(5,:), repmat(res_array(5), [1, size(width_mat, 2)]), sz, auc_matt3_3d(5,:,noise_level),'filled');
p1 = scatter(width_mat(6,:), repmat(res_array(6), [1, size(width_mat, 2)]), sz, auc_matt3_3d(6,:,noise_level),'filled');
p1 = scatter(width_mat(7,:), repmat(res_array(7), [1, size(width_mat, 2)]), sz, auc_matt3_3d(7,:,noise_level),'filled');
caxis([0.5 1]);


[max_lvl2 idx_lvl2] = max(auc_matt2_3d(:,:,noise_level));
hold on; 
pline = plot(width_mat(1,:), res_array(idx_lvl2),  'k', 'LineWidth', 2);

%[X,Y] = meshgrid(width_mat(1,:),res_array);
%[C,h] = contour(X,Y,auc_matt3_3d(:,:,noise_level), 'LineWidth', 2, 'ShowText','on');
xlim([0 3]); ylim([0.3 2]);
% set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
%clabel(C,h,'FontSize',18,'Color','black');
% axis off;
colormap(brewermap([],'*RdYlBu'));
%%
sz = 72;
noise_level = 5;

[max_lvl2 idx_lvl2] = max(auc_matt2_3d(:,:,noise_level));


noise_level = 10;
[max_lvl3 idx_lvl3] = max(auc_matt2_3d(:,:,noise_level));


figure(); 
p1 = scatter(width_mat(1,:), max_lvl2, sz,'filled');
hold on;

p1 = plot(width_mat(1,:), max_lvl2);

p2 = scatter(width_mat(1,:), max_lvl3, sz,'filled');
p2 = plot(width_mat(1,:), max_lvl3);


%% Voxel vs auc
noise_level = 10;
figure('Position', [0 100 400 400]);
px_array = res_array;
plot(px_array, auc_matt3_04(:,noise_level), 'k', 'LineWidth', 2);
hold on;
plot(px_array, auc_matt3_06(:,noise_level), 'k', 'LineWidth', 2);
plot(px_array, auc_matt3_08(:,noise_level), 'k', 'LineWidth', 2);
plot(px_array, auc_matt3_10(:,noise_level), 'k', 'LineWidth', 2);
plot(px_array, auc_matt3_15(:,noise_level), 'k', 'LineWidth', 2);
plot(px_array, auc_matt3_20(:,noise_level), 'k', 'LineWidth', 2);
plot(px_array, auc_matt3_25(:,noise_level), 'k', 'LineWidth', 2);
plot(px_array, auc_matt3_30(:,noise_level), 'k', 'LineWidth', 2);

ylim([0.5 1]);

%%
noise_level = 5;
figure('Position', [0 100 400 400]);
px_array = res_array;
plot(px_array, auc_matt3_04(:,noise_level), 'k', 'LineWidth', 2);
hold on;
plot(px_array, auc_matt3_06(:,noise_level), 'k', 'LineWidth', 2);
plot(px_array, auc_matt3_08(:,noise_level), 'k', 'LineWidth', 2);
plot(px_array, auc_matt3_10(:,noise_level), 'k', 'LineWidth', 2);
plot(px_array, auc_matt3_15(:,noise_level), 'k', 'LineWidth', 2);
plot(px_array, auc_matt3_20(:,noise_level), 'k', 'LineWidth', 2);
plot(px_array, auc_matt3_25(:,noise_level), 'k', 'LineWidth', 2);
plot(px_array, auc_matt3_30(:,noise_level), 'k', 'LineWidth', 2);
ylim([0.5 1]);

%%
noise_level = 1;
figure('Position', [0 100 400 400]);
px_array = res_array;
plot(px_array, auc_matt3_04(:,noise_level), 'k', 'LineWidth', 2);
hold on;
plot(px_array, auc_matt3_06(:,noise_level), 'k', 'LineWidth', 2);
plot(px_array, auc_matt3_08(:,noise_level), 'k', 'LineWidth', 2);
plot(px_array, auc_matt3_10(:,noise_level), 'k', 'LineWidth', 2);
plot(px_array, auc_matt3_15(:,noise_level), 'k', 'LineWidth', 2);
plot(px_array, auc_matt3_20(:,noise_level), 'k', 'LineWidth', 2);
plot(px_array, auc_matt3_25(:,noise_level), 'k', 'LineWidth', 2);
plot(px_array, auc_matt3_30(:,noise_level), 'k', 'LineWidth', 2);
ylim([0.5 1]);

