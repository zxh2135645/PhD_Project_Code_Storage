clear all;
close all;

addpath('../function/')
base_dir = uigetdir;
nonhemo_fu = load(cat(2, base_dir, '/FF_Data/484060000041/FU/Heart_Localizer - ZS0023213085_T2_star_mapping_GRE_2D_8echo_WB_25.mat'));
nonhemo_bl = load(cat(2, base_dir, '/FF_Data/484060000041/BL/Heart_Localizer - ZS0021574796_T2_star_mapping_GRE_2D_8echo_WB_24.mat'));
hemo_bl = load(cat(2, base_dir, '/FF_Data/484060000009/BL/T2star_8echo_db_1recovHB_new_31.mat'));
hemo_fu = load(cat(2, base_dir, '/FF_Data/484060000009/FU/T2star_8echo_wb_1recovHB_new_29.mat'));

%% colorbar
f = figure();
n_colors = 256;
h = axes( ...
    'Parent',f, ...
    'Units','normalized', ...
    'Position',[0.91 0.34 0.022 0.475], ...
    'XTick',[], ...
    'YTick',[], ...
    'YAxisLocation','right', ...
    'XLim',[0 1], ...
    'YLim',[0 n_colors-1], ...
    'Visible','off');
surface( ...
    'Parent',h, ...
    'XData',[0 1], ...
    'YData',0:n_colors-1, ...
    'ZData',zeros(n_colors,2), ...
    'FaceColor','flat', ...
    'CData',permute(jet(n_colors),[1 3 2]), ...
    'EdgeColor','none');

h = axes( ...
    'Parent',f, ...
    'Units','normalized', ...
    'Position',[0.5 0.34 0.022 0.475], ...
    'XTick',[], ...
    'YTick',[], ...
    'YAxisLocation','right', ...
    'XLim',[0 1], ...
    'YLim',[0 n_colors-1], ...
    'Visible','off');
values_brewermap = brewermap([],'*RdYlBu');
surface( ...
    'Parent',h, ...
    'XData',[0 1], ...
    'YData',0:n_colors-1, ...
    'ZData',zeros(n_colors,2), ...
    'FaceColor','flat', ...
    'CData',permute(values_brewermap,[1 3 2]), ...
    'EdgeColor','none');

%% Hemo-
figure(); 
imagesc(nonhemo_bl.fwmc_ff); caxis([0 40]); axis image; colorbar('TickLabels',{''}); colormap jet;
axis off; 

figure(); 
imagesc(nonhemo_fu.fwmc_ff); caxis([0 40]); axis image; colorbar('TickLabels',{''}); colormap jet;
axis off;

figure(); 
imagesc(nonhemo_bl.fwmc_r2star); caxis([0 100]); axis image; colorbar('TickLabels',{''}); colormap jet;
axis off; 

figure(); 
imagesc(nonhemo_fu.fwmc_r2star); caxis([0 100]); axis image; colorbar('TickLabels',{''}); colormap jet;
axis off;

%% Hemo- brewer
figure(); 
imagesc(nonhemo_bl.fwmc_ff); caxis([0 50]); axis image; colorbar('TickLabels',{''}); 
axis off; 
colormap(brewermap([],'*RdYlBu'));

figure(); 
imagesc(nonhemo_fu.fwmc_ff); caxis([0 50]); axis image; colorbar('TickLabels',{''}); 
axis off;
colormap(brewermap([],'*RdYlBu'));

figure(); 
imagesc(nonhemo_bl.fwmc_r2star); caxis([0 200]); axis image; colorbar('TickLabels',{''});
axis off; 
colormap(brewermap([],'*RdYlBu'));

figure(); 
imagesc(nonhemo_fu.fwmc_r2star); caxis([0 200]); axis image; colorbar('TickLabels',{''});
axis off;
colormap(brewermap([],'*RdYlBu'));

%%  Raw
figure(); 
imagesc(nonhemo_bl.fwmc_ff); caxis([0 50]); axis image; colorbar('TickLabels',{''}); colormap gray;
axis off; 

figure(); 
imagesc(nonhemo_fu.fwmc_ff); caxis([0 50]); axis image; colorbar('TickLabels',{''}); colormap gray;
axis off;

figure(); 
imagesc(nonhemo_bl.fwmc_r2star); caxis([0 200]); axis image; colorbar('TickLabels',{''}); colormap gray;
axis off; 

figure(); 
imagesc(nonhemo_fu.fwmc_r2star); caxis([0 200]); axis image; colorbar('TickLabels',{''}); colormap gray;
axis off;

%% Hemo+
figure(); 
imagesc(hemo_bl.fwmc_ff); caxis([0 40]); axis image; colorbar('TickLabels',{''}); colormap jet;
axis off; 

figure(); 
imagesc(hemo_fu.fwmc_ff); caxis([0 40]); axis image; colorbar('TickLabels',{''}); colormap jet;
axis off;

figure(); 
imagesc(hemo_bl.fwmc_r2star); caxis([0 100]); axis image; colorbar('TickLabels',{''}); colormap jet;
axis off; 

figure(); 
imagesc(hemo_fu.fwmc_r2star); caxis([0 100]); axis image; colorbar('TickLabels',{''}); colormap jet;
axis off;
%% Hemo+ brewer
figure(); 
imagesc(hemo_bl.fwmc_ff); caxis([0 50]); axis image; colorbar('TickLabels',{''}); colormap jet;
axis off; 
colormap(brewermap([],'*RdYlBu'));

figure(); 
imagesc(hemo_fu.fwmc_ff); caxis([0 50]); axis image; colorbar('TickLabels',{''}); colormap jet;
axis off;
colormap(brewermap([],'*RdYlBu'));

figure(); 
imagesc(hemo_bl.fwmc_r2star); caxis([0 200]); axis image; colorbar('TickLabels',{''}); colormap jet;
axis off; 
colormap(brewermap([],'*RdYlBu'));

figure(); 
imagesc(hemo_fu.fwmc_r2star); caxis([0 200]); axis image; colorbar('TickLabels',{''}); colormap jet;
axis off;
colormap(brewermap([],'*RdYlBu'));

%% Hemo+ raw
figure(); 
imagesc(hemo_bl.fwmc_ff); caxis([0 50]); axis image; colorbar('TickLabels',{''}); colormap gray;
axis off; 

figure(); 
imagesc(hemo_fu.fwmc_ff); caxis([0 50]); axis image; colorbar('TickLabels',{''}); colormap gray;
axis off;

figure(); 
imagesc(hemo_bl.fwmc_r2star); caxis([0 200]); axis image; colorbar('TickLabels',{''}); colormap gray;
axis off; 

figure(); 
imagesc(hemo_fu.fwmc_r2star); caxis([0 200]); axis image; colorbar('TickLabels',{''}); colormap gray;
axis off;
