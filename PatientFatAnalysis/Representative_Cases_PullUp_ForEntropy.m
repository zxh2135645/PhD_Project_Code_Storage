clear all;
close all;

addpath('../function/')
base_dir = uigetdir;
nonhemo_fu = load(cat(2, base_dir, '/FF_Data/484060000018/FU/Heart_Cardiac_Dot_Engine - ZS0022208988_T2star_8echo_wb_1recovHB_new_31.mat'));
nonhemo_bl = load(cat(2, base_dir, '/FF_Data/484060000018/BL/Heart_Cardiac_Dot_Engine - ZS0020044304_T2star_8echo_db_1recovHB_new_23.mat'));
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


%% T1 entropy map
figure();
%roi_in_myo_t1_nan = roi_in_myo_t1;
%roi_in_myo_t1_nan = double(roi_rimmed_t1_new);
%roi_in_myo_t1_nan(roi_in_myo_t1 == 0) = nan;
%roi_in_myo_t1_nan(roi_rimmed_t1_new == 0) = nan;

roi_in_myo_t1_nan = double(roi_in_myo_t1);
roi_in_myo_t1_nan(roi_in_myo_t1 == 0) = nan;

myo_t1_nan = double(myo_t1);
myo_t1_nan(roi_in_myo_t1 == 0) = nan;

myo_t1_nan = myo_t1;
myo_t1_nan(myo_t1 == 0) = nan;
%for slc = 1:size(roi_in_myo_t1_nan,3)
for slc = 2:2
    t1_slc = t1(:,:,slc);

    [x,y] = find(remote_in_myo_t1(:,:,slc) == 1);

    t1_array = zeros(length(x),1);
    for i = 1:length(x)
        t1_array(i) = t1(x(i),y(i),slc);
    end

    t1_array(t1_array < 0) = 0;
    mean_t1_remote = mean(t1_array);

    ax1 = axes;
    imagesc(t1_slc); pbaspect([size(t1_slc, 2), size(t1_slc, 1) 1]);
    colormap(ax1, 'gray');
    set(ax1, 'xticklabel', []); set(ax1, 'yticklabel', []);

    ax2 = axes;
    imagesc(ax2, myo_t1(:,:,slc), 'AlphaData', myo_t1(:,:,slc));
    pbaspect([size(t1_slc, 2), size(t1_slc, 1) 1]); colormap(ax2, 'cool');
    ax2.Visible = 'off';

    ax3 = axes;
    imagesc(ax3, myo_t1_nan(:,:,slc) .* t1_slc - mean_t1_remote, 'AlphaData', myo_t1(:,:,slc));
    pbaspect([size(t1_slc, 2), size(t1_slc, 1) 1]); colormap(ax3, 'cool');
    caxis(ax1, [0 2000]); caxis(ax2, [0 2]); caxis(ax3, [-500 500]); linkprop([ax1 ax2 ax3], 'Position');
    ax3.Visible = 'off';

end