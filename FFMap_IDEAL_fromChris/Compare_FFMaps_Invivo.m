clear all;
%close all;

addpath('../function/');
Scan_Type='Cardiac';

[HBfile, HBpath] = uigetfile('*.mat');
load(cat(2, HBpath, HBfile));

name = 'Sahara';
%name = 'Mojave'
%name = 'Merry'
time_point = '1YR'
label_t2star = 'T2star';
base_dir = uigetdir;
load(cat(2, base_dir, '/ContourData/', name, '/', name, '_', time_point, '/', label_t2star, '/', label_t2star, '_Index.mat'));

ff_map = cell(1, length(glob_names));
for f = 1:length(ff_map)
    ff_map{f} = load(cat(2,  base_dir, '/FF_Data/',  name, '/', name, '_', time_point, '_', glob_names{f}, '.mat'), 'fwmc_ff');
end

% convert ff_map to matrix
ff = zeros(size(ff_map{1}.fwmc_ff,1), size(ff_map{1}.fwmc_ff, 2), length(ff_map));
for f = 1:length(ff_map)
    ff(:,:,f) = ff_map{f}.fwmc_ff;
end

%% Comparison
ff_hb = abs(fat) ./ (abs(fat) + abs(water));
% ff_hb = real(fat) ./ (real(fat) + real(water));
filename = cat(2, base_dir, '/img/', name,  '/', name, '_', time_point, '/', 'FFMap_Comparison_QualGuided.gif');

h = figure();
for slc = 1:size(fat, 3)
    drawnow;
    subplot(1,2,1);
    imagesc(ff(:,:,slc)); caxis([0 20]); colormap jet; colorbar; axis image;
    subplot(1,2,2);
    imagesc(abs(ff_hb(:,:,slc))*100); caxis([0 20]); colormap jet; colorbar; axis image;
    pause(1);

    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    % Write to the GIF File
    if slc == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end

%% Masked
%ff_hb = real(fat) ./ (real(fat) + real(water));
filename = cat(2, base_dir, '/img/', name,  '/', name, '_', time_point, '/', 'FFMap_Comparison_masked_QualGuided_real.gif');

h = figure();
for slc = 1:size(fat, 3)
    drawnow;
    subplot(1,2,1);
    imagesc(ff(:,:,slc).*AllPhasemap(slc).Mask); caxis([0 20]); colormap jet; colorbar; axis image;
    subplot(1,2,2);
    imagesc(abs(ff_hb(:,:,slc))*100.*AllPhasemap(slc).Mask); caxis([0 20]); colormap jet; colorbar; axis image;
    pause(1);

    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    % Write to the GIF File
    if slc == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end

%%
anatomy_label = {'BloodPool', 'excludeArea', 'freeROI', 'Heart', 'Myocardium', 'MyoReference', 'noReflowArea'}; 

ff_hb = abs(fat) ./ (abs(fat) + abs(water));
save_dir = cat(2, base_dir, '/img/');
name_save_dir = cat(2, save_dir, name);
if ~exist(name_save_dir, 'dir')
    mkdir(name_save_dir);
end
tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', name, '_', time_point,  '/');

myo_glob = glob(cat(2, tp_dir, label_t2star, '/', anatomy_label{5}, '/*'));
roi_glob = glob(cat(2, tp_dir, label_t2star, '/',anatomy_label{3}, '/*'));
remote_glob = glob(cat(2, tp_dir, label_t2star, '/',anatomy_label{6}, '/*'));

load(cat(2, tp_dir, label_t2star, '/', label_t2star, '_vol_img_3D.mat'));
load(myo_glob{1});
load(roi_glob{1});
load(remote_glob{1});
load(cat(2, tp_dir, label_t2star, '/', label_t2star, '_SliceLoc.mat'));
[slc_array_ff, idx_reordered] = sort(slc_array);

roi_in_myo_ff = mask_myocardium_3D .* freeROIMask_3D;
remote_in_myo_ff = mask_myocardium_3D .* myoRefMask_3D;
roi_ff = roi_in_myo_ff .* vol_img_3D;
remote_ff = remote_in_myo_ff .* vol_img_3D;
myo_ff = mask_myocardium_3D;

roi_in_myo_ff = roi_in_myo_ff(:,:,idx_reordered);
remote_in_myo_ff = remote_in_myo_ff(:,:,idx_reordered);
roi_ff = roi_ff(:,:,idx_reordered);
remote_ff = remote_ff(:,:,idx_reordered);
ff = ff(:,:,idx_reordered);
myo_ff = myo_ff(:,:,idx_reordered);


Imgin = ff;
Maskin = myo_ff;
Segn = 50;
Groove = 0;

[Segmentpix, stats, Mask_Segn] =AHASegmentation(Imgin,Maskin,Segn,Groove);
figure();
for i = 1:size(Mask_Segn, 3)
    subplot(2,2,i);
    imagesc(Mask_Segn(:,:,i)); axis image;
end

%% AHA Segmentation and Plot BullsEye
% ========================== To simplify, name the reference point at (100, 100)
% ========================== To simplify, centroid is centroid of the first slice 
% Preprocess - 

ff(ff < 0) = 0;
ff(ff > 100) = 100;
se = strel('disk', 1);
myo_ff_eroded = imerode(myo_ff, se);

x = 50;
y = 50;
C = regionprops(myo_ff_eroded(:,:,1));
x_centroid = C.Centroid(2);
y_centroid = C.Centroid(1);
BaseGroove = atan2(x - x_centroid, y - y_centroid) * 180 / pi;
[~, idx_array] = sort(slc_array_ff); % assume minimal value is basal slice location
LocPixCount = AHA_16Seg(ff, myo_ff_eroded, BaseGroove, idx_array);
LocPixCount_hb = AHA_16Seg(ff_hb*100, myo_ff_eroded, BaseGroove, idx_array);




% =================== Create the AHA 17-segment bullseye =================
% ========================================================================
figure('Position', [400 400 800 800]);
PlotBullsEye(LocPixCount);
figure('Position', [400 400 800 800]);
PlotBullsEye(LocPixCount_hb);   

%% Analysis - (Linear Regression and Bland-Altman)
color_cell_x = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell_y = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};

mdl = fitlm(LocPixCount, LocPixCount_hb);

figure(); hold on;
Y = LocPixCount .* mdl.Coefficients.Estimate(2) + mdl.Coefficients.Estimate(1);
xlabel('FF (Without Phase  Map, %)');
ylabel('FF (With Phase  Map, %)');
%xlim([0 50]); ylim([0 150]);

scatter(LocPixCount, LocPixCount_hb,  64, 'MarkerEdgeColor', color_cell_y{5}, 'MarkerFaceColor', color_cell_x{3});
plot(LocPixCount, Y, 'b', 'LineWidth', 1);
yl = ylim;
xl = xlim;
text(0.3*xl(2), yl(1)*0.3, cat(2,'Y = ', num2str(mdl.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl.Rsquared.Ordinary,3)), 'FontSize', 12)

%% Bland-Altman
addpath('../function/BlandAltman/');
% Heart muscle territories per patient
% territories = {'9MO'};
territories = {time_point};
nterritories = length(territories);

% Patient states during measurement
states = {'Myocardium'};
nstates = length(states);

data1 = reshape(LocPixCount, [], 1, 1);
data2 = reshape(LocPixCount_hb, [], 1, 1);

% BA plot paramters
tit = 'FF Comparison'; % figure title
gnames = {territories, states}; % names of groups in data {dimension 1 and 2}
label = {'Fat Fraction (Diego)','Exvivo Fat Fraction (Homebrew)','%'}; % Names of data sets
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
[cr, fig, statsStruct] = BlandAltman(data1, data2,label,tit,gnames,'corrInfo',corrinfo,'baInfo',BAinfo,'axesLimits',limits,'colors',colors, 'showFitCI',' on');

