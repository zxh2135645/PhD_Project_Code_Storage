clear all;
close all;

addpath('../function/');
Scan_Type='Cardiac';

% [HBfile, HBpath] = uigetfile('*.mat');
% load(cat(2, HBpath, HBfile));


base_dir = uigetdir;
ff_hb = load(cat(2, base_dir, '/Result_FromWorkstation_07052022/AllPhasemap_SAHARA.mat'));
ff_hb = load(cat(2, base_dir, '/Result_FromWorkstation_07052022/AllPhasemap_RYN.mat'));
%ff_phase = real(ff_hb.fat) ./ (real(ff_hb.fat) + real(ff_hb.water));
ff_diego = load(cat(2, base_dir, '/FF_Data_06232022/7777-00170_SAHARA_EXVIVO_06222022/3D_IDEAL_Avg4_Bipolar.mat'));
ff_diego = load(cat(2, base_dir, '/FF_Data_06232022/RYN_EXVIVO3_06092022/3D_IDEAL_Avg4_Bipolar_0051.mat'));


%% Show Image
ff_phase = abs(ff_hb.fat_r2s) ./ (abs(ff_hb.fat_r2s) + abs(ff_hb.water_r2s));
ff_hd = ff_diego.fwmc_ff;

ff_hd_trunc = ff_hd(:,:,49:52);
ff_phase = ff_phase(:,:,49:52);
figure();
for i = 1:size(ff_phase, 3)
    subplot(2,2,i);
    imagesc(transpose(ff_phase(:,:,i))); caxis([0 0.2]); axis image;
end

figure();
for i = 1:size(ff_hd_trunc, 3)
    subplot(2,2,i);
    imagesc(ff_hd_trunc(:,:,i)); caxis([0 20]); axis image;
end

%% Save image
name = 'Sahara';
name = 'Ryn';
time_point = '1YR'

filename = cat(2, base_dir, '/../MATLAB/T1_Fat_Project/img/', name,  '/', name, '_', time_point, '/', 'FFMap_Comparison_QualGuided_Exvivo_r2s.gif');

h = figure();
for slc = 1:size(ff_phase, 3)
    drawnow;
    subplot(1,2,1);
    imagesc(ff_hd_trunc(:,:,slc)); caxis([0 20]); colormap jet; colorbar; axis image;
    subplot(1,2,2);
    imagesc(abs(transpose(ff_phase(:,:,slc)))*100); caxis([0 20]); colormap jet; colorbar; axis image;
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

%% Draw mask
figure();
for slc = 1:size(ff_phase, 3)
    ff_phase_trans(:,:,slc) = transpose(ff_phase(:,:,slc));
    imagesc(ff_phase_trans(:,:,slc)); axis image; axis off; colormap gray; caxis([0 0.2]); title(cat(2, 'ROI, slice ', num2str(slc)));
    roi = drawpolygon;
    mask_roi(:,:,slc) = createMask(roi);
end

for slc = 1:size(ff_phase, 3)
    ff_phase_trans(:,:,slc) = transpose(ff_phase(:,:,slc));
    imagesc(ff_phase_trans(:,:,slc)); axis image; axis off; colormap gray; caxis([0 0.2]); title(cat(2, 'ENDO, slice ', num2str(slc)));
    roi = drawpolygon;
    mask_endo(:,:,slc) = createMask(roi);
end

%% AHA Segmentation and Plot BullsEye
addpath('../AHA16Segment/');

myo_ff = mask_roi - mask_endo;
ff_hd_trunc(ff_hd_trunc < 0) = 0;
ff_hd_trunc(ff_hd_trunc > 100) = 100;
se = strel('disk', 2);
myo_ff_eroded = imerode(myo_ff, se);

idx_array = 1:size(ff_phase_trans, 3);
BaseGroove = 0;
LocPixCount = AHA_16Seg(ff_hd_trunc, myo_ff_eroded, BaseGroove, idx_array);
LocPixCount_hb = AHA_16Seg(ff_phase_trans*100, myo_ff_eroded, BaseGroove, idx_array);

% =================== Create the AHA 17-segment bullseye =================
% ========================================================================
figure('Position', [400 400 800 800]);
PlotBullsEye(LocPixCount);
figure('Position', [400 400 800 800]);
PlotBullsEye(LocPixCount_hb);  

%% 50 chords
Nseg = 50;
[Segmentpix_hd, stats_hd, Mask_Segn_hd] = AHASegmentation(ff_hd_trunc, myo_ff_eroded, Nseg, BaseGroove);
[Segmentpix_hb, stats_hb, Mask_Segn_hb] = AHASegmentation(ff_phase_trans*100, myo_ff_eroded, Nseg, BaseGroove);

mean_hd = reshape(permute(stats_hd(1,:,:), [2,1,3]), Nseg, []);
mean_hb = reshape(permute(stats_hb(1,:,:), [2,1,3]), Nseg, []);
std_hd = reshape(permute(stats_hd(2,:,:), [2,1,3]), Nseg, []);
std_hb = reshape(permute(stats_hb(2,:,:), [2,1,3]), Nseg, []);
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



% data1 = reshape(LocPixCount, [], 1, 1);
% data2 = reshape(LocPixCount_hb, [], 1, 1);

% idx = find(myo_ff_eroded ~= 0);
% data1 = reshape(ff_hd_trunc(idx), [], 1, 1);
% data2 = reshape(ff_phase_trans(idx)*100, [], 1, 1);

data1 = reshape(mean_hd(:), [], 1, 1);
data2 = reshape(mean_hb(:), [], 1, 1);

% data1 = reshape(std_hd(:), [], 1, 1);
% data2 = reshape(std_hb(:), [], 1, 1);

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


