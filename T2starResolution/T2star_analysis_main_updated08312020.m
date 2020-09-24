% This script is specifically used for workflow updated in 08/28/2020
% The updated detail please refer to logbook in benchling

clear all;
close all;

addpath('../function/');

%%%%%%%%%%%%%%%%%%%%%%%%%% output file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% roi
% mask_struct
% aha_anlysis
% T2star_meanSD_table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 20P10
% The starting with 0014
base_dir = uigetdir;
folder_glob = glob(cat(2, base_dir, '\*'));

labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};

label = labels{5};
idx_array = contains(folder_glob, label);

[list_to_read, order_to_read] = NamePicker(folder_glob);

proj_dir = GetFullPath(cat(2, pwd, '/../../T2star_Resolution_Project/'));
if ~exist(proj_dir, 'dir')
    mkdir(proj_dir)
end

save_dir = GetFullPath(cat(2, proj_dir, 'img/'));
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

data_dir = GetFullPath(cat(2, proj_dir, 'data/'));
if ~exist(data_dir, 'dir')
    mkdir(data_dir)
end

subject_name = input('Please type subject name here:  ', 's');
subject_dir = GetFullPath(cat(2, save_dir, subject_name, '/'));
if ~exist(subject_dir, 'dir')
    mkdir(subject_dir)
end
subject_data_dir = GetFullPath(cat(2, data_dir, subject_name, '/'));
if ~exist(subject_data_dir, 'dir')
    mkdir(subject_data_dir)
end

disp('Avg 0016 starts here: ');
avg_num = input('Please type average number here:  ');
avg_name = cat(2, 'Avg', num2str(avg_num, '%04.f'));

%% Read T2* DICOM files
whatsinit = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    whatsinit{i} = dicom23D(f);
end

% Display images
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    subplot(4,7,i);
    imagesc(whatsinit{i}); axis image;
    caxis([0 100])
end

%% Draw contours @ epi, endo, MI, remote, fluid
img = whatsinit{1};
myo_coords_cell = cell(size(img, 3), 2);
roi_save = cat(2, subject_data_dir, 'roi.mat');

if ~exist(roi_save, 'file')
for i = 1:(size(img, 3))
    disp('Epicardium: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 100])
    epi = drawpolygon(gca);
    epi_coords = epi.Position;
    
    disp('Endocardium: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 100])
    endo = drawpolygon(gca);
    endo_coords = endo.Position;
    
    disp('MI roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 100])
    mi = drawpolygon(gca);
    mi_coords = mi.Position;
    
    disp('remote roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 100])
    remote = drawpolygon(gca);
    remote_coords = remote.Position;
    
    disp('fluid roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 100])
    fluid = drawpolygon(gca);
    fluid_coords = fluid.Position;
    
    disp('Center line roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 100])
    center_line = drawpolygon(gca);
    center_coords = center_line.Position;
    
    disp('air roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 100])
    air = drawpolygon(gca);
    air_coords = air.Position;
    
    myo_coords_cell{i, 1} = epi.Position;
    myo_coords_cell{i, 2} = endo.Position;
    epi_mask = createMask(epi);
    endo_mask = createMask(endo);
    center_mask = createMask(center_line);
    
    close all;
end

roi.myo_coords_cell = myo_coords_cell;
roi.mi_coords = mi_coords;
roi.remote_coords = remote_coords;
roi.fluid_coords = fluid_coords;
roi.center_coords = center_coords;
roi.air_coords = air_coords;

save(roi_save, 'roi');

else
    load(roi_save);
    myo_coords_cell = roi.myo_coords_cell;
    mi_coords = roi.mi_coords;
    remote_coords = roi.remote_coords;
    fluid_coords = roi.fluid_coords;
    center_coords = roi.center_coords;
    air_coords = roi.air_coords;
end

% 28 different set of parameters
% Convert coords to masks for 28 images
img_size_truth = size(img);
mask_save = cat(2, subject_data_dir, 'mask.mat');

if ~exist(mask_save, 'file')
    figure();
    mask_struct = struct;
    for i = 1:length(whatsinit)
        img2 = whatsinit{i};
        img2_size = size(whatsinit{i});
        ratio = img_size_truth(1) ./ img2_size(1);
        imagesc(img2); caxis([0 100]);
        epi = drawpolygon(gca,'Position', [myo_coords_cell{1}(:,1)/ratio + (ratio-1)/ratio, myo_coords_cell{1}(:,2)/ratio + (ratio-1)/ratio]);
        endo = drawpolygon(gca,'Position', [myo_coords_cell{2}(:,1)/ratio + (ratio-1)/ratio, myo_coords_cell{2}(:,2)/ratio + (ratio-1)/ratio]);
        mi = drawpolygon(gca,'Position', [mi_coords(:,1)/ratio + (ratio-1)/ratio, mi_coords(:,2)/ratio + (ratio-1)/ratio]);
        remote = drawpolygon(gca,'Position', [remote_coords(:,1)/ratio + (ratio-1)/ratio, remote_coords(:,2)/ratio + (ratio-1)/ratio]);
        fluid = drawpolygon(gca,'Position', [fluid_coords(:,1)/ratio + (ratio-1)/ratio, fluid_coords(:,2)/ratio + (ratio-1)/ratio]);
        air = drawpolygon(gca,'Position', [air_coords(:,1)/ratio + (ratio-1)/ratio, air_coords(:,2)/ratio + (ratio-1)/ratio]);
        center_line = drawpolygon(gca,'Position', [center_coords(:,1)/ratio + (ratio-1)/ratio, center_coords(:,2)/ratio + (ratio-1)/ratio]);
        
        epi_mask = createMask(epi);
        endo_mask = createMask(endo);
        myo_mask = epi_mask - endo_mask;
        
        mi_mask = createMask(mi);
        remote_mask = createMask(remote);
        fluid_mask = createMask(fluid);
        air_mask = createMask(air);
        center_mask = createMask(center_line);
        
        mask_struct(i).myo_mask = myo_mask;
        mask_struct(i).mi_mask = mi_mask;
        mask_struct(i).remote_mask = remote_mask;
        mask_struct(i).fluid_mask = fluid_mask;
        mask_struct(i).air_mask = air_mask;
        
        mask_struct(i).epi_mask = epi_mask;
        mask_struct(i).endo_mask = endo_mask;
        
        myo_mask_endo = myo_mask .* center_mask;
        myo_mask_epi = myo_mask - myo_mask_endo;
        mask_struct(i).myo_mask_endo = myo_mask_endo;
        mask_struct(i).myo_mask_epi = myo_mask_epi;
    end
    
    save(mask_save, 'mask_struct');
else
    load(mask_save);
end

%% Display image overlay with Mean-2SD mask
save_array = 1:1:length(whatsinit);
for i = 1:length(whatsinit)
    save_idx = save_array(i);
    figure();
    img2 = whatsinit{save_idx};
    thresh = mean(nonzeros(img2 .* mask_struct(save_idx).remote_mask)) - 2*std(nonzeros(img2 .* mask_struct(save_idx).remote_mask));
    hemo_mask = img2 < thresh;
    subplot(1,2,1);
    imagesc(img2); caxis([0 50]); axis image; colormap default; %colorbar;
    
    %img_rgb = zeros(size(img2,1), size(img2, 2), 3);
    %img2_norm = img2 ./ max(img2(:));
    %img_rgb(:,:,1) = img2_norm;
    %img_rgb(:,:,2) = img2_norm;
    %img_rgb(:,:,3) = img2_norm;
    
    %map = hsv(256);
    map = gray(256);
    rgbImage = ind2rgb(img2, map);
    rgbImage(:,:,1) = rgbImage(:,:,1) + hemo_mask.*mask_struct(save_idx).mi_mask;
    subplot(1,2,2);
    imagesc(rgbImage); axis image;
    
    
    mean2sd_dir = cat(2, subject_dir, 'Mean2SD_', avg_name, '/');
    if ~exist(mean2sd_dir, 'dir')
        mkdir(mean2sd_dir)
    end
    saveas(gcf, cat(2, mean2sd_dir, num2str(save_idx), '.png'));
end

close all;
%% Otsu segmentation
num_cluster = 3;

for i = 1:length(whatsinit)
    Img = whatsinit{i};
    img_to_otsu = Img .* mask_struct(i).mi_mask;
    img_to_otsu(img_to_otsu == 0) = NaN;
    img_to_otsu(img_to_otsu > 50) = 50;
    
    [IDX, seg] = otsu(img_to_otsu, num_cluster);
    % Cluster1 - water; Cluster2 - myocardium; Cluster3 - hemorrhage
    
    otsu_mean = zeros(num_cluster, 1);
    for c = 1:num_cluster
        otsu_mean(c) = mean(Img(IDX == c));
    end
    
    [M,I] = min(otsu_mean);
    
    mask_struct(i).hemo_mask = IDX == I;
    
    figure();
    subplot(1,2,1);
    imagesc(Img); caxis([0 50]); axis image; % colorbar;
    
    map = gray(256);
    rgbImage = ind2rgb(Img, map);
    rgbImage(:,:,1) = rgbImage(:,:,1) + mask_struct(i).hemo_mask;
    subplot(1,2,2);
    imagesc(rgbImage); axis image;
    
    otsu_dir = cat(2, subject_dir, 'Otsu_', num2str(num_cluster), '_', avg_name, '/');
    if ~exist(otsu_dir, 'dir')
        mkdir(otsu_dir)
    end
    saveas(gcf, cat(2, otsu_dir, num2str(i), '.png'));
end

close all;

%% Mean + SD plot 
t2star_mean_array = zeros(1, length(whatsinit));
t2star_sd_array = zeros(1, length(whatsinit));
for i = 1:length(whatsinit)
    t2star_mean_array(i) = mean(nonzeros(mask_struct(i).remote_mask .* whatsinit{i}));
    t2star_sd_array(i) = std(nonzeros(mask_struct(i).remote_mask .* whatsinit{i}));
end

d = length(whatsinit)/4;
t2star_mean_reshape = reshape(t2star_mean_array, d, length(whatsinit)/d).';
t2star_sd_reshape = reshape(t2star_sd_array, d, length(whatsinit)/d).';

x = [0, d, 2*d, 3*d, 4*d] + [0 , 0.5, 0.5, 0.5, 1];
inplane_res = 1:d;
res = [inplane_res; inplane_res + d; inplane_res + 2*d; inplane_res + 3*d];
figure('Position', [100 0 1600 1600]);
errorbar(t2star_mean_array, t2star_sd_array, 'LineStyle', 'none' );
hold on;
ylim_lb = min(ylim); ylim_ub = max(ylim);
patch([x(1) x(2) x(2) x(1)], [max(ylim) max(ylim) 0 0], [241 194 151]/255, 'FaceAlpha',.5)
patch([x(2) x(3) x(3) x(2)], [max(ylim) max(ylim) 0 0], [199 213 161]/255, 'FaceAlpha',.5)
patch([x(3) x(4) x(4) x(3)], [max(ylim) max(ylim) 0 0], [159 203 219]/255, 'FaceAlpha',.5)
patch([x(4) x(5) x(5) x(4)], [max(ylim) max(ylim) 0 0], [98 141 207]/255, 'FaceAlpha',.5)
errorbar(res', t2star_mean_reshape', t2star_sd_reshape', '-o', 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980]); 
xticks([1.2 4 6.5 8.5 11 13.5 15.5 18 20.5 22.5 25 27.5]);
xticklabels({'0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1'})
xlim([0 x(5)]);ylim([ylim_lb, ylim_ub])

text(2,ylim_ub-2, 'Slice Thickness = 2 mm', 'FontSize', 16);
text(9,ylim_ub-2, 'Slice Thickness = 4 mm', 'FontSize', 16);
text(16,ylim_ub-2, 'Slice Thickness = 6 mm', 'FontSize', 16);
text(23,ylim_ub-2, 'Slice Thickness = 8 mm', 'FontSize', 16)
set(gca, 'FontSize', 16);
xlabel('Resolution (mm^2)', 'FontSize', 24); ylabel('T2^* (ms)', 'FontSize', 24);
hold off

%% Mean + SD Heatmap
figure('Position', [100 0 1600 1600]);
subplot(2,1,1);
imagesc(t2star_mean_reshape); colorbar; axis off;
for i = 1:size(t2star_mean_reshape, 1)
    for j = 1:size(t2star_mean_reshape, 2)
        text(0.9+(j-1),1+(i-1),num2str(round(t2star_mean_reshape(i,j),2)), 'Color', [0.8500, 0.3250, 0.0980], 'FontSize', 16);
    end
end
title('Mean of Remote T2* value', 'FontSize', 24);

subplot(2,1,2);
imagesc(t2star_sd_reshape); colorbar; axis off;
for i = 1:size(t2star_sd_reshape, 1)
    for j = 1:size(t2star_sd_reshape, 2)
        text(0.9+(j-1),1+(i-1),num2str(round(t2star_sd_reshape(i,j),2)), 'Color', [0.8500, 0.3250, 0.0980], 'FontSize', 16);
    end
end
title('Standard Deviation of Remote T2* value', 'FontSize', 24);

%% Mean + SD plot (boxplot)
clear t2star_remote_array
t2star_remote_array = [];
len_array = zeros(length(whatsinit), 1);
res_array = {'03', '06', '08', '10', '13', '16', '21'};
slc_array = {'2', '4', '6', '8'};
g = {};
for i = 1:length(whatsinit)
    temp = nonzeros(mask_struct(i).remote_mask .* whatsinit{i});
    len_array(i) = length(temp);
    t2star_remote_array = [t2star_remote_array; temp];
    q = fix((i-1)/7)+1;
    r = mod(i, 7);
    if r == 0
        real_r = 7 - r;
    else
        real_r = r;
    end
    g_temp = repmat({cat(2, slc_array{q}, '_', res_array{real_r})}, len_array(i), 1);
    g = [g; g_temp];
end

y_lim = [min(t2star_remote_array) max(t2star_remote_array)];
x = [0, d, 2*d, 3*d, 4*d] + [0 , 0.5, 0.5, 0.5, 1];
inplane_res = 1:d;
res = [inplane_res; inplane_res + d; inplane_res + 2*d; inplane_res + 3*d];
figure('Position', [100 0 1600 1600]);
errorbar(t2star_mean_array, t2star_sd_array, 'LineStyle', 'none' );
hold on;
ylim_lb = min(ylim); ylim_ub = max(ylim);
patch([x(1) x(2) x(2) x(1)], [max(y_lim) max(y_lim) 0 0], [241 194 151]/255, 'FaceAlpha',.5)
patch([x(2) x(3) x(3) x(2)], [max(y_lim) max(y_lim) 0 0], [199 213 161]/255, 'FaceAlpha',.5)
patch([x(3) x(4) x(4) x(3)], [max(y_lim) max(y_lim) 0 0], [159 203 219]/255, 'FaceAlpha',.5)
patch([x(4) x(5) x(5) x(4)], [max(y_lim) max(y_lim) 0 0], [98 141 207]/255, 'FaceAlpha',.5)
h = boxplot(t2star_remote_array,g);
set(h,'LineWidth',2); grid on;
ylim([min(t2star_remote_array), max(t2star_remote_array)]);

%% Table
VarNames = {' ', '0.3x0.3 mm^2', '0.6x0.6 mm^2',  '0.8x0.8 mm^2', '1.0x1.0 mm^2', '1.3x1.3 mm^2', '1.6x1.6 mm^2', '2.1x2.1 mm^2'};
mean_table = table({'2 mm'; '4 mm'; '6 mm'; '8 mm'},t2star_mean_reshape(:,1), t2star_mean_reshape(:,2), t2star_mean_reshape(:,3),...
    t2star_mean_reshape(:,4), t2star_mean_reshape(:,5), t2star_mean_reshape(:,6), t2star_mean_reshape(:,7), 'VariableNames',VarNames)

sd_table = table({'2 mm'; '4 mm'; '6 mm'; '8 mm'},t2star_sd_reshape(:,1), t2star_sd_reshape(:,2), t2star_sd_reshape(:,3),...
    t2star_sd_reshape(:,4), t2star_sd_reshape(:,5), t2star_sd_reshape(:,6), t2star_sd_reshape(:,7), 'VariableNames',VarNames)

t2star_table = struct;
t2star_table.mean_table = mean_table;
t2star_table.sd_table = sd_table;
save_table_f = cat(2, subject_data_dir, 'T2star_meanSD_table_', avg_name, '.mat');
save(save_table_f, 'save_table_f');

%% AHA
addpath('../AHA16Segment/');
Segn = 50;
Groove = 0;
aha50 = struct;
for i = 1:length(whatsinit)
    [Segmentpix, stats, Mask_Segn] = AHASegmentation(whatsinit{i}, mask_struct(i).myo_mask, Segn, Groove, mask_struct(i).endo_mask);
    aha50(i).Segmentpix = Segmentpix;
    aha50(i).Mask_Segn = Mask_Segn;
end

%figure(); imagesc(aha50(1).Mask_Segn.* mask_struct(1).myo_mask);

% AHA analysis
perc_array = zeros(length(whatsinit),Segn);

for i = 1:length(whatsinit)
    for j = 1:Segn
        thresh = mean(nonzeros(whatsinit{i} .* mask_struct(i).remote_mask)) - 2*std(nonzeros(whatsinit{i} .* mask_struct(i).remote_mask));
        perc_array(i,j) = sum(aha50(i).Segmentpix{j}<thresh) / length(aha50(i).Segmentpix{j});
    end
end

res = perc_array > 0.1;

%% 50 Segments in MI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aha_mi2 = struct;
for i = 1:length(whatsinit)
    Mipix = cell(Segn, 1);
    Img = whatsinit{i};
    for j = 1:Segn
        Mipix{j,1} = Img(aha50(i).Mask_Segn .* mask_struct(i).mi_mask == j);
    end
    aha_mi2(i).Mipix = Mipix;
end

Mipix_flat = {};
for i = 1:length(whatsinit)
    idx = 1;
        for j = 1:Segn
            if ~isempty(aha_mi2(1).Mipix{j,1})
                Mipix_flat{idx} = aha_mi2(i).Mipix{j,1};
                idx = idx + 1;
            end
        end
    aha_mi2(i).Mipix_flat = Mipix_flat;
end

% AHA analysis
perc_array_mi2 = zeros(length(whatsinit),length(Mipix_flat));

for i = 1:length(whatsinit)
    for j = 1:length(Mipix_flat)
        thresh = mean(nonzeros(whatsinit{i} .* mask_struct(i).remote_mask)) - 2*std(nonzeros(whatsinit{i} .* mask_struct(i).remote_mask));
        perc_array_mi2(i,j) = sum(aha_mi2(i).Mipix_flat{j}<thresh) / length(aha_mi2(i).Mipix_flat{j});
    end
end


res_mi2 = perc_array_mi2 > 0.1;

%%
figure(); imagesc(aha50(1).Mask_Segn.* mask_struct(1).mi_mask);
%% Confusion Matrix
figure('Position', [100 0 1600 1600]);
sens_mi2 = zeros(length(whatsinit),1);
spec_mi2 = zeros(length(whatsinit),1);
for i = 1:length(whatsinit)
    [cm, order] = confusionmat(res_mi2(1,:),res_mi2(i,:));
    subplot(4,7,i);confusionchart(cm, order); 
    set(gca, 'FontSize', 18);
    sens_mi2(i) = cm(2,2) / (cm(2,2) + cm(2,1));
    spec_mi2(i) = cm(1,1) / (cm(1,1) + cm(1,2));
end

% Plot sensitivity and specificity
sens_reshape = reshape(sens_mi2, [7, 4])';
spec_reshape = reshape(spec_mi2, [7, 4])';
figure();
subplot(2,1,1);
imagesc(sens_reshape); axis image
title('Sensitivity');
set(gca, 'FontSize', 18); colorbar;
subplot(2,1,2);
imagesc(spec_reshape); axis image
title('Specificity');
set(gca, 'FontSize', 18); colorbar;

aha_analysis.perc_array_mi = perc_array_mi2;
aha_analysis.sens_mi = sens_mi2;
aha_analysis.spec_mi = spec_mi2;
aha_analysis.aha_mi = aha_mi2;

% ROC analysis %%
auc_array_mi2 = zeros(length(whatsinit),1);
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    [X,Y,T,AUC,OPTROCPT] = perfcurve(res_mi2(1,:),perc_array_mi2(i,:), 1);
    subplot(4,7,i);
    plot(X,Y, 'LineWidth', 2);
    xlabel('FPR');
    ylabel('TPR');
    text(0.5,0.5,num2str(round(AUC, 2)),'FontSize',24);
    % title('ROC for Classification')
    set(gca, 'FontSize', 16);
    auc_array_mi2(i) = AUC;
end

aha_analysis.auc_array_mi = auc_array_mi2;

% Plot AUC %%
auc_reshape = reshape(auc_array_mi2, [7, 4])';
figure();
imagesc(auc_reshape); axis image
title('Area Under Curve');
set(gca, 'FontSize', 18); colorbar;

save(cat(2, subject_data_dir, 'aha_analysis_', avg_name, '.mat'), 'aha_analysis');
%% 2. For Avg0001
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%#########################################################################%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Avg0001 Starts here: ');
[list_to_read, order_to_read] = NamePicker(folder_glob);

avg_num = input('Please type average number here:  ');
avg_name = cat(2, 'Avg', num2str(avg_num, '%04.f'));

% Read T2* DICOM files
whatsinit = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    whatsinit{i} = dicom23D(f);
end

%% Display images
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    subplot(4,7,i);
    imagesc(whatsinit{i}); axis image;
    caxis([0 100])
end

%% Display image overlay with Mean-2SD mask
save_array = 1:1:length(whatsinit);
for i = 1:length(whatsinit)
    save_idx = save_array(i);
    figure();
    img2 = whatsinit{save_idx};
    thresh = mean(nonzeros(img2 .* mask_struct(save_idx).remote_mask)) - 2*std(nonzeros(img2 .* mask_struct(save_idx).remote_mask));
    hemo_mask = img2 < thresh;
    subplot(1,2,1);
    imagesc(img2); caxis([0 50]); axis image; %colorbar;
    map = gray(256);
    rgbImage = ind2rgb(img2, map);
    rgbImage(:,:,1) = rgbImage(:,:,1) + hemo_mask.*mask_struct(save_idx).mi_mask;
    subplot(1,2,2);
    imagesc(rgbImage); axis image;
    
    
    mean2sd_dir = cat(2, subject_dir, 'Mean2SD_', avg_name, '/');
    if ~exist(mean2sd_dir, 'dir')
        mkdir(mean2sd_dir)
    end
    saveas(gcf, cat(2, mean2sd_dir, num2str(save_idx), '.png'));
end

close all;
%% Otsu segmentation
num_cluster = 3;

for i = 1:length(whatsinit)
    Img = whatsinit{i};
    img_to_otsu = Img .* mask_struct(i).mi_mask;
    img_to_otsu(img_to_otsu == 0) = NaN;
    img_to_otsu(img_to_otsu > 50) = 50;
    
    [IDX, seg] = otsu(img_to_otsu, num_cluster);
    % Cluster1 - water; Cluster2 - myocardium; Cluster3 - hemorrhage
    
    otsu_mean = zeros(num_cluster, 1);
    for c = 1:num_cluster
        otsu_mean(c) = mean(Img(IDX == c));
    end
    
    [M,I] = min(otsu_mean);
    
    mask_struct(i).hemo_mask = IDX == I;
    
    figure();
    subplot(1,2,1);
    imagesc(Img); caxis([0 50]); axis image; %colorbar;
    map = gray(256);
    rgbImage = ind2rgb(Img, map);
    rgbImage(:,:,1) = rgbImage(:,:,1) + mask_struct(i).hemo_mask;
    subplot(1,2,2);
    imagesc(rgbImage); axis image;
    
    otsu_dir = cat(2, subject_dir, 'Otsu_', num2str(num_cluster), '_', avg_name, '/');
    if ~exist(otsu_dir, 'dir')
        mkdir(otsu_dir)
    end
    saveas(gcf, cat(2, otsu_dir, num2str(i), '.png'));
end

close all;

%% Mean + SD plot 
t2star_mean_array = zeros(1, length(whatsinit));
t2star_sd_array = zeros(1, length(whatsinit));
for i = 1:length(whatsinit)
    t2star_mean_array(i) = mean(nonzeros(mask_struct(i).remote_mask .* whatsinit{i}));
    t2star_sd_array(i) = std(nonzeros(mask_struct(i).remote_mask .* whatsinit{i}));
end

d = length(whatsinit)/4;
t2star_mean_reshape = reshape(t2star_mean_array, d, length(whatsinit)/d).';
t2star_sd_reshape = reshape(t2star_sd_array, d, length(whatsinit)/d).';

x = [0, d, 2*d, 3*d, 4*d] + [0 , 0.5, 0.5, 0.5, 1];
inplane_res = 1:d;
res = [inplane_res; inplane_res + d; inplane_res + 2*d; inplane_res + 3*d];
figure('Position', [100 0 1600 1600]);
errorbar(t2star_mean_array, t2star_sd_array, 'LineStyle', 'none' );
hold on;
ylim_lb = min(ylim); ylim_ub = max(ylim);
patch([x(1) x(2) x(2) x(1)], [max(ylim) max(ylim) 0 0], [241 194 151]/255, 'FaceAlpha',.5)
patch([x(2) x(3) x(3) x(2)], [max(ylim) max(ylim) 0 0], [199 213 161]/255, 'FaceAlpha',.5)
patch([x(3) x(4) x(4) x(3)], [max(ylim) max(ylim) 0 0], [159 203 219]/255, 'FaceAlpha',.5)
patch([x(4) x(5) x(5) x(4)], [max(ylim) max(ylim) 0 0], [98 141 207]/255, 'FaceAlpha',.5)
errorbar(res', t2star_mean_reshape', t2star_sd_reshape', '-o', 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980]); 
xticks([1.2 4 6.5 8.5 11 13.5 15.5 18 20.5 22.5 25 27.5]);
xticklabels({'0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1'})
xlim([0 x(5)]);ylim([ylim_lb, ylim_ub])

text(2,ylim_ub-2, 'Slice Thickness = 2 mm', 'FontSize', 16);
text(9,ylim_ub-2, 'Slice Thickness = 4 mm', 'FontSize', 16);
text(16,ylim_ub-2, 'Slice Thickness = 6 mm', 'FontSize', 16);
text(23,ylim_ub-2, 'Slice Thickness = 8 mm', 'FontSize', 16)
set(gca, 'FontSize', 16);
xlabel('Resolution (mm^2)', 'FontSize', 24); ylabel('T2^* (ms)', 'FontSize', 24);
hold off

%% Mean + SD Heatmap
figure('Position', [100 0 1600 1600]);
subplot(2,1,1);
imagesc(t2star_mean_reshape); colorbar; axis off;
for i = 1:size(t2star_mean_reshape, 1)
    for j = 1:size(t2star_mean_reshape, 2)
        text(0.9+(j-1),1+(i-1),num2str(round(t2star_mean_reshape(i,j),2)), 'Color', [0.8500, 0.3250, 0.0980], 'FontSize', 16);
    end
end
title('Mean of Remote T2* value', 'FontSize', 24);

subplot(2,1,2);
imagesc(t2star_sd_reshape); colorbar; axis off;
for i = 1:size(t2star_sd_reshape, 1)
    for j = 1:size(t2star_sd_reshape, 2)
        text(0.9+(j-1),1+(i-1),num2str(round(t2star_sd_reshape(i,j),2)), 'Color', [0.8500, 0.3250, 0.0980], 'FontSize', 16);
    end
end
title('Standard Deviation of Remote T2* value', 'FontSize', 24);

%% Mean + SD plot (boxplot)
clear t2star_remote_array
t2star_remote_array = [];
len_array = zeros(length(whatsinit), 1);
res_array = {'03', '06', '08', '10', '13', '16', '21'};
slc_array = {'2', '4', '6', '8'};
g = {};
for i = 1:length(whatsinit)
    temp = nonzeros(mask_struct(i).remote_mask .* whatsinit{i});
    len_array(i) = length(temp);
    t2star_remote_array = [t2star_remote_array; temp];
    q = fix((i-1)/7)+1;
    r = mod(i, 7);
    if r == 0
        real_r = 7 - r;
    else
        real_r = r;
    end
    g_temp = repmat({cat(2, slc_array{q}, '_', res_array{real_r})}, len_array(i), 1);
    g = [g; g_temp];
end

y_lim = [min(t2star_remote_array) max(t2star_remote_array)];
x = [0, d, 2*d, 3*d, 4*d] + [0 , 0.5, 0.5, 0.5, 1];
inplane_res = 1:d;
res = [inplane_res; inplane_res + d; inplane_res + 2*d; inplane_res + 3*d];
figure('Position', [100 0 1600 1600]);
errorbar(t2star_mean_array, t2star_sd_array, 'LineStyle', 'none' );
hold on;
ylim_lb = min(ylim); ylim_ub = max(ylim);
patch([x(1) x(2) x(2) x(1)], [max(y_lim) max(y_lim) 0 0], [241 194 151]/255, 'FaceAlpha',.5)
patch([x(2) x(3) x(3) x(2)], [max(y_lim) max(y_lim) 0 0], [199 213 161]/255, 'FaceAlpha',.5)
patch([x(3) x(4) x(4) x(3)], [max(y_lim) max(y_lim) 0 0], [159 203 219]/255, 'FaceAlpha',.5)
patch([x(4) x(5) x(5) x(4)], [max(y_lim) max(y_lim) 0 0], [98 141 207]/255, 'FaceAlpha',.5)
h = boxplot(t2star_remote_array,g);
set(h,'LineWidth',2); grid on;
ylim([min(t2star_remote_array), max(t2star_remote_array)]);

%% Table
VarNames = {' ', '0.3x0.3 mm^2', '0.6x0.6 mm^2',  '0.8x0.8 mm^2', '1.0x1.0 mm^2', '1.3x1.3 mm^2', '1.6x1.6 mm^2', '2.1x2.1 mm^2'};
mean_table = table({'2 mm'; '4 mm'; '6 mm'; '8 mm'},t2star_mean_reshape(:,1), t2star_mean_reshape(:,2), t2star_mean_reshape(:,3),...
    t2star_mean_reshape(:,4), t2star_mean_reshape(:,5), t2star_mean_reshape(:,6), t2star_mean_reshape(:,7), 'VariableNames',VarNames)

sd_table = table({'2 mm'; '4 mm'; '6 mm'; '8 mm'},t2star_sd_reshape(:,1), t2star_sd_reshape(:,2), t2star_sd_reshape(:,3),...
    t2star_sd_reshape(:,4), t2star_sd_reshape(:,5), t2star_sd_reshape(:,6), t2star_sd_reshape(:,7), 'VariableNames',VarNames)

t2star_table = struct;
t2star_table.mean_table = mean_table;
t2star_table.sd_table = sd_table;
save_table_f = cat(2, subject_data_dir, 'T2star_meanSD_table_', avg_name, '.mat');
save(save_table_f, 'save_table_f');

%% AHA
Segn = 50;
Groove = 0;
aha50 = struct;
for i = 1:length(whatsinit)
    [Segmentpix, stats, Mask_Segn] = AHASegmentation(whatsinit{i}, mask_struct(i).myo_mask, Segn, Groove, mask_struct(i).endo_mask);
    aha50(i).Segmentpix = Segmentpix;
    aha50(i).Mask_Segn = Mask_Segn;
end

% AHA analysis
perc_array = zeros(length(whatsinit),Segn);

for i = 1:length(whatsinit)
    for j = 1:Segn
        thresh = mean(nonzeros(whatsinit{i} .* mask_struct(i).remote_mask)) - 2*std(nonzeros(whatsinit{i} .* mask_struct(i).remote_mask));
        perc_array(i,j) = sum(aha50(i).Segmentpix{j}<thresh) / length(aha50(i).Segmentpix{j});
    end
end

res = perc_array > 0.1;
%% 50 Segments in MI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aha_mi2 = struct;
for i = 1:length(whatsinit)
    Mipix = cell(Segn, 1);
    Img = whatsinit{i};
    for j = 1:Segn
        Mipix{j,1} = Img(aha50(i).Mask_Segn .* mask_struct(i).mi_mask == j);
    end
    aha_mi2(i).Mipix = Mipix;
end

Mipix_flat = {};
for i = 1:length(whatsinit)
    idx = 1;
        for j = 1:Segn
            if ~isempty(aha_mi2(1).Mipix{j,1})
                Mipix_flat{idx} = aha_mi2(i).Mipix{j,1};
                idx = idx + 1;
            end
        end
    aha_mi2(i).Mipix_flat = Mipix_flat;
end

% AHA analysis
perc_array_mi2 = zeros(length(whatsinit),length(Mipix_flat));

for i = 1:length(whatsinit)
    for j = 1:length(Mipix_flat)
        thresh = mean(nonzeros(whatsinit{i} .* mask_struct(i).remote_mask)) - 2*std(nonzeros(whatsinit{i} .* mask_struct(i).remote_mask));
        perc_array_mi2(i,j) = sum(aha_mi2(i).Mipix_flat{j}<thresh) / length(aha_mi2(i).Mipix_flat{j});
    end
end


res_mi2 = perc_array_mi2 > 0.1;
res_mi_avg16 = aha_analysis.perc_array_mi > 0.1;

%% Confusion Matrix
figure('Position', [100 0 1600 1600]);
sens_mi2 = zeros(length(whatsinit),1);
spec_mi2 = zeros(length(whatsinit),1);
for i = 1:length(whatsinit)
    [cm, order] = confusionmat(res_mi_avg16(1,:),res_mi2(i,:));
    subplot(4,7,i);confusionchart(cm, order); 
    set(gca, 'FontSize', 18);
    sens_mi2(i) = cm(2,2) / (cm(2,2) + cm(2,1));
    spec_mi2(i) = cm(1,1) / (cm(1,1) + cm(1,2));
end

%% Plot sensitivity and specificity
sens_reshape = reshape(sens_mi2, [7, 4])';
spec_reshape = reshape(spec_mi2, [7, 4])';
figure();
subplot(2,1,1);
imagesc(sens_reshape); axis image
title('Sensitivity');
set(gca, 'FontSize', 18); colorbar;
subplot(2,1,2);
imagesc(spec_reshape); axis image
title('Specificity');
set(gca, 'FontSize', 18); colorbar;

aha_analysis2.perc_array_mi = perc_array_mi2;
aha_analysis2.sens_mi = sens_mi2;
aha_analysis2.spec_mi = spec_mi2;
aha_analysis2.aha_mi = aha_mi2;

% ROC analysis %%
auc_array_mi2 = zeros(length(whatsinit),1);
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    [X,Y,T,AUC,OPTROCPT] = perfcurve(res_mi_avg16(1,:),perc_array_mi2(i,:), 1);
    subplot(4,7,i);
    plot(X,Y, 'LineWidth', 2);
    xlabel('FPR');
    ylabel('TPR');
    text(0.5,0.5,num2str(round(AUC, 2)),'FontSize',24);
    % title('ROC for Classification')
    set(gca, 'FontSize', 16);
    auc_array_mi2(i) = AUC;
end

aha_analysis2.auc_array_mi = auc_array_mi2;

% Plot AUC %%
auc_reshape = reshape(auc_array_mi2, [7, 4])';
figure();
imagesc(auc_reshape); axis image
title('Area Under Curve');
set(gca, 'FontSize', 18); colorbar;

save(cat(2, subject_data_dir, 'aha_analysis_', avg_name, '.mat'), 'aha_analysis2');

%% Invivo Invivo parameters
disp('Invivo Starts here: ');
[list_to_read, order_to_read] = NamePicker(folder_glob);

invivo_label = 'Invivo';

% Read T2* DICOM files
whatsinit = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    whatsinit{i} = dicom23D(f);
end

% Display images
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    subplot(4,5,i);
    imagesc(whatsinit{i}); axis image;
    caxis([0 100])
end

%% Display image overlay with Mean-2SD mask
save_array = 1:1:length(whatsinit);
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
for i = 1:length(whatsinit)
    save_idx = save_array(i);
    figure();
    img2 = whatsinit{save_idx};
    mask_idx = mask_idx_array(save_idx);
    thresh = mean(nonzeros(img2 .* mask_struct(mask_idx).remote_mask)) - 2*std(nonzeros(img2 .* mask_struct(mask_idx).remote_mask));
    hemo_mask = img2 < thresh;
    subplot(1,2,1);
    imagesc(img2); caxis([0 50]); axis image; % colorbar;
    
    map = gray(256);
    rgbImage = ind2rgb(img2, map);
    rgbImage(:,:,1) = rgbImage(:,:,1) + hemo_mask.*mask_struct(mask_idx).mi_mask;
    subplot(1,2,2);
    imagesc(rgbImage); axis image;
    
    %subplot(1,2,2);
    %imagesc(img2 + hemo_mask*50.*mask_struct(mask_idx).mi_mask); caxis([0 100]); axis image; colorbar;
    
    
    mean2sd_dir = cat(2, subject_dir, 'Mean2SD_', invivo_label, '/');
    if ~exist(mean2sd_dir, 'dir')
        mkdir(mean2sd_dir)
    end
    saveas(gcf, cat(2, mean2sd_dir, num2str(save_idx), '.png'));
end

close all;
%% Otsu segmentation
num_cluster = 3;

for i = 1:length(whatsinit)
    Img = whatsinit{i};
    mask_idx = mask_idx_array(i);
    img_to_otsu = Img .* mask_struct(mask_idx).mi_mask;
    img_to_otsu(img_to_otsu == 0) = NaN;
    img_to_otsu(img_to_otsu > 50) = 50;
    
    [IDX, seg] = otsu(img_to_otsu, num_cluster);
    % Cluster1 - water; Cluster2 - myocardium; Cluster3 - hemorrhage
    
    otsu_mean = zeros(num_cluster, 1);
    for c = 1:num_cluster
        otsu_mean(c) = mean(Img(IDX == c));
    end
    
    [M,I] = min(otsu_mean);
    
    mask_struct(mask_idx).hemo_mask = IDX == I;
    
    figure();
    subplot(1,2,1);
    imagesc(Img); caxis([0 50]); axis image; %colorbar;
    
    map = gray(256);
    rgbImage = ind2rgb(Img, map);
    rgbImage(:,:,1) = rgbImage(:,:,1) + mask_struct(mask_idx).hemo_mask;
    subplot(1,2,2);
    imagesc(rgbImage); axis image;
    
    %subplot(1,2,2);
    %imagesc(Img + 50.*mask_struct(mask_idx).hemo_mask); caxis([0 100]); axis image; colorbar;
    
    otsu_dir = cat(2, subject_dir, 'Otsu_', num2str(num_cluster), '_', invivo_label, '/');
    if ~exist(otsu_dir, 'dir')
        mkdir(otsu_dir)
    end
    saveas(gcf, cat(2, otsu_dir, num2str(i), '.png'));
end

close all;

%% Mean + SD plot 
t2star_mean_array = zeros(1, length(whatsinit));
t2star_sd_array = zeros(1, length(whatsinit));
for i = 1:length(whatsinit)
    mask_idx = mask_idx_array(i);
    t2star_mean_array(i) = mean(nonzeros(mask_struct(mask_idx).remote_mask .* whatsinit{i}));
    t2star_sd_array(i) = std(nonzeros(mask_struct(mask_idx).remote_mask .* whatsinit{i}));
end

d = length(whatsinit)/4;
t2star_mean_reshape = reshape(t2star_mean_array, d, length(whatsinit)/d).';
t2star_sd_reshape = reshape(t2star_sd_array, d, length(whatsinit)/d).';

x = [0, d, 2*d, 3*d, 4*d] + [0 , 0.5, 0.5, 0.5, 1];
inplane_res = 1:d;
res = [inplane_res; inplane_res + d; inplane_res + 2*d; inplane_res + 3*d];
figure('Position', [100 0 1600 1600]);
errorbar(t2star_mean_array, t2star_sd_array, 'LineStyle', 'none' );
hold on;
ylim_lb = min(ylim); ylim_ub = max(ylim);
patch([x(1) x(2) x(2) x(1)], [max(ylim) max(ylim) 0 0], [241 194 151]/255, 'FaceAlpha',.5)
patch([x(2) x(3) x(3) x(2)], [max(ylim) max(ylim) 0 0], [199 213 161]/255, 'FaceAlpha',.5)
patch([x(3) x(4) x(4) x(3)], [max(ylim) max(ylim) 0 0], [159 203 219]/255, 'FaceAlpha',.5)
patch([x(4) x(5) x(5) x(4)], [max(ylim) max(ylim) 0 0], [98 141 207]/255, 'FaceAlpha',.5)
errorbar(res', t2star_mean_reshape', t2star_sd_reshape', '-o', 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980]); 
xticks([1.2 4 6.5 8.5 11 13.5 15.5 18 20.5 22.5 25 27.5]);
xticklabels({'0.8x0.8','---->','2.1x2.1','0.8x0.8','---->','2.1x2.1','0.8x0.8','---->','2.1x2.1','0.8x0.8','---->','2.1x2.1'})
xlim([0 x(5)]);ylim([ylim_lb, ylim_ub])

text(1,ylim_ub-2, 'Slice Thickness = 2 mm', 'FontSize', 16);
text(6.5,ylim_ub-2, 'Slice Thickness = 4 mm', 'FontSize', 16);
text(11.5,ylim_ub-2, 'Slice Thickness = 6 mm', 'FontSize', 16);
text(16.5,ylim_ub-2, 'Slice Thickness = 8 mm', 'FontSize', 16);
set(gca, 'FontSize', 16);
xlabel('Resolution (mm^2)', 'FontSize', 24); ylabel('T2^* (ms)', 'FontSize', 24);
hold off

%% Mean + SD Heatmap
figure('Position', [100 0 1600 1600]);
subplot(2,1,1);
imagesc(t2star_mean_reshape); colorbar; axis off;
for i = 1:size(t2star_mean_reshape, 1)
    for j = 1:size(t2star_mean_reshape, 2)
        text(0.9+(j-1),1+(i-1),num2str(round(t2star_mean_reshape(i,j),2)), 'Color', [0.8500, 0.3250, 0.0980], 'FontSize', 16);
    end
end
title('Mean of Remote T2* value', 'FontSize', 24);

subplot(2,1,2);
imagesc(t2star_sd_reshape); colorbar; axis off;
for i = 1:size(t2star_sd_reshape, 1)
    for j = 1:size(t2star_sd_reshape, 2)
        text(0.9+(j-1),1+(i-1),num2str(round(t2star_sd_reshape(i,j),2)), 'Color', [0.8500, 0.3250, 0.0980], 'FontSize', 16);
    end
end
title('Standard Deviation of Remote T2* value', 'FontSize', 24);

%% Mean + SD plot (boxplot)
clear t2star_remote_array
t2star_remote_array = [];
len_array = zeros(length(whatsinit), 1);
res_array = {'08', '10', '13', '16', '21'};
slc_array = {'2', '4', '6', '8'};
g = {};
for i = 1:length(whatsinit)
    mask_idx = mask_idx_array(i);
    temp = nonzeros(mask_struct(mask_idx).remote_mask .* whatsinit{i});
    len_array(i) = length(temp);
    t2star_remote_array = [t2star_remote_array; temp];
    q = fix((i-1)/5)+1;
    r = mod(i, 5);
    if r == 0
        real_r = 5 - r;
    else
        real_r = r;
    end
    g_temp = repmat({cat(2, slc_array{q}, '_', res_array{real_r})}, len_array(i), 1);
    g = [g; g_temp];
end

y_lim = [min(t2star_remote_array) max(t2star_remote_array)];
x = [0, d, 2*d, 3*d, 4*d] + [0 , 0.5, 0.5, 0.5, 1];
inplane_res = 1:d;
res = [inplane_res; inplane_res + d; inplane_res + 2*d; inplane_res + 3*d];
figure('Position', [100 0 1600 1600]);
errorbar(t2star_mean_array, t2star_sd_array, 'LineStyle', 'none' );
hold on;
ylim_lb = min(ylim); ylim_ub = max(ylim);
patch([x(1) x(2) x(2) x(1)], [max(y_lim) max(y_lim) 0 0], [241 194 151]/255, 'FaceAlpha',.5)
patch([x(2) x(3) x(3) x(2)], [max(y_lim) max(y_lim) 0 0], [199 213 161]/255, 'FaceAlpha',.5)
patch([x(3) x(4) x(4) x(3)], [max(y_lim) max(y_lim) 0 0], [159 203 219]/255, 'FaceAlpha',.5)
patch([x(4) x(5) x(5) x(4)], [max(y_lim) max(y_lim) 0 0], [98 141 207]/255, 'FaceAlpha',.5)
h = boxplot(t2star_remote_array,g);
set(h,'LineWidth',2); grid on;
ylim([min(t2star_remote_array), max(t2star_remote_array)]);

text(1.5,y_lim(2)-5, 'Slice Thickness = 2 mm', 'FontSize', 16);
text(6.5,y_lim(2)-5, 'Slice Thickness = 4 mm', 'FontSize', 16);
text(11.5,y_lim(2)-5, 'Slice Thickness = 6 mm', 'FontSize', 16);
text(16.5,y_lim(2)-5, 'Slice Thickness = 8 mm', 'FontSize', 16);
set(gca, 'FontSize', 16);
%% Table
VarNames = {' ', '0.8x0.8 mm^2', '1.0x1.0 mm^2', '1.3x1.3 mm^2', '1.6x1.6 mm^2', '2.1x2.1 mm^2'};
mean_table = table({'2 mm'; '4 mm'; '6 mm'; '8 mm'},t2star_mean_reshape(:,1), t2star_mean_reshape(:,2), t2star_mean_reshape(:,3),...
    t2star_mean_reshape(:,4), t2star_mean_reshape(:,5), 'VariableNames',VarNames)

sd_table = table({'2 mm'; '4 mm'; '6 mm'; '8 mm'},t2star_sd_reshape(:,1), t2star_sd_reshape(:,2), t2star_sd_reshape(:,3),...
    t2star_sd_reshape(:,4), t2star_sd_reshape(:,5), 'VariableNames',VarNames)

t2star_table = struct;
t2star_table.mean_table = mean_table;
t2star_table.sd_table = sd_table;
save_table_f = cat(2, subject_data_dir, 'T2star_meanSD_table_', invivo_label, '.mat');
save(save_table_f, 'save_table_f');


%% AHA
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
Segn = 50;
Groove = 0;
aha50 = struct;
for i = 1:length(whatsinit)
    mask_idx = mask_idx_array(i);
    [Segmentpix, stats, Mask_Segn] = AHASegmentation(whatsinit{i}, mask_struct(mask_idx).myo_mask, Segn, Groove, mask_struct(mask_idx).endo_mask);
    aha50(i).Segmentpix = Segmentpix;
    aha50(i).Mask_Segn = Mask_Segn;
end

% AHA analysis
perc_array = zeros(length(whatsinit),Segn);

for i = 1:length(whatsinit)
    mask_idx = mask_idx_array(i);
    for j = 1:Segn
        thresh = mean(nonzeros(whatsinit{i} .* mask_struct(mask_idx).remote_mask)) - 2*std(nonzeros(whatsinit{i} .* mask_struct(mask_idx).remote_mask));
        perc_array(i,j) = sum(aha50(i).Segmentpix{j}<thresh) / length(aha50(i).Segmentpix{j});
    end
end

res = perc_array > 0.1;
%%
figure(); imagesc(aha50(1).Mask_Segn.* mask_struct(3).mi_mask);
%% 50 Segments in MI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aha_mi2 = struct;
for i = 1:length(whatsinit)
    Mipix = cell(Segn, 1);
    Img = whatsinit{i};
    mask_idx = mask_idx_array(i);
    for j = 1:Segn
        Mipix{j,1} = Img(aha50(i).Mask_Segn .* mask_struct(mask_idx).mi_mask == j);
    end
    aha_mi2(i).Mipix = Mipix;
end

Mipix_flat = {};
for i = 1:length(whatsinit)
    idx = 1;
        for j = 1:Segn
            if ~isempty(aha_mi2(1).Mipix{j,1})
                Mipix_flat{idx} = aha_mi2(i).Mipix{j,1};
                idx = idx + 1;
            end
        end
    aha_mi2(i).Mipix_flat = Mipix_flat;
end

% AHA analysis
perc_array_mi2 = zeros(length(whatsinit),length(Mipix_flat));

for i = 1:length(whatsinit)
    mask_idx = mask_idx_array(i);
    for j = 1:length(Mipix_flat)
        thresh = mean(nonzeros(whatsinit{i} .* mask_struct(mask_idx).remote_mask)) - 2*std(nonzeros(whatsinit{i} .* mask_struct(mask_idx).remote_mask));
        perc_array_mi2(i,j) = sum(aha_mi2(i).Mipix_flat{j}<thresh) / length(aha_mi2(i).Mipix_flat{j});
    end
end


res_mi2 = perc_array_mi2 > 0.1;
res_mi_avg16 = aha_analysis.perc_array_mi > 0.1;

%% Confusion Matrix
figure('Position', [100 0 1600 1600]);
sens_mi2 = zeros(length(whatsinit),1);
spec_mi2 = zeros(length(whatsinit),1);
for i = 1:length(whatsinit)
    [cm, order] = confusionmat(res_mi_avg16(1,:),res_mi2(i,:));
    subplot(4,5,i);confusionchart(cm, order); 
    set(gca, 'FontSize', 18);
    sens_mi2(i) = cm(2,2) / (cm(2,2) + cm(2,1));
    spec_mi2(i) = cm(1,1) / (cm(1,1) + cm(1,2));
end

%% Plot sensitivity and specificity
sens_reshape = reshape(sens_mi2, [5, 4])';
spec_reshape = reshape(spec_mi2, [5, 4])';
figure();
subplot(2,1,1);
imagesc(sens_reshape); axis image
title('Sensitivity');
set(gca, 'FontSize', 18); colorbar;
subplot(2,1,2);
imagesc(spec_reshape); axis image
title('Specificity');
set(gca, 'FontSize', 18); colorbar;

aha_analysis2.perc_array_mi = perc_array_mi2;
aha_analysis2.sens_mi = sens_mi2;
aha_analysis2.spec_mi = spec_mi2;
aha_analysis2.aha_mi = aha_mi2;

% ROC analysis %%
auc_array_mi2 = zeros(length(whatsinit),1);
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    [X,Y,T,AUC,OPTROCPT] = perfcurve(res_mi_avg16(1,:),perc_array_mi2(i,:), 1);
    subplot(4,5,i);
    plot(X,Y, 'LineWidth', 2);
    xlabel('FPR');
    ylabel('TPR');
    text(0.5,0.5,num2str(round(AUC, 2)),'FontSize',24);
    % title('ROC for Classification')
    set(gca, 'FontSize', 16);
    auc_array_mi2(i) = AUC;
end

aha_analysis2.auc_array_mi = auc_array_mi2;

% Plot AUC %%
auc_reshape = reshape(auc_array_mi2, [5, 4])';
figure();
imagesc(auc_reshape); axis image
title('Area Under Curve');
set(gca, 'FontSize', 18); colorbar;

save(cat(2, subject_data_dir, 'aha_analysis_', invivo_label, '.mat'), 'aha_analysis2');
