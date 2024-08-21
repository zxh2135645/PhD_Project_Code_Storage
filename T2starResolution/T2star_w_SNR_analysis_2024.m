clear all;
close all;

% the main body for T2* SNR analysis

addpath('../function/');
%%%%%%%%%%%%%%%%%%%%%%%%%% input file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% roi 
% mask_struct
% Both from T2star_analysis_main.m
%%%%%%%%%%%%%%%%%%%%%%%%%% output file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SNR - struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% {'18P90', '18P93', '20P40', '20P10_Exvivo7', '20P11_Exvivo6', '18P92', '18P94_Exvivo3', '18P95', '17P73', '20P48'};
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

avg_num = input('Please type average number here:  ');
if isnumeric(avg_num)
    avg_name = cat(2, 'Avg', num2str(avg_num, '%04.f'));
else
    avg_name = avg_num;
end
%% Read T2* weighted DICOM files
whatsinit = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    whatsinit{i} = dicom23D(f);
end


%% Display images
figure('Position', [100 0 1600 1600]);
row = 4;
col = length(whatsinit) / row;
for i = 1:length(whatsinit)
    subplot(row,col,i);
    imagesc(whatsinit{i}(:,:,5)); axis image;
    %caxis([0 100])
end

% For scheme figures
% figure('Position', [100 0 1600 1600]);
% row = 2;
% col = 3;
% for i = 1:size(whatsinit{1}, 3)
%     subplot(row,col,i);
%     imagesc(whatsinit{1}); axis image;
%     caxis([0 100])
%     colormap gray;
% end
%% Draw contours @ epi, endo, MI, remote, fluid
% Coords and Masks should be generated already in data folder
% But I'm still not assuming that
% images size is 3-dimensional
img = whatsinit{1};
myo_coords_cell = cell(size(img, 4), 2);
roi_save = cat(2, subject_data_dir, 'roi.mat');

if ~exist(roi_save, 'file')
for i = 1:(size(img, 4))
    disp('Epicardium: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,1,i)); axis image;
    %caxis([0 100])
    epi = drawpolygon(gca);
    epi_coords = epi.Position;
    
    disp('Endocardium: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,1,i)); axis image;
    %caxis([0 100])
    endo = drawpolygon(gca);
    endo_coords = endo.Position;
    
    disp('MI roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,1,i)); axis image;
    %caxis([0 100])
    mi = drawpolygon(gca);
    mi_coords = mi.Position;
    
    disp('remote roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,1,i)); axis image;
    %caxis([0 100])
    remote = drawpolygon(gca);
    remote_coords = remote.Position;
    
    disp('fluid roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,1,i)); axis image;
    %caxis([0 100])
    fluid = drawpolygon(gca);
    fluid_coords = fluid.Position;
    
    disp('Center line roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,1,i)); axis image;
    %caxis([0 100])
    center_line = drawpolygon(gca);
    center_coords = center_line.Position;
    
    disp('air roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,1,i)); axis image;
    %caxis([0 100])
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
end

% Convert coords to masks for 28 images
img_size_truth = size(img);
mask_save = cat(2, subject_data_dir, 'mask.mat');


if ~exist(mask_save, 'file')
    figure();
    mask_struct = struct;
    for i = 1:length(whatsinit)
        img2 = whatsinit{i};
        img2_size = size(img2);
        ratio = img_size_truth(1) ./ img2_size(1);
        imagesc(img2(:,:,1)); % caxis([0 100]);
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

%% SNR Analysis starts here
% 
img_size = size(whatsinit{1});
snr_remote = zeros(length(whatsinit), size(whatsinit{1},3));
snr_air = zeros(length(whatsinit), size(whatsinit{1},3));
sig_remote_mean = zeros(length(whatsinit), size(whatsinit{1},3));
sig_remote_sd = zeros(length(whatsinit), size(whatsinit{1},3));
sig_air_mean = zeros(length(whatsinit), size(whatsinit{1},3));
sig_air_sd = zeros(length(whatsinit), size(whatsinit{1},3));
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];

for i = 1:length(whatsinit)
    img = whatsinit{i};
    
    for j = 1:size(whatsinit{1},3)
        
        if ~strcmp('Invivo', avg_name)
            idx = find(mask_struct(i).remote_mask == 1);
            remote = mask_struct(i).remote_mask .* img(:,:,j);
            idx_air = find(mask_struct(i).air_mask == 1);
            air = mask_struct(i).air_mask .* img(:,:,j);

            thresh = mean(nonzeros(remote)) - 2*std(nonzeros(remote));
            hemo_mask = (img(:,:,j) < thresh) .* mask_struct(i).mi_mask .* mask_struct(i).myo_mask;
        else
            mask_idx = mask_idx_array(i);
            idx = find(mask_struct(mask_idx).remote_mask == 1);
            remote = mask_struct(mask_idx).remote_mask .* img(:,:,j);
            idx_air = find(mask_struct(mask_idx).air_mask == 1);
            air = mask_struct(mask_idx).air_mask .* img(:,:,j);

            thresh = mean(nonzeros(remote)) - 2*std(nonzeros(remote));
            hemo_mask = (img(:,:,j) < thresh) .* mask_struct(mask_idx).mi_mask .* mask_struct(mask_idx).myo_mask;
        end
        
        sig_remote_mean(i,j) = mean(remote(idx));
        sig_remote_sd(i,j) = std(remote(idx));
        sig_air_mean(i,j) = mean(air(idx_air));
        sig_air_sd(i,j) = std(air(idx_air));
        sig_hemo_mean(i,j) = mean(nonzeros(hemo_mask.*img(:,:,j)));
        snr_remote(i,j) = mean(remote(idx)) / std(remote(idx));
        snr_air(i,j) = mean(remote(idx)) / std(air(idx_air));
        cnr_remote(i,j) = (sig_remote_mean(i,j)-sig_hemo_mean(i,j)) ./ std(remote(idx));
    end
end

%% Read CNR
label = labels{5};
idx_array = contains(folder_glob, label);
[list_to_read, order_to_read] = NamePicker(folder_glob);
groundtruth = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    groundtruth{i} = dicom23D(f);
end

cnr_remote = zeros(length(whatsinit), size(whatsinit{1},3));
remote_gt = (mask_struct(1).remote_mask) .* groundtruth{1};
thresh = mean(nonzeros(remote_gt)) - 2*std(nonzeros(remote_gt));
hemo_mask_gt = (groundtruth{1} < thresh) .* mask_struct(1).mi_mask .* mask_struct(1).myo_mask;

for i = 1:length(whatsinit)
    img = whatsinit{i};
    img_size = size(img);
    hemo_mask_resized = imresize(hemo_mask_gt, img_size)>0.25;
    for j = 1:size(whatsinit{1},3)
        
        if ~strcmp('Invivo', avg_name)
            idx = find(mask_struct(i).remote_mask == 1);
            remote = mask_struct(i).remote_mask .* img(:,:,j);
            
        else
            mask_idx = mask_idx_array(i);
            idx = find(mask_struct(mask_idx).remote_mask == 1);
            remote = mask_struct(mask_idx).remote_mask .* img(:,:,j);      
        end
        
        sig_remote_mean(i,j) = mean(remote(idx));
        sig_remote_sd(i,j) = std(remote(idx));
        sig_hemo_mean(i,j) = mean(nonzeros(hemo_mask_resized.*img(:,:,j)));
        snr_remote(i,j) = mean(remote(idx)) / std(remote(idx));
        cnr_remote(i,j) = (sig_remote_mean(i,j)-sig_hemo_mean(i,j)) ./ std(remote(idx));
    end
end

%% SNR average (Different color scheme)
mean_snr_remote_invivo_1d = [3.889101759, 6.796342678, 8.90610959, 11.69427981, 15.14968152,...
    7.989617412, 12.37448587, 15.04384232, 18.1386487, 24.18634346, 11.57426404, 16.37096552,...
    19.903348, 21.82762437, 27.80253859, 14.41929372, 18.08539358, 21.71116137, 24.08643505, 28.5705392];
sd_snr_remote_invivo_1d = [1.321445202, 0.910185045, 1.487865209, 1.737043979, 2.000581889,...
    1.303934476, 1.654983191, 2.516943858, 3.099431548, 3.655129308, 1.921314911, 1.944751457,...
    3.142859944, 3.135530109, 8.800968977, 2.160788375, 2.878735521, 4.896095446, 6.407345975, 8.745891485];

d = 7;
x = [0, d, 2*d, 3*d, 4*d] + [0 , 0.5, 0.5, 0.5, 1];
inplane_res = 1:d;
res = [inplane_res; inplane_res + d; inplane_res + 2*d; inplane_res + 3*d];
figure('Position', [100 0 800 400]);
hax = axes;
plotHandles = zeros(4,2);
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
color_cell = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]};

avg_temp = mean_snr_remote_invivo_1d;
avg_sd_temp = sd_snr_remote_invivo_1d;
avg_x = reshape(mask_idx_array, [], 4).';
avg_compare = reshape(avg_temp, [], 4).';
avg_err = reshape(avg_sd_temp, [], 4).';

hold on;
ylim([0 40]);
ylim_lb = 0; ylim_ub = max(ylim);

patch([x(1) x(2) x(2) x(1)], [max(ylim) max(ylim) 0 0], [247 247 247]/255, 'FaceAlpha',.5)
patch([x(2) x(3) x(3) x(2)], [max(ylim) max(ylim) 0 0], [204 204 204]/255, 'FaceAlpha',.5)
patch([x(3) x(4) x(4) x(3)], [max(ylim) max(ylim) 0 0], [150 150 150]/255, 'FaceAlpha',.5)
patch([x(4) x(5) x(5) x(4)], [max(ylim) max(ylim) 0 0], [99 99 99]/255, 'FaceAlpha',.5)
plotHandles(:,2) = plot(avg_x.', avg_compare.', 'LineWidth', 2, 'Color', color_cell{2}); hold on;
plotHandles(:,1) = errorbar(avg_x.', avg_compare.', avg_err.', 'LineWidth', 2, 'Color', color_cell{2});

set(gca, 'FontSize', 18);
%grid on;
%xticks([1.2 4 6.5 8.5 11 13.5 15.5 18 20.5 22.5 25 27.5]);
%xticklabels({'0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1'})
xlim([0 x(5)]);ylim([ylim_lb, ylim_ub])

color_cell_invivo = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
% color_cell_invivo = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};
set(plotHandles(:,1), 'LineStyle', 'none', 'Marker', 's', 'Color', color_cell_invivo{4})
set(plotHandles(:,1), 'LineWidth', 1, 'Marker', 's', 'MarkerSize', 8, ...
    'MarkerEdgeColor', color_cell_invivo{5}, 'MarkerFaceColor' , color_cell_invivo{2})
set(plotHandles(:,2), 'LineWidth', 1, 'Color', color_cell_invivo{4});

set(gca, 'Xcolor', 'w', 'Ycolor', 'w')
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'XTickLabels', []);
%hax.YAxis(1).Visible='off';
set(gca, 'YTickLabels', []);
set(gca,'LineWidth', 1.5,'TickLength',[0.02 0.02]);
set(gca,'TickDir','out');
%% Display results
figure('Position', [100 0 400 1600]);
subplot(1,2,1);
imagesc(snr_air); axis image; colorbar;
subplot(1,2,2);
imagesc(snr_remote); axis image; colorbar;

%%
row = 4;
col = length(whatsinit) / row;
snr_air_max = round(max(snr_air(:)),-2);
snr_remote_max = round(max(snr_remote(:)),-1);
snr_air_reshape = permute(reshape(snr_air, col, row, []), [2,1,3]);
snr_remote_reshape = permute(reshape(snr_remote, col, row, []), [2,1,3]);
figure('Position', [100 0 800 1600]);
subplot(5,2,1);
imagesc(snr_air_reshape(:,:,1)); axis image; caxis([0 snr_air_max]); colorbar;
subplot(5,2,2);
imagesc(snr_remote_reshape(:,:,1)); axis image; caxis([0 snr_remote_max]);colorbar;
subplot(5,2,3);
imagesc(snr_air_reshape(:,:,2)); axis image; caxis([0 snr_air_max]);colorbar;
subplot(5,2,4);
imagesc(snr_remote_reshape(:,:,2)); axis image; caxis([0 snr_remote_max]);colorbar;
subplot(5,2,5);
imagesc(snr_air_reshape(:,:,3)); axis image; caxis([0 snr_air_max]);colorbar;
subplot(5,2,6);
imagesc(snr_remote_reshape(:,:,3)); axis image; caxis([0 snr_remote_max]);colorbar;
subplot(5,2,7);
imagesc(snr_air_reshape(:,:,4)); axis image; caxis([0 snr_air_max]);colorbar;
subplot(5,2,8);
imagesc(snr_remote_reshape(:,:,4)); axis image; caxis([0 snr_remote_max]);colorbar;
subplot(5,2,9);
imagesc(snr_air_reshape(:,:,5)); axis image; caxis([0 snr_air_max]);colorbar;
subplot(5,2,10);
imagesc(snr_remote_reshape(:,:,5)); axis image; caxis([0 snr_remote_max]);colorbar;

%% Save SNR  (Didn't save)
SNR = struct;
SNR.snr_remote = snr_remote;
SNR.snr_air = snr_air;
SNR.sig.sig_remote_mean = sig_remote_mean;
SNR.sig.sig_air_mean = sig_air_mean;
SNR.sig.sig_remote_sd = sig_remote_sd;
SNR.sig.sig_air_sd = sig_air_sd;

save_f = cat(2, subject_data_dir, 'SNR_weighted_', avg_name, '.mat');
save(save_f, 'SNR');

%% R^2 analysis of fitting residual might be a good analysis?
