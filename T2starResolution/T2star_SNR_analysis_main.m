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
avg_name = cat(2, 'Avg', num2str(avg_num, '%04.f'));
%% Read T2* weighted DICOM files
whatsinit = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    whatsinit{i} = dicom23D(f);
end

% Display images
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    subplot(4,7,i);
    imagesc(whatsinit{i}(:,:,5)); axis image;
    %caxis([0 100])
end

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
snr_remote = zeros(length(whatsinit), img_size(3));
snr_air = zeros(length(whatsinit), img_size(3));

for i = 1:length(whatsinit)
    img = whatsinit{i};
    
    for j = 1:img_size(3)
        idx = find(mask_struct(i).remote_mask == 1);
        
        remote = mask_struct(i).remote_mask .* img(:,:,j);
        snr_remote(i,j) = mean(remote(idx)) / std(remote(idx));
        
        idx_air = find(mask_struct(i).air_mask == 1);
        air = mask_struct(i).air_mask .* img(:,:,j);
        snr_air(i,j) = mean(remote(idx)) / std(air(idx_air));
    end
end
%% DIsplay results
figure('Position', [100 0 400 1600]);
subplot(1,2,1);
imagesc(snr_air); axis image; colorbar;
subplot(1,2,2);
imagesc(snr_remote); axis image; colorbar;

%%
snr_air_max = round(max(snr_air(:)),-2);
snr_remote_max = round(max(snr_remote(:)),-2);
snr_air_reshape = permute(reshape(snr_air, 7, 4, []), [2,1,3]);
snr_remote_reshape = permute(reshape(snr_remote, 7, 4, []), [2,1,3]);
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
%% Line shape
figure();
for j = 1:size(snr_air_reshape, 3)
    subplot(3,2,j);
    for i = 1:size(snr_air_reshape, 1)

        plot(snr_air_reshape(i,:,j))
        hold on;
    end
    
end

figure();
for j = 1:size(snr_remote_reshape, 3)
    subplot(3,2,j);
    for i = 1:size(snr_remote_reshape, 1)

        plot(snr_remote_reshape(i,:,j))
        hold on;
    end
end

%% Save SNR 
SNR = struct;
SNR.snr_remote = snr_remote;
SNR.snr_air = snr_air;
save_f = cat(2, subject_data_dir, 'SNR_', avg_name, '.mat');
save(save_f, 'SNR');

%% R^2 analysis of fitting residual might be a good analysis?
