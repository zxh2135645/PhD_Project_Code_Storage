% T2star analysis for 20P03 Ersnest angle analysis specifically
clear all;
close all;

addpath('../function/');

%%%%%%%%%%%%%%%%%%%%%%%%%% output file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% roi
% mask_struct
% T2star_meanSD_table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% Read T2* DICOM files
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
end



%% 28 different set of parameters
%% Convert coords to masks for 28 images
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

%% Mean + SD plot (errorbar)
t2star_mean_array = zeros(1, length(whatsinit));
t2star_sd_array = zeros(1, length(whatsinit));
for i = 1:length(whatsinit)
    t2star_mean_array(i) = mean(nonzeros(mask_struct(i).remote_mask .* whatsinit{i}));
    t2star_sd_array(i) = std(nonzeros(mask_struct(i).remote_mask .* whatsinit{i}));
end

figure('Position', [100 0 1600 1600]);
errorbar(t2star_mean_array, t2star_sd_array, '-o', 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980]); 

set(gca, 'FontSize', 16);
xlabel('Flip Angel (Degree)', 'FontSize', 24); ylabel('T2^* (ms)', 'FontSize', 24);
xticks([1 2 3 4 5 6]);
xticklabels({'15','30','45','60','75','90'})
grid on;
%% Mean + SD plot (boxplot)
clear t2star_remote_array
t2star_remote_array = [];
for i = 1:length(whatsinit)
    temp = nonzeros(mask_struct(i).remote_mask .* whatsinit{i});
    t2star_remote_array = [t2star_remote_array; temp];
end
g1 = repmat({'15'},length(temp),1);
g2 = repmat({'30'},length(temp),1);
g3 = repmat({'45'},length(temp),1);
g4 = repmat({'60'},length(temp),1);
g5 = repmat({'75'},length(temp),1);
g6 = repmat({'90'},length(temp),1);
g = [g1; g2; g3; g4; g5; g6];
figure('Position', [100 0 1600 1600]);
h = boxplot(t2star_remote_array,g);
set(h,'LineWidth',2); grid on;
%% Read T2star weighted images
[list_to_read, order_to_read] = NamePicker(folder_glob);
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

%% SNR Analysis starts here
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

%% Display results
figure('Position', [100 0 400 1600]);
subplot(1,2,1);
imagesc(snr_air); axis image; colorbar;
subplot(1,2,2);
imagesc(snr_remote); axis image; colorbar;

%% Line shape
figure();
for j = 1:size(snr_air, 3)
    for i = 1:size(snr_air, 1)

        plot(snr_air(i,:,j), 'LineWidth', 2)
        hold on;
    end
    legend({'15 degree', '30 degree', '45 degree', '60 degree', '75 degree', '90 degree'});
    ylabel('SNR');
    xticks([1 2 3 4 5]);
    xticklabels({'TE1','TE2','TE3','TE4','TE5'})
    set(gca, 'FontSize', 16);
    grid on;
end

figure();
for j = 1:size(snr_remote, 3)
    for i = 1:size(snr_remote, 1)
        
        plot(snr_remote(i,:,j),  'LineWidth', 2)
        hold on;
    end
    legend({'15 degree', '30 degree', '45 degree', '60 degree', '75 degree', '90 degree'});
    ylabel('SNR');
    xticks([1 2 3 4 5]);
    xticklabels({'TE1','TE2','TE3','TE4','TE5'})
    set(gca, 'FontSize', 16);
    grid on;
end