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
dicom_fields = {'EchoTime'};
[list_to_read, order_to_read] = NamePicker(folder_glob);
whatsinit = cell(length(list_to_read), 1);
slice_data = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit{i}, slice_data{i}] = dicom23D(f, dicom_fields);
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

%% Save SNR 
SNR = struct;
SNR.snr_remote = snr_remote;
SNR.snr_air = snr_air;
save_f = cat(2, subject_data_dir, 'SNR.mat');
save(save_f, 'SNR');

%% T2* fitting
qMRinfo('mono_t2');

TE_array = zeros(length(slice_data{1,1}), 1);

for i = 1:length(TE_array)
    TE_array(i, 1) = slice_data{1,1}(i).EchoTime;
end

FitResults_struct = struct;
for i = 1:length(whatsinit)
    % Reshape data and mask
    % Reshape matrix as [Width x Height x #Slice x #TE]
    dicom_size = size(whatsinit{i});
    dicom_reshape = reshape(whatsinit{i}, dicom_size(1), dicom_size(2), 1, []);
    %
    Model = mono_t2;  % Create class from model
    %Model = Custom_OptionsGUI(Model);
    Model.Prot.SEdata.Mat = TE_array; %
    Model.st = [100 2000];
    Model.lb = [1 2000];
    Model.fx = [0 0];
    Model.voxelwise = 1;
    Model.options.FitType = 'Linear';
    data = struct;  % Create data structure
    data.SEdata = dicom_reshape;
    data.Mask = mask_struct(i).myo_mask;
    FitResults = FitData(data, Model); %fit data
    % FitResultsSave_mat(FitResults);
    FitResults_struct(i).FitResults = FitResults;
end

%% T2* map
[list_to_read, order_to_read] = NamePicker(folder_glob);
T2star_map = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    T2star_map{i} = dicom23D(f);
end

%% Compare between T2* map and fitted T2* map
figure();
subplot(2,2,1);
imagesc(T2star_map{1}.*mask_struct(1).myo_mask); caxis([0 100]);axis image;
colorbar;title('T2* Map from console');
subplot(2,2,2);
imagesc(FitResults_struct(1).FitResults.T2); caxis([0 100]);axis image;colorbar;
title('T2* map from qMRLab')
subplot(2,2,3);
diff_img = abs(T2star_map{1}.*mask_struct(1).myo_mask - FitResults_struct(1).FitResults.T2);
imagesc(diff_img);caxis([0 20]);axis image;colorbar;
title('Difference Map');
subplot(2,2,4);
imagesc(100*abs(diff_img)./(T2star_map{1}.*mask_struct(1).myo_mask));caxis([0 50]);axis image;colorbar;
title('Percentage of Difference');
%% Plot residual images
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    subplot(2,3,i);
    imagesc(FitResults_struct(i).FitResults.res);
    colorbar;caxis([0 100]);
end
%% Save FittedResults
save_f = cat(2, subject_data_dir, 'FitResults.mat');
save(save_f, 'FitResults_struct');