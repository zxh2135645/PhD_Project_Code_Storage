% To pull up images for presentation use
% colormapped original image with epicardium masks
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

%% Display original colormap image
% Cropped center of the endo_mask 03/21/2021
save_array = 1:1:length(whatsinit);
len = 36;
base_x = size(mask_struct(7).myo_mask, 1);
for i = 1:length(whatsinit)
    save_idx = save_array(i);
    figure();
    img2 = whatsinit{save_idx};
    img2 = img2 .* mask_struct(i).epi_mask;
    imagesc(img2); caxis([0 100]); axis image; %colormap default; %colorbar;
    axis off;
    colormap(brewermap([],'RdBu'));
    c = colorbar;
    w = c.FontSize;
    c.FontSize = 20;
    % '*YlGnBu'
        
    colormap_dir = cat(2, subject_dir, 'Colormap_', avg_name, '/');
    if ~exist(colormap_dir, 'dir')
        mkdir(colormap_dir)
    end
    saveas(gcf, cat(2, colormap_dir, num2str(save_idx), '.png'));
    
    stats = regionprops(mask_struct(i).myo_mask);
    centroid = round(stats.Centroid);
    
    ioi_x = size(mask_struct(i).myo_mask, 1);
    multiplier = ioi_x / base_x;
    len_mul = multiplier * len;
    img2_cropped = imcrop(img2, [centroid(1) - len_mul/2, centroid(2) - len_mul/2, len_mul, len_mul]);
    figure(); 
    imagesc(img2_cropped);
    caxis([0 100]); axis image;
    axis off;
    colormap(brewermap([],'RdBu'));
    saveas(gcf, cat(2, colormap_dir, num2str(save_idx), '_cropped.png'));
end

close all;


%% Get intensity from a pixel
n = 28;
figure('Position', [100 0 1600 1600]);
img2 = whatsinit{n};
img2 = img2 .* mask_struct(n).epi_mask;
hIm = imagesc(img2); caxis([0 100]); axis image;
colormap(brewermap([],'RdBu'));
roi = drawpolygon(gca);

roi_mask = createMask(roi);

mean_v = mean(nonzeros(img2 .* roi_mask));
% 20P10_Exvivo7 
% T2* is 13.08/15.09 ms at core @ 0.3x0.3x2
% For hemo- MI region, T2* = 59.73
% @ 2.1x2.1x8
% T2* is 20.50 ms and 52.50

% 18P93
% T2* is 27.50/24.67/26.40 ms @ 0.3x0.3x2
% T2* is 37 ms @ 2.1x2.1x8

%% Measure Distance of thinned wall and hemo+ (Optional)
figure('Position', [100 0 1600 1600]);
img2 = whatsinit{1};
img2 = img2 .* mask_struct(1).epi_mask;
hIm = imagesc(img2); caxis([0 100]); axis image; %colormap default; %colorbar;
sz = size(img2);
colormap(brewermap([],'RdBu'));

myData = MeasureDist_Image(hIm, sz);

% 18P93
% (24.03 + 16.00 + 19.61 + 14.65 + 27.98 + 20.43 + 16.80 + 21.37 + 18.30 + 28.36)/10 * (200/sz(1))
% 5.6394 mm 
% MI thinnest is ~0.3 mm 
% Thickest is ~ 1.36 mm
% From transmurality -> mean thickness is 0.31 mm 
% 20P10_Exvivo7
% (20.63 + 12.58 + 21.00 + 21.57 + 23.76 + 25.54 + 19.82 + 25.16) / 8 * (200/sz(1))
% 5.5358 mm
% MI thinnest is ~ 0.6-0.7 mm
% MI thickest is ~ 3.26 mm
% From transmurality -> mean thickness is 1.1145 