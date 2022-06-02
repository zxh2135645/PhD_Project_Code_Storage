clear all;
close all;

%% Read Ex-vivo FF map (Manual)
figure();
n = ceil(sqrt(size(fwmc_ff, 3)));
for i = 1:size(fwmc_ff, 3)
   subplot(n,n,i);
   imagesc(fwmc_ff(:,:,i)); axis image; axis off; caxis([0 20]);
end
%%
addpath('../function/');
dicom_dir = uigetdir;
folder_glob = glob(cat(2, dicom_dir, '\*'));

labels = {'T2STAR'};
label = labels{1};

idx_array = contains(folder_glob, label);
[list_to_read, order_to_read] = NamePicker(folder_glob(idx_array));

dicom_fields = {...
    'Filename',...
    'Height', ...
    'Width', ...
    'Rows',...
    'Columns', ...
    'PixelSpacing',...
    'SliceThickness',...
    'SliceLocation',...
    'ImagePositionPatient',...
    'ImageOrientationPatient',...
    'MediaStorageSOPInstanceUID',...
    'TriggerTime',...
    'RepetitionTime',...
    };

whatsinit = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit{i} slice_data] = dicom23D(f, dicom_fields);
end
%%
fh = figure('Position', [100 100 300 400]);
axis tight manual

img = whatsinit{1};
mask_epi = zeros(size(img));
mask_endo = zeros(size(img));
mask_roi = zeros(size(img));

for i = 1:size(img, 3)
    imagesc(img(:,:,i)); axis image; axis off; colormap gray; caxis([0 50]); title(cat(2, 'Epi, slice ', num2str(i)));
    roi_epi = drawpolygon;
    mask_epi(:,:,i) = createMask(roi_epi); % mask of
    imagesc(img(:,:,i)); axis image; axis off; colormap gray; caxis([0 50]); title(cat(2, 'Endo, slice ', num2str(i)));
    roi_endo = drawpolygon;
    mask_endo(:,:,i) = createMask(roi_endo); 
end

for i = 1:size(img, 3)
    imagesc(img(:,:,i)); axis image; axis off; colormap gray; caxis([0 50]); title(cat(2, 'ROI, slice ', num2str(i)));
    roi = drawpolygon;
    mask_roi(:,:,i) = createMask(roi); 
end
%%
se = strel('disk', 1);
myo = mask_epi - mask_endo;
myo_eroded = imerode(myo, se);

figure();
for i = 1:size(img, 3)
    subplot(n,n,i);
    imagesc(myo_eroded(:,:,i).*img(:,:,i)); axis image; axis off; caxis([0 50]);
end


figure();
for i = 1:size(img, 3)
    subplot(n,n,i);
    imagesc(myo_eroded(:,:,i)); axis image; axis off;
end


mean_ff_array = zeros(size(fwmc_ff,3), 1);
std_ff_array = zeros(size(fwmc_ff,3), 1);
max_ff_array = zeros(size(fwmc_ff,3), 1);
min_ff_array = zeros(size(fwmc_ff,3), 1);
mean_roi_ff_array = zeros(size(fwmc_ff,3), 1);
std_roi_ff_array = zeros(size(fwmc_ff,3), 1);

figure();
for i = 1:size(img, 3)
    subplot(n,n,i);
    imagesc(myo_eroded(:,:,i).*fwmc_ff(:,:,i)); axis image; axis off; caxis([0 20]);
    
    fwmc_ff(fwmc_ff < 0) = 0;
    fwmc_ff(fwmc_ff > 100) = 100;
    mean_ff_array(i) = mean(nonzeros(myo_eroded(:,:,i).*fwmc_ff(:,:,i)));
    std_ff_array(i) = std(nonzeros(myo_eroded(:,:,i).*fwmc_ff(:,:,i)));
    mean_roi_ff_array(i) = mean(nonzeros(myo_eroded(:,:,i).*fwmc_ff(:,:,i).*mask_roi(:,:,i)));
    std_roi_ff_array(i) = std(nonzeros(myo_eroded(:,:,i).*fwmc_ff(:,:,i).*mask_roi(:,:,i)));
    
    if ~isempty(nonzeros(myo_eroded(:,:,i).*fwmc_ff(:,:,i)))
        max_ff_array(i) = max(nonzeros(myo_eroded(:,:,i).*fwmc_ff(:,:,i)));
        min_ff_array(i) = min(nonzeros(myo_eroded(:,:,i).*fwmc_ff(:,:,i)));
    else
        max_ff_array(i) = nan;
        min_ff_array(i) = nan;
    end
end

%% Manually load masks
figure();
for i = 1:size(vol_img_3D_te8,3)
    subplot(2,2,i); imagesc(vol_img_3D_te8(:,:,i).*myoRefMask_3D(:,:,i)); axis image;
end

figure();
for i = 1:size(vol_img_3D_te8,3)
    subplot(2,2,i); imagesc(vol_img_3D_te8(:,:,i)); axis image;
end

%% 
base_dir = uigetdir;
name = '18D16';
time_point = '6MO';
ff_map = cell(1, length(glob_names));
for f = 1:length(ff_map)
    ff_map{f} = load(cat(2,  base_dir, '/FF_Data/',  name, '/', name, '_', time_point, '_', glob_names{f}, '.mat'), 'fwmc_ff');
end

% convert ff_map to matrix
ff = zeros(size(ff_map{1}.fwmc_ff,1), size(ff_map{1}.fwmc_ff, 2), length(ff_map));
for f = 1:length(ff_map)
    ff(:,:,f) = ff_map{f}.fwmc_ff;
end

mean_ff_array = zeros(size(ff,3), 1);
std_ff_array = zeros(size(ff,3), 1);
max_ff_array = zeros(size(ff,3), 1);
min_ff_array = zeros(size(ff,3), 1);

figure();
for i = 1:size(ff, 3)
    subplot(2,2,i);
    %imagesc(myoRefMask_3D(:,:,i).*ff(:,:,i)); axis image; axis off; caxis([0 20]);
    imagesc(ff(:,:,i)); axis image; axis off; caxis([0 20]);
    
    ff(ff < 0) = 0;
    ff(ff > 100) = 100;
    mean_ff_array(i) = mean(nonzeros(myoRefMask_3D(:,:,i).*ff(:,:,i)));
    std_ff_array(i) = std(nonzeros(myoRefMask_3D(:,:,i).*ff(:,:,i)));
    if ~isempty(nonzeros(myoRefMask_3D(:,:,i).*ff(:,:,i)))
        max_ff_array(i) = max(nonzeros(myoRefMask_3D(:,:,i).*ff(:,:,i)));
        min_ff_array(i) = min(nonzeros(myoRefMask_3D(:,:,i).*ff(:,:,i)));
    end
end
