clear all;
close all;
%% Read DICOM file (merry)
addpath('../function/');

base_dir = uigetdir;
folder_glob = glob(cat(2, base_dir, '\*'));

labels = {'T1MAP'};

label = labels{1};
idx_array = contains(folder_glob, label);

[list_to_read, order_to_read] = NamePicker(folder_glob(idx_array));

whatsinit = cell(length(list_to_read), 1);
slice_data = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit{i}, slice_data{i}] = dicom23D(f);
end

%% T1 Map
figure();
t1map_cropped = imcrop(whatsinit{1}(:,:,3), [20 40 45 60]);
imagesc(t1map_cropped); caxis([800 1800]); axis image;
colorbar;
colormap(brewermap([],'*RdYlBu'))

%% T1 MOCO
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
    'InversionTime',...
    };

whatsinit = cell(length(list_to_read), 1);
slice_data = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit{i}, slice_data{i}] = dicom23D(f, dicom_fields);
end
%%
l = 17:24;
t1_reshaped = whatsinit{1};
t1moco_3d = t1_reshaped(:,:,l);
t1moco_3d_norm = t1moco_3d ./ max(t1moco_3d(:));

for i = 1:length(l)
    figure()
    t1moco_cropped = imcrop(t1moco_3d_norm(:,:,i), [20 40 45 60]);
    imagesc(t1moco_cropped); axis image;
    colormap gray;
    axis off;
end

figure();
implay(t1moco_3d_norm);
%%
ff_cropped = imcrop(fwmc_ff(:,:,4), [20 40 45 60]);
figure(); imagesc(ff_cropped); caxis([0 40]);
axis image; colorbar;
colormap(brewermap([],'*RdYlBu'));

r2star_cropped = imcrop(fwmc_r2star(:,:,4), [20 40 45 60]);
figure(); imagesc(r2star_cropped); caxis([0 150]);
axis image; colorbar;
colormap(brewermap([],'*RdYlBu'));
%%
inv_t = [172, 252, 1172, 1252, 2172, 2252, 3172, 4172];
t1moco_3d = t1_reshaped(:,:,l);
y_1 =  squeeze(t1moco_3d(33,17,:)); % 554 ms , %33.10, 29.67
y_2 =  squeeze(t1moco_3d(30,15,:)); % 1771 ms  %15.33, 31.93
y_3 =  squeeze(t1moco_3d(27,13,:)); % 1231 ms  %24.51,  4.30

figure();
plot(inv_t, y_1, 'LineWidth', 1.5);
hold on;
plot(inv_t, y_2, 'LineWidth', 1.5);
plot(inv_t, y_3, 'LineWidth', 1.5);
%% Load contours
load(GetFullPath(cat(2, pwd, '/../../T1_Fat_Project/ContourData_old03012021/Merry/Merry_Exvivo/T1_CMR/freeROI/freeROI.mat')))
roi_cropped = imcrop(squeeze(freeROIMask_3D), [20 40 45 60]);

figure();
subplot(2,2,1);
imagesc(t1map_cropped .* roi_cropped); axis image; colorbar;
colormap(brewermap([],'*RdYlBu'));
subplot(2,2,2);
imagesc(ff_cropped .* roi_cropped);axis image; colorbar;
colormap(brewermap([],'*RdYlBu'));
subplot(2,2,3);
imagesc(r2star_cropped .* roi_cropped);axis image; colorbar;
colormap(brewermap([],'*RdYlBu'));

%% Version 1 (no pixel shift)
se = strel('disk',1);
t1map_temp = t1map_cropped .* roi_cropped;
t1map_temp_eroded = imerode(t1map_temp, se);
nonzero_mask = ones(size(t1map_temp));
nonzero_mask(t1map_temp_eroded == 0) = nan;

t1_array = nonzero_mask.*t1map_cropped;
t1_array = t1_array(~isnan(t1_array));

ff_array = nonzero_mask.*ff_cropped;
ff_array = ff_array(~isnan(ff_array));
ff_array(ff_array<0) = 0;

figure();
subplot(1,2,1);
imagesc(nonzero_mask.*t1map_cropped); caxis([800 1800]); axis image;
colorbar;
colormap(brewermap([],'*RdYlBu'));
subplot(1,2,2);
scatter(ff_array, t1_array, 'LineWidth', 1.5);

%% Version 2 (no pixel shift)
se = strel('disk',1);
t1map_temp = t1map_cropped .* roi_cropped;
t1map_temp_eroded = imerode(t1map_temp, se);
ind = find(t1map_temp_eroded);
sz = size(t1map_temp_eroded);
[row, col] = ind2sub(sz, ind)

row_shift = row + 1;
col_shift = col + 1;
mask_shift = zeros(sz);
ind_shift = sub2ind(sz, row_shift, col_shift);
mask_shift(ind_shift) = 1;

nonzero_mask = ones(size(t1map_temp_eroded));
nonzero_mask_shift = ones(size(t1map_temp_eroded));
nonzero_mask(t1map_temp_eroded == 0) = nan;
nonzero_mask_shift(mask_shift == 0) = nan;

t1_array = nonzero_mask.*t1map_cropped;
t1_array = t1_array(~isnan(t1_array));

ff_array = nonzero_mask_shift.*ff_cropped;
ff_array = ff_array(~isnan(ff_array));
ff_array(ff_array<0) = 0;

r2star_array = nonzero_mask_shift.*r2star_cropped;
r2star_array = r2star_array(~isnan(r2star_array));
r2star_array(r2star_array<0) = 0;

figure();
subplot(1,2,1);
imagesc(nonzero_mask.*t1map_cropped); caxis([800 1800]); axis image;
colorbar;
colormap(brewermap([],'*RdYlBu'));
subplot(2,2,2);
scatter(ff_array, t1_array, 'LineWidth', 1.5);
xlabel('FF (%)'); ylabel('T1 (ms)');
subplot(2,2,4);
scatter(r2star_array, t1_array, 'LineWidth', 1.5);
xlabel('r2star (s-1)'); ylabel('T1 (ms)');

%% Load contours (Remote)
load(GetFullPath(cat(2, pwd, '/../../T1_Fat_Project/ContourData_old03012021/Merry/Merry_Exvivo/T1_CMR/MyoReference/myoRef.mat')))
ref_cropped = imcrop(squeeze(myoRefMask_3D), [20 40 45 60]);

se = strel('disk',1);
t1map_temp = t1map_cropped .* ref_cropped;
t1map_temp_eroded = imerode(t1map_temp, se);
ind = find(t1map_temp_eroded);
sz = size(t1map_temp_eroded);
[row, col] = ind2sub(sz, ind)

row_shift = row + 1;
col_shift = col + 1;
mask_shift = zeros(sz);
ind_shift = sub2ind(sz, row_shift, col_shift);
mask_shift(ind_shift) = 1;

nonzero_mask = ones(size(t1map_temp_eroded));
nonzero_mask_shift = ones(size(t1map_temp_eroded));
nonzero_mask(t1map_temp_eroded == 0) = nan;
nonzero_mask_shift(mask_shift == 0) = nan;

t1_array = nonzero_mask.*t1map_cropped;
t1_array = t1_array(~isnan(t1_array));

ff_array = nonzero_mask_shift.*ff_cropped;
ff_array = ff_array(~isnan(ff_array));
ff_array(ff_array<0) = 0;

r2star_array = nonzero_mask_shift.*r2star_cropped;
r2star_array = r2star_array(~isnan(r2star_array));
r2star_array(r2star_array<0) = 0;

figure();
subplot(1,2,1);
imagesc(nonzero_mask.*t1map_cropped); caxis([800 1800]); axis image;
colorbar;
colormap(brewermap([],'*RdYlBu'));
subplot(2,2,2);
scatter(ff_array, t1_array, 'LineWidth', 1.5);
xlabel('FF (%)'); ylabel('T1 (ms)');
subplot(2,2,4);
scatter(r2star_array, t1_array, 'LineWidth', 1.5);
xlabel('r2star (s-1)'); ylabel('T1 (ms)');

%% ------------------------------------------------------------------------
%  ************************************************************************
%  ------------------------------------------------------------------------
%% Read DICOM file (Ryn)
addpath('../function/');

base_dir = uigetdir;
folder_glob = glob(cat(2, base_dir, '\*'));

labels = {'T1MAP'};

label = labels{1};
idx_array = contains(folder_glob, label);

[list_to_read, order_to_read] = NamePicker(folder_glob(idx_array));

whatsinit = cell(length(list_to_read), 1);
slice_data = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit{i}, slice_data{i}] = dicom23D(f);
end

%% T1 Map
figure();
t1map_cropped = imcrop(whatsinit{1}(:,:,11), [20 38 56 60]);
imagesc(t1map_cropped); caxis([800 1800]); axis image;
colorbar;
colormap(brewermap([],'*RdYlBu'));

%% T1 MOCO
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
    'InversionTime',...
    };

whatsinit = cell(length(list_to_read), 1);
slice_data = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit{i}, slice_data{i}] = dicom23D(f, dicom_fields);
end

%%
l = 81:88;
t1_reshaped = whatsinit{1};
t1moco_3d = t1_reshaped(:,:,l);
t1moco_3d_norm = t1moco_3d ./ max(t1moco_3d(:));

for i = 1:length(l)
    figure();
    t1moco_cropped = imcrop(t1moco_3d_norm(:,:,i), [20 38 56 60]);
    imagesc(t1moco_cropped); axis image;
    colormap gray;
    axis off;
end

figure();
implay(t1moco_3d_norm);

%%
fwmc_ff = (abs(fat) ./ (abs(fat) + abs(water)))*100;
fwmc_r2star = R2s;

% ff_cropped = imcrop(fwmc_ff(:,:,3), [20 38 56 60]);
% figure(); imagesc(ff_cropped); caxis([0 40]);
% axis image; colorbar;
% colormap(brewermap([],'*RdYlBu'));
% 
% r2star_cropped = imcrop(fwmc_r2star(:,:,3), [20 38 56 60]);
% figure(); imagesc(r2star_cropped); caxis([0 150]);
% axis image; colorbar;
% colormap(brewermap([],'*RdYlBu'));

figure(); imagesc(fwmc_ff(:,:,42)); caxis([0 40]);
axis image; colorbar;
colormap(brewermap([],'*RdYlBu'));

figure(); imagesc(fwmc_r2star(:,:,42)); caxis([0 150]);
axis image; colorbar;
colormap(brewermap([],'*RdYlBu'));
%%
inv_t = [172, 252, 1172, 1252, 2172, 2252, 3172, 4172];
t1moco_3d = t1_reshaped(:,:,l);
y_1 =  squeeze(t1moco_3d(33,17,:)); % 554 ms , %33.10, 29.67
y_2 =  squeeze(t1moco_3d(30,15,:)); % 1771 ms  %15.33, 31.93
y_3 =  squeeze(t1moco_3d(27,13,:)); % 1231 ms  %24.51,  4.30

figure();
plot(inv_t, y_1, 'LineWidth', 1.5);
hold on;
plot(inv_t, y_2, 'LineWidth', 1.5);
plot(inv_t, y_3, 'LineWidth', 1.5);

%% Load contours
load(GetFullPath(cat(2, pwd, '/../../T1_Fat_Project/ContourData_old03012021/Ryn/Ryn_Exvivo/T1_CMR/freeROI/freeROI.mat')))
roi_cropped = imcrop(squeeze(freeROIMask_3D(:,:,11)), [20 38 56 60]);

figure();
subplot(2,2,1);
imagesc(t1map_cropped .* roi_cropped); axis image; colorbar;
colormap(brewermap([],'*RdYlBu'));
subplot(2,2,2);
imagesc(ff_cropped .* roi_cropped);axis image; colorbar;
colormap(brewermap([],'*RdYlBu'));
subplot(2,2,3);
imagesc(r2star_cropped .* roi_cropped);axis image; colorbar;
colormap(brewermap([],'*RdYlBu'));

%% Version 1 (no pixel shift)
se = strel('disk',1);
t1map_temp = t1map_cropped .* roi_cropped;
t1map_temp_eroded = imerode(t1map_temp, se);
nonzero_mask = ones(size(t1map_temp));
nonzero_mask(t1map_temp_eroded == 0) = nan;

t1_array = nonzero_mask.*t1map_cropped;
t1_array = t1_array(~isnan(t1_array));

ff_array = nonzero_mask.*ff_cropped;
ff_array = ff_array(~isnan(ff_array));
ff_array(ff_array<0) = 0;

figure();
subplot(1,2,1);
imagesc(nonzero_mask.*t1map_cropped); caxis([800 1800]); axis image;
colorbar;
colormap(brewermap([],'*RdYlBu'));
subplot(1,2,2);
scatter(ff_array, t1_array, 'LineWidth', 1.5);

%% Version 2 (no pixel shift)
se = strel('disk',1);
t1map_temp = t1map_cropped .* roi_cropped;
t1map_temp_eroded = imerode(t1map_temp, se);
ind = find(t1map_temp_eroded);
sz = size(t1map_temp_eroded);
[row, col] = ind2sub(sz, ind)

row_shift = row + 1;
col_shift = col + 1;
mask_shift = zeros(sz);
ind_shift = sub2ind(sz, row_shift, col_shift);
mask_shift(ind_shift) = 1;

nonzero_mask = ones(size(t1map_temp_eroded));
nonzero_mask_shift = ones(size(t1map_temp_eroded));
nonzero_mask(t1map_temp_eroded == 0) = nan;
nonzero_mask_shift(mask_shift == 0) = nan;

t1_array = nonzero_mask.*t1map_cropped;
t1_array = t1_array(~isnan(t1_array));

ff_array = nonzero_mask_shift.*ff_cropped;
ff_array = ff_array(~isnan(ff_array));
ff_array(ff_array<0) = 0;

r2star_array = nonzero_mask_shift.*r2star_cropped;
r2star_array = r2star_array(~isnan(r2star_array));
r2star_array(r2star_array<0) = 0;

figure();
subplot(1,2,1);
imagesc(nonzero_mask.*t1map_cropped); caxis([800 1800]); axis image;
colorbar;
colormap(brewermap([],'*RdYlBu'));
subplot(2,2,2);
scatter(ff_array, t1_array, 'LineWidth', 1.5);
xlabel('FF (%)'); ylabel('T1 (ms)');
subplot(2,2,4);
scatter(r2star_array, t1_array, 'LineWidth', 1.5);
xlabel('r2star (s-1)'); ylabel('T1 (ms)');

%% Load contours (Remote)
load(GetFullPath(cat(2, pwd, '/../../T1_Fat_Project/ContourData_old03012021/Ryn/Ryn_Exvivo/T1_CMR/MyoReference/myoRef.mat')))
ref_cropped = imcrop(squeeze(myoRefMask_3D(:,:,10)), [20 38 56 60]);

se = strel('disk',1);
t1map_temp = t1map_cropped .* ref_cropped;
t1map_temp_eroded = imerode(t1map_temp, se);
ind = find(t1map_temp_eroded);
sz = size(t1map_temp_eroded);
[row, col] = ind2sub(sz, ind)

row_shift = row + 1;
col_shift = col + 1;
mask_shift = zeros(sz);
ind_shift = sub2ind(sz, row_shift, col_shift);
mask_shift(ind_shift) = 1;

nonzero_mask = ones(size(t1map_temp_eroded));
nonzero_mask_shift = ones(size(t1map_temp_eroded));
nonzero_mask(t1map_temp_eroded == 0) = nan;
nonzero_mask_shift(mask_shift == 0) = nan;

t1_array = nonzero_mask.*t1map_cropped;
t1_array = t1_array(~isnan(t1_array));

ff_array = nonzero_mask_shift.*ff_cropped;
ff_array = ff_array(~isnan(ff_array));
ff_array(ff_array<0) = 0;

r2star_array = nonzero_mask_shift.*r2star_cropped;
r2star_array = r2star_array(~isnan(r2star_array));
r2star_array(r2star_array<0) = 0;

figure();
subplot(1,2,1);
imagesc(nonzero_mask.*t1map_cropped); caxis([800 1800]); axis image;
colorbar;
colormap(brewermap([],'*RdYlBu'));
subplot(2,2,2);
scatter(ff_array, t1_array, 'LineWidth', 1.5);
xlabel('FF (%)'); ylabel('T1 (ms)');
subplot(2,2,4);
scatter(r2star_array, t1_array, 'LineWidth', 1.5);
xlabel('r2star (s-1)'); ylabel('T1 (ms)');
