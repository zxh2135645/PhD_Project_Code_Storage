% Analysis of complete phantom version 1: 48 vials with iron conc. 
% 0-50 ug/mL and fat volume percentage from 0% to 50% (6x8)

%% This Is Not Main for in-vivo analysis
clear all;
close all;

addpath('../function/');

% The reading of T2 mapping and T2* mapping CMR does not work
%%%%%%%%%%%%%%%%%%%%%%%%%% output file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
% 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_dir = uigetdir;
folder_glob = glob(cat(2, base_dir, '\*'));

labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};

label = labels{5};
idx_array = contains(folder_glob, label);

disp('Read CMR single slice quantitative mapping first: ');
[list_to_read, order_to_read] = NamePicker(folder_glob);

proj_dir = GetFullPath(cat(2, pwd, '/../../T1_Fat_Project/'));
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

%% Read CMR DICOM files (Single slice)
whatsinit = cell(length(list_to_read), 1);
slice_data = cell(length(list_to_read), 1);
dicom_fields = {'RescaleSlope', ...
                'SliceLocation'};


for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit{i}, slice_data{i}] = dicom23D(f, dicom_fields);
    if ~isempty(slice_data{i}.RescaleSlope)
        whatsinit{i} = whatsinit{i} .* slice_data{i}.RescaleSlope;
    end
end

Slc_Loc_cmr = slice_data{1}.SliceLocation;

%% Display images
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    subplot(2,2,i);
    imagesc(whatsinit{i}); axis image;
    
    if i == 1
        caxis([0 3000]);
    else
        caxis([0 100]);
    end
    colormap(brewermap([],'*RdBu'));
    colorbar;
end

%% Draw contours of each vial
figure(); imagesc(whatsinit{3}); axis image; 
epi = drawpolygon(gca); % You may need to use roipoly and its equivalent if MATLAB version is earlier than 2020a
mask_epi = createMask(epi);

dim = input('Dimension of vials (1 or 2): ');
roi_save = cat(2, subject_data_dir, 'roi_cmr.mat');
caxis_rg = [0 200];
img = whatsinit{3};

[roi_cmr, roi_row, roi_col, N] = Func_DrawROI_inPhantom(img, mask_epi, roi_save, caxis_rg, dim);
save(roi_save, 'roi_cmr');

%% Read SIEMENS DICOM files (Multi-slice)
disp('Read SIEMENS multi-slice quantitative mapping then: ');
[list_to_read, order_to_read] = NamePicker(folder_glob);
dicom_fields = {'RescaleSlope', ...
                'SliceLocation'};
            
whatsinit2 = cell(length(list_to_read), 1);
slice_data = cell(length(list_to_read), 1);

for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit2{i}, slice_data{i}] = dicom23D(f, dicom_fields);
    if ~isempty(slice_data{i}(1).RescaleSlope)
        whatsinit{i} = whatsinit{i} .* slice_data{i}(1).RescaleSlope;
    end
end

%% Find corresponding slice location
idx_array = zeros(length(whatsinit2), 1);
for j = 1:length(whatsinit2)
    Slc_Loc = zeros(length(slice_data{j}), 1);
    for i = 1:length(slice_data{j})
        Slc_Loc(i) = slice_data{j}(i).SliceLocation;
    end
    
    [min_val,idx] = min(abs(Slc_Loc - Slc_Loc_cmr));
    idx_array(j) = idx;
end

% Draw contours of each vial (Different resolution in multi-slice)
img = whatsinit2{3}(:,:,idx_array(3));
roi_save = cat(2, subject_data_dir, 'roi.mat');
caxis_rg = [0 100];

figure(); imagesc(img); axis image; caxis(caxis_rg);
epi = drawpolygon(gca); % You may need to use roipoly and its equivalent if MATLAB version is earlier than 2020a
mask_epi = createMask(epi);



[roi, roi_row, roi_col, N] = Func_DrawROI_inPhantom(img, mask_epi, roi_save, caxis_rg, dim);
save(roi_save, 'roi');

%% Get Quantitative mapping
vial_mask_cell_cmr = roi_cmr.vial_mask_cell;
vial_mask_cell_siemens = roi.vial_mask_cell;

mask_composite = zeros(size(vial_mask_cell_cmr{1}));
t1_cmr = zeros(roi_row, roi_col);
t2_cmr = zeros(roi_row, roi_col);
t2star_cmr = zeros(roi_row, roi_col);

t1_siemens = zeros(roi_row, roi_col);
t2_siemens = zeros(roi_row, roi_col);
t2star_siemens = zeros(roi_row, roi_col);

for j = 1:roi_row
    for k = 1:roi_col
        ind = k+(j-1)*roi_col;
        mask_composite = mask_composite + ind*vial_mask_cell_cmr{j,k};
        
        t1_cmr(j,k) = mean(nonzeros(whatsinit{1} .* vial_mask_cell_cmr{j,k}));
        t2_cmr(j,k) = mean(nonzeros(whatsinit{2} .* vial_mask_cell_cmr{j,k}));
        t2star_cmr(j,k) = mean(nonzeros(whatsinit{3} .* vial_mask_cell_cmr{j,k}));
        
        t1_siemens(j,k) = mean(nonzeros(whatsinit2{1} .* vial_mask_cell_siemens{j,k}));
        t2_siemens(j,k) = mean(nonzeros(whatsinit2{2} .* vial_mask_cell_siemens{j,k}));
        t2star_siemens(j,k) = mean(nonzeros(whatsinit2{3} .* vial_mask_cell_siemens{j,k}));
    end
end

figure();
imagesc(mask_composite); title('Composite image'); colormap(brewermap([],'YlGnBu'));


%% Heatmaps
figure();
subplot(2,2,1);
imagesc(t1_cmr);title('T1 MOLLI');colormap(brewermap([],'*YlGnBu'));
colorbar;
subplot(2,2,2);
imagesc(t2_cmr); title('T2 CMR');% Unable to read T2 CMR correctly
colorbar;
subplot(2,2,3);
imagesc(t2star_cmr); title('T2* CMR'); 
colorbar;
subplot(2,2,4);
r = t1_cmr / max(t1_cmr(:));
g = t2_cmr / max(t2_cmr(:));
b = t2star_cmr / max(t2star_cmr(:));
rgb_cmr = zeros(size(r,1), size(r,2), 3);
rgb_cmr(:,:,1) = r;
rgb_cmr(:,:,2) = g;
rgb_cmr(:,:,3) = b;
imagesc(rgb_cmr); title('RGB'); 

figure();
subplot(2,2,1);
imagesc(t1_siemens); title('T1 Mapping');colormap(brewermap([],'*YlGnBu'));
colorbar;
subplot(2,2,2);
imagesc(t2_siemens); title('T2 Mapping');
colorbar;
subplot(2,2,3);
imagesc(t2star_siemens); title('T2* Mapping');
colorbar;
subplot(2,2,4);
r = t1_siemens / max(t1_siemens(:));
g = t2_siemens / max(t2_siemens(:));
b = t2star_siemens / max(t2star_siemens(:));
rgb_siemens = zeros(size(r,1), size(r,2), 3);
rgb_siemens(:,:,1) = r;
rgb_siemens(:,:,2) = g;
rgb_siemens(:,:,3) = b;
imagesc(rgb_siemens); title('RGB');

%% Save mean value of CMR and Siemens maps
map_save = cat(2, subject_data_dir, 'maps.mat');
maps = struct;
maps.t1_cmr = t1_cmr;
maps.t2_cmr = t2_cmr;
maps.t2star_cmr = t2star_cmr;

maps.t1_siemens = t1_siemens;
maps.t2_siemens = t2_siemens;
maps.t2star_siemens = t2star_siemens;

save(map_save, 'maps');