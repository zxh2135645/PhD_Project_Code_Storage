% Analysis of complete phantom version 1: 48 vials with iron conc. 
% 0-50 ug/mL and fat volume percentage from 0% to 50% (6x8)

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

roi_row = input('Please input number of rows in your phantom: ');
roi_col = input('Please input number of columns in your phantom: ');

%% Read CMR DICOM files (Single slice)
whatsinit = cell(length(list_to_read), 1);
slice_data = cell(length(list_to_read), 1);

for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit{i}, slice_data{i}] = dicom23D(f);
end

Slc_Loc_cmr = slice_data{1}.SliceLocation;

% Display images
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    subplot(2,2,i);
    imagesc(whatsinit{i}); axis image;
    %caxis([0 100]) 
end

%% Draw contours of each vial

img = whatsinit{3};
roi_save = cat(2, subject_data_dir, 'roi_cmr.mat');
vial_coords_cell = cell(roi_row, roi_col);
vial_mask_cell = cell(roi_row, roi_col);

if ~exist(roi_save, 'file')
    for i = 1:(size(img, 3))
        for j = 1:roi_row
            for k = 1:roi_col
                disp(['Vial #', num2str(k+(j-1)*roi_col)]);
                figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
                vial = drawpolygon(gca);
                vial_coords_cell{j,k} = vial.Position;
                
                
                vial_mask_cell{j,k} = createMask(vial);
                close all;
            end
        end
    end
    
    roi_cmr.vial_coords_cell = vial_coords_cell;
    roi_cmr.vial_mask_cell = vial_mask_cell;
    save(roi_save, 'roi_cmr');
    
else
    load(roi_save);
end

%% Read SIEMENS DICOM files (Multi-slice)
disp('Read SIEMENS multi-slice quantitative mapping then: ');
[list_to_read, order_to_read] = NamePicker(folder_glob);

whatsinit2 = cell(length(list_to_read), 1);
slice_data = cell(length(list_to_read), 1);

for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit2{i}, slice_data{i}] = dicom23D(f);
end

%% Find corresponding slice location
Slc_Loc = zeros(length(slice_data{1}), 1);
for i = 1:length(slice_data{1})
    Slc_Loc(i) = slice_data{1}(i).SliceLocation;
end

[min_val,idx] = min(abs(Slc_Loc - Slc_Loc_cmr));

% Draw contours of each vial (Different resolution in multi-slice)
img = whatsinit2{3}(:,:,idx);
roi_save = cat(2, subject_data_dir, 'roi.mat');
vial_coords_cell = cell(roi_row, roi_col);
vial_mask_cell = cell(roi_row, roi_col);

if ~exist(roi_save, 'file')
    for i = 1:(size(img, 3))
        for j = 1:roi_row
            for k = 1:roi_col
                disp(['Vial #', num2str(k+(j-1)*roi_col)]);
                figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;caxis([0 100])
                vial = drawpolygon(gca);
                vial_coords_cell{j,k} = vial.Position;
                
                
                vial_mask_cell{j,k} = createMask(vial);
                close all;
            end
        end
    end
    
    roi.vial_coords_cell = vial_coords_cell;
    roi.vial_mask_cell = vial_mask_cell;
    roi.idx = idx;
    save(roi_save, 'roi');
    
else
    load(roi_save);
end

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
imagesc(mask_composite); title('Composite image');

%% Heatmaps
figure();
subplot(2,2,1);
imagesc(t1_cmr);title('T1 MOLLI');
subplot(2,2,2);
imagesc(t2_cmr); title('T2 CMR');% Unable to read T2 CMR correctly
subplot(2,2,3);
imagesc(t2star_cmr); title('T2* CMR');
subplot(2,2,4);
r = t1_cmr / max(t1_cmr(:));
g = t2_cmr / max(t2_cmr(:));
b = t2star_cmr / max(t2star_cmr(:));
rgb_cmr = zeros(size(r,1), size(r,2), 3);
rgb_cmr(:,:,1) = r;
rgb_cmr(:,:,2) = g;
rgb_cmr(:,:,3) = b;
imagesc(rgb_cmr);

figure();
subplot(2,2,1);
imagesc(t1_siemens); title('T1 Mapping');
subplot(2,2,2);
imagesc(t2_siemens); title('T2 Mapping');
subplot(2,2,3);
imagesc(t2star_siemens); title('T2* Mapping');
subplot(2,2,4);
r = t1_siemens / max(t1_siemens(:));
g = t2_siemens / max(t2_siemens(:));
b = t2star_siemens / max(t2star_siemens(:));
rgb_siemens = zeros(size(r,1), size(r,2), 3);
rgb_siemens(:,:,1) = r;
rgb_siemens(:,:,2) = g;
rgb_siemens(:,:,3) = b;
imagesc(rgb_siemens); title('RGB');