% Analysis of complete phantoms with various iron conc. 
% 0-50 ug/mL and fat volume percentage from 0% to 50% (6x8)
% Specifically for weighted image


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

%% Display last weighted images
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    subplot(2,2,i);
    imagesc(whatsinit{i}(:,:,end)); axis image;
    caxis([0 100]) 
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

%% Get Weighted images
vial_mask_cell_cmr = roi_cmr.vial_mask_cell;
vial_mask_cell_siemens = roi.vial_mask_cell;

mask_composite = zeros(size(vial_mask_cell_cmr{1}));
t1_cmr_mean = zeros(roi_row, roi_col, size(whatsinit{1}, 3));
t2_cmr_mean = zeros(roi_row, roi_col, size(whatsinit{2}, 3));
t2star_cmr_mean = zeros(roi_row, roi_col, size(whatsinit{3}, 3));
t1_cmr_sd = zeros(roi_row, roi_col, size(whatsinit{1}, 3));
t2_cmr_sd = zeros(roi_row, roi_col, size(whatsinit{2}, 3));
t2star_cmr_sd = zeros(roi_row, roi_col, size(whatsinit{3}, 3));

t1_siemens1_mean = zeros(roi_row, roi_col, size(whatsinit2{1}, 3)/size(whatsinit2{1}, 3));
t1_siemens2_mean = zeros(roi_row, roi_col, size(whatsinit2{2}, 3)/size(whatsinit2{1}, 3));
t2_siemens_mean = zeros(roi_row, roi_col, size(whatsinit2{3}, 3)/size(whatsinit2{1}, 3));
t2star_siemens_mean = zeros(roi_row, roi_col, size(whatsinit2{4}, 3)/size(whatsinit2{1}, 3));
t1_siemens1_sd = zeros(roi_row, roi_col, size(whatsinit2{1}, 3)/size(whatsinit2{1}, 3));
t1_siemens2_sd = zeros(roi_row, roi_col, size(whatsinit2{2}, 3)/size(whatsinit2{1}, 3));
t2_siemens_sd = zeros(roi_row, roi_col, size(whatsinit2{3}, 3)/size(whatsinit2{1}, 3));
t2star_siemens_sd = zeros(roi_row, roi_col, size(whatsinit2{4}, 3)/size(whatsinit2{1}, 3));

t1_siemens1 = zeros(size(whatsinit2{1},1), size(whatsinit2{1},2), 1);
t1_siemens2 = zeros(size(whatsinit2{1},1), size(whatsinit2{1},2), 1);
t2_siemens = zeros(size(whatsinit2{3},1), size(whatsinit2{3},2), size(whatsinit2{3}, 3)/size(whatsinit2{1}, 3));
t2star_siemens = zeros(size(whatsinit2{4},1), size(whatsinit2{4},2), size(whatsinit2{4}, 3)/size(whatsinit2{1}, 3));

for j = 1:roi_row
    for k = 1:roi_col
        ind = k+(j-1)*roi_col;
        mask_composite = mask_composite + ind*vial_mask_cell_cmr{j,k};
        % T1 MOLLI
        t1_cmr = whatsinit{1};
        for m = 1:size(whatsinit{1}, 3)
            t1_cmr_mean(j,k,m) = mean(nonzeros(whatsinit{1}(:,:,m) .* vial_mask_cell_cmr{j,k}));
            t1_cmr_sd(j,k,m) = std(nonzeros(whatsinit{1}(:,:,m) .* vial_mask_cell_cmr{j,k}));
        end
        
        % T2 CMR
        t2_cmr = whatsinit{2};
        for m = 1:size(whatsinit{2}, 3)
            t2_cmr_mean(j,k,m) = mean(nonzeros(whatsinit{2}(:,:,m) .* vial_mask_cell_cmr{j,k}));
            t2_cmr_sd(j,k,m) = std(nonzeros(whatsinit{2}(:,:,m) .* vial_mask_cell_cmr{j,k}));
        end
        
        % T2* CMR
        t2star_cmr = whatsinit{3};
        for m = 1:size(whatsinit{2}, 3)
            t2star_cmr_mean(j,k,m) = mean(nonzeros(whatsinit{3}(:,:,m) .* vial_mask_cell_cmr{j,k}));
            t2star_cmr_sd(j,k,m) = std(nonzeros(whatsinit{3}(:,:,m) .* vial_mask_cell_cmr{j,k}));
            
        end
        
        % T1 Siemens
        n = size(whatsinit2{1}, 3)/size(whatsinit2{1}, 3);
        for m = 1:size(whatsinit2{1}, 3)/size(whatsinit2{1}, 3)
            ind = n*(idx-1)+m;
            t1_siemens1(:,:,m) = whatsinit2{1}(:,:,ind);
            t1_siemens2(:,:,m) = whatsinit2{2}(:,:,ind);
            t1_siemens1_mean(j,k,m) = mean(nonzeros(whatsinit2{1}(:,:,ind) .* vial_mask_cell_siemens{j,k}));
            t1_siemens2_mean(j,k,m) = mean(nonzeros(whatsinit2{2}(:,:,ind) .* vial_mask_cell_siemens{j,k}));
            t1_siemens1_sd(j,k,m) = std(nonzeros(whatsinit2{1}(:,:,ind) .* vial_mask_cell_siemens{j,k}));
            t1_siemens2_sd(j,k,m) = std(nonzeros(whatsinit2{2}(:,:,ind) .* vial_mask_cell_siemens{j,k}));
            %figure(); imagesc(whatsinit2{2}(:,:,ind));
        end
        
        % T2 Siemens
        n = size(whatsinit2{3}, 3)/size(whatsinit2{1}, 3);
        for m = 1:size(whatsinit2{3}, 3)/size(whatsinit2{1}, 3)
            ind = n*(idx-1)+m;
            t2_siemens(:,:,m) = whatsinit2{3}(:,:,ind);
            t2_siemens(j,k,m) = mean(nonzeros(whatsinit2{3}(:,:,ind) .* vial_mask_cell_siemens{j,k}));
            t2_siemens_sd(j,k,m) = std(nonzeros(whatsinit2{3}(:,:,ind) .* vial_mask_cell_siemens{j,k}));
            %figure(); imagesc(whatsinit2{3}(:,:,ind));
        end
        
        % T2 star Siemens
        n = size(whatsinit2{4}, 3)/size(whatsinit2{1}, 3);
        for m = 1:size(whatsinit2{4}, 3)/size(whatsinit2{1}, 3)
            ind = n*(idx-1)+m;
            t2star_siemens(:,:,m) = whatsinit2{4}(:,:,ind);
            t2star_siemens(j,k) = mean(nonzeros(whatsinit2{4}(:,:,ind) .* vial_mask_cell_siemens{j,k}));
            t2star_siemens_sd(j,k) = std(nonzeros(whatsinit2{4}(:,:,ind) .* vial_mask_cell_siemens{j,k}));
            %figure(); imagesc(whatsinit2{4}(:,:,ind));
        end
    end
end

figure();
imagesc(mask_composite); title('Composite image');

% Save as mat
weighted_intensity = struct;
weighted_intensity.t1_cmr = t1_cmr;
weighted_intensity.t1_cmr_mean = t1_cmr_mean;
weighted_intensity.t1_cmr_sd = t1_cmr_sd;

weighted_intensity.t2_cmr = t2_cmr;
weighted_intensity.t2_cmr_mean = t2_cmr_mean;
weighted_intensity.t2_cmr_sd = t2_cmr_sd;

weighted_intensity.t2star_cmr = t2star_cmr;
weighted_intensity.t2star_cmr_mean = t2star_cmr_mean;
weighted_intensity.t2star_cmr_sd = t2star_cmr_sd;

weighted_intensity.t1_siemens1 = t1_siemens1;
weighted_intensity.t1_siemens1_mean = t1_siemens1_mean;
weighted_intensity.t1_siemens1_sd = t1_siemens1_sd;
weighted_intensity.t1_siemens2 = t1_siemens2;
weighted_intensity.t1_siemens2_mean = t1_siemens2_mean;
weighted_intensity.t1_siemens1_sd = t1_siemens2_sd;

weighted_intensity.t2_siemens = t2_siemens;
weighted_intensity.t2_siemens_mean = t2_siemens_mean;
weighted_intensity.t2_siemens_sd = t2_siemens_sd;

weighted_intensity.t2star_siemens = t2star_siemens;
weighted_intensity.t2star_siemens_mean = t2star_siemens_mean;
weighted_intensity.t2star_siemens_sd = t2star_siemens_sd;

save(cat(2, subject_data_dir, 'weighted_intensity.mat'), 'weighted_intensity');

%% Try to analyse it
figure();
for i = 1:size(t1_cmr_mean, 3)
    subplot(3,3,i);
    imagesc(t1_cmr_mean(:,:,i)); colorbar;
end

%%
figure();
% Go over ff
for i = 1:size(t1_cmr_mean, 1)
    p = plot(squeeze(t1_cmr_mean(i,1,:)), 'LineWidth', 2); 
    hold on;
end
legend({'1', '2', '3', '4', '5'});

figure();
% Go over iron
for i = 1:size(t1_cmr_mean, 2)
    p = plot(squeeze(t1_cmr_mean(1,i,:)), 'LineWidth', 2); 
    hold on;
end
legend({'1', '2', '3'});