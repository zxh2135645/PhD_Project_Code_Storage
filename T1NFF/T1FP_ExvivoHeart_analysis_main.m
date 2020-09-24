% Analysis of non-hemorrhagic exvivo heart
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
%% Read IR-SE DICOM files
dicom_fields = {'InversionTime'};
whatsinit = cell(length(list_to_read), 1);
slice_data = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit{i}, slice_data{i}] = dicom23D(f, dicom_fields);
end

% Get IR_array
IR_array = zeros(1, length(slice_data));
for i = 1:length(slice_data)
    IR_array(i) = slice_data{i}.InversionTime;
end

%% I - DESCRIPTION
qMRinfo('inversion_recovery'); % Describe the model

% II - MODEL PARAMETERS
% a - create object
Model = inversion_recovery;

% b - modifiy options
Model = Custom_OptionsGUI(Model); % You need to close GUI to move on.

%% III - FIT EXPERIMENTAL DATASET
%Draw Masks
figure(); imagesc(whatsinit{1}); axis image; 
epi = drawpolygon(gca); % You may need to use roipoly and its equivalent if MATLAB version is earlier than 2020a
mask_epi = createMask(epi); % You may need to use roipoly and its equivalent if MATLAB version is earlier than 2020a

% Reshape matrix as [Width x Height x #Slice x #IR]
dicom_size = size(whatsinit{1});
dicom_reshape = zeros(dicom_size(1), dicom_size(2), 1, length(whatsinit));
for i = 1:length(whatsinit)
    dicom_reshape(:,:,1,i) = whatsinit{i};
end

% a - load experimental data
data = struct();

% IRData.mat contains [128  128    1    9] data.
%load('inversion_recovery_data/IRData.mat');
% Mask.mat contains [128  128] data.
%load('inversion_recovery_data/Mask.mat');
data.IRData= double(dicom_reshape);
data.Mask= double(mask_epi);

% Reset IRs
Model.Prot.IRData.Mat = IR_array.';
Model.voxelwise = 1;

% b- fit dataset
FitResults = FitData(data,Model,0);

% c- show fitting results
%qMRshowOutput(FitResults,data,Model);

%% Draw ROIs in N vials
dim = input('Dimension of vials (1 or 2): ');

mask_save = cat(2, subject_data_dir, 'roi_cmr.mat');
if ~exist(mask_save, 'file')
    if dim == 1
        N = input('Number of vials: ');
        roi = cell(N, 1);
        for i = 1:size(roi, 3)
            figure(); imagesc(FitResults.T1 .* mask_epi); axis image; caxis([0 2000]); colorbar;
            temp = drawpolygon(gca); % You may need to use roipoly and its equivalent if MATLAB version is earlier than 2020a
            roi{i} = createMask(temp);
        end
    elseif dim == 2
        row = input('Number of rows: ');
        col = input('Number of cols: ');
        N = row * col;
        roi = cell(row, col);
        for k = 1:row
            for j = 1:col
                figure(); imagesc(FitResults.T1 .* mask_epi); axis image; caxis([0 2000]); colorbar;
                temp = drawpolygon(gca); % You may need to use roipoly and its equivalent if MATLAB version is earlier than 2020a
                roi{k, j} = createMask(temp);
            end
        end
    end
else
    load(mask_save);
    roi = roi_cmr.vial_mask_cell;
    row = size(roi, 1);
    col = size(roi, 2);
    N = row * col;
end

composite = zeros(dicom_size(1), dicom_size(2));

if dim == 1
    mean_t1 = zeros(N, 1);
    for i = 1:N
        
        composite = composite + FitResults.T1 .* roi{i};
        t1_masked = nonzeros(FitResults.T1 .* roi{i});
        t1_real = t1_masked(~isnan(t1_masked));
        mean_t1(i) = mean(t1_real);
    end
elseif dim == 2
    mean_t1 = zeros(row, col);
    for k = 1:row
        for j = 1:col
            composite = composite + FitResults.T1 .* roi{k,j};
            t1_masked = nonzeros(FitResults.T1 .* roi{k,j});
            t1_real = t1_masked(~isnan(t1_masked));
            mean_t1(k,j) = mean(t1_real);
        end
    end
end

%% Display Phantom
figure(); imagesc(FitResults.T1 .* mask_epi); axis image; caxis([0 2000]); colorbar; axis off;

%% For Comparing T1 MOLLI to IRSE 

[list_to_read, order_to_read] = NamePicker(folder_glob);

dicom_fields = {'InversionTime'};
T1_molli = cell(length(list_to_read), 1);
slice_data = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [T1_molli{i}, slice_data{i}] = dicom23D(f, dicom_fields);
end

%% Draw ROIs in N vials for T1 MOLLI

if ~exist(mask_save, 'file')
    if dim == 1
        roi_molli = cell(N,1);
        for i = 1:N
            figure(); imagesc(T1_molli{1}); axis image; caxis([0 3000]); colorbar;
            temp = drawpolygon(gca); % You may need to use roipoly and its equivalent if MATLAB version is earlier than 2020a
            roi_molli{i} = createMask(temp);
        end
    elseif dim == 2
        roi_molli = cell(row, col);
        for k = 1:row
            for j = 1:col
                figure(); imagesc(T1_molli{1}); axis image; caxis([0 3000]); colorbar;
                temp = drawpolygon(gca); % You may need to use roipoly and its equivalent if MATLAB version is earlier than 2020a
                roi_molli{k, j} = createMask(temp);
            end
        end
    end
else
    load(mask_save);
    roi_molli = roi_cmr.vial_mask_cell;
    row = size(roi_molli, 1);
    col = size(roi_molli, 2);
    N = row * col;
end


if dim == 1
    % mean value of T1 MOLLI
    mean_t1_molli = zeros(N, 1);
    for i = 1:N
        t1_masked = nonzeros(T1_molli{1} .* roi_molli{i});
        t1_real = t1_masked(~isnan(t1_masked));
        mean_t1_molli(i) = mean(t1_real);
    end
elseif dim == 2
    % mean value of T1 MOLLI
    mean_t1_molli = zeros(row, col);
    for k = 1:row
        for j = 1:col
            t1_masked = nonzeros(T1_molli{1} .* roi_molli{k,j});
            t1_real = t1_masked(~isnan(t1_masked));
            mean_t1_molli(k,j) = mean(t1_real);
        end
    end
end

%% Display as comparison (create a circle to include everything)
figure(); imagesc(T1_molli{1}); axis image; 
epi = drawpolygon(gca); % You may need to use roipoly and its equivalent if MATLAB version is earlier than 2020a
mask_epi_molli = createMask(epi); % You may need to use roipoly and its equivalent if MATLAB version is earlier than 2020a

%% IRSE vs MOLLI
figure('Position', [0 100 1600 800]); 
subplot(1,2,1);
imagesc(FitResults.T1 .* mask_epi); title('T1-IRSE', 'FontSize', 24);
axis image; caxis([0 3000]); colorbar; axis off;

subplot(1,2,2);
imagesc(T1_molli{1} .* mask_epi_molli); title('T1-MOLLI', 'FontSize', 24)
axis image; caxis([0 3000]); colorbar; axis off;
%% May not be useful in terms of xticklabels
figure();
plot(mean_t1, 'LineWidth', 2); hold on;
plot(mean_t1_molli, 'LineWidth', 2); grid on;
ylabel('T1 value (ms)'); xlabel('Fat Fraction (%)');
xticks(1:8);
xticklabels({'0','2.5','5','10','20','30','40','50'})
set(gca, 'FontSize', 16);
legend({'T1 IRSE', 'T1 MOLLI'})

%% Save the saved results as data
ir_fitting = struct;
ir_fitting.FitResults = FitResults;
ir_fitting.composite = composite;
ir_fitting.mean_t1 = mean_t1;
save(cat(2, subject_data_dir, 'ir_fitting.mat'), 'ir_fitting');

%% Draw contours in heart
img = T1_molli{1};
myo_coords_cell = cell(size(img, 3), 2);
roi_save = cat(2, subject_data_dir, 'roi_cmr.mat');

if ~exist(roi_save, 'file')
for i = 1:(size(img, 3))
    disp('Epicardium: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    %caxis([0 100])
    epi = drawpolygon(gca);
    epi_coords = epi.Position;
    
    disp('Endocardium: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    %caxis([0 100])
    endo = drawpolygon(gca);
    endo_coords = endo.Position;
    
    disp('MI roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    %caxis([0 100])
    mi = drawpolygon(gca);
    mi_coords = mi.Position;
    
    disp('remote roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    %caxis([0 100])
    remote = drawpolygon(gca);
    remote_coords = remote.Position;
    
    disp('fluid roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    %caxis([0 100])
    fluid = drawpolygon(gca);
    fluid_coords = fluid.Position;
    
    disp('Center line roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    %caxis([0 100])
    center_line = drawpolygon(gca);
    center_coords = center_line.Position;
    
    disp('air roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
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

roi_cmr = struct;
roi_cmr.myo_coords_cell = myo_coords_cell;
roi_cmr.mi_coords = mi_coords;
roi_cmr.remote_coords = remote_coords;
roi_cmr.fluid_coords = fluid_coords;
roi_cmr.center_coords = center_coords;
roi_cmr.air_coords = air_coords;

save(roi_save, 'roi_cmr');

else
    load(roi_save);
    myo_coords_cell = roi.myo_coords_cell;
    mi_coords = roi.mi_coords;
    remote_coords = roi.remote_coords;
    fluid_coords = roi.fluid_coords;
    center_coords = roi.center_coords;
    air_coords = roi.air_coords;
end

img_size_truth = size(img);
mask_save = cat(2, subject_data_dir, 'mask_cmr.mat');

if ~exist(mask_save, 'file')
    figure();
    mask_cmr_struct = struct;
    for i = 1:size(img, 3)
        img2 = whatsinit{i};
        img2_size = size(whatsinit{i});
        ratio = img_size_truth(1) ./ img2_size(1);
        imagesc(img2); %caxis([0 100]);
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
        
        mask_cmr_struct(i).myo_mask = myo_mask;
        mask_cmr_struct(i).mi_mask = mi_mask;
        mask_cmr_struct(i).remote_mask = remote_mask;
        mask_cmr_struct(i).fluid_mask = fluid_mask;
        mask_cmr_struct(i).air_mask = air_mask;
        
        mask_cmr_struct(i).epi_mask = epi_mask;
        mask_cmr_struct(i).endo_mask = endo_mask;
        
        myo_mask_endo = myo_mask .* center_mask;
        myo_mask_epi = myo_mask - myo_mask_endo;
        mask_cmr_struct(i).myo_mask_endo = myo_mask_endo;
        mask_cmr_struct(i).myo_mask_epi = myo_mask_epi;
    end
    
    save(mask_save, 'mask_cmr_struct');
else
    load(mask_save);
end

%% Read CMR quantitative mapping
[list_to_read, order_to_read] = NamePicker(folder_glob);

cmr = cell(length(list_to_read), 1);
slice_data_cmr = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [cmr{i}, slice_data_cmr{i}] = dicom23D(f);
end
%% Get MI region values
cmr_parameter = zeros(length(cmr), 2);
for i = 1:length(cmr)
    cmr_parameter(i, 1) = mean(nonzeros(mask_cmr_struct.mi_mask .* cmr{i}));
    cmr_parameter(i, 2) = mean(nonzeros(mask_cmr_struct.remote_mask .* cmr{i})); 
end

%% Read CMR quantitative mapping
[list_to_read, order_to_read] = NamePicker(folder_glob);
siemens = cell(length(list_to_read), 1);
slice_data_siemens = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [siemens{i}, slice_data_siemens{i}] = dicom23D(f);
end

%% Get MI region values
Slc_Loc = zeros(length(slice_data_siemens{1}), 1);
for i = 1:length(slice_data_siemens{1})
    Slc_Loc(i) = slice_data_siemens{1}(i).SliceLocation;
end
Slc_Loc_cmr = slice_data_cmr{1}.SliceLocation;
[min_val,idx] = min(abs(Slc_Loc - Slc_Loc_cmr));


%% Draw contours in heart (Siemens)
img = siemens{1};
myo_coords_cell = cell(1, 2);
roi_save = cat(2, subject_data_dir, 'roi_siemens.mat');

if ~exist(roi_save, 'file')
for i = 1:(size(myo_coords_cell, 1))
    disp('Epicardium: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,idx)); axis image;
    %caxis([0 100])
    epi = drawpolygon(gca);
    epi_coords = epi.Position;
    
    disp('Endocardium: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,idx)); axis image;
    %caxis([0 100])
    endo = drawpolygon(gca);
    endo_coords = endo.Position;
    
    disp('MI roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,idx)); axis image;
    %caxis([0 100])
    mi = drawpolygon(gca);
    mi_coords = mi.Position;
    
    disp('remote roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,idx)); axis image;
    %caxis([0 100])
    remote = drawpolygon(gca);
    remote_coords = remote.Position;
    
    disp('fluid roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,idx)); axis image;
    %caxis([0 100])
    fluid = drawpolygon(gca);
    fluid_coords = fluid.Position;
    
    disp('Center line roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,idx)); axis image;
    %caxis([0 100])
    center_line = drawpolygon(gca);
    center_coords = center_line.Position;
    
    disp('air roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,idx)); axis image;
    %caxis([0 100])
    air = drawpolygon(gca);
    air_coords = air.Position;
    
    myo_coords_cell{1, 1} = epi.Position;
    myo_coords_cell{1, 2} = endo.Position;
    epi_mask = createMask(epi);
    endo_mask = createMask(endo);
    center_mask = createMask(center_line);
    
    close all;
end

roi_siemens = struct;
roi_siemens.myo_coords_cell = myo_coords_cell;
roi_siemens.mi_coords = mi_coords;
roi_siemens.remote_coords = remote_coords;
roi_siemens.fluid_coords = fluid_coords;
roi_siemens.center_coords = center_coords;
roi_siemens.air_coords = air_coords;
roi_siemens.idx = idx;

save(roi_save, 'roi_siemens');

else
    load(roi_save);
    myo_coords_cell = roi_siemens.myo_coords_cell;
    mi_coords = roi_siemens.mi_coords;
    remote_coords = roi_siemens.remote_coords;
    fluid_coords = roi_siemens.fluid_coords;
    center_coords = roi_siemens.center_coords;
    air_coords = roi_siemens.air_coords;
end

img_size_truth = size(img);
mask_save = cat(2, subject_data_dir, 'mask_siemens.mat');

if ~exist(mask_save, 'file')
    figure();
    mask_siemens_struct = struct;
    for i = 1:size(myo_coords_cell, 1)
        img2 = siemens{i};
        img2_size = size(siemens{i});
        ratio = img_size_truth(1) ./ img2_size(1);
        imagesc(img2(:,:,idx)); %caxis([0 100]);
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
        
        mask_siemens_struct(i).myo_mask = myo_mask;
        mask_siemens_struct(i).mi_mask = mi_mask;
        mask_siemens_struct(i).remote_mask = remote_mask;
        mask_siemens_struct(i).fluid_mask = fluid_mask;
        mask_siemens_struct(i).air_mask = air_mask;
        
        mask_siemens_struct(i).epi_mask = epi_mask;
        mask_siemens_struct(i).endo_mask = endo_mask;
        
        myo_mask_endo = myo_mask .* center_mask;
        myo_mask_epi = myo_mask - myo_mask_endo;
        mask_siemens_struct(i).myo_mask_endo = myo_mask_endo;
        mask_siemens_struct(i).myo_mask_epi = myo_mask_epi;
    end
    
    save(mask_save, 'mask_siemens_struct');
else
    load(mask_save);
end
%% Siemens parameter metrics
siemens_parameter = zeros(length(siemens), 2);
for i = 1:length(siemens)
    siemens_parameter(i, 1) = mean(nonzeros(mask_siemens_struct.mi_mask .* siemens{i}(:,:,idx)));
    siemens_parameter(i, 2) = mean(nonzeros(mask_siemens_struct.remote_mask .* siemens{i}(:,:,idx)));
end

siemens_parameter_save = cat(2, subject_data_dir, 'siemens_parameter.mat');
save(siemens_parameter_save, 'siemens_parameter');

%%
T1 = mask_cmr_struct.mi_mask .* FitResults.T1;
T1(isnan(T1)) = 0;
mean(nonzeros(mask_cmr_struct.mi_mask .* T1))

figure();
bar(siemens_parameter);
xticklabels({'T1', 'T2', 'T2star'}); grid on;