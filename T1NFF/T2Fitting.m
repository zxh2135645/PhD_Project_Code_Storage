% T2 SE
% This is incorperated into T1FP analysis workflow
% Need masks from T1FP analysis 
%%%%%%%%%%%%%%%%%%%%%%%%%% input  file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%% output file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;

% Example for generating T1 Map MOLLI from qMRLab
% addpath('/Users/jameszhang/Documents/qMRLab');
% This line needs to be modified

addpath('../function/');
base_dir = uigetdir;
folder_glob = glob(cat(2, base_dir, '/*'));

labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT', 'IR'};

label = labels{3};
idx_array = contains(folder_glob, label);

proj_dir = GetFullPath(cat(2, pwd, '/../../T1_Fat_Project/'));
if ~exist(proj_dir, 'dir')
    mkdir(proj_dir)
end

data_dir = GetFullPath(cat(2, proj_dir, 'data/'));
if ~exist(data_dir, 'dir')
    mkdir(data_dir)
end

subject_name = input('Please type subject name here:  ', 's');
subject_data_dir = GetFullPath(cat(2, data_dir, subject_name, '/'));
if ~exist(subject_data_dir, 'dir')
    mkdir(subject_data_dir)
end

[list_to_read, order_to_read] = NamePicker(folder_glob);

%% Read IRSE data
dicom_fields = {'EchoTime'};
whatsinit = cell(length(list_to_read), 1);
slice_data = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit{i}, slice_data{i}] = dicom23D(f, dicom_fields);
end

% Get IR_array
TE_array = zeros(1, length(slice_data));
for i = 1:length(slice_data)
    TE_array(i) = slice_data{i}.EchoTime;
end

%% I - DESCRIPTION
qMRinfo('mono_t2'); % Describe the model

% II - MODEL PARAMETERS
% a - create object
Model = mono_t2;

% b - modifiy options
Model = Custom_OptionsGUI(Model);

%% III - FIT EXPERIMENTAL DATASET
%Draw Masks
figure(); imagesc(whatsinit{1}); axis image; 
epi = drawpolygon(gca); % You may need to use roipoly and its equivalent if MATLAB version is earlier than 2020a
mask_epi = createMask(epi); % You may need to use roipoly and its equivalent if MATLAB version is earlier than 2020a
%%
% Reshape matrix as [Width x Height x #Slice x #IR]
dicom_size = size(whatsinit{1});
dicom_reshape = zeros(dicom_size(1), dicom_size(2), 1, length(whatsinit));
for i = 1:length(whatsinit)
    dicom_size = size(whatsinit{i});
    dicom_reshape(:,:,:,i) = reshape(whatsinit{i}, dicom_size(1), dicom_size(2), 1, []);
end

% a - load experimental data
data = struct();

% IRData.mat contains [128  128    1    9] data.
%load('inversion_recovery_data/IRData.mat');
% Mask.mat contains [128  128] data.
%load('inversion_recovery_data/Mask.mat');
data.SEdata= double(dicom_reshape);
data.Mask= double(mask_epi);

% Reset IRs
Model.Prot.SEdata.Mat = TE_array.';
Model.st = [100 2000];
Model.lb = [1 1];
Model.ub = [300 1e4];
Model.fx = [0 1];
Model.voxelwise = 1;
%Model.options.FitType = 'Linear';
Model.options.FitType = 'Exponential';

% b- fit dataset
FitResults = FitData(data,Model,0);

%% Show image
figure();
imagesc(FitResults.T2); caxis([0 200]);
colormap(brewermap([],'*YlGnBu'));
title('T2 SE', 'FontSize', 24); axis off;
%% Read T2SE-WE data
[list_to_read, order_to_read] = NamePicker(folder_glob);

dicom_fields = {'EchoTime'};
whatsinit2 = cell(length(list_to_read), 1);
slice_data = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit2{i}, slice_data{i}] = dicom23D(f, dicom_fields);
end

%% III - FIT EXPERIMENTAL DATASET

% Reshape matrix as [Width x Height x #Slice x #IR]
dicom_size = size(whatsinit2{1});
dicom_reshape = zeros(dicom_size(1), dicom_size(2), 1, length(whatsinit2));
for i = 1:length(whatsinit2)
    dicom_size = size(whatsinit2{i});
    dicom_reshape(:,:,:,i) = reshape(whatsinit2{i}, dicom_size(1), dicom_size(2), 1, []);
end

% a - load experimental data
data = struct();

% IRData.mat contains [128  128    1    9] data.
%load('inversion_recovery_data/IRData.mat');
% Mask.mat contains [128  128] data.
%load('inversion_recovery_data/Mask.mat');
data.SEdata= double(dicom_reshape);
data.Mask= double(mask_epi);

% Reset IRs
Model.Prot.SEdata.Mat = TE_array.';
Model.st = [100 2000];
Model.lb = [1 1];
Model.ub = [300 1e4];
Model.fx = [0 1];
Model.voxelwise = 1;
%Model.options.FitType = 'Linear';
Model.options.FitType = 'Exponential';

% b- fit dataset
FitResults_WE = FitData(data,Model,0);

%% Show image
figure('Position', [100 0 1600 1600]);
subplot(1,2,1);
imagesc(FitResults.T2); caxis([0 100]);axis image;
title('T2SE', 'FontSize', 20); axis off;
subplot(1,2,2);
imagesc(FitResults_WE.T2); caxis([0 100]); axis image;
title('T2SE-WE', 'FontSize', 20); axis off;
colormap(brewermap([],'*YlGnBu'));

%% Plot mean values
dim = input('Dimension of vials (1 or 2): ');
mask_save = cat(2, subject_data_dir, 'roi_cmr.mat');
img = FitResults.T2;
caxis_rg = [0 200];

[roi, row, col, N] = Func_DrawROI_inPhantom(img, mask_epi, mask_save, caxis_rg, dim);

composite = zeros(dicom_size(1), dicom_size(2));
if dim == 1
    mean_t2 = zeros(N, 1);
    for i = 1:N
        
        composite = composite + FitResults.T2 .* roi.vial_mask_cell{i};
        t2_masked = nonzeros(FitResults.T2 .* roi.vial_mask_cell{i});
        t2_real = t2_masked(~isnan(t2_masked));
        mean_t2(i) = mean(t2_real);
    end
elseif dim == 2
    mean_t2 = zeros(row, col);
    for k = 1:row
        for j = 1:col
            composite = composite + FitResults.T2 .* roi.vial_mask_cell{k,j};
            t2_masked = nonzeros(FitResults.T2 .* roi.vial_mask_cell{k,j});
            t2_real = t2_masked(~isnan(t2_masked));
            mean_t2(k,j) = mean(t2_real);
        end
    end
end

%% mean value of se-we
composite_we = zeros(dicom_size(1), dicom_size(2));
if dim == 1
    mean_t2_we = zeros(N, 1);
    for i = 1:N
        composite_we = composite_we + FitResults_WE.T2 .* roi.vial_mask_cell{i};
        t2_masked = nonzeros(FitResults_WE.T2 .* roi.vial_mask_cell{i});
        t2_real = t1_masked(~isnan(t2_masked));
        mean_t2_we(i) = mean(t2_real);
    end
elseif dim == 2
    mean_t2_we = zeros(row, col);
    for k = 1:row
        for j = 1:col
            composite_we = composite_we + FitResults_WE.T2 .* roi.vial_mask_cell{k,j};
            t2_masked = nonzeros(FitResults_WE.T2 .* roi.vial_mask_cell{k,j});
            t2_real = t2_masked(~isnan(t2_masked));
            mean_t2_we(k,j) = mean(t2_real);
        end
    end
end
%% FF (optional)
ff = {{'2.5','5','10','20','30','40','100'}, {'2.5','5','10','20','30','40','100'}};
row_info = {'W&F', 'W Only'};
xlabel_info = {'Fat Fraction (%)', 'Fat Fraction (%)'};
figure('Position', [100 0 800 800]);
for i = 1:size(mean_t2, 1)
    subplot(1,2,i);
    plot(mean_t2(i,:), 'LineWidth', 2); hold on;
    plot(mean_t2_we(i,:), 'LineWidth', 2); grid on;
    ylabel('T1 value (ms)'); xlabel(xlabel_info{i});
    xticks(1:size(mean_t2, 2));
    xticklabels(ff{i})
    set(gca, 'FontSize', 16);
    legend({'T2 SE', 'T2 SE-WE'}, 'Location', 'SouthWest')
    title(cat(2, row_info{i}));
    ylim([0 100]);
end
%% Save the saved results as data
se_fitting = struct;
se_fitting.FitResults = FitResults;
se_fitting.composite = composite;
se_fitting.mean_t2 = mean_t2;
save(cat(2, subject_data_dir, 'se_fitting.mat'), 'se_fitting');

se_fitting_we = struct;
se_fitting_we.FitResults = FitResults_WE;
se_fitting_we.composite = composite_we;
se_fitting_we.mean_t1 = mean_t2_we;
save(cat(2, subject_data_dir, 'se_fitting_we.mat'), 'se_fitting_we');