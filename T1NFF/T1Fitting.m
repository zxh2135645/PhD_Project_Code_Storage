% T1 MOLLI
% This is incorperated into T1FP analysis workflow
% Need masks from T1FP analysis 
%%%%%%%%%%%%%%%%%%%%%%%%%% input  file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%% output file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ir_weighted_metrics.mat
% ir_fitting.mat
% ir_fitting_we.mat (OPTIONAL)
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

label = labels{11};
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
img = FitResults.T1;
caxis_rg = [0 2000];

[roi, row, col, N] = Func_DrawROI_inPhantom(img, mask_epi, mask_save, caxis_rg, dim);
roi_cmr = roi;
save(mask_save, 'roi_cmr');

composite = zeros(dicom_size(1), dicom_size(2));
if dim == 1
    mean_t1 = zeros(N, 1);
    for i = 1:N
        
        composite = composite + FitResults.T1 .* roi.vial_mask_cell{i};
        t1_masked = nonzeros(FitResults.T1 .* roi.vial_mask_cell{i});
        t1_real = t1_masked(~isnan(t1_masked));
        mean_t1(i) = mean(t1_real);
    end
elseif dim == 2
    mean_t1 = zeros(row, col);
    for k = 1:row
        for j = 1:col
            composite = composite + FitResults.T1 .* roi.vial_mask_cell{k,j};
            t1_masked = nonzeros(FitResults.T1 .* roi.vial_mask_cell{k,j});
            t1_real = t1_masked(~isnan(t1_masked));
            mean_t1(k,j) = mean(t1_real);
        end
    end
end
%% Display Phantom
figure(); imagesc(FitResults.T1 .* mask_epi); axis image; caxis([0 3000]); colorbar; axis off;
colormap(brewermap([],'*YlGnBu'));
%% For Comparing T1 MOLLI to IRSE 
folder_glob = glob(cat(2, base_dir, '/*'));
[list_to_read, order_to_read] = NamePicker(folder_glob);

T1_molli = cell(length(list_to_read), 1);
slice_data = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [T1_molli{i}, slice_data{i}] = dicom23D(f);
end
%% Display as comparison (create a circle to include everything)
figure(); imagesc(T1_molli{1}); axis image; 
epi = drawpolygon(gca); % You may need to use roipoly and its equivalent if MATLAB version is earlier than 2020a
mask_epi_molli = createMask(epi); % You may need to use roipoly and its equivalent if MATLAB version is earlier than 2020a
%% Draw ROIs in N vials for T1 MOLLI
mask_save = cat(2, subject_data_dir, 'roi_cmr.mat');
img = T1_molli{1};
caxis_rg = [0 2000];
[roi, row, col, N] = Func_DrawROI_inPhantom(img, mask_epi_molli, mask_save, caxis_rg, dim);
fn = fieldnames(roi);
roi_molli = roi;

if dim == 1
    % mean value of T1 MOLLI
    mean_t1_molli = zeros(N, 1);
    for i = 1:N
        t1_masked = nonzeros(T1_molli{1} .* roi_molli.vial_mask_cell{i});
        t1_real = t1_masked(~isnan(t1_masked));
        mean_t1_molli(i) = mean(t1_real);
    end
elseif dim == 2
    % mean value of T1 MOLLI
    mean_t1_molli = zeros(row, col);
    for k = 1:row
        for j = 1:col
            t1_masked = nonzeros(T1_molli{1} .* roi_molli.vial_mask_cell{k,j});
            t1_real = t1_masked(~isnan(t1_masked));
            mean_t1_molli(k,j) = mean(t1_real);
        end
    end
end
    


%% IRSE vs MOLLI
figure('Position', [0 100 1600 800]); 
subplot(1,2,1);
imagesc(FitResults.T1 .* mask_epi); title('T1-IRSE', 'FontSize', 24);
axis image; caxis([0 3000]); colorbar; axis off;
colormap(brewermap([],'*YlGnBu'));


subplot(1,2,2);
imagesc(T1_molli{1} .* mask_epi_molli); title('T1-MOLLI', 'FontSize', 24)
axis image; caxis([0 3000]); colorbar; axis off;
%% May not be useful in terms of xticklabels
% Iron concentration
iron_conc = {'0', '10', '20', '30', '40', '50'};
figure();
for i = 1:size(mean_t1, 2)
    subplot(2,3,i)
    plot(mean_t1(:,i), 'LineWidth', 2); hold on;
    plot(mean_t1_molli(:,i), 'LineWidth', 2); grid on;
    ylabel('T1 value (ms)'); xlabel('Fat Fraction (%)');
    xticks(1:8);
    xticklabels({'0','2.5','5','10','20','30','40','100'})
    set(gca, 'FontSize', 16);
    legend({'T1 IRSE', 'T1 MOLLI'})
    title(cat(2, 'c(Fe) ', iron_conc{i}, ' ug/mL'))
end

%% FF (optional)
ff = {{'2.5','5','10','20','30', '40', '100'}, {'2.5','5','10','20','30', '40', '100'}};
row_info = {'W&F', 'W Only'};
xlabel_info = {'Fat Fraction (%)', 'Gd Concentration (mM)'};
figure('Position', [100 0 800 800]);
for i = 1:size(mean_t1, 1)
    subplot(1,2,i);
    plot(mean_t1(i,:), 'LineWidth', 2); hold on;
    plot(mean_t1_molli(i,:), 'LineWidth', 2); grid on;
    ylabel('T1 value (ms)'); xlabel(xlabel_info{i});
    xticks(1:size(mean_t1, 2));
    xticklabels(ff{i})
    set(gca, 'FontSize', 16);
    legend({'T1 IRSE', 'T1 MOLLI'}, 'Location', 'SouthWest')
    title(cat(2, row_info{i}));
    ylim([0 3000]);
end
%% Save the saved results as data
ir_fitting = struct;
ir_fitting.FitResults = FitResults;
ir_fitting.composite = composite;
ir_fitting.mean_t1 = mean_t1;
save(cat(2, subject_data_dir, 'ir_fitting.mat'), 'ir_fitting');
%% Read IRSE data (WE) (OPTIONAL)

folder_glob = glob(cat(2, base_dir, '/*'));
[list_to_read, order_to_read] = NamePicker(folder_glob);

dicom_fields = {'InversionTime'};
whatsinit_we = cell(length(list_to_read), 1);
slice_data = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit_we{i}, slice_data{i}] = dicom23D(f, dicom_fields);
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
figure(); imagesc(whatsinit_we{1}); axis image; 
epi = drawpolygon(gca); % You may need to use roipoly and its equivalent if MATLAB version is earlier than 2020a
mask_epi = createMask(epi); % You may need to use roipoly and its equivalent if MATLAB version is earlier than 2020a

% Reshape matrix as [Width x Height x #Slice x #IR]
dicom_size = size(whatsinit_we{1});
dicom_reshape = zeros(dicom_size(1), dicom_size(2), 1, length(whatsinit));
for i = 1:length(whatsinit_we)
    dicom_reshape(:,:,1,i) = whatsinit_we{i};
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

%% Draw ROIs in N vials (WE)
mask_save = cat(2, subject_data_dir, 'roi_cmr.mat');
img = FitResults.T1;
caxis_rg = [0 2000];
[roi, row, col, N] = Func_DrawROI_inPhantom(img, mask_epi, mask_save, caxis_rg, dim);
fn = fieldnames(roi);
roi_cmr_we = roi;

composite = zeros(dicom_size(1), dicom_size(2));

if dim == 1
    mean_t1_we = zeros(N, 1);
    for i = 1:N
        
        composite = composite + FitResults.T1 .* roi_cmr_we.vial_mask_cell{i};
        t1_masked = nonzeros(FitResults.T1 .* roi_cmr_we.vial_mask_cell{i});
        t1_real = t1_masked(~isnan(t1_masked));
        mean_t1_we(i) = mean(t1_real);
    end
elseif dim == 2
    mean_t1_we = zeros(row, col);
    for k = 1:row
        for j = 1:col
            composite = composite + FitResults.T1 .* roi_cmr_we.vial_mask_cell{k,j};
            t1_masked = nonzeros(FitResults.T1 .* roi_cmr_we.vial_mask_cell{k,j});
            t1_real = t1_masked(~isnan(t1_masked));
            mean_t1_we(k,j) = mean(t1_real);
        end
    end
end

%% Save the saved results as data
ir_fitting_we = struct;
ir_fitting_we.FitResults = FitResults;
ir_fitting_we.composite = composite;
ir_fitting_we.mean_t1 = mean_t1_we;
save(cat(2, subject_data_dir, 'ir_fitting_we.mat'), 'ir_fitting_we');

%% display IRSE and IRSE-WE
figure();
for i = 1:length(whatsinit)
    subplot(3,3,i)
    imagesc(whatsinit{i});
end

figure();
for i = 1:length(whatsinit_we)
    subplot(3,3,i)
    imagesc(whatsinit_we{i});
end

%% display IRSE and IRSE-WE

figure();
subplot(1,2,1);
imagesc(ir_fitting.FitResults.T1); axis image; caxis([0 3000])
title('IRSE', 'FontSize', 24); colormap(brewermap([],'*YlGnBu'));

subplot(1,2,2);
imagesc(ir_fitting_we.FitResults.T1); axis image; caxis([0 3000])
title('IRSE-WE', 'FontSize', 24); colormap(brewermap([],'*YlGnBu'));

%% display IRSE vs T1 MOLLI weighted
%% For Comparing T1 MOLLI weighted image
folder_glob = glob(cat(2, base_dir, '/*'));
[list_to_read, order_to_read] = NamePicker(folder_glob);

dicom_fields = {'InversionTime'};
T1_molli_w = cell(length(list_to_read), 1);
slice_data = cell(length(list_to_read), 1);


for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [T1_molli_w{i}, slice_data{i}] = dicom23D(f, dicom_fields);
end

% Get IR_array
IR_array_molli = zeros(1, length(slice_data{1}));
temp = slice_data{1};
for i = 1:length(slice_data{1})
    IR_array_molli(i) = temp(i).InversionTime;
end

%% Calculate T1 MOLLI
if dim == 1
    % mean value of T1 MOLLI
    mean_t1_molli_w = zeros(N, 8);
    sd_t1_molli_w = zeros(N, 8);
    for i = 1:N
        for t = 1:size(mean_t1_molli_w, 2)
            t1_masked = nonzeros(T1_molli_w{1}(:,:,t) .* roi_molli.vial_mask_cell{i});
            t1_real = t1_masked(~isnan(t1_masked));
            mean_t1_molli_w(i, t) = mean(t1_real);
            sd_t1_molli_w(i, t) = std(t1_real);
        end
    end
elseif dim == 2
    % mean value of T1 MOLLI
    mean_t1_molli_w = zeros(row, col, 8);
    sd_t1_molli_w = zeros(row, col, 8);
    for k = 1:row
        for j = 1:col
            for t = 1:size(mean_t1_molli_w, 3)
                t1_masked = nonzeros(T1_molli_w{1}(:,:,t) .* roi_molli.vial_mask_cell{k,j});
                t1_real = t1_masked(~isnan(t1_masked));
                mean_t1_molli_w(k,j,t) = mean(t1_real);
                sd_t1_molli_w(k,j,t) = std(t1_real);
            end
        end
    end
end
%% IRSE weighted image

if dim == 1
    % mean value of IRSE
    mean_t1_irse_w = zeros(N, length(whatsinit));
    sd_t1_irse_w = zeros(N, length(whatsinit));
    for i = 1:N
        for t = 1:size(mean_t1_irse_w, 2)
            t1_masked = nonzeros(whatsinit{t} .* roi_cmr.vial_mask_cell{i});
            t1_real = t1_masked(~isnan(t1_masked));
            mean_t1_irse_w(i, t) = mean(t1_real);
            sd_t1_irse_w(i, t) = std(t1_real);
        end
    end
elseif dim == 2
    % mean value of T1 MOLLI
    mean_t1_irse_w = zeros(row, col, length(whatsinit));
    sd_t1_irse_w = zeros(row, col, length(whatsinit));
    for k = 1:row
        for j = 1:col
            for t = 1:size(mean_t1_irse_w, 3)
                t1_masked = nonzeros(whatsinit{t} .* roi_cmr.vial_mask_cell{k,j});
                t1_real = t1_masked(~isnan(t1_masked));
                mean_t1_irse_w(k,j,t) = mean(t1_real);
                sd_t1_irse_w(k,j,t) = std(t1_real);
            end
        end
    end
end

%% IRSE-WE weighted image (OPTIONAL)
if dim == 1
    % mean value of T1 MOLLI
    mean_t1_irsewe_w = zeros(N, length(whatsinit_we));
    sd_t1_irsewe_w = zeros(N, length(whatsinit_we));
    for i = 1:N
        for t = 1:size(mean_t1_irsewe_w, 2)
            t1_masked = nonzeros(whatsinit_we{t} .* roi_cmr_we.vial_mask_cell{i});
            t1_real = t1_masked(~isnan(t1_masked));
            mean_t1_irsewe_w(i, t) = mean(t1_real);
            sd_t1_irsewe_w(i, t) = std(t1_real);
        end
    end
elseif dim == 2
    % mean value of T1 MOLLI
    mean_t1_irsewe_w = zeros(row, col, length(whatsinit_we));
    sd_t1_irsewe_w = zeros(row, col, length(whatsinit_we));
    for k = 1:row
        for j = 1:col
            for t = 1:size(mean_t1_irsewe_w, 3)
                t1_masked = nonzeros(whatsinit_we{t} .* roi_cmr_we.vial_mask_cell{k,j});
                t1_real = t1_masked(~isnan(t1_masked));
                mean_t1_irsewe_w(k,j,t) = mean(t1_real);
                sd_t1_irsewe_w(k,j,t) = std(t1_real);
            end
        end
    end
end
%% Plot the raw IR curves
ff = {'2.5', '5', '10', '20', '30', '40', '100'};
figure();
for j = 1:size(mean_t1_molli_w, 2)
    subplot(3,3,j);
    %for i = 1:size(mean_t1_molli_w, 1)
     
        plot(IR_array_molli, squeeze(mean_t1_molli_w(1, j, :)), 'LineWidth', 2);ylim([0 800]);
        hold on;
        yyaxis right;
        plot(IR_array, squeeze(mean_t1_irse_w(1, j, :)), 'LineWidth', 2);ylim([0 3000]);
        plot(IR_array, squeeze(mean_t1_irsewe_w(1, j, :)), 'LineWidth', 2);
        title(cat(2, 'Fat Fraction: ', ff{j}, '%'));
        legend({'MOLLI', 'IRSE'});
    %end
end


gg = {'2.5', '5', '10', '20', '30', '40', '100'};
figure();
for j = 1:size(mean_t1_molli_w, 2)
    subplot(3,3,j);
    %for i = 1:size(mean_t1_molli_w, 1)
     
        plot(IR_array_molli, squeeze(mean_t1_molli_w(2, j, :)), 'LineWidth', 2);ylim([0 1500]);
        hold on;
        yyaxis right;
        plot(IR_array, squeeze(mean_t1_irse_w(2, j, :)), 'LineWidth', 2);ylim([0 3000]);
        plot(IR_array, squeeze(mean_t1_irsewe_w(2, j, :)), 'LineWidth', 2);
        title(cat(2, 'Fat Fraction: ', gg{j}, '%'));
        legend({'MOLLI', 'IRSE'});
    %end
end
%% Save mean T1, IRSE and IRSE-WE
ir_weighted_save = cat(2, subject_data_dir, 'ir_weighted_metrics.mat');
% load(ir_weighted_save);

ir_weighted_metrics = struct;
ir_weighted_metrics.mean_t1_molli_w = mean_t1_molli_w;
ir_weighted_metrics.mean_t1_irse_w = mean_t1_irse_w;

ir_weighted_metrics.sd_t1_molli_w = sd_t1_molli_w;
ir_weighted_metrics.sd_t1_irse_w = sd_t1_irse_w;

ir_weighted_metrics.IR_array = IR_array;
ir_weighted_metrics.IR_array_molli = IR_array_molli;

if exist('mean_t1_irsewe_w', 'var')
    ir_weighted_metrics.mean_t1_irsewe_w = mean_t1_irsewe_w;
    ir_weighted_metrics.sd_t1_irsewe_w = sd_t1_irsewe_w;
end

save(ir_weighted_save, 'ir_weighted_metrics');