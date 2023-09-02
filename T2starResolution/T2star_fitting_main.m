% Fitting T2* map using qMRLab, 

clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%% input file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DICOM files
% fitting results 
%%%%%%%%%%%%%%%%%%%%%%%%%% output file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FitResults_<avg_name>.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../function/'); 
addpath('../../../qMRLab/');

% Read files
base_dir = uigetdir;
folder_glob = glob(cat(2, base_dir, '/*'));

dicom_fields = {'EchoTime'};
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

avg_num = input('Please type average number here:  ');
if isnumeric(avg_num)
    avg_name = cat(2, 'Avg', num2str(avg_num, '%04.f'));
else
    avg_name = avg_num;
end
%% Read T2* DICOM files
whatsinit = cell(length(list_to_read), 1);
slice_data = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit{i}, slice_data{i}] = dicom23D(f, dicom_fields);
end

% Start qMRinfo and TE array
qMRinfo('mono_t2');

TE_array = zeros(length(slice_data{1,1}), 1);

for i = 1:length(TE_array)
    TE_array(i, 1) = slice_data{1,1}(i).EchoTime;
end

%% Load mask
load(cat(2, subject_data_dir, 'mask.mat'));
FitResults_struct = struct;
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];

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
    if length(whatsinit) == 20
        data.Mask = mask_struct(mask_idx_array(i)).myo_mask;
    else
        data.Mask = mask_struct(i).myo_mask;
    end
    FitResults = FitData(data, Model); %fit data
    % FitResultsSave_mat(FitResults);
    FitResults_struct(i).FitResults = FitResults;
end
%% T2* map (Console generated)
[list_to_read, order_to_read] = NamePicker(folder_glob);
T2star_map = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    T2star_map{i} = dicom23D(f);
end

%% Compare between T2* map and fitted T2* map
i = 1;
if length(whatsinit) == 20
    idx = mask_idx_array(i);
else
    idx = 1;
end

figure();
subplot(2,2,1);
imagesc(T2star_map{i}.*mask_struct(idx).myo_mask); caxis([0 100]);axis image;
colorbar;title('T2* Map from console');
subplot(2,2,2);
imagesc(FitResults_struct(i).FitResults.T2); caxis([0 100]);axis image;colorbar;
title('T2* map from qMRLab')
subplot(2,2,3);
diff_img = abs(T2star_map{i}.*mask_struct(idx).myo_mask - FitResults_struct(i).FitResults.T2);
imagesc(diff_img);caxis([0 20]);axis image;colorbar;
title('Difference Map');
subplot(2,2,4);
imagesc(100*abs(diff_img)./(T2star_map{i}.*mask_struct(idx).myo_mask));caxis([0 50]);axis image;colorbar;
title('Percentage of Difference');

%% Plot residual images
figure('Position', [100 0 1600 1600]);
row = 4;
col = length(whatsinit) / row;
for i = 1:length(whatsinit)
    subplot(row,col,i);
    imagesc(FitResults_struct(i).FitResults.res);
    colorbar;caxis([0 100]);
end

% Save FittedResults
save_f = cat(2, subject_data_dir, 'FitResults_', avg_name, '.mat');
save(save_f, 'FitResults_struct');
