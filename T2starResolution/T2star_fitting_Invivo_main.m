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
dicom_dir = uigetdir;

name_glob = glob(cat(2, dicom_dir, '/*'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'/');
    Names{i} = strings{end-1};
end

time_label = {'8WK'};

names_to_rule_out = {'20P10', '20P40_1Month'};
RuleOutLabel = NameRuleOutFunc(Names, names_to_rule_out);
Names = Names(RuleOutLabel == 0);
name_glob = name_glob(RuleOutLabel == 0);

sequence_label = {'T2STAR'};
label = labels{5};

dicom_fields = {'EchoTime'};

%% Load T2star weighted image
% 17P17, 18P90, 18P92, 18P93, 18P94, 18P95, 20P10, 20P11, 20P40, 20P48
sel_array = [147, 152, 241001, 169, 256001, 171, 269, 263, 201, 200];
% 17P17, 18P90, 18P92, 18P93, 18P94, 18P95, 20P10C, 20P11, 20P40, 20P48
sel_array = [147, 152, 241001, 169, 256001, 171, 266, 263, 201, 200];

whatsinit = cell(length(sel_array), 1);
slice_data = cell(length(sel_array), 8);
count = 0;

for i = 1:length(name_glob)
        count = count + 1;
        idx_array = contains(name_glob, label);

        ya_glob = glob(cat(2, name_glob{i}, 'DICOM/*'));
        dst_names = ExtractNames(ya_glob);
        sub_name = num2str(sel_array(count), '%04.f');

        ind_array = zeros(size(dst_names, 1), 1);
        ind_array = ind_array + contains(dst_names, sub_name);
        ind_array2 = find(ind_array == 1);
        list_to_read = ya_glob(ind_array2);
        yya_glob = glob(cat(2, list_to_read{1}, '*'));

        for j = 1:length(yya_glob)
            % [whatsinit{count}, slice_data{count}] = dicom23D(list_to_read{j});
            whatsinit{count}(:,:,j) = dicomread(yya_glob{j});
            slice_data{count,j} = dicominfo(yya_glob{j});
        end
end

TE_array = zeros(size(slice_data));

for i = 1:size(TE_array,1)
    for j = 1:size(TE_array, 2)
        if ~isempty(slice_data{i,j})
            TE_array(i, j) = slice_data{i,j}.EchoTime;
        end
    end
end

%% Load mask
mask_cell = cell(length(name_glob),1);
for i = 1:length(name_glob)
     load(cat(2, base_dir, '/ContourData_Invivo/', Names{i}, '/',...
        Names{i}, '_', time_label{1}, '/Myocardium/mask_myocardium.mat'));
     mask_cell{i} = mask_myocardium_3D{1};
end
%%

% Start qMRinfo and TE array
qMRinfo('mono_t2');
FitResults_struct = struct;

for i = 1:length(whatsinit)
    % Reshape data and mask
    % Reshape matrix as [Width x Height x #Slice x #TE]
    dicom_size = size(whatsinit{i});
    dicom_reshape = double(reshape(whatsinit{i}, dicom_size(1), dicom_size(2), 1, []));
    %
    Model = mono_t2;  % Create class from model
    %Model = Custom_OptionsGUI(Model);
    Model.Prot.SEdata.Mat = nonzeros(TE_array(i,:)); %
    Model.st = [100 2000];
    Model.lb = [1 2000];
    Model.fx = [0 0];
    Model.voxelwise = 1;
    Model.options.FitType = 'Linear';
    data = struct;  % Create data structure
    data.SEdata = dicom_reshape;
    data.Mask = mask_cell{i};
    FitResults = FitData(data, Model); %fit data
    % FitResultsSave_mat(FitResults);
    FitResults_struct(i).FitResults = FitResults;
end
%% R2 map
for i = 1:length(whatsinit)
    te_array = nonzeros(TE_array(i,:));
    dicom_size = size(whatsinit{i});
    dicom_reshape_2d = double(reshape(whatsinit{i},[],dicom_size(3)));
    S0 = zeros(dicom_size);
    for j = 1:length(te_array)
        S0(:,:,j) = FitResults_struct(i).FitResults.M0 .* exp(-te_array(j) ./ FitResults_struct(i).FitResults.T2);
    end
    S0_reshape_2d = reshape(S0,[],dicom_size(3));
    
    dicom_reshape_2d_mean = mean(dicom_reshape_2d, 2);
    SS_res = sum((dicom_reshape_2d - S0_reshape_2d).^2, 2);
    SS_tot = sum((dicom_reshape_2d - dicom_reshape_2d_mean).^2, 2);
    R2 = reshape(1 - SS_res ./ SS_tot, dicom_size(1), dicom_size(2));
    FitResults_struct(i).FitResults.R2 = R2;
end

% Why R2 is negative???

%% T2* map (Console generated)
sel_array = [501, 153, 242001, 170, 257001, 172, 267, 264, 202, 201];
T2star = cell(length(sel_array), 1);
slice_data_t2star = cell(length(sel_array), 1);
count = 0;

for i = 1:length(name_glob)
        count = count + 1;
        idx_array = contains(name_glob, label);

        ya_glob = glob(cat(2, name_glob{i}, 'DICOM/*'));
        dst_names = ExtractNames(ya_glob);
        sub_name = num2str(sel_array(count), '%04.f');

        ind_array = zeros(size(dst_names, 1), 1);
        ind_array = ind_array + contains(dst_names, sub_name);
        ind_array2 = find(ind_array == 1);
        list_to_read = ya_glob(ind_array2);
        yya_glob = glob(cat(2, list_to_read{1}, '*'));

        for j = 1:length(yya_glob)
            % [whatsinit{count}, slice_data{count}] = dicom23D(list_to_read{j});
            T2star{count} = dicomread(yya_glob{j});
            slice_data_t2star{count} = dicominfo(yya_glob{j});
        end
end

rescaled_t2star = cell(length(sel_array), 1);
for i = 1:length(sel_array)
    rescaled_t2star{i} = double(T2star{i} .* slice_data_t2star{i}.RescaleSlope);
end

%% Compare between T2* map and fitted T2* map
i = 1;

figure();
subplot(2,2,1);
imagesc(rescaled_t2star{i}.*mask_cell{i}); caxis([0 100]);axis image;
colorbar;title('T2* Map from console');
subplot(2,2,2);
imagesc(FitResults_struct(i).FitResults.T2); caxis([0 100]);axis image;colorbar;
title('T2* map from qMRLab')
subplot(2,2,3);
diff_img = abs(rescaled_t2star{i}.*mask_cell{i} - FitResults_struct(i).FitResults.T2);
imagesc(diff_img);caxis([0 20]);axis image;colorbar;
title('Difference Map');
subplot(2,2,4);
imagesc(100*abs(diff_img)./(rescaled_t2star{i}.*mask_cell{i}));caxis([0 50]);axis image;colorbar;
title('Percentage of Difference');

%% Plot residual images
figure('Position', [100 0 1600 1600]);
row = 2;
col = length(whatsinit) / row;
for i = 1:length(whatsinit)
    subplot(row,col,i);
    imagesc(FitResults_struct(i).FitResults.res);
    colorbar;caxis([0 100]);
end
%% Plot R2 images
figure('Position', [100 0 1600 1600]);
row = 2;
col = length(whatsinit) / row;
for i = 1:length(whatsinit)
    subplot(row,col,i);
    imagesc(FitResults_struct(i).FitResults.R2);
    colorbar;caxis([0 1]);
end

%% Save FittedResults
save_dir = cat(2, base_dir, '/data/Invivo_Fitting/');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
save_f = cat(2, save_dir, 'FitResults_Invivo.mat');
save(save_f, 'FitResults_struct');
