clear all;
close all;
clc;
current_dir = pwd;
% Dog data configuration for Ting 06/18/2021

addpath('../function/');
addpath('../T1NFF/');
base_dir = uigetdir; % D:\Data\Ting

name_glob = glob(cat(2, base_dir, '/*'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'\');
    Names{i} = strings{end-1};
end

sequence_label = {'T1'};


names_to_rule_out = {};
RuleOutLabel = NameRuleOutFunc(Names, names_to_rule_out);
Names = Names(RuleOutLabel == 0);

name_check = {'AXEL_Day0'};
name_idx_list = linspace(1, length(Names), length(Names)); % initialize with incremental add

if length(name_check) == 1
    starting_point = find(strcmp(name_check, Names),1);
else
    name_idx_list = zeros(1, length(name_check));
    for n = 1:length(name_check)
        % Check an array of names
        name_idxo = find(strcmp(name_check(n), Names),1);
        name_idx_list(n) = name_idxo;
    end
end


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
    };

t1w_dicom_fields = {...
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

output_label = {'T1'};

%% All are in cell format
time_points = {'Before', '10min', '20min', '30min', '40min', '50min', '60min', '90min'};
Ting_Data_Storage = struct;
for n = starting_point:length(Names)
    name = Names{n};
    label = sequence_label{1};
    contour_dir = cat(2, base_dir, '/', name, '/ContourData/', label, '/' );
    load(cat(2, contour_dir, label, '_vol_img_3D.mat'));
    load(cat(2, contour_dir, 'freeROI/freeROI.mat'));
    load(cat(2, contour_dir, 'MyoReference/myoRef.mat'));
    
    %     figure();
    %     for i = 1:length(vol_img_3D)
    %         subplot(3,3,i);
    %         imagesc(vol_img_3D{i});
    %     end
    roi_cell = {};
    remote_cell = {};
    roi_mean = zeros(length(vol_img_3D), 1);
    roi_sd = zeros(length(vol_img_3D), 1);
    remote_mean = zeros(length(vol_img_3D), 1);
    remote_sd = zeros(length(vol_img_3D), 1);
    for i = 1:length(vol_img_3D)
        roi_cell{i} = nonzeros(vol_img_3D{i} .* freeROIMask_3D{i});
        remote_cell{i} = nonzeros(vol_img_3D{i} .* myoRefMask_3D{i});
        roi_mean(i) = mean(roi_cell{i});
        remote_mean(i) = mean(remote_cell{i});
        roi_sd(i) = std(roi_cell{i});
        remote_sd(i) = std(remote_cell{i});
    end
    
    x = 1:1:length(time_points);
    figure();
    plotHandles(:,1) = errorbar(x, roi_mean, roi_sd, 'LineWidth', 2);
    hold on;
    plotHandles(:,2) = errorbar(x, remote_mean, remote_sd, 'LineWidth', 2);
    grid on;
    
    Ting_Data_Storage.Name = name;
    Ting_Data_Storage.metrics = struct;
    Ting_Data_Storage.metrics.roi_cell = roi_cell;
    Ting_Data_Storage.metrics.remote_cell = remote_cell;
    Ting_Data_Storage.metrics.roi_mean = roi_mean;
    Ting_Data_Storage.metrics.remote_mean = remote_mean;
    Ting_Data_Storage.metrics.roi_sd = roi_sd;
    Ting_Data_Storage.metrics.remote_sd = remote_sd;
end

save_dir = cat(2, base_dir, '/', name, '/Results/');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

save(cat(2, save_dir, label, '_data_storage.mat'), 'Ting_Data_Storage');

