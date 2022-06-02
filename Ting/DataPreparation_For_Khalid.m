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
%name_check = {'Evelyn_Day0'};
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
%% Endo-Epi
time_points = {'Baseline', '10min', '20min', '30min', '40min', '50min', '60min', '90min'}; % AXEL
%time_points = {'Baseline', '15min', '20min', '30min', '35min', '55min', 'Occlu_1hr', 'Occlu_2hr'}; % Evelyn

data = struct;
%[optimizer, metric] = imregconfig('multimodal');
[optimizer, metric] = imregconfig('monomodal');
fix_slc = 2;
for n = starting_point:starting_point
    name = Names{n};
    label = sequence_label{1};
    contour_dir = cat(2, base_dir, '/', name, '/ContourData/', label, '/' );
    load(cat(2, contour_dir, label, '_vol_img_3D.mat'));
    load(cat(2, contour_dir, 'freeROI/freeROI.mat'));
    load(cat(2, contour_dir, 'Myocardium/mask_myocardium.mat'));
    
    fixed = mask_myocardium_3D{fix_slc};
    myocardium_3D = zeros(size(fixed,1), size(fixed,2), length(mask_myocardium_3D));
    t1map_3D = zeros(size(fixed,1), size(fixed,2), length(mask_myocardium_3D));
    roi_3D = zeros(size(fixed,1), size(fixed,2), length(mask_myocardium_3D));
    
    
        
    for i = 1:length(mask_myocardium_3D)
        if any(strcmp(name, {'Evelyn_Day0', 'AXEL_Day0'}))
            wy0 = size(vol_img_3D{2},2);
            wyi = size(vol_img_3D{i},2);
            vol_img_3D{i} = vol_img_3D{i}(:,((wyi-wy0)+1)/2:(end-(wyi-wy0)/2));
            mask_myocardium_3D{i} = mask_myocardium_3D{i}(:,((wyi-wy0)+1)/2:(end-(wyi-wy0)/2));
            freeROIMask_3D{i} = freeROIMask_3D{i}(:,((wyi-wy0)+1)/2:(end-(wyi-wy0)/2));
        end
        moving = mask_myocardium_3D{i};
        moving_t1 = vol_img_3D{i};
        moving_roi = freeROIMask_3D{i};
        if i ~= fix_slc
            % tform = imregtform(moving, fixed, 'rigid', optimizer, metric);
            tform = imregtform(moving, fixed, 'affine', optimizer, metric);
            % tform = imregtform(moving_t1, vol_img_3D{fix_slc}, 'affine', optimizer, metric);
            
            movingRegistered_myo = imwarp(moving,tform,'OutputView',imref2d(size(fixed)))>0.5;
            movingRegistered_t1 = imwarp(moving_t1,tform,'OutputView',imref2d(size(fixed)));
            movingRegistered_roi = imwarp(moving_roi,tform,'OutputView',imref2d(size(fixed)));
            myocardium_3D(:,:,i) = movingRegistered_myo;
            t1map_3D(:,:,i) = movingRegistered_t1;
            roi_3D(:,:,i) = movingRegistered_roi;
        else
            myocardium_3D(:,:,i) = fixed;
            t1map_3D(:,:,i) = moving_t1;
            roi_3D(:,:,i) = moving_roi;
        end
        
    end
    
    data.myocardium_3D = myocardium_3D;
    data.t1map_3D = t1map_3D;
    data.roi_3D = roi_3D;
    data.time_points = time_points;
    
    save_dir = cat(2, base_dir, '/', name, '/Results/');
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    
    save(cat(2, save_dir, name, '_affine_mono.mat'), '-struct', 'data');
end
    
%% Test how well are they registered

