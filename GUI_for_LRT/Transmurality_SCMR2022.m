clear all;
close all;
clc;
current_dir = pwd;
% Dog data configuration for Ting 06/18/2021

addpath('../function/');
addpath('../T1NFF/');
base_dir = uigetdir;

name_glob = glob(cat(2, base_dir, '/*'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'/');
    Names{i} = strings{end-1};
end

% sequence_label = {'MAG'};
% sequence_label = {'DICOM_E', 'DICOM_L', 'PSIR', 'MAG'};
sequence_label = {'LRT_LGE'};
time_label = {'D6', 'D8', 'WK8', 'WK8+2'};

names_to_rule_out = {};
RuleOutLabel = NameRuleOutFunc(Names, names_to_rule_out);
Names = Names(RuleOutLabel == 0);
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
%output_label = {'MAG'};
output_label = {'LRT_LGE'};
%% Draw insertion point (Can be skipped if there is no data added)
for n = 10:length(Names)
% for n = 6:6
    name = Names{n}
    for t = 1:length(time_label)
        %mvo_perc = zeros(length(sequence_label), 1);
        %infarct_perc = zeros(length(sequence_label), 1);
        time = time_label{t}
        dst_dir = cat(2, base_dir, '/', name,'/', name, '_', time_label{t});

        if exist(cat(2, dst_dir, '/ContourData/'))
            save_dir = cat(2, dst_dir, '/Results/');
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end

            for i = 1:length(sequence_label)
                label = sequence_label{i};
                contour_dir = cat(2, dst_dir, '/ContourData/', label, '/');
                load(cat(2, contour_dir, label, '_vol_img_3D.mat'));
                load(cat(2, contour_dir, label, '_SliceLoc.mat'));

                [slc_array_sorted, idx_sorted] = sort(slc_array);
                if ~strcmp(label, 'LRT_LGE')
                    vol_img_tight = zeros(size(vol_img_3D{1},1), size(vol_img_3D{1},2), length(idx_sorted));
                else
                    vol_img_tight = zeros(size(vol_img_3D,1), size(vol_img_3D,2), length(idx_sorted));
                end

                if ~strcmp(label, 'LRT_LGE')
                    for ii = 1:length(vol_img_3D)
                        vol_img_tight(:,:,ii) = vol_img_3D{idx_sorted(ii)};
                    end
                else
                    for ii = 1:size(vol_img_3D,3)
                        vol_img_tight(:,:,ii) = vol_img_3D(:,:,idx_sorted(ii));
                    end
                end

                ref_pts_save = cat(2, contour_dir, '/ref_pts.mat');
                if ~exist(ref_pts_save, 'file')
                    figure();

                    x = zeros(size(vol_img_tight, 3),1);
                    y = zeros(size(vol_img_tight, 3),1);
                    for slc = 1:size(vol_img_tight, 3)
                        imagesc(vol_img_tight(:,:,slc)); % caxis([0 100]);
                        [x(slc),y(slc)] = getpts(gca);
                    end
                    ref_pts.x = x;
                    ref_pts.y = y;
                    save(ref_pts_save, 'ref_pts');
                else
                    load(ref_pts_save);
                    x = ref_pts.x;
                    y = ref_pts.y;
                end
            end
        end
    end
end

%% Transmurality analysis
addpath('../AHA16Segment/');
for n = 1:length(Names)
    % for n = 6:6
    name = Names{n}
    for t = 1:length(time_label)
        %mvo_perc = zeros(length(sequence_label), 1);
        %infarct_perc = zeros(length(sequence_label), 1);
        time = time_label{t}
        dst_dir = cat(2, base_dir, '/', name,'/', name, '_', time_label{t});

        if exist(cat(2, dst_dir, '/ContourData/'))
            save_dir = cat(2, dst_dir, '/Results/');
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end

            for i = 1:length(sequence_label)
                
                label = sequence_label{i};
                contour_dir = cat(2, dst_dir, '/ContourData/', label, '/');
                load(cat(2, contour_dir, label, '_vol_img_3D.mat'));
                load(cat(2, contour_dir, label, '_SliceLoc.mat'));
                load(cat(2, contour_dir, 'Myocardium/mask_myocardium.mat'));
                load(cat(2, contour_dir, 'Heart/mask_heart.mat'));
                load(cat(2, contour_dir, 'freeROI/freeROI.mat'));
                load(cat(2, contour_dir, 'BloodPool/mask_blood.mat'));

                [slc_array_sorted, idx_sorted] = sort(slc_array);
                if ~strcmp(label, 'LRT_LGE')
                    vol_img_tight = zeros(size(vol_img_3D{1},1), size(vol_img_3D{1},2), length(idx_sorted));
                    mask_myocardium_tight = zeros(size(mask_myocardium_3D{1},1), size(mask_myocardium_3D{1},2), length(idx_sorted));
                    mask_heart_tight = zeros(size(mask_heart_3D{1},1), size(mask_heart_3D{1},2), length(idx_sorted));
                    freeROIMask_tight = zeros(size(mask_myocardium_3D{1},1), size(mask_myocardium_3D{1},2), length(idx_sorted));
                    mask_blood_tight = zeros(size(mask_blood_3D{1},1), size(mask_blood_3D{1},2), length(idx_sorted));
                else
                    vol_img_tight = zeros(size(vol_img_3D,1), size(vol_img_3D,2), length(idx_sorted));
                    mask_myocardium_tight = zeros(size(mask_myocardium_3D,1), size(mask_myocardium_3D,2), length(idx_sorted));
                    mask_heart_tight = zeros(size(mask_heart_3D,1), size(mask_heart_3D,2), length(idx_sorted));
                    freeROIMask_tight = zeros(size(mask_myocardium_3D,1), size(mask_myocardium_3D,2), length(idx_sorted));
                    mask_blood_tight = zeros(size(mask_blood_3D,1), size(mask_blood_3D,2), length(idx_sorted));
                end

                if ~strcmp(label, 'LRT_LGE')
                    for ii = 1:length(vol_img_3D)
                        vol_img_tight(:,:,ii) = vol_img_3D{idx_sorted(ii)};
                        mask_myocardium_tight(:,:,ii) = mask_myocardium_3D{idx_sorted(ii)};
                        mask_heart_tight(:,:,ii) = mask_heart_3D{idx_sorted(ii)};
                        freeROIMask_tight(:,:,ii) = freeROIMask_3D{idx_sorted(ii)};
                        mask_blood_tight(:,:,ii) = mask_blood_3D{idx_sorted(ii)};
                    end
                else
                    for ii = 1:size(vol_img_3D, 3)
                        vol_img_tight(:,:,ii) = vol_img_3D(:,:,idx_sorted(ii));
                        mask_myocardium_tight(:,:,ii) = mask_myocardium_3D(:,:,idx_sorted(ii));
                        mask_heart_tight(:,:,ii) = mask_heart_3D(:,:,idx_sorted(ii));
                        freeROIMask_tight(:,:,ii) = freeROIMask_3D(:,:,idx_sorted(ii));
                        mask_blood_tight(:,:,ii) = mask_blood_3D(:,:,idx_sorted(ii));
                    end
                end
                
                freeROIMask_tight = mask_myocardium_tight .* freeROIMask_tight;
                ref_pts_save = cat(2, contour_dir, '/ref_pts.mat');

                load(ref_pts_save);
                x = ref_pts.x;
                y = ref_pts.y;
                
                Segn = 50;
                BaseGroove = zeros(size(vol_img_tight, 3), 1);
                transmurality_array = zeros(size(vol_img_tight, 3), 1);
                
                if any(mask_heart_tight(:))
                    blood_pool_tight = (mask_heart_tight - mask_myocardium_tight)>0;
                else
                    blood_pool_tight = mask_blood_tight;
                end

                for slc = 1:size(vol_img_tight, 3)
                    
                    C = regionprops(blood_pool_tight(:,:,slc));
                    if ~isempty(C)
                        x_centroid = C.Centroid(2);
                        y_centroid = C.Centroid(1);
                        BaseGroove(slc) = atan2(x(slc) - x_centroid, y(slc) - y_centroid) * 180 / pi;

                        [Segmentpix_roi, stats, Mask_Segn] = AHASegmentation(freeROIMask_tight(:,:,slc), mask_myocardium_tight(:,:,slc)>0, Segn, BaseGroove(slc));
                        [Segmentpix_myo, stats, Mask_Segn] = AHASegmentation(vol_img_tight(:,:,slc), mask_myocardium_tight(:,:,slc)>0, Segn, BaseGroove(slc));

                        mi_numel = zeros(Segn, 1);
                        myo_numel = zeros(Segn, 1);
                        for seg = 1:Segn
                            mi_numel(seg) = sum(Segmentpix_roi{seg});
                            myo_numel(seg) = length(Segmentpix_myo{seg});
                        end
                        transmurality_temp = (mi_numel ./ myo_numel)*100;
                        transmurality_temp(transmurality_temp < 1) = [];
                    else
                        transmurality_temp = [];
                    end

                    if isempty(transmurality_temp)
                        transmurality_array(slc) = 0;
                    else
                        transmurality_array(slc) = mean(transmurality_temp);
                    end
                end

                transmurality_array
                transmurality_array_total = mean(transmurality_array(transmurality_array>1))
                pause;
            end
        end
    end
end
%% AHA 16 segment
BaseGroove = zeros(size(t2star_map, 3), 1);
for i = 1:size(t2star_map, 3)
    C = regionprops(myo_mask(:,:,i));
    x_centroid = C.Centroid(2);
    y_centroid = C.Centroid(1);
    BaseGroove(i) = atan2(x(i) - x_centroid, y(i) - y_centroid) * 180 / pi;
end

% tease out margin cases
m = 1;
BaseGroove_truc = BaseGroove((1+m):(end-m));
t2star_map_truc = t2star_map(:,:,(1+m):(end-m));
edg = zeros(size(t2star_map_truc));
myo_mask_truc = myo_mask(:,:,(1+m):(end-m));
for i = 1:size(t2star_map_truc, 3)
    edg(:,:,i)  = edge(squeeze(myo_mask_truc(:,:,i)),'Canny');
end

myo_mask_peel = myo_mask_truc - edg;


n = size(t2star_map_truc, 3);
mode = mod(n,3);
integ = fix(n/3);
if n >= 3
    switch mode
        case {0}
            aha_slice = cat(2, repmat([1], [1, integ]), repmat([2], [1, integ]), repmat([3], [1, integ]));
        case {1}
            aha_slice = cat(2, repmat([1], [1, integ+1]), repmat([2], [1, integ]), repmat([3], [1, integ]));
        case {2}
            aha_slice = cat(2, repmat([1], [1, integ+1]), repmat([2], [1, integ+1]), repmat([3], [1, integ]));
    end
else
    error("Available slice numbers are smaller than 3.");
end