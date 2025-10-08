clear all;
close all;
%% Load target dir and check files
addpath('/Users/jameszhang/Documents/MATLAB/PhD_Project_Code_Storage/function/');
addpath('/Users/jameszhang/Documents/MATLAB/PhD_Project_Code_Storage/AHA16Segment/');

% Write a base directory
baseDir = '/Users/jameszhang/Documents/Data/LIUXIN_Project/';

contourDir = fullfile(baseDir, 'ContourData');
folder_glob = glob(cat(2, contourDir, '/*'));

Names = ExtractNames(folder_glob);

hasLGE = false(length(Names),1);
hasBOOST2 = false(length(Names),1);
hasT2 = false(length(Names),1);

for i = 1:length(Names)
    subfolder = fullfile(contourDir, Names{i});
    items = dir(subfolder);
    itemNames = {items([items.isdir]).name};
    hasLGE(i) = any(strcmp(itemNames, 'LGE'));
    hasBOOST2(i) = any(strcmp(itemNames, 'BOOST2'));
    hasT2(i) = any(strcmp(itemNames, 'T2'));
end

%% 

target_subject_label = hasLGE & hasBOOST2 & hasT2;
target_subjects = Names(target_subject_label);


%for i = 1:length(target_subjects)
for i = 1:1
    subfolder = fullfile(contourDir, target_subjects{i});
    items = dir(subfolder);
    itemNames = {items([items.isdir]).name};
    
    if any(strcmp(itemNames, 'LGE')) && any(strcmp(itemNames, 'BOOST2')) && any(strcmp(itemNames, 'T2'))
        fprintf('Processing subject: %s\n', target_subjects{i});
        % Call the processing function here
        % Example: ProcessSubject(subfolder);
        % name = target_subjects{i};
        load(fullfile(subfolder, 'LGE', 'Myocardium', 'mask_myocardium.mat'));
        load(fullfile(subfolder, 'LGE', 'freeROI', 'freeROI.mat'));
        load(fullfile(subfolder, 'LGE', 'LGE_SliceLoc.mat'));
        load(fullfile(subfolder, 'LGE', 'LGE_vol_img_3D.mat'));
        load(fullfile(subfolder, 'LGE', 'BloodPool', 'mask_blood.mat'));
        freeROIMask_3D_lge = freeROIMask_3D;
        mask_myocardium_3D_lge = mask_myocardium_3D;
        slc_array_lge = slc_array;
        vol_img_3D_lge = vol_img_3D;
        mask_blood_3D_lge = mask_blood_3D;


        load(fullfile(subfolder, 'BOOST2', 'Myocardium', 'mask_myocardium.mat'));
        load(fullfile(subfolder, 'BOOST2', 'freeROI', 'freeROI.mat'));
        load(fullfile(subfolder, 'BOOST2', 'BOOST2_SliceLoc.mat'));
        load(fullfile(subfolder, 'BOOST2', 'BOOST2_vol_img_3D.mat'));
        load(fullfile(subfolder, 'BOOST2', 'BloodPool', 'mask_blood.mat'));
        load(fullfile(subfolder, 'BOOST2', 'MyoReference', 'myoRef.mat'));
        load(fullfile(subfolder, 'BOOST2', 'noReflowArea', 'noReflow.mat'));
        load(fullfile(subfolder, 'BOOST2', 'excludeArea', 'excludeArea.mat'));
        freeROIMask_3D_boost2 = freeROIMask_3D;
        mask_myocardium_3D_boost2 = mask_myocardium_3D;
        slc_array_boost2 = slc_array;
        vol_img_3D_boost2 = vol_img_3D;
        mask_blood_3D_boost2 = mask_blood_3D;

        load(fullfile(subfolder, 'T2', 'Myocardium', 'mask_myocardium.mat'));
        load(fullfile(subfolder, 'T2', 'freeROI', 'freeROI.mat'));
        load(fullfile(subfolder, 'T2', 'T2_SliceLoc.mat'));
        load(fullfile(subfolder, 'T2', 'T2_vol_img_3D.mat'));
        load(fullfile(subfolder, 'T2', 'BloodPool', 'mask_blood.mat'));
        freeROIMask_3D_t2 = freeROIMask_3D;
        mask_myocardium_3D_t2 = mask_myocardium_3D;
        slc_array_t2 = slc_array;
        vol_img_3D_t2 = vol_img_3D;
        mask_blood_3D_t2 = mask_blood_3D;

        % Check if the masks are consistent across modalitie
        ub_lge  = 0.9*max(vol_img_3D_lge(:));
        ub_boost2 = 0.9*max(vol_img_3D_boost2(:));
        ub_t2 = 0.9*max(vol_img_3D_t2(:));
        lb_lge = 1.4*min(vol_img_3D_lge(:));
        lb_boost2 = 1.4*min(vol_img_3D_boost2(:));
        lb_t2 = 1.4*min(vol_img_3D_t2(:));
        
        fname_lge = fullfile(subfolder, 'LGE', 'Groove.mat');
        fname_boost2 = fullfile(subfolder, 'BOOST2', 'Groove.mat');
        fname_t2 = fullfile(subfolder, 'T2', 'Groove.mat');
        % disp('Draw Groove points for each LGE slice:');
        if ~exist(fullfile(subfolder, 'LGE', 'Groove.mat'), 'file')
            fprintf('Drawing Groove points for LGE...\n');
            coords_lge = GroovePickGeneric(vol_img_3D_lge, mask_blood_3D_lge, lb_lge, ub_lge);
            save(fname_lge, 'coords_lge');
        else
            fprintf('LGE Groove points already exist. Skipping...\n');
            load(fname_lge, 'coords_lge');
        end
        if ~exist(fullfile(subfolder, 'BOOST2', 'Groove.mat'), 'file')
            fprintf('Drawing Groove points for BOOST2...\n');
            coords_boost2 = GroovePickGeneric(vol_img_3D_boost2, mask_blood_3D_boost2, lb_boost2, ub_boost2);
            save(fname_boost2, 'coords_boost2');
        else
            fprintf('BOOST2 Groove points already exist. Skipping...\n');
            load(fname_boost2, 'coords_boost2');
        end
        if ~exist(fullfile(subfolder, 'T2', 'Groove.mat'), 'file')
            fprintf('Drawing Groove points for T2...\n');
            coords_t2 = GroovePickGeneric(vol_img_3D_t2, mask_blood_3D_t2, lb_boost2, ub_t2);
            save(fname_t2, 'coords_t2');
        else
            fprintf('T2 Groove points already exist. Skipping...\n');
            load(fname_t2, 'coords_t2');
        end


        x = coords_lge.x_array;
        y = coords_lge.y_array;
        x_centroid = coords_lge.x_centroid_array;
        y_centroid = coords_lge.y_centroid_array;

        BaseGroove_lge = zeros(size(vol_img_3D_lge,3), 1);
        for i = 1:size(vol_img_3D_lge,3)
            if ~isnan(x(i))
                BaseGroove_lge(i) = atan2(x(i) - x_centroid(i), y(i) - y_centroid(i)) * 180 / pi;
            end
        end

        x = coords_boost2.x_array;
        y = coords_boost2.y_array;
        x_centroid = coords_boost2.x_centroid_array;
        y_centroid = coords_boost2.y_centroid_array;

        BaseGroove_boost2 = zeros(size(vol_img_3D_boost2,3), 1);
        for i = 1:size(vol_img_3D_boost2,3)
            if ~isnan(x(i))
                BaseGroove_boost2(i) = atan2(x(i) - x_centroid(i), y(i) - y_centroid(i)) * 180 / pi;
            end
        end

        x = coords_t2.x_array;
        y = coords_t2.y_array;
        x_centroid = coords_t2.x_centroid_array;
        y_centroid = coords_t2.y_centroid_array;

        BaseGroove_t2 = zeros(size(vol_img_3D_t2,3), 1);
        for i = 1:size(vol_img_3D_t2,3)
            if ~isnan(x(i))
                BaseGroove_t2(i) = atan2(x(i) - x_centroid(i), y(i) - y_centroid(i)) * 180 / pi;
            end
        end


        % keyboard;
        % Check if it is apex -> base

        % vol_img_3D_lge = vol_img_3D_lge(:,:,idx_array_lge);
        % freeROIMask_3D_lge = freeROIMask_3D_lge(:,:,idx_array_lge);
        % mask_myocardium_3D_lge = mask_myocardium_3D_lge(:,:,idx_array_lge);
        % BaseGroove = BaseGroove(idx_array_lge);
        % mask_blood_3D_lge  = mask_blood_3D_lge (:,:,idx_array_lge);
        % x_centroid_lge = x_centroid(idx_array_lge);
        % y_centroid_lge = y_centroid(idx_array_lge);
        freeROI_check = sum(reshape(freeROIMask_3D_lge, [], size(freeROIMask_3D_lge,3)),1) > 0;

        [~, idx_array_lge] = sort(slc_array_lge(freeROI_check));
        [~, idx_array_boost2] = sort(slc_array_boost2(freeROI_check));
        [~, idx_array_t2] = sort(slc_array_t2(freeROI_check));

        LocPixCount_lge = AHA_16Seg(freeROIMask_3D_lge(:,:,freeROI_check).*mask_myocardium_3D_lge(:,:,freeROI_check)>0, ...
            mask_myocardium_3D_lge(:,:,freeROI_check)>0, BaseGroove_lge(freeROI_check), flip(idx_array_lge));


        LocPixCount_t2 = AHA_16Seg(freeROIMask_3D_t2(:,:,freeROI_check).*mask_myocardium_3D_t2(:,:,freeROI_check)>0, ...
            mask_myocardium_3D_t2(:,:,freeROI_check)>0, BaseGroove_t2(freeROI_check), flip(idx_array_t2));

        aha_gt_lge = LocPixCount_lge(:,1) > 0.01;
        aha_gt_t2 = LocPixCount_t2(:,1) > 0.01;

        mean_boost2 = mean(nonzeros(vol_img_3D_boost2(:,:,freeROI_check) .* myoRefMask_3D(:,:,freeROI_check)));
        sd_boost2 = std(nonzeros(vol_img_3D_boost2(:,:,freeROI_check) .* myoRefMask_3D(:,:,freeROI_check)));

        thresh_array = mean_boost2 + [2:8] * sd_boost2;
        LocPixCount_boost2 = zeros(16, 2, length(thresh_array));
        ROIMask_3D_boost2 = zeros(size(vol_img_3D_boost2,1), size(vol_img_3D_boost2,2), size(vol_img_3D_boost2,3), length(thresh_array));
        for v = 1:length(thresh_array)
            ROIMask_3D_boost2(:,:,:,v) = vol_img_3D_boost2 > thresh_array(v);
            ROIMask_3D_boost2(:,:,:,v) = (ROIMask_3D_boost2(:,:,:,v) .* mask_myocardium_3D_boost2 .* ~excludeMask_3D > 0 + ...
                noReflowMask_3D > 0) > 0;
            LocPixCount_boost2(:,:,v) = AHA_16Seg(ROIMask_3D_boost2(:,:,freeROI_check,v), mask_myocardium_3D_boost2(:,:,freeROI_check)>0, ...
                BaseGroove_boost2(freeROI_check), flip(idx_array_boost2));
        end

        LocPixCount_boost2_aha = LocPixCount_boost2 > 0.01;

        figure;
        numSlices = size(vol_img_3D_lge, 3);
        for s = 1:numSlices
            subplot(ceil(sqrt(numSlices)), ceil(sqrt(numSlices)), s);
            imagesc(vol_img_3D_lge(:,:,s) .* mask_myocardium_3D_lge(:,:,s));
            axis image off;
            title(sprintf('Slice %d', s));
            colormap gray;
        end
        sgtitle('LGE Volume Slices');

        figure;
        numSlices = size(vol_img_3D_boost2, 3);
        for s = 1:numSlices
            subplot(ceil(sqrt(numSlices)), ceil(sqrt(numSlices)), s);
            imagesc(vol_img_3D_boost2(:,:,s) .* mask_myocardium_3D_boost2(:,:,s));
            axis image off;
            title(sprintf('Slice %d', s));
            colormap gray;
        end
        sgtitle('BOOST2 Volume Slices');

        figure;
        numSlices = size(vol_img_3D_lge, 3);
        for s = 1:numSlices
            subplot(ceil(sqrt(numSlices)), ceil(sqrt(numSlices)), s);
            imagesc(vol_img_3D_lge(:,:,s) .* freeROIMask_3D_lge(:,:,s));
            axis image off;
            title(sprintf('Slice %d', s));
            colormap gray;
        end
        sgtitle('LGE Volume Slices Masked by freeROI');
    else
        fprintf('Skipping subject: %s due to missing data.\n', target_subjects{i});
    end
end



%% 
centroid_mask_lge = zeros(1, 2); % Initialize centroid array
for s = 1:size(mask_myocardium_3D_lge, 3)
    [y_indices, x_indices] = find(mask_myocardium_3D_lge(:,:,s)); % Get indices of non-zero elements
    if ~isempty(x_indices) && ~isempty(y_indices)
        centroid_mask_lge(1) = mean(x_indices); % Calculate mean x-coordinate
        centroid_mask_lge(2) = mean(y_indices); % Calculate mean y-coordinate

        centroid_masks_lge(s, :) = centroid_mask_lge; % Store centroid for the current slice

       
    end
end

% Find slices where centroid_masks_lge is smaller than 32
small_centroid_indices = find(centroid_masks_lge(:, 1) < 32 & centroid_masks_lge(:, 2) < 32);
if ~isempty(small_centroid_indices)
    fprintf('Slices with centroids smaller than 32 found at indices: %s\n', num2str(small_centroid_indices'));
else
    fprintf('No slices with centroids smaller than 32 found.\n');
end

centroid_masks_lge = zeros(size(mask_myocardium_3D_lge, 3), 2); % Initialize array for centroids
for s = 1:size(mask_myocardium_3D_lge, 3)
    centroid_mask_lge = zeros(1, 2); % Initialize centroid array for current slice
    [y_indices, x_indices] = find(mask_myocardium_3D_lge(:,:,s)); % Get indices of non-zero elements
    if ~isempty(x_indices) && ~isempty(y_indices)
        centroid_mask_lge(1) = mean(x_indices); % Calculate mean x-coordinate
        centroid_mask_lge(2) = mean(y_indices); % Calculate mean y-coordinate

        centroid_masks_lge(s, :) = centroid_mask_lge; % Store centroid for the current slice

        % Multiply coordinates by 8 if the slice index is in small_centroid_indices
        if ismember(s, small_centroid_indices)
            centroid_masks_lge(s, :) = centroid_masks_lge(s, :) * 8;
            % Multiply y_indices and x_indices by 8
            x_indices = x_indices * 8;
            y_indices = y_indices * 8;

            % Create a new mask for the modified indices
            new_mask = zeros(size(mask_myocardium_3D_lge(:,:,s)));
            new_mask(sub2ind(size(new_mask), y_indices, x_indices)) = 1;

            % Update the original mask_myocardium_3D_lge with the new mask
            mask_myocardium_3D_lge(:,:,s) = new_mask;
        end
    end
end