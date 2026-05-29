clear all;
close all;
%% Load target dir and check files
addpath('/Users/jameszhang/Documents/MATLAB/PhD_Project_Code_Storage/function/');
addpath('/Users/jameszhang/Documents/MATLAB/PhD_Project_Code_Storage/AHA16Segment/');

% Write a base directory
baseDir = '/Volumes/SSD_Data/LIUXIN_Project/';
saveDir = '/Volumes/SSD_Data/LIUXIN_Project/Results/';
if ~exist(saveDir)
    mkdir(saveDir);
end
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

cutoff_array = [0.01, 0.05, 0.1 0.2 0.3 0.4 0.5];

%for i = 1:length(target_subjects)
ROC_struct = struct;
ROC_struct_persub = struct;
ROC_struct_persub_v2 = struct;
for nn = 1:length(target_subjects)
    subfolder = fullfile(contourDir, target_subjects{nn});
    items = dir(subfolder);
    itemNames = {items([items.isdir]).name};
    
    saveSubDir = fullfile(saveDir, target_subjects{nn});
    if ~exist(saveSubDir)
        mkdir(saveSubDir);
    end

    if any(strcmp(itemNames, 'LGE')) && any(strcmp(itemNames, 'BOOST2')) && any(strcmp(itemNames, 'T2'))
        fprintf('Processing subject: %s\n', target_subjects{nn});

        ROC_struct.(target_subjects{nn}) = struct;

        % Call the processing function here
        % Example: ProcessSubject(subfolder);
        % name = target_subjects{i};
        load(fullfile(subfolder, 'LGE', 'Myocardium', 'mask_myocardium.mat'));
        load(fullfile(subfolder, 'LGE', 'freeROI', 'freeROI.mat'));
        load(fullfile(subfolder, 'LGE', 'LGE_SliceLoc.mat'));
        load(fullfile(subfolder, 'LGE', 'LGE_vol_img_3D.mat'));
        load(fullfile(subfolder, 'LGE', 'BloodPool', 'mask_blood.mat'));
        load(fullfile(subfolder, 'LGE', 'LGE_refPoints.mat'));
        freeROIMask_3D_lge = freeROIMask_3D;
        mask_myocardium_3D_lge = mask_myocardium_3D;
        slc_array_lge = slc_array;
        vol_img_3D_lge = vol_img_3D;
        mask_blood_3D_lge = mask_blood_3D;

        x_array_lge = zeros(length(refPoints.refPoint), 1);
        y_array_lge = zeros(length(refPoints.refPoint), 1);
        for ref_p = 1:length(refPoints.refPoint)
            if ~isempty(refPoints.refPoint{ref_p})
                x_array_lge(ref_p) = refPoints.refPoint{ref_p}(1);
                y_array_lge(ref_p) = refPoints.refPoint{ref_p}(2);
            end
        end
        x_array_lge(x_array_lge == 0) = nan;
        y_array_lge(y_array_lge == 0) = nan;

        load(fullfile(subfolder, 'BOOST2', 'Myocardium', 'mask_myocardium.mat'));
        load(fullfile(subfolder, 'BOOST2', 'freeROI', 'freeROI.mat'));
        load(fullfile(subfolder, 'BOOST2', 'BOOST2_SliceLoc.mat'));
        load(fullfile(subfolder, 'BOOST2', 'BOOST2_vol_img_3D.mat'));
        load(fullfile(subfolder, 'BOOST2', 'BloodPool', 'mask_blood.mat'));
        load(fullfile(subfolder, 'BOOST2', 'MyoReference', 'myoRef.mat'));
        load(fullfile(subfolder, 'BOOST2', 'noReflowArea', 'noReflow.mat'));
        load(fullfile(subfolder, 'BOOST2', 'excludeArea', 'excludeArea.mat'));
        load(fullfile(subfolder, 'BOOST2', 'BOOST2_refPoints.mat'));
        freeROIMask_3D_boost2 = freeROIMask_3D;
        mask_myocardium_3D_boost2 = mask_myocardium_3D;
        slc_array_boost2 = slc_array;
        vol_img_3D_boost2 = vol_img_3D;
        mask_blood_3D_boost2 = mask_blood_3D;

        x_array_boost2 = zeros(length(refPoints.refPoint), 1);
        y_array_boost2 = zeros(length(refPoints.refPoint), 1);
        for ref_p = 1:length(refPoints.refPoint)
            if ~isempty(refPoints.refPoint{ref_p})
                x_array_boost2(ref_p) = refPoints.refPoint{ref_p}(1);
                y_array_boost2(ref_p) = refPoints.refPoint{ref_p}(2);
            end
        end
        x_array_boost2(x_array_boost2 == 0) = nan;
        y_array_boost2(y_array_boost2 == 0) = nan;

        load(fullfile(subfolder, 'T2', 'Myocardium', 'mask_myocardium.mat'));
        load(fullfile(subfolder, 'T2', 'freeROI', 'freeROI.mat'));
        load(fullfile(subfolder, 'T2', 'T2_SliceLoc.mat'));
        load(fullfile(subfolder, 'T2', 'T2_vol_img_3D.mat'));
        load(fullfile(subfolder, 'T2', 'BloodPool', 'mask_blood.mat'));
        load(fullfile(subfolder, 'T2', 'T2_refPoints.mat'));
        freeROIMask_3D_t2 = freeROIMask_3D;
        mask_myocardium_3D_t2 = mask_myocardium_3D;
        slc_array_t2 = slc_array;
        vol_img_3D_t2 = vol_img_3D;
        mask_blood_3D_t2 = mask_blood_3D;

        x_array_t2 = zeros(length(refPoints.refPoint), 1);
        y_array_t2 = zeros(length(refPoints.refPoint), 1);
        for ref_p = 1:length(refPoints.refPoint)
            if ~isempty(refPoints.refPoint{ref_p})
                x_array_t2(ref_p) = refPoints.refPoint{ref_p}(1);
                y_array_t2(ref_p) = refPoints.refPoint{ref_p}(2);
            end
        end

        x_array_t2(x_array_t2 == 0) = nan;
        y_array_t2(y_array_t2 == 0) = nan;

        %
        % Fill NaNs by nearest non-NaN neighbour along the slice dimension
        vars = {'x_array_lge','y_array_lge','x_array_boost2','y_array_boost2','x_array_t2','y_array_t2'};
        for v = 1:numel(vars)
            arr = eval(vars{v});
            n = numel(arr);
            if n == 0
                continue;
            end
            nanIdx = isnan(arr);
            if any(nanIdx)
                goodIdx = find(~nanIdx);
                if isempty(goodIdx)
                    warning('All elements are NaN for %s; leaving as-is.', vars{v});
                    continue;
                end
                % Fill NaNs using nearest neighbour interpolation (with extrapolation)
                arr(nanIdx) = interp1(goodIdx, arr(goodIdx), find(nanIdx), 'nearest', 'extrap');
                eval([vars{v} ' = arr;']);
                fprintf('Filled %d NaNs in %s\n', sum(nanIdx), vars{v});
            end
        end

        
        % Check if the masks are consistent across modalitie
        ub_lge  = 0.9*max(vol_img_3D_lge(:));
        ub_boost2 = 0.9*max(vol_img_3D_boost2(:));
        ub_t2 = 0.9*max(vol_img_3D_t2(:));
        lb_lge = 1.4*min(vol_img_3D_lge(:));
        lb_boost2 = 1.4*min(vol_img_3D_boost2(:));
        lb_t2 = 1.4*min(vol_img_3D_t2(:));
        
        % Deprecated 11/03/2025 XZ
        % fname_lge = fullfile(subfolder, 'LGE', 'Groove.mat');
        % fname_boost2 = fullfile(subfolder, 'BOOST2', 'Groove.mat');
        % fname_t2 = fullfile(subfolder, 'T2', 'Groove.mat');
        % % disp('Draw Groove points for each LGE slice:');
        % if ~exist(fullfile(subfolder, 'LGE', 'Groove.mat'), 'file')
        %     fprintf('Drawing Groove points for LGE...\n');
        %     coords_lge = GroovePickGeneric(vol_img_3D_lge, mask_blood_3D_lge, lb_lge, ub_lge);
        %     save(fname_lge, 'coords_lge');
        % else
        %     fprintf('LGE Groove points already exist. Skipping...\n');
        %     load(fname_lge, 'coords_lge');
        % end
        % if ~exist(fullfile(subfolder, 'BOOST2', 'Groove.mat'), 'file')
        %     fprintf('Drawing Groove points for BOOST2...\n');
        %     coords_boost2 = GroovePickGeneric(vol_img_3D_boost2, mask_blood_3D_boost2, lb_boost2, ub_boost2);
        %     save(fname_boost2, 'coords_boost2');
        % else
        %     fprintf('BOOST2 Groove points already exist. Skipping...\n');
        %     load(fname_boost2, 'coords_boost2');
        % end
        % if ~exist(fullfile(subfolder, 'T2', 'Groove.mat'), 'file')
        %     fprintf('Drawing Groove points for T2...\n');
        %     coords_t2 = GroovePickGeneric(vol_img_3D_t2, mask_blood_3D_t2, lb_boost2, ub_t2);
        %     save(fname_t2, 'coords_t2');
        % else
        %     fprintf('T2 Groove points already exist. Skipping...\n');
        %     load(fname_t2, 'coords_t2');
        % end
        
        % Get centroid of three images
        x_centroid_array_lge = zeros(size(mask_blood_3D_lge,3),1);
        y_centroid_array_lge = zeros(size(mask_blood_3D_lge,3),1);
        for slice_num = 1:size(mask_blood_3D_lge,3)

            [x_heart, y_heart] = find(mask_blood_3D_lge(:,:,slice_num) ~= 0);

            x_centroid_array_lge(slice_num) = round(mean(x_heart),1);
            y_centroid_array_lge(slice_num) = round(mean(y_heart),1);
        end

        x_centroid_array_t2 = zeros(size(mask_blood_3D_t2,3),1);
        y_centroid_array_t2 = zeros(size(mask_blood_3D_t2,3),1);
        for slice_num = 1:size(mask_blood_3D_t2,3)

            [x_heart, y_heart] = find(mask_blood_3D_t2(:,:,slice_num) ~= 0);

            x_centroid_array_t2(slice_num) = round(mean(x_heart),1);
            y_centroid_array_t2(slice_num) = round(mean(y_heart),1);
        end

        x_centroid_array_boost2 = zeros(size(mask_blood_3D_boost2,3),1);
        y_centroid_array_boost2 = zeros(size(mask_blood_3D_boost2,3),1);
        for slice_num = 1:size(mask_blood_3D_boost2,3)

            [x_heart, y_heart] = find(mask_blood_3D_boost2(:,:,slice_num) ~= 0);

            x_centroid_array_boost2(slice_num) = round(mean(x_heart),1);
            y_centroid_array_boost2(slice_num) = round(mean(y_heart),1);
        end


        [x_array_lge, y_array_lge] = check_coords(vol_img_3D_lge, x_array_lge, y_array_lge);
        [x_array_boost2, y_array_boost2] = check_coords(vol_img_3D_boost2, x_array_boost2, y_array_boost2);
        [x_array_t2, y_array_t2] = check_coords(vol_img_3D_t2, x_array_t2, y_array_t2);


        x = x_array_lge;
        y = y_array_lge;
        x_centroid = x_centroid_array_lge;
        y_centroid = y_centroid_array_lge;

        % To visualize the LGE slices
        % Overlay coordinates on each LGE slice
        % figure;
        % numSlices = size(vol_img_3D_lge, 3);
        % for s = 1:numSlices
        %     imagesc(vol_img_3D_lge(:,:,s)); colormap gray; axis image off; hold on;
        %     if s <= numel(x) && ~isnan(x(s)) && ~isnan(y(s))
        %         % plot marker (y,x) because imagesc displays columns on x-axis and rows on y-axis
        %         plot(y(s), x(s), 'r+','MarkerSize',12,'LineWidth',2);
        %         % small text label next to the marker
        %         txt = sprintf('x=%.1f, y=%.1f', x(s), y(s));
        %         text(y(s)+3, x(s), txt, 'Color','y','FontSize',9, 'FontWeight','bold');
        % 
        %         plot(y_centroid(s), x_centroid(s), 'g+','MarkerSize',12,'LineWidth',2);
        %     end
        %     title(sprintf('LGE Slice %d', s));
        %     hold off;
        %     drawnow;
        %     pause(1); % adjust pause or remove for faster display
        % end


        BaseGroove_lge = zeros(size(vol_img_3D_lge,3), 1);
        for i = 1:size(vol_img_3D_lge,3)
            if ~isnan(x(i))
                BaseGroove_lge(i) = atan2(x(i) - x_centroid(i), y(i) - y_centroid(i)) * 180 / pi;
            end
        end

        x = x_array_boost2;
        y = y_array_boost2;
        x_centroid = x_centroid_array_boost2;
        y_centroid = y_centroid_array_boost2;

        BaseGroove_boost2 = zeros(size(vol_img_3D_boost2,3), 1);
        for i = 1:size(vol_img_3D_boost2,3)
            if ~isnan(x(i))
                BaseGroove_boost2(i) = atan2(x(i) - x_centroid(i), y(i) - y_centroid(i)) * 180 / pi;
            end
        end

        x = x_array_t2;
        y = y_array_t2;
        x_centroid = x_centroid_array_t2;
        y_centroid = y_centroid_array_t2;

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

        if sum(freeROI_check) > 0
            myo_check = sum(reshape(mask_myocardium_3D_boost2, [], size(mask_myocardium_3D_boost2,3)),1) > 0;

            [~, idx_array_lge] = sort(slc_array_lge(freeROI_check));
            [~, idx_array_boost2] = sort(abs(slc_array_boost2(myo_check)));
            % idx_array_boost2 = idx_array_boost2(freeROI_check);
            [~, idx_array_t2] = sort(slc_array_t2(freeROI_check));

            LocPixCount_lge = AHA_16Seg(freeROIMask_3D_lge(:,:,freeROI_check).*mask_myocardium_3D_lge(:,:,freeROI_check)>0, ...
                mask_myocardium_3D_lge(:,:,freeROI_check)>0, BaseGroove_lge(freeROI_check), flip(idx_array_lge));


            LocPixCount_t2 = AHA_16Seg(freeROIMask_3D_t2(:,:,freeROI_check).*mask_myocardium_3D_t2(:,:,freeROI_check)>0, ...
                mask_myocardium_3D_t2(:,:,freeROI_check)>0, BaseGroove_t2(freeROI_check), flip(idx_array_t2));



            mean_boost2 = mean(nonzeros(vol_img_3D_boost2(:,:,myo_check) .* myoRefMask_3D(:,:,myo_check)));
            sd_boost2 = std(nonzeros(vol_img_3D_boost2(:,:,myo_check) .* myoRefMask_3D(:,:,myo_check)));

            thresh_array = mean_boost2 + [1:8] * sd_boost2;
            LocPixCount_boost2 = zeros(16, 2, length(thresh_array));
            LocPixCount_boost2_v2 = zeros(16, 2, length(thresh_array));
            ROIMask_3D_boost2 = zeros(size(vol_img_3D_boost2,1), size(vol_img_3D_boost2,2), size(vol_img_3D_boost2,3), length(thresh_array));
            for v = 1:length(thresh_array)
                ROIMask_3D_boost2(:,:,:,v) = vol_img_3D_boost2 > thresh_array(v);
                ROIMask_3D_boost2(:,:,:,v) = (ROIMask_3D_boost2(:,:,:,v) .* mask_myocardium_3D_boost2 .* ~excludeMask_3D > 0 + ...
                    noReflowMask_3D > 0) > 0;
                LocPixCount_boost2(:,:,v) = AHA_16Seg(ROIMask_3D_boost2(:,:,myo_check,v), mask_myocardium_3D_boost2(:,:,myo_check)>0, ...
                    BaseGroove_boost2(myo_check), flip(idx_array_boost2));

                LocPixCount_boost2_v2(:,:,v) = AHA_16Seg(ROIMask_3D_boost2(:,:,myo_check,v), mask_myocardium_3D_boost2(:,:,myo_check)>0, ...
                    BaseGroove_boost2(myo_check), idx_array_boost2);
            end

            auc_roc_array_lge = zeros(length(cutoff_array), 1);
            auc_roc_matrix_lge = zeros(length(cutoff_array), length(thresh_array));
            auc_roc_array_t2 = zeros(length(cutoff_array), 1);
            auc_roc_matrix_t2 = zeros(length(cutoff_array), length(thresh_array));

            auc_roc_array_lge_v2 = zeros(length(cutoff_array), 1);
            auc_roc_matrix_lge_v2 = zeros(length(cutoff_array), length(thresh_array));
            auc_roc_array_t2_v2 = zeros(length(cutoff_array), 1);
            auc_roc_matrix_t2_v2 = zeros(length(cutoff_array), length(thresh_array));

            figure('Visible', 'off');
            for cf = 1:length(cutoff_array)
                cutoff = cutoff_array(cf);

                aha_gt_lge = LocPixCount_lge(:,1) > cutoff;
                aha_gt_t2 = LocPixCount_t2(:,1) > cutoff;
                LocPixCount_boost2_aha = LocPixCount_boost2 > cutoff;
                LocPixCount_boost2_aha_v2 = LocPixCount_boost2_v2 > cutoff;

                % Build ROC points robustly and plot individual points
                %figure();
                %hold on;
                %plot(0,0,'o', 'LineWidth',2);
                %plot(1,1,'o', 'LineWidth',2);

                roc_curve_x_lge = [];
                roc_curve_y_lge = [];

                roc_curve_x_lge(1) = 1;
                roc_curve_y_lge(1) = 1;

                roc_curve_x_t2 = [];
                roc_curve_y_t2 = [];

                roc_curve_x_t2(1) = 1;
                roc_curve_y_t2(1) = 1;

                roc_curve_x_lge_v2 = [];
                roc_curve_y_lge_v2 = [];

                roc_curve_x_lge_v2(1) = 1;
                roc_curve_y_lge_v2(1) = 1;

                roc_curve_x_t2_v2 = [];
                roc_curve_y_t2_v2 = [];

                roc_curve_x_t2_v2(1) = 1;
                roc_curve_y_t2_v2(1) = 1;

                nPos = sum(aha_gt_lge == 1);
                nNeg = sum(aha_gt_lge == 0);
                if nPos == 0
                    warning('No positive examples (aha_gt_lge==1); TPR will be set to 0 for all points.');
                end
                if nNeg == 0
                    warning('No negative examples (aha_gt_lge==0); FPR will be set to 0 for all points.');
                end

                for vv = 1:length(thresh_array)
                    y_pred = LocPixCount_boost2_aha(:,1,vv);
                    % Confusion matrix elements
                    TP = sum((aha_gt_lge == 1) & (y_pred == 1));
                    FP = sum((aha_gt_lge == 0) & (y_pred == 1));
                    TN = sum((aha_gt_lge == 0) & (y_pred == 0));
                    FN = sum((aha_gt_lge == 1) & (y_pred == 0));

                    % Compute rates with guards against zero denominators
                    if (TP + FN) > 0
                        TPR = TP / (TP + FN);   % Sensitivity
                    else
                        TPR = 0; % no positives -> define TPR as 0 for plotting/aggregation
                    end
                    if (FP + TN) > 0
                        FPR = FP / (FP + TN);   % False Positive Rate
                    else
                        FPR = 0; % no negatives -> define FPR as 0
                    end

                    roc_curve_x_lge(end+1) = FPR;
                    roc_curve_y_lge(end+1) = TPR;
                    %fprintf('TPR = %.3f\n', TPR);
                    %fprintf('FPR = %.3f\n', FPR);
                    %plot(FPR, TPR, 'o', 'LineWidth',2);


                    % Another version, AUC of each thresholds
                    if sum(aha_gt_lge) > 0
                        [~,~,~,AUC] = perfcurve(aha_gt_lge,LocPixCount_boost2(:,1,vv),1);
                    else
                        AUC = 0.5;
                    end
                    auc_roc_matrix_lge(cf,vv) = AUC; % Row cutoff, Col threshold Mean+NSD



                    TP = sum((aha_gt_t2 == 1) & (y_pred == 1));
                    FP = sum((aha_gt_t2 == 0) & (y_pred == 1));
                    TN = sum((aha_gt_t2 == 0) & (y_pred == 0));
                    FN = sum((aha_gt_t2 == 1) & (y_pred == 0));

                    % Compute rates with guards against zero denominators
                    if (TP + FN) > 0
                        TPR = TP / (TP + FN);   % Sensitivity
                    else
                        TPR = 0; % no positives -> define TPR as 0 for plotting/aggregation
                    end
                    if (FP + TN) > 0
                        FPR = FP / (FP + TN);   % False Positive Rate
                    else
                        FPR = 0; % no negatives -> define FPR as 0
                    end

                    roc_curve_x_t2(end+1) = FPR;
                    roc_curve_y_t2(end+1) = TPR;

                    % Another version, AUC of each thresholds
                    %[~,~,~,AUC] = perfcurve(aha_gt_t2,LocPixCount_boost2(:,1,vv),1);
                    %auc_roc_matrix_t2(cf,vv) = AUC; % Row cutoff, Col threshold Mean+NSD




                    y_pred = LocPixCount_boost2_aha_v2(:,1,vv);
                    % Confusion matrix elements
                    TP = sum((aha_gt_lge == 1) & (y_pred == 1));
                    FP = sum((aha_gt_lge == 0) & (y_pred == 1));
                    TN = sum((aha_gt_lge == 0) & (y_pred == 0));
                    FN = sum((aha_gt_lge == 1) & (y_pred == 0));

                    % Compute rates with guards against zero denominators
                    if (TP + FN) > 0
                        TPR = TP / (TP + FN);   % Sensitivity
                    else
                        TPR = 0; % no positives -> define TPR as 0 for plotting/aggregation
                    end
                    if (FP + TN) > 0
                        FPR = FP / (FP + TN);   % False Positive Rate
                    else
                        FPR = 0; % no negatives -> define FPR as 0
                    end

                    roc_curve_x_lge_v2(end+1) = FPR;
                    roc_curve_y_lge_v2(end+1) = TPR;

                    % Another version, AUC of each thresholds
                    if sum(aha_gt_lge) > 0
                        [~,~,~,AUC] = perfcurve(aha_gt_lge,LocPixCount_boost2_v2(:,1,vv),1);
                    else
                        AUC = 0.5;
                    end
                    auc_roc_matrix_lge_v2(cf,vv) = AUC; % Row cutoff, Col threshold Mean+NSD
                end

                % Add endpoints (0,0) and (1,1) to ensure full curve coverage
                roc_curve_x_lge(end+1) = 0;
                roc_curve_y_lge(end+1) = 0;
                roc_curve_x_t2(end+1) = 0;
                roc_curve_y_t2(end+1) = 0;
                roc_curve_x_lge_v2(end+1) = 0;
                roc_curve_y_lge_v2(end+1) = 0;

                % Plot the raw points
                % subplot(3,3,cf)
                % plot(roc_curve_x_t2, roc_curve_y_t2, 'o-', 'LineWidth', 1.5);
                % hold on;
                % plot([0 1], [0 1], '--', 'LineWidth', 1.5);
                % axis image;
                % saveas(gcf, fullfile(saveSubDir, 'T2_AUC.png'));

                subplot(3,3,cf)
                plot(roc_curve_x_lge, roc_curve_y_lge, 'o-', 'LineWidth', 1.5);
                hold on;
                plot([0 1], [0 1], '--', 'LineWidth', 1.5);
                axis image;
                saveas(gcf, fullfile(saveSubDir, 'LGE_AUC.png'));

                subplot(3,3,cf)
                plot(roc_curve_x_lge_v2, roc_curve_y_lge_v2, 'o-', 'LineWidth', 1.5);
                hold on;
                plot([0 1], [0 1], '--', 'LineWidth', 1.5);
                axis image;
                saveas(gcf, fullfile(saveSubDir, 'LGE_AUC_v2.png'));

                auc_roc_lge = auc_from_fpr_tpr(roc_curve_x_lge, roc_curve_y_lge);
                fprintf('LGE ROC AUC = %.4f\n', auc_roc_lge);
                auc_roc_array_lge(cf) = auc_roc_lge;

                %auc_roc_t2 = auc_from_fpr_tpr(roc_curve_x_t2, roc_curve_y_t2);
                %fprintf('T2 ROC AUC = %.4f\n', auc_roc_t2);
                %auc_roc_array_t2(cf) = auc_roc_t2;

                auc_roc_lge_v2 = auc_from_fpr_tpr(roc_curve_x_lge_v2, roc_curve_y_lge_v2);
                fprintf('LGE ROC AUC V2 = %.4f\n', auc_roc_lge_v2);
                auc_roc_array_lge_v2(cf) = auc_roc_lge_v2;
            end

            ROC_struct.(target_subjects{nn}).auc_roc_array_t2 = auc_roc_array_t2;
            ROC_struct.(target_subjects{nn}).auc_roc_array_lge = auc_roc_array_lge;
            ROC_struct.(target_subjects{nn}).auc_roc_matrix_lge = auc_roc_matrix_lge;
            ROC_struct.(target_subjects{nn}).auc_roc_matrix_t2 = auc_roc_matrix_t2;
            ROC_struct.(target_subjects{nn}).auc_roc_array_t2_v2 = auc_roc_array_t2_v2;
            ROC_struct.(target_subjects{nn}).auc_roc_array_lge_v2 = auc_roc_array_lge_v2;
            ROC_struct.(target_subjects{nn}).auc_roc_matrix_lge_v2 = auc_roc_matrix_lge_v2;
            ROC_struct.(target_subjects{nn}).auc_roc_matrix_t2_v2 = auc_roc_matrix_t2_v2;


            ROC_struct_persub.auc_roc_array_t2 = auc_roc_array_t2;
            ROC_struct_persub.auc_roc_array_lge = auc_roc_array_lge;
            ROC_struct_persub.auc_roc_matrix_lge = auc_roc_matrix_lge;
            ROC_struct_persub.auc_roc_matrix_t2 = auc_roc_matrix_t2;

            ROC_struct_persub_v2.auc_roc_array_t2 = auc_roc_array_t2_v2;
            ROC_struct_persub_v2.auc_roc_array_lge = auc_roc_array_lge_v2;
            ROC_struct_persub_v2.auc_roc_matrix_lge = auc_roc_matrix_lge_v2;
            ROC_struct_persub_v2.auc_roc_matrix_t2 = auc_roc_matrix_t2_v2;

            save(fullfile(saveSubDir, 'ROC_analysis_persub.mat'), 'ROC_struct_persub');
            save(fullfile(saveSubDir, 'ROC_analysis_persub_v2.mat'), 'ROC_struct_persub_v2');
        end
        % Robust AUC calc handling duplicate FPRs and non-monotonic TPRs
        % xx = roc_curve_x(:);
        % yy = roc_curve_y(:);
        % ok = ~isnan(xx) & ~isnan(yy);
        % xx = xx(ok); yy = yy(ok);
        % 
        % % Sort by FPR
        % [xxs, ord] = sort(xx);
        % yys = yy(ord);
        % 
        % % Collapse duplicate FPR by taking maximum TPR (upper envelope)
        % [xu, ~, ic] = unique(xxs);
        % yu = accumarray(ic, yys, [], @max);
        % 
        % % Enforce non-decreasing TPR (monotonic ROC curve)
        % yu = cummax(yu);
        % 
        % % Ensure endpoints at 0 and 1 exist
        % if xu(1) > 0
        %     xu = [0; xu];
        %     yu = [0; yu];
        % end
        % if xu(end) < 1
        %     xu = [xu; 1];
        %     yu = [yu; 1];
        % end
        % 
        % % Final sort just in case
        % [xu, sidx] = sort(xu);
        % yu = yu(sidx);
        % 
        % % Compute AUC via trapezoidal integration
        % auc_roc = trapz(xu, yu);
        % fprintf('ROC AUC = %.6f\n', auc_roc);
        

        % % Compute AUC for ROC curve
        % [xs, idx] = sort(roc_curve_x);
        % ys = roc_curve_y(idx);
        % % Collapse duplicate x by taking max y for each unique x (ensures valid integration)
        % [xu, ~, ic] = unique(xs);
        % yu = accumarray(ic, ys, [], @max);
        % auc_roc = trapz(xu, yu);
        % fprintf('ROC AUC = %.4f\n', auc_roc);
        


        % figure;
        % numSlices = size(vol_img_3D_lge, 3);
        % for s = 1:numSlices
        %     subplot(ceil(sqrt(numSlices)), ceil(sqrt(numSlices)), s);
        %     imagesc(vol_img_3D_lge(:,:,s) .* mask_myocardium_3D_lge(:,:,s));
        %     axis image off;
        %     title(sprintf('Slice %d', s));
        %     colormap gray;
        % end
        % sgtitle('LGE Volume Slices');

        % figure;
        % numSlices = size(vol_img_3D_boost2, 3);
        % for s = 1:numSlices
        %     subplot(ceil(sqrt(numSlices)), ceil(sqrt(numSlices)), s);
        %     imagesc(vol_img_3D_boost2(:,:,s) .* mask_myocardium_3D_boost2(:,:,s));
        %     %imagesc(vol_img_3D_boost2(:,:,s));
        %     axis image off;
        %     title(sprintf('Slice %d', s));
        %     colormap gray;
        % end
        % sgtitle('BOOST2 Volume Slices');

        % figure;
        % numSlices = size(vol_img_3D_lge, 3);
        % for s = 1:numSlices
        %     subplot(ceil(sqrt(numSlices)), ceil(sqrt(numSlices)), s);
        %     imagesc(vol_img_3D_lge(:,:,s) .* freeROIMask_3D_lge(:,:,s));
        %     axis image off;
        %     title(sprintf('Slice %d', s));
        %     colormap gray;
        % end
        % sgtitle('LGE Volume Slices Masked by freeROI');



    else
        fprintf('Skipping subject: %s due to missing data.\n', target_subjects{i});
    end
end

save(fullfile(saveDir, 'ROC_analysis.mat'), 'ROC_struct');

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

function [x_array_lge, y_array_lge] = check_coords(vol_img_3D_lge, x_array_lge, y_array_lge)

        if any(x_array_lge < size(vol_img_3D_lge,1)/8)
            x_array_lge(x_array_lge < size(vol_img_3D_lge,1)/8) = x_array_lge(x_array_lge < size(vol_img_3D_lge,1)/8) * 8;
        end

        if any(y_array_lge < size(vol_img_3D_lge,2)/8)
            y_array_lge(y_array_lge < size(vol_img_3D_lge,2)/8) = y_array_lge(y_array_lge < size(vol_img_3D_lge,2)/8) * 8;
        end

        if any(x_array_lge > size(vol_img_3D_lge,1)*7/8)
            x_array_lge(x_array_lge > size(vol_img_3D_lge,1)*7/8) = x_array_lge(x_array_lge > size(vol_img_3D_lge,1)*7/8) / 8;
        end

        if any(y_array_lge > size(vol_img_3D_lge,2)*7/8)
            y_array_lge(y_array_lge > size(vol_img_3D_lge,2)*7/8) = y_array_lge(y_array_lge > size(vol_img_3D_lge,2)*7/8) / 8;
        end
end

function auc = auc_from_fpr_tpr(fpr, tpr)
% AUC under ROC curve with duplicate-FPR handling.
% - For repeated FPR values, keeps the maximum TPR
% - Enforces nondecreasing TPR vs FPR
% - Adds (0,0) and (1,1) if missing
%
% Usage:
%   auc = auc_from_fpr_tpr(fpr, tpr)

    % column-ize and clean
    fpr = fpr(:); tpr = tpr(:);
    m = ~isnan(fpr) & ~isnan(tpr);
    fpr = fpr(m); tpr = tpr(m);

    % clip to [0,1] just in case
    fpr = max(0, min(1, fpr));
    tpr = max(0, min(1, tpr));

    % ensure endpoints present
    if isempty(fpr) || fpr(1) ~= 0 || tpr(1) ~= 0
        fpr = [fpr; 0]; tpr = [tpr; 0];
    end
    if ~any(fpr == 1)
        fpr = [fpr; 1]; tpr = [tpr; 1];
    end

    % sort by FPR, then by TPR
    [~, ord] = sortrows([fpr tpr], [1 2]);
    fpr = fpr(ord); tpr = tpr(ord);

    % group identical FPRs -> take max TPR for each FPR
    % [ufpr, ~, g] = unique(tpr, 'stable');
    % unique on (fpr,tpr) pairs (treat coordinates as complete)
    [pairs, ~, g] = unique([fpr, tpr], 'rows', 'stable');
    ufpr = pairs(:,1);
    tpr_mono = pairs(:,2);
    % For each unique (fpr,tpr) pair take the maximum tpr (identical pairs -> same tpr)
    %tpr_max = accumarray(g, tpr, [], @max);

    % enforce nondecreasing ROC (optional but common):
    % if numerical noise drops a point, lift it up
    %tpr_mono = cummax(tpr_max);

    % trapezoidal integration
    auc = trapz(ufpr, tpr_mono);
end