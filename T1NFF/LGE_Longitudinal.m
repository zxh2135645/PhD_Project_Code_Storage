clear all;
close all;

addpath('../function/');
base_dir = uigetdir;
contour_glob = glob(cat(2, base_dir, '/ContourData/*'));
%Names = ExtractNames(contour_glob);
Names = {'Merry', 'Ryn', 'Mojave', 'Sahara', 'ZZ', 'Tina', 'Sunny', 'Queenie', 'Hope', 'Gobi', 'Felicity', 'Evelyn', '18D15', '18D16', '11D05', '11D26', '11D33'};
%time_points = {'0D_baseline','1D', '7D', '28D', '8WK', '6MO', '9MO', '1YR', '15YR'};
time_points = {'7D', '8WK', '12WK', '14WK' '4MO', '6MO', '9MO', '1YR', '15YR'};

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

sequence_label = {'T1', 'T2star', 'LGE'};
anatomy_label = {'BloodPool', 'excludeArea', 'freeROI', 'Heart', 'Myocardium', 'MyoReference', 'noReflowArea'}; 
%name_check = 'Evelyn';
%starting_point = find(strcmp(name_check, Names),1);

save_dir = cat(2, base_dir, '/img/');
data_save_dir = cat(2, base_dir, '/data/');

label_t1 = sequence_label{1};
label_lge = sequence_label{3};
label_t2star = sequence_label{2};

metrics_save_dir = cat(2, base_dir, '/Results/');
if ~exist(metrics_save_dir, 'dir')
   mkdir(metrics_save_dir); 
end

% Before analysis, parse pre_QualControl
load(cat(2, metrics_save_dir, 'pre_QualControl.mat'));

status_check = struct;
for n = 1:length(Names)
    name = Names{n};
    status_check(n).Name = name;
    status_check(n).status = [];
    status_check(n).status_final = [];
    tp_count = 1;
    for tp = 1:length(time_points)
        
        time_point = time_points{end-tp+1};
        tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', name, '_', time_point,  '/');
        
        if ~exist(tp_dir, 'dir')
            disp(cat(2, 'No folder at: ', name, ' ', time_point));
        else
            for i = 1:(length(fieldnames(pre_QualControl(n).status))-1)
                slc_loc = cat(2, 'Slice', num2str(i));
                if pre_QualControl(n).status(end-tp+1).(slc_loc) == 1
                    status_check(n).status(tp_count, i) = 1;
                elseif pre_QualControl(n).status(end-tp+1).(slc_loc) == 0
                    status_check(n).status(tp_count, i) = 0;
                end
            end
            
            tp_count = tp_count + 1;
        end
    end
    for i = 1:(length(fieldnames(pre_QualControl(n).status))-1)
        if i <= size(status_check(n).status,2)
            status_check(n).status_final(1,i) = all(status_check(n).status(:,i));
            % The timepoints of status_check goes from end to beginning
        end
    end
end

%% Pull LGE
% 09/03/2023 modified for LGE
se = strel('disk', 1);
metrics = struct;
L_cell = {};
count = 1;
for n = 1:length(Names)
%for n = 1:1
    % for n = starting_point:starting_point
    % Do not need to pull up images for baseline
    name = Names{n};
    name_save_dir = cat(2, save_dir, name);
    if ~exist(name_save_dir, 'dir')
        mkdir(name_save_dir);
    end
    name_data_save_dir = cat(2, data_save_dir, name);
    if ~exist(name_data_save_dir, 'dir')
        mkdir(name_data_save_dir);
    end
    tp_count = 0;

    
    metrics(n).name = name;
    metrics(n).TimePoints = struct;
    for tp = 1:length(time_points)
    %for tp = 2:2
        time_point = time_points{end-tp+1};
        tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', name, '_', time_point,  '/');
        if ~exist(tp_dir, 'dir')
            disp(cat(2, 'No folder at: ', name, ' ', time_point));
        else
            
            tp_dir2 = cat(2, name_save_dir, '/', name, '_', time_point, '/');
            if ~exist(tp_dir2, 'dir')
                mkdir(tp_dir2);
            end
            %LR_mdl_fname = cat(2, name_data_save_dir, '/LinearRegression_', name, '_', time_point, '.mat');

            tp_count = tp_count+1;
            status = status_check(n).status(tp_count,:);
            % AHA Segment

            if (strcmp(name, '18D16') && strcmp(time_point, '9MO'))
                % WHY 
                % Some issue with status?
                disp(cat(2, 'Skipped: ', name, ' ', time_point))
            else
                % T1 slice location
                load(cat(2, tp_dir, label_t1, '/', label_t1, '_SliceLoc.mat'));
                [slc_array_t1, idx_reordered] = sort(slc_array);

                % LGE

                myo_glob = glob(cat(2, tp_dir, label_lge, '/', anatomy_label{5}, '/*'));
                roi_glob = glob(cat(2, tp_dir, label_lge, '/',anatomy_label{3}, '/*'));
                remote_glob = glob(cat(2, tp_dir, label_lge, '/',anatomy_label{6}, '/*'));

                load(cat(2, tp_dir, label_lge, '/', label_lge, '_vol_img_3D.mat'));
                load(myo_glob{1});
                load(roi_glob{1});
                load(remote_glob{1});
                load(cat(2, tp_dir, label_lge, '/', label_lge, '_SliceLoc.mat'));
                slc_array_lge = slc_array;

                clear myo_lge_eroded
                for i = 1:size(mask_myocardium_3D, 3)
                    myo_lge_eroded(:,:,i) = imerode(mask_myocardium_3D(:,:,i), se);
                end

                roi_in_myo_lge = myo_lge_eroded .* freeROIMask_3D;
                remote_in_myo_lge = myo_lge_eroded.* myoRefMask_3D;
                roi_lge = roi_in_myo_lge .* vol_img_3D;
                remote_lge = remote_in_myo_lge .* vol_img_3D;
                lge = vol_img_3D;
                myo_lge = myo_lge_eroded;

                idx_reordered = Func_AlignSliceLoc(slc_array_t1, slc_array_lge);
                lge = lge(:,:,idx_reordered);
                myo_lge = myo_lge(:,:,idx_reordered);
                remote_lge = remote_lge(:,:,idx_reordered);
                roi_lge = roi_lge(:,:,idx_reordered);
                remote_in_myo_lge = remote_in_myo_lge(:,:,idx_reordered);
                roi_in_myo_lge = roi_in_myo_lge(:,:,idx_reordered);

                if (strcmp(name, 'Tina') && strcmp(time_point, '6MO'))

                    % Do nothing
                else
                    exclude_idx = find(status == 0);
                    lge(:,:,exclude_idx) = [];
                    roi_in_myo_lge(:,:,exclude_idx) = [];
                    remote_in_myo_lge(:,:,exclude_idx) = [];
                    myo_lge(:,:,exclude_idx) = [];

                end
                metrics(n).TimePoints(tp_count).time_point = time_point;
                metrics(n).TimePoints(tp_count).mean_lge_roi = mean(nonzeros(lge .* roi_in_myo_lge));
                metrics(n).TimePoints(tp_count).mean_lge_remote = mean(nonzeros(lge .* remote_in_myo_lge));
                metrics(n).TimePoints(tp_count).std_lge_roi = std(nonzeros(lge .* roi_in_myo_lge));
                metrics(n).TimePoints(tp_count).std_lge_remote = std(nonzeros(lge .* remote_in_myo_lge));
                metrics(n).TimePoints(tp_count).skewness_lge_roi = skewness(nonzeros(lge .* roi_in_myo_lge));
                metrics(n).TimePoints(tp_count).skewness_lge_remote = skewness(nonzeros(lge .* remote_in_myo_lge));
                metrics(n).TimePoints(tp_count).kurtosis_lge_roi = kurtosis(nonzeros(lge .* roi_in_myo_lge))-3;
                metrics(n).TimePoints(tp_count).kurtosis_lge_remote = kurtosis(nonzeros(lge .* remote_in_myo_lge))-3;

                metrics(n).TimePoints(tp_count).SliceAnalysis = struct;
                metrics(n).TimePoints(tp_count).HeteroAnalysis = struct;

                % lge_norm = (lge - min(lge(:)))./ max(lge(:));
                temp_roi = (lge .* roi_in_myo_lge);
                temp_remote = (lge .* remote_in_myo_lge);
                temp_myo = (lge .* myo_lge);
                lge_norm_remote = zeros(size(lge));
                lge_norm_roi = zeros(size(lge));


                for slc = 1:size(roi_in_myo_lge,3)
                    %temp_norm_roi = temp_roi(:,:,slc);
                    %temp_norm_roi(temp_norm_roi == 0) = nan;
                    % temp_norm_roi = uint8((temp_roi(:,:,slc) - min(min((nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc))))))./ (max(max(nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc))))-min(min((nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc)))))) * 256);
                    % temp_norm_remote = uint8((temp_remote(:,:,slc) - min(min((nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc))))))./ (max(max(nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc))))-min(min((nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc)))))) * 256);

                    temp_norm_roi = (temp_roi(:,:,slc) - min(min((nonzeros(temp_myo(:,:,slc))))))./ (max(max(nonzeros(temp_myo(:,:,slc))))-min(min((nonzeros(temp_myo(:,:,slc))))));
                    temp_norm_remote = (temp_remote(:,:,slc) - min(min((nonzeros(temp_myo(:,:,slc))))))./ (max(max(nonzeros(temp_myo(:,:,slc))))-min(min((nonzeros(temp_myo(:,:,slc))))));

                    % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_lge_roi = mean(nonzeros(lge_norm(:,:,slc) .* roi_in_myo_lge(:,:,slc)));
                    % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_lge_remote = mean(nonzeros(lge_norm(:,:,slc) .* remote_in_myo_lge(:,:,slc)));
                    % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_lge_roi = std(nonzeros(lge_norm(:,:,slc) .* roi_in_myo_lge(:,:,slc)));
                    % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_lge_remote = std(nonzeros(lge_norm(:,:,slc) .* remote_in_myo_lge(:,:,slc)));
                    % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_lge_roi = skewness(nonzeros(lge_norm(:,:,slc) .* roi_in_myo_lge(:,:,slc)));
                    % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_lge_remote = skewness(nonzeros(lge_norm(:,:,slc) .* remote_in_myo_lge(:,:,slc)));
                    % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_lge_roi = kurtosis(nonzeros(lge_norm(:,:,slc) .* roi_in_myo_lge(:,:,slc)))-3;
                    % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_lge_remote = kurtosis(nonzeros(lge_norm(:,:,slc) .* remote_in_myo_lge(:,:,slc)))-3;
                    %
                    % metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).hetero_roi = Func_Hetero_Analysis(slc, roi_in_myo_lge, lge_norm);
                    % metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).hetero_remote = Func_Hetero_Analysis(slc, remote_in_myo_lge, lge_norm);

                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_lge_roi = mean(nonzeros(temp_norm_roi .* roi_in_myo_lge(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_lge_remote = mean(nonzeros(temp_norm_remote .* remote_in_myo_lge(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_lge_roi = std(nonzeros(temp_norm_roi .* roi_in_myo_lge(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_lge_remote = std(nonzeros(temp_norm_remote .* remote_in_myo_lge(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_lge_roi = skewness(nonzeros(temp_norm_roi .* roi_in_myo_lge(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_lge_remote = skewness(nonzeros(temp_norm_remote .* remote_in_myo_lge(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_lge_roi = kurtosis(nonzeros(temp_norm_roi .* roi_in_myo_lge(:,:,slc)))-3;
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_lge_remote = kurtosis(nonzeros(temp_norm_remote .* remote_in_myo_lge(:,:,slc)))-3;

                    lge_norm_roi(:,:,slc) = temp_norm_roi;
                    lge_norm_remote(:,:,slc) = temp_norm_remote;
                    metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).hetero_roi = Func_Hetero_Analysis(slc, roi_in_myo_lge, lge_norm_roi);
                    metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).hetero_remote = Func_Hetero_Analysis(slc, remote_in_myo_lge, lge_norm_remote);

                    temp_norm_roi = uint8((temp_roi(:,:,slc) - min(min((nonzeros(temp_myo(:,:,slc))))))./ (max(max(nonzeros(temp_myo(:,:,slc))))-min(min((nonzeros(temp_myo(:,:,slc)))))) * 256);
                    temp_norm_remote = uint8((temp_remote(:,:,slc) - min(min((nonzeros(temp_myo(:,:,slc))))))./ (max(max(nonzeros(temp_myo(:,:,slc))))-min(min((nonzeros(temp_myo(:,:,slc)))))) * 256);


                    bipolar = temp_roi(:,:,slc) - mean(nonzeros(temp_remote(:,:,slc)));
                    weighted_map = double(bipolar < 0) .* 2 .* roi_in_myo_lge(:,:,slc) .* bipolar + double(bipolar >= 0) .* roi_in_myo_lge(:,:,slc) .* bipolar;

                    bipolar_remote = temp_remote(:,:,slc) - mean(nonzeros(temp_remote(:,:,slc)));
                    weighted_map_remote = double(bipolar_remote < 0) .* 2 .* remote_in_myo_lge(:,:,slc) .* bipolar_remote + double(bipolar_remote >= 0) .* remote_in_myo_lge(:,:,slc) .* bipolar_remote;

                    lb = min(nonzeros(weighted_map(:)));
                    ub = max(nonzeros(weighted_map(:)));
                    temp_norm_roi = uint8((weighted_map - lb)./ (ub-lb) * 256);

                    weighted_map_remote(weighted_map_remote < lb) = lb;
                    weighted_map_remote(weighted_map_remote > ub) = ub;
                    temp_norm_remote = uint8((weighted_map_remote - lb)./ (ub-lb) * 256);

                    non_roi_value = uint8((0 - lb) ./ (ub - lb) * 256);

                    lge_norm_roi_slc = temp_norm_roi;
                    lge_norm_roi_slc(lge_norm_roi_slc == non_roi_value) = [];
                    %figure(); imhist(lge_norm_roi_slc,20);
                    p_roi = imhist(lge_norm_roi_slc,20);
                    nonZeros = find(p_roi);
                    len = length((nonZeros));
                    pNonZeros = zeros(1,len);

                    for i = 1:len
                        pNonZeros(i) = p_roi(nonZeros(i));
                    end

                    % normalize pNonZeros so that sum(p) is one.
                    pNonZeros = pNonZeros ./ sum(p_roi);
                    E_roi = -sum(pNonZeros.*log2(pNonZeros));
                    

                    non_remote_value = uint8((0 - lb) ./ (ub - lb) * 256);

                    lge_norm_remote_slc = temp_norm_remote;
                    lge_norm_remote_slc(lge_norm_remote_slc == non_remote_value) = [];
                    %figure(); imhist(lge_norm_remote_slc,20);
                    p_remote = imhist(lge_norm_remote_slc,20);

                    nonZeros = find(p_remote);
                    len = length((nonZeros));
                    pNonZeros = zeros(1,len);

                    for i = 1:len
                        pNonZeros(i) = p_remote(nonZeros(i));
                    end

                    % normalize pNonZeros so that sum(p) is one.
                    pNonZeros = pNonZeros ./ sum(p_remote);
                    E_remote = -sum(pNonZeros.*log2(pNonZeros));


                    metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).entropy_roi = E_roi;
                    metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).entropy_remote = E_remote;

                    L = bwlabel(roi_in_myo_lge(:,:,slc));
                    
                    if length(unique(L)) > 2
                        L_cell{count} = cat(2, name, ' ', time_point, ' Slice ', num2str(slc));
                        count = count + 1;
                    end
                end
                


%                 figure();
%                 histogram(nonzeros(t1 .* remote_in_myo_t1), 'Normalization', 'Probability');
%                 hold on;
%                 histogram(nonzeros(t1 .* roi_in_myo_t1), 'Normalization', 'Probability');

            end
        end
    end
    close all;
end

%% Save data
fname = 'Demographic_Metrics_rim_LGE_weightedEntropy_20bins';
save(cat(2, data_save_dir, '/', fname), 'metrics');
