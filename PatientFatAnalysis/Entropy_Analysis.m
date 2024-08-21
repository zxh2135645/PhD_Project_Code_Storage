clear all;
close all;
clc;
current_dir = pwd;
% Patient data configuration for Khalid
%% 
addpath('../function/');
addpath('../AHA16Segment/');
addpath('../T1NFF/');
base_dir = uigetdir;

contour_glob = glob(cat(2, base_dir, '/ContourData/*'));
Names = cell(length(contour_glob), 1); 
for i = 1:length(contour_glob)
    strings = strsplit(contour_glob{i},'/');
    Names{i} = strings{end-1};
end

sequence_label = {'LGE', 'T1'};
anatomy_label = {'BloodPool', 'excludeArea', 'freeROI', 'Heart', 'Myocardium', 'MyoReference', 'noReflowArea'}; 

names_to_rule_out = {};
RuleOutLabel = NameRuleOutFunc(Names, names_to_rule_out);
Names = Names(RuleOutLabel == 0);

name_check = {'484060000009'};
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


output_label = {'LGE', 'T2star'};
save_dir = GetFullPath(cat(2, base_dir, '/Analysis/'));
data_save_dir = cat(2, base_dir, '/data/');

time_points = {'FU'};
time_points = {'BL'};
%time_points = {'BL2'};

label_lge = sequence_label{1};
label_t1 = sequence_label{2};
label_t2 = 'T2';
label_mag = 'MAG';
label_psir = 'PSIR';
label_t1molli = 'T1MOLLI';
label_t2star = 'T2star';

metrics_save_dir = cat(2, base_dir, '/Results/');
if ~exist(metrics_save_dir, 'dir')
   mkdir(metrics_save_dir); 
end

% Read excel file
% T = readtable(cat(2, base_dir, '/STEMI_with_IMH-1.xlsx'), 'VariableNamingRule', 'preserve');
T = readtable(cat(2, base_dir, '/STEMI_with_IMH-1.xlsx'));
id_array = T.AnonymizationID;
hemo_array = T.withOrWithout;


%% Main Body of Entropy Analysis
se = strel('disk', 1);
% metrics = struct;
% metrics_t1 = struct;
L_cell = {};
count = 1;

min_roi_value = 450;
max_roi_value = 1950;

for n = starting_point:starting_point
%for n = starting_point:length(Names)
%for n = 29:29
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

    %metrics_t1(n).name = name;
    %metrics_t1(n).TimePoints = struct;

    for tp = 1:length(time_points)
    %for tp = 1:1
        time_point = time_points{end-tp+1};
        tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', time_point,  '/');
        if ~exist(cat(2, tp_dir, label_t1, '/'), 'dir')
            disp(cat(2, 'No folder at: ', name, ' ', time_point));
        else

            tp_dir2 = cat(2, name_save_dir, '/', name, '/', time_point, '/');
            if ~exist(tp_dir2, 'dir')
                mkdir(tp_dir2);
            end
            %LR_mdl_fname = cat(2, name_data_save_dir, '/LinearRegression_', name, '_', time_point, '.mat');

            tp_count = tp_count+1;
            % status = status_check(n).status(tp_count,:);
            % AHA Segment


            % T1 slice location
            load(cat(2, tp_dir, label_t1, '/', label_t1, '_SliceLoc.mat'));
            % [slc_array_t1, idx_reordered] = sort(slc_array);
            slc_array_t1 = slc_array;

            % LGE
            % myo_glob = glob(cat(2, tp_dir, label_lge, '/', anatomy_label{5}, '/*'));
            % roi_glob = glob(cat(2, tp_dir, label_lge, '/',anatomy_label{3}, '/*'));
            % remote_glob = glob(cat(2, tp_dir, label_lge, '/',anatomy_label{6}, '/*'));
            % 
            % load(cat(2, tp_dir, label_lge, '/', label_mag, '_vol_img_3D.mat'));
            % load(myo_glob{1});
            % load(roi_glob{1});
            % load(remote_glob{1});
            % load(cat(2, tp_dir, label_lge, '/', label_lge, '_SliceLoc.mat'));
            % slc_array_lge = slc_array;
            % 
            % clear myo_lge_eroded
            % for i = 1:size(mask_myocardium_3D, 3)
            %     myo_lge_eroded(:,:,i) = imerode(mask_myocardium_3D(:,:,i), se);
            % end
            % 
            % vol_img_3D = svi3.sig_vol_img_3D;
            % roi_in_myo_lge = myo_lge_eroded .* freeROIMask_3D;
            % remote_in_myo_lge = myo_lge_eroded.* myoRefMask_3D;
            % roi_lge = roi_in_myo_lge .* vol_img_3D;
            % remote_lge = remote_in_myo_lge .* vol_img_3D;
            % lge = vol_img_3D;
            % myo_lge = myo_lge_eroded;
            % 
            % idx_reordered = Func_AlignSliceLoc(slc_array_t1, slc_array_lge);
            % lge = lge(:,:,idx_reordered);
            % myo_lge = myo_lge(:,:,idx_reordered);
            % remote_lge = remote_lge(:,:,idx_reordered);
            % roi_lge = roi_lge(:,:,idx_reordered);
            % remote_in_myo_lge = remote_in_myo_lge(:,:,idx_reordered);
            % roi_in_myo_lge = roi_in_myo_lge(:,:,idx_reordered);

            %exclude_idx = find(status == 0);
            %lge(:,:,exclude_idx) = [];
            %roi_in_myo_lge(:,:,exclude_idx) = [];
            %remote_in_myo_lge(:,:,exclude_idx) = [];
            %myo_lge(:,:,exclude_idx) = [];

            % metrics(n).TimePoints(tp_count).time_point = time_point;
            % metrics(n).TimePoints(tp_count).mean_lge_roi = mean(nonzeros(lge .* roi_in_myo_lge));
            % metrics(n).TimePoints(tp_count).mean_lge_remote = mean(nonzeros(lge .* remote_in_myo_lge));
            % metrics(n).TimePoints(tp_count).std_lge_roi = std(nonzeros(lge .* roi_in_myo_lge));
            % metrics(n).TimePoints(tp_count).std_lge_remote = std(nonzeros(lge .* remote_in_myo_lge));
            % metrics(n).TimePoints(tp_count).skewness_lge_roi = skewness(nonzeros(lge .* roi_in_myo_lge));
            % metrics(n).TimePoints(tp_count).skewness_lge_remote = skewness(nonzeros(lge .* remote_in_myo_lge));
            % metrics(n).TimePoints(tp_count).kurtosis_lge_roi = kurtosis(nonzeros(lge .* roi_in_myo_lge))-3;
            % metrics(n).TimePoints(tp_count).kurtosis_lge_remote = kurtosis(nonzeros(lge .* remote_in_myo_lge))-3;


            % metrics(n).TimePoints(tp_count).SliceAnalysis = struct;
            % metrics(n).TimePoints(tp_count).HeteroAnalysis = struct;
            % 
            % % lge_norm = (lge - min(lge(:)))./ max(lge(:));
            % temp_roi = (lge .* roi_in_myo_lge);
            % temp_remote = (lge .* remote_in_myo_lge);
            % temp_myo = (lge .* myo_lge);
            % lge_norm_remote = zeros(size(lge));
            % lge_norm_roi = zeros(size(lge));
            % 
            % 
            % for slc = 1:size(roi_in_myo_lge,3)
            %     %temp_norm_roi = temp_roi(:,:,slc);
            %     %temp_norm_roi(temp_norm_roi == 0) = nan;
            %     % temp_norm_roi = uint8((temp_roi(:,:,slc) - min(min((nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc))))))./ (max(max(nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc))))-min(min((nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc)))))) * 256);
            %     % temp_norm_remote = uint8((temp_remote(:,:,slc) - min(min((nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc))))))./ (max(max(nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc))))-min(min((nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc)))))) * 256);
            % 
            %     temp_norm_roi = (temp_roi(:,:,slc) - min(min((nonzeros(temp_myo(:,:,slc))))))./ (max(max(nonzeros(temp_myo(:,:,slc))))-min(min((nonzeros(temp_myo(:,:,slc))))));
            %     temp_norm_remote = (temp_remote(:,:,slc) - min(min((nonzeros(temp_myo(:,:,slc))))))./ (max(max(nonzeros(temp_myo(:,:,slc))))-min(min((nonzeros(temp_myo(:,:,slc))))));
            % 
            %     % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_lge_roi = mean(nonzeros(lge_norm(:,:,slc) .* roi_in_myo_lge(:,:,slc)));
            %     % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_lge_remote = mean(nonzeros(lge_norm(:,:,slc) .* remote_in_myo_lge(:,:,slc)));
            %     % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_lge_roi = std(nonzeros(lge_norm(:,:,slc) .* roi_in_myo_lge(:,:,slc)));
            %     % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_lge_remote = std(nonzeros(lge_norm(:,:,slc) .* remote_in_myo_lge(:,:,slc)));
            %     % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_lge_roi = skewness(nonzeros(lge_norm(:,:,slc) .* roi_in_myo_lge(:,:,slc)));
            %     % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_lge_remote = skewness(nonzeros(lge_norm(:,:,slc) .* remote_in_myo_lge(:,:,slc)));
            %     % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_lge_roi = kurtosis(nonzeros(lge_norm(:,:,slc) .* roi_in_myo_lge(:,:,slc)))-3;
            %     % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_lge_remote = kurtosis(nonzeros(lge_norm(:,:,slc) .* remote_in_myo_lge(:,:,slc)))-3;
            %     %
            %     % metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).hetero_roi = Func_Hetero_Analysis(slc, roi_in_myo_lge, lge_norm);
            %     % metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).hetero_remote = Func_Hetero_Analysis(slc, remote_in_myo_lge, lge_norm);
            % 
            %     metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_lge_roi = mean(nonzeros(temp_norm_roi .* roi_in_myo_lge(:,:,slc)));
            %     metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_lge_remote = mean(nonzeros(temp_norm_remote .* remote_in_myo_lge(:,:,slc)));
            %     metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_lge_roi = std(nonzeros(temp_norm_roi .* roi_in_myo_lge(:,:,slc)));
            %     metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_lge_remote = std(nonzeros(temp_norm_remote .* remote_in_myo_lge(:,:,slc)));
            %     metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_lge_roi = skewness(nonzeros(temp_norm_roi .* roi_in_myo_lge(:,:,slc)));
            %     metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_lge_remote = skewness(nonzeros(temp_norm_remote .* remote_in_myo_lge(:,:,slc)));
            %     metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_lge_roi = kurtosis(nonzeros(temp_norm_roi .* roi_in_myo_lge(:,:,slc)))-3;
            %     metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_lge_remote = kurtosis(nonzeros(temp_norm_remote .* remote_in_myo_lge(:,:,slc)))-3;
            % 
            %     lge_norm_roi(:,:,slc) = temp_norm_roi;
            %     lge_norm_remote(:,:,slc) = temp_norm_remote;
            %     %metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).hetero_roi = Func_Hetero_Analysis(slc, roi_in_myo_lge, lge_norm_roi);
            %     %metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).hetero_remote = Func_Hetero_Analysis(slc, remote_in_myo_lge, lge_norm_remote);
            % 
            %     temp_norm_roi = uint8((temp_roi(:,:,slc) - min(min((nonzeros(temp_myo(:,:,slc))))))./ (max(max(nonzeros(temp_myo(:,:,slc))))-min(min((nonzeros(temp_myo(:,:,slc)))))) * 256);
            %     temp_norm_remote = uint8((temp_remote(:,:,slc) - min(min((nonzeros(temp_myo(:,:,slc))))))./ (max(max(nonzeros(temp_myo(:,:,slc))))-min(min((nonzeros(temp_myo(:,:,slc)))))) * 256);
            % 
            %     bipolar = temp_roi(:,:,slc) - mean(nonzeros(temp_remote(:,:,slc)));
            %     weighted_map = double(bipolar < 0) .* 2 .* roi_in_myo_lge(:,:,slc) .* bipolar + double(bipolar >= 0) .* roi_in_myo_lge(:,:,slc) .* bipolar;
            % 
            %     bipolar_remote = temp_remote(:,:,slc) - mean(nonzeros(temp_remote(:,:,slc)));
            %     weighted_map_remote = double(bipolar_remote < 0) .* 2 .* remote_in_myo_lge(:,:,slc) .* bipolar_remote + double(bipolar_remote >= 0) .* remote_in_myo_lge(:,:,slc) .* bipolar_remote;
            % 
            %     lb = min(nonzeros(weighted_map(:)));
            %     ub = max(nonzeros(weighted_map(:)));
            %     temp_norm_roi = uint8((weighted_map - lb)./ (ub-lb) * 256);
            % 
            %     weighted_map_remote(weighted_map_remote < lb) = lb;
            %     weighted_map_remote(weighted_map_remote > ub) = ub;
            %     temp_norm_remote = uint8((weighted_map_remote - lb)./ (ub-lb) * 256);
            % 
            %     non_roi_value = uint8((0 - lb) ./ (ub - lb) * 256);
            % 
            % 
            %     bins = 20;
            %     lge_norm_roi_slc = temp_norm_roi;
            %     lge_norm_roi_slc(lge_norm_roi_slc == non_roi_value) = [];
            %     %figure(); imhist(lge_norm_roi_slc,20);
            %     p_roi = imhist(lge_norm_roi_slc, bins);
            %     nonZeros = find(p_roi);
            %     len = length((nonZeros));
            %     pNonZeros = zeros(1,len);
            % 
            %     for i = 1:len
            %         pNonZeros(i) = p_roi(nonZeros(i));
            %     end
            % 
            %     % normalize pNonZeros so that sum(p) is one.
            %     pNonZeros = pNonZeros ./ sum(p_roi);
            %     E_roi = -sum(pNonZeros.*log2(pNonZeros));
            % 
            %     non_remote_value = uint8((0 - lb) ./ (ub - lb) * 256);
            % 
            %     lge_norm_remote_slc = temp_norm_remote;
            %     lge_norm_remote_slc(lge_norm_remote_slc == non_remote_value ) = [];
            %     %figure(); imhist(lge_norm_remote_slc,20);
            %     p_remote = imhist(lge_norm_remote_slc,bins);
            % 
            %     nonZeros = find(p_remote);
            %     len = length((nonZeros));
            %     pNonZeros = zeros(1,len);
            % 
            %     for i = 1:len
            %         pNonZeros(i) = p_remote(nonZeros(i));
            %     end
            % 
            %     % normalize pNonZeros so that sum(p) is one.
            %     pNonZeros = pNonZeros ./ sum(p_remote);
            %     E_remote = -sum(pNonZeros.*log2(pNonZeros));
            % 
            % 
            %     metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).entropy_roi = E_roi;
            %     metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).entropy_remote = E_remote;
            % 
            %     L = bwlabel(roi_in_myo_lge(:,:,slc));
            % 
            %     if length(unique(L)) > 2
            %         L_cell{count} = cat(2, name, ' ', time_point, ' Slice ', num2str(slc));
            %         count = count + 1;
            %     end
            % 
            % end

                %                 figure();
                %                 histogram(nonzeros(t1 .* remote_in_myo_t1), 'Normalization', 'Probability');
                %                 hold on;
                %                 histogram(nonzeros(t1 .* roi_in_myo_t1), 'Normalization', 'Probability');

                % T1
                myo_glob = glob(cat(2, tp_dir, label_t1, '/', anatomy_label{5}, '/*'));
                roi_glob = glob(cat(2, tp_dir, label_t1, '/',anatomy_label{3}, '/*'));
                remote_glob = glob(cat(2, tp_dir, label_t1, '/',anatomy_label{6}, '/*'));

                load(cat(2, tp_dir, label_t1, '/', label_t1molli, '_vol_img_3D.mat'));
                load(myo_glob{1});
                load(roi_glob{1});
                load(remote_glob{1});
                load(cat(2, tp_dir, label_t1, '/', label_t1, '_SliceLoc.mat'));
                slc_array_t1 = slc_array;

                
                clear myo_t1_eroded
                for i = 1:size(mask_myocardium_3D, 3)
                    myo_t1_eroded(:,:,i) = imerode(mask_myocardium_3D(:,:,i), se);
                end

                vol_img_3D = T1_struct.vol_img_3D;
                roi_in_myo_t1 = myo_t1_eroded .* freeROIMask_3D;
                remote_in_myo_t1 = myo_t1_eroded .* myoRefMask_3D;
                roi_t1 = roi_in_myo_t1 .* vol_img_3D;
                remote_t1 = remote_in_myo_t1 .* vol_img_3D;
                t1 = vol_img_3D;
                myo_t1 = myo_t1_eroded;

                metrics_t1(n).TimePoints(tp_count).time_point = time_point;
                metrics_t1(n).TimePoints(tp_count).mean_t1_roi = mean(nonzeros(t1 .* roi_in_myo_t1));
                metrics_t1(n).TimePoints(tp_count).mean_t1_remote = mean(nonzeros(t1 .* remote_in_myo_t1));
                metrics_t1(n).TimePoints(tp_count).std_t1_roi = std(nonzeros(t1 .* roi_in_myo_t1));
                metrics_t1(n).TimePoints(tp_count).std_t1_remote = std(nonzeros(t1 .* remote_in_myo_t1));
                metrics_t1(n).TimePoints(tp_count).skewness_t1_roi = skewness(nonzeros(t1 .* roi_in_myo_t1));
                metrics_t1(n).TimePoints(tp_count).skewness_t1_remote = skewness(nonzeros(t1 .* remote_in_myo_t1));
                metrics_t1(n).TimePoints(tp_count).kurtosis_t1_roi = kurtosis(nonzeros(t1 .* roi_in_myo_t1))-3;
                metrics_t1(n).TimePoints(tp_count).kurtosis_t1_remote = kurtosis(nonzeros(t1 .* remote_in_myo_t1))-3;

                metrics_t1(n).TimePoints(tp_count).min_t1_roi = min(min(nonzeros(t1 .* roi_in_myo_t1)));
                metrics_t1(n).TimePoints(tp_count).max_t1_roi = max(max(nonzeros(t1 .* roi_in_myo_t1)));
                metrics_t1(n).TimePoints(tp_count).min_t1_remote = min(min(nonzeros(t1 .* remote_in_myo_t1)));
                metrics_t1(n).TimePoints(tp_count).max_t1_remote = max(max(nonzeros(t1 .* remote_in_myo_t1)));

                CI_roi = ConfidenceInterval(nonzeros(t1 .* roi_in_myo_t1));
                CI_remote = ConfidenceInterval(nonzeros(t1 .* remote_in_myo_t1));
                metrics_t1(n).TimePoints(tp_count).ci_roi_lower = metrics_t1(n).TimePoints(tp_count).mean_t1_roi + CI_roi(1);
                metrics_t1(n).TimePoints(tp_count).ci_roi_upper = metrics_t1(n).TimePoints(tp_count).mean_t1_roi + CI_roi(2);
                metrics_t1(n).TimePoints(tp_count).ci_remote_lower = metrics_t1(n).TimePoints(tp_count).mean_t1_remote + CI_remote(1);
                metrics_t1(n).TimePoints(tp_count).ci_remote_upper = metrics_t1(n).TimePoints(tp_count).mean_t1_remote + CI_remote(2);

                metrics_t1(n).TimePoints(tp_count).roi_1_percentile = prctile(nonzeros(t1 .* roi_in_myo_t1), 1);
                metrics_t1(n).TimePoints(tp_count).roi_99_percentile = prctile(nonzeros(t1 .* roi_in_myo_t1), 99);
                metrics_t1(n).TimePoints(tp_count).remote_1_percentile = prctile(nonzeros(t1 .* remote_in_myo_t1), 1);
                metrics_t1(n).TimePoints(tp_count).remote_99_percentile = prctile(nonzeros(t1 .* remote_in_myo_t1), 99);


                metrics_t1(n).TimePoints(tp_count).SliceAnalysis = struct;
                metrics_t1(n).TimePoints(tp_count).HeteroAnalysis = struct;

                temp_roi = (t1 .* roi_in_myo_t1);
                temp_remote = (t1 .* remote_in_myo_t1);
                temp_myo = (t1 .* myo_t1);
                t1_norm_remote = zeros(size(t1));
                t1_norm_roi = zeros(size(t1));

            %for slc = 1:size(roi_in_myo_t1,3)
            for slc = 3:3
                %temp_norm_roi = temp_roi(:,:,slc);
                %temp_norm_roi(temp_norm_roi == 0) = nan;
                % temp_norm_roi = uint8((temp_roi(:,:,slc) - min(min((nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc))))))./ (max(max(nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc))))-min(min((nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc)))))) * 256);
                % temp_norm_remote = uint8((temp_remote(:,:,slc) - min(min((nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc))))))./ (max(max(nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc))))-min(min((nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc)))))) * 256);

                % Hard-coded for Patient data T1 mapping
                temp_norm_roi = (temp_roi(:,:,slc) - min_roi_value) ./ (max_roi_value - min_roi_value);
                temp_norm_remote = (temp_remote(:,:,slc) - min_roi_value) ./ (max_roi_value - min_roi_value);

                metrics_t1(n).TimePoints(tp_count).SliceAnalysis(slc).mean_t1_roi = mean(nonzeros(temp_norm_roi .* roi_in_myo_t1(:,:,slc)));
                metrics_t1(n).TimePoints(tp_count).SliceAnalysis(slc).mean_t1_remote = mean(nonzeros(temp_norm_remote .* remote_in_myo_t1(:,:,slc)));
                metrics_t1(n).TimePoints(tp_count).SliceAnalysis(slc).std_t1_roi = std(nonzeros(temp_norm_roi .* roi_in_myo_t1(:,:,slc)));
                metrics_t1(n).TimePoints(tp_count).SliceAnalysis(slc).std_t1_remote = std(nonzeros(temp_norm_remote .* remote_in_myo_t1(:,:,slc)));
                metrics_t1(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_t1_roi = skewness(nonzeros(temp_norm_roi .* roi_in_myo_t1(:,:,slc)));
                metrics_t1(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_t1_remote = skewness(nonzeros(temp_norm_remote .* remote_in_myo_t1(:,:,slc)));
                metrics_t1(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_t1_roi = kurtosis(nonzeros(temp_norm_roi .* roi_in_myo_t1(:,:,slc)))-3;
                metrics_t1(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_t1_remote = kurtosis(nonzeros(temp_norm_remote .* remote_in_myo_t1(:,:,slc)))-3;

                t1_norm_roi(:,:,slc) = temp_norm_roi;
                t1_norm_remote(:,:,slc) = temp_norm_remote;
                %metrics_t1(n).TimePoints(tp_count).HeteroAnalysis(slc).hetero_roi = Func_Hetero_Analysis(slc, roi_in_myo_t1, t1_norm_roi);
                %metrics_t1(n).TimePoints(tp_count).HeteroAnalysis(slc).hetero_remote = Func_Hetero_Analysis(slc, remote_in_myo_t1, t1_norm_remote);

                %bipolar = temp_roi(:,:,slc) - mean(nonzeros(temp_remote(:,:,slc)));
                %bipolar(bipolar < 0) * 2

                %roi_in_myo_t1;

                % bipolar = temp_roi(:,:,slc) - mean(nonzeros(temp_remote(:,:,slc)));
                % weighted_map = double(bipolar < 0) .* 2 .* roi_in_myo_t1(:,:,slc) .* bipolar + double(bipolar >= 0) .* roi_in_myo_t1(:,:,slc) .* bipolar;
                % 
                % bipolar_remote = temp_remote(:,:,slc) - mean(nonzeros(temp_remote(:,:,slc)));
                % weighted_map_remote = double(bipolar_remote < 0) .* 2 .* remote_in_myo_t1(:,:,slc) .* bipolar_remote + double(bipolar_remote >= 0) .* remote_in_myo_t1(:,:,slc) .* bipolar_remote;

                bipolar = temp_roi(:,:,slc) - mean(nonzeros(temp_remote(:,:,slc)));
                % weighted_map = double(bipolar < 0) .* 2 .* roi_in_myo_t1(:,:,slc) .* bipolar + double(bipolar >= 0) .* roi_in_myo_t1(:,:,slc) .* bipolar;
                weighted_map = double(bipolar < 0) .* roi_in_myo_t1(:,:,slc) .* bipolar + double(bipolar >= 0) .* roi_in_myo_t1(:,:,slc) .* bipolar;

                bipolar_remote = temp_remote(:,:,slc) - mean(nonzeros(temp_remote(:,:,slc)));
                % weighted_map_remote = double(bipolar_remote < 0) .* 2 .* remote_in_myo_t1(:,:,slc) .* bipolar_remote + double(bipolar_remote >= 0) .* remote_in_myo_t1(:,:,slc) .* bipolar_remote;
                weighted_map_remote = double(bipolar_remote < 0) .* remote_in_myo_t1(:,:,slc) .* bipolar_remote + double(bipolar_remote >= 0) .* remote_in_myo_t1(:,:,slc) .* bipolar_remote;


                % lb = 2*(400 - mean(nonzeros(temp_remote(:,:,slc))));
                lb = (min_roi_value - mean(nonzeros(temp_remote(:,:,slc))));
                ub = max_roi_value - mean(nonzeros(temp_remote(:,:,slc)));
                weighted_map(weighted_map < lb) = lb;
                weighted_map(weighted_map > ub) = ub;

                weighted_map_remote(weighted_map_remote < lb) = lb;
                weighted_map_remote(weighted_map_remote > ub) = ub;

                temp_norm_roi = uint8((weighted_map - lb)./ (ub-lb) * 256);
                temp_norm_remote = uint8((weighted_map_remote - lb)./ (ub-lb) * 256);


                non_roi_value = uint8((0 - lb) ./ (ub - lb) * 256);

                % temp_norm_roi = uint8((temp_roi(:,:,slc) - 500) ./ (1500 - 500) * 256);
                % temp_norm_remote = uint8((temp_remote(:,:,slc) - 500) ./ (1500 - 500) * 256);

                bins = 20;
                t1_norm_roi_slc = temp_norm_roi;
                t1_norm_roi_slc(t1_norm_roi_slc == non_roi_value) = [];


                figure(); imhist(t1_norm_roi_slc, bins);
                p_roi = imhist(t1_norm_roi_slc, bins);
                nonZeros = find(p_roi);
                len = length((nonZeros));
                pNonZeros = zeros(1,len);

                for i = 1:len
                    pNonZeros(i) = p_roi(nonZeros(i));
                end

                % normalize pNonZeros so that sum(p) is one.
                pNonZeros = pNonZeros ./ sum(p_roi);
                E_roi = -sum(pNonZeros.*log2(pNonZeros));

                text(50, 10, cat(2, 'Entropy = ', num2str(E_roi)));
                fname = cat(2, 'T1entropy_Histo_ROI_', time_point, '_', num2str(slc), 'unweighted_original_', num2str(bins), 'bins.png');
                saveas(gcf,cat(2, name_save_dir, '/', time_point, '/', fname));


                non_remote_value = uint8((0 - lb) ./ (ub - lb) * 256);

                t1_norm_remote_slc = temp_norm_remote;
                t1_norm_remote_slc(t1_norm_remote_slc == non_remote_value) = [];
                figure(); imhist(t1_norm_remote_slc,bins);
                p_remote = imhist(t1_norm_remote_slc,bins);

                nonZeros = find(p_remote);
                len = length((nonZeros));
                pNonZeros = zeros(1,len);

                for i = 1:len
                    pNonZeros(i) = p_remote(nonZeros(i));
                end

                % normalize pNonZeros so that sum(p) is one.
                pNonZeros = pNonZeros ./ sum(p_remote);
                E_remote = -sum(pNonZeros.*log2(pNonZeros));

                text(50, 4, cat(2, 'Entropy = ', num2str(E_remote)));
                fname = cat(2, 'T1entropy_Histo_Remote_', time_point, '_', num2str(slc), 'unweighted_original_', num2str(bins), 'bins.png');
                saveas(gcf,cat(2, name_save_dir, '/', time_point, '/', fname));


                metrics_t1(n).TimePoints(tp_count).HeteroAnalysis(slc).entropy_roi = E_roi;
                metrics_t1(n).TimePoints(tp_count).HeteroAnalysis(slc).entropy_remote = E_remote;

                L = bwlabel(roi_in_myo_t1(:,:,slc));
                if length(unique(L)) > 2
                    L_cell{count} = cat(2, name, ' ', time_point, ' Slice ', num2str(slc));
                    count = count + 1;
                end
            end
            


            
        end
    end
    close all;
end


%% FF and R2star
all_names = {'484060000001', '484060000003', '484060000004', '484060000008', '484060000009', '484060000010', '484060000011',...
    '484060000012', '484060000013', '484060000015', '484060000016', '484060000023', '484060000018', '484060000020', '484060000017',...
    '484060000021', '484060000030', '484060000029', '484060000032', '484060000035', '484060000036', '484060000038', '484060000039', ...
    '484060000040', '484060000031', '484060000037'};

% % For BL2
% all_names = {'484060000001', '484060000003', '484060000004', '484060000008', '484060000009', '484060000010', '484060000011',...
%     '484060000012', '484060000013', '484060000015', '484060000016', '484060000023', '484060000018', '484060000020', '484060000017',...
%     '484060000021', '484060000029', '484060000032', '484060000035', '484060000036', '484060000038', '484060000039', ...
%     '484060000040', '484060000031', '484060000037'};

%time_points = {'BL'};
%time_points = {'BL2'};
time_points = {'FU'};

mean_ff_roi_array = [];
sd_ff_roi_array = [];
mean_ff_remote_array = [];
sd_ff_remote_array = [];

mean_r2star_roi_array = [];
sd_r2star_roi_array = [];
mean_r2star_remote_array = [];
sd_r2star_remote_array = [];

ff_pixel_roi_array = [];
r2star_pixel_roi_array = [];
ff_pixel_remote_array = [];
r2star_pixel_remote_array = [];

name_label = {};
slice_count = 1;
vec = @(x) x(:);
slc_loc_label = {};


% For dichotomize hemo+ and hemo-
mean_ff_roi_array_hemo_n = [];
sd_ff_roi_array_hemo_n = [];
mean_ff_remote_array_hemo_n = [];
sd_ff_remote_array_hemo_n = [];

mean_r2star_roi_array_hemo_n = [];
sd_r2star_roi_array_hemo_n = [];
mean_r2star_remote_array_hemo_n = [];
sd_r2star_remote_array_hemo_n = [];

% Positive
mean_ff_roi_array_hemo_p = [];
sd_ff_roi_array_hemo_p = [];
mean_ff_remote_array_hemo_p = [];
sd_ff_remote_array_hemo_p = [];

mean_r2star_roi_array_hemo_p = [];
sd_r2star_roi_array_hemo_p = [];
mean_r2star_remote_array_hemo_p = [];
sd_r2star_remote_array_hemo_p = [];

name_label_hemo_n = {};
slice_count_hemo_n = 1;
name_label_hemo_p = {};
slice_count_hemo_p = 1;

se = strel('disk', 1);
nhood = [1 1 1; 1 1 1; 1 1 1];

exclude_pixels_flag = 1;
exclude_pixels_flag = 0;

metrics_ff = struct;

%for n = 1:length(Names)
for n = starting_point:starting_point
%for n = 13:13
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

    metrics_ff(n).name = name;
    metrics_ff(n).TimePoints = struct;

    %if any(contains(excel_names, name))
    if any(contains(all_names, name))
        name_for_table_searching = insertAfter(name, 6, '-');
        row = find(contains(id_array,name_for_table_searching));

        IMH_cell = table2cell(T(row, 13)); % IMH
        IMH = IMH_cell{1};
        for tp = 1:length(time_points)
        %for tp = 1:1
            time_point = time_points{end-tp+1};
            tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', time_point,  '/');
            if ~exist(tp_dir, 'dir')
                disp(cat(2, 'No folder at: ', name, ' ', time_point));
            else
                % T1
                tp_count = tp_count+1;

                % T1 slice location
                load(cat(2, tp_dir, label_t1, '/', label_t1, '_SliceLoc.mat'));
                %[slc_array_t1, idx_reordered] = sort(slc_array);
                slc_array_t1 = slc_array;


                % FF
                myo_glob = glob(cat(2, tp_dir, label_t2star, '/', anatomy_label{5}, '/*'));
                roi_glob = glob(cat(2, tp_dir, label_t2star, '/',anatomy_label{3}, '/*'));
                remote_glob = glob(cat(2, tp_dir, label_t2star, '/',anatomy_label{6}, '/*'));

                load(myo_glob{1});
                load(roi_glob{1});
                load(remote_glob{1});
                load(cat(2, tp_dir, label_t2star, '/', label_t2star, '_Index.mat')); % glob_names
                load(cat(2, tp_dir, label_t2star, '/', label_t2star, '_SliceLoc.mat')); % slc_array
                slc_array_t2star = slc_array;

                %idx_reordered = Func_AlignSliceLoc(slc_array_t1, slc_array_t2star);


                ff_map = cell(1, length(glob_names));
                num_array = ExtractNum(glob_names);
                for f = 1:length(ff_map)
                    % series_name = glob_names{f};
                    ff_glob = glob(cat(2,  base_dir, '/FF_Data/',  name, '/', time_point, '/*_', num2str(num_array(f)), '.mat'));
                    ff_map{f} = load(ff_glob{1}, 'fwmc_ff');
                end

                % convert ff_map to matrix
                ff = zeros(size(ff_map{1}.fwmc_ff,1), size(ff_map{1}.fwmc_ff, 2), length(ff_map));
                for f = 1:length(ff_map)
                    ff(:,:,f) = ff_map{f}.fwmc_ff;
                end

                % R2star Map
                r2star_map = cell(1, length(glob_names));
                for f = 1:length(r2star_map)
                    r2star_glob = glob(cat(2,  base_dir, '/FF_Data/',  name, '/', time_point, '/*_', num2str(num_array(f)), '.mat'));
                    r2star_map{f} = load(r2star_glob{1}, 'fwmc_r2star');
                end

                % convert ff_map to matrix
                r2star = zeros(size(r2star_map{1}.fwmc_r2star,1), size(r2star_map{1}.fwmc_r2star, 2), length(r2star_map));
                for f = 1:length(r2star_map)
                    r2star(:,:,f) = r2star_map{f}.fwmc_r2star;
                end

                mask_myocardium_3D = imerode(mask_myocardium_3D, se);

                roi_in_myo_ff = mask_myocardium_3D .* freeROIMask_3D;
                remote_in_myo_ff = mask_myocardium_3D .* myoRefMask_3D;
                roi_ff = roi_in_myo_ff .* ff;
                remote_ff = remote_in_myo_ff .* ff;
                myo_ff = mask_myocardium_3D;


                % TODO add QC for myocardium map here
                if exclude_pixels_flag == 1
                    f_exclusion_pts = cat(2, name_data_save_dir, '/ExclusionPts_', name, '_', time_point, '.mat');
                    figure('Position', [100 100 800 600]);

                    if exist(f_exclusion_pts, 'file')
                        load(f_exclusion_pts);


                        for slc = 1:size(roi_in_myo_ff,3)
                            exclude_mask = zeros(size(mask_myocardium_3D,1),size(mask_myocardium_3D,2));

                            x_array = exclusion_points{slc, 1};
                            y_array = exclusion_points{slc, 2};

                            for pts = 1:length(y_array)
                                exclude_mask(round(x_array(pts)),round(y_array(pts))) = 1;
                            end

                            exclude_mask = imdilate(exclude_mask, nhood);

                            freeROIMask_3D(:,:,slc) = freeROIMask_3D(:,:,slc) .* ~exclude_mask;

                            roi_in_myo_ff(:,:,slc) = mask_myocardium_3D(:,:,slc) .* freeROIMask_3D(:,:,slc);
                            roi_ff(:,:,slc) = roi_in_myo_ff(:,:,slc) .* ff(:,:,slc);
                            % myo_ff(:,:,slc) = mask_myocardium_3D(:,:,slc);
                        end

                    else
                        exclusion_points = cell(size(roi_in_myo_ff,3), 2);


                        for slc = 1:size(roi_in_myo_ff,3)
                            if any(any(roi_in_myo_ff(:,:,slc)))
                                exclude_mask = zeros(size(mask_myocardium_3D,1),size(mask_myocardium_3D,2));
                                flag = 1;
                                y_array = [];
                                x_array = [];
                                count = 0;
                                while flag
                                    % crop image
                                    stats = regionprops(mask_myocardium_3D(:,:,slc));
                                    centroid = stats.Centroid;
                                    X = round(centroid(1)-size(roi_in_myo_ff,2)/8);
                                    Y = round(centroid(2)-size(roi_in_myo_ff,1)/8);
                                    w = size(roi_in_myo_ff,2)/4;
                                    h = size(roi_in_myo_ff,1)/4;

                                    roi_in_myo_ff_crop = imcrop(roi_in_myo_ff(:,:,slc), [X,Y,w,h]);
                                    r2star_crop = imcrop(r2star(:,:,slc), [X,Y,w,h]);
                                    roi_ff_crop = imcrop(roi_ff(:,:,slc), [X,Y,w,h]);
                                    myo_ff_crop = imcrop(myo_ff(:,:,slc), [X,Y,w,h]);

                                    roi_in_myo_ff_nan = roi_in_myo_ff_crop;
                                    roi_in_myo_ff_nan(roi_in_myo_ff_crop == 0) = nan;

                                    ax1 = axes;
                                    imagesc(r2star_crop); pbaspect([size(r2star_crop, 2), size(r2star_crop, 1) 1]);
                                    colormap(ax1, 'gray');
                                    set(ax1, 'xticklabel', []); set(ax1, 'yticklabel', []);

                                    ax2 = axes;
                                    imagesc(ax2, roi_in_myo_ff_nan .* roi_ff_crop, 'AlphaData', myo_ff_crop);
                                    pbaspect([size(r2star_crop, 2), size(r2star_crop, 1) 1]); colormap(ax2, 'cool');
                                    caxis(ax1, [0 100]); caxis(ax2, [-2 20]); linkprop([ax1 ax2], 'Position');
                                    ax2.Visible = 'off';
                                    if count == 0
                                        cb = colorbar; title(cb, 'FF (%)');
                                    end

                                    [y,x] = getpts;

                                    if length(y) == 1
                                        txt = 'y';
                                    else
                                        y = y + X - 1;
                                        x = x + Y - 1;

                                        y_array = [y_array; y];
                                        x_array = [x_array; x];

                                        for pts = 1:length(y)
                                            exclude_mask(round(x_array(pts)),round(y_array(pts))) = 1;
                                        end

                                        exclude_mask = imdilate(exclude_mask, nhood);

                                        freeROIMask_3D(:,:,slc) = freeROIMask_3D(:,:,slc) .* ~exclude_mask;

                                        roi_in_myo_ff(:,:,slc) = mask_myocardium_3D(:,:,slc) .* freeROIMask_3D(:,:,slc);
                                        roi_ff(:,:,slc) = roi_in_myo_ff(:,:,slc) .* ff(:,:,slc);
                                        % myo_ff(:,:,slc) = mask_myocardium_3D(:,:,slc);

                                        r2star_crop = imcrop(r2star(:,:,slc), [X,Y,w,h]);
                                        roi_ff_crop = imcrop(roi_ff(:,:,slc), [X,Y,w,h]);

                                        imagesc(roi_in_myo_ff_nan .* roi_ff_crop);
                                        pbaspect([size(r2star_crop, 2), size(r2star_crop, 1) 1]); caxis([-2 10]);

                                        prompt = 'Are you happy??? (y/n): ';
                                        txt = input(prompt,"s");
                                    end

                                    if strcmp(txt, 'y')
                                        flag = 0;
                                    elseif strcmp(txt, 'n')
                                        flag = 1;
                                        count = count + 1;
                                    end
                                end
                                exclusion_points{slc, 1} = x_array;
                                exclusion_points{slc, 2} = y_array;
                            end
                        end

                        save(f_exclusion_pts, 'exclusion_points');
                    end
                end

                %  need to reorder for analysis
                % addpath('../T1NFF/');

                % Special cases: 23, 
                idx_reordered = Func_AlignSliceLoc(slc_array_t1, slc_array_t2star);
                ff = ff(:,:,idx_reordered);
                myo_ff = myo_ff(:,:,idx_reordered);
                remote_ff = remote_ff(:,:,idx_reordered);
                roi_ff = roi_ff(:,:,idx_reordered);
                remote_in_myo_ff = remote_in_myo_ff(:,:,idx_reordered);
                roi_in_myo_ff = roi_in_myo_ff(:,:,idx_reordered);


                roi_in_myo_r2star = mask_myocardium_3D .* freeROIMask_3D;
                remote_in_myo_r2star = mask_myocardium_3D .* myoRefMask_3D;
                roi_r2star = roi_in_myo_r2star .* r2star;
                remote_r2star = remote_in_myo_r2star .* r2star;
                myo_r2star = mask_myocardium_3D;


                r2star = r2star(:,:,idx_reordered);
                myo_r2star = myo_r2star(:,:,idx_reordered);
                remote_r2star = remote_r2star(:,:,idx_reordered);
                roi_r2star = roi_r2star(:,:,idx_reordered);
                remote_in_myo_r2star = remote_in_myo_r2star(:,:,idx_reordered);
                roi_in_myo_r2star = roi_in_myo_r2star(:,:,idx_reordered);


                tp_dir2 = cat(2, name_save_dir, '/', time_point, '/');
                if ~exist(tp_dir2, 'dir')
                    mkdir(tp_dir2);
                end


                center_mask_fname = cat(2, name_data_save_dir, '/CenterLine_', name, '_', time_point, '.mat');

                strat_fname = cat(2, name_data_save_dir, '/FF_Stratify_', name, '_', time_point, '.mat');

                if ~exist(center_mask_fname, 'file')
                    center_mask_ff = zeros(size(roi_in_myo_ff));
                    BW_skel = zeros(size(roi_in_myo_ff));
                    %caxis_rg = [0 1];
                    %for i = 1:size(roi_in_myo_ff, 3)
                    % Draw center line on myocardium ff
                    %disp(cat(2, 'Center line roi: ', name, '  ', time_point, 'Slice', num2str(i)));

                    %moving = myo_ff(:,:,i);
                    %center_mask_ff(:,:,i) = Func_DrawCenterLine(moving, caxis_rg);
                    %BW_skel(:,:,i) = bwmorph(moving, 'skel', Inf);
                    %center_mask_ff(:,:,i) = imfill(BW_skel(:,:,i), 'hole');
                    %end

                    figure();
                    sz = ceil(sqrt(size(roi_in_myo_ff, 3)));
                    se = strel('disk', 1);
                    for i = 1:size(roi_in_myo_ff, 3)
                        subplot(sz,sz,i);
                        myo_ff_eroded = imerode(myo_ff(:,:,i), se);
                        BW_skel(:,:,i) = bwmorph(myo_ff_eroded, 'skel', Inf);
                        center_mask_ff(:,:,i) = imfill(BW_skel(:,:,i), 'hole');

                        % Found 18D16 8WK is not a closed shape (Hard-Coded)
                        % Felicity 6MO
                        %                     if (n == 14 && tp == 8 && i == 3) || (n == 11 && tp == 4 && i == 3)
                        %                         center_mask_ff(:,:,3) = bwconvhull(center_mask_ff(:,:,3));
                        %                     end
                        epi = myo_ff_eroded - center_mask_ff(:,:,i) > 0;
                        endo = center_mask_ff(:,:,i) + myo_ff_eroded > 1;
                        imagesc(endo*2 + epi);
                        colormap(brewermap([],'*RdYlBu'));
                    end
                    save(center_mask_fname, 'center_mask_ff');
                    saveas(gcf, cat(2, tp_dir2, 'CenterLineMask.png'))
                else
                    load(center_mask_fname);
                end

                % status = status_check(n).status(tp_count,:);
                % AHA Segment
                Segn = 50;
                Groove = 0;

                r2star(r2star > 200) = 200;
                r2star(r2star < 0) = 0;
                ff(ff > 100) = 100;
                ff(ff < 0) = 0;
                %roi_in_myo_r2star(roi_in_myo_r2star == 0) = nan;
                %roi_in_myo_ff(roi_in_myo_ff == 0) = nan;
                %remote_in_myo_r2star(remote_in_myo_r2star == 0) = nan;
                %remote_in_myo_ff(remote_in_myo_ff == 0) = nan; % ???

                % ROI-wise
                r2star_roi_masked = roi_in_myo_r2star .* r2star;
                ff_roi_masked = roi_in_myo_ff .* ff;
                r2star_remote_masked = remote_in_myo_r2star .* r2star;
                ff_remote_masked = remote_in_myo_ff .* ff;

                % remove 0 and 100s for ff map
                ff_roi_masked_px = ff_roi_masked;
                r2star_roi_masked_px = r2star_roi_masked;
                ff_roi_masked_px(ff_roi_masked == 0) = nan;
                ff_roi_masked_px(ff_roi_masked == 100) = nan;
                r2star_roi_masked_px(ff_roi_masked == 0) = nan;
                r2star_roi_masked_px(ff_roi_masked == 100) = nan;


                for slc = 1:size(roi_in_myo_ff,3)

                    % Hard-coded for Patient data
                    metrics_ff(n).TimePoints(tp_count).SliceAnalysis(slc).mean_ff_roi = mean(vec(ff_roi_masked_px(:,:,slc)), 'omitnan');
                    metrics_ff(n).TimePoints(tp_count).SliceAnalysis(slc).mean_ff_remote = mean(vec(nonzeros(ff_remote_masked(:,:,slc))), 'omitnan');
                    metrics_ff(n).TimePoints(tp_count).SliceAnalysis(slc).std_ff_roi = std(vec(ff_roi_masked_px(:,:,slc)), 'omitnan');
                    metrics_ff(n).TimePoints(tp_count).SliceAnalysis(slc).std_ff_remote = std(vec(nonzeros(ff_remote_masked(:,:,slc))), 'omitnan');

                    metrics_ff(n).TimePoints(tp_count).SliceAnalysis(slc).mean_r2star_roi = mean(vec(r2star_roi_masked_px(:,:,slc)), 'omitnan');
                    metrics_ff(n).TimePoints(tp_count).SliceAnalysis(slc).mean_r2star_remote = mean(vec(nonzeros(r2star_remote_masked(:,:,slc))), 'omitnan');
                    metrics_ff(n).TimePoints(tp_count).SliceAnalysis(slc).std_r2star_roi = std(vec(r2star_roi_masked_px(:,:,slc)), 'omitnan');
                    metrics_ff(n).TimePoints(tp_count).SliceAnalysis(slc).std_r2star_remote = std(vec(nonzeros(r2star_remote_masked(:,:,slc))), 'omitnan');

                end

                figure();
                for slc = 1:size(r2star_roi_masked, 3)
                    mean_r2star_roi_array = [mean_r2star_roi_array, mean(vec(r2star_roi_masked_px(:,:,slc)), 'omitnan')];
                    sd_r2star_roi_array = [sd_r2star_roi_array, std(vec(r2star_roi_masked_px(:,:,slc)), 'omitnan')];
                    mean_ff_roi_array = [mean_ff_roi_array, mean(vec(ff_roi_masked_px(:,:,slc)), 'omitnan')];
                    sd_ff_roi_array = [sd_ff_roi_array, std(vec(ff_roi_masked_px(:,:,slc)), 'omitnan')];
                    mean_r2star_remote_array = [mean_r2star_remote_array, mean(vec(r2star_remote_masked(:,:,slc)), 'omitnan')];
                    sd_r2star_remote_array = [sd_r2star_remote_array, std(vec(r2star_remote_masked(:,:,slc)), 'omitnan')];
                    mean_ff_remote_array = [mean_ff_remote_array, mean(vec(ff_remote_masked(:,:,slc)), 'omitnan')];
                    sd_ff_remote_array = [sd_ff_remote_array, std(vec(ff_remote_masked(:,:,slc)), 'omitnan')];

                    name_label{slice_count} = cat(2, name, '_', time_point, '_Slice', num2str(slc));
                    slc_loc_label{slice_count} = cat(2, num2str(slc_array_t2star(slc)));
                    slice_count = slice_count + 1;

                    subplot(1,2,1); imagesc(ff_roi_masked_px(:,:,slc)); axis image; caxis([0 20]); colorbar;
                    subplot(1,2,2); imagesc(r2star_roi_masked_px(:,:,slc)); axis image; caxis([0 100]); colorbar;
                end


                % Dichotomize into hemo+ and hemo-
                if strcmp(IMH, '-')
                    % ROI
                    for slc = 1:size(r2star_roi_masked, 3)
                        mean_r2star_roi_array_hemo_n = [mean_r2star_roi_array_hemo_n, mean(vec(r2star_roi_masked_px(:,:,slc)), 'omitnan')];
                        sd_r2star_roi_array_hemo_n = [sd_r2star_roi_array_hemo_n, std(vec(r2star_roi_masked_px(:,:,slc)), 'omitnan')];
                        mean_ff_roi_array_hemo_n = [mean_ff_roi_array_hemo_n, mean(vec(ff_roi_masked_px(:,:,slc)), 'omitnan')];
                        sd_ff_roi_array_hemo_n = [sd_ff_roi_array_hemo_n, std(vec(ff_roi_masked_px(:,:,slc)), 'omitnan')];
                        mean_r2star_remote_array_hemo_n = [mean_r2star_remote_array_hemo_n, mean(vec(r2star_remote_masked(:,:,slc)), 'omitnan')];
                        sd_r2star_remote_array_hemo_n = [sd_r2star_remote_array_hemo_n, std(vec(r2star_remote_masked(:,:,slc)), 'omitnan')];
                        mean_ff_remote_array_hemo_n = [mean_ff_remote_array_hemo_n, mean(vec(ff_remote_masked(:,:,slc)), 'omitnan')];
                        sd_ff_remote_array_hemo_n = [sd_ff_remote_array_hemo_n, std(vec(ff_remote_masked(:,:,slc)), 'omitnan')];

                        name_label_hemo_n{slice_count_hemo_n} = cat(2, name, '_', time_point, '_Slice', num2str(slc));
                        slice_count_hemo_n = slice_count_hemo_n + 1;
                    end

                elseif strcmp(IMH, '+')
                    % ROI
                    for slc = 1:size(r2star_roi_masked, 3)
                        mean_r2star_roi_array_hemo_p = [mean_r2star_roi_array_hemo_p, mean(vec(r2star_roi_masked_px(:,:,slc)), 'omitnan')];
                        sd_r2star_roi_array_hemo_p = [sd_r2star_roi_array_hemo_p, std(vec(r2star_roi_masked_px(:,:,slc)), 'omitnan')];
                        mean_ff_roi_array_hemo_p = [mean_ff_roi_array_hemo_p, mean(vec(ff_roi_masked_px(:,:,slc)), 'omitnan')];
                        sd_ff_roi_array_hemo_p = [sd_ff_roi_array_hemo_p, std(vec(ff_roi_masked_px(:,:,slc)), 'omitnan')];
                        mean_r2star_remote_array_hemo_p = [mean_r2star_remote_array_hemo_p, mean(vec(r2star_remote_masked(:,:,slc)), 'omitnan')];
                        sd_r2star_remote_array_hemo_p = [sd_r2star_remote_array_hemo_p, std(vec(r2star_remote_masked(:,:,slc)), 'omitnan')];
                        mean_ff_remote_array_hemo_p = [mean_ff_remote_array_hemo_p, mean(vec(ff_remote_masked(:,:,slc)), 'omitnan')];
                        sd_ff_remote_array_hemo_p = [sd_ff_remote_array_hemo_p, std(vec(ff_remote_masked(:,:,slc)), 'omitnan')];

                        name_label_hemo_p{slice_count_hemo_p} = cat(2, name, '_', time_point, '_Slice', num2str(slc));
                        slice_count_hemo_p = slice_count_hemo_p + 1;
                    end
                end
            end

            close all;
        end
    end % find the excellent subjects
end


%% T2

all_names = {'484060000001', '484060000010', '484060000012', '484060000016', '484060000018', '484060000030', '484060000029', '484060000036', '484060000031', '484060000037'};
mean_t2_roi_array = [];
sd_t2_roi_array = [];
mean_t2_remote_array = [];
sd_t2_remote_array = [];

t2_pixel_roi_array = [];
t2_pixel_remote_array = [];

name_label = {};
slice_count = 1;
vec = @(x) x(:);
slc_loc_label = {};

% For dichotomize hemo+ and hemo-
mean_t2_roi_array_hemo_n = [];
sd_t2_roi_array_hemo_n = [];
mean_t2_remote_array_hemo_n = [];
sd_t2_remote_array_hemo_n = [];

% Positive
mean_t2_roi_array_hemo_p = [];
sd_t2_roi_array_hemo_p = [];
mean_t2_remote_array_hemo_p = [];
sd_t2_remote_array_hemo_p = [];

name_label_hemo_n = {};
slice_count_hemo_n = 1;
name_label_hemo_p = {};
slice_count_hemo_p = 1;

se = strel('disk', 1);
nhood = [1 1 1; 1 1 1; 1 1 1];

exclude_pixels_flag = 1;

for n = 1:length(Names)
%for n = starting_point:starting_point
%for n = 3:3
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
    % tp_count = 0;

    %if any(contains(excel_names, name))
    if any(contains(all_names, name))
        name_for_table_searching = insertAfter(name, 6, '-');
        row = find(contains(id_array,name_for_table_searching));

        IMH_cell = table2cell(T(row, 13)); % IMH
        IMH = IMH_cell{1};
        for tp = 1:length(time_points)
        %for tp = 1:1
            time_point = time_points{end-tp+1};
            tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', time_point,  '/');
            if ~exist(tp_dir, 'dir')
                disp(cat(2, 'No folder at: ', name, ' ', time_point));
            else
                % T2
                myo_glob = glob(cat(2, tp_dir, label_t2, '/', anatomy_label{5}, '/*'));
                roi_glob = glob(cat(2, tp_dir, label_t2, '/',anatomy_label{3}, '/*'));
                remote_glob = glob(cat(2, tp_dir, label_t2, '/',anatomy_label{6}, '/*'));

                load(myo_glob{1});
                load(roi_glob{1});
                load(remote_glob{1});
                % load(cat(2, tp_dir, label_t2, '/', label_t2, '_Index.mat')); % glob_names
                load(cat(2, tp_dir, label_t2, '/', label_t2, '_SliceLoc.mat')); % slc_array
                load(cat(2, tp_dir, label_t2, '/', label_t2, '_vol_img_3D.mat'));
                slc_array_t2star = slc_array;

                t2 = T1_struct.vol_img_3D * 0.1;
                mask_myocardium_3D = imerode(mask_myocardium_3D, se);

                roi_in_myo_t2 = mask_myocardium_3D .* freeROIMask_3D;
                remote_in_myo_t2 = mask_myocardium_3D .* myoRefMask_3D;
                roi_t2 = roi_in_myo_t2 .* t2;
                remote_t2 = remote_in_myo_t2 .* t2;
                myo_t2 = mask_myocardium_3D;
              

                % % No need to reorder for analysis
                % addpath('../T1NFF/');
                % idx_reordered = Func_AlignSliceLoc(slc_array_t1, slc_array);
                % t2 = t2(:,:,idx_reordered);
                % myo_t2 = myo_t2(:,:,idx_reordered);
                % remote_t2 = remote_t2(:,:,idx_reordered);
                % roi_t2 = roi_t2(:,:,idx_reordered);
                % remote_in_myo_t2 = remote_in_myo_t2(:,:,idx_reordered);
                % roi_in_myo_t2 = roi_in_myo_t2(:,:,idx_reordered);


                tp_dir2 = cat(2, name_save_dir, '/', time_point, '/');
                if ~exist(tp_dir2, 'dir')
                    mkdir(tp_dir2);
                end


                center_mask_fname = cat(2, name_data_save_dir, '/CenterLine_', name, '_', time_point, '.mat');
                strat_fname = cat(2, name_data_save_dir, '/FF_Stratify_', name, '_', time_point, '.mat');


                % status = status_check(n).status(tp_count,:);
                % AHA Segment
                Segn = 50;
                Groove = 0;

                t2(t2 < 0) = 0;
                %roi_in_myo_r2star(roi_in_myo_r2star == 0) = nan;
                %roi_in_myo_ff(roi_in_myo_ff == 0) = nan;
                %remote_in_myo_r2star(remote_in_myo_r2star == 0) = nan;
                %remote_in_myo_ff(remote_in_myo_ff == 0) = nan; % ???

                % ROI-wise
                t2_roi_masked = roi_in_myo_t2 .* t2;
                t2_remote_masked = remote_in_myo_t2 .* t2;

                % remove 0 and 100s for ff map
                t2_roi_masked_px = t2_roi_masked;
                t2_roi_masked_px(t2_roi_masked == 0) = nan;


                figure();
                for slc = 1:size(t2_roi_masked, 3)
                    mean_t2_roi_array = [mean_t2_roi_array, mean(vec(t2_roi_masked_px(:,:,slc)), 'omitnan')];
                    sd_t2_roi_array = [sd_t2_roi_array, std(vec(t2_roi_masked_px(:,:,slc)), 'omitnan')];
                    mean_t2_remote_array = [mean_t2_remote_array, mean(vec(t2_remote_masked(:,:,slc)), 'omitnan')];
                    sd_t2_remote_array = [sd_t2_remote_array, std(vec(t2_remote_masked(:,:,slc)), 'omitnan')];

                    name_label{slice_count} = cat(2, name, '_', time_point, '_Slice', num2str(slc));
                    slc_loc_label{slice_count} = cat(2, num2str(slc_array(slc)));
                    slice_count = slice_count + 1;

                    imagesc(t2_roi_masked_px(:,:,slc)); axis image; caxis([0 100]); colorbar;
                end


                % Dichotomize into hemo+ and hemo-
                if strcmp(IMH, '-')
                    % ROI
                    for slc = 1:size(t2_roi_masked, 3)
                       
                        mean_t2_roi_array_hemo_n = [mean_t2_roi_array_hemo_n, mean(vec(t2_roi_masked_px(:,:,slc)), 'omitnan')];
                        sd_t2_roi_array_hemo_n = [sd_t2_roi_array_hemo_n, std(vec(t2_roi_masked_px(:,:,slc)), 'omitnan')];
                        mean_t2_remote_array_hemo_n = [mean_t2_remote_array_hemo_n, mean(vec(t2_remote_masked(:,:,slc)), 'omitnan')];
                        sd_t2_remote_array_hemo_n = [sd_t2_remote_array_hemo_n, std(vec(t2_remote_masked(:,:,slc)), 'omitnan')];

                        name_label_hemo_n{slice_count_hemo_n} = cat(2, name, '_', time_point, '_Slice', num2str(slc));
                        slice_count_hemo_n = slice_count_hemo_n + 1;
                    end

                elseif strcmp(IMH, '+')
                    % ROI
                    for slc = 1:size(t2_roi_masked, 3)

                        mean_t2_roi_array_hemo_p = [mean_t2_roi_array_hemo_p, mean(vec(t2_roi_masked_px(:,:,slc)), 'omitnan')];
                        sd_t2_roi_array_hemo_p = [sd_t2_roi_array_hemo_p, std(vec(t2_roi_masked_px(:,:,slc)), 'omitnan')];
                        mean_t2_remote_array_hemo_p = [mean_t2_remote_array_hemo_p, mean(vec(t2_remote_masked(:,:,slc)), 'omitnan')];
                        sd_t2_remote_array_hemo_p = [sd_t2_remote_array_hemo_p, std(vec(t2_remote_masked(:,:,slc)), 'omitnan')];

                        name_label_hemo_p{slice_count_hemo_p} = cat(2, name, '_', time_point, '_Slice', num2str(slc));
                        slice_count_hemo_p = slice_count_hemo_p + 1;
                    end
                end
            end

            close all;
        end
    end % find the excellent subjects
end

%% Save as struct (FU) T2
metrics_to_save = struct;
metrics_to_save.name_label = name_label;
metrics_to_save.slc_loc_label = slc_loc_label;
metrics_to_save.name_label_hemo_p = name_label_hemo_p;
metrics_to_save.name_label_hemo_n = name_label_hemo_n;

metrics_to_save.ROI = struct;
metrics_to_save.ROI.mean_t2_roi_array = mean_t2_roi_array;
metrics_to_save.ROI.sd_t2_roi_array = sd_t2_roi_array;
metrics_to_save.ROI.mean_t2_remote_array = mean_t2_remote_array;
metrics_to_save.ROI.sd_t2_remote_array = sd_t2_remote_array;

% For dichotomize hemo+ and hemo-
metrics_to_save.ROI.mean_t2_roi_array_hemo_n = mean_t2_roi_array_hemo_n;
metrics_to_save.ROI.sd_t2_roi_array_hemo_n = sd_t2_roi_array_hemo_n;
metrics_to_save.ROI.mean_t2_remote_array_hemo_n = mean_t2_remote_array_hemo_n;
metrics_to_save.ROI.sd_t2_remote_array_hemo_n = sd_t2_remote_array_hemo_n;

% Positive
metrics_to_save.ROI.mean_t2_roi_array_hemo_p = mean_t2_roi_array_hemo_p;
metrics_to_save.ROI.sd_t2_roi_array_hemo_p = sd_t2_roi_array_hemo_p;
metrics_to_save.ROI.mean_t2_remote_array_hemo_p = mean_t2_remote_array_hemo_p;
metrics_to_save.ROI.sd_t2_remote_array_hemo_p = sd_t2_remote_array_hemo_p;


%% Save as struct (FU) (Go back to main body and reiterate BL)
metrics_to_save = struct;
metrics_to_save.name_label = name_label;
metrics_to_save.slc_loc_label = slc_loc_label;
metrics_to_save.name_label_hemo_p = name_label_hemo_p;
metrics_to_save.name_label_hemo_n = name_label_hemo_n;

metrics_to_save.ROI = struct;
metrics_to_save.ROI.mean_ff_roi_array = mean_ff_roi_array;
metrics_to_save.ROI.sd_ff_roi_array = sd_ff_roi_array;
metrics_to_save.ROI.mean_ff_remote_array = mean_ff_remote_array;
metrics_to_save.ROI.sd_ff_remote_array = sd_ff_remote_array;

metrics_to_save.ROI.mean_r2star_roi_array = mean_r2star_roi_array;
metrics_to_save.ROI.sd_r2star_roi_array = sd_r2star_roi_array;
metrics_to_save.ROI.mean_r2star_remote_array = mean_r2star_remote_array;
metrics_to_save.ROI.sd_r2star_remote_array = sd_r2star_remote_array;

% For dichotomize hemo+ and hemo-
metrics_to_save.ROI.mean_ff_roi_array_hemo_n = mean_ff_roi_array_hemo_n;
metrics_to_save.ROI.sd_ff_roi_array_hemo_n = sd_ff_roi_array_hemo_n;
metrics_to_save.ROI.mean_ff_remote_array_hemo_n = mean_ff_remote_array_hemo_n;
metrics_to_save.ROI.sd_ff_remote_array_hemo_n = sd_ff_remote_array_hemo_n;

metrics_to_save.ROI.mean_r2star_roi_array_hemo_n = mean_r2star_roi_array_hemo_n;
metrics_to_save.ROI.sd_r2star_roi_array_hemo_n = sd_r2star_roi_array_hemo_n;
metrics_to_save.ROI.mean_r2star_remote_array_hemo_n = mean_r2star_remote_array_hemo_n;
metrics_to_save.ROI.sd_r2star_remote_array_hemo_n = sd_r2star_remote_array_hemo_n;

% Positive
metrics_to_save.ROI.mean_ff_roi_array_hemo_p = mean_ff_roi_array_hemo_p;
metrics_to_save.ROI.sd_ff_roi_array_hemo_p = sd_ff_roi_array_hemo_p;
metrics_to_save.ROI.mean_ff_remote_array_hemo_p = mean_ff_remote_array_hemo_p;
metrics_to_save.ROI.sd_ff_remote_array_hemo_p = sd_ff_remote_array_hemo_p;

metrics_to_save.ROI.mean_r2star_roi_array_hemo_p = mean_r2star_roi_array_hemo_p;
metrics_to_save.ROI.sd_r2star_roi_array_hemo_p = sd_r2star_roi_array_hemo_p;
metrics_to_save.ROI.mean_r2star_remote_array_hemo_p = mean_r2star_remote_array_hemo_p;
metrics_to_save.ROI.sd_r2star_remote_array_hemo_p = sd_r2star_remote_array_hemo_p;

% save(cat(2, metrics_save_dir, 'Metrics_FatNIron_Analysis.mat'), '-struct', 'metrics_to_save');


%% T2
se = strel('disk', 1);
metrics_t2 = struct;
L_cell = {};
count = 1;
%for n = 1:length(Names)
for n = starting_point:starting_point
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

    metrics_t1(n).name = name;
    metrics_t1(n).TimePoints = struct;

    for tp = 1:length(time_points)
        %for tp = 2:2
        time_point = time_points{end-tp+1};
        tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', time_point,  '/');
        if ~exist(cat(2, tp_dir, label_t1, '/'), 'dir')
            disp(cat(2, 'No folder at: ', name, ' ', time_point));
        else

            tp_dir2 = cat(2, name_save_dir, '/', name, '/', time_point, '/');
            if ~exist(tp_dir2, 'dir')
                mkdir(tp_dir2);
            end

                % T2
                myo_glob = glob(cat(2, tp_dir, label_t2, '/', anatomy_label{5}, '/*'));
                roi_glob = glob(cat(2, tp_dir, label_t2, '/',anatomy_label{3}, '/*'));
                remote_glob = glob(cat(2, tp_dir, label_t2, '/',anatomy_label{6}, '/*'));

                load(cat(2, tp_dir, label_t1, '/', label_t2, '_vol_img_3D.mat'));
                load(myo_glob{1});
                load(roi_glob{1});
                load(remote_glob{1});
                load(cat(2, tp_dir, label_t2, '/', label_t2, '_SliceLoc.mat'));
                slc_array_t1 = slc_array;

                clear myo_t2_eroded
                for i = 1:size(mask_myocardium_3D, 3)
                    myo_t2_eroded(:,:,i) = imerode(mask_myocardium_3D(:,:,i), se);
                end

                
                roi_in_myo_t2 = myo_t2_eroded .* freeROIMask_3D;
                remote_in_myo_t2 = myo_t2_eroded .* myoRefMask_3D;
                roi_t2 = roi_in_myo_t2 .* vol_img_3D;
                remote_t2 = remote_in_myo_t2 .* vol_img_3D;
                t2 = vol_img_3D;
                myo_t2 = myo_t2_eroded;

                metrics_t2(n).TimePoints(tp_count).time_point = time_point;
                metrics_t2(n).TimePoints(tp_count).mean_t2_roi = mean(nonzeros(t2 .* roi_in_myo_t2));
                metrics_t2(n).TimePoints(tp_count).mean_t2_remote = mean(nonzeros(t2 .* remote_in_myo_t2));
                metrics_t2(n).TimePoints(tp_count).std_t2_roi = std(nonzeros(t2 .* roi_in_myo_t2));
                metrics_t2(n).TimePoints(tp_count).std_t2_remote = std(nonzeros(t2 .* remote_in_myo_t2));
                metrics_t2(n).TimePoints(tp_count).skewness_t2_roi = skewness(nonzeros(t2 .* roi_in_myo_t2));
                metrics_t2(n).TimePoints(tp_count).skewness_t2_remote = skewness(nonzeros(t2 .* remote_in_myo_t2));
                metrics_t2(n).TimePoints(tp_count).kurtosis_t2_roi = kurtosis(nonzeros(t2 .* roi_in_myo_t2))-3;
                metrics_t2(n).TimePoints(tp_count).kurtosis_t2_remote = kurtosis(nonzeros(t2 .* remote_in_myo_t2))-3;

                metrics_t2(n).TimePoints(tp_count).SliceAnalysis = struct;
                metrics_t2(n).TimePoints(tp_count).HeteroAnalysis = struct;

                temp_roi = (t2 .* roi_in_myo_t2);
                temp_remote = (t2 .* remote_in_myo_t2);
                temp_myo = (t2 .* myo_t2);
                t2_norm_remote = zeros(size(t2));
                t2_norm_roi = zeros(size(t2));

                for slc = 1:size(roi_in_myo_t2,3)
                    %temp_norm_roi = temp_roi(:,:,slc);
                    %temp_norm_roi(temp_norm_roi == 0) = nan;
                    % temp_norm_roi = uint8((temp_roi(:,:,slc) - min(min((nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc))))))./ (max(max(nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc))))-min(min((nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc)))))) * 256);
                    % temp_norm_remote = uint8((temp_remote(:,:,slc) - min(min((nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc))))))./ (max(max(nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc))))-min(min((nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc)))))) * 256);

                    % Hard-coded for Patient data T1 mapping
                    temp_norm_roi = (temp_roi(:,:,slc) - 100) ./ (100 - 0);
                    temp_norm_remote = (temp_remote(:,:,slc) - 100) ./ (100 - 0);

                    metrics_t2(n).TimePoints(tp_count).SliceAnalysis(slc).mean_t2_roi = mean(nonzeros(temp_norm_roi .* roi_in_myo_t2(:,:,slc)));
                    metrics_t2(n).TimePoints(tp_count).SliceAnalysis(slc).mean_t2_remote = mean(nonzeros(temp_norm_remote .* remote_in_myo_t2(:,:,slc)));
                    metrics_t2(n).TimePoints(tp_count).SliceAnalysis(slc).std_t2_roi = std(nonzeros(temp_norm_roi .* roi_in_myo_t2(:,:,slc)));
                    metrics_t2(n).TimePoints(tp_count).SliceAnalysis(slc).std_t2_remote = std(nonzeros(temp_norm_remote .* remote_in_myo_t2(:,:,slc)));
                    metrics_t2(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_t2_roi = skewness(nonzeros(temp_norm_roi .* roi_in_myo_t2(:,:,slc)));
                    metrics_t2(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_t2_remote = skewness(nonzeros(temp_norm_remote .* remote_in_myo_t2(:,:,slc)));
                    metrics_t2(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_t2_roi = kurtosis(nonzeros(temp_norm_roi .* roi_in_myo_t2(:,:,slc)))-3;
                    metrics_t2(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_t2_remote = kurtosis(nonzeros(temp_norm_remote .* remote_in_myo_t2(:,:,slc)))-3;

                    t2_norm_roi(:,:,slc) = temp_norm_roi;
                    t2_norm_remote(:,:,slc) = temp_norm_remote;
                    metrics_t2(n).TimePoints(tp_count).HeteroAnalysis(slc).hetero_roi = Func_Hetero_Analysis(slc, roi_in_myo_t2, t2_norm_roi);
                    metrics_t2(n).TimePoints(tp_count).HeteroAnalysis(slc).hetero_remote = Func_Hetero_Analysis(slc, remote_in_myo_t2, t2_norm_remote);

                    temp_norm_roi = uint8((temp_roi(:,:,slc) - 100) ./ (100 - 0) * 256);
                    temp_norm_remote = uint8((temp_remote(:,:,slc) - 100) ./ (100 - 0) * 256);

                    t1_norm_roi_slc = temp_norm_roi;
                    t1_norm_roi_slc(t1_norm_roi_slc == 0) = [];
                    %figure(); imhist(lge_norm_roi_slc,20);
                    p_roi = imhist(t1_norm_roi_slc,20);
                    nonZeros = find(p_roi);
                    len = length((nonZeros));
                    pNonZeros = zeros(1,len);

                    for i = 1:len
                        pNonZeros(i) = p_roi(nonZeros(i));
                    end

                    % normalize pNonZeros so that sum(p) is one.
                    pNonZeros = pNonZeros ./ sum(p_roi);
                    E_roi = -sum(pNonZeros.*log2(pNonZeros));


                    t1_norm_remote_slc = temp_norm_remote;
                    t1_norm_remote_slc(t1_norm_remote_slc == 0) = [];
                    %figure(); imhist(lge_norm_remote_slc,20);
                    p_remote = imhist(t1_norm_remote_slc,20);

                    nonZeros = find(p_remote);
                    len = length((nonZeros));
                    pNonZeros = zeros(1,len);

                    for i = 1:len
                        pNonZeros(i) = p_remote(nonZeros(i));
                    end

                    % normalize pNonZeros so that sum(p) is one.
                    pNonZeros = pNonZeros ./ sum(p_remote);
                    E_remote = -sum(pNonZeros.*log2(pNonZeros));


                    metrics_t2(n).TimePoints(tp_count).HeteroAnalysis(slc).entropy_roi = E_roi;
                    metrics_t2(n).TimePoints(tp_count).HeteroAnalysis(slc).entropy_remote = E_remote;

                    L = bwlabel(roi_in_myo_t1(:,:,slc));
                    if length(unique(L)) > 2
                        L_cell{count} = cat(2, name, ' ', time_point, ' Slice ', num2str(slc));
                        count = count + 1;
                    end
                end

        end
    end
end


%% Hypothetical 
r = normrnd(3,60,[1,100]) + 127;
R = uint8(r);
figure(); imhist(R, bins);
p_roi = imhist(R, bins);
nonZeros = find(p_roi);
len = length((nonZeros));
pNonZeros = zeros(1,len);

for i = 1:len
    pNonZeros(i) = p_roi(nonZeros(i));
end

% normalize pNonZeros so that sum(p) is one.
pNonZeros = pNonZeros ./ sum(p_roi);
E_roi = -sum(pNonZeros.*log2(pNonZeros))