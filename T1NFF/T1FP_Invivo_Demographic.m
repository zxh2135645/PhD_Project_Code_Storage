clear all;
close all;
%% 
%% Demographic in the MI region (SD,Skewness,Kurtosis)
addpath('../function/');
addpath('../AHA16Segment/');
addpath('../function/demon_registration_version_8f_winOS/');
base_dir = uigetdir;
contour_glob = glob(cat(2, base_dir, '/ContourData/*'));
%Names = ExtractNames(contour_glob);
Names = {'Merry', 'Ryn', 'Mojave', 'Sahara', 'ZZ', 'Tina', 'Sunny', 'Queenie', 'Hope', 'Gobi', 'Felicity', 'Evelyn', '18D15', '18D16', '11D05', '11D26', '11D33'};
%time_points = {'6MO', '9MO', '1YR', '15YR'};
time_points = {'7D', '8WK', '12WK', '14WK' '4MO', '6MO', '9MO', '1YR', '15YR'};
%time_points = {'7D', '8WK'};
%time_points = {'3D'};


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
% sequence_label = {'T1_PostCon', 'T2star', 'LGE'}; % PostCon T1 entropy
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

%%
% 07/27/2023 modified
se = strel('disk', 1);
metrics = struct;
L_cell = {};
count = 1;

min_roi_value = 1000;
max_roi_value = 1000;

vec=@(x) x(:);
qc_flag = 0;

% This is defined for postCon T1 map

% min_roi = 100;
% max_roi = 1500;

min_roi = 750;
max_roi = 1750;

%for n = 1:length(Names)
for n = 14:14
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
    %for tp = 1:length(time_points)
    for tp = 9:9
        time_point = time_points{end-tp+1};
        tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', name, '_', time_point,  '/');
        if ~exist(cat(2, tp_dir, label_t1, '/'), 'dir')
            disp(cat(2, 'No folder at: ', name, ' ', time_point));
        else
            % T1
            tp_count = tp_count+1;
            myo_glob = glob(cat(2, tp_dir, label_t1, '/', anatomy_label{5}, '/*'));
            roi_glob = glob(cat(2, tp_dir, label_t1, '/',anatomy_label{3}, '/*'));
            remote_glob = glob(cat(2, tp_dir, label_t1, '/',anatomy_label{6}, '/*'));
            
            load(cat(2, tp_dir, label_t1, '/', label_t1, '_vol_img_3D.mat'));
            load(myo_glob{1});
            load(roi_glob{1});
            load(remote_glob{1});
            load(cat(2, tp_dir, label_t1, '/', label_t1, '_SliceLoc.mat'));
            [slc_array_t1, idx_reordered] = sort(slc_array);

            clear myo_t1_eroded
            for i = 1:size(mask_myocardium_3D, 3)
                myo_t1_eroded(:,:,i) = imerode(mask_myocardium_3D(:,:,i), se);
            end

            % roi_in_myo_t1 = mask_myocardium_3D .* freeROIMask_3D;
            % remote_in_myo_t1 = mask_myocardium_3D .* myoRefMask_3D;
            roi_in_myo_t1 = myo_t1_eroded .* freeROIMask_3D;
            remote_in_myo_t1 = myo_t1_eroded .* myoRefMask_3D;
            roi_t1 = roi_in_myo_t1 .* vol_img_3D;
            remote_t1 = remote_in_myo_t1 .* vol_img_3D;
            t1 = vol_img_3D;
            % myo_t1 = mask_myocardium_3D;
            myo_t1 = myo_t1_eroded;
            
            roi_in_myo_t1 = roi_in_myo_t1(:,:,idx_reordered);
            remote_in_myo_t1 = remote_in_myo_t1(:,:,idx_reordered);
            roi_t1 = roi_t1(:,:,idx_reordered);
            remote_t1 = remote_t1(:,:,idx_reordered);
            t1 = t1(:,:,idx_reordered);
            myo_t1 = myo_t1(:,:,idx_reordered);
            
            % FF
            myo_glob = glob(cat(2, tp_dir, label_t2star, '/', anatomy_label{5}, '/*'));
            roi_glob = glob(cat(2, tp_dir, label_t2star, '/',anatomy_label{3}, '/*'));
            remote_glob = glob(cat(2, tp_dir, label_t2star, '/',anatomy_label{6}, '/*'));
            
            hemo_core_flag = 0;

            load(myo_glob{1});
            load(roi_glob{1});
            if length(roi_glob) == 2
                load(roi_glob{2});
                hemo_core_flag = 1;
            end
            load(remote_glob{1});
            load(cat(2, tp_dir, label_t2star, '/', label_t2star, '_Index.mat'));
            load(cat(2, tp_dir, label_t2star, '/', label_t2star, '_SliceLoc.mat'));
            slc_array_ff = slc_array;
            
            ff_map = cell(1, length(glob_names));
            for f = 1:length(ff_map)
                ff_map{f} = load(cat(2,  base_dir, '/FF_Data/',  name, '/', name, '_', time_point, '_', glob_names{f}, '.mat'), 'fwmc_ff');
            end
            
            % convert ff_map to matrix
            ff = zeros(size(ff_map{1}.fwmc_ff,1), size(ff_map{1}.fwmc_ff, 2), length(ff_map));
            for f = 1:length(ff_map)
                ff(:,:,f) = ff_map{f}.fwmc_ff;
            end
            
            clear myo_ff_eroded
            for i = 1:size(mask_myocardium_3D, 3)
                myo_ff_eroded(:,:,i) = imerode(mask_myocardium_3D(:,:,i), se);
            end

            % roi_in_myo_ff = mask_myocardium_3D .* freeROIMask_3D;
            % remote_in_myo_ff = mask_myocardium_3D .* myoRefMask_3D;
            roi_in_myo_ff = myo_ff_eroded .* freeROIMask_3D;
            remote_in_myo_ff = myo_ff_eroded .* myoRefMask_3D;
            roi_ff = roi_in_myo_ff .* ff;
            remote_ff = remote_in_myo_ff .* ff;
            % myo_ff = mask_myocardium_3D;
            myo_ff = myo_ff_eroded;

            if hemo_core_flag == 1
                roi_in_myo_ff_te8 = myo_ff_eroded .* freeROIMask_3D_te8;
                roi_ff_te8 = roi_in_myo_ff_te8 .* ff;
            end

            idx_reordered = Func_AlignSliceLoc(slc_array_t1, slc_array_ff);
            % Special cases:
            if strcmp(name, '18D15') && strcmp(time_point, '9MO')
                idx_reordered = [3 2 1];
            end

            ff = ff(:,:,idx_reordered);
            myo_ff = myo_ff(:,:,idx_reordered);
            remote_ff = remote_ff(:,:,idx_reordered);
            roi_ff = roi_ff(:,:,idx_reordered);
            remote_in_myo_ff = remote_in_myo_ff(:,:,idx_reordered);
            roi_in_myo_ff = roi_in_myo_ff(:,:,idx_reordered);
            
            if hemo_core_flag == 1
                roi_ff_te8 = roi_ff_te8(:,:,idx_reordered);
                roi_in_myo_ff_te8 = roi_in_myo_ff_te8(:,:,idx_reordered);
            end

            % R2star Map
            r2star_map = cell(1, length(glob_names));
            for f = 1:length(r2star_map)
                r2star_map{f} = load(cat(2,  base_dir, '/FF_Data/',  name, '/', name, '_', time_point, '_', glob_names{f}, '.mat'), 'fwmc_r2star');
            end
            
            % convert ff_map to matrix
            r2star = zeros(size(r2star_map{1}.fwmc_r2star,1), size(r2star_map{2}.fwmc_r2star, 2), length(r2star_map));
            for f = 1:length(r2star_map)
                r2star(:,:,f) = r2star_map{f}.fwmc_r2star;
            end
            
            
            % roi_in_myo_r2star = mask_myocardium_3D .* freeROIMask_3D;
            % remote_in_myo_r2star = mask_myocardium_3D .* myoRefMask_3D;
            roi_in_myo_r2star = myo_ff_eroded .* freeROIMask_3D;
            remote_in_myo_r2star = myo_ff_eroded .* myoRefMask_3D;
            roi_r2star = roi_in_myo_r2star .* r2star;
            remote_r2star = remote_in_myo_r2star .* r2star;
            % myo_r2star = mask_myocardium_3D;
            myo_r2star = myo_ff_eroded;
            
            if hemo_core_flag == 1
                roi_in_myo_r2star_te8 = myo_ff_eroded .* freeROIMask_3D_te8;
                roi_r2star_te8 = roi_in_myo_r2star_te8 .* r2star;
            end

            r2star = r2star(:,:,idx_reordered);
            myo_r2star = myo_r2star(:,:,idx_reordered);
            remote_r2star = remote_r2star(:,:,idx_reordered);
            roi_r2star = roi_r2star(:,:,idx_reordered);
            remote_in_myo_r2star = remote_in_myo_r2star(:,:,idx_reordered);
            roi_in_myo_r2star = roi_in_myo_r2star(:,:,idx_reordered);
            

            if hemo_core_flag == 1
                roi_r2star_te8 = roi_r2star_te8(:,:,idx_reordered);
                roi_in_myo_r2star_te8 = roi_in_myo_r2star_te8(:,:,idx_reordered);
            end

            tp_dir2 = cat(2, name_save_dir, '/', name, '_', time_point, '/');
            if ~exist(tp_dir2, 'dir')
                mkdir(tp_dir2);
            end
            %LR_mdl_fname = cat(2, name_data_save_dir, '/LinearRegression_', name, '_', time_point, '.mat');
           
            %if ~exist(center_mask_fname, 'file')
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
                    if (n == 14 && tp == 8 && i == 3) || (n == 11 && tp == 4 && i == 3)
                        center_mask_ff(:,:,3) = bwconvhull(center_mask_ff(:,:,3));
                    end
                    epi = myo_ff_eroded - center_mask_ff(:,:,i) > 0;
                    endo = center_mask_ff(:,:,i) + myo_ff_eroded > 1;
                    imagesc(endo*2 + epi);
                    colormap(brewermap([],'*RdYlBu'));
                end
                %save(center_mask_fname, 'center_mask_ff');
                %saveas(gcf, cat(2, tp_dir2, 'CenterLineMask.png'))

                if qc_flag == 1
                    status = status_check(n).status(tp_count,:);
                end
               
            % AHA Segment


            if (strcmp(name, '18D18') && strcmp(time_point, '9MO'))
                % WHY 
                % Some issue with status?
                disp(cat(2, 'Skipped: ', name, ' ', time_point))
            else

                if qc_flag == 1 % See if applies QC or not; Only works for Day3
                    exclude_idx = find(status == 0); % Some how previous version doesn't work for Tina
                else
                    exclude_idx = [];
                end

                if (strcmp(name, 'Tina') && strcmp(time_point, '6MO'))
                    exclude_idx = [];
                elseif (strcmp(name, 'Ryn') && strcmp(time_point, '8WK'))
                    exclude_idx = [];
                elseif (strcmp(name, 'Gobi') && strcmp(time_point, '1YR'))
                    exclude_idx = [];
                elseif (strcmp(name, 'Gobi') && strcmp(time_point, '9MO'))
                    exclude_idx = [];
                elseif (strcmp(name, 'Gobi') && strcmp(time_point, '7D'))
                    exclude_idx = [];
                elseif (strcmp(name, 'Felicity') && strcmp(time_point, '7D'))
                    exclude_idx = [];
                elseif (strcmp(name, '18D16') && strcmp(time_point, '9MO'))
                    exclude_idx = [];
                elseif (strcmp(name, '18D16') && strcmp(time_point, '7D'))
                    exclude_idx = [];
                end

                t1(:,:,exclude_idx) = [];
                ff(:,:,exclude_idx) = [];
                r2star(:,:,exclude_idx) = [];
                roi_in_myo_t1(:,:,exclude_idx) = [];
                remote_in_myo_t1(:,:,exclude_idx) = [];
                roi_in_myo_ff(:,:,exclude_idx) = [];
                remote_in_myo_ff(:,:,exclude_idx) = [];
                roi_in_myo_r2star(:,:,exclude_idx) = [];
                remote_in_myo_r2star(:,:,exclude_idx) = [];
                myo_t1(:,:,exclude_idx) = [];
                
                if hemo_core_flag == 1
                    roi_in_myo_ff_te8(:,:,exclude_idx) = [];
                    roi_in_myo_r2star_te8(:,:,exclude_idx) = [];
                end

                ff(ff < 0) = 1e-8;
                ff(ff > 100) = 100;
                r2star(r2star<0)=1e-8;
                r2star(r2star>200)=200;
                
                
                metrics(n).TimePoints(tp_count).time_point = time_point;
                metrics(n).TimePoints(tp_count).mean_t1_roi = mean(nonzeros(t1 .* roi_in_myo_t1));
                metrics(n).TimePoints(tp_count).mean_t1_remote = mean(nonzeros(t1 .* remote_in_myo_t1));
                metrics(n).TimePoints(tp_count).std_t1_roi = std(nonzeros(t1 .* roi_in_myo_t1));
                metrics(n).TimePoints(tp_count).std_t1_remote = std(nonzeros(t1 .* remote_in_myo_t1));
                metrics(n).TimePoints(tp_count).skewness_t1_roi = skewness(nonzeros(t1 .* roi_in_myo_t1));
                metrics(n).TimePoints(tp_count).skewness_t1_remote = skewness(nonzeros(t1 .* remote_in_myo_t1));
                metrics(n).TimePoints(tp_count).kurtosis_t1_roi = kurtosis(nonzeros(t1 .* roi_in_myo_t1))-3;
                metrics(n).TimePoints(tp_count).kurtosis_t1_remote = kurtosis(nonzeros(t1 .* remote_in_myo_t1))-3;


                metrics(n).TimePoints(tp_count).min_t1_roi = min(min(nonzeros(t1 .* roi_in_myo_t1)));
                metrics(n).TimePoints(tp_count).max_t1_roi = max(max(nonzeros(t1 .* roi_in_myo_t1)));
                metrics(n).TimePoints(tp_count).min_t1_remote = min(min(nonzeros(t1 .* remote_in_myo_t1)));
                metrics(n).TimePoints(tp_count).max_t1_remote = max(max(nonzeros(t1 .* remote_in_myo_t1)));
                
                CI_roi = ConfidenceInterval(nonzeros(t1 .* roi_in_myo_t1));
                CI_remote = ConfidenceInterval(nonzeros(t1 .* remote_in_myo_t1));
                metrics(n).TimePoints(tp_count).ci_roi_lower = metrics(n).TimePoints(tp_count).mean_t1_roi + CI_roi(1);
                metrics(n).TimePoints(tp_count).ci_roi_upper = metrics(n).TimePoints(tp_count).mean_t1_roi + CI_roi(2);
                metrics(n).TimePoints(tp_count).ci_remote_lower = metrics(n).TimePoints(tp_count).mean_t1_remote + CI_remote(1);
                metrics(n).TimePoints(tp_count).ci_remote_upper = metrics(n).TimePoints(tp_count).mean_t1_remote + CI_remote(2);

                metrics(n).TimePoints(tp_count).roi_1_percentile = prctile(nonzeros(t1 .* roi_in_myo_t1), 1);
                metrics(n).TimePoints(tp_count).roi_99_percentile = prctile(nonzeros(t1 .* roi_in_myo_t1), 99);
                metrics(n).TimePoints(tp_count).remote_1_percentile = prctile(nonzeros(t1 .* remote_in_myo_t1), 1);
                metrics(n).TimePoints(tp_count).remote_99_percentile = prctile(nonzeros(t1 .* remote_in_myo_t1), 99);

                metrics(n).TimePoints(tp_count).mean_ff_roi = mean(nonzeros(ff .* roi_in_myo_ff));
                metrics(n).TimePoints(tp_count).mean_ff_remote = mean(nonzeros(ff .* remote_in_myo_ff));
                metrics(n).TimePoints(tp_count).std_ff_roi = std(nonzeros(ff .* roi_in_myo_ff));
                metrics(n).TimePoints(tp_count).std_ff_remote = std(nonzeros(ff .* remote_in_myo_ff));
                metrics(n).TimePoints(tp_count).skewness_ff_roi = skewness(nonzeros(ff .* roi_in_myo_ff));
                metrics(n).TimePoints(tp_count).skewness_ff_remote = skewness(nonzeros(ff .* remote_in_myo_ff));
                metrics(n).TimePoints(tp_count).kurtosis_ff_roi = kurtosis(nonzeros(ff .* roi_in_myo_ff))-3;
                metrics(n).TimePoints(tp_count).kurtosis_ff_remote = kurtosis(nonzeros(ff .* remote_in_myo_ff))-3;

                metrics(n).TimePoints(tp_count).mean_r2star_roi = mean(nonzeros(r2star .* roi_in_myo_r2star));
                metrics(n).TimePoints(tp_count).mean_r2star_remote = mean(nonzeros(r2star .* remote_in_myo_r2star));
                metrics(n).TimePoints(tp_count).std_r2star_roi = std(nonzeros(r2star .* roi_in_myo_r2star));
                metrics(n).TimePoints(tp_count).std_r2star_remote = std(nonzeros(r2star .* remote_in_myo_r2star));
                metrics(n).TimePoints(tp_count).skewness_r2star_roi = skewness(nonzeros(r2star .* roi_in_myo_r2star));
                metrics(n).TimePoints(tp_count).skewness_r2star_remote = skewness(nonzeros(r2star .* remote_in_myo_r2star));
                metrics(n).TimePoints(tp_count).kurtosis_r2star_roi = kurtosis(nonzeros(r2star .* roi_in_myo_r2star))-3;
                metrics(n).TimePoints(tp_count).kurtosis_r2star_remote = kurtosis(nonzeros(r2star .* remote_in_myo_r2star))-3;

                metrics(n).TimePoints(tp_count).SliceAnalysis = struct;
                metrics(n).TimePoints(tp_count).HeteroAnalysis = struct;

                temp_roi = (t1 .* roi_in_myo_t1);
                temp_remote = (t1 .* remote_in_myo_t1);
                temp_myo = (t1 .* myo_t1);
                t1_norm_roi = (t1 - min(t1(:))) / max(t1(:));
                t1_norm_remote = (t1 - min(t1(:))) / max(t1(:));

               

                for slc = 1:size(roi_in_myo_ff,3)
                    min_temp = min(min((nonzeros(temp_roi(:,:,slc)))));
                    max_temp = max(max((nonzeros(temp_roi(:,:,slc)))));

                    if min_temp < min_roi_value
                        min_roi_value = min_temp;
                    end

                    if max_temp > max_roi_value
                        max_roi_value = max_temp;
                    end
                end

                roi_in_myo_ff_nan = roi_in_myo_ff;
                remote_in_myo_ff_nan = remote_in_myo_ff;
                roi_in_myo_ff_nan(roi_in_myo_ff_nan == 0) = nan;
                remote_in_myo_ff_nan(remote_in_myo_ff_nan == 0) = nan;

                %for slc = 1:size(roi_in_myo_ff,3)
                for slc = 3:3

                    temp_norm_roi = (temp_roi(:,:,slc) - min(min((nonzeros(temp_roi(:,:,slc))))))./ max(max(nonzeros(temp_roi(:,:,slc))));
                    temp_norm_remote = (temp_remote(:,:,slc) - min(min(nonzeros(temp_remote(:,:,slc)))))./ max(max(nonzeros(temp_remote(:,:,slc))));

                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_t1_roi = mean(nonzeros(temp_roi(:,:,slc) .* roi_in_myo_t1(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_t1_remote = mean(nonzeros(temp_remote(:,:,slc) .* remote_in_myo_t1(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_t1_roi = std(nonzeros(temp_roi(:,:,slc) .* roi_in_myo_t1(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_t1_remote = std(nonzeros(temp_remote(:,:,slc) .* remote_in_myo_t1(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_t1_roi = skewness(nonzeros(temp_roi(:,:,slc) .* roi_in_myo_t1(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_t1_remote = skewness(nonzeros(temp_remote(:,:,slc) .* remote_in_myo_t1(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_t1_roi = kurtosis(nonzeros(temp_roi(:,:,slc) .* roi_in_myo_t1(:,:,slc)))-3;
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_t1_remote = kurtosis(nonzeros(temp_remote(:,:,slc) .* remote_in_myo_t1(:,:,slc)))-3;

                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_ff_roi = mean(vec(ff(:,:,slc) .* roi_in_myo_ff_nan(:,:,slc)), 'omitnan');
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_ff_remote = mean(vec(ff(:,:,slc) .* remote_in_myo_ff_nan(:,:,slc)), 'omitnan');
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_ff_roi = std(vec(ff(:,:,slc) .* roi_in_myo_ff_nan(:,:,slc)), 'omitnan');
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_ff_remote = std(vec(ff(:,:,slc) .* remote_in_myo_ff_nan(:,:,slc)), 'omitnan');

                    % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_ff_roi = mean(nonzeros(ff(:,:,slc) .* roi_in_myo_ff(:,:,slc)));
                    % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_ff_remote = mean(nonzeros(ff(:,:,slc) .* remote_in_myo_ff(:,:,slc)));
                    % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_ff_roi = std(nonzeros(ff(:,:,slc) .* roi_in_myo_ff(:,:,slc)));
                    % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_ff_remote = std(nonzeros(ff(:,:,slc) .* remote_in_myo_ff(:,:,slc)));

                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_ff_roi = skewness(nonzeros(ff(:,:,slc) .* roi_in_myo_ff(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_ff_remote = skewness(nonzeros(ff(:,:,slc) .* remote_in_myo_ff(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_ff_roi = kurtosis(nonzeros(ff(:,:,slc) .* roi_in_myo_ff(:,:,slc)))-3;
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_ff_remote = kurtosis(nonzeros(ff(:,:,slc) .* remote_in_myo_ff(:,:,slc)))-3;


                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_r2star_roi = mean(nonzeros(r2star(:,:,slc) .* roi_in_myo_r2star(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_r2star_remote = mean(nonzeros(r2star(:,:,slc) .* remote_in_myo_r2star(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_r2star_roi = std(nonzeros(r2star(:,:,slc) .* roi_in_myo_r2star(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_r2star_remote = std(nonzeros(r2star(:,:,slc) .* remote_in_myo_r2star(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_r2star_roi = skewness(nonzeros(r2star(:,:,slc) .* roi_in_myo_r2star(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_r2star_remote = skewness(nonzeros(r2star(:,:,slc) .* remote_in_myo_r2star(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_r2star_roi = kurtosis(nonzeros(r2star(:,:,slc) .* roi_in_myo_r2star(:,:,slc)))-3;
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_r2star_remote = kurtosis(nonzeros(r2star(:,:,slc) .* remote_in_myo_r2star(:,:,slc)))-3;

                    if hemo_core_flag == 1
                        metrics(n).TimePoints(tp_count).SliceAnalysis(slc).hemo_status = sum(vec(roi_in_myo_ff_te8(:,:,slc))) > 0;

                        metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_ff_roi_te8 = mean(nonzeros(ff(:,:,slc) .* roi_in_myo_ff_te8(:,:,slc)));
                        metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_ff_roi_te8 = std(nonzeros(ff(:,:,slc) .* roi_in_myo_ff_te8(:,:,slc)));
                        metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_ff_roi_te8 = skewness(nonzeros(ff(:,:,slc) .* roi_in_myo_ff_te8(:,:,slc)));
                        metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_ff_roi_te8 = kurtosis(nonzeros(ff(:,:,slc) .* roi_in_myo_ff_te8(:,:,slc)))-3;

                        metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_r2star_roi_te8 = mean(nonzeros(r2star(:,:,slc) .* roi_in_myo_r2star_te8(:,:,slc)));
                        metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_r2star_roi_te8 = std(nonzeros(r2star(:,:,slc) .* roi_in_myo_r2star_te8(:,:,slc)));
                        metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_r2star_roi_te8 = skewness(nonzeros(r2star(:,:,slc) .* roi_in_myo_r2star_te8(:,:,slc)));
                        metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_r2star_roi_te8 = kurtosis(nonzeros(r2star(:,:,slc) .* roi_in_myo_r2star_te8(:,:,slc)))-3;

                        metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_ff_roi_te8_peri = mean(nonzeros(ff(:,:,slc) .* (roi_in_myo_ff(:,:,slc) - roi_in_myo_ff_te8(:,:,slc))));
                        metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_ff_roi_te8_peri = std(nonzeros(ff(:,:,slc) .* (roi_in_myo_ff(:,:,slc) - roi_in_myo_ff_te8(:,:,slc))));
                        metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_ff_roi_te8_peri = skewness(nonzeros(ff(:,:,slc) .* (roi_in_myo_ff(:,:,slc) - roi_in_myo_ff_te8(:,:,slc))));
                        metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_ff_roi_te8_peri = kurtosis(nonzeros(ff(:,:,slc) .* (roi_in_myo_ff(:,:,slc) - roi_in_myo_ff_te8(:,:,slc))))-3;

                        metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_r2star_roi_te8_peri = mean(nonzeros(r2star(:,:,slc) .* (roi_in_myo_ff(:,:,slc) - roi_in_myo_r2star_te8(:,:,slc))));
                        metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_r2star_roi_te8_peri = std(nonzeros(r2star(:,:,slc) .* (roi_in_myo_ff(:,:,slc) - roi_in_myo_r2star_te8(:,:,slc))));
                        metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_r2star_roi_te8_peri = skewness(nonzeros(r2star(:,:,slc) .* (roi_in_myo_ff(:,:,slc) - roi_in_myo_r2star_te8(:,:,slc))));
                        metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_r2star_roi_te8_peri = kurtosis(nonzeros(r2star(:,:,slc) .* (roi_in_myo_ff(:,:,slc) - roi_in_myo_r2star_te8(:,:,slc))))-3;
                    end

                    % t1_norm_roi(:,:,slc) = temp_norm_roi;
                    % t1_norm_remote(:,:,slc) = temp_norm_remote;
                    % metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).hetero_roi = Func_Hetero_Analysis(slc, roi_in_myo_t1, t1_norm_roi);
                    % metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).hetero_remote = Func_Hetero_Analysis(slc, remote_in_myo_t1, t1_norm_remote);


                    % temp_norm_roi = uint8((temp_roi(:,:,slc) - min(min((nonzeros(temp_myo(:,:,slc))))))./ (max(max(nonzeros(temp_myo(:,:,slc))))-min(min((nonzeros(temp_myo(:,:,slc)))))) * 256);
                    % temp_norm_remote = uint8((temp_remote(:,:,slc) - min(min((nonzeros(temp_myo(:,:,slc))))))./ (max(max(nonzeros(temp_myo(:,:,slc))))-min(min((nonzeros(temp_myo(:,:,slc)))))) * 256);

                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).max_t1_roi = max(max(nonzeros(temp_roi(:,:,slc) .* roi_in_myo_t1(:,:,slc))));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).min_t1_roi = min(min(nonzeros(temp_roi(:,:,slc) .* roi_in_myo_t1(:,:,slc))));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).max_t1_remote = max(max(nonzeros(temp_remote(:,:,slc) .* remote_in_myo_t1(:,:,slc))));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).min_t1_remote = min(min(nonzeros(temp_remote(:,:,slc) .* remote_in_myo_t1(:,:,slc))));


                    % Hard-coded for T1 mapping
                    bipolar = temp_roi(:,:,slc) - mean(nonzeros(temp_remote(:,:,slc)));
                    % weighted_map = double(bipolar < 0) .* 2 .* roi_in_myo_t1(:,:,slc) .* bipolar + double(bipolar >= 0) .* roi_in_myo_t1(:,:,slc) .* bipolar;
                    weighted_map = double(bipolar < 0) .* roi_in_myo_t1(:,:,slc) .* bipolar + double(bipolar >= 0) .* roi_in_myo_t1(:,:,slc) .* bipolar;

                    bipolar_remote = temp_remote(:,:,slc) - mean(nonzeros(temp_remote(:,:,slc)));
                    % weighted_map_remote = double(bipolar_remote < 0) .* 2 .* remote_in_myo_t1(:,:,slc) .* bipolar_remote + double(bipolar_remote >= 0) .* remote_in_myo_t1(:,:,slc) .* bipolar_remote;
                    weighted_map_remote = double(bipolar_remote < 0) .* remote_in_myo_t1(:,:,slc) .* bipolar_remote + double(bipolar_remote >= 0) .* remote_in_myo_t1(:,:,slc) .* bipolar_remote;

                    % lb = 2*(400 - mean(nonzeros(temp_remote(:,:,slc))));
                    lb = (min_roi - mean(nonzeros(temp_remote(:,:,slc))));
                    ub = max_roi - mean(nonzeros(temp_remote(:,:,slc)));
                    weighted_map(weighted_map < lb) = lb;
                    weighted_map(weighted_map > ub) = ub;

                    weighted_map_remote(weighted_map_remote < lb) = lb;
                    weighted_map_remote(weighted_map_remote > ub) = ub;

                    %temp_norm_roi = uint8((temp_roi(:,:,slc) - 800)./ (1800-800) * 256);
                    %temp_norm_remote = uint8((temp_remote(:,:,slc) - 800)./ (1800-800) * 256);
      
                    temp_norm_roi = uint8((weighted_map - lb)./ (ub-lb) * 256);
                    % temp_norm_remote = uint8((temp_remote(:,:,slc) - 800)./ (1800-800) * 256);
                    temp_norm_remote = uint8((weighted_map_remote - lb)./ (ub-lb) * 256);


                    non_roi_value = uint8((0 - lb) ./ (ub - lb) * 256);
                    t1_norm_roi_slc = temp_norm_roi;
                    t1_norm_roi_slc(t1_norm_roi_slc ==  non_roi_value) = [];
                    
                    bins = 20;
                   
                    figure(); imhist(t1_norm_roi_slc,bins);
                    p_roi = imhist(t1_norm_roi_slc,bins);
                    nonZeros = find(p_roi);
                    len = length((nonZeros));
                    pNonZeros = zeros(1,len);
                    
                    
                    for i = 1:len
                        if p_roi(nonZeros(i)) > 0 % No removal
                            pNonZeros(i) = p_roi(nonZeros(i));
                        end
                    end

                    % normalize pNonZeros so that sum(p) is one.
                    %pNonZeros = pNonZeros ./ sum(p_roi);
                    pNonZeros = nonzeros(pNonZeros ./ sum(pNonZeros));
                    E_roi = -sum(pNonZeros.*log2(pNonZeros));
                    text(50, 10, cat(2, 'Entropy = ', num2str(E_roi)));
                    fname = cat(2, 'T1entropy_Histo_ROI_', time_point, '_', num2str(slc), 'unweighted_original_', num2str(bins), 'bins.png');
                    saveas(gcf,cat(2, name_save_dir, '/', name, '_', time_point, '/', fname));


                    non_remote_value = uint8((0 - lb) ./ (ub - lb) * 256);
                    t1_norm_remote_slc = temp_norm_remote;
                    t1_norm_remote_slc(t1_norm_remote_slc == non_remote_value) = [];

                    
                    figure(); imhist(t1_norm_remote_slc,bins);
                    p_remote = imhist(t1_norm_remote_slc,bins);

                    nonZeros = find(p_remote);
                    len = length((nonZeros));
                    pNonZeros = zeros(1,len);

                    for i = 1:len
                        if p_remote(nonZeros(i)) > 0
                            pNonZeros(i) = p_remote(nonZeros(i));
                        end
                    end

                    % normalize pNonZeros so that sum(p) is one.
                    %pNonZeros = pNonZeros ./ sum(p_remote);
                    pNonZeros = nonzeros(pNonZeros ./ sum(pNonZeros));
                    E_remote = -sum(pNonZeros.*log2(pNonZeros));

                    text(50, 4, cat(2, 'Entropy = ', num2str(E_remote)));
                    fname = cat(2, 'T1entropy_Histo_Remote_', time_point, '_', num2str(slc), 'unweighted_original_', num2str(bins), 'bins.png');
                    saveas(gcf,cat(2, name_save_dir, '/', name, '_', time_point, '/', fname));

                    metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).entropy_roi = E_roi;
                    metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).entropy_remote = E_remote;

                    L = bwlabel(roi_in_myo_t1(:,:,slc));
                    
                    if length(unique(L)) > 2
                        L_cell{count} = cat(2, name, ' ', time_point, ' Slice ', num2str(slc));
                        count = count + 1;
                    end

                    bin_interval = 256/bins;
                    x = (bin_interval/2 + 0:bin_interval:255).';
                    y = p_roi;
                    f = fit(x, y, 'gauss1');
                    figure();
                    plot(x,y,'o');
                    hold on; hold on; plot(f, x, y);
                    fname = cat(2, 'T1entropy_Histo_ROI_', time_point, '_', num2str(slc), 'Gaussian_', num2str(bins), 'bins.png');
                    saveas(gcf,cat(2, name_save_dir, '/', name, '_', time_point, '/', fname));

                    y_eval = f.a1 * exp(-((x-f.b1)/f.c1).^2);
                    y_roi = round(y - y_eval);
                    p_roi_removal = y_roi(round(y - y_eval)>0);
                    
                    len = length(p_roi_removal);
                    pNonZeros = zeros(1,len);                                      
                    for i = 1:len
                        if p_roi_removal(i) > 0 % No removal
                            pNonZeros(i) = p_roi_removal(i);
                        end
                    end

                    pNonZeros = nonzeros(pNonZeros ./ sum(pNonZeros));
                    E_roi_removal = -sum(pNonZeros.*log2(pNonZeros));

                    bin_interval = 256/bins;
                    x = (bin_interval/2 + 0:bin_interval:255).';
                    y = p_remote;
                    f = fit(x, y, 'gauss1');
                    figure();
                    plot(x,y,'o');
                    hold on; hold on; plot(f, x, y);
                    fname = cat(2, 'T1entropy_Histo_Remote_', time_point, '_', num2str(slc), 'Gaussian_', num2str(bins), 'bins.png');
                    saveas(gcf,cat(2, name_save_dir, '/', name, '_', time_point, '/', fname));


                    y_eval = f.a1 * exp(-((x-f.b1)/f.c1).^2);
                    y_remote = round(y - y_eval);
                    p_remote_removal = y_remote(round(y - y_eval)>0);

                    len = length(p_remote_removal);
                    pNonZeros = zeros(1,len);
                    for i = 1:len
                        if p_remote_removal(i) > 0 % No removal
                            pNonZeros(i) = p_remote_removal(i);
                        end
                    end

                    pNonZeros = nonzeros(pNonZeros ./ sum(pNonZeros));
                    E_remote_removal = -sum(pNonZeros.*log2(pNonZeros));

                    metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).entropy_roi_removal = E_roi_removal;
                    metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).entropy_remote_removal = E_remote_removal;

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
%fname = 'Demographic_Metrics_rim_unaware_normalizedinROIMean_N_Entropy_bin40_weighted_rm2px';
fname = 'Demographic_Metrics_rim_unaware_normalizedinROIMean_N_Entropy_bin40_original';

save(cat(2, data_save_dir, '/', fname), 'metrics');

%% For example pull up

figure(); 
for i = 1:4
    subplot(2,2,i);
    imagesc(t1(:,:,i) .* roi_in_myo_t1(:,:,i));
    %imagesc(t1(:,:,i) .* remote_in_myo_t1(:,:,i));
    std(nonzeros(t1(:,:,i) .* roi_in_myo_t1(:,:,i))) / mean(nonzeros(t1(:,:,i) .* roi_in_myo_t1(:,:,i)));
end

%% Merry 1.5 YR, slice 2,
%% To DO: what's wrong with Sahara 1YR slice 1: need to do shape unaware

slc = 1;

[X,Y] = find(roi_in_myo_t1(:,:,slc));
C_array = 1:length(X);
C = nchoosek(C_array, 2);

%Xm = 59; Ym = 52; Im = 1221;
%Xn = 62; Yn = 67; In = 1383;
Delta_I_ens = zeros(size(C,1),2);

% example of i = 7534
for i = 1:size(C,1)
%for i = 7534:7534
%for i = 3695:3695
%for i = 1288:1288
    Xm = X(C(i,1)); Ym = Y(C(i,1));
    Im = t1(Xm, Ym, slc);
    Xn = X(C(i,2)); Yn = Y(C(i,2));
    In = t1(Xn, Yn, slc);
    %if Xm ~= Xn && Ym ~= Yn
        rmn = sqrt((Xm-Xn).^2 + (Ym-Yn).^2);
        Imn = In - Im;

        [XL, YL] = bresenham(Xm, Ym, Xn, Yn);
        %figure();
        %plot(XL, YL, 'or');

        % I_rmL = zeros(length(XL), 1);
        Delta_I_array = zeros(length(XL), 1);
        L = length(XL);
        count = 0;
        for l = 1:L
            Xl = XL(l); Yl = YL(l);
            Il = t1(Xl, Yl, slc) .* (roi_in_myo_t1(Xl,Yl,slc));

            if Il ~= 0 % To make it shape unaware
                rml = sqrt((Xm-Xl).^2 + (Ym-Yl).^2);
                I_rml = Im + (In - Im)/rmn * rml;
                Delta_I_array(l) = abs(I_rml - Il);
            else
                count = count + 1;
            end
        end

        L = L - count; % To make it shape unaware

        Delta_I = sum(Delta_I_array) / L;
        Delta_I_ens(i,1) = L;
        Delta_I_ens(i,2) = Delta_I;
    %end
end


L_array = unique(Delta_I_ens(:,1));
L_max = max(L_array);
L_array_norm = L_array / L_max;
Delta_I_ens_avg = zeros(length(L_array), 1);
for l = 1:length(L_array)
    L = L_array(l);
    [idx] = find(Delta_I_ens(:,1) == L);
    Delta_I_ens_avg(l) = mean(Delta_I_ens(idx,2));
end
figure();
plot(L_array_norm, Delta_I_ens_avg);
Q = trapz(L_array_norm,Delta_I_ens_avg);
Q
