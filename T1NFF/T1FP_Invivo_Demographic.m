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
time_points = {'6MO', '9MO', '1YR', '15YR'};

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

%%
metrics = struct;
for n = 1:length(Names)
%for n = 6:6
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
    %for tp = 4:4
        time_point = time_points{end-tp+1};
        tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', name, '_', time_point,  '/');
        if ~exist(tp_dir, 'dir')
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
            
            roi_in_myo_t1 = mask_myocardium_3D .* freeROIMask_3D;
            remote_in_myo_t1 = mask_myocardium_3D .* myoRefMask_3D;
            roi_t1 = roi_in_myo_t1 .* vol_img_3D;
            remote_t1 = remote_in_myo_t1 .* vol_img_3D;
            t1 = vol_img_3D;
            myo_t1 = mask_myocardium_3D;
            
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
            
            load(myo_glob{1});
            load(roi_glob{1});
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
            
            roi_in_myo_ff = mask_myocardium_3D .* freeROIMask_3D;
            remote_in_myo_ff = mask_myocardium_3D .* myoRefMask_3D;
            roi_ff = roi_in_myo_ff .* ff;
            remote_ff = remote_in_myo_ff .* ff;
            myo_ff = mask_myocardium_3D;
            
            idx_reordered = Func_AlignSliceLoc(slc_array_t1, slc_array_ff);
            ff = ff(:,:,idx_reordered);
            myo_ff = myo_ff(:,:,idx_reordered);
            remote_ff = remote_ff(:,:,idx_reordered);
            roi_ff = roi_ff(:,:,idx_reordered);
            remote_in_myo_ff = remote_in_myo_ff(:,:,idx_reordered);
            roi_in_myo_ff = roi_in_myo_ff(:,:,idx_reordered);
            
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

            status = status_check(n).status(tp_count,:);
            % AHA Segment

            if (strcmp(name, '18D16') && strcmp(time_point, '9MO'))
                % WHY 
                % Some issue with status?
                disp(cat(2, 'Skipped: ', name, ' ', time_point))
            else
                exclude_idx = find(status == 0);
                t1(:,:,exclude_idx) = [];
                ff(:,:,exclude_idx) = [];
                r2star(:,:,exclude_idx) = [];
                roi_in_myo_t1(:,:,exclude_idx) = [];
                remote_in_myo_t1(:,:,exclude_idx) = [];
                roi_in_myo_ff(:,:,exclude_idx) = [];
                remote_in_myo_ff(:,:,exclude_idx) = [];
                roi_in_myo_r2star(:,:,exclude_idx) = [];
                remote_in_myo_r2star(:,:,exclude_idx) = [];


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
                for slc = 1:size(roi_in_myo_ff,3)
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_t1_roi = mean(nonzeros(t1(:,:,slc) .* roi_in_myo_t1(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_t1_remote = mean(nonzeros(t1(:,:,slc) .* remote_in_myo_t1(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_t1_roi = std(nonzeros(t1(:,:,slc) .* roi_in_myo_t1(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_t1_remote = std(nonzeros(t1(:,:,slc) .* remote_in_myo_t1(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_t1_roi = skewness(nonzeros(t1(:,:,slc) .* roi_in_myo_t1(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_t1_remote = skewness(nonzeros(t1(:,:,slc) .* remote_in_myo_t1(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_t1_roi = kurtosis(nonzeros(t1(:,:,slc) .* roi_in_myo_t1(:,:,slc)))-3;
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_t1_remote = kurtosis(nonzeros(t1(:,:,slc) .* remote_in_myo_t1(:,:,slc)))-3;

                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_ff_roi = mean(nonzeros(ff(:,:,slc) .* roi_in_myo_ff(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_ff_remote = mean(nonzeros(ff(:,:,slc) .* remote_in_myo_ff(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_ff_roi = std(nonzeros(ff(:,:,slc) .* roi_in_myo_ff(:,:,slc)));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_ff_remote = std(nonzeros(ff(:,:,slc) .* remote_in_myo_ff(:,:,slc)));
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

fname = 'Demographic_Metrics';
save(cat(2, data_save_dir, '/', fname), 'metrics');
