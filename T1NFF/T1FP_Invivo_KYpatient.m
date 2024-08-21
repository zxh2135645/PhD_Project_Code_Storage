clear all;
close all;

%% Demographic in the MI region (SD,Skewness,Kurtosis)
addpath('../function/');
addpath('../AHA16Segment/');
addpath('../function/demon_registration_version_8f_winOS/');
base_dir = uigetdir;
contour_glob = glob(cat(2, base_dir, '/ContourData/*'));
Names = ExtractNames(contour_glob);

%%
time_points = {'FU'};
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

sequence_label = {'T1MOLLI', 'T2star', 'T2'};
anatomy_label = {'BloodPool', 'excludeArea', 'freeROI', 'Heart', 'Myocardium', 'MyoReference', 'noReflowArea'};
%name_check = 'Evelyn';
%starting_point = find(strcmp(name_check, Names),1);

save_dir = cat(2, base_dir, '/img/');
data_save_dir = cat(2, base_dir, '/data/');

label_t1 = sequence_label{1};
label_t2 = sequence_label{3};
label_t2star = sequence_label{2};

metrics_save_dir = cat(2, base_dir, '/Results/');
if ~exist(metrics_save_dir, 'dir')
    mkdir(metrics_save_dir);
end

%% Before analysis, parse pre_QualControl
% load(cat(2, metrics_save_dir, 'pre_QualControl.mat'));
%
% status_check = struct;
% for n = 1:length(Names)
%     name = Names{n};
%     status_check(n).Name = name;
%     status_check(n).status = [];
%     status_check(n).status_final = [];
%     tp_count = 1;
%     for tp = 1:length(time_points)
%
%         time_point = time_points{end-tp+1};
%         tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', name, '_', time_point,  '/');
%
%         if ~exist(tp_dir, 'dir')
%             disp(cat(2, 'No folder at: ', name, ' ', time_point));
%         else
%             for i = 1:(length(fieldnames(pre_QualControl(n).status))-1)
%                 slc_loc = cat(2, 'Slice', num2str(i));
%                 if pre_QualControl(n).status(end-tp+1).(slc_loc) == 1
%                     status_check(n).status(tp_count, i) = 1;
%                 elseif pre_QualControl(n).status(end-tp+1).(slc_loc) == 0
%                     status_check(n).status(tp_count, i) = 0;
%                 end
%             end
%
%             tp_count = tp_count + 1;
%         end
%     end
%     for i = 1:(length(fieldnames(pre_QualControl(n).status))-1)
%         if i <= size(status_check(n).status,2)
%             status_check(n).status_final(1,i) = all(status_check(n).status(:,i));
%             % The timepoints of status_check goes from end to beginning
%         end
%     end
% end

%%
% 07/27/2023 modified
se = strel('disk', 1);
metrics = struct;
L_cell = {};
count = 1;

min_roi_value = 1000;
max_roi_value = 1000;

vec=@(x) x(:);
for n = 1:length(Names)
    %for n = 1:2
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
            t1 = vol_img_3D;
            load(myo_glob{1});
            load(roi_glob{1});
            load(remote_glob{1});

            load(cat(2, tp_dir, label_t1, '/', label_t1, '_SliceLoc.mat'));
            [slc_array_t1, idx_reordered] = sort(slc_array);

            clear myo_t1_eroded roi_in_myo_t1 remote_in_myo_t1 roi_t1 remote_t1 t1 myo_t1
            myo_t1_eroded = {};
            for i = 1:length(idx_reordered)
                myo_t1_eroded{i} = imerode(mask_myocardium_3D{i}, se);
            end

            for slc = 1:length(mask_myocardium_3D)
                slc_reordered = idx_reordered(slc);
                roi_in_myo_t1{slc} = myo_t1_eroded{slc_reordered} .* freeROIMask_3D{slc_reordered};
                remote_in_myo_t1{slc} = myo_t1_eroded{slc_reordered} .* myoRefMask_3D{slc_reordered};
                roi_t1{slc} = roi_in_myo_t1{slc} .* vol_img_3D{slc_reordered};
                remote_t1{slc} = remote_in_myo_t1{slc} .* vol_img_3D{slc_reordered};
                t1{slc} = vol_img_3D{slc_reordered};
                % myo_t1 = mask_myocardium_3D;
                myo_t1{slc} = myo_t1_eroded{slc_reordered};
            end

            % T2
            tp_count = tp_count+1;
            myo_glob = glob(cat(2, tp_dir, label_t2, '/', anatomy_label{5}, '/*'));
            roi_glob = glob(cat(2, tp_dir, label_t2, '/',anatomy_label{3}, '/*'));
            remote_glob = glob(cat(2, tp_dir, label_t2, '/',anatomy_label{6}, '/*'));

            load(cat(2, tp_dir, label_t2, '/', label_t2, '_vol_img_3D.mat'));
            t2 = vol_img_3D;
            load(myo_glob{1});
            load(roi_glob{1});
            load(remote_glob{1});

            load(cat(2, tp_dir, label_t2, '/', label_t2, '_SliceLoc.mat'));
            slc_array_t2 = slc_array;

            idx_reordered = Func_AlignSliceLoc(slc_array_t1, slc_array_t2);

            clear myo_t2_eroded roi_in_myo_t2 remote_in_myo_t2 roi_t2 remote_t2 t2 myo_t2
            myo_t2_eroded = {};
            for i = 1:length(mask_myocardium_3D)
                myo_t2_eroded{i} = imerode(mask_myocardium_3D{i}, se);
            end

            for slc = 1:length(mask_myocardium_3D)
                slc_reordered = idx_reordered(slc);
                roi_in_myo_t2{slc} = myo_t2_eroded{slc_reordered} .* freeROIMask_3D{slc_reordered};
                remote_in_myo_t2{slc} = myo_t2_eroded{slc_reordered} .* myoRefMask_3D{slc_reordered};
                roi_t2{slc} = roi_in_myo_t2{slc} .* vol_img_3D{slc_reordered};
                remote_t2{slc} = remote_in_myo_t2{slc} .* vol_img_3D{slc_reordered};
                t2{slc} = vol_img_3D{slc_reordered};
                % myo_t1 = mask_myocardium_3D;
                myo_t2{slc} = myo_t2_eroded{slc_reordered};
            end


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

            idx_reordered = Func_AlignSliceLoc(slc_array_t1, slc_array_ff);

            ff_map = cell(1, length(glob_names));
            r2star_map = cell(1, length(glob_names));
            for f = 1:length(ff_map)
                % f_glob = glob(cat(2,  base_dir, '/FF_Data/',  name, '/', time_point, '/', glob_names{f}, '*.mat'));
                f_glob = glob(cat(2,  base_dir, '/FF_Data_Diego/',  name, '/', time_point, '/*', glob_names{f}, '*.mat'));
                load(f_glob{1}, 'fwmc_ff');
                % ff_map{f} = permute(ff, [2,1]);
                ff_map{f} = fwmc_ff;
                load(f_glob{1}, 'fwmc_r2star')
                % r2star_map{f} = permute(R2s, [2,1]);
                r2star_map{f} = fwmc_r2star;
            end

            % % convert ff_map to matrix
            % ff = zeros(size(ff_map{1}.fwmc_ff,1), size(ff_map{1}.fwmc_ff, 2), length(ff_map));
            % for f = 1:length(ff_map)
            %     ff(:,:,f) = ff_map{f}.fwmc_ff;
            % end



            clear myo_ff_eroded myo_r2star_eroded roi_in_myo_ff remote_in_myo_ff roi_ff remote_ff ff myo_ff ...
                roi_in_myo_r2star remote_in_myo_r2star roi_r2star remote_r2star r2star myo_r2star
            for i = 1:length(mask_myocardium_3D)
                myo_ff_eroded{i} = imerode(mask_myocardium_3D{i}, se);
                myo_r2star_eroded{i} = imerode(mask_myocardium_3D{i}, se);
            end


            % R2star Map

            for slc = 1:length(mask_myocardium_3D)
                slc_reordered = idx_reordered(slc);
                roi_in_myo_ff{slc} = myo_ff_eroded{slc_reordered} .* freeROIMask_3D{slc_reordered};
                remote_in_myo_ff{slc} = myo_ff_eroded{slc_reordered} .* myoRefMask_3D{slc_reordered};
                roi_ff{slc} = roi_in_myo_ff{slc} .* ff_map{slc_reordered};
                remote_ff{slc} = remote_in_myo_ff{slc} .* ff_map{slc_reordered};
                ff{slc} = ff_map{slc_reordered};
                % myo_t1 = mask_myocardium_3D;
                myo_ff{slc} = myo_ff_eroded{slc_reordered};

                roi_in_myo_r2star{slc} = myo_r2star_eroded{slc_reordered} .* freeROIMask_3D{slc_reordered};
                remote_in_myo_r2star{slc} = myo_r2star_eroded{slc_reordered} .* myoRefMask_3D{slc_reordered};
                roi_r2star{slc} = roi_in_myo_r2star{slc} .* r2star_map{slc_reordered};
                remote_r2star{slc} = remote_in_myo_r2star{slc} .* r2star_map{slc_reordered};
                r2star{slc} = r2star_map{slc_reordered};
                % myo_t1 = mask_myocardium_3D;
                myo_r2star{slc} = myo_r2star_eroded{slc_reordered};
            end

            % % convert ff_map to matrix
            % r2star = zeros(size(r2star_map{1}.fwmc_r2star,1), size(r2star_map{2}.fwmc_r2star, 2), length(r2star_map));
            % for f = 1:length(r2star_map)
            %     r2star(:,:,f) = r2star_map{f}.fwmc_r2star;
            % end

            tp_dir2 = cat(2, name_save_dir, '/', name, '_', time_point, '/');
            if ~exist(tp_dir2, 'dir')
                mkdir(tp_dir2);
            end

            % status = status_check(n).status(tp_count,:);
            % AHA Segment

            if 0
                % WHY
                % Some issue with status?
                disp(cat(2, 'Skipped: ', name, ' ', time_point))
            else

                % % exclude_idx = find(status == 0); % Some how previous version doesn't work for Tina
                % if (strcmp(name, 'Tina') && strcmp(time_point, '6MO'))
                %     exclude_idx = [];
                % elseif (strcmp(name, 'Ryn') && strcmp(time_point, '8WK'))
                %     exclude_idx = [];
                % elseif (strcmp(name, 'Gobi') && strcmp(time_point, '1YR'))
                %     exclude_idx = [];
                % elseif (strcmp(name, 'Gobi') && strcmp(time_point, '9MO'))
                %     exclude_idx = [];
                % elseif (strcmp(name, 'Gobi') && strcmp(time_point, '7D'))
                %     exclude_idx = [];
                % elseif (strcmp(name, 'Felicity') && strcmp(time_point, '7D'))
                %     exclude_idx = [];
                % end

                exclude_idx = [];
                % t1{exclude_idx} = [];
                % t2{exclude_idx} = [];
                % ff{exclude_idx) = [];
                % r2star{exclude_idx} = [];
                % roi_in_myo_t1{exclude_idx} = [];
                % remote_in_myo_t1{exclude_idx} = [];
                % roi_in_myo_ff{exclude_idx} = [];
                % remote_in_myo_ff{exclude_idx} = [];
                % roi_in_myo_r2star{exclude_idx} = [];
                % remote_in_myo_r2star{exclude_idx} = [];
                % myo_t1{exclude_idx} = [];
                % roi_in_myo_t2{exclude_idx} = [];
                % remote_in_myo_t2{exclude_idx} = [];
                % myo_t2{exclude_idx} = [];

                metrics(n).TimePoints(tp_count).time_point = time_point;
                metrics(n).TimePoints(tp_count).SliceAnalysis = struct;
                metrics(n).TimePoints(tp_count).HeteroAnalysis = struct;

                
                %
                % for slc = 1:size(roi_in_myo_ff,3)
                %     min_temp = min(min((nonzeros(temp_roi(:,:,slc)))));
                %     max_temp = max(max((nonzeros(temp_roi(:,:,slc)))));
                %
                %     if min_temp < min_roi_value
                %         min_roi_value = min_temp;
                %     end
                %
                %     if max_temp > max_roi_value
                %         max_roi_value = max_temp;
                %     end
                % end

                for slc = 1:length(roi_in_myo_ff)

                    ff{slc}(ff{slc} < 0) = 1e-8;
                    ff{slc}(ff{slc} > 100) = 100;
                    r2star{slc}(r2star{slc}<0)=1e-8;
                    r2star{slc}(r2star{slc}>200)=200;

                    temp_roi{slc} = (t1{slc} .* roi_in_myo_t1{slc});
                    temp_remote{slc} = (t1{slc} .* remote_in_myo_t1{slc});
                    temp_myo{slc} = (t1{slc} .* myo_t1{slc});

                    temp_roi_t2{slc} = (t2{slc} .* roi_in_myo_t2{slc});
                    temp_remote_t2{slc} = (t2{slc} .* remote_in_myo_t2{slc});
                    temp_myo_t2{slc} = (t2{slc} .* myo_t2{slc});

                end

                for slc = 1:length(roi_in_myo_ff)
                    
                    temp_norm_roi = (temp_roi{slc} - min(min((nonzeros(temp_roi{slc})))))./ max(max(nonzeros(temp_roi{slc})));
                    temp_norm_remote = (temp_remote{slc} - min(min((nonzeros(temp_remote{slc})))))./ max(max(nonzeros(temp_remote{slc})));

                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_t1_roi = mean(nonzeros(temp_roi{slc} .* roi_in_myo_t1{slc}));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_t1_remote = mean(nonzeros(temp_remote{slc} .* remote_in_myo_t1{slc}));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_t1_roi = std(nonzeros(temp_roi{slc} .* roi_in_myo_t1{slc}));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_t1_remote = std(nonzeros(temp_remote{slc} .* remote_in_myo_t1{slc}));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_t1_roi = skewness(nonzeros(temp_roi{slc} .* roi_in_myo_t1{slc}));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_t1_remote = skewness(nonzeros(temp_remote{slc} .* remote_in_myo_t1{slc}));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_t1_roi = kurtosis(nonzeros(temp_roi{slc} .* roi_in_myo_t1{slc}))-3;
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_t1_remote = kurtosis(nonzeros(temp_remote{slc} .* remote_in_myo_t1{slc}))-3;

                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_t2_roi = mean(nonzeros(temp_roi_t2{slc} .* roi_in_myo_t2{slc}));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_t2_remote = mean(nonzeros(temp_remote_t2{slc} .* remote_in_myo_t2{slc}));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_t2_roi = std(nonzeros(temp_roi_t2{slc} .* roi_in_myo_t2{slc}));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_t2_remote = std(nonzeros(temp_remote_t2{slc} .* remote_in_myo_t2{slc}));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_t2_roi = skewness(nonzeros(temp_roi_t2{slc} .* roi_in_myo_t2{slc}));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_t2_remote = skewness(nonzeros(temp_remote_t2{slc} .* remote_in_myo_t2{slc}));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_t2_roi = kurtosis(nonzeros(temp_roi_t2{slc} .* roi_in_myo_t2{slc}))-3;
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_t2_remote = kurtosis(nonzeros(temp_remote_t2{slc} .* remote_in_myo_t2{slc}))-3;

                    % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_ff_roi = mean(nonzeros(ff{slc} .* roi_in_myo_ff{slc}));
                    % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_ff_remote = mean(nonzeros(ff{slc} .* remote_in_myo_ff{slc}));
                    % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_ff_roi = std(nonzeros(ff{slc} .* roi_in_myo_ff{slc}));
                    % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_ff_remote = std(nonzeros(ff{slc} .* remote_in_myo_ff{slc}));

                    % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_ff_roi = mean(nonzeros(vec(ff{slc} .* roi_in_myo_ff{slc})), 'omitnan');
                    % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_ff_remote = mean(nonzeros(vec(ff{slc} .* remote_in_myo_ff{slc})), 'omitnan');
                    % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_ff_roi = std(nonzeros(vec(ff{slc} .* roi_in_myo_ff{slc})), 'omitnan');
                    % metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_ff_remote = std(nonzeros(vec(ff{slc} .* remote_in_myo_ff{slc})), 'omitnan');


                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_ff_roi = mean(vec(ff{slc} .* roi_in_myo_ff{slc}), 'omitnan');
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_ff_remote = mean(vec(ff{slc} .* remote_in_myo_ff{slc}), 'omitnan');
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_ff_roi = std(vec(ff{slc} .* roi_in_myo_ff{slc}), 'omitnan');
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_ff_remote = std(vec(ff{slc} .* remote_in_myo_ff{slc}), 'omitnan');

                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_ff_roi = skewness(nonzeros(ff{slc} .* roi_in_myo_ff{slc}));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_ff_remote = skewness(nonzeros(ff{slc} .* remote_in_myo_ff{slc}));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_ff_roi = kurtosis(nonzeros(ff{slc} .* roi_in_myo_ff{slc}))-3;
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_ff_remote = kurtosis(nonzeros(ff{slc} .* remote_in_myo_ff{slc}))-3;

                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_r2star_roi = mean(nonzeros(r2star{slc} .* roi_in_myo_r2star{slc}));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).mean_r2star_remote = mean(nonzeros(r2star{slc} .* remote_in_myo_r2star{slc}));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_r2star_roi = std(nonzeros(r2star{slc} .* roi_in_myo_r2star{slc}));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).std_r2star_remote = std(nonzeros(r2star{slc} .* remote_in_myo_r2star{slc}));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_r2star_roi = skewness(nonzeros(r2star{slc} .* roi_in_myo_r2star{slc}));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).skewness_r2star_remote = skewness(nonzeros(r2star{slc} .* remote_in_myo_r2star{slc}));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_r2star_roi = kurtosis(nonzeros(r2star{slc} .* roi_in_myo_r2star{slc}))-3;
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).kurtosis_r2star_remote = kurtosis(nonzeros(r2star{slc} .* remote_in_myo_r2star{slc}))-3;

                    % temp_norm_roi = uint8((temp_roi(:,:,slc) - min(min((nonzeros(temp_myo(:,:,slc))))))./ (max(max(nonzeros(temp_myo(:,:,slc))))-min(min((nonzeros(temp_myo(:,:,slc)))))) * 256);
                    % temp_norm_remote = uint8((temp_remote(:,:,slc) - min(min((nonzeros(temp_myo(:,:,slc))))))./ (max(max(nonzeros(temp_myo(:,:,slc))))-min(min((nonzeros(temp_myo(:,:,slc)))))) * 256);

                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).max_t1_roi = max(max(nonzeros(temp_roi{slc} .* roi_in_myo_t1{slc})));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).min_t1_roi = min(min(nonzeros(temp_roi{slc} .* roi_in_myo_t1{slc})));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).max_t1_remote = max(max(nonzeros(temp_remote{slc} .* remote_in_myo_t1{slc})));
                    metrics(n).TimePoints(tp_count).SliceAnalysis(slc).min_t1_remote = min(min(nonzeros(temp_remote{slc} .* remote_in_myo_t1{slc})));


                    % Hard-coded for T1 mapping

                    bipolar = temp_roi{slc} - mean(nonzeros(temp_remote{slc}));
                    % weighted_map = double(bipolar < 0) .* 2 .* roi_in_myo_t1(:,:,slc) .* bipolar + double(bipolar >= 0) .* roi_in_myo_t1(:,:,slc) .* bipolar;
                    weighted_map = double(bipolar < 0) .* roi_in_myo_t1{slc} .* bipolar + double(bipolar >= 0) .* roi_in_myo_t1{slc} .* bipolar;

                    bipolar_remote = temp_remote{slc} - mean(nonzeros(temp_remote{slc}));
                    % weighted_map_remote = double(bipolar_remote < 0) .* 2 .* remote_in_myo_t1(:,:,slc) .* bipolar_remote + double(bipolar_remote >= 0) .* remote_in_myo_t1(:,:,slc) .* bipolar_remote;
                    weighted_map_remote = double(bipolar_remote < 0) .* remote_in_myo_t1{slc} .* bipolar_remote + double(bipolar_remote >= 0) .* remote_in_myo_t1{slc} .* bipolar_remote;

                    % lb = 2*(400 - mean(nonzeros(temp_remote(:,:,slc))));
                    lb = (400 - mean(nonzeros(temp_remote{slc})));
                    ub = 1800 - mean(nonzeros(temp_remote{slc}));
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

                    figure(); imhist(t1_norm_roi_slc,40);
                    p_roi = imhist(t1_norm_roi_slc,40);
                    nonZeros = find(p_roi);
                    len = length((nonZeros));
                    pNonZeros = zeros(1,len);

                    for i = 1:len
                        if p_roi(nonZeros(i)) > 0
                            pNonZeros(i) = p_roi(nonZeros(i));
                        end
                    end

                    % normalize pNonZeros so that sum(p) is one.
                    %pNonZeros = pNonZeros ./ sum(p_roi);
                    pNonZeros = nonzeros(pNonZeros ./ sum(pNonZeros));
                    E_roi = -sum(pNonZeros.*log2(pNonZeros));
                    text(50, 10, cat(2, 'Entropy = ', num2str(E_roi)));
                    fname = cat(2, 'T1entropy_Histo_ROI_', time_point, '_', num2str(slc), 'unweighted_original.png');
                    saveas(gcf,cat(2, name_save_dir, '/', name, '_', time_point, '/', fname));


                    non_remote_value = uint8((0 - lb) ./ (ub - lb) * 256);
                    t1_norm_remote_slc = temp_norm_remote;
                    t1_norm_remote_slc(t1_norm_remote_slc == non_remote_value) = [];
                    figure(); imhist(t1_norm_remote_slc,40);

                    p_remote = imhist(t1_norm_remote_slc,40);

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

                    fname = cat(2, 'T1entropy_Histo_Remote_', time_point, '_', num2str(slc), 'unweighted_original.png');
                    saveas(gcf,cat(2, name_save_dir, '/', name, '_', time_point, '/', fname));

                    metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).entropy_roi = E_roi;
                    metrics(n).TimePoints(tp_count).HeteroAnalysis(slc).entropy_remote = E_remote;

                    L = bwlabel(roi_in_myo_t1{slc});

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