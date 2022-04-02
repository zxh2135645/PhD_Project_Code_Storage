clear all;
close all;
clc;
current_dir = pwd;
% Patient data configuration for Khalid
%% 
addpath('../function/');
addpath('../AHA16Segment/');
base_dir = uigetdir;

contour_glob = glob(cat(2, base_dir, '/ContourData/*'));
Names = cell(length(contour_glob), 1); 
for i = 1:length(contour_glob)
    strings = strsplit(contour_glob{i},'/');
    Names{i} = strings{end-1};
end

sequence_label = {'LGE', 'T2star'};
anatomy_label = {'BloodPool', 'excludeArea', 'freeROI', 'Heart', 'Myocardium', 'MyoReference', 'noReflowArea'}; 

names_to_rule_out = {};
RuleOutLabel = NameRuleOutFunc(Names, names_to_rule_out);
Names = Names(RuleOutLabel == 0);

name_check = {'484060000001'};
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

time_points = {'BL', 'BL2', 'FU'};

label_lge = sequence_label{1};
label_t2star = sequence_label{2};

metrics_save_dir = cat(2, base_dir, '/Results/');
if ~exist(metrics_save_dir, 'dir')
   mkdir(metrics_save_dir); 
end

% Read excel file
T = readtable(cat(2, base_dir, '/STEMI_with_IMH-1.xlsx'));
id_array = T.AnonymizationID;
hemo_array = T.withOrWithout;
%% May be used for future used
% Before analysis, parse pre_QualControl
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
%         status_check(n).status_final(1,i) = all(status_check(n).status(:,i));
%         % The timepoints of status_check goes from end to beginning
%     end
% end

%% Main Body

% Array initialization
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

mean_ff_array_50chord = [];
mean_r2star_array_50chord = [];
mean_ff_array_100chord = [];
mean_r2star_array_100chord = [];

% For dichotomize hemo+ and hemo-
mean_ff_roi_array_hemo_n = [];
sd_ff_roi_array_hemo_n = [];
mean_ff_remote_array_hemo_n = [];
sd_ff_remote_array_hemo_n = [];

mean_r2star_roi_array_hemo_n = [];
sd_r2star_roi_array_hemo_n = [];
mean_r2star_remote_array_hemo_n = [];
sd_r2star_remote_array_hemo_n = [];

ff_pixel_roi_array_hemo_n = [];
r2star_pixel_roi_array_hemo_n = [];
ff_pixel_remote_array_hemo_n = [];
r2star_pixel_remote_array_hemo_n = [];

mean_ff_array_50chord_hemo_n = [];
mean_r2star_array_50chord_hemo_n = [];
mean_ff_array_100chord_hemo_n = [];
mean_r2star_array_100chord_hemo_n = [];

% Positive
mean_ff_roi_array_hemo_p = [];
sd_ff_roi_array_hemo_p = [];
mean_ff_remote_array_hemo_p = [];
sd_ff_remote_array_hemo_p = [];

mean_r2star_roi_array_hemo_p = [];
sd_r2star_roi_array_hemo_p = [];
mean_r2star_remote_array_hemo_p = [];
sd_r2star_remote_array_hemo_p = [];

ff_pixel_roi_array_hemo_p = [];
r2star_pixel_roi_array_hemo_p = [];
ff_pixel_remote_array_hemo_p = [];
r2star_pixel_remote_array_hemo_p = [];

mean_ff_array_50chord_hemo_p = [];
mean_r2star_array_50chord_hemo_p = [];
mean_ff_array_100chord_hemo_p = [];
mean_r2star_array_100chord_hemo_p = [];

name_label_hemo_n = {};
slice_count_hemo_n = 1;
name_label_hemo_p = {};
slice_count_hemo_p = 1;

for n = 1:length(Names)
%for n = 7:7
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
    
    name_for_table_searching = insertAfter(name, 6, '-');
    row = find(contains(id_array,name_for_table_searching));
    
    IMH_cell = table2cell(T(row, 13)); % IMH
    IMH = IMH_cell{1};
    %for tp = 1:length(time_points)
    for tp = 2:3
        time_point = time_points{end-tp+1};
        tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', time_point,  '/');
        if ~exist(tp_dir, 'dir')
            disp(cat(2, 'No folder at: ', name, ' ', time_point));
        else
            % T1
            % tp_count = tp_count+1;
            
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
            
            roi_in_myo_ff = mask_myocardium_3D .* freeROIMask_3D;
            remote_in_myo_ff = mask_myocardium_3D .* myoRefMask_3D;
            roi_ff = roi_in_myo_ff .* ff;
            remote_ff = remote_in_myo_ff .* ff;
            myo_ff = mask_myocardium_3D;
            
            % No need to reorder for analysis
%             addpath('../T1NFF/');
%             idx_reordered = Func_AlignSliceLoc(slc_array_t1, slc_array_ff);
%             ff = ff(:,:,idx_reordered);
%             myo_ff = myo_ff(:,:,idx_reordered);
%             remote_ff = remote_ff(:,:,idx_reordered);
%             roi_ff = roi_ff(:,:,idx_reordered);
%             remote_in_myo_ff = remote_in_myo_ff(:,:,idx_reordered);
%             roi_in_myo_ff = roi_in_myo_ff(:,:,idx_reordered);
            
            % R2star Map
            r2star_map = cell(1, length(glob_names));
            for f = 1:length(r2star_map)
                r2star_glob = glob(cat(2,  base_dir, '/FF_Data/',  name, '/', time_point, '/*_', num2str(num_array(f)), '.mat'));
                r2star_map{f} = load(r2star_glob{1}, 'fwmc_r2star');
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
            
%             r2star = r2star(:,:,idx_reordered);
%             myo_r2star = myo_r2star(:,:,idx_reordered);
%             remote_r2star = remote_r2star(:,:,idx_reordered);
%             roi_r2star = roi_r2star(:,:,idx_reordered);
%             remote_in_myo_r2star = remote_in_myo_r2star(:,:,idx_reordered);
%             roi_in_myo_r2star = roi_in_myo_r2star(:,:,idx_reordered);
            
            
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
            roi_in_myo_r2star(roi_in_myo_r2star == 0) = nan;
            roi_in_myo_ff(roi_in_myo_ff == 0) = nan;
            remote_in_myo_r2star(remote_in_myo_r2star == 0) = nan;
            remote_in_myo_ff(roi_in_myo_ff == 0) = nan;
            
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
            
            for slc = 1:size(r2star_roi_masked, 3)
                mean_r2star_roi_array = [mean_r2star_roi_array, nanmean(vec(r2star_roi_masked_px(:,:,slc)))];
                sd_r2star_roi_array = [sd_r2star_roi_array, nanstd(vec(r2star_roi_masked_px(:,:,slc)))];
                mean_ff_roi_array = [mean_ff_roi_array, nanmean(vec(ff_roi_masked_px(:,:,slc)))];
                sd_ff_roi_array = [sd_ff_roi_array, nanstd(vec(ff_roi_masked_px(:,:,slc)))];
                mean_r2star_remote_array = [mean_r2star_remote_array, nanmean(vec(r2star_remote_masked(:,:,slc)))];
                sd_r2star_remote_array = [sd_r2star_remote_array, nanstd(vec(r2star_remote_masked(:,:,slc)))];
                mean_ff_remote_array = [mean_ff_remote_array, nanmean(vec(ff_remote_masked(:,:,slc)))];
                sd_ff_remote_array = [sd_ff_remote_array, nanstd(vec(ff_remote_masked(:,:,slc)))];
                
                name_label{slice_count} = cat(2, name, '_', time_point, '_Slice', num2str(slc));
                slice_count = slice_count + 1;
            end
            
            % Pixel-wise            
            r2star_pixel_roi_array = [r2star_pixel_roi_array; vec(r2star_roi_masked_px(~isnan(r2star_roi_masked_px)))];
            ff_pixel_roi_array = [ff_pixel_roi_array; vec(ff_roi_masked_px(~isnan(ff_roi_masked_px)))];
            r2star_pixel_remote_array = [r2star_pixel_remote_array; vec(r2star_remote_masked(~isnan(r2star_remote_masked)))];
            ff_pixel_remote_array = [ff_pixel_remote_array; vec(ff_remote_masked(~isnan(ff_remote_masked)))];
            
            % AHA 50 chords
            % AHA 100 chords
            % BW_skel = zeros(size(roi_in_myo_ff));
            for slc = 1:size(r2star_roi_masked, 3)
                epi = myo_ff(:,:,slc) - center_mask_ff(:,:,slc) > 0;
                endo = center_mask_ff(:,:,slc) + myo_ff(:,:,slc) > 1;
                
                ff_single_slc = ff(:,:,slc);
                r2star_single_slc = r2star(:,:,slc);
                ff_single_slc(ff_single_slc == 0) = nan;
                ff_single_slc(ff_single_slc == 100) = nan;
                r2star_single_slc(ff_single_slc == 0) = nan;
                r2star_single_slc(ff_single_slc == 100) = nan;
                
                [Segmentpix, stats, Mask_Segn] = AHASegmentation(ff_single_slc, myo_ff, Segn, Groove);
                [Segmentpix, stats, Mask_Segn_epi] = AHASegmentation(ff_single_slc, epi, Segn, Groove);
                [Segmentpix, stats, Mask_Segn_endo] = AHASegmentation(ff_single_slc, endo, Segn, Groove);
                
                for j = 1:Segn
                    seg_mask_ff = Mask_Segn .* roi_in_myo_ff(:,:,slc) == j;
                    seg_mask_ff_epi = Mask_Segn_epi .* roi_in_myo_ff(:,:,slc) == j;
                    seg_mask_ff_endo = Mask_Segn_endo .* roi_in_myo_ff(:,:,slc) == j;
                    
                    if any(seg_mask_ff(:))
                        mean_ff_array_50chord = [mean_ff_array_50chord, nanmean(nonzeros(vec(ff_single_slc .* seg_mask_ff)))];
                        mean_r2star_array_50chord = [mean_r2star_array_50chord, nanmean(nonzeros(vec(r2star_single_slc .* seg_mask_ff)))];
                    end
                    
                    if any(seg_mask_ff_epi(:))
                        mean_ff_array_100chord = [mean_ff_array_100chord, nanmean(nonzeros(vec(ff_single_slc .* seg_mask_ff_epi)))];
                        mean_r2star_array_100chord = [mean_r2star_array_100chord, nanmean(nonzeros(vec(r2star_single_slc .* seg_mask_ff_epi)))];
                    end
                    if any(seg_mask_ff_endo(:))
                        mean_ff_array_100chord = [mean_ff_array_100chord, nanmean(nonzeros(vec(ff_single_slc .* seg_mask_ff_endo)))];
                        mean_r2star_array_100chord = [mean_r2star_array_100chord, nanmean(nonzeros(vec(r2star_single_slc .* seg_mask_ff_endo)))];
                    end
                    
                    if strcmp(IMH, '-')
                        if any(seg_mask_ff(:))
                            mean_ff_array_50chord_hemo_n = [mean_ff_array_50chord_hemo_n, nanmean(nonzeros(vec(ff_single_slc .* seg_mask_ff)))];
                            mean_r2star_array_50chord_hemo_n = [mean_r2star_array_50chord_hemo_n, nanmean(nonzeros(vec(r2star_single_slc .* seg_mask_ff)))];
                        end
                        
                        if any(seg_mask_ff_epi(:))
                            mean_ff_array_100chord_hemo_n = [mean_ff_array_100chord_hemo_n, nanmean(nonzeros(vec(ff_single_slc .* seg_mask_ff_epi)))];
                            mean_r2star_array_100chord_hemo_n = [mean_r2star_array_100chord_hemo_n, nanmean(nonzeros(vec(r2star_single_slc .* seg_mask_ff_epi)))];
                        end
                        if any(seg_mask_ff_endo(:))
                            mean_ff_array_100chord_hemo_n = [mean_ff_array_100chord_hemo_n, nanmean(nonzeros(vec(ff_single_slc .* seg_mask_ff_endo)))];
                            mean_r2star_array_100chord_hemo_n = [mean_r2star_array_100chord_hemo_n, nanmean(nonzeros(vec(r2star_single_slc .* seg_mask_ff_endo)))];
                        end
                        
                    elseif strcmp(IMH, '+')
                        if any(seg_mask_ff(:))
                            mean_ff_array_50chord_hemo_p = [mean_ff_array_50chord_hemo_p, nanmean(nonzeros(vec(ff_single_slc .* seg_mask_ff)))];
                            mean_r2star_array_50chord_hemo_p = [mean_r2star_array_50chord_hemo_p, nanmean(nonzeros(vec(r2star_single_slc .* seg_mask_ff)))];
                        end
                        
                        if any(seg_mask_ff_epi(:))
                            mean_ff_array_100chord_hemo_p = [mean_ff_array_100chord_hemo_p, nanmean(nonzeros(vec(ff_single_slc .* seg_mask_ff_epi)))];
                            mean_r2star_array_100chord_hemo_p = [mean_r2star_array_100chord_hemo_p, nanmean(nonzeros(vec(r2star_single_slc .* seg_mask_ff_epi)))];
                        end
                        if any(seg_mask_ff_endo(:))
                            mean_ff_array_100chord_hemo_p = [mean_ff_array_100chord_hemo_p, nanmean(nonzeros(vec(ff_single_slc .* seg_mask_ff_endo)))];
                            mean_r2star_array_100chord_hemo_p = [mean_r2star_array_100chord_hemo_p, nanmean(nonzeros(vec(r2star_single_slc .* seg_mask_ff_endo)))];
                        end
                    end
                end
            end
            
            % Dichotomize into hemo+ and hemo-
            if strcmp(IMH, '-')
                % ROI
                for slc = 1:size(r2star_roi_masked, 3)
                    mean_r2star_roi_array_hemo_n = [mean_r2star_roi_array_hemo_n, nanmean(vec(r2star_roi_masked_px(:,:,slc)))];
                    sd_r2star_roi_array_hemo_n = [sd_r2star_roi_array_hemo_n, nanstd(vec(r2star_roi_masked_px(:,:,slc)))];
                    mean_ff_roi_array_hemo_n = [mean_ff_roi_array_hemo_n, nanmean(vec(ff_roi_masked_px(:,:,slc)))];
                    sd_ff_roi_array_hemo_n = [sd_ff_roi_array_hemo_n, nanstd(vec(ff_roi_masked_px(:,:,slc)))];
                    mean_r2star_remote_array_hemo_n = [mean_r2star_remote_array_hemo_n, nanmean(vec(r2star_remote_masked(:,:,slc)))];
                    sd_r2star_remote_array_hemo_n = [sd_r2star_remote_array_hemo_n, nanstd(vec(r2star_remote_masked(:,:,slc)))];
                    mean_ff_remote_array_hemo_n = [mean_ff_remote_array_hemo_n, nanmean(vec(ff_remote_masked(:,:,slc)))];
                    sd_ff_remote_array_hemo_n = [sd_ff_remote_array_hemo_n, nanstd(vec(ff_remote_masked(:,:,slc)))];
                    
                    name_label_hemo_n{slice_count_hemo_n} = cat(2, name, '_', time_point, '_Slice', num2str(slc));
                    slice_count_hemo_n = slice_count_hemo_n + 1;
                end
                
                % Pixel-wise
                r2star_pixel_roi_array_hemo_n = [r2star_pixel_roi_array_hemo_n; vec(r2star_roi_masked_px(~isnan(r2star_roi_masked_px)))];
                ff_pixel_roi_array_hemo_n = [ff_pixel_roi_array_hemo_n; vec(ff_roi_masked_px(~isnan(ff_roi_masked_px)))];
                r2star_pixel_remote_array_hemo_n = [r2star_pixel_remote_array_hemo_n; vec(r2star_remote_masked(~isnan(r2star_remote_masked)))];
                ff_pixel_remote_array_hemo_n = [ff_pixel_remote_array_hemo_n; vec(ff_remote_masked(~isnan(ff_remote_masked)))];
                
            elseif strcmp(IMH, '+')
                % ROI
                for slc = 1:size(r2star_roi_masked, 3)
                    mean_r2star_roi_array_hemo_p = [mean_r2star_roi_array_hemo_p, nanmean(vec(r2star_roi_masked_px(:,:,slc)))];
                    sd_r2star_roi_array_hemo_p = [sd_r2star_roi_array_hemo_p, nanstd(vec(r2star_roi_masked_px(:,:,slc)))];
                    mean_ff_roi_array_hemo_p = [mean_ff_roi_array_hemo_p, nanmean(vec(ff_roi_masked_px(:,:,slc)))];
                    sd_ff_roi_array_hemo_p = [sd_ff_roi_array_hemo_p, nanstd(vec(ff_roi_masked_px(:,:,slc)))];
                    mean_r2star_remote_array_hemo_p = [mean_r2star_remote_array_hemo_p, nanmean(vec(r2star_remote_masked(:,:,slc)))];
                    sd_r2star_remote_array_hemo_p = [sd_r2star_remote_array_hemo_p, nanstd(vec(r2star_remote_masked(:,:,slc)))];
                    mean_ff_remote_array_hemo_p = [mean_ff_remote_array_hemo_p, nanmean(vec(ff_remote_masked(:,:,slc)))];
                    sd_ff_remote_array_hemo_p = [sd_ff_remote_array_hemo_p, nanstd(vec(ff_remote_masked(:,:,slc)))];
                    
                    name_label_hemo_p{slice_count_hemo_p} = cat(2, name, '_', time_point, '_Slice', num2str(slc));
                    slice_count_hemo_p = slice_count_hemo_p + 1;
                end
                
                % Pixel-wise
                r2star_pixel_roi_array_hemo_p = [r2star_pixel_roi_array_hemo_p; vec(r2star_roi_masked_px(~isnan(r2star_roi_masked_px)))];
                ff_pixel_roi_array_hemo_p = [ff_pixel_roi_array_hemo_p; vec(ff_roi_masked_px(~isnan(ff_roi_masked_px)))];
                r2star_pixel_remote_array_hemo_p = [r2star_pixel_remote_array_hemo_p; vec(r2star_remote_masked(~isnan(r2star_remote_masked)))];
                ff_pixel_remote_array_hemo_p = [ff_pixel_remote_array_hemo_p; vec(ff_remote_masked(~isnan(ff_remote_masked)))];
            end
        end
        close all;
    end
end

%% Plot
color_cell_roi = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell_remote = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};

figure(); 
subplot(2,2,1);
% plot(mean_ff_roi_array, mean_r2star_roi_array, 'o')
mdl = fitlm(mean_ff_roi_array, mean_r2star_roi_array);
scatter(mean_ff_roi_array, mean_r2star_roi_array, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
Y = mean_ff_roi_array .* mdl.Coefficients.Estimate(2) + mdl.Coefficients.Estimate(1);
hold on;
plot(mean_ff_roi_array, Y, 'k', 'LineWidth', 1);
xlabel('FF (%)');
ylabel('R2star (s^{-1})');
%xlim([0 50]); ylim([0 150]);
yl = ylim;
xl = xlim;
text(0.3*xl(2), yl(1)+10, cat(2,'Y = ', num2str(mdl.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl.Rsquared.Ordinary,3)), 'FontSize', 12)
title('ROI');

subplot(2,2,2); % plot(ff_pixel_roi_array, r2star_pixel_roi_array, 'o')
mdl_px = fitlm(ff_pixel_roi_array, r2star_pixel_roi_array);
scatter(ff_pixel_roi_array, r2star_pixel_roi_array, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
Y = ff_pixel_roi_array .* mdl_px.Coefficients.Estimate(2) + mdl_px.Coefficients.Estimate(1);
hold on;
plot(ff_pixel_roi_array, Y, 'k', 'LineWidth', 1);
xlabel('FF (%)');
ylabel('R2star (s^{-1})');
%xlim([0 50]); ylim([0 150]);
yl = ylim;
xl = xlim;
text(0.3*xl(2), yl(1)+10, cat(2,'Y = ', num2str(mdl_px.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_px.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_px.Rsquared.Ordinary,3)), 'FontSize', 12)
title('Pixel-wise');

subplot(2,2,3);
%plot(mean_ff_array_50chord, mean_r2star_array_50chord, 'o');
mdl_50chord = fitlm(mean_ff_array_50chord, mean_r2star_array_50chord);
scatter(mean_ff_array_50chord, mean_r2star_array_50chord, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
Y = mean_ff_array_50chord .* mdl_50chord.Coefficients.Estimate(2) + mdl_50chord.Coefficients.Estimate(1);
hold on;
plot(mean_ff_array_50chord, Y, 'k', 'LineWidth', 1);
xlabel('FF (%)');
ylabel('R2star (s^{-1})');
%xlim([0 50]); ylim([0 150]);
yl = ylim;
xl = xlim;
text(0.3*xl(2), yl(1)+10, cat(2,'Y = ', num2str(mdl_50chord.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_50chord.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_50chord.Rsquared.Ordinary,3)), 'FontSize', 12)
title('50 chords');

subplot(2,2,4);
% plot(mean_ff_array_100chord, mean_r2star_array_100chord, 'o');
mdl_100chord = fitlm(mean_ff_array_100chord, mean_r2star_array_100chord);
scatter(mean_ff_array_100chord, mean_r2star_array_100chord, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
Y = mean_ff_array_100chord .* mdl_100chord.Coefficients.Estimate(2) + mdl_100chord.Coefficients.Estimate(1);
hold on;
plot(mean_ff_array_100chord, Y, 'k', 'LineWidth', 1);
xlabel('FF (%)');
ylabel('R2star (s^{-1})');
%xlim([0 50]); ylim([0 150]);
yl = ylim;
xl = xlim;
text(0.3*xl(2), yl(1)+10, cat(2,'Y = ', num2str(mdl_100chord.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_100chord.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_100chord.Rsquared.Ordinary,3)), 'FontSize', 12)
title('100 chords');

%% Plot hemo+ and hemo-
%% hemo-
color_cell_roi = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell_remote = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};

figure(); 
subplot(2,2,1);
% plot(mean_ff_roi_array, mean_r2star_roi_array, 'o')
mdl_hemo_n = fitlm(mean_ff_roi_array_hemo_n, mean_r2star_roi_array_hemo_n);
scatter(mean_ff_roi_array_hemo_n, mean_r2star_roi_array_hemo_n, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
Y = mean_ff_roi_array_hemo_n .* mdl_hemo_n.Coefficients.Estimate(2) + mdl_hemo_n.Coefficients.Estimate(1);
hold on;
plot(mean_ff_roi_array_hemo_n, Y, 'k', 'LineWidth', 1);
xlabel('FF (%)');
ylabel('R2star (s^{-1})');
xlim([0 50]); ylim([0 150]);
yl = ylim;
xl = xlim;
text(0.3*xl(2), yl(1)+10, cat(2,'Y = ', num2str(mdl_hemo_n.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_hemo_n.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_hemo_n.Rsquared.Ordinary,3)), 'FontSize', 12)
title('ROI');

subplot(2,2,2); % plot(ff_pixel_roi_array, r2star_pixel_roi_array, 'o')
mdl_px_hemo_n = fitlm(ff_pixel_roi_array_hemo_n, r2star_pixel_roi_array_hemo_n);
scatter(ff_pixel_roi_array_hemo_n, r2star_pixel_roi_array_hemo_n, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
Y = ff_pixel_roi_array_hemo_n .* mdl_px_hemo_n.Coefficients.Estimate(2) + mdl_px_hemo_n.Coefficients.Estimate(1);
hold on;
plot(ff_pixel_roi_array_hemo_n, Y, 'k', 'LineWidth', 1);
xlabel('FF (%)');
ylabel('R2star (s^{-1})');
xlim([0 50]); ylim([0 150]);
yl = ylim;
xl = xlim;
text(0.3*xl(2), yl(1)+10, cat(2,'Y = ', num2str(mdl_px_hemo_n.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_px_hemo_n.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_px_hemo_n.Rsquared.Ordinary,3)), 'FontSize', 12)
title('Pixel-wise');

subplot(2,2,3);
%plot(mean_ff_array_50chord, mean_r2star_array_50chord, 'o');
mdl_50chord_hemo_n = fitlm(mean_ff_array_50chord_hemo_n, mean_r2star_array_50chord_hemo_n);
scatter(mean_ff_array_50chord_hemo_n, mean_r2star_array_50chord_hemo_n, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
Y = mean_ff_array_50chord_hemo_n .* mdl_50chord_hemo_n.Coefficients.Estimate(2) + mdl_50chord_hemo_n.Coefficients.Estimate(1);
hold on;
plot(mean_ff_array_50chord_hemo_n, Y, 'k', 'LineWidth', 1);
xlabel('FF (%)');
ylabel('R2star (s^{-1})');
xlim([0 50]); ylim([0 150]);
yl = ylim;
xl = xlim;
text(0.3*xl(2), yl(1)+10, cat(2,'Y = ', num2str(mdl_50chord_hemo_n.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_50chord_hemo_n.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_50chord_hemo_n.Rsquared.Ordinary,3)), 'FontSize', 12)
title('50 chords');

subplot(2,2,4);
% plot(mean_ff_array_100chord, mean_r2star_array_100chord, 'o');
mdl_100chord_hemo_n = fitlm(mean_ff_array_100chord_hemo_n, mean_r2star_array_100chord_hemo_n);
scatter(mean_ff_array_100chord_hemo_n, mean_r2star_array_100chord_hemo_n, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
Y = mean_ff_array_100chord_hemo_n .* mdl_100chord_hemo_n.Coefficients.Estimate(2) + mdl_100chord_hemo_n.Coefficients.Estimate(1);
hold on;
plot(mean_ff_array_100chord_hemo_n, Y, 'k', 'LineWidth', 1);
xlabel('FF (%)');
ylabel('R2star (s^{-1})');
xlim([0 50]); ylim([0 150]);
yl = ylim;
xl = xlim;
text(0.3*xl(2), yl(1)+10, cat(2,'Y = ', num2str(mdl_100chord_hemo_n.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_100chord_hemo_n.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_100chord_hemo_n.Rsquared.Ordinary,3)), 'FontSize', 12)
title('100 chords');

%% hemo+
color_cell_roi = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell_remote = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};

figure(); 
subplot(2,2,1);
% plot(mean_ff_roi_array, mean_r2star_roi_array, 'o')
mdl_hemo_p = fitlm(mean_ff_roi_array_hemo_p, mean_r2star_roi_array_hemo_p);
scatter(mean_ff_roi_array_hemo_p, mean_r2star_roi_array_hemo_p, 64, 'MarkerEdgeColor', color_cell_remote{5}, 'MarkerFaceColor', color_cell_remote{3});
Y = mean_ff_roi_array_hemo_p .* mdl_hemo_p.Coefficients.Estimate(2) + mdl_hemo_p.Coefficients.Estimate(1);
hold on;
plot(mean_ff_roi_array_hemo_p, Y, 'k', 'LineWidth', 1);
xlabel('FF (%)');
ylabel('R2star (s^{-1})');
xlim([0 50]); ylim([0 150]);
yl = ylim;
xl = xlim;
text(0.3*xl(2), yl(1)+10, cat(2,'Y = ', num2str(mdl_hemo_p.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_hemo_p.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_hemo_p.Rsquared.Ordinary,3)), 'FontSize', 12)
title('ROI');

subplot(2,2,2); % plot(ff_pixel_roi_array, r2star_pixel_roi_array, 'o')
mdl_px_hemo_p = fitlm(ff_pixel_roi_array_hemo_p, r2star_pixel_roi_array_hemo_p);
scatter(ff_pixel_roi_array_hemo_p, r2star_pixel_roi_array_hemo_p, 64, 'MarkerEdgeColor', color_cell_remote{5}, 'MarkerFaceColor', color_cell_remote{3});
Y = ff_pixel_roi_array_hemo_p .* mdl_px_hemo_p.Coefficients.Estimate(2) + mdl_px_hemo_p.Coefficients.Estimate(1);
hold on;
plot(ff_pixel_roi_array_hemo_p, Y, 'k', 'LineWidth', 1);
xlabel('FF (%)');
ylabel('R2star (s^{-1})');
xlim([0 50]); ylim([0 150]);
yl = ylim;
xl = xlim;
text(0.3*xl(2), yl(1)+10, cat(2,'Y = ', num2str(mdl_px_hemo_p.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_px_hemo_p.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_px_hemo_p.Rsquared.Ordinary,3)), 'FontSize', 12)
title('Pixel-wise');

subplot(2,2,3);
%plot(mean_ff_array_50chord, mean_r2star_array_50chord, 'o');
mdl_50chord_hemo_p = fitlm(mean_ff_array_50chord_hemo_p, mean_r2star_array_50chord_hemo_p);
scatter(mean_ff_array_50chord_hemo_p, mean_r2star_array_50chord_hemo_p, 64, 'MarkerEdgeColor', color_cell_remote{5}, 'MarkerFaceColor', color_cell_remote{3});
Y = mean_ff_array_50chord_hemo_p .* mdl_50chord_hemo_p.Coefficients.Estimate(2) + mdl_50chord_hemo_p.Coefficients.Estimate(1);
hold on;
plot(mean_ff_array_50chord_hemo_p, Y, 'k', 'LineWidth', 1);
xlabel('FF (%)');
ylabel('R2star (s^{-1})');
xlim([0 50]); ylim([0 150]);
yl = ylim;
xl = xlim;
text(0.3*xl(2), yl(1)+10, cat(2,'Y = ', num2str(mdl_50chord_hemo_p.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_50chord_hemo_p.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_50chord_hemo_p.Rsquared.Ordinary,3)), 'FontSize', 12)
title('50 chords');

subplot(2,2,4);
% plot(mean_ff_array_100chord, mean_r2star_array_100chord, 'o');
mdl_100chord_hemo_p = fitlm(mean_ff_array_100chord_hemo_p, mean_r2star_array_100chord_hemo_p);
scatter(mean_ff_array_100chord_hemo_p, mean_r2star_array_100chord_hemo_p, 64, 'MarkerEdgeColor', color_cell_remote{5}, 'MarkerFaceColor', color_cell_remote{3});
Y = mean_ff_array_100chord_hemo_p .* mdl_100chord_hemo_p.Coefficients.Estimate(2) + mdl_100chord_hemo_p.Coefficients.Estimate(1);
hold on;
plot(mean_ff_array_100chord_hemo_p, Y, 'k', 'LineWidth', 1);
xlabel('FF (%)');
ylabel('R2star (s^{-1})');
xlim([0 50]); ylim([0 150]);
yl = ylim;
xl = xlim;
text(0.3*xl(2), yl(1)+10, cat(2,'Y = ', num2str(mdl_100chord_hemo_p.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_100chord_hemo_p.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_100chord_hemo_p.Rsquared.Ordinary,3)), 'FontSize', 12)
title('100 chords');

%% Save as struct (FU)
metrics_to_save = struct;
metrics_to_save.name_label = name_label;
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

metrics_to_save.Pixel = struct;
metrics_to_save.Pixel.ff_pixel_roi_array = ff_pixel_roi_array;
metrics_to_save.Pixel.r2star_pixel_roi_array = r2star_pixel_roi_array;
metrics_to_save.Pixel.ff_pixel_remote_array = ff_pixel_remote_array;
metrics_to_save.Pixel.r2star_pixel_remote_array = r2star_pixel_remote_array;

metrics_to_save.Chord50 = struct;
metrics_to_save.Chord50.mean_ff_array_50chord = mean_ff_array_50chord;
metrics_to_save.Chord50.mean_r2star_array_50chord = mean_r2star_array_50chord;

metrics_to_save.Chord100 = struct;
metrics_to_save.Chord100.mean_ff_array_100chord = mean_ff_array_100chord;
metrics_to_save.Chord100.mean_r2star_array_100chord = mean_r2star_array_100chord;

% For dichotomize hemo+ and hemo-
metrics_to_save.ROI.mean_ff_roi_array_hemo_n = mean_ff_roi_array_hemo_n;
metrics_to_save.ROI.sd_ff_roi_array_hemo_n = sd_ff_roi_array_hemo_n;
metrics_to_save.ROI.mean_ff_remote_array_hemo_n = mean_ff_remote_array_hemo_n;
metrics_to_save.ROI.sd_ff_remote_array_hemo_n = sd_ff_remote_array_hemo_n;

metrics_to_save.ROI.mean_r2star_roi_array_hemo_n = mean_r2star_roi_array_hemo_n;
metrics_to_save.ROI.sd_r2star_roi_array_hemo_n = sd_r2star_roi_array_hemo_n;
metrics_to_save.ROI.mean_r2star_remote_array_hemo_n = mean_r2star_remote_array_hemo_n;
metrics_to_save.ROI.sd_r2star_remote_array_hemo_n = sd_r2star_remote_array_hemo_n;

metrics_to_save.Pixel.ff_pixel_roi_array_hemo_n = ff_pixel_roi_array_hemo_n;
metrics_to_save.Pixel.r2star_pixel_roi_array_hemo_n = r2star_pixel_roi_array_hemo_n;
metrics_to_save.Pixel.ff_pixel_remote_array_hemo_n = ff_pixel_remote_array_hemo_n;
metrics_to_save.Pixel.r2star_pixel_remote_array_hemo_n = r2star_pixel_remote_array_hemo_n;

metrics_to_save.Chord50.mean_ff_array_50chord_hemo_n = mean_ff_array_50chord_hemo_n;
metrics_to_save.Chord50.mean_r2star_array_50chord_hemo_n = mean_r2star_array_50chord_hemo_n;
metrics_to_save.Chord100.mean_ff_array_100chord_hemo_n = mean_ff_array_100chord_hemo_n;
metrics_to_save.Chord100.mean_r2star_array_100chord_hemo_n = mean_r2star_array_100chord_hemo_n;

% Positive
metrics_to_save.ROI.mean_ff_roi_array_hemo_p = mean_ff_roi_array_hemo_p;
metrics_to_save.ROI.sd_ff_roi_array_hemo_p = sd_ff_roi_array_hemo_p;
metrics_to_save.ROI.mean_ff_remote_array_hemo_p = mean_ff_remote_array_hemo_p;
metrics_to_save.ROI.sd_ff_remote_array_hemo_p = sd_ff_remote_array_hemo_p;

metrics_to_save.ROI.mean_r2star_roi_array_hemo_p = mean_r2star_roi_array_hemo_p;
metrics_to_save.ROI.sd_r2star_roi_array_hemo_p = sd_r2star_roi_array_hemo_p;
metrics_to_save.ROI.mean_r2star_remote_array_hemo_p = mean_r2star_remote_array_hemo_p;
metrics_to_save.ROI.sd_r2star_remote_array_hemo_p = sd_r2star_remote_array_hemo_p;

metrics_to_save.Pixel.ff_pixel_roi_array_hemo_p = ff_pixel_roi_array_hemo_p;
metrics_to_save.Pixel.r2star_pixel_roi_array_hemo_p = r2star_pixel_roi_array_hemo_p;
metrics_to_save.Pixel.ff_pixel_remote_array_hemo_p = ff_pixel_remote_array_hemo_p;
metrics_to_save.Pixel.r2star_pixel_remote_array_hemo_p = r2star_pixel_remote_array_hemo_p;

metrics_to_save.Chord50.mean_ff_array_50chord_hemo_p = mean_ff_array_50chord_hemo_p;
metrics_to_save.Chord50.mean_r2star_array_50chord_hemo_p = mean_r2star_array_50chord_hemo_p;
metrics_to_save.Chord100.mean_ff_array_100chord_hemo_p = mean_ff_array_100chord_hemo_p;
metrics_to_save.Chord100.mean_r2star_array_100chord_hemo_p = mean_r2star_array_100chord_hemo_p;

metrics_to_save.ROI.mdl = mdl;
metrics_to_save.ROI.mdl_hemo_n = mdl_hemo_n;
metrics_to_save.ROI.mdl_hemo_p = mdl_hemo_p;
metrics_to_save.Pixel.mdl_px = mdl_px;
metrics_to_save.Pixel.mdl_px_hemo_n = mdl_px_hemo_n;
metrics_to_save.Pixel.mdl_px_hemo_p = mdl_px_hemo_p;
metrics_to_save.Chord50.mdl_50chord = mdl_50chord;
metrics_to_save.Chord50.mdl_50chord_hemo_n = mdl_50chord_hemo_n;
metrics_to_save.Chord50.mdl_50chord_hemo_p = mdl_50chord_hemo_p;
metrics_to_save.Chord100.mdl_50chord = mdl_100chord;
metrics_to_save.Chord100.mdl_50chord_hemo_n = mdl_100chord_hemo_n;
metrics_to_save.Chord100.mdl_50chord_hemo_p = mdl_100chord_hemo_p;

save(cat(2, metrics_save_dir, 'Metrics_FatNIron_Analysis.mat'), '-struct', 'metrics_to_save');

%% BL struct
metrics_to_save = struct;
metrics_to_save.name_label_BL = name_label;
metrics_to_save.name_label_hemo_p_BL = name_label_hemo_p;
metrics_to_save.name_label_hemo_n_BL = name_label_hemo_n;

metrics_to_save.ROI_BL = struct;
metrics_to_save.ROI_BL.mean_ff_roi_array = mean_ff_roi_array;
metrics_to_save.ROI_BL.sd_ff_roi_array = sd_ff_roi_array;
metrics_to_save.ROI_BL.mean_ff_remote_array = mean_ff_remote_array;
metrics_to_save.ROI_BL.sd_ff_remote_array = sd_ff_remote_array;

metrics_to_save.ROI_BL.mean_r2star_roi_array = mean_r2star_roi_array;
metrics_to_save.ROI_BL.sd_r2star_roi_array = sd_r2star_roi_array;
metrics_to_save.ROI_BL.mean_r2star_remote_array = mean_r2star_remote_array;
metrics_to_save.ROI_BL.sd_r2star_remote_array = sd_r2star_remote_array;

metrics_to_save.Pixel_BL = struct;
metrics_to_save.Pixel_BL.ff_pixel_roi_array = ff_pixel_roi_array;
metrics_to_save.Pixel_BL.r2star_pixel_roi_array = r2star_pixel_roi_array;
metrics_to_save.Pixel_BL.ff_pixel_remote_array = ff_pixel_remote_array;
metrics_to_save.Pixel_BL.r2star_pixel_remote_array = r2star_pixel_remote_array;

metrics_to_save.Chord50_BL = struct;
metrics_to_save.Chord50_BL.mean_ff_array_50chord = mean_ff_array_50chord;
metrics_to_save.Chord50_BL.mean_r2star_array_50chord = mean_r2star_array_50chord;

metrics_to_save.Chord100_BL = struct;
metrics_to_save.Chord100_BL.mean_ff_array_100chord = mean_ff_array_100chord;
metrics_to_save.Chord100_BL.mean_r2star_array_100chord = mean_r2star_array_100chord;

% For dichotomize hemo+ and hemo-
metrics_to_save.ROI_BL.mean_ff_roi_array_hemo_n = mean_ff_roi_array_hemo_n;
metrics_to_save.ROI_BL.sd_ff_roi_array_hemo_n = sd_ff_roi_array_hemo_n;
metrics_to_save.ROI_BL.mean_ff_remote_array_hemo_n = mean_ff_remote_array_hemo_n;
metrics_to_save.ROI_BL.sd_ff_remote_array_hemo_n = sd_ff_remote_array_hemo_n;

metrics_to_save.ROI_BL.mean_r2star_roi_array_hemo_n = mean_r2star_roi_array_hemo_n;
metrics_to_save.ROI_BL.sd_r2star_roi_array_hemo_n = sd_r2star_roi_array_hemo_n;
metrics_to_save.ROI_BL.mean_r2star_remote_array_hemo_n = mean_r2star_remote_array_hemo_n;
metrics_to_save.ROI_BL.sd_r2star_remote_array_hemo_n = sd_r2star_remote_array_hemo_n;

metrics_to_save.Pixel_BL.ff_pixel_roi_array_hemo_n = ff_pixel_roi_array_hemo_n;
metrics_to_save.Pixel_BL.r2star_pixel_roi_array_hemo_n = r2star_pixel_roi_array_hemo_n;
metrics_to_save.Pixel_BL.ff_pixel_remote_array_hemo_n = ff_pixel_remote_array_hemo_n;
metrics_to_save.Pixel_BL.r2star_pixel_remote_array_hemo_n = r2star_pixel_remote_array_hemo_n;

metrics_to_save.Chord50_BL.mean_ff_array_50chord_hemo_n = mean_ff_array_50chord_hemo_n;
metrics_to_save.Chord50_BL.mean_r2star_array_50chord_hemo_n = mean_r2star_array_50chord_hemo_n;
metrics_to_save.Chord100_BL.mean_ff_array_100chord_hemo_n = mean_ff_array_100chord_hemo_n;
metrics_to_save.Chord100_BL.mean_r2star_array_100chord_hemo_n = mean_r2star_array_100chord_hemo_n;

% Positive
metrics_to_save.ROI_BL.mean_ff_roi_array_hemo_p = mean_ff_roi_array_hemo_p;
metrics_to_save.ROI_BL.sd_ff_roi_array_hemo_p = sd_ff_roi_array_hemo_p;
metrics_to_save.ROI_BL.mean_ff_remote_array_hemo_p = mean_ff_remote_array_hemo_p;
metrics_to_save.ROI_BL.sd_ff_remote_array_hemo_p = sd_ff_remote_array_hemo_p;

metrics_to_save.ROI_BL.mean_r2star_roi_array_hemo_p = mean_r2star_roi_array_hemo_p;
metrics_to_save.ROI_BL.sd_r2star_roi_array_hemo_p = sd_r2star_roi_array_hemo_p;
metrics_to_save.ROI_BL.mean_r2star_remote_array_hemo_p = mean_r2star_remote_array_hemo_p;
metrics_to_save.ROI_BL.sd_r2star_remote_array_hemo_p = sd_r2star_remote_array_hemo_p;

metrics_to_save.Pixel_BL.ff_pixel_roi_array_hemo_p = ff_pixel_roi_array_hemo_p;
metrics_to_save.Pixel_BL.r2star_pixel_roi_array_hemo_p = r2star_pixel_roi_array_hemo_p;
metrics_to_save.Pixel_BL.ff_pixel_remote_array_hemo_p = ff_pixel_remote_array_hemo_p;
metrics_to_save.Pixel_BL.r2star_pixel_remote_array_hemo_p = r2star_pixel_remote_array_hemo_p;

metrics_to_save.Chord50_BL.mean_ff_array_50chord_hemo_p = mean_ff_array_50chord_hemo_p;
metrics_to_save.Chord50_BL.mean_r2star_array_50chord_hemo_p = mean_r2star_array_50chord_hemo_p;
metrics_to_save.Chord100_BL.mean_ff_array_100chord_hemo_p = mean_ff_array_100chord_hemo_p;
metrics_to_save.Chord100_BL.mean_r2star_array_100chord_hemo_p = mean_r2star_array_100chord_hemo_p;

metrics_to_save.ROI_BL.mdl = mdl;
metrics_to_save.ROI_BL.mdl_hemo_n = mdl_hemo_n;
metrics_to_save.ROI_BL.mdl_hemo_p = mdl_hemo_p;
metrics_to_save.Pixel_BL.mdl_px = mdl_px;
metrics_to_save.Pixel_BL.mdl_px_hemo_n = mdl_px_hemo_n;
metrics_to_save.Pixel_BL.mdl_px_hemo_p = mdl_px_hemo_p;
metrics_to_save.Chord50_BL.mdl_50chord = mdl_50chord;
metrics_to_save.Chord50_BL.mdl_50chord_hemo_n = mdl_50chord_hemo_n;
metrics_to_save.Chord50_BL.mdl_50chord_hemo_p = mdl_50chord_hemo_p;
metrics_to_save.Chord100_BL.mdl_50chord = mdl_100chord;
metrics_to_save.Chord100_BL.mdl_50chord_hemo_n = mdl_100chord_hemo_n;
metrics_to_save.Chord100_BL.mdl_50chord_hemo_p = mdl_100chord_hemo_p;

save(cat(2, metrics_save_dir, 'Metrics_FatNIron_Analysis_BL.mat'), '-struct', 'metrics_to_save');

%% Before run the code below, need to load both Metrics_FatNIron_Analysis_BL.mat and Metrics_FatNIron_Analysis.mat
%% BL + FU
color_cell_roi = {[254,240,217]/255, [253,204,138]/255, [252,141,89]/255, [227,74,51]/255, [179,0,0]/255};
color_cell_remote = {[241, 238, 246]/255, [189, 201, 225]/255, [116, 169, 207]/255, [43, 140, 190]/255, [4, 90, 141]/255};
figure();
mdl_50chord = fitlm(Chord50.mean_ff_array_50chord, Chord50.mean_r2star_array_50chord);
scatter(Chord50.mean_ff_array_50chord, Chord50.mean_r2star_array_50chord, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
Y = Chord50.mean_ff_array_50chord .* mdl_50chord.Coefficients.Estimate(2) + mdl_50chord.Coefficients.Estimate(1);
hold on;
plot(Chord50.mean_ff_array_50chord, Y, 'r', 'LineWidth', 1);
xlabel('FF (%)');
ylabel('R2star (s^{-1})');
%xlim([0 50]); ylim([0 150]);

scatter(Chord50_BL.mean_ff_array_50chord, Chord50_BL.mean_r2star_array_50chord, 64, 'MarkerEdgeColor', color_cell_remote{5}, 'MarkerFaceColor', color_cell_remote{3});

yl = ylim;
xl = xlim;
text(0.3*xl(2), yl(1)+10, cat(2,'Y = ', num2str(mdl_50chord.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_50chord.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_50chord.Rsquared.Ordinary,3)), 'FontSize', 12)
title('50 chords');
mdl_50chord_BL = fitlm(Chord50_BL.mean_ff_array_50chord, Chord50_BL.mean_r2star_array_50chord);
Y_BL = Chord50_BL.mean_ff_array_50chord .* mdl_50chord_BL.Coefficients.Estimate(2) + mdl_50chord_BL.Coefficients.Estimate(1);
plot(Chord50_BL.mean_ff_array_50chord, Y_BL, 'b', 'LineWidth', 1);
text(0.3*xl(2), yl(1)+5, cat(2,'Y = ', num2str(mdl_50chord_BL.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_50chord_BL.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_50chord_BL.Rsquared.Ordinary,3)), 'FontSize', 12)
legend({'FU','', 'BL',''});

%% BL + FU (Hemo-)
figure();
mdl_50chord_hemo_n = fitlm(Chord50.mean_ff_array_50chord_hemo_n, Chord50.mean_r2star_array_50chord_hemo_n);
scatter(Chord50.mean_ff_array_50chord_hemo_n, Chord50.mean_r2star_array_50chord_hemo_n, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
Y = Chord50.mean_ff_array_50chord_hemo_n .* mdl_50chord_hemo_n.Coefficients.Estimate(2) + mdl_50chord_hemo_n.Coefficients.Estimate(1);
hold on;
plot(Chord50.mean_ff_array_50chord_hemo_n, Y, 'r', 'LineWidth', 1);
xlabel('FF (%)');
ylabel('R2star (s^{-1})');
%xlim([0 50]); ylim([0 150]);

scatter(Chord50_BL.mean_ff_array_50chord_hemo_n, Chord50_BL.mean_r2star_array_50chord_hemo_n, 64, 'MarkerEdgeColor', color_cell_remote{5}, 'MarkerFaceColor', color_cell_remote{3});

yl = ylim;
xl = xlim;
text(0.3*xl(2), yl(1)+10, cat(2,'Y = ', num2str(mdl_50chord_hemo_n.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_50chord_hemo_n.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_50chord_hemo_n.Rsquared.Ordinary,3)), 'FontSize', 12)
title('50 chords');
mdl_50chord_BL_hemo_n = fitlm(Chord50_BL.mean_ff_array_50chord_hemo_n, Chord50_BL.mean_r2star_array_50chord_hemo_n);
Y_BL = Chord50_BL.mean_ff_array_50chord_hemo_n .* mdl_50chord_BL_hemo_n.Coefficients.Estimate(2) + mdl_50chord_BL_hemo_n.Coefficients.Estimate(1);
plot(Chord50_BL.mean_ff_array_50chord_hemo_n, Y_BL, 'b', 'LineWidth', 1);
text(0.3*xl(2), yl(1)+5, cat(2,'Y = ', num2str(mdl_50chord_BL_hemo_n.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_50chord_BL_hemo_n.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_50chord_BL_hemo_n.Rsquared.Ordinary,3)), 'FontSize', 12)
legend({'FU','', 'BL',''});

%% BL + FU (Hemo+)
figure();
mdl_50chord_hemo_p = fitlm(Chord50.mean_ff_array_50chord_hemo_p, Chord50.mean_r2star_array_50chord_hemo_p);
scatter(Chord50.mean_ff_array_50chord_hemo_p, Chord50.mean_r2star_array_50chord_hemo_p, 64, 'MarkerEdgeColor', color_cell_roi{5}, 'MarkerFaceColor', color_cell_roi{3});
Y = Chord50.mean_ff_array_50chord_hemo_p .* mdl_50chord_hemo_p.Coefficients.Estimate(2) + mdl_50chord_hemo_p.Coefficients.Estimate(1);
hold on;
plot(Chord50.mean_ff_array_50chord_hemo_p, Y, 'r', 'LineWidth', 1);
xlabel('FF (%)');
ylabel('R2star (s^{-1})');
%xlim([0 50]); ylim([0 150]);

scatter(Chord50_BL.mean_ff_array_50chord_hemo_p, Chord50_BL.mean_r2star_array_50chord_hemo_p, 64, 'MarkerEdgeColor', color_cell_remote{5}, 'MarkerFaceColor', color_cell_remote{3});

yl = ylim;
xl = xlim;
text(0.3*xl(2), yl(1)+10, cat(2,'Y = ', num2str(mdl_50chord_hemo_p.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_50chord_hemo_p.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_50chord_hemo_p.Rsquared.Ordinary,3)), 'FontSize', 12)
title('50 chords');
mdl_50chord_BL_hemo_p = fitlm(Chord50_BL.mean_ff_array_50chord_hemo_p, Chord50_BL.mean_r2star_array_50chord_hemo_p);
Y_BL = Chord50_BL.mean_ff_array_50chord_hemo_p .* mdl_50chord_BL_hemo_p.Coefficients.Estimate(2) + mdl_50chord_BL_hemo_p.Coefficients.Estimate(1);
plot(Chord50_BL.mean_ff_array_50chord_hemo_p, Y_BL, 'b', 'LineWidth', 1);
text(0.3*xl(2), yl(1)+5, cat(2,'Y = ', num2str(mdl_50chord_BL_hemo_p.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl_50chord_BL_hemo_p.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl_50chord_BL_hemo_p.Rsquared.Ordinary,3)), 'FontSize', 12)
legend({'FU','', 'BL',''});