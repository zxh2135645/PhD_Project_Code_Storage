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
FF_glob = glob(cat(2, base_dir, '/FF_Data/*'));
count = 1;
Names = cell(length(FF_glob)-1, 1); 
for i = 1:length(FF_glob)
    strings = strsplit(FF_glob{i},'/');
    if ~isempty(regexp(strings{end-1}, '\d', 'once'))
        Names{count} = strings{end-1};
        count = count + 1;
    end
end

sequence_label = {'T2_MULTIECHO'};
anatomy_label = {'BloodPool', 'excludeArea', 'freeROI', 'Heart', 'Myocardium', 'MyoReference', 'noReflowArea'}; 

names_to_rule_out = {'SS10113', '117898'};
RuleOutLabel = NameRuleOutFunc(Names, names_to_rule_out);
Names = Names(RuleOutLabel == 0);

%
name_check = {'108516'};
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


output_label = {'T2_MULTIECHO'};
save_dir = GetFullPath(cat(2, base_dir, '/Analysis/'));
data_save_dir = cat(2, base_dir, '/data/');

time_points = {'Acute', 'Chronic'};
label_t2star = sequence_label{1};

metrics_save_dir = cat(2, base_dir, '/Results/');
if ~exist(metrics_save_dir, 'dir')
   mkdir(metrics_save_dir); 
end

% IMH+
Names = {'SS4887', '108516', 'SS4570', '121381', '139460', '110868', '122341', 'SS4542', 'SS5559'};

% IMH-
Names = {'SS4922', 'SS5357', '111525', '114627', 'SS4510', '112111', '131534', '139461', '140381'};
%% Main Body (Need to run twice, one with FU, the other with BL and BL2)
% Array initialization
mean_ff_roi_array = [];
sd_ff_roi_array = [];
mean_ff_remote_array = [];
sd_ff_remote_array = [];

mean_r2star_roi_array = [];
sd_r2star_roi_array = [];
mean_r2star_remote_array = [];
sd_r2star_remote_array = [];

name_label = {};
name_label_remote = {};
slice_count = 1;
slice_count_remote = 1;
vec = @(x) x(:);

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

hemo_label_array = [];

for n = 1:length(Names)
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
    
    %for tp = 1:length(time_points)
    for tp = 2:2
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
            exclude_glob = glob(cat(2, tp_dir, label_t2star, '/',anatomy_label{2}, '/*'));

            load(myo_glob{1});
            load(roi_glob{2});
            load(remote_glob{2});
            load(exclude_glob{2});
            load(cat(2, tp_dir, label_t2star, '/', label_t2star, '_Index.mat')); % glob_names
            load(cat(2, tp_dir, label_t2star, '/', label_t2star, '_SliceLoc.mat')); % slc_array
            load(cat(2, tp_dir, label_t2star, '/', label_t2star, '_vol_img_3D_te.mat'));

           

            slc_array_t2star = slc_array;
            
            num_array = regexp(glob_names,'\d*','Match');

            ff_glob = glob(cat(2,  base_dir, '/FF_Data/',  name, '/', time_point, '/*_', num2str(str2num(num2str(num_array{1}))), '.mat'));
            ff_map = load(ff_glob{1}, 'fwmc_ff');
            ff = ff_map.fwmc_ff;
            
            r2star_glob = glob(cat(2,  base_dir, '/FF_Data/',  name, '/', time_point, '/*_', num2str(str2num(num2str(num_array{1}))), '.mat'));
            r2star_map = load(r2star_glob{1}, 'fwmc_r2star');
            r2star = r2star_map.fwmc_r2star;
            
            
            mask_myocardium_3D = imerode(mask_myocardium_3D, se);
            
            roi_in_myo_ff = mask_myocardium_3D .* freeROIMask_3D_te;
            remote_in_myo_ff = mask_myocardium_3D .* myoRef_3D_te;
            roi_ff = roi_in_myo_ff .* ff;
            remote_ff = remote_in_myo_ff .* ff;
            myo_ff = mask_myocardium_3D;
            
            
            % TODO add QC for myocardium map here
            % f_exclusion_pts = cat(2, name_data_save_dir, '/ExclusionPts_', name, '_', time_point, '.mat');
            % figure('Position', [100 100 800 600]);
            % 
            % if exist(f_exclusion_pts, 'file')
            %     load(f_exclusion_pts);
            % 
            % 
            %     for slc = 1:size(roi_in_myo_ff,3)
            %         exclude_mask = zeros(size(mask_myocardium_3D,1),size(mask_myocardium_3D,2));
            % 
            %         x_array = exclusion_points{slc, 1};
            %         y_array = exclusion_points{slc, 2};
            % 
            %         for pts = 1:length(y_array)
            %             exclude_mask(round(x_array(pts)),round(y_array(pts))) = 1;
            %         end
            % 
            %         exclude_mask = imdilate(exclude_mask, nhood);
            % 
            %         freeROIMask_3D(:,:,slc) = freeROIMask_3D(:,:,slc) .* ~exclude_mask;
            % 
            %         roi_in_myo_ff(:,:,slc) = mask_myocardium_3D(:,:,slc) .* freeROIMask_3D(:,:,slc);
            %         roi_ff(:,:,slc) = roi_in_myo_ff(:,:,slc) .* ff(:,:,slc);
            %         % myo_ff(:,:,slc) = mask_myocardium_3D(:,:,slc);
            %     end
            % 
            % else
            %     exclusion_points = cell(size(roi_in_myo_ff,3), 2);
            % 
            % 
            %     for slc = 1:size(roi_in_myo_ff,3)
            % 
            %         if any(any(roi_in_myo_ff(:,:,slc)))
            %             exclude_mask = zeros(size(mask_myocardium_3D,1),size(mask_myocardium_3D,2));
            %             flag = 1;
            %             y_array = [];
            %             x_array = [];
            %             count = 0;
            %             while flag
            %                 % crop image
            %                 stats = regionprops(mask_myocardium_3D(:,:,slc));
            %                 centroid = stats.Centroid;
            %                 X = round(centroid(1)-size(roi_in_myo_ff,2)/8);
            %                 Y = round(centroid(2)-size(roi_in_myo_ff,1)/8);
            %                 w = size(roi_in_myo_ff,2)/4;
            %                 h = size(roi_in_myo_ff,1)/4;
            % 
            %                 roi_in_myo_ff_crop = imcrop(roi_in_myo_ff(:,:,slc), [X,Y,w,h]);
            %                 r2star_crop = imcrop(r2star(:,:,slc), [X,Y,w,h]);
            %                 roi_ff_crop = imcrop(roi_ff(:,:,slc), [X,Y,w,h]);
            %                 myo_ff_crop = imcrop(myo_ff(:,:,slc), [X,Y,w,h]);
            % 
            %                 roi_in_myo_ff_nan = roi_in_myo_ff_crop;
            %                 roi_in_myo_ff_nan(roi_in_myo_ff_crop == 0) = nan;
            % 
            %                 ax1 = axes;
            %                 imagesc(r2star_crop); pbaspect([size(r2star_crop, 2), size(r2star_crop, 1) 1]);
            %                 colormap(ax1, 'gray');
            %                 set(ax1, 'xticklabel', []); set(ax1, 'yticklabel', []);
            %                 ax2 = axes;
            % 
            %                 imagesc(ax2, roi_in_myo_ff_nan .* roi_ff_crop, 'AlphaData', myo_ff_crop);
            %                 pbaspect([size(r2star_crop, 2), size(r2star_crop, 1) 1]); colormap(ax2, 'cool');
            %                 caxis(ax1, [0 100]); caxis(ax2, [-2 10]); linkprop([ax1 ax2], 'Position');
            %                 ax2.Visible = 'off';
            %                 if count == 0
            %                     cb = colorbar; title(cb, 'FF (%)');
            %                 end
            % 
            %                 [y,x] = getpts;
            % 
            %                 if length(y) == 1
            %                     txt = 'y';
            %                 else
            %                     y = y + X - 1;
            %                     x = x + Y - 1;
            % 
            %                     y_array = [y_array; y];
            %                     x_array = [x_array; x];
            % 
            %                     for pts = 1:length(y)
            %                         exclude_mask(round(x_array(pts)),round(y_array(pts))) = 1;
            %                     end
            % 
            %                     exclude_mask = imdilate(exclude_mask, nhood);
            % 
            %                     freeROIMask_3D(:,:,slc) = freeROIMask_3D(:,:,slc) .* ~exclude_mask;
            % 
            %                     roi_in_myo_ff(:,:,slc) = mask_myocardium_3D(:,:,slc) .* freeROIMask_3D(:,:,slc);
            %                     roi_ff(:,:,slc) = roi_in_myo_ff(:,:,slc) .* ff(:,:,slc);
            %                     % myo_ff(:,:,slc) = mask_myocardium_3D(:,:,slc);
            % 
            %                     r2star_crop = imcrop(r2star(:,:,slc), [X,Y,w,h]);
            %                     roi_ff_crop = imcrop(roi_ff(:,:,slc), [X,Y,w,h]);
            % 
            %                     imagesc(roi_in_myo_ff_nan .* roi_ff_crop); 
            %                     pbaspect([size(r2star_crop, 2), size(r2star_crop, 1) 1]); caxis([-2 10]);
            % 
            %                     prompt = 'Are you happy??? (y/n): ';
            %                     txt = input(prompt,"s");
            %                 end
            % 
            %                 if strcmp(txt, 'y')
            %                     flag = 0;
            %                 elseif strcmp(txt, 'n')
            %                     flag = 1;
            %                     count = count + 1;
            %                 end
            %             end
            %             exclusion_points{slc, 1} = x_array;
            %             exclusion_points{slc, 2} = y_array;
            %         end
            %     end
            % 
            %     save(f_exclusion_pts, 'exclusion_points');
            % end
                
            % No need to reorder for analysis
            
            roi_in_myo_r2star = mask_myocardium_3D .* freeROIMask_3D_te;
            remote_in_myo_r2star = mask_myocardium_3D .* myoRef_3D_te;
            roi_r2star = roi_in_myo_r2star .* r2star;
            remote_r2star = remote_in_myo_r2star .* r2star;
            myo_r2star = mask_myocardium_3D;
            
            tp_dir2 = cat(2, name_save_dir, '/', time_point, '/');
            if ~exist(tp_dir2, 'dir')
                mkdir(tp_dir2);
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
            remote_in_myo_ff(remote_in_myo_ff == 0) = nan;
            
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

            % Mean - 2SD to see if it's above 1%
            hemo_label = zeros(1, size(r2star_roi_masked, 3));
            for slc = 1:size(r2star_roi_masked, 3)
                thresh = mean(nonzeros(vol_img_3D_te(:,:,slc) .* myoRef_3D_te(:,:,slc))) - 2*std(nonzeros(vol_img_3D_te(:,:,slc) .* myoRef_3D_te(:,:,slc)));
                hemo = (vol_img_3D_te(:,:,slc) <= thresh) .* roi_in_myo_ff(:,:,slc) .* ~excludeMask_3D_te(:,:,slc);
                hemo_perc = sum(vec(hemo==1)) ./ sum(vec(mask_myocardium_3D(:,:,slc)==1));
                if hemo_perc > 0.01
                    hemo_label(slc) = 1;
                end
            end
            

            % figure();
            for slc = 1:size(r2star_roi_masked, 3)

                if any(vec(ff_roi_masked_px(:,:,slc)))

                    mean_r2star_roi_array = [mean_r2star_roi_array, mean(vec(r2star_roi_masked_px(:,:,slc)), 'omitnan')];
                    sd_r2star_roi_array = [sd_r2star_roi_array, std(vec(r2star_roi_masked_px(:,:,slc)), 'omitnan')];
                    mean_ff_roi_array = [mean_ff_roi_array, mean(vec(ff_roi_masked_px(:,:,slc)), 'omitnan')];
                    sd_ff_roi_array = [sd_ff_roi_array, std(vec(ff_roi_masked_px(:,:,slc)), 'omitnan')];
                    
                    hemo_label_array = [hemo_label_array, hemo_label];
                    name_label{slice_count} = cat(2, name, '_', time_point, '_Slice', num2str(slc));
                    slice_count = slice_count + 1;

                    if (hemo_label(slc)>0)
                        fprintf('Hemo+ Name: %s, Timepoint: %s\n', name, time_point);
                        mean_r2star_roi_array_hemo_p = [mean_r2star_roi_array_hemo_p, mean(vec(r2star_roi_masked_px(:,:,slc)), 'omitnan')];
                        sd_r2star_roi_array_hemo_p = [sd_r2star_roi_array_hemo_p, std(vec(r2star_roi_masked_px(:,:,slc)), 'omitnan')];
                        mean_ff_roi_array_hemo_p = [mean_ff_roi_array_hemo_p, mean(vec(ff_roi_masked_px(:,:,slc)), 'omitnan')];
                        sd_ff_roi_array_hemo_p = [sd_ff_roi_array_hemo_p, std(vec(ff_roi_masked_px(:,:,slc)), 'omitnan')];
                    else
                        mean_r2star_roi_array_hemo_n = [mean_r2star_roi_array_hemo_n, mean(vec(r2star_roi_masked_px(:,:,slc)), 'omitnan')];
                        sd_r2star_roi_array_hemo_n = [sd_r2star_roi_array_hemo_n, std(vec(r2star_roi_masked_px(:,:,slc)), 'omitnan')];
                        mean_ff_roi_array_hemo_n = [mean_ff_roi_array_hemo_n, mean(vec(ff_roi_masked_px(:,:,slc)), 'omitnan')];
                        sd_ff_roi_array_hemo_n = [sd_ff_roi_array_hemo_n, std(vec(ff_roi_masked_px(:,:,slc)), 'omitnan')];
                    end
                end

                if any(vec(ff_remote_masked(:,:,slc)))
                    mean_r2star_remote_array = [mean_r2star_remote_array, mean(vec(r2star_remote_masked(:,:,slc)), 'omitnan')];
                    sd_r2star_remote_array = [sd_r2star_remote_array, std(vec(r2star_remote_masked(:,:,slc)), 'omitnan')];
                    mean_ff_remote_array = [mean_ff_remote_array, mean(vec(ff_remote_masked(:,:,slc)), 'omitnan')];
                    sd_ff_remote_array = [sd_ff_remote_array, std(vec(ff_remote_masked(:,:,slc)), 'omitnan')];

                    name_label_remote{slice_count_remote} = cat(2, name, '_', time_point, '_Slice', num2str(slc));
                    slice_count_remote = slice_count_remote + 1;
                end

                %subplot(1,2,1); imagesc(ff_roi_masked_px(:,:,slc)); axis image; caxis([0 20]); colorbar;
                %subplot(1,2,2); imagesc(r2star_roi_masked_px(:,:,slc)); axis image; caxis([0 100]); colorbar;
                
            end
        end
    end
    close all;
end
