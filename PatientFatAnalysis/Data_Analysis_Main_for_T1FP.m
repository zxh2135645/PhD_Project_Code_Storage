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

time_points = {'FU'};

label_lge = sequence_label{1};
label_t2star = sequence_label{2};

metrics_save_dir = cat(2, base_dir, '/Results/');
if ~exist(metrics_save_dir, 'dir')
   mkdir(metrics_save_dir); 
end

% Read excel file
% T = readtable(cat(2, base_dir, '/STEMI_with_IMH-1.xlsx'), 'VariableNamingRule', 'preserve');
T = readtable(cat(2, base_dir, '/STEMI_with_IMH-1.xlsx'));
id_array = T.AnonymizationID;
hemo_array = T.withOrWithout;

excel_names = {'484060000008', '484060000009', '484060000045', '484060000010', '484060000018', '484060000052',...
    '484060000054', '484060000056', '484060000029', '484060000060',...
    '484060000030', '484060000031', '484060000033', '484060000041'};
excel_names = {'484060000001'};
good_names = {'484060000003', '484060000012', '484060000015', '484060000016', '484060000020', '484060000021',...
    '484060000022', '484060000023', '484060000028', '484060000032',...
    '484060000037', '484060000040', '484060000055', '484060000063','484060000036'};
ok_names = {'484060000002', '484060000006', '484060000017', '484060000058', '484060000027', '484060000035', ...
    '484060000039'};
bad_names = {'484060000001', '484060000011', '484060000013'};

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

se = strel('disk', 1);
nhood = [1 1 1; 1 1 1; 1 1 1];

for n = 1:(length(Names)-1)
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
    
    if any(contains(excel_names, name))
    %if any(contains(good_names, name))
    name_for_table_searching = insertAfter(name, 6, '-');
    row = find(contains(id_array,name_for_table_searching));
    
    IMH_cell = table2cell(T(row, 13)); % IMH
    IMH = IMH_cell{1};
    for tp = 1:length(time_points)
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
                            caxis(ax1, [0 100]); caxis(ax2, [-2 10]); linkprop([ax1 ax2], 'Position');
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
                
            % No need to reorder for analysis
            
            roi_in_myo_r2star = mask_myocardium_3D .* freeROIMask_3D;
            remote_in_myo_r2star = mask_myocardium_3D .* myoRefMask_3D;
            roi_r2star = roi_in_myo_r2star .* r2star;
            remote_r2star = remote_in_myo_r2star .* r2star;
            myo_r2star = mask_myocardium_3D;
            
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
