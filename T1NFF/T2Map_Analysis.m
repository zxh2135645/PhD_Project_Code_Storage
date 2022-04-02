clear all;
close all;

%% T2 mapping analysis for probing edema
addpath('../function/');
addpath('../AHA16Segment/');
addpath('../function/demon_registration_version_8f_winOS/');
base_dir = uigetdir;
contour_glob = glob(cat(2, base_dir, '/ContourData/*'));
%Names = ExtractNames(contour_glob);
Names = {'Merry', 'Ryn', 'Mojave', 'Sahara', 'ZZ', 'Tina', 'Sunny', 'Queenie', 'Hope', 'Gobi', 'Felicity', 'Evelyn', '18D15', '18D16'};
time_points = {'7D', '8WK', '12WK', '14WK', '4MO', '6MO', '9MO', '1YR', '15YR'};
time_points = {'6MO', '9MO', '1YR', '15YR'};
% time_points = {'1YR', '15YR'};
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
 
sequence_label = {'T1', 'T2star', 'LGE', 'T2'};
anatomy_label = {'BloodPool', 'excludeArea', 'freeROI', 'Heart', 'Myocardium', 'MyoReference', 'noReflowArea'}; 
%name_check = 'Evelyn';
%starting_point = find(strcmp(name_check, Names),1);
 
save_dir = cat(2, base_dir, '/img/');
data_save_dir = cat(2, base_dir, '/data/');
 
label_t2 = sequence_label{4};
label_t2star = sequence_label{2};
 
metrics_save_dir = cat(2, base_dir, '/Results/');
if ~exist(metrics_save_dir, 'dir')
   mkdir(metrics_save_dir); 
end
 
% Before analysis, parse pre_QualControl
load(cat(2, metrics_save_dir, 'pre_QualControl.mat'));
 
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
 
%%
mean_T2_array = [];
sd_T2_array = [];
mean_T2_array_remote = [];
sd_T2_array_remote = [];
name_label = {};
count = 1;

mean_T2_array_hemo = [];
sd_T2_array_hemo = [];
mean_T2_array_peri = [];
sd_T2_array_peri = [];

mean_T2_array_nonhemo = [];
sd_T2_array_nonhemo = [];
mean_T2_array_remote_nonhemo = [];
sd_T2_array_remote_nonhemo = [];
name_label_nonhemo = {};
count_nonhemo = 1;

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
    tp_count = 0;
    for tp = 1:length(time_points)
    %for tp = 4:4
        time_point = time_points{end-tp+1};
        tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', name, '_', time_point,  '/');
        if ~exist(tp_dir, 'dir') || ~exist(cat(2, tp_dir, label_t2, '/', label_t2, '_vol_img_3D.mat'))
            disp(cat(2, 'No folder at: ', name, ' ', time_point));
        else
            % T2
            tp_count = tp_count+1;
            myo_glob = glob(cat(2, tp_dir, label_t2, '/', anatomy_label{5}, '/*'));
            roi_glob = glob(cat(2, tp_dir, label_t2, '/',anatomy_label{3}, '/*'));
            remote_glob = glob(cat(2, tp_dir, label_t2, '/',anatomy_label{6}, '/*'));
            hemo_glob = glob(cat(2, tp_dir, label_t2, '/',anatomy_label{7}, '/*'));
            
            
            load(cat(2, tp_dir, label_t2, '/', label_t2, '_vol_img_3D.mat'));
            load(myo_glob{1});
            load(roi_glob{1});
            load(remote_glob{1});
            load(hemo_glob{1});
            load(cat(2, tp_dir, label_t2, '/', label_t2, '_SliceLoc.mat'));
            [slc_array_t2, idx_reordered] = sort(slc_array);
            
            se = strel('disk', 1);
            mask_myocardium_3D_eroded = mask_myocardium_3D;
            for ii = 1:size(mask_myocardium_3D,3)
                mask_myocardium_3D_eroded(:,:,ii) = imerode(mask_myocardium_3D(:,:,ii), se);
                %mask_myocardium_3D
            end
            roi_in_myo_t2 = mask_myocardium_3D_eroded .* freeROIMask_3D;
            remote_in_myo_t2 = mask_myocardium_3D_eroded .* myoRefMask_3D;
            roi_t2 = roi_in_myo_t2 .* vol_img_3D;
            hemo_core = noReflowMask_3D .* roi_in_myo_t2;
%             hemo_in_roi_t2 = roi_t2 .* noReflowMask_3D; % Hemo Core
%             hemo_peri_roi_t2 = roi_in_myo_t2 .* ~noReflowMask_3D; 
            remote_t2 = remote_in_myo_t2 .* vol_img_3D;
            t2 = vol_img_3D;
            myo_t2 = mask_myocardium_3D_eroded;
            
%             roi_in_myo_t2 = roi_in_myo_t2(:,:,idx_reordered);
%             remote_in_myo_t2 = remote_in_myo_t2(:,:,idx_reordered);
%             roi_t2 = roi_t2(:,:,idx_reordered);
%             remote_t2 = remote_t2(:,:,idx_reordered);
%             t2 = t2(:,:,idx_reordered);
%             myo_t2 = myo_t2(:,:,idx_reordered);
                        
            
            tp_dir2 = cat(2, name_save_dir, '/', name, '_', time_point, '/');
            if ~exist(tp_dir2, 'dir')
                mkdir(tp_dir2);
            end

            % if condition
            center_mask_fname = cat(2, name_data_save_dir, '/CenterLine_', name, '_', time_point, '.mat');
                        
%             %if ~exist(center_mask_fname, 'file')
%                 center_mask_ff = zeros(size(roi_in_myo_ff));
%                 BW_skel = zeros(size(roi_in_myo_ff));
%                 %caxis_rg = [0 1];
%                 %for i = 1:size(roi_in_myo_ff, 3)
%                     % Draw center line on myocardium ff
%                     %disp(cat(2, 'Center line roi: ', name, '  ', time_point, 'Slice', num2str(i)));
%  
%                     %moving = myo_ff(:,:,i);
%                     %center_mask_ff(:,:,i) = Func_DrawCenterLine(moving, caxis_rg);
%                     %BW_skel(:,:,i) = bwmorph(moving, 'skel', Inf);
%                     %center_mask_ff(:,:,i) = imfill(BW_skel(:,:,i), 'hole');
%                 %end
%                 
%                 figure();
%                 sz = ceil(sqrt(size(roi_in_myo_ff, 3)));
%                 se = strel('disk', 1);
%                 for i = 1:size(roi_in_myo_ff, 3)
%                     subplot(sz,sz,i);
%                     myo_ff_eroded = imerode(myo_ff(:,:,i), se);
%                     BW_skel(:,:,i) = bwmorph(myo_ff_eroded, 'skel', Inf);
%                     center_mask_ff(:,:,i) = imfill(BW_skel(:,:,i), 'hole');
%                     
%                     % Found 18D16 8WK is not a closed shape (Hard-Coded)
%                     % Felicity 6MO
%                     if (n == 14 && tp == 8 && i == 3) || (n == 11 && tp == 4 && i == 3)
%                         center_mask_ff(:,:,3) = bwconvhull(center_mask_ff(:,:,3));
%                     end
%                     epi = myo_ff_eroded - center_mask_ff(:,:,i) > 0;
%                     endo = center_mask_ff(:,:,i) + myo_ff_eroded > 1;
%                     imagesc(endo*2 + epi);
%                     colormap(brewermap([],'*RdYlBu'));
%                 end
%                 save(center_mask_fname, 'center_mask_ff');
%                 saveas(gcf, cat(2, tp_dir2, 'CenterLineMask.png'))
            %else
            %    load(center_mask_fname);
            %end
            %status = status_check(n).status(tp_count,:);
            % AHA Segment
            Segn = 50;
            Groove = 0;
            %figure();
            disp(name);
            disp(sum(hemo_core(:)) / sum(roi_in_myo_t2(:)));
            disp(sum(hemo_core(:)) / sum(mask_myocardium_3D_eroded(:)));
            hemo_perc  = sum(hemo_core(:)) / sum(mask_myocardium_3D_eroded(:));
            for slc = 1:size(t2, 3)
                temp = noReflowMask_3D(:,:,slc);
                %subplot(2,3,slc);
                %imagesc(temp); axis off;
                %title(cat(2, name, '_', time_point, '_Slice', num2str(slc)));
                
                
                % if any(temp(:))
                if hemo_perc > 0.03 || any(temp(:)) % Hard-coded here 
                    mean_T2_array = [mean_T2_array, mean(nonzeros(roi_in_myo_t2(:,:,slc) .* t2(:,:,slc)))];
                    sd_T2_array = [sd_T2_array, std(nonzeros(roi_in_myo_t2(:,:,slc) .* t2(:,:,slc)))];
                    mean_T2_array_remote = [mean_T2_array_remote, mean(nonzeros(remote_in_myo_t2(:,:,slc) .* t2(:,:,slc)))];
                    sd_T2_array_remote = [sd_T2_array_remote, std(nonzeros(remote_in_myo_t2(:,:,slc) .* t2(:,:,slc)))];
                    name_label{count} = cat(2, name, '_', time_point, '_Slice', num2str(slc));
                    count = count + 1;
                    
                    % hemo core
                    mean_T2_array_hemo = [mean_T2_array_hemo, mean(nonzeros(roi_in_myo_t2(:,:,slc) .* t2(:,:,slc) .* hemo_core(:,:,slc)))];
                    sd_T2_array_hemo = [sd_T2_array_hemo, std(nonzeros(roi_in_myo_t2(:,:,slc) .* t2(:,:,slc) .* hemo_core(:,:,slc)))];
                    
                    mean_T2_array_peri = [mean_T2_array_peri, mean(nonzeros(roi_in_myo_t2(:,:,slc) .* t2(:,:,slc) .* ~hemo_core(:,:,slc)))];
                    sd_T2_array_peri = [sd_T2_array_peri, std(nonzeros(roi_in_myo_t2(:,:,slc) .* t2(:,:,slc) .* ~hemo_core(:,:,slc)))];
                    
                else
                    mean_T2_array_nonhemo = [mean_T2_array_nonhemo, mean(nonzeros(roi_in_myo_t2(:,:,slc) .* t2(:,:,slc)))];
                    sd_T2_array_nonhemo = [sd_T2_array_nonhemo, std(nonzeros(roi_in_myo_t2(:,:,slc) .* t2(:,:,slc)))];
                    mean_T2_array_remote_nonhemo = [mean_T2_array_remote_nonhemo, mean(nonzeros(remote_in_myo_t2(:,:,slc) .* t2(:,:,slc)))];
                    sd_T2_array_remote_nonhemo = [sd_T2_array_remote_nonhemo, std(nonzeros(remote_in_myo_t2(:,:,slc) .* t2(:,:,slc)))];
                    name_label_nonhemo{count_nonhemo} = cat(2, name, '_', time_point, '_Slice', num2str(slc));
                    count_nonhemo = count_nonhemo + 1;
                end
            end
            
        end
    end
    %close all;
end

%% 
figure();
x = [1, 2, 3, 4, 5, 6];
y = [0.1*mean(mean_T2_array); 0.1*nanmean(mean_T2_array_hemo); 0.1*mean(mean_T2_array_peri); 0.1*mean(mean_T2_array_remote); 0.1*mean(mean_T2_array_nonhemo); 0.1*mean(mean_T2_array_remote_nonhemo)];
bar(x, y); ylim([20 50]);
xticklabels({'MI (Hemo+)', 'Hemo Core', 'Hemo Peri', 'Remote (Hemo+)', 'MI (Hemo-)', 'Remote (Hemo-)'});

err = [0.1*std(mean_T2_array); 0.1*nanstd(mean_T2_array_hemo); 0.1*std(mean_T2_array_peri); 0.1*std(mean_T2_array_remote); 0.1*std(mean_T2_array_nonhemo); 0.1*std(mean_T2_array_remote_nonhemo)];
hold on;
er = errorbar(x,y, err);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
%% R2star and FF
mean_R2star_array = [];
sd_R2star_array = [];
mean_R2star_array_remote = [];
sd_R2star_array_remote = [];
name_label = {};
count = 1;

mean_R2star_array_hemo = [];
sd_R2star_array_hemo = [];
mean_R2star_array_peri = [];
sd_R2star_array_peri = [];

mean_R2star_array_nonhemo = [];
sd_R2star_array_nonhemo = [];
mean_R2star_array_remote_nonhemo = [];
sd_R2star_array_remote_nonhemo = [];
name_label_nonhemo = {};
count_nonhemo = 1;

mean_ff_array = [];
sd_ff_array = [];
mean_ff_array_remote = [];
sd_ff_array_remote = [];

mean_ff_array_hemo = [];
sd_ff_array_hemo = [];
mean_ff_array_peri = [];
sd_ff_array_peri = [];

mean_ff_array_nonhemo = [];
sd_ff_array_nonhemo = [];
mean_ff_array_remote_nonhemo = [];
sd_ff_array_remote_nonhemo = [];

name_label_fat_nonhemo = {};
name_label_fat_hemo = {};
name_label_nonfat_nonhemo = {};
name_label_nonfat_hemo = {};

hemo_perc_array = [];
ff_array = [];
ct = 1;
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
    tp_count = 0;
    for tp = 1:length(time_points)
    %for tp = 4:4
        time_point = time_points{end-tp+1};
        tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', name, '_', time_point,  '/');
        if ~exist(tp_dir, 'dir') || ~exist(cat(2, tp_dir, label_t2star, '/', label_t2star, '_vol_img_3D.mat'))
            disp(cat(2, 'No folder at: ', name, ' ', time_point));
        else
            % T2
            tp_count = tp_count+1;
            myo_glob = glob(cat(2, tp_dir, label_t2star, '/', anatomy_label{5}, '/*'));
            roi_glob = glob(cat(2, tp_dir, label_t2star, '/',anatomy_label{3}, '/*'));
            remote_glob = glob(cat(2, tp_dir, label_t2star, '/',anatomy_label{6}, '/*'));
            hemo_glob = glob(cat(2, tp_dir, label_t2star, '/',anatomy_label{7}, '/*'));
            
            
            load(cat(2, tp_dir, label_t2star, '/', label_t2star, '_vol_img_3D.mat'));
            load(myo_glob{1});
            load(roi_glob{1});
            load(roi_glob{2});
            load(remote_glob{1});
            load(hemo_glob{1});
            load(cat(2, tp_dir, label_t2star, '/', label_t2star, '_Index.mat'));
            load(cat(2, tp_dir, label_t2star, '/', label_t2star, '_SliceLoc.mat'));
            [slc_array_t2, idx_reordered] = sort(slc_array);
            
            ff_map = cell(1, length(glob_names));
            r2star_map = cell(1, length(glob_names));
            for f = 1:length(ff_map)
                ff_map{f} = load(cat(2,  base_dir, '/FF_Data/',  name, '/', name, '_', time_point, '_', glob_names{f}, '.mat'), 'fwmc_ff');
                r2star_map{f} = load(cat(2,  base_dir, '/FF_Data/',  name, '/', name, '_', time_point, '_', glob_names{f}, '.mat'), 'fwmc_r2star');                
            end
            
            % convert ff_map to matrix
            ff = zeros(size(ff_map{1}.fwmc_ff,1), size(ff_map{1}.fwmc_ff, 2), length(ff_map));
            r2star = zeros(size(r2star_map{1}.fwmc_r2star,1), size(r2star_map{2}.fwmc_r2star, 2), length(r2star_map));
            for f = 1:length(ff_map)
                ff(:,:,f) = ff_map{f}.fwmc_ff;
                r2star(:,:,f) = r2star_map{f}.fwmc_r2star;
            end  
            
            ff(ff > 100) = 100;
            ff(ff < 0) = 0;
            r2star(r2star>200) = 200;
            
            se = strel('disk', 1);
            mask_myocardium_3D_eroded = mask_myocardium_3D;
            for ii = 1:size(mask_myocardium_3D,3)
                mask_myocardium_3D_eroded(:,:,ii) = imerode(mask_myocardium_3D(:,:,ii), se);
                %mask_myocardium_3D
            end
            
            roi_in_myo_t2star = mask_myocardium_3D_eroded .* freeROIMask_3D;
            remote_in_myo_t2star = mask_myocardium_3D_eroded .* myoRefMask_3D;
            roi_t2star = roi_in_myo_t2star .* vol_img_3D;
            hemo_core = freeROIMask_3D_te8 .* roi_in_myo_t2star;
%             hemo_in_roi_t2 = roi_t2 .* noReflowMask_3D; % Hemo Core
%             hemo_peri_roi_t2 = roi_in_myo_t2 .* ~noReflowMask_3D; 
            remote_t2star = remote_in_myo_t2star .* vol_img_3D;
            t2star = vol_img_3D;
            myo_t2star = mask_myocardium_3D_eroded;
            ff_core = roi_in_myo_t2star .* ff > 10;
%             roi_in_myo_t2 = roi_in_myo_t2(:,:,idx_reordered);
%             remote_in_myo_t2 = remote_in_myo_t2(:,:,idx_reordered);
%             roi_t2 = roi_t2(:,:,idx_reordered);
%             remote_t2 = remote_t2(:,:,idx_reordered);
%             t2 = t2(:,:,idx_reordered);
%             myo_t2 = myo_t2(:,:,idx_reordered);
                        
            
            tp_dir2 = cat(2, name_save_dir, '/', name, '_', time_point, '/');
            if ~exist(tp_dir2, 'dir')
                mkdir(tp_dir2);
            end

            % if condition
            center_mask_fname = cat(2, name_data_save_dir, '/CenterLine_', name, '_', time_point, '.mat');
                        
            % AHA Segment
            Segn = 50;
            Groove = 0;
            %figure();
            disp(name);
            disp(sum(hemo_core(:)) / sum(roi_in_myo_t2star(:)));
            disp(sum(hemo_core(:)) / sum(mask_myocardium_3D_eroded(:)));
            hemo_perc  = sum(hemo_core(:)) / sum(mask_myocardium_3D_eroded(:));
            
            for slc = 1:size(r2star, 3)
                temp = freeROIMask_3D_te8(:,:,slc);
                %subplot(2,3,slc);
                %imagesc(temp); axis off;
                %title(cat(2, name, '_', time_point, '_Slice', num2str(slc)));
                
                
                % if any(temp(:))
                if hemo_perc > 0.03 || any(temp(:)) % Hard-coded here 
                    mean_R2star_array = [mean_R2star_array, mean(nonzeros(roi_in_myo_t2star(:,:,slc) .* r2star(:,:,slc)))];
                    sd_R2star_array = [sd_R2star_array, std(nonzeros(roi_in_myo_t2star(:,:,slc) .* r2star(:,:,slc)))];
                    mean_R2star_array_remote = [mean_R2star_array_remote, mean(nonzeros(remote_in_myo_t2star(:,:,slc) .* r2star(:,:,slc)))];
                    sd_R2star_array_remote = [sd_R2star_array_remote, std(nonzeros(remote_in_myo_t2star(:,:,slc) .* r2star(:,:,slc)))];
                    name_label{count} = cat(2, name, '_', time_point, '_Slice', num2str(slc));
                    count = count + 1;
                    
                    % hemo core
                    mean_R2star_array_hemo = [mean_R2star_array_hemo, mean(nonzeros(roi_in_myo_t2star(:,:,slc) .* r2star(:,:,slc) .* hemo_core(:,:,slc)))];
                    sd_R2star_array_hemo = [sd_R2star_array_hemo, std(nonzeros(roi_in_myo_t2star(:,:,slc) .* r2star(:,:,slc) .* hemo_core(:,:,slc)))];
                    
                    mean_R2star_array_peri = [mean_R2star_array_peri, mean(nonzeros(roi_in_myo_t2star(:,:,slc) .* r2star(:,:,slc) .* ~hemo_core(:,:,slc)))];
                    sd_R2star_array_peri = [sd_R2star_array_peri, std(nonzeros(roi_in_myo_t2star(:,:,slc) .* r2star(:,:,slc) .* ~hemo_core(:,:,slc)))];
                    
                else
                    mean_R2star_array_nonhemo = [mean_R2star_array_nonhemo, mean(nonzeros(roi_in_myo_t2star(:,:,slc) .* r2star(:,:,slc)))];
                    sd_R2star_array_nonhemo = [sd_R2star_array_nonhemo, std(nonzeros(roi_in_myo_t2star(:,:,slc) .* r2star(:,:,slc)))];
                    mean_R2star_array_remote_nonhemo = [mean_R2star_array_remote_nonhemo, mean(nonzeros(remote_in_myo_t2star(:,:,slc) .* r2star(:,:,slc)))];
                    sd_R2star_array_remote_nonhemo = [sd_R2star_array_remote_nonhemo, std(nonzeros(remote_in_myo_t2star(:,:,slc) .* r2star(:,:,slc)))];
                    name_label_nonhemo{count_nonhemo} = cat(2, name, '_', time_point, '_Slice', num2str(slc));
                    count_nonhemo = count_nonhemo + 1;
                end
                
                if hemo_perc > 0.03 || any(temp(:)) % Hard-coded here 
                    mean_ff_array = [mean_ff_array, mean(nonzeros(roi_in_myo_t2star(:,:,slc) .* ff(:,:,slc)))];
                    sd_ff_array = [sd_ff_array, std(nonzeros(roi_in_myo_t2star(:,:,slc) .* ff(:,:,slc)))];
                    mean_ff_array_remote = [mean_ff_array_remote, mean(nonzeros(remote_in_myo_t2star(:,:,slc) .* ff(:,:,slc)))];
                    sd_ff_array_remote = [sd_ff_array_remote, std(nonzeros(remote_in_myo_t2star(:,:,slc) .* ff(:,:,slc)))];
                    %name_label{count} = cat(2, name, '_', time_point, '_Slice', num2str(slc));
                    %count = count + 1;
                    
                    % hemo core
                    mean_ff_array_hemo = [mean_ff_array_hemo, mean(nonzeros(roi_in_myo_t2star(:,:,slc) .* ff(:,:,slc) .* hemo_core(:,:,slc)))];
                    sd_ff_array_hemo = [sd_ff_array_hemo, std(nonzeros(roi_in_myo_t2star(:,:,slc) .* ff(:,:,slc) .* hemo_core(:,:,slc)))];
                    
                    mean_ff_array_peri = [mean_ff_array_peri, mean(nonzeros(roi_in_myo_t2star(:,:,slc) .* ff(:,:,slc) .* ~hemo_core(:,:,slc)))];
                    sd_ff_array_peri = [sd_ff_array_peri, std(nonzeros(roi_in_myo_t2star(:,:,slc) .* ff(:,:,slc) .* ~hemo_core(:,:,slc)))];
                    
                else
                    mean_ff_array_nonhemo = [mean_ff_array_nonhemo, mean(nonzeros(roi_in_myo_t2star(:,:,slc) .* ff(:,:,slc)))];
                    sd_ff_array_nonhemo = [sd_ff_array_nonhemo, std(nonzeros(roi_in_myo_t2star(:,:,slc) .* ff(:,:,slc)))];
                    mean_ff_array_remote_nonhemo = [mean_ff_array_remote_nonhemo, mean(nonzeros(remote_in_myo_t2star(:,:,slc) .* ff(:,:,slc)))];
                    sd_ff_array_remote_nonhemo = [sd_ff_array_remote_nonhemo, std(nonzeros(remote_in_myo_t2star(:,:,slc) .* ff(:,:,slc)))];
                    %name_label_nonhemo{count_nonhemo} = cat(2, name, '_', time_point, '_Slice', num2str(slc));
                    %count_nonhemo = count_nonhemo + 1;
                end
                
                ff_perc = sum(sum(ff_core(:,:,slc))) / sum(roi_in_myo_t2star(:));
                if any(temp(:)) && ff_perc > 0.05
                    name_label_fat_hemo{end+1} = cat(2, name, '_', time_point, '_Slice', num2str(slc));
                elseif any(temp(:)) && ff_perc <= 0.05
                    name_label_nonfat_hemo{end+1} = cat(2, name, '_', time_point, '_Slice', num2str(slc));
                elseif ~any(temp(:)) && ff_perc > 0.05
                    name_label_fat_nonhemo{end+1} = cat(2, name, '_', time_point, '_Slice', num2str(slc));
                elseif ~any(temp(:)) && ff_perc <= 0.05
                    name_label_nonfat_nonhemo{end+1} = cat(2, name, '_', time_point, '_Slice', num2str(slc));
                end
                
                hemo_perc_slc_base = sum(sum(hemo_core(:,:,slc))) / sum(mask_myocardium_3D_eroded(:));
                hemo_perc_array(end+1) = hemo_perc_slc_base;
                ff_array(end+1) = ff_perc;
                
            end
            
        end
    end
    %close all;
end

%% 
figure();
x = [1, 2, 3, 4, 5, 6];
y = [mean(mean_R2star_array); nanmean(mean_R2star_array_hemo); nanmean(mean_R2star_array_peri); mean(mean_R2star_array_remote); mean(mean_R2star_array_nonhemo); mean(mean_R2star_array_remote_nonhemo)];
bar(x, y); ylim([0 100]);
xticklabels({'MI (Hemo+)', 'Hemo Core', 'Hemo Peri', 'Remote (Hemo+)', 'MI (Hemo-)', 'Remote (Hemo-)'});

err = [std(mean_R2star_array); nanstd(mean_R2star_array_hemo); nanstd(mean_R2star_array_peri); std(mean_R2star_array_remote); std(mean_R2star_array_nonhemo); std(mean_R2star_array_remote_nonhemo)];
hold on;
er = errorbar(x,y, err);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
ylabel('R2* (s^{-1})');

figure();
x = [1, 2, 3, 4, 5, 6];
y = [mean(mean_ff_array); nanmean(mean_ff_array_hemo); nanmean(mean_ff_array_peri); mean(mean_ff_array_remote); mean(mean_ff_array_nonhemo); mean(mean_ff_array_remote_nonhemo)];
bar(x, y); ylim([0 20]);
xticklabels({'MI (Hemo+)', 'Hemo Core', 'Hemo Peri', 'Remote (Hemo+)', 'MI (Hemo-)', 'Remote (Hemo-)'});

err = [std(mean_ff_array); nanstd(mean_ff_array_hemo); nanstd(mean_ff_array_peri); std(mean_ff_array_remote); std(mean_ff_array_nonhemo); std(mean_ff_array_remote_nonhemo)];
hold on;
er = errorbar(x,y, err);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
ylabel('FF (%)');

%%
figure();
plot(mean_R2star_array, mean_ff_array, 'o');
hold on;
plot(mean_R2star_array_nonhemo,mean_ff_array_nonhemo, 'o');
mdl = fitlm([mean_R2star_array,mean_R2star_array_nonhemo], [mean_ff_array,mean_ff_array_nonhemo]);
Y = [mean_R2star_array,mean_R2star_array_nonhemo] .* mdl.Coefficients.Estimate(2) + mdl.Coefficients.Estimate(1);
plot([mean_R2star_array,mean_R2star_array_nonhemo], Y, 'k', 'LineWidth', 1);
yl = ylim;
xl = xlim;
text(0.6*xl(2), yl(1)+5, cat(2,'Y = ', num2str(mdl.Coefficients.Estimate(2), 2), 'X + ', num2str(mdl.Coefficients.Estimate(1),2), ', R^2 = ', num2str(mdl.Rsquared.Ordinary,3)), 'FontSize', 12)
xlabel('R2star'); ylabel('FF');

%% 

figure();
plot(hemo_perc_array, ff_array, 'o');
xlabel('Hemo Percentage'); ylabel('Fat Percentage');

