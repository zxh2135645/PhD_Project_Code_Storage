close all;
clear all;
% T1 and FF analysis 12142020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% pre_QualControl.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% data_storage_rim.mat
% Time_Evolution_rim.png
% Time_Evolution.png
% for_analysis.mat
% for_analysis_rim.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

metrics_save_dir = cat(2, base_dir, '/Results/');
if ~exist(metrics_save_dir, 'dir')
   mkdir(metrics_save_dir); 
end
load(cat(2, metrics_save_dir, 'pre_QualControl.mat'));

%% Pull LGE, T1 Map, FF Map, and True R2* Map
% This can be skipped for the second time

label_t1 = sequence_label{1};
label_lge = sequence_label{3};
label_t2star = sequence_label{2};
for_analysis = struct;


%for n = 9:length(Names)
for n = 1:1
% for n = starting_point:starting_point
% Do not need to pull up images for baseline
    name = Names{n};
    name_save_dir = cat(2, save_dir, name);
    if ~exist(name_save_dir, 'dir')
        mkdir(name_save_dir);
    end
    
    % for tp = 1:length(time_points)
    for tp = 9:9
        time_point = time_points{end-tp+1};
        tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', name, '_', time_point,  '/');
        if ~exist(tp_dir, 'dir')
            disp(cat(2, 'No folder at: ', name, ' ', time_point));
        else
            % T1
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
            
            roi_in_myo_lge = mask_myocardium_3D .* freeROIMask_3D;
            remote_in_myo_lge = mask_myocardium_3D .* myoRefMask_3D;
            roi_lge = roi_in_myo_lge .* vol_img_3D;
            remote_lge = remote_in_myo_lge .* vol_img_3D;
            lge = vol_img_3D;
            myo_lge = mask_myocardium_3D;
            
            idx_reordered = Func_AlignSliceLoc(slc_array_t1, slc_array_lge);
            lge = lge(:,:,idx_reordered);
            myo_lge = myo_lge(:,:,idx_reordered);
            remote_lge = remote_lge(:,:,idx_reordered);
            roi_lge = roi_lge(:,:,idx_reordered);
            remote_in_myo_lge = remote_in_myo_lge(:,:,idx_reordered);
            roi_in_myo_lge = roi_in_myo_lge(:,:,idx_reordered);
            
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
            r2star = zeros(size(r2star_map{1}.fwmc_r2star,1), size(r2star_map{1}.fwmc_r2star, 2), length(r2star_map));
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
            
            for i = 1:size(roi_t1, 3)
                %for i = 1:1
                figure();
                subplot(2,2,1);
                imagesc(t1(:,:,i)); axis image; colorbar;
                title(['T1 Map: ', name, ' ', time_point, '  Slice = ', num2str(i)]);
                colormap(brewermap([],'*RdYlBu'));
                subplot(2,2,2);
                imagesc(lge(:,:,i)); axis image; colorbar;
                title(['LGE: ', name, ' ', time_point, '  Slice = ', num2str(i)]);
                
                subplot(2,2,3);
                imagesc(ff(:,:,i)); axis image; caxis([0 40]); colorbar;
                title(['FF Map: ', name, ' ', time_point, '  Slice = ', num2str(i)]);
                subplot(2,2,4);
                imagesc(r2star(:,:,i)); axis image; caxis([0 200]); colorbar;
                title(['R2star Map: ', name, ' ', time_point, '  Slice = ', num2str(i)]);
                
                overview_dir = cat(2, name_save_dir, '/overview/');
                if ~exist(overview_dir, 'dir')
                    mkdir(overview_dir);
                end
                saveas(gcf, cat(2, overview_dir, name, '_', time_point, '_Slice', num2str(i), '.png'));
            end

            for i = 1:size(roi_t1, 3)
            %for i = 1:1
                
                
                figure();
                subplot(2,2,1);
                roi_edg_t1  = edge(squeeze(roi_in_myo_t1(:,:,i)),'Canny');
                remote_edg_t1  = edge(squeeze(remote_in_myo_t1(:,:,i)),'Canny');
                myo_edg_t1 = edge(squeeze(myo_t1(:,:,i)), 'Canny');
                RGB_t1 = Func_Display_As_RGB(t1(:,:,i), roi_edg_t1, remote_edg_t1, myo_edg_t1);
                
                imagesc(RGB_t1); axis image;
                title(['T1 Map: ', name, ' ', time_point, '  Slice = ', num2str(i)]);
                %colormap(brewermap([],'*RdYlBu'));
                
                subplot(2,2,2);
                roi_edg_lge  = edge(squeeze(roi_in_myo_lge(:,:,i)),'Canny');
                remote_edg_lge  = edge(squeeze(remote_in_myo_lge(:,:,i)),'Canny');
                myo_edg_lge = edge(squeeze(myo_lge(:,:,i)), 'Canny');
                RGB_lge = Func_Display_As_RGB(lge(:,:,i), roi_edg_lge, remote_edg_lge, myo_edg_lge);
                
                imagesc(RGB_lge); axis image; 
                title(['LGE: ', name, ' ', time_point, '  Slice = ', num2str(i)]);
                
                subplot(2,2,3);
                roi_edg_ff  = edge(squeeze(roi_in_myo_ff(:,:,i)),'Canny');
                remote_edg_ff  = edge(squeeze(remote_in_myo_ff(:,:,i)),'Canny');
                myo_edg_ff = edge(squeeze(myo_ff(:,:,i)), 'Canny');
                ff_scaled = ff(:,:,i);
                ff_scaled(ff_scaled > 40) = 40;
                ff_scaled(ff_scaled < 0) = 0;
                RGB_ff = Func_Display_As_RGB(ff_scaled, roi_edg_ff, remote_edg_ff, myo_edg_ff);
                
                
                imagesc(RGB_ff); axis image;
                title(['FF Map: ', name, ' ', time_point, '  Slice = ', num2str(i)]);
                
                subplot(2,2,4);
                roi_edg_r2star  = edge(squeeze(roi_in_myo_r2star(:,:,i)),'Canny');
                remote_edg_r2star  = edge(squeeze(remote_in_myo_r2star(:,:,i)),'Canny');
                myo_edg_r2star = edge(squeeze(myo_r2star(:,:,i)), 'Canny');
                r2star_scaled = r2star(:,:,i);
                r2star_scaled(r2star_scaled > 200) = 200;
                r2star_scaled(r2star_scaled < 0) = 0;
                
                RGB_r2star = Func_Display_As_RGB(r2star_scaled, roi_edg_r2star, remote_edg_r2star, myo_edg_r2star);
                
                imagesc(RGB_r2star); axis image;
                title(['R2star Map: ', name, ' ', time_point, '  Slice = ', num2str(i)]);
                
                overview_dir = cat(2, name_save_dir, '/overview/');
                if ~exist(overview_dir, 'dir')
                    mkdir(overview_dir);
                end
                saveas(gcf, cat(2, overview_dir, name, '_', time_point, '_Contour_Slice', num2str(i), '.png'));
            end
            
            %saveas(gcf, cat(2, name_save_dir, '\', name, '_', time_point, '_MI.png'));
            close all;
            
            name_tp = cat(2, name, '_', time_point);
            for_analysis(n).Name = name;
            for_analysis(n).metrics = struct;
            for_analysis(n).metrics.time_point = time_point;
            for_analysis(n).metrics.mean_roi_t1 = mean(nonzeros(roi_t1));
            for_analysis(n).metrics.sd_roi_t1 = std(nonzeros(roi_t1));
            for_analysis(n).metrics.mean_remote_t1 = mean(nonzeros(remote_t1));
            for_analysis(n).metrics.sd_remote_t1 = std(nonzeros(remote_t1));
            
            for_analysis(n).metrics.mean_roi_ff = mean(nonzeros(roi_ff));
            for_analysis(n).metrics.sd_roi_ff = std(nonzeros(roi_ff));
            for_analysis(n).metrics.mean_remote_ff = mean(nonzeros(remote_ff));
            for_analysis(n).metrics.sd_remote_ff = std(nonzeros(remote_ff));
            
            for_analysis(n).metrics.mean_roi_r2star = mean(nonzeros(roi_r2star));
            for_analysis(n).metrics.sd_roi_r2star = std(nonzeros(roi_r2star));
            for_analysis(n).metrics.mean_remote_r2star = mean(nonzeros(remote_r2star));
            for_analysis(n).metrics.sd_remote_r2star = std(nonzeros(remote_r2star));
        end
        
        
%         figure('Position', [1000 500 900 600]);
%         errorbar(metrics, metrics_sd, '-o', 'LineWidth', 2); xlabel('Time point'); ylabel('mean T1 (ms)');
%         hold on;  errorbar(metrics_remote, metrics_sd_remote, '-o', 'LineWidth', 2);
%         xticklabels(time_points);
%         set(gca, 'FontSize', 16);
%         grid on;
%         legend({'MI', 'Remote'}, 'Location', 'SouthEast');
        
    end
end

%% Here is where analysis starts
% First looking at mean+sd value time evolution
%% The pre_QualControl is generated in T1FP_GenDict_FromQC.m
% after running the section above
%% Before analysis, parse pre_QualControl
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

%% Skip to main body
%% Get the matrics of T1, FF and R2star (Without peeling off edge pixels)
%for ll = 1:length(sequence_label)
% T1 Map
% Images will be stored at img/<name>/overview/
% To exclude certain slices that has bad image quality

label_t1 = sequence_label{1};
label_lge = sequence_label{3};
label_t2star = sequence_label{2};
for_analysis = struct;

for n = 1:length(Names)
%for n = starting_point:starting_point
    name = Names{n};
    name_save_dir = cat(2, save_dir, name);
    if ~exist(name_save_dir, 'dir')
        mkdir(name_save_dir);
    end
    for_analysis(n).Name = name;
    for_analysis(n).metrics = struct;
    tp_count = 1;
    for tp = 1:length(time_points)
        %for tp = 1:length(time_points)
        time_point = time_points{end-tp+1};
        tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', name, '_', time_point,  '/');
        if ~exist(tp_dir, 'dir')
            disp(cat(2, 'No folder at: ', name, ' ', time_point));
        else
            % T1
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
            
            % Pre-QC
            [roi_in_myo_t1_new, remote_in_myo_t1_new, roi_t1_new, remote_t1_new, t1_new, myo_t1_new] = ...
                Func_status_check(status_check, n, name, tp_count, roi_in_myo_t1, remote_in_myo_t1,...
                roi_t1, remote_t1, t1, myo_t1);

            [row_roi_t1, col_roi_t1, v_roi_t1] = find(roi_in_myo_t1_new);
            [row_remote_t1, col_remote_t1, v_remote_t1] = find(remote_in_myo_t1_new);
            
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
            ff = zeros(size(ff_map{1}.fwmc_ff,1), size(ff_map{2}.fwmc_ff, 2), length(ff_map));
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
            
            % Pre-QC
            [roi_in_myo_ff_new, remote_in_myo_ff_new, roi_ff_new, remote_ff_new, ff_new, myo_ff_new] = ...
                Func_status_check(status_check, n, name, tp_count, roi_in_myo_ff, remote_in_myo_ff,...
                roi_ff, remote_ff, ff, myo_ff);
            
            [row_roi, col_roi, v_roi] = find(roi_in_myo_ff_new);
            [row_remote, col_remote, v_remote] = find(remote_in_myo_ff_new);
            
            % R2star Map
            r2star_map = cell(1, length(glob_names));
            for f = 1:length(r2star_map)
                r2star_map{f} = load(cat(2,  base_dir, '/FF_Data/',  name, '/', name, '_', time_point, '_', glob_names{f}, '.mat'), 'fwmc_r2star');
            end
            
            % convert ff_map to matrix
            r2star = zeros(size(r2star_map{1}.fwmc_r2star,1), size(r2star_map{1}.fwmc_r2star, 2), length(r2star_map));
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
            
            % Pre-QC
            [roi_in_myo_r2star_new, remote_in_myo_r2star_new, roi_r2star_new, remote_r2star_new, r2star_new, myo_r2star_new] = ...
                Func_status_check(status_check, n, name, tp_count, roi_in_myo_r2star, remote_in_myo_r2star,...
                roi_r2star, remote_r2star, r2star, myo_r2star);
            
            [row_roi_r2star, col_roi_r2star, v_roi_r2star] = find(roi_in_myo_r2star_new);
            [row_remote_r2star, col_remote_r2star, v_remote_r2star] = find(remote_in_myo_r2star_new);
            
            
            roi_ff_array = zeros(1, length(row_roi));
            remote_ff_array = zeros(1, length(row_remote));
            roi_r2star_array = zeros(1, length(row_roi_r2star));
            remote_r2star_array = zeros(1, length(row_remote_r2star));
            roi_t1_array = zeros(1, length(row_roi_t1));
            remote_t1_array = zeros(1, length(row_remote_t1));
            
            for fff = 1:length(row_roi)
               roi_ff_array(fff) = roi_ff_new(row_roi(fff), col_roi(fff));
            end
            
            for fff = 1:length(row_remote)
                remote_ff_array(fff) = remote_ff_new(row_remote(fff), col_remote(fff));
            end
            
            for fff = 1:length(row_roi_r2star)
                roi_r2star_array(fff) = roi_r2star_new(row_roi_r2star(fff), col_roi_r2star(fff));
            end
            
            for fff = 1:length(row_remote_r2star)
                remote_r2star_array(fff) = remote_r2star_new(row_remote_r2star(fff), col_remote_r2star(fff));
            end
            
            for fff = 1:length(row_roi_t1)
                roi_t1_array(fff) = roi_t1_new(row_roi_t1(fff), col_roi_t1(fff));
            end
            
            for fff = 1:length(row_remote_t1)
                remote_t1_array(fff) = remote_t1_new(row_remote_t1(fff), col_remote_t1(fff));
            end
            
            roi_ff_array(roi_ff_array < 0) = 0;
            roi_ff_array(roi_ff_array > 100) = 100;
            remote_ff_array(remote_ff_array < 0) = 0;
            remote_ff_array(remote_ff_array > 100) = 100;
            
            roi_r2star_array(roi_r2star_array > 100) = 100;
            remote_r2star_array(remote_r2star_array > 100) = 100;
                        
            % Plot Histogram of ROI vs Remote (Gross view)
            figure();
            h_ff = histogram(roi_ff_array, 'Normalization', 'probability');xlabel('Fat Fraction (%)'); ylabel('Frequency');
            NumBins = h_ff.NumBins;
            BinWidth = h_ff.BinWidth;
            BinEdges = h_ff.BinEdges;
            hold on;
            h_remote = histogram(remote_ff_array, 'Normalization', 'probability');
            h_remote.BinEdges = BinEdges;
            set(gca, 'FontSize', 16); title([name, ' ', time_point]);
            xlim([0 100]);
            legend({'MI', 'Remote'});
            name_tp_dir = cat(2, name_save_dir, '/', name, '_', time_point);
            if ~exist(name_tp_dir, 'dir')
                mkdir(name_tp_dir);
            end
            saveas(gcf, cat(2, name_save_dir, '/', name, '_', time_point, '/Histogram_FF_Whole.png'));
            
            % R2star
            figure();
            h_r2star = histogram(roi_r2star_array, 'Normalization', 'probability');xlabel('R2* (Hz)'); ylabel('Frequency');
            NumBins = h_r2star.NumBins;
            BinWidth = h_r2star.BinWidth;
            BinEdges = h_r2star.BinEdges;
            hold on;
            h_r2star_remote = histogram(remote_r2star_array, 'Normalization', 'probability');
            h_r2star_remote.BinEdges = BinEdges;
            set(gca, 'FontSize', 16); title([name, ' ', time_point]);
            xlim([0 100]);
            legend({'MI', 'Remote'});
            name_tp_dir = cat(2, name_save_dir, '/', name, '_', time_point);
            if ~exist(name_tp_dir, 'dir')
                mkdir(name_tp_dir);
            end
            saveas(gcf, cat(2, name_save_dir, '/', name, '_', time_point, '/Histogram_R2star_Whole.png'));
            
            % T1
            figure();
            h_t1 = histogram(roi_t1_array, 'Normalization', 'probability');xlabel('T1 (ms)'); ylabel('Frequency');
            NumBins = h_t1.NumBins;
            BinWidth = h_t1.BinWidth;
            BinEdges = h_t1.BinEdges;
            hold on;
            h_t1_remote = histogram(remote_t1_array, 'Normalization', 'probability');
            h_t1_remote.BinEdges = BinEdges;
            set(gca, 'FontSize', 16); title([name, ' ', time_point]);
            %xlim([0 100]);
            legend({'MI', 'Remote'});

            name_tp_dir = cat(2, name_save_dir, '/', name, '_', time_point);
            if ~exist(name_tp_dir, 'dir')
                mkdir(name_tp_dir);
            end
            saveas(gcf, cat(2, name_save_dir, '/', name, '_', time_point, '/Histogram_T1_Whole.png'));
            
            % Histogram and curve, not necessary
%             figure();
%             hh_ff = histfit(roi_ff_array, 20, 'exponential');
%             yt = get(gca, 'YTick');
%             set(gca, 'YTick', yt, 'YTickLabel', round(yt/numel(roi_ff_array), 2))
%             left_yt = round(yt/numel(roi_ff_array), 2);
%             hold on;
%             yyaxis right; histfit(remote_ff_array, 20, 'exponential');
%             yt = get(gca, 'YTick');
%             set(gca, 'YTick', yt, 'YTickLabel', round(yt/numel(remote_ff_array), 2))

            name_tp = cat(2, name, '_', time_point);
            
            for_analysis(n).metrics(tp).time_point = time_point;
            for_analysis(n).metrics(tp).mean_roi_t1 = mean(nonzeros(roi_t1_array));
            for_analysis(n).metrics(tp).sd_roi_t1 = std(nonzeros(roi_t1_array));
            for_analysis(n).metrics(tp).mean_remote_t1 = mean(nonzeros(remote_t1_array));
            for_analysis(n).metrics(tp).sd_remote_t1 = std(nonzeros(remote_t1_array));
            
            for_analysis(n).metrics(tp).mean_roi_ff = mean(roi_ff_array);
            for_analysis(n).metrics(tp).sd_roi_ff = std(roi_ff_array);
            for_analysis(n).metrics(tp).mean_remote_ff = mean(remote_ff_array);
            for_analysis(n).metrics(tp).sd_remote_ff = std(remote_ff_array);
            
            for_analysis(n).metrics(tp).mean_roi_r2star = mean(nonzeros(roi_r2star_array));
            for_analysis(n).metrics(tp).sd_roi_r2star = std(nonzeros(roi_r2star_array));
            for_analysis(n).metrics(tp).mean_remote_r2star = mean(nonzeros(remote_r2star_array));
            for_analysis(n).metrics(tp).sd_remote_r2star = std(nonzeros(remote_r2star_array));
            
            close all;
            %break;
            tp_count = tp_count + 1;
        end
        
        
%                     figure('Position', [1000 500 900 600]);
%                     errorbar(metrics, metrics_sd, '-o', 'LineWidth', 2); xlabel('Time point'); ylabel('mean T1 (ms)');
%                     hold on;  errorbar(metrics_remote, metrics_sd_remote, '-o', 'LineWidth', 2);
%                     xticklabels(time_points);
%                     set(gca, 'FontSize', 16);
%                     grid on;
%                     legend({'MI', 'Remote'}, 'Location', 'SouthEast');

    end
    close all;
end
%% The following needs to be deprecated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Longitudinal time evolution
time_points_lr = fliplr(time_points);
xtcks = 1:length(time_points);
xtcks_lr = fliplr(xtcks);
for n = 1:length(Names)
%for n = starting_point:starting_point
    name = Names{n};
    name_save_dir = cat(2, save_dir, name);
    long_name_save_dir = cat(2, name_save_dir, '/Longitudinal/');
    if ~exist(long_name_save_dir, 'dir')
        mkdir(long_name_save_dir);
    end
    
    mean_roi_t1 = zeros(1, length(time_points));
    sd_roi_t1 = zeros(1, length(time_points));
    mean_remote_t1 = zeros(1, length(time_points));
    sd_remote_t1 = zeros(1, length(time_points));
    
    mean_roi_ff = zeros(1, length(time_points));
    sd_roi_ff = zeros(1, length(time_points));
    mean_remote_ff = zeros(1, length(time_points));
    sd_remote_ff = zeros(1, length(time_points));
    
    mean_roi_r2star = zeros(1, length(time_points));
    sd_roi_r2star = zeros(1, length(time_points));
    mean_remote_r2star = zeros(1, length(time_points));
    sd_remote_r2star = zeros(1, length(time_points));
    
    for tp = 1:length(time_points)
        
        roi_t1_temp = for_analysis(n).metrics(tp).mean_roi_t1;
        if isempty(roi_t1_temp)
            mean_roi_t1(tp) = NaN;
            sd_roi_t1(tp) = NaN;
            mean_remote_t1(tp) = NaN;
            sd_remote_t1(tp) = NaN;
            
            mean_roi_ff(tp) = NaN;
            sd_roi_ff(tp) = NaN;
            mean_remote_ff(tp) = NaN;
            sd_remote_ff(tp) = NaN;
            
            mean_roi_r2star(tp) = NaN;
            sd_roi_r2star(tp) = NaN;
            mean_remote_r2star(tp) = NaN;
            sd_remote_r2star(tp) = NaN;
            
        else
            sd_roi_t1(tp) = for_analysis(n).metrics(tp).sd_roi_t1;
            mean_roi_t1(tp) = for_analysis(n).metrics(tp).mean_roi_t1;
            mean_remote_t1(tp) = for_analysis(n).metrics(tp).mean_remote_t1;
            sd_remote_t1(tp) = for_analysis(n).metrics(tp).sd_remote_t1;
            
            mean_roi_ff(tp) = for_analysis(n).metrics(tp).mean_roi_ff;
            sd_roi_ff(tp) = for_analysis(n).metrics(tp).sd_roi_ff;
            mean_remote_ff(tp) = for_analysis(n).metrics(tp).mean_remote_ff;
            sd_remote_ff(tp) = for_analysis(n).metrics(tp).sd_remote_ff;
            
            mean_roi_r2star(tp) = for_analysis(n).metrics(tp).mean_roi_r2star;
            sd_roi_r2star(tp) = for_analysis(n).metrics(tp).sd_roi_r2star;
            mean_remote_r2star(tp) = for_analysis(n).metrics(tp).mean_remote_r2star;
            sd_remote_r2star(tp) = for_analysis(n).metrics(tp).sd_remote_r2star;
        end
    end
    
    I = ~isnan(mean_roi_t1);
    
    figure('Position', [1000 500 900 600]);
    subplot(3,1,1);
    errorbar(xtcks_lr(I), mean_roi_t1(I), sd_roi_t1(I), '-o', 'LineWidth', 2); xlabel('Time point'); ylabel('mean T1 (ms)');
    hold on;  errorbar(xtcks_lr(I), mean_remote_t1(I), sd_remote_t1(I), '-o', 'LineWidth', 2);
    xticks(fliplr(xtcks_lr(I)));
    xticklabels(fliplr(time_points_lr(I)));
    set(gca, 'FontSize', 16);
    grid on;
    ylim([800 1600]); xlim([xtcks(1) xtcks(end)]);
    legend({'MI', 'Remote'}, 'Location', 'SouthEast');
    title(cat(2, name, '  T1 evolution'));
    
    
    subplot(3,1,2);
    errorbar(xtcks_lr(I), mean_roi_ff(I), sd_roi_ff(I), '-o', 'LineWidth', 2); xlabel('Time point'); ylabel('Fat Fraction (%)');
    hold on;  errorbar(xtcks_lr(I), mean_remote_ff(I), sd_remote_ff(I), '-o', 'LineWidth', 2);
    xticks(fliplr(xtcks_lr(I)));
    xticklabels(fliplr(time_points_lr(I)));
    set(gca, 'FontSize', 16);
    grid on;
    ylim([-5 20]); xlim([xtcks(1) xtcks(end)]);
    legend({'MI', 'Remote'}, 'Location', 'SouthEast');
    title(cat(2, name, '  FF evolution'));
    
    subplot(3,1,3);
    errorbar(xtcks_lr(I), mean_roi_r2star(I), sd_roi_r2star(I), '-o', 'LineWidth', 2); xlabel('Time point'); ylabel('R2 star (Hz)');
    hold on;  errorbar(xtcks_lr(I), mean_remote_r2star(I), sd_remote_r2star(I), '-o', 'LineWidth', 2);
    xticks(fliplr(xtcks_lr(I)));
    xticklabels(fliplr(time_points_lr(I)));
    set(gca, 'FontSize', 16);
    grid on;
    xlim([xtcks(1) xtcks(end)]);
    %ylim([-30 30])
    legend({'MI', 'Remote'}, 'Location', 'SouthEast');
    title(cat(2, name, '  R2 star evolution'));
    
    saveas(gcf, cat(2, long_name_save_dir, '/Time_Evolution.png'));
end

close all;
%% The below is main body of this script
% Because it peels off the edge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get the matrics of T1, FF and R2star ( exclude edges of MI region)
%for ll = 1:length(sequence_label)
% T1 Map
% Images will be stored at img/<name>/overview/
% To exclude certain slices that has bad image quality

label_t1 = sequence_label{1};
label_lge = sequence_label{3};
label_t2star = sequence_label{2};
for_analysis_rim = struct;
data_storage_rim = struct;

%time_points = {'6MO', '9MO', '1YR', '15YR'};

for n = 1:length(Names)
%for n = 2:2
    name = Names{n};
    name_save_dir = cat(2, save_dir, name);
    if ~exist(name_save_dir, 'dir')
        mkdir(name_save_dir);
    end
    for_analysis_rim(n).Name = name;
    for_analysis_rim(n).metrics = struct;
    tp_count = 1;
    data_storage_rim(n).Name = name;
    data_storage_rim(n).data = struct;
        
    for tp = 1:length(time_points)
    %for tp = 2:2
        time_point = time_points{end-tp+1};
        tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', name, '_', time_point,  '/');
        if ~exist(tp_dir, 'dir')
            disp(cat(2, 'No folder at: ', name, ' ', time_point));
        else
            % T1
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
            
            % Pre-QC
            [roi_in_myo_t1_new, remote_in_myo_t1_new, roi_t1_new, remote_t1_new, t1_new, myo_t1_new] = ...
                Func_status_check(status_check, n, name, tp_count, roi_in_myo_t1, remote_in_myo_t1,...
                roi_t1, remote_t1, t1, myo_t1);
            
            % remove edges of MI region
            roi_edg_t1_new = zeros(size(roi_in_myo_t1_new));
            for i = 1:size(roi_in_myo_t1_new, 3)
                roi_edg_t1_new(:,:,i)  = edge(squeeze(roi_in_myo_t1_new(:,:,i)),'Canny');
            end
            roi_rimmed_t1_new = (roi_in_myo_t1_new - roi_edg_t1_new)>0;
            [row_roi_t1, col_roi_t1, v_roi_t1] = find(roi_rimmed_t1_new);
            [row_remote_t1, col_remote_t1, v_remote_t1] = find(remote_in_myo_t1_new);
            
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
            ff = zeros(size(ff_map{1}.fwmc_ff,1), size(ff_map{2}.fwmc_ff, 2), length(ff_map));
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
            
            % Pre-QC
            [roi_in_myo_ff_new, remote_in_myo_ff_new, roi_ff_new, remote_ff_new, ff_new, myo_ff_new] = ...
                Func_status_check(status_check, n, name, tp_count, roi_in_myo_ff, remote_in_myo_ff,...
                roi_ff, remote_ff, ff, myo_ff);
            
            % remove edges of MI region
            roi_edg_ff_new = zeros(size(roi_in_myo_ff_new));
            for i = 1:size(roi_in_myo_ff_new, 3)
                roi_edg_ff_new(:,:,i)  = edge(squeeze(roi_in_myo_ff_new(:,:,i)),'Canny');
            end
            roi_rimmed_ff_new = (roi_in_myo_ff_new - roi_edg_ff_new)>0;
            roi_ff_new = roi_ff_new .* roi_rimmed_ff_new;


            
            [row_roi, col_roi, v_roi] = find(roi_rimmed_ff_new);
            [row_remote, col_remote, v_remote] = find(remote_in_myo_ff_new);
            
            % R2star Map
            r2star_map = cell(1, length(glob_names));
            for f = 1:length(r2star_map)
                r2star_map{f} = load(cat(2,  base_dir, '/FF_Data/',  name, '/', name, '_', time_point, '_', glob_names{f}, '.mat'), 'fwmc_r2star');
            end
            
            % convert ff_map to matrix
            r2star = zeros(size(r2star_map{1}.fwmc_r2star,1), size(r2star_map{1}.fwmc_r2star, 2), length(r2star_map));
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
            
            % Pre-QC
            [roi_in_myo_r2star_new, remote_in_myo_r2star_new, roi_r2star_new, remote_r2star_new, r2star_new, myo_r2star_new] = ...
                Func_status_check(status_check, n, name, tp_count, roi_in_myo_r2star, remote_in_myo_r2star,...
                roi_r2star, remote_r2star, r2star, myo_r2star);
            
            % remove edges of MI region
            roi_edg_r2star_new = zeros(size(roi_in_myo_r2star_new));
            for i = 1:size(roi_in_myo_r2star_new, 3)
                roi_edg_r2star_new(:,:,i)  = edge(squeeze(roi_in_myo_r2star_new(:,:,i)),'Canny');
            end
            roi_rimmed_r2star_new = (roi_in_myo_r2star_new - roi_edg_r2star_new)>0;
            roi_r2star_new = roi_rimmed_r2star_new .* roi_r2star_new;

            [row_roi_r2star, col_roi_r2star, v_roi_r2star] = find(roi_rimmed_r2star_new);
            [row_remote_r2star, col_remote_r2star, v_remote_r2star] = find(remote_in_myo_r2star_new);

           
            
            roi_ff_array = zeros(1, length(row_roi));
            remote_ff_array = zeros(1, length(row_remote));
            roi_r2star_array = zeros(1, length(row_roi_r2star));
            remote_r2star_array = zeros(1, length(row_remote_r2star));
            roi_t1_array = zeros(1, length(row_roi_t1));
            remote_t1_array = zeros(1, length(row_remote_t1));
            
            for fff = 1:length(row_roi)
               roi_ff_array(fff) = roi_ff_new(row_roi(fff), col_roi(fff));
            end
            
            for fff = 1:length(row_remote)
                remote_ff_array(fff) = remote_ff_new(row_remote(fff), col_remote(fff));
            end
            
            for fff = 1:length(row_roi_r2star)
                roi_r2star_array(fff) = roi_r2star_new(row_roi_r2star(fff), col_roi_r2star(fff));
            end
            
            for fff = 1:length(row_remote_r2star)
                remote_r2star_array(fff) = remote_r2star_new(row_remote_r2star(fff), col_remote_r2star(fff));
            end
            
            for fff = 1:length(row_roi_t1)
                roi_t1_array(fff) = roi_t1_new(row_roi_t1(fff), col_roi_t1(fff));
            end
            
            for fff = 1:length(row_remote_t1)
                remote_t1_array(fff) = remote_t1_new(row_remote_t1(fff), col_remote_t1(fff));
            end
            
            roi_ff_array(roi_ff_array < 0) = 0;
            roi_ff_array(roi_ff_array > 100) = 100;
            remote_ff_array(remote_ff_array < 0) = 0;
            remote_ff_array(remote_ff_array > 100) = 100;
            
            roi_r2star_array(roi_r2star_array > 100) = 100;
            remote_r2star_array(remote_r2star_array > 100) = 100;

            
             % XZ 09/20/2023
            roi_ff_new(roi_ff_new < 0) = 0;
            roi_ff_new(roi_ff_new > 100) = 100;
            roi_r2star_new(roi_r2star_new > 100) = 100;

            thresh = mean(remote_r2star_array) + 2*std(remote_r2star_array);
            fib_perc = zeros(1, size(roi_ff_new, 3));
            mi_perc = zeros(1, size(roi_ff_new, 3));
            mi_pix = zeros(1, size(roi_ff_new, 3));
            ff_pix = zeros(1, size(roi_ff_new, 3));
            r2star_pix = zeros(1, size(roi_ff_new, 3));
            union_pix = zeros(1, size(roi_ff_new, 3));
            intercept_pix = zeros(1, size(roi_ff_new, 3));

            for xx = 1:size(roi_ff_new, 3)
                roi_ff_thresh = roi_ff_new(:,:,xx) > 6;
                ff_pix(xx) = sum(sum(roi_ff_thresh));

                roi_r2star_thresh = roi_r2star_new(:,:,xx) > thresh;
                r2star_pix(xx) = sum(sum(roi_r2star_thresh));

                roi_union = roi_ff_thresh | roi_r2star_thresh;
                union_pix(xx) = sum(sum(roi_union));

                roi_intercept = roi_ff_thresh & roi_r2star_thresh;
                intercept_pix(xx) = sum(sum(roi_intercept));

                fib_perc(xx) = (1 - sum(roi_union(:)) / sum(sum(roi_rimmed_ff_new(:,:,xx))))*100;

                mi_perc(xx) = sum(sum(roi_in_myo_ff(:,:,xx))) / sum(sum(myo_ff_new(:,:,xx)));

                mi_pix(xx) = sum(sum(roi_in_myo_ff(:,:,xx)));
            end


            % Plot Histogram of ROI vs Remote (Gross view)
            figure();
            h_ff = histogram(roi_ff_array, 'Normalization', 'probability');xlabel('Fat Fraction (%)'); ylabel('Frequency');
            NumBins = h_ff.NumBins;
            BinWidth = h_ff.BinWidth;
            BinEdges = h_ff.BinEdges;
            hold on;
            h_remote = histogram(remote_ff_array, 'Normalization', 'probability');
            h_remote.BinEdges = BinEdges;
            set(gca, 'FontSize', 16); title([name, ' ', time_point]);
            xlim([0 100]);
            legend({'MI', 'Remote'});
            name_tp_dir = cat(2, name_save_dir, '/', name, '_', time_point);
            if ~exist(name_tp_dir, 'dir')
                mkdir(name_tp_dir);
            end
            saveas(gcf, cat(2, name_save_dir, '/', name, '_', time_point, '/Histogram_FF_Whole_rim.png'));
            
            % R2star
            figure();
            h_r2star = histogram(roi_r2star_array, 'Normalization', 'probability');xlabel('R2* (Hz)'); ylabel('Frequency');
            NumBins = h_r2star.NumBins;
            BinWidth = h_r2star.BinWidth;
            BinEdges = h_r2star.BinEdges;
            hold on;
            h_r2star_remote = histogram(remote_r2star_array, 'Normalization', 'probability');
            h_r2star_remote.BinEdges = BinEdges;
            set(gca, 'FontSize', 16); title([name, ' ', time_point]);
            xlim([0 100]);
            legend({'MI', 'Remote'});
            name_tp_dir = cat(2, name_save_dir, '/', name, '_', time_point);
            if ~exist(name_tp_dir, 'dir')
                mkdir(name_tp_dir);
            end
            saveas(gcf, cat(2, name_save_dir, '/', name, '_', time_point, '/Histogram_R2star_Whole_rim.png'));
            
            % T1
            figure();
            h_t1 = histogram(roi_t1_array, 'Normalization', 'probability');xlabel('T1 (ms)'); ylabel('Frequency');
            NumBins = h_t1.NumBins;
            BinWidth = h_t1.BinWidth;
            BinEdges = h_t1.BinEdges;
            hold on;
            h_t1_remote = histogram(remote_t1_array, 'Normalization', 'probability');
            h_t1_remote.BinEdges = BinEdges;
            set(gca, 'FontSize', 16); title([name, ' ', time_point]);
            %xlim([0 100]);
            legend({'MI', 'Remote'});

            name_tp_dir = cat(2, name_save_dir, '/', name, '_', time_point);
            if ~exist(name_tp_dir, 'dir')
                mkdir(name_tp_dir);
            end
            saveas(gcf, cat(2, name_save_dir, '/', name, '_', time_point, '/Histogram_T1_Whole_rim.png'));
            
            
            name_tp = cat(2, name, '_', time_point);
            
            for_analysis_rim(n).metrics(tp).time_point = time_point;
            for_analysis_rim(n).metrics(tp).mean_roi_t1 = mean(nonzeros(roi_t1_array));
            for_analysis_rim(n).metrics(tp).sd_roi_t1 = std(nonzeros(roi_t1_array));
            for_analysis_rim(n).metrics(tp).mean_remote_t1 = mean(nonzeros(remote_t1_array));
            for_analysis_rim(n).metrics(tp).sd_remote_t1 = std(nonzeros(remote_t1_array));
            
            for_analysis_rim(n).metrics(tp).mean_roi_ff = mean(roi_ff_array);
            for_analysis_rim(n).metrics(tp).sd_roi_ff = std(roi_ff_array);
            for_analysis_rim(n).metrics(tp).mean_remote_ff = mean(remote_ff_array);
            for_analysis_rim(n).metrics(tp).sd_remote_ff = std(remote_ff_array);
            
            for_analysis_rim(n).metrics(tp).mean_roi_r2star = mean(nonzeros(roi_r2star_array));
            for_analysis_rim(n).metrics(tp).sd_roi_r2star = std(nonzeros(roi_r2star_array));
            for_analysis_rim(n).metrics(tp).mean_remote_r2star = mean(nonzeros(remote_r2star_array));
            for_analysis_rim(n).metrics(tp).sd_remote_r2star = std(nonzeros(remote_r2star_array));

                        
            data_storage_rim(n).data(tp).time_point = time_point;
            data_storage_rim(n).data(tp).roi_ff_array = roi_ff_array;
            data_storage_rim(n).data(tp).remote_ff_array = remote_ff_array;
            data_storage_rim(n).data(tp).roi_r2star_array = roi_r2star_array;
            data_storage_rim(n).data(tp).remote_r2star_array = remote_r2star_array;
            data_storage_rim(n).data(tp).roi_t1_array = roi_t1_array;
            data_storage_rim(n).data(tp).remote_t1_array = remote_t1_array;

            data_storage_rim(n).data(tp).fib_perc = fib_perc;
            data_storage_rim(n).data(tp).mi_perc = mi_perc;
            data_storage_rim(n).data(tp).mi_pix = mi_pix;
            
            data_storage_rim(n).data(tp).ff_pix = ff_pix;
            data_storage_rim(n).data(tp).r2star_pix = r2star_pix;
            data_storage_rim(n).data(tp).union_pix = union_pix;
            data_storage_rim(n).data(tp).intercept_pix = intercept_pix;

            tp_count = tp_count + 1;
        end
    end
    close all;
end

%%
metrics_save_dir = cat(2, base_dir, '/Results/');
if ~exist(metrics_save_dir, 'dir')
   mkdir(metrics_save_dir); 
end
save(cat(2, metrics_save_dir, 'data_storage_rim.mat'), 'data_storage_rim');


%%

%% T2
output_label = {'LGE', 'T2star'};
save_dir = GetFullPath(cat(2, base_dir, '/Analysis/'));
data_save_dir = cat(2, base_dir, '/data/');

% Names = {'Merry', 'Ryn', 'Mojave', 'Sahara', 'ZZ', 'Tina', 'Sunny', 'Queenie', 'Hope', 'Gobi', 'Felicity', 'Evelyn', '18D15', '18D16', '11D05', '11D26', '11D33'};
%time_points = {'0D_baseline','1D', '7D', '28D', '8WK', '6MO', '9MO', '1YR', '15YR'};
time_points = {'7D', '8WK', '12WK', '14WK', '4MO', '6MO', '9MO', '1YR', '15YR'};

label_lge = sequence_label{1};
label_t1 = sequence_label{2};
label_t2 = 'T2';
label_mag = 'MAG';
label_psir = 'PSIR';
label_t1molli = 'T1MOLLI';

metrics_save_dir = cat(2, base_dir, '/Results/');
if ~exist(metrics_save_dir, 'dir')
   mkdir(metrics_save_dir); 
end

se = strel('disk', 1);
metrics_t2 = struct;
L_cell = {};
count = 1;

min_roi_value = 50;
max_roi_value = 50;

for n = 3:3
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

    metrics_t2(n).name = name;
    metrics_t2(n).TimePoints = struct;

    %for tp = 1:length(time_points)
    for tp = 9:9
        time_point = time_points{end-tp+1};
        tp_dir = cat(2, base_dir, '/ContourData/', name, '/',  name, '_', time_point,  '/');
        if ~exist(cat(2, tp_dir, label_t1, '/'), 'dir')
            disp(cat(2, 'No folder at: ', name, ' ', time_point));
        elseif (strcmp(name, '18D16') && strcmp(time_point, '9MO'))
            % skip
        else

            tp_dir2 = cat(2, name_save_dir, '/', name, '/', time_point, '/');
            if ~exist(tp_dir2, 'dir')
                mkdir(tp_dir2);
            end

            tp_count = tp_count+1;

                % T2
                myo_glob = glob(cat(2, tp_dir, label_t2, '/', anatomy_label{5}, '/*'));
                roi_glob = glob(cat(2, tp_dir, label_t2, '/',anatomy_label{3}, '/*'));
                remote_glob = glob(cat(2, tp_dir, label_t2, '/',anatomy_label{6}, '/*'));

                load(cat(2, tp_dir, label_t2, '/', label_t2, '_vol_img_3D.mat'));
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
                t2 = vol_img_3D * 0.1; % SCALE FACTOR
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
                    min_temp = min(min((nonzeros(temp_roi(:,:,slc)))));
                    max_temp = max(max((nonzeros(temp_roi(:,:,slc)))));

                    if min_temp < min_roi_value
                        min_roi_value = min_temp;
                    end

                    if max_temp > max_roi_value
                        max_roi_value = max_temp;
                    end
                end

                for slc = 1:size(roi_in_myo_t2,3)
                    %temp_norm_roi = temp_roi(:,:,slc);
                    %temp_norm_roi(temp_norm_roi == 0) = nan;
                    % temp_norm_roi = uint8((temp_roi(:,:,slc) - min(min((nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc))))))./ (max(max(nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc))))-min(min((nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc)))))) * 256);
                    % temp_norm_remote = uint8((temp_remote(:,:,slc) - min(min((nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc))))))./ (max(max(nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc))))-min(min((nonzeros(temp_roi(:,:,slc)+temp_remote(:,:,slc)))))) * 256);

                    % Hard-coded for Patient data T1 mapping
                    temp_norm_roi = temp_roi(:,:,slc);
                    temp_norm_remote = temp_remote(:,:,slc);

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


                    bipolar = temp_roi(:,:,slc) - mean(nonzeros(temp_remote(:,:,slc)));
                    weighted_map = double(bipolar < 0) .* 2 .* roi_in_myo_t2(:,:,slc) .* bipolar + double(bipolar >= 0) .* roi_in_myo_t2(:,:,slc) .* bipolar;
                    
                    lb = 2*(10 - mean(nonzeros(temp_remote(:,:,slc))));
                    ub = 60 - mean(nonzeros(temp_remote(:,:,slc)));
                    weighted_map(weighted_map < lb) = lb;
                    weighted_map(weighted_map > ub) = ub;

                    %temp_norm_roi = uint8((temp_roi(:,:,slc) - 800)./ (1800-800) * 256);
                    %temp_norm_remote = uint8((temp_remote(:,:,slc) - 800)./ (1800-800) * 256);
      
                    bipolar = temp_remote(:,:,slc) - mean(nonzeros(temp_remote(:,:,slc)));
                    weighted_map_remote = double(bipolar < 0) .* 2 .* roi_in_myo_t2(:,:,slc) .* bipolar + double(bipolar >= 0) .* roi_in_myo_t2(:,:,slc) .* bipolar;
                    
                    weighted_map_remote(weighted_map_remote < lb) = lb;
                    weighted_map_remote(weighted_map_remote > ub) = ub;

                    temp_norm_roi = uint8((weighted_map - lb)./ (ub-lb) * 256);
                    % temp_norm_remote = uint8((temp_remote(:,:,slc) - 800)./ (1800-800) * 256);

                    non_roi_value = uint8((0 - lb) ./ (ub - lb) * 256);

                    %temp_norm_roi = uint8((temp_roi(:,:,slc) - 0) ./ (60 - 0) * 256);
                    % temp_norm_remote = uint8((temp_remote(:,:,slc) - 0) ./ (60 - 10) * 256);

                    temp_norm_remote = uint8((weighted_map_remote - lb) ./ (ub - lb) * 256);

                    t1_norm_roi_slc = temp_norm_roi;
                    t1_norm_roi_slc(t1_norm_roi_slc == non_roi_value) = [];
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


                    non_remote_value = uint8((0 - lb) ./ (ub - lb) * 256);

                    t1_norm_remote_slc = temp_norm_remote;
                    t1_norm_remote_slc(t1_norm_remote_slc == non_remote_value) = [];
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

                    L = bwlabel(roi_in_myo_t2(:,:,slc));
                    if length(unique(L)) > 2
                        L_cell{count} = cat(2, name, ' ', time_point, ' Slice ', num2str(slc));
                        count = count + 1;
                    end
                end

        end
    end
end
%% Plot (Longitudinal plot)
time_points_lr = fliplr(time_points);
xtcks = 1:length(time_points);
xtcks_lr = fliplr(xtcks);
for n = 1:length(Names)
%for n = starting_point:starting_point
    name = Names{n};
    name_save_dir = cat(2, save_dir, name);
    long_name_save_dir = cat(2, name_save_dir, '/Longitudinal/');
    if ~exist(long_name_save_dir, 'dir')
        mkdir(long_name_save_dir);
    end
    
    mean_roi_t1 = zeros(1, size(data_storage_rim(n).data,2));
    sd_roi_t1 = zeros(1, size(data_storage_rim(n).data,2));
    mean_remote_t1 = zeros(1, size(data_storage_rim(n).data,2));
    sd_remote_t1 = zeros(1, size(data_storage_rim(n).data,2));
    
    mean_roi_ff = zeros(1, size(data_storage_rim(n).data,2));
    sd_roi_ff = zeros(1, size(data_storage_rim(n).data,2));
    mean_remote_ff = zeros(1, size(data_storage_rim(n).data,2));
    sd_remote_ff = zeros(1, size(data_storage_rim(n).data,2));
    
    mean_roi_r2star = zeros(1, size(data_storage_rim(n).data,2));
    sd_roi_r2star = zeros(1, size(data_storage_rim(n).data,2));
    mean_remote_r2star = zeros(1, size(data_storage_rim(n).data,2));
    sd_remote_r2star = zeros(1, size(data_storage_rim(n).data,2));
    
    for tp = 1:size(data_storage_rim(n).data,2)
        
        roi_t1_temp = for_analysis_rim(n).metrics(tp).mean_roi_t1;
        if isempty(roi_t1_temp)
            mean_roi_t1(tp) = NaN;
            sd_roi_t1(tp) = NaN;
            mean_remote_t1(tp) = NaN;
            sd_remote_t1(tp) = NaN;
            
            mean_roi_ff(tp) = NaN;
            sd_roi_ff(tp) = NaN;
            mean_remote_ff(tp) = NaN;
            sd_remote_ff(tp) = NaN;
            
            mean_roi_r2star(tp) = NaN;
            sd_roi_r2star(tp) = NaN;
            mean_remote_r2star(tp) = NaN;
            sd_remote_r2star(tp) = NaN;
            
        else
            sd_roi_t1(tp) = for_analysis_rim(n).metrics(tp).sd_roi_t1;
            mean_roi_t1(tp) = for_analysis_rim(n).metrics(tp).mean_roi_t1;
            mean_remote_t1(tp) = for_analysis_rim(n).metrics(tp).mean_remote_t1;
            sd_remote_t1(tp) = for_analysis_rim(n).metrics(tp).sd_remote_t1;
            
            mean_roi_ff(tp) = for_analysis_rim(n).metrics(tp).mean_roi_ff;
            sd_roi_ff(tp) = for_analysis_rim(n).metrics(tp).sd_roi_ff;
            mean_remote_ff(tp) = for_analysis_rim(n).metrics(tp).mean_remote_ff;
            sd_remote_ff(tp) = for_analysis_rim(n).metrics(tp).sd_remote_ff;
            
            mean_roi_r2star(tp) = for_analysis_rim(n).metrics(tp).mean_roi_r2star;
            sd_roi_r2star(tp) = for_analysis_rim(n).metrics(tp).sd_roi_r2star;
            mean_remote_r2star(tp) = for_analysis_rim(n).metrics(tp).mean_remote_r2star;
            sd_remote_r2star(tp) = for_analysis_rim(n).metrics(tp).sd_remote_r2star;
        end
    end
    
    I = ~isnan(mean_roi_t1);
    
    figure('Position', [1000 500 900 600]);
    subplot(3,1,1);
    errorbar(xtcks_lr(I), mean_roi_t1(I), sd_roi_t1(I), '-o', 'LineWidth', 2); xlabel('Time point'); ylabel('mean T1 (ms)');
    hold on;  errorbar(xtcks_lr(I), mean_remote_t1(I), sd_remote_t1(I), '-o', 'LineWidth', 2);
    xticks(fliplr(xtcks_lr(I)));
    xticklabels(fliplr(time_points_lr(I)));
    set(gca, 'FontSize', 16);
    grid on;
    ylim([800 1600]); xlim([xtcks(1) xtcks(end)]);
    legend({'MI', 'Remote'}, 'Location', 'SouthEast');
    title(cat(2, name, '  T1 evolution'));
    
    
    subplot(3,1,2);
    errorbar(xtcks_lr(I), mean_roi_ff(I), sd_roi_ff(I), '-o', 'LineWidth', 2); xlabel('Time point'); ylabel('Fat Fraction (%)');
    hold on;  errorbar(xtcks_lr(I), mean_remote_ff(I), sd_remote_ff(I), '-o', 'LineWidth', 2);
    xticks(fliplr(xtcks_lr(I)));
    xticklabels(fliplr(time_points_lr(I)));
    set(gca, 'FontSize', 16);
    grid on;
    ylim([-5 20]); xlim([xtcks(1) xtcks(end)]);
    legend({'MI', 'Remote'}, 'Location', 'SouthEast');
    title(cat(2, name, '  FF evolution'));
    
    subplot(3,1,3);
    errorbar(xtcks_lr(I), mean_roi_r2star(I), sd_roi_r2star(I), '-o', 'LineWidth', 2); xlabel('Time point'); ylabel('R2 star (Hz)');
    hold on;  errorbar(xtcks_lr(I), mean_remote_r2star(I), sd_remote_r2star(I), '-o', 'LineWidth', 2);
    xticks(fliplr(xtcks_lr(I)));
    xticklabels(fliplr(time_points_lr(I)));
    set(gca, 'FontSize', 16);
    grid on;
    xlim([xtcks(1) xtcks(end)]);
    %ylim([-30 30])
    legend({'MI', 'Remote'}, 'Location', 'SouthEast');
    title(cat(2, name, '  R2 star evolution'));
    
    saveas(gcf, cat(2, long_name_save_dir, '/Time_Evolution_rim_03012021.png'));
end

close all;

%% At least save for_analysis metrics
metrics_save_dir = cat(2, base_dir, '/Results/');
if ~exist(metrics_save_dir, 'dir')
   mkdir(metrics_save_dir); 
end
%save(cat(2, metrics_save_dir, 'for_analysis.mat'), 'for_analysis');
save(cat(2, metrics_save_dir, 'for_analysis_rim.mat'), 'for_analysis_rim');

%% Make one tight-subplot for longitudinal image
metrics_save_dir = cat(2, base_dir, '/Results/');
%if ~exist('for_analysis', 'var')
%   load(cat(2, metrics_save_dir, 'for_analysis.mat'));
%end

if ~exist(cat(2, metrics_save_dir, 'for_analysis_rim.mat'), 'var')
    load(cat(2, metrics_save_dir, 'for_analysis_rim.mat'));
end

figure('Position', [200 500 1200 800]);
ha = tight_subplot(length(Names),3);
%for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end

time_points_lr = fliplr(time_points);
xtcks = 1:length(time_points);
xtcks_lr = fliplr(xtcks);
for n = 1:length(Names)
%for n = starting_point:starting_point
    name = Names{n};

    mean_roi_t1 = zeros(1, size(data_storage_rim(n).data,2));
    sd_roi_t1 = zeros(1, size(data_storage_rim(n).data,2));
    mean_remote_t1 = zeros(1, size(data_storage_rim(n).data,2));
    sd_remote_t1 = zeros(1, size(data_storage_rim(n).data,2));
    
    mean_roi_ff = zeros(1, size(data_storage_rim(n).data,2));
    sd_roi_ff = zeros(1, size(data_storage_rim(n).data,2));
    mean_remote_ff = zeros(1, size(data_storage_rim(n).data,2));
    sd_remote_ff = zeros(1, size(data_storage_rim(n).data,2));
    
    mean_roi_r2star = zeros(1, size(data_storage_rim(n).data,2));
    sd_roi_r2star = zeros(1, size(data_storage_rim(n).data,2));
    mean_remote_r2star = zeros(1, size(data_storage_rim(n).data,2));
    sd_remote_r2star = zeros(1, size(data_storage_rim(n).data,2));
    
    for tp = 1:size(data_storage_rim(n).data,2)
        
        roi_t1_temp = for_analysis_rim(n).metrics(tp).mean_roi_t1;
        if isempty(roi_t1_temp)
            mean_roi_t1(tp) = NaN;
            sd_roi_t1(tp) = NaN;
            mean_remote_t1(tp) = NaN;
            sd_remote_t1(tp) = NaN;
            
            mean_roi_ff(tp) = NaN;
            sd_roi_ff(tp) = NaN;
            mean_remote_ff(tp) = NaN;
            sd_remote_ff(tp) = NaN;
            
            mean_roi_r2star(tp) = NaN;
            sd_roi_r2star(tp) = NaN;
            mean_remote_r2star(tp) = NaN;
            sd_remote_r2star(tp) = NaN;
            
        else
            sd_roi_t1(tp) = for_analysis_rim(n).metrics(tp).sd_roi_t1;
            mean_roi_t1(tp) = for_analysis_rim(n).metrics(tp).mean_roi_t1;
            mean_remote_t1(tp) = for_analysis_rim(n).metrics(tp).mean_remote_t1;
            sd_remote_t1(tp) = for_analysis_rim(n).metrics(tp).sd_remote_t1;
            
            mean_roi_ff(tp) = for_analysis_rim(n).metrics(tp).mean_roi_ff;
            sd_roi_ff(tp) = for_analysis_rim(n).metrics(tp).sd_roi_ff;
            mean_remote_ff(tp) = for_analysis_rim(n).metrics(tp).mean_remote_ff;
            sd_remote_ff(tp) = for_analysis_rim(n).metrics(tp).sd_remote_ff;
            
            mean_roi_r2star(tp) = for_analysis_rim(n).metrics(tp).mean_roi_r2star;
            sd_roi_r2star(tp) = for_analysis_rim(n).metrics(tp).sd_roi_r2star;
            mean_remote_r2star(tp) = for_analysis_rim(n).metrics(tp).mean_remote_r2star;
            sd_remote_r2star(tp) = for_analysis_rim(n).metrics(tp).sd_remote_r2star;
        end
    end
    
    I = ~isnan(mean_roi_t1);
    axes(ha(3*(n-1)+1));
    errorbar(xtcks_lr(I), mean_roi_t1(I), sd_roi_t1(I), '-o', 'LineWidth', 2); 
    ylabel(name, 'FontSize', 12);
    hold on;  errorbar(xtcks_lr(I), mean_remote_t1(I), sd_remote_t1(I), '-o', 'LineWidth', 2);
    xticks(fliplr(xtcks_lr));
    xticklabels(fliplr(time_points_lr));
    grid on;
    ylim([800 1600]); xlim([xtcks(1) xtcks(end)+1]);
    ha(3*(n-1)+1).XAxis.FontSize = 16;
    if n == 1
        title(' T1 evolution', 'FontSize', 16);
    end
    if n == length(Names)
        legend({'MI', 'Remote'}, 'Location', 'SouthEast', 'FontSize', 14);
    end
    
    axes(ha(3*(n-1)+2));
    errorbar(xtcks_lr(I), mean_roi_ff(I), sd_roi_ff(I), '-o', 'LineWidth', 2); 
    %ylabel(name, 'FontSize', 12);
    hold on;  errorbar(xtcks_lr(I), mean_remote_ff(I), sd_remote_ff(I), '-o', 'LineWidth', 2);
    xticks(fliplr(xtcks_lr));
    xticklabels(fliplr(time_points_lr));
    grid on;
    ylim([-5 20]); xlim([xtcks(1) xtcks(end)+1]);
    ha(3*(n-1)+2).XAxis.FontSize = 16;
    if n == 1
        title(' FF evolution', 'FontSize', 16);
    end
    if n == length(Names)
        legend({'MI', 'Remote'}, 'Location', 'SouthEast', 'FontSize', 14);
    end
    
    axes(ha(3*(n-1)+3));
    errorbar(xtcks_lr(I), mean_roi_r2star(I), sd_roi_r2star(I), '-o', 'LineWidth', 2); 
    %ylabel(name, 'FontSize', 12);
    hold on;  errorbar(xtcks_lr(I), mean_remote_r2star(I), sd_remote_r2star(I), '-o', 'LineWidth', 2);
    xticks(fliplr(xtcks_lr));
    xticklabels(fliplr(time_points_lr));
    grid on;
    ylim([0 100]); xlim([xtcks(1) xtcks(end)+1]);
    ha(3*(n-1)+3).XAxis.FontSize = 16;
    if n == 1
        title(' R2* evolution', 'FontSize', 16);
    end
    if n == length(Names)
        legend({'MI', 'Remote'}, 'Location', 'SouthEast', 'FontSize', 14);
    end
end

set(ha(1:(3*length(Names)-3)),'XTickLabel',''); set(ha,'YTickLabel','');

t1fp_save_dir = cat(2, save_dir, '/T1FP_Compilation/');
if ~exist(t1fp_save_dir, 'dir')
    mkdir(t1fp_save_dir);
end
    
saveas(gcf, cat(2, t1fp_save_dir, '/Time_Evolution_rim_new_03012021.png'));

%% Plot the table & statistical analysis
% Please go to T1FP_Stats_Analysis.m

%% Get the matrics of T1, FF and R2star ( exclude edges of MI region)
%% Slice-by-Slice
% for ll = 1:length(sequence_label)
% T1 Map
% Images will be stored at img/<name>/overview/
% To exclude certain slices that has bad image quality

label_t1 = sequence_label{1};
label_lge = sequence_label{3};
label_t2star = sequence_label{2};
for_analysis_rim = struct;
data_storage_rim = struct;
get_2d_slice = @(x,s)x(:,:,s); 
for n = 1:length(Names)
%for n = length(Names):length(Names)
    name = Names{n};
    name_save_dir = cat(2, save_dir, name);
    if ~exist(name_save_dir, 'dir')
        mkdir(name_save_dir);
    end
    for_analysis_rim(n).Name = name;
    for_analysis_rim(n).metrics = struct;
    
    tp_count = 1;
    data_storage_rim(n).Name = name;
    data_storage_rim(n).data = struct;
    
    for tp = 1:length(time_points)
        time_point = time_points{end-tp+1};
        for_analysis_rim(n).metrics(tp).time_point = time_point;
        for_analysis_rim(n).metrics(tp).slices = struct;
        
        data_storage_rim(n).data(tp).time_point = time_point;
        data_storage_rim(n).data(tp).slices = struct;
        
        tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', name, '_', time_point,  '/');
        if ~exist(tp_dir, 'dir')
            disp(cat(2, 'No folder at: ', name, ' ', time_point));
        else
            % T1
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
            
            % Pre-QC
            [roi_in_myo_t1_new, remote_in_myo_t1_new, roi_t1_new, remote_t1_new, t1_new, myo_t1_new] = ...
                Func_status_check(status_check, n, name, tp_count, roi_in_myo_t1, remote_in_myo_t1,...
                roi_t1, remote_t1, t1, myo_t1);
            
            % remove edges of MI region
            roi_edg_t1_new = zeros(size(roi_in_myo_t1_new));
            for i = 1:size(roi_in_myo_t1_new, 3)
                roi_edg_t1_new(:,:,i)  = edge(squeeze(roi_in_myo_t1_new(:,:,i)),'Canny');
            end
            roi_rimmed_t1_new = (roi_in_myo_t1_new - roi_edg_t1_new)>0;
            
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
            ff = zeros(size(ff_map{1}.fwmc_ff,1), size(ff_map{2}.fwmc_ff, 2), length(ff_map));
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
            
            % Pre-QC
            [roi_in_myo_ff_new, remote_in_myo_ff_new, roi_ff_new, remote_ff_new, ff_new, myo_ff_new] = ...
                Func_status_check(status_check, n, name, tp_count, roi_in_myo_ff, remote_in_myo_ff,...
                roi_ff, remote_ff, ff, myo_ff);
            
            % remove edges of MI region
            roi_edg_ff_new = zeros(size(roi_in_myo_ff_new));
            for i = 1:size(roi_in_myo_ff_new, 3)
                roi_edg_ff_new(:,:,i)  = edge(squeeze(roi_in_myo_ff_new(:,:,i)),'Canny');
            end
            roi_rimmed_ff_new = (roi_in_myo_ff_new - roi_edg_ff_new)>0;
            
            
            % R2star Map
            r2star_map = cell(1, length(glob_names));
            for f = 1:length(r2star_map)
                r2star_map{f} = load(cat(2,  base_dir, '/FF_Data/',  name, '/', name, '_', time_point, '_', glob_names{f}, '.mat'), 'fwmc_r2star');
            end
            
            % convert ff_map to matrix
            r2star = zeros(size(r2star_map{1}.fwmc_r2star,1), size(r2star_map{1}.fwmc_r2star, 2), length(r2star_map));
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
            
            % Pre-QC
            [roi_in_myo_r2star_new, remote_in_myo_r2star_new, roi_r2star_new, remote_r2star_new, r2star_new, myo_r2star_new] = ...
                Func_status_check(status_check, n, name, tp_count, roi_in_myo_r2star, remote_in_myo_r2star,...
                roi_r2star, remote_r2star, r2star, myo_r2star);
            
            % remove edges of MI region
            roi_edg_r2star_new = zeros(size(roi_in_myo_r2star_new));
            for i = 1:size(roi_in_myo_r2star_new, 3)
                roi_edg_r2star_new(:,:,i)  = edge(squeeze(roi_in_myo_r2star_new(:,:,i)),'Canny');
            end
            roi_rimmed_r2star_new = (roi_in_myo_r2star_new - roi_edg_r2star_new)>0;
            
            for slc = 1:size(roi_rimmed_ff_new, 3)
                [row_roi, col_roi, v_roi] = find(roi_rimmed_ff_new(:,:,slc));
                [row_remote, col_remote, v_remote] = find(remote_in_myo_ff_new(:,:,slc));
                
                [row_roi_t1, col_roi_t1, v_roi_t1] = find(roi_rimmed_t1_new(:,:,slc));
                [row_remote_t1, col_remote_t1, v_remote_t1] = find(remote_in_myo_t1_new(:,:,slc));
                
                [row_roi_r2star, col_roi_r2star, v_roi_r2star] = find(roi_rimmed_r2star_new(:,:,slc));
                [row_remote_r2star, col_remote_r2star, v_remote_r2star] = find(remote_in_myo_r2star_new(:,:,slc));
                
                roi_ff_array = zeros(1, length(row_roi));
                remote_ff_array = zeros(1, length(row_remote));
                roi_r2star_array = zeros(1, length(row_roi_r2star));
                remote_r2star_array = zeros(1, length(row_remote_r2star));
                roi_t1_array = zeros(1, length(row_roi_t1));
                remote_t1_array = zeros(1, length(row_remote_t1));
                
                roi_ff_new_temp = get_2d_slice(roi_ff_new, slc);
                remote_ff_new_temp = get_2d_slice(remote_ff_new, slc);
                roi_r2star_new_temp = get_2d_slice(roi_r2star_new, slc);
                remote_r2star_new_temp = get_2d_slice(remote_r2star_new, slc);
                roi_t1_new_temp = get_2d_slice(roi_t1_new, slc);
                remote_t1_new_temp = get_2d_slice(remote_t1_new, slc);
                
                for fff = 1:length(row_roi)
                    roi_ff_array(fff) = roi_ff_new_temp(row_roi(fff), col_roi(fff));
                end
                
                for fff = 1:length(row_remote)
                    remote_ff_array(fff) = remote_ff_new_temp(row_remote(fff), col_remote(fff));
                end
                
                for fff = 1:length(row_roi_r2star)
                    roi_r2star_array(fff) = roi_r2star_new_temp(row_roi_r2star(fff), col_roi_r2star(fff));
                end
                
                for fff = 1:length(row_remote_r2star)
                    remote_r2star_array(fff) = remote_r2star_new_temp(row_remote_r2star(fff), col_remote_r2star(fff));
                end
                
                for fff = 1:length(row_roi_t1)
                    roi_t1_array(fff) = roi_t1_new_temp(row_roi_t1(fff), col_roi_t1(fff));
                end
                
                for fff = 1:length(row_remote_t1)
                    remote_t1_array(fff) = remote_t1_new_temp(row_remote_t1(fff), col_remote_t1(fff));
                end
                
                roi_ff_array(roi_ff_array < 0) = 0;
                roi_ff_array(roi_ff_array > 100) = 100;
                remote_ff_array(remote_ff_array < 0) = 0;
                remote_ff_array(remote_ff_array > 100) = 100;
                
                roi_r2star_array(roi_r2star_array > 100) = 100;
                remote_r2star_array(remote_r2star_array > 100) = 100;
                
                
%                 % Plot Histogram of ROI vs Remote (Gross view)
%                 figure();
%                 h_ff = histogram(roi_ff_array, 'Normalization', 'probability');xlabel('Fat Fraction (%)'); ylabel('Frequency');
%                 NumBins = h_ff.NumBins;
%                 BinWidth = h_ff.BinWidth;
%                 BinEdges = h_ff.BinEdges;
%                 hold on;
%                 h_remote = histogram(remote_ff_array, 'Normalization', 'probability');
%                 h_remote.BinEdges = BinEdges;
%                 set(gca, 'FontSize', 16); title([name, ' ', time_point, ' Slice', num2str(slc)]);
%                 xlim([0 100]);
%                 legend({'MI', 'Remote'});
%                 name_tp_dir = cat(2, name_save_dir, '/', name, '_', time_point);
%                 if ~exist(name_tp_dir, 'dir')
%                     mkdir(name_tp_dir);
%                 end
%                 %saveas(gcf, cat(2, name_save_dir, '/', name, '_', time_point, '/Histogram_FF_Whole_rim.png'));
%                 
%                 % R2star
%                 figure();
%                 h_r2star = histogram(roi_r2star_array, 'Normalization', 'probability');xlabel('R2* (Hz)'); ylabel('Frequency');
%                 NumBins = h_r2star.NumBins;
%                 BinWidth = h_r2star.BinWidth;
%                 BinEdges = h_r2star.BinEdges;
%                 hold on;
%                 h_r2star_remote = histogram(remote_r2star_array, 'Normalization', 'probability');
%                 h_r2star_remote.BinEdges = BinEdges;
%                 set(gca, 'FontSize', 16); title([name, ' ', time_point, ' Slice', num2str(slc)]);
%                 xlim([0 100]);
%                 legend({'MI', 'Remote'});
%                 name_tp_dir = cat(2, name_save_dir, '/', name, '_', time_point);
%                 if ~exist(name_tp_dir, 'dir')
%                     mkdir(name_tp_dir);
%                 end
%                 %saveas(gcf, cat(2, name_save_dir, '/', name, '_', time_point, '/Histogram_R2star_Whole_rim.png'));
%                 
%                 % T1
%                 figure();
%                 h_t1 = histogram(roi_t1_array, 'Normalization', 'probability');xlabel('T1 (ms)'); ylabel('Frequency');
%                 NumBins = h_t1.NumBins;
%                 BinWidth = h_t1.BinWidth;
%                 BinEdges = h_t1.BinEdges;
%                 hold on;
%                 h_t1_remote = histogram(remote_t1_array, 'Normalization', 'probability');
%                 h_t1_remote.BinEdges = BinEdges;
%                 set(gca, 'FontSize', 16); title([name, ' ', time_point, ' Slice', num2str(slc)]);
%                 %xlim([0 100]);
%                 legend({'MI', 'Remote'});
%                 
%                 name_tp_dir = cat(2, name_save_dir, '/', name, '_', time_point);
%                 if ~exist(name_tp_dir, 'dir')
%                     mkdir(name_tp_dir);
%                 end
                %saveas(gcf, cat(2, name_save_dir, '/', name, '_', time_point, '/Histogram_T1_Whole_rim.png'));
                
                
                name_tp = cat(2, name, '_', time_point);
                
                for_analysis_rim(n).metrics(tp).slices(slc).slice = ['Slice', num2str(slc)];
                for_analysis_rim(n).metrics(tp).slices(slc).mean_roi_t1 = mean(nonzeros(roi_t1_array));
                for_analysis_rim(n).metrics(tp).slices(slc).sd_roi_t1 = std(nonzeros(roi_t1_array));
                for_analysis_rim(n).metrics(tp).slices(slc).mean_remote_t1 = mean(nonzeros(remote_t1_array));
                for_analysis_rim(n).metrics(tp).slices(slc).sd_remote_t1 = std(nonzeros(remote_t1_array));
                
                for_analysis_rim(n).metrics(tp).slices(slc).mean_roi_ff = mean(roi_ff_array);
                for_analysis_rim(n).metrics(tp).slices(slc).sd_roi_ff = std(roi_ff_array);
                for_analysis_rim(n).metrics(tp).slices(slc).mean_remote_ff = mean(remote_ff_array);
                for_analysis_rim(n).metrics(tp).slices(slc).sd_remote_ff = std(remote_ff_array);
                
                for_analysis_rim(n).metrics(tp).slices(slc).mean_roi_r2star = mean(nonzeros(roi_r2star_array));
                for_analysis_rim(n).metrics(tp).slices(slc).sd_roi_r2star = std(nonzeros(roi_r2star_array));
                for_analysis_rim(n).metrics(tp).slices(slc).mean_remote_r2star = mean(nonzeros(remote_r2star_array));
                for_analysis_rim(n).metrics(tp).slices(slc).sd_remote_r2star = std(nonzeros(remote_r2star_array));
                
                
                data_storage_rim(n).data(tp).slices(slc).slice = ['Slice', num2str(slc)];
                data_storage_rim(n).data(tp).slices(slc).roi_ff_array = roi_ff_array;
                data_storage_rim(n).data(tp).slices(slc).remote_ff_array = remote_ff_array;
                data_storage_rim(n).data(tp).slices(slc).roi_r2star_array = roi_r2star_array;
                data_storage_rim(n).data(tp).slices(slc).remote_r2star_array = remote_r2star_array;
                data_storage_rim(n).data(tp).slices(slc).roi_t1_array = roi_t1_array;
                data_storage_rim(n).data(tp).slices(slc).remote_t1_array = remote_t1_array;
            end
            tp_count = tp_count + 1;
        end
    end
    close all;
end

metrics_save_dir = cat(2, base_dir, '/Results/');
if ~exist(metrics_save_dir, 'dir')
   mkdir(metrics_save_dir); 
end
save(cat(2, metrics_save_dir, 'for_analysis_rim_slices.mat'), 'for_analysis_rim');
save(cat(2, metrics_save_dir, 'data_storage_rim_slices.mat'), 'data_storage_rim');
