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
Names = {'Merry', 'Ryn', 'Hope'};
time_points = {'Exvivo'};

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

sequence_label = {'T1_CMR', 'T2star_CMR'};
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
label_t2star = sequence_label{2};
for_analysis = struct;


%for n = 1:length(Names)
for n = 3:3
% for n = starting_point:starting_point
% Do not need to pull up images for baseline
    name = Names{n};
    name_save_dir = cat(2, save_dir, name);
    if ~exist(name_save_dir, 'dir')
        mkdir(name_save_dir);
    end
    
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
            ff = zeros(size(ff_map{1}.fwmc_ff,1), size(ff_map{1}.fwmc_ff, 2),size(ff_map{1}.fwmc_ff,3),  length(ff_map));
            
            for f = 1:length(ff_map)
                   ff(:,:,:,f) = ff_map{f}.fwmc_ff;
            end
            idx_reordered = Func_AlignSliceLoc(slc_array_t1, slc_array_ff);
            mask_myocardium_3D = mask_myocardium_3D(:,:,idx_reordered);
            freeROIMask_3D = freeROIMask_3D(:,:,idx_reordered);
            myoRefMask_3D = myoRefMask_3D(:,:,idx_reordered);
            ff = squeeze(ff(:,:,idx_reordered,:));
            roi_in_myo_ff = mask_myocardium_3D .* freeROIMask_3D;
            remote_in_myo_ff = mask_myocardium_3D .* myoRefMask_3D;
            roi_ff = roi_in_myo_ff .* ff;
            remote_ff = remote_in_myo_ff .* ff;
            myo_ff = mask_myocardium_3D;
            
         
            % R2star Map
            r2star_map = cell(1, length(glob_names));
            for f = 1:length(r2star_map)
                r2star_map{f} = load(cat(2,  base_dir, '/FF_Data/',  name, '/', name, '_', time_point, '_', glob_names{f}, '.mat'), 'fwmc_r2star');
            end
            
            % convert ff_map to matrix
            r2star = zeros(size(r2star_map{1}.fwmc_r2star,1), size(r2star_map{1}.fwmc_r2star, 2), size(r2star_map{1}.fwmc_r2star, 3), length(r2star_map));
            for f = 1:length(ff_map)
                   r2star(:,:,:,f) = r2star_map{f}.fwmc_r2star;
            end
            
            idx_reordered = Func_AlignSliceLoc(slc_array_t1, slc_array_ff);
            
            r2star = squeeze(r2star(:,:,idx_reordered,:));
            roi_in_myo_r2star = mask_myocardium_3D .* freeROIMask_3D;
            remote_in_myo_r2star = mask_myocardium_3D .* myoRefMask_3D;
            roi_r2star = roi_in_myo_r2star .* r2star;
            remote_r2star = remote_in_myo_r2star .* r2star;
            myo_r2star = mask_myocardium_3D;
            
            for i = 1:size(roi_t1, 3)
                %for i = 1:1
                figure();
                subplot(2,2,1);
                imagesc(t1(:,:,i)); axis image; colorbar;
                title(['T1 Map: ', name, ' ', time_point, '  Slice = ', num2str(i)]);
                colormap(brewermap([],'*RdYlBu'));

                
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

%% Get the matrics of T1, FF and R2star ( exclude edges of MI region)
% T1 Map
% Images will be stored at img/<name>/overview/
% To exclude certain slices that has bad image quality

label_t1 = sequence_label{1};
label_t2star = sequence_label{2};
% for_analysis = struct;
for_analysis_rim = struct;
data_storage_rim = struct;

%for n = 1:length(Names)
for n = 2:2
    % for n = starting_point:starting_point
    % Do not need to pull up images for baseline
    name = Names{n};
    name_save_dir = cat(2, save_dir, name);
    if ~exist(name_save_dir, 'dir')
        mkdir(name_save_dir);
    end
    
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
            %[roi_in_myo_t1_new, remote_in_myo_t1_new, roi_t1_new, remote_t1_new, t1_new, myo_t1_new] = ...
            %    Func_status_check(status_check, n, name, tp_count, roi_in_myo_t1, remote_in_myo_t1,...
            %    roi_t1, remote_t1, t1, myo_t1);
            
            % remove edges of MI region
            roi_edg_t1 = zeros(size(roi_in_myo_t1));
            for i = 1:size(roi_in_myo_t1, 3)
                roi_edg_t1(:,:,i)  = edge(squeeze(roi_in_myo_t1(:,:,i)),'Canny');
            end
            roi_rimmed_t1 = (roi_in_myo_t1 - roi_edg_t1)>0;
            [row_roi_t1, col_roi_t1, v_roi_t1] = find(roi_rimmed_t1);
            [row_remote_t1, col_remote_t1, v_remote_t1] = find(remote_in_myo_t1);
            
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
            ff = zeros(size(ff_map{1}.fwmc_ff,1), size(ff_map{1}.fwmc_ff, 2),size(ff_map{1}.fwmc_ff,3),  length(ff_map));
            
            for f = 1:length(ff_map)
                ff(:,:,:,f) = ff_map{f}.fwmc_ff;
            end
            idx_reordered = Func_AlignSliceLoc(slc_array_t1, slc_array_ff);
            mask_myocardium_3D = mask_myocardium_3D(:,:,idx_reordered);
            freeROIMask_3D = freeROIMask_3D(:,:,idx_reordered);
            myoRefMask_3D = myoRefMask_3D(:,:,idx_reordered);
            ff = squeeze(ff(:,:,idx_reordered,:));
            roi_in_myo_ff = mask_myocardium_3D .* freeROIMask_3D;
            remote_in_myo_ff = mask_myocardium_3D .* myoRefMask_3D;
            roi_ff = roi_in_myo_ff .* ff;
            remote_ff = remote_in_myo_ff .* ff;
            myo_ff = mask_myocardium_3D;
            
            % Pre-QC
            %[roi_in_myo_ff_new, remote_in_myo_ff_new, roi_ff_new, remote_ff_new, ff_new, myo_ff_new] = ...
            %    Func_status_check(status_check, n, name, tp_count, roi_in_myo_ff, remote_in_myo_ff,...
            %    roi_ff, remote_ff, ff, myo_ff);
            
            % remove edges of MI region
            roi_edg_ff = zeros(size(roi_in_myo_ff));
            for i = 1:size(roi_in_myo_ff, 3)
                roi_edg_ff(:,:,i)  = edge(squeeze(roi_in_myo_ff(:,:,i)),'Canny');
            end
            roi_rimmed_ff = (roi_in_myo_ff - roi_edg_ff)>0;
            
            [row_roi, col_roi, v_roi] = find(roi_rimmed_ff);
            [row_remote, col_remote, v_remote] = find(remote_in_myo_ff);
            
            
            % R2star Map
            r2star_map = cell(1, length(glob_names));
            for f = 1:length(r2star_map)
                r2star_map{f} = load(cat(2,  base_dir, '/FF_Data/',  name, '/', name, '_', time_point, '_', glob_names{f}, '.mat'), 'fwmc_r2star');
            end
            
            % convert ff_map to matrix
            r2star = zeros(size(r2star_map{1}.fwmc_r2star,1), size(r2star_map{1}.fwmc_r2star, 2), size(r2star_map{1}.fwmc_r2star, 3), length(r2star_map));
            for f = 1:length(ff_map)
                r2star(:,:,:,f) = r2star_map{f}.fwmc_r2star;
            end
            
            idx_reordered = Func_AlignSliceLoc(slc_array_t1, slc_array_ff);
            
            r2star = squeeze(r2star(:,:,idx_reordered,:));
            roi_in_myo_r2star = mask_myocardium_3D .* freeROIMask_3D;
            remote_in_myo_r2star = mask_myocardium_3D .* myoRefMask_3D;
            roi_r2star = roi_in_myo_r2star .* r2star;
            remote_r2star = remote_in_myo_r2star .* r2star;
            myo_r2star = mask_myocardium_3D;
            
            
            % Pre-QC
            %[roi_in_myo_r2star_new, remote_in_myo_r2star_new, roi_r2star_new, remote_r2star_new, r2star_new, myo_r2star_new] = ...
            %    Func_status_check(status_check, n, name, tp_count, roi_in_myo_r2star, remote_in_myo_r2star,...
            %    roi_r2star, remote_r2star, r2star, myo_r2star);
            
            % remove edges of MI region
            roi_edg_r2star = zeros(size(roi_in_myo_r2star));
            for i = 1:size(roi_in_myo_r2star, 3)
                roi_edg_r2star(:,:,i)  = edge(squeeze(roi_in_myo_r2star(:,:,i)),'Canny');
            end
            roi_rimmed_r2star = (roi_in_myo_r2star - roi_edg_r2star)>0;
            [row_roi_r2star, col_roi_r2star, v_roi_r2star] = find(roi_rimmed_r2star);
            [row_remote_r2star, col_remote_r2star, v_remote_r2star] = find(remote_in_myo_r2star);
            
            
            
            roi_ff_array = zeros(1, length(row_roi));
            remote_ff_array = zeros(1, length(row_remote));
            roi_r2star_array = zeros(1, length(row_roi_r2star));
            remote_r2star_array = zeros(1, length(row_remote_r2star));
            roi_t1_array = zeros(1, length(row_roi_t1));
            remote_t1_array = zeros(1, length(row_remote_t1));
            
            for fff = 1:length(row_roi)
                roi_ff_array(fff) = roi_ff(row_roi(fff), col_roi(fff));
            end
            
            for fff = 1:length(row_remote)
                remote_ff_array(fff) = remote_ff(row_remote(fff), col_remote(fff));
            end
            
            for fff = 1:length(row_roi_r2star)
                roi_r2star_array(fff) = roi_r2star(row_roi_r2star(fff), col_roi_r2star(fff));
            end
            
            for fff = 1:length(row_remote_r2star)
                remote_r2star_array(fff) = remote_r2star(row_remote_r2star(fff), col_remote_r2star(fff));
            end
            
            for fff = 1:length(row_roi_t1)
                roi_t1_array(fff) = roi_t1(row_roi_t1(fff), col_roi_t1(fff));
            end
            
            for fff = 1:length(row_remote_t1)
                remote_t1_array(fff) = remote_t1(row_remote_t1(fff), col_remote_t1(fff));
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
            
            tp_count = tp_count + 1;
        end
    end
    close all;
end

metrics_save_dir = cat(2, base_dir, '/Results/');
if ~exist(metrics_save_dir, 'dir')
   mkdir(metrics_save_dir); 
end
save(cat(2, metrics_save_dir, 'data_storage_rim.mat'), 'data_storage_rim');

