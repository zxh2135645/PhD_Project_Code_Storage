close all;
clear all;

% T1 and FF and true T2* mapping analysis
% This version is developed in 10212020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% img/<name>/overview/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('../function/');
base_dir = uigetdir;
contour_glob = glob(cat(2, base_dir, '/ContourData/*'));

Names = ExtractNames(contour_glob);

% labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};
time_points = {'0D_baseline', '0D_occl', '0D', '1D', '3D', '7D', '21D', '28D', '8WK', '12WK', '14WK', '6MO', '9MO', '1YR', '15YR'};

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
name_check = 'Evelyn';
starting_point = find(strcmp(name_check, Names),1);

save_dir = cat(2, base_dir, '/img/');

%% What is the metrics? (Old version used for longitudinal analysis of Merry)
%for ll = 1:length(sequence_label)
% T1 Map
for ll = 1:1
    label = sequence_label{ll};
    T1_analysis = struct;
    for n = starting_point:length(Names)
    %for n = starting_point:starting_point
        name = Names{n};
        name_save_dir = cat(2, save_dir, name);
        if ~exist(name_save_dir, 'dir')
            mkdir(name_save_dir);
        end
        
        T1_analysis(n).Name = name;
        
        metrics = zeros(1, length(time_points));
        metrics_remote = zeros(1, length(time_points));
        metrics_sd = zeros(1, length(time_points));
        metrics_sd_remote = zeros(1, length(time_points));
        for tp = 1:length(time_points)
            %for tp = 1:length(time_points)
            time_point = time_points{end-tp+1};
            tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', name, '_', time_point,  '/');
            if ~exist(tp_dir, 'dir')
                disp(cat(2, 'No folder at: ', name, ' ', time_point));
            else
                myo_glob = glob(cat(2, tp_dir, label, '/', anatomy_label{5}, '/*'));
                roi_glob = glob(cat(2, tp_dir, label, '/',anatomy_label{3}, '/*'));
                remote_glob = glob(cat(2, tp_dir, label, '/',anatomy_label{6}, '/*'));
                
                load(cat(2, tp_dir, label, '/', label, '_vol_img_3D.mat'));
                load(myo_glob{1});
                load(roi_glob{1});
                load(remote_glob{1});
                
                roi_in_myo = mask_myocardium_3D .* freeROIMask_3D;
                remote_in_myo = mask_myocardium_3D .* myoRefMask_3D;
                roi = roi_in_myo .* vol_img_3D;
                remote = remote_in_myo .* vol_img_3D;
                
                figure();
                if size(roi, 3) > 4
                    for i = 1:size(roi, 3)
                        subplot(2,3,i);
                        imagesc(roi(:,:,i)); axis image;
                        title([name, ' ', time_point, '  Slice = ', num2str(i)]);
                    end
                else
                    for i = 1:size(roi, 3)
                        subplot(2,2,i);
                        imagesc(roi(:,:,i)); axis image;
                        title([name, ' ', time_point, '  Slice = ', num2str(i)]);
                    end
                end
                saveas(gcf, cat(2, name_save_dir, '\', name, '_', time_point, '_MI.png'));
                %             nel = unique(nonzeros(roi));
                %             count = zeros(1, numel(nel));
                %             for nn = 1:length(nel)
                %                 count(nn) = sum(nonzeros(roi) == nel(nn));
                %             end
                
                figure();
                histogram(nonzeros(roi),'Normalization', 'probability');xlabel('T1 (ms)'); ylabel('Frequency');
                hold on;
                histogram(nonzeros(remote),'Normalization', 'probability');
                set(gca, 'FontSize', 16); title([name, ' ', time_point]);
                xlim([0 2000]);
                legend({'MI', 'Remote'});
                
                saveas(gcf, cat(2, name_save_dir, '\', name, '_', time_point, '_histogram.png'));
                %yyaxis right;
                %plot(nel, count);
                
                metrics(end-tp+1) = mean(nonzeros(roi));
                metrics_remote(end-tp+1) = mean(nonzeros(remote));
                metrics_sd(end-tp+1) = std(nonzeros(roi));
                metrics_sd_remote(end-tp+1) = std(nonzeros(remote));
                break;
            end
            %
            %             figure('Position', [1000 500 900 600]);
            %             errorbar(metrics, metrics_sd, '-o', 'LineWidth', 2); xlabel('Time point'); ylabel('mean T1 (ms)');
            %             hold on;  errorbar(metrics_remote, metrics_sd_remote, '-o', 'LineWidth', 2);
            %             xticklabels(time_points);
            %             set(gca, 'FontSize', 16);
            %             grid on;
            %             legend({'MI', 'Remote'}, 'Location', 'SouthEast');
            
        end
        T1_analysis(n).Name = name;
        T1_analysis(n).roi_mean = metrics;
        T1_analysis(n).roi_sd = metrics_sd;
        T1_analysis(n).remote_mean = metrics_remote;
        T1_analysis(n).remote_sd = metrics_sd_remote;
    end
end

%% Pull LGE, T1 Map, FF map and True T2* map
%for ll = 1:length(sequence_label)
% T1 Map
% Images will be stored at img/<name>/overview/

label_t1 = sequence_label{1};
label_lge = sequence_label{3};
label_t2star = sequence_label{2};
for_analysis = struct;

for n = starting_point:length(Names)
%for n = starting_point:starting_point
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
            break;
        end
        
        %
        %             figure('Position', [1000 500 900 600]);
        %             errorbar(metrics, metrics_sd, '-o', 'LineWidth', 2); xlabel('Time point'); ylabel('mean T1 (ms)');
        %             hold on;  errorbar(metrics_remote, metrics_sd_remote, '-o', 'LineWidth', 2);
        %             xticklabels(time_points);
        %             set(gca, 'FontSize', 16);
        %             grid on;
        %             legend({'MI', 'Remote'}, 'Location', 'SouthEast');

    end
end

%% Analysis

mean_t1 = zeros(1, length(for_analysis));
mean_t1_remote = zeros(1, length(for_analysis));
mean_ff = zeros(1, length(for_analysis));
mean_r2star = zeros(1, length(for_analysis));
name_cell = cell(1, length(for_analysis));

sd_t1 = zeros(1, length(for_analysis));
sd_t1_remote = zeros(1, length(for_analysis));
sd_ff = zeros(1, length(for_analysis));
sd_r2star = zeros(1, length(for_analysis));

for i = 1:length(for_analysis)
    mean_t1(i) = for_analysis(i).metrics.mean_roi_t1;
    mean_t1_remote(i) = for_analysis(i).metrics.mean_remote_t1;
    mean_ff(i) = for_analysis(i).metrics.mean_roi_ff;
    mean_r2star(i) = for_analysis(i).metrics.mean_roi_r2star;
    name_cell{i} = for_analysis(i).Name;
    
    sd_t1(i) = for_analysis(i).metrics.sd_roi_t1;
    sd_t1_remote(i) = for_analysis(i).metrics.sd_remote_t1;
    sd_ff(i) = for_analysis(i).metrics.sd_roi_ff;
    sd_r2star(i) = for_analysis(i).metrics.sd_roi_r2star;
end

figure();
X = 1:length(for_analysis);
errorbar(X, mean_t1, sd_t1, 'LineWidth', 2);
hold on;
%plot(X, mean_t1_remote);
xticks(X);
xticklabels(name_cell);
yyaxis right; errorbar(X, mean_ff, sd_ff, 'LineWidth', 2); 
legend({'T1', 'FF'});
%plot(mean_r2star);

%% Get the matrics of T1, FF and R2star
%for ll = 1:length(sequence_label)
% T1 Map
% Images will be stored at img/<name>/overview/
% To exclude certain slices that has bad image quality

label_t1 = sequence_label{1};
label_lge = sequence_label{3};
label_t2star = sequence_label{2};
for_analysis = struct;

for n = starting_point:length(Names)
%for n = starting_point:starting_point
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
            [row_roi_t1, col_roi_t1, v_roi_t1] = find(roi_in_myo_t1);
            [row_remote_t1, col_remote_t1, v_remote_t1] = find(remote_in_myo_t1);
            
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
            [row_roi, col_roi, v_roi] = find(roi_in_myo_ff);
            [row_remote, col_remote, v_remote] = find(remote_in_myo_ff);
            
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
            [row_roi_r2star, col_roi_r2star, v_roi_r2star] = find(roi_in_myo_r2star);
            [row_remote_r2star, col_remote_r2star, v_remote_r2star] = find(remote_in_myo_r2star);
            
            name_tp = cat(2, name, '_', time_point);
            for_analysis(n).Name = name;
            for_analysis(n).metrics = struct;
            
            for_analysis(n).metrics.time_point = time_point;
            for_analysis(n).metrics.mean_roi_t1 = mean(nonzeros(roi_t1));
            for_analysis(n).metrics.sd_roi_t1 = std(nonzeros(roi_t1));
            for_analysis(n).metrics.mean_remote_t1 = mean(nonzeros(remote_t1));
            for_analysis(n).metrics.sd_remote_t1 = std(nonzeros(remote_t1));
            
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
%            x_ff = zeros(1, NumBins);
%             for nb = 1:NumBins
%                 x_ff(nb) = BinEdges(nb) + BinWidth/2;
%             end
%             y_ff = h_ff.Values;
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

            % Plot Histogram of ROI vs Remote (Slice by slice)
            % Missing remote in ZZ slice 3
            for slc = 1:size(roi_ff, 3)
                % FF
                roi_ff_slc = roi_ff(:,:,slc);
                remote_ff_slc = remote_ff(:,:,slc);
                [row_roi, col_roi, v_roi] = find(roi_in_myo_ff(:,:,slc));
                roi_ff_array = zeros(1, length(row_roi));
                for fff = 1:length(row_roi)
                    roi_ff_array(fff) = roi_ff_slc(row_roi(fff), col_roi(fff));
                end
                
                [row_remote, col_remote, v_remote] = find(remote_in_myo_ff(:,:,slc));
                remote_ff_array = zeros(1, length(row_remote));
                for fff = 1:length(row_remote)
                    remote_ff_array(fff) = remote_ff_slc(row_remote(fff), col_remote(fff));
                end
                roi_ff_array(roi_ff_array < 0) = 0;
                roi_ff_array(roi_ff_array > 100) = 100;
                remote_ff_array(remote_ff_array < 0) = 0;
                remote_ff_array(remote_ff_array > 100) = 100;
                
                figure();
                h_ff = histogram(roi_ff_array, 'Normalization', 'probability');xlabel('Fat Fraction (%)'); ylabel('Frequency');
                NumBins = h_ff.NumBins;
                BinWidth = h_ff.BinWidth;
                BinEdges = h_ff.BinEdges;
                x_ff = zeros(1, NumBins);
                hold on;
                h_remote = histogram(remote_ff_array, 'Normalization', 'probability');
                h_remote.BinEdges = BinEdges;
                set(gca, 'FontSize', 16); title([name, ' ', time_point, ' Slice=', num2str(slc)]);
                xlim([0 100]);
                legend({'MI', 'Remote'});
                name_tp_dir = cat(2, name_save_dir, '/', name, '_', time_point);
                if ~exist(name_tp_dir, 'dir')
                    mkdir(name_tp_dir);
                end
                saveas(gcf, cat(2, name_save_dir, '/', name, '_', time_point, '/Histogram_FF_Slice', num2str(slc),'.png'));
                
                % R2star
                roi_r2star_slc = roi_r2star(:,:,slc);
                remote_r2star_slc = remote_r2star(:,:,slc);
                [row_roi_r2star, col_roi_r2star, v_roi_r2star] = find(roi_in_myo_r2star(:,:,slc));
                roi_r2star_array = zeros(1, length(row_roi_r2star));
                for fff = 1:length(row_roi_r2star)
                    roi_r2star_array(fff) = roi_r2star_slc(row_roi_r2star(fff), col_roi_r2star(fff));
                end
                
                [row_remote_r2star, col_remote_r2star, v_remote_r2star] = find(remote_in_myo_r2star(:,:,slc));
                remote_r2star_array = zeros(1, length(row_remote_r2star));
                for fff = 1:length(row_remote_r2star)
                    remote_r2star_array(fff) = remote_r2star_slc(row_remote_r2star(fff), col_remote_r2star(fff));
                end
                roi_r2star_array(roi_r2star_array > 100) = 100;
                remote_r2star_array(remote_r2star_array > 100) = 100;
                
                figure();
                h_r2star = histogram(roi_r2star_array, 'Normalization', 'probability');xlabel('R2* (Hz)'); ylabel('Frequency');
                NumBins = h_r2star.NumBins;
                BinWidth = h_r2star.BinWidth;
                BinEdges = h_r2star.BinEdges;
                x_r2star = zeros(1, NumBins);
                hold on;
                h_r2star_remote = histogram(remote_r2star_array, 'Normalization', 'probability');
                h_r2star_remote.BinEdges = BinEdges;
                set(gca, 'FontSize', 16); title([name, ' ', time_point, ' Slice=', num2str(slc)]);
                xlim([0 100]);
                legend({'MI', 'Remote'});
                name_tp_dir = cat(2, name_save_dir, '/', name, '_', time_point);
                if ~exist(name_tp_dir, 'dir')
                    mkdir(name_tp_dir);
                end
                saveas(gcf, cat(2, name_save_dir, '/', name, '_', time_point, '/Histogram_R2star_Slice', num2str(slc),'.png'));
                
                % T1
                roi_t1_slc = roi_t1(:,:,slc);
                remote_t1_slc = remote_t1(:,:,slc);
                [row_roi_t1, col_roi_t1, v_roi_t1] = find(roi_in_myo_t1(:,:,slc));
                roi_t1_array = zeros(1, length(row_roi_t1));
                for fff = 1:length(row_roi_t1)
                    roi_t1_array(fff) = roi_t1_slc(row_roi_t1(fff), col_roi_t1(fff));
                end
                
                [row_remote_t1, col_remote_t1, v_remote_t1] = find(remote_in_myo_t1(:,:,slc));
                remote_t1_array = zeros(1, length(row_remote_t1));
                for fff = 1:length(row_remote_t1)
                    remote_t1_array(fff) = remote_t1_slc(row_remote_t1(fff), col_remote_t1(fff));
                end
                
                figure();
                h_t1 = histogram(roi_t1_array, 'Normalization', 'probability');xlabel('T1 (ms)'); ylabel('Frequency');
                NumBins = h_t1.NumBins;
                BinWidth = h_t1.BinWidth;
                BinEdges = h_t1.BinEdges;
                x_t1 = zeros(1, NumBins);
                hold on;
                h_t1_remote = histogram(remote_t1_array, 'Normalization', 'probability');
                h_t1_remote.BinEdges = BinEdges;
                set(gca, 'FontSize', 16); title([name, ' ', time_point, ' Slice=', num2str(slc)]);
                %xlim([0 100]);
                legend({'MI', 'Remote'});
                name_tp_dir = cat(2, name_save_dir, '/', name, '_', time_point);
                if ~exist(name_tp_dir, 'dir')
                    mkdir(name_tp_dir);
                end
                saveas(gcf, cat(2, name_save_dir, '/', name, '_', time_point, '/Histogram_T1_Slice', num2str(slc),'.png'));
            end
            
            for_analysis(n).metrics.mean_roi_ff = mean(roi_ff_array);
            for_analysis(n).metrics.sd_roi_ff = std(roi_ff_array);
            for_analysis(n).metrics.mean_remote_ff = mean(remote_ff_array);
            for_analysis(n).metrics.sd_remote_ff = std(remote_ff_array);
            
            for_analysis(n).metrics.mean_roi_r2star = mean(nonzeros(roi_r2star));
            for_analysis(n).metrics.sd_roi_r2star = std(nonzeros(roi_r2star));
            for_analysis(n).metrics.mean_remote_r2star = mean(nonzeros(remote_r2star));
            for_analysis(n).metrics.sd_remote_r2star = std(nonzeros(remote_r2star));
            
            close all;
            break;
        end
        
        %
        %             figure('Position', [1000 500 900 600]);
        %             errorbar(metrics, metrics_sd, '-o', 'LineWidth', 2); xlabel('Time point'); ylabel('mean T1 (ms)');
        %             hold on;  errorbar(metrics_remote, metrics_sd_remote, '-o', 'LineWidth', 2);
        %             xticklabels(time_points);
        %             set(gca, 'FontSize', 16);
        %             grid on;
        %             legend({'MI', 'Remote'}, 'Location', 'SouthEast');

    end
end

%% Analysis

mean_t1 = zeros(1, length(for_analysis));
mean_t1_remote = zeros(1, length(for_analysis));
mean_ff = zeros(1, length(for_analysis));
mean_r2star = zeros(1, length(for_analysis));
name_cell = cell(1, length(for_analysis));

sd_t1 = zeros(1, length(for_analysis));
sd_t1_remote = zeros(1, length(for_analysis));
sd_ff = zeros(1, length(for_analysis));
sd_r2star = zeros(1, length(for_analysis));

for i = 1:length(for_analysis)
    mean_t1(i) = for_analysis(i).metrics.mean_roi_t1;
    mean_t1_remote(i) = for_analysis(i).metrics.mean_remote_t1;
    mean_ff(i) = for_analysis(i).metrics.mean_roi_ff;
    mean_r2star(i) = for_analysis(i).metrics.mean_roi_r2star;
    name_cell{i} = for_analysis(i).Name;
    
    sd_t1(i) = for_analysis(i).metrics.sd_roi_t1;
    sd_t1_remote(i) = for_analysis(i).metrics.sd_remote_t1;
    sd_ff(i) = for_analysis(i).metrics.sd_roi_ff;
    sd_r2star(i) = for_analysis(i).metrics.sd_roi_r2star;
end

figure();
X = 1:length(for_analysis);
errorbar(X, mean_t1, sd_t1);
hold on;
%plot(X, mean_t1_remote);
xticks(X);
xticklabels(name_cell);
yyaxis right; errorbar(X, mean_ff, sd_ff);