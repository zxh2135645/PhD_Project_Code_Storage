close all;
clear all;

% T1 and FF analysis


addpath('../function/');
base_dir = '../T1_Fat_Project/';
contour_glob = glob(cat(2, base_dir, 'ContourData\*'));

Names = ExtractNames(contour_glob);

% labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};
time_points = {'0D_occl', '0D', '1D', '3D', '7D', '21D', '28D', '8WK', '6MO', '1YR', '15YR'};

% OutputPath = GetFullPath(cat(2, base_dir, 'ContourData\'));
% if ~exist(OutputPath, 'dir')
%     mkdir(OutputPath);
% end

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


sequence_label = {'T1', 'T2star'};
anatomy_label = {'BloodPool', 'excludeArea', 'freeROI', 'Heart', 'Myocardium', 'MyoReference', 'noReflowArea'}; 
name_check = 'Merry';
starting_point = find(strcmp(name_check, Names),1);

save_dir = cat(2, base_dir, 'img\');
%%
% What is the metrics? 
%for ll = 1:length(sequence_label)
% T1 Map
time_points = {'0D_occl', '0D', '1D', '3D', '7D', '21D', '28D', '8WK', '6MO', '1YR', '15YR'};
for ll = 1:1
    label = sequence_label{ll};
    
    for n = 1:starting_point
        name = Names{n};
        name_save_dir = cat(2, save_dir, name);
        if ~exist(name_save_dir, 'dir')
            mkdir(name_save_dir);
        end
        
        metrics = zeros(1, length(time_points));
        metrics_remote = zeros(1, length(time_points));
        metrics_sd = zeros(1, length(time_points));
        metrics_sd_remote = zeros(1, length(time_points));
        for tp = 1:length(time_points)
            %for tp = 1:length(time_points)
            time_point = time_points{tp};
            
            myo_glob = glob(cat(2, base_dir, 'ContourData\',  name, '\', name, '_', time_point,  '\', label, '\', anatomy_label{5}, '\*'));
            roi_glob = glob(cat(2, base_dir, 'ContourData\',  name, '\', name, '_', time_point,  '\', label, '\',anatomy_label{3}, '\*'));
            remote_glob = glob(cat(2, base_dir, 'ContourData\',  name, '\', name, '_', time_point,  '\', label, '\',anatomy_label{6}, '\*'));
            
            load(cat(2, base_dir, 'ContourData\',  name, '\', name, '_', time_point,  '\', label, '\', label, '_vol_img_3D.mat'));
            load(myo_glob{1});
            load(roi_glob{1});
            load(remote_glob{1});
            
            roi_in_myo = mask_myocardium_3D .* freeROIMask_3D;
            remote_in_myo = mask_myocardium_3D .* myoRefMask_3D;
            roi = roi_in_myo .* vol_img_3D;
            remote = remote_in_myo .* vol_img_3D;
            
            figure();
            for i = 1:size(roi, 3)
                subplot(2,2,i);
                imagesc(roi(:,:,i)); axis image;
                title([name, ' ', time_point, '  Slice = ', num2str(i)]);
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

            metrics(tp) = mean(nonzeros(roi));
            metrics_remote(tp) = mean(nonzeros(remote));
            metrics_sd(tp) = std(nonzeros(roi));
            metrics_sd_remote(tp) = std(nonzeros(remote));
        end
        
        figure('Position', [1000 500 900 600]);
        errorbar(metrics, metrics_sd, '-o', 'LineWidth', 2); xlabel('Time point'); ylabel('mean T1 (ms)');
        hold on;  errorbar(metrics_remote, metrics_sd_remote, '-o', 'LineWidth', 2);
        xticklabels(time_points);
        set(gca, 'FontSize', 16);
        grid on;
        legend({'MI', 'Remote'}, 'Location', 'SouthEast');
        
    end
end

%%
% FF Map
time_points = {'0D_occl', '0D', '1D', '3D', '7D', '21D', '28D', '8WK', '6MO', '1YR', '15YR'};
for ll = 2:2
    label = sequence_label{ll};
    
    for n = 1:starting_point
        name = Names{n};
                name_save_dir = cat(2, save_dir, name);
        if ~exist(name_save_dir, 'dir')
            mkdir(name_save_dir);
        end
        
        metrics = zeros(1, length(time_points));
        metrics_remote = zeros(1, length(time_points));
        metrics_sd = zeros(1, length(time_points));
        metrics_sd_remote = zeros(1, length(time_points));
        
        for tp = 1:length(time_points)
            %for tp = 1:length(time_points)
            time_point = time_points{tp};
            
            myo_glob = glob(cat(2, base_dir, 'ContourData\',  name, '\', name, '_', time_point,  '\', label, '\', anatomy_label{5}, '\*'));
            roi_glob = glob(cat(2, base_dir, 'ContourData\',  name, '\', name, '_', time_point,  '\', label, '\',anatomy_label{3}, '\*'));
            remote_glob = glob(cat(2, base_dir, 'ContourData\',  name, '\', name, '_', time_point,  '\', label, '\',anatomy_label{6}, '\*'));
            
%            ff_glob = glob(cat(2, base_dir, 'Data\', name, '\', name, '_', time_point, '\ff_map.mat'));
            
%             [list_to_read, order_to_read] = NamePicker(ff_glob, 1);
%            
%             clear ff
%             for i = 1:length(order_to_read)
%                 f = list_to_read{order_to_read(i)};
%                 load(f, 'fwmc_ff');
%                 ff(:,:,i) = fwmc_ff;
%             end
            load(cat(2, base_dir, 'Data\', name, '\', name, '_', time_point, '\ff_map.mat'));
            load(myo_glob{1});
            load(roi_glob{1});
            load(remote_glob{1});
            
            roi_in_myo = mask_myocardium_3D .* freeROIMask_3D;
            remote_in_myo = mask_myocardium_3D .* myoRefMask_3D;
            
            roi_ind = find(roi_in_myo == 1);
            remote_ind = find(remote_in_myo == 1);
            roi = roi_in_myo .* ff;
            remote = remote_in_myo .* ff;
            
            figure();
            for i = 1:size(roi, 3)
                subplot(2,2,i);
                imagesc(roi(:,:,i)); axis image; caxis([0 50]); colorbar;
                title([name, ' ', time_point, '  Slice = ', num2str(i)]);
            end
            saveas(gcf, cat(2, name_save_dir, '\', name, '_', time_point, '_MI_FF.png'));
            
            roi(roi < 0) = 0;
            roi(roi > 100) = 100;
            remote(remote < 0) = 0;
            remote(remote > 100) = 100;
            
            
            figure();
            histogram(roi(roi_ind),'Normalization', 'probability');xlabel('Fat Fraction (%)'); ylabel('Frequency');
            hold on;
            histogram(remote(remote_ind),'Normalization', 'probability');
            set(gca, 'FontSize', 16); title([name, ' ', time_point]);
            xlim([0 100]);
            legend({'MI', 'Remote'});
            saveas(gcf, cat(2, name_save_dir, '\', name, '_', time_point, '_Histogram_FF.png'));

            
            metrics(tp) = mean(roi(roi_ind));
            metrics_remote(tp) = mean(remote(remote_ind));
            metrics_sd(tp) = std(roi(roi_ind));
            metrics_sd_remote(tp) = std(remote(remote_ind));
        end
        
        figure('Position', [1000 500 900 600]);
        errorbar(metrics, metrics_sd, '-o', 'LineWidth', 2); xlabel('Time point'); ylabel('mean T1 (ms)');
        hold on; errorbar(metrics_remote, metrics_sd_remote, '-o', 'LineWidth', 2);
        xticklabels(time_points);
        set(gca, 'FontSize', 16);
        grid on;
        legend({'MI', 'Remote'}, 'Location', 'SouthEast');
        
    end
end