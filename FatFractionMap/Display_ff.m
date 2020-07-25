clear all;
close all;
clc
current_dir = pwd;
%% 
addpath('..\function\');
base_dir = GetFullPath(cat(2, current_dir, '\..\..\T1_Fat_Project\Data\'));

name_glob = glob(cat(2, base_dir, '*'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'\');
    Names{i} = strings{end-1};
end

sequence_label = {'T1', 'LGE', 'T2star', 'FF'};
anatomy_label = {'Heart', 'Myocardium', 'excludeArea', 'MyoReference', 'noReflowArea', 'freeROI'};
output_label = {'Mask'};

time_label = {'0D_occl', '0D', '1D', '3D', '7D', '21D', '28D', '8WK', '6MO', '1YR', '15YR'};


for n = 1:length(Names)
    name = Names{n};
    for tp = 1:length(time_label)
        
        time_point = time_label{tp};
        
        % where is orig_img.mat generated? From drawMyo.m
        img = load(cat(2, base_dir, name, '\', name, '_', time_point, '\orig_img.mat'));
        
        ff_glob = glob(cat(2, base_dir, name, '\', name, '_', time_point, '\FF\*'));

        if ~isempty(ff_glob)
            [list_to_read, order_to_read] = NamePicker(ff_glob, 1);
        end
        
        clear ff
        for i = 1:length(order_to_read)
            f = list_to_read{order_to_read(i)};
            load(f, 'fwmc_ff');
            ff(:,:,i) = fwmc_ff;
        end
        
        
        t1_label = sequence_label{1};
        ff_label = sequence_label{3};
        t1_img = img.orig_img.(t1_label);
        t1_mask_myo = load(cat(2, base_dir, name, '\', name, '_', time_point, '\Mask\', t1_label, '_Myocardium.mat'));
        ff_mask_myo = load(cat(2, base_dir, name, '\', name, '_', time_point, '\Mask\', ff_label, '_Myocardium.mat'));
        
        for i = 1:length(order_to_read)
            t1_edge = edge(t1_mask_myo.myo_mask(:,:,i), 'Canny');
            C = regionprops(t1_mask_myo.myo_mask(:,:,i));
            t1_cropped = CropAroundHeart(C.Centroid, t1_img(:,:,i));
            t1_myo_cropped = CropAroundHeart(C.Centroid, t1_edge);
            
            C = regionprops(ff_mask_myo.myo_mask(:,:,i));
            ff_cropped = CropAroundHeart(C.Centroid, ff(:,:,i));
            ff_myo_cropped = CropAroundHeart(C.Centroid, ff_mask_myo.myo_mask(:,:,i));

            figure('Position', [100 0 1600 1600]);
            t1_rgb = zeros(size(t1_cropped, 1), size(t1_cropped, 2), 3, 'uint8');
            t1_cropped = t1_cropped / max(t1_cropped(:)) * 255;

            t1_rgb(:,:,1) = t1_cropped + t1_myo_cropped*255;
            t1_rgb(:,:,2) = t1_cropped;
            t1_rgb(:,:,3) = t1_cropped;

            subplot(1,2,1); imagesc(t1_rgb); axis image;
            title('T1 MAP', 'FontSize', 24);
            
            ff_img = ff_cropped .* ff_myo_cropped;
            
            subplot(1,2,2); imagesc(ff_img); axis image;
            title('FF MAP', 'FontSize', 24);
            caxis([0 50]); colorbar;
            
            save_dir = GetFullPath(cat(2, base_dir, name, '\img\'));
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end
            
            fsave = GetFullPath(cat(2, save_dir, name, '_', time_point, '_slc', num2str(i), '.png'));
            saveas(gcf, fsave);
        end
        
        close all;
        ffsave = GetFullPath(cat(2, save_dir, '..\', name, '_', time_point, '\ff_map.mat'));
        save(ffsave, 'ff');
        
    end
end