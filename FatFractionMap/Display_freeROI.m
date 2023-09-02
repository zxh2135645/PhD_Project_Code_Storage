clear all;
close all;
clc;
current_dir = pwd;
%% 
addpath('..\function\');
base_dir = GetFullPath(cat(2, current_dir, '\..\..\Data\Diane\ContourData\'));
img_dir = GetFullPath(cat(2, base_dir, '../../Diane/ResultsFolder_180718/'));

name_glob = glob(cat(2, base_dir, '*'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'\');
    Names{i} = strings{end-1};
end

sequence_label = {'MultiEcho'};
anatomy_label = {'Heart', 'Myocardium', 'excludeArea', 'MyoReference', 'noReflowArea', 'freeROI'};
output_label = {'MultiEcho'};

name_check = 'Queenie_Final_28May2018';
te_array = [2.52, 4.46, 6.34, 8.22, 10.1, 11.98, 13.86, 15.74];

ff_glob = glob(cat(2, img_dir, name, '\', 'T2_multiecho_sax*_bright*.mat'));
idx_array = [4 5 6 7 8 9 10 1];
for la = 1:length(sequence_label)
    for i = 1:length(Names)
        name = Names{i};
        disp(name)
        label = sequence_label{la};
        labelo = label;
        
        % Assume images must be present
        load(cat(2, base_dir, name, '/', labelo, '/', label, '_vol_img_4D.mat'));
        if exist('vol_img_4D', 'var')
            img = vol_img_4D;
        else
            img = sig_vol_img_4D;
        end
        
        mask_myocardium_4D = [];
        freeROIMask_4D = [];
        freeROIMask_3D = [];
        
        if length(ls(cat(2, base_dir, name, '/', labelo, '/', anatomy_label{2}))) > 2
            load(cat(2, base_dir, name, '/', labelo, '/', anatomy_label{2}, '/', 'mask_myocardium.mat'));
        end
        if length(ls(cat(2, base_dir, name, '/', labelo, '/', anatomy_label{6}))) > 2
            load(cat(2, base_dir, name, '/', labelo, '/', anatomy_label{6}, '/', 'freeROI.mat'));
        end
        

        img_size = size(img);
        
        fig = figure('Renderer', 'painters', 'Position', [100 100 1800 1200]);
        %dim = [.2 .5 .3 .3];
        for tp = 1:img_size(3)
            %str = 'Straight Line Plot from 1 to 10';
            %annotation('textbox',dim,'String',str,'FitBoxToText','on');
            for slc = 1:img_size(4)
                n = ceil(sqrt(img_size(3)));
                subplot(n,n,slc);
                myo =  mask_myocardium_4D(:,:,tp,slc);
                C = regionprops(myo);
                img_cropped = CropAroundHeart(C.Centroid, img(:,:,tp,slc));
                myo_cropped = CropAroundHeart(C.Centroid, myo);
                freeROIMask_cropped = CropAroundHeart(C.Centroid, freeROIMask_4D(:,:,1,slc));
                myo_edge = edge(myo_cropped, 'Canny');
                %imagesc(img_cropped); axis image; colormap bone;
                title(cat(2, 'TE = ', num2str(te_array(tp)), ' ms'))
                %imagesc(img(:,:,1,slc) .* mask_myocardium_4D(:,:,1,slc)); axis image;
                
                img_cropped = img_cropped / max(img_cropped(:)) * 255;
                img_rgb = zeros(size(img_cropped, 1), size(img_cropped,2), 3, 'uint8');
                img_rgb(:,:,1) = img_cropped(:,:) + 255*myo_edge;
                img_rgb(:,:,2) = img_cropped(:,:) + 255*freeROIMask_cropped;
                img_rgb(:,:,3) = img_cropped(:,:);
                imagesc(img_rgb); axis image;
                
            end
            %movieVector(tp) = getframe(fig);
        end
        
        %fname = cat(2, base_dir, name, '/multiecho.avi');
        %FrameRate = 1;
        %RecordVideo(movieVector, fname, FrameRate)
        fig = figure('Renderer', 'painters', 'Position', [100 100 1800 1200]);
        %dim = [.2 .5 .3 .3];
            %str = 'Straight Line Plot from 1 to 10';
            %annotation('textbox',dim,'String',str,'FitBoxToText','on');
            tp = 1;
            for slc = 1:img_size(4)
                ff = load(ff_glob{idx_array(slc)});
                n = ceil(sqrt(img_size(3)));
                subplot(n,n,slc);
                myo =  mask_myocardium_4D(:,:,tp,slc);
                C = regionprops(myo);
                ff_cropped = CropAroundHeart(C.Centroid, ff.fwmc_ff);
                myo_cropped = CropAroundHeart(C.Centroid, myo);
                freeROIMask_cropped = CropAroundHeart(C.Centroid, freeROIMask_4D(:,:,1,slc));
                myo_edge = edge(myo_cropped, 'Canny');
                %imagesc(img_cropped); axis image; colormap bone;
                ff_cropped = ff_cropped .* myo_cropped;
                %imagesc(img(:,:,1,slc) .* mask_myocardium_4D(:,:,1,slc)); axis image;
                imagesc(ff_cropped); axis image; caxis([-5 25]); colorbar;
                title(cat(2, 'Slice ', num2str(slc)));
%               
%                 ff_cropped(ff_cropped > 100) = 100;
%                 ff_cropped(ff_cropped < 0) = 0;
%                 %ff_cropped = ff_cropped .* myo_cropped;
%                 ff_cropped = ff_cropped / max(ff_cropped(:)) * 255;
%                 img_rgb = zeros(size(ff_cropped, 1), size(ff_cropped,2), 3, 'uint8');
%                 img_rgb(:,:,1) = ff_cropped(:,:) + 255*myo_edge;
%                 img_rgb(:,:,2) = ff_cropped(:,:) + 255*freeROIMask_cropped;
%                 img_rgb(:,:,3) = ff_cropped(:,:);
%                 imagesc(img_rgb); axis image;
                
            end
            
            % Assume images must be present
            % LGE MAG
            lge = load(cat(2, base_dir, name, '/LGE/MAG_vol_img_3D.mat'));
            if length(ls(cat(2, base_dir, name, '/LGE/', anatomy_label{2}))) > 2
                lge_myo = load(cat(2, base_dir, name, '/LGE/', anatomy_label{2}, '/', 'mask_myocardium.mat'));
            end
            t1 = load(cat(2, base_dir, name, '/T1/T1Map_vol_img_3D.mat'));
            if length(ls(cat(2, base_dir, name, '/T1/', anatomy_label{2}))) > 2
                t1_myo = load(cat(2, base_dir, name, '/T1/', anatomy_label{2}, '/', 'mask_myocardium.mat'));
            end
            
            if length(ls(cat(2, base_dir, name, '/T1/', anatomy_label{6}))) > 2
                load(cat(2, base_dir, name, '/T1/', anatomy_label{6}, '/', 'freeROI.mat'));
            end
            %fig_t1 = figure('Renderer', 'painters', 'Position', [100 100 1800 1200]);
            fig_lge = figure('Renderer', 'painters', 'Position', [100 100 1800 1200]);
            
            for slc = 1:img_size(4)
                subplot(n,n,slc)
                myo = lge_myo.mask_myocardium_3D(:,:,slc);
                C = regionprops(myo);
                myo_cropped = CropAroundHeart(C.Centroid, myo);
                myo_edge = edge(myo_cropped,'Canny');
                
                img_cropped = CropAroundHeart(C.Centroid, lge.vol_img_3D(:,:,slc));
                img_cropped = img_cropped / max(img_cropped(:));
                img_rgb = zeros(size(img_cropped, 1), size(img_cropped, 2), 3);
                img_rgb(:,:,1) = img_cropped + myo_edge/2;
                img_rgb(:,:,2) = img_cropped;
                img_rgb(:,:,3) = img_cropped;
                imagesc(img_rgb); axis image;
                
            end
            
            % T1 Map            
            fig_t1 = figure('Renderer', 'painters', 'Position', [100 100 1800 1200]);
            
            for slc = 1:img_size(4)
                subplot(n,n,slc)
                myo = t1_myo.mask_myocardium_3D(:,:,slc);
                C = regionprops(myo);
                myo_cropped = CropAroundHeart(C.Centroid, myo);
                myo_edge = edge(myo_cropped,'Canny');
                
                img_cropped = CropAroundHeart(C.Centroid, t1.vol_img_3D(:,:,slc));
                img_cropped = img_cropped / max(img_cropped(:));
                img_rgb = zeros(size(img_cropped, 1), size(img_cropped, 2), 3);
                img_rgb(:,:,1) = img_cropped + myo_edge/2;
                img_rgb(:,:,2) = img_cropped;
                img_rgb(:,:,3) = img_cropped;
                imagesc(img_rgb); axis image;
                
            end
            
            % T1 Map + freeROI        
            fig_t1 = figure('Renderer', 'painters', 'Position', [100 100 1800 1200]);
            
            for slc = 1:img_size(4)
                subplot(n,n,slc)
                myo = t1_myo.mask_myocardium_3D(:,:,slc);
                C = regionprops(myo);
                myo_cropped = CropAroundHeart(C.Centroid, myo);
                myo_edge = edge(myo_cropped,'Canny');
                
                freeROIMask_cropped = CropAroundHeart(C.Centroid, freeROIMask_3D(:,:,slc));
                
                img_cropped = CropAroundHeart(C.Centroid, t1.vol_img_3D(:,:,slc));
                img_cropped = img_cropped / max(img_cropped(:));
                img_rgb = zeros(size(img_cropped, 1), size(img_cropped, 2), 3);
                img_rgb(:,:,1) = img_cropped + myo_edge/2;
                img_rgb(:,:,2) = img_cropped + freeROIMask_cropped/2;
                img_rgb(:,:,3) = img_cropped;
                imagesc(img_rgb); axis image;
                
            end
    end
end