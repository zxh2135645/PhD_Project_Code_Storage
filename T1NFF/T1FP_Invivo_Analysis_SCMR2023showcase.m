clear all;
close all;
addpath('../function/');
addpath('../AHA16Segment/');
base_dir = uigetdir;
contour_glob = glob(cat(2, base_dir, '/ContourData/*'));
%Names = ExtractNames(contour_glob);
Names = {'Merry', 'Ryn', 'Mojave', 'Sahara', 'ZZ', 'Tina', 'Sunny', 'Queenie', 'Hope', 'Gobi', 'Felicity', 'Evelyn', '18D15', '18D16', '11D05', '11D26', '11D33'};
time_points = {'7D', '8WK', '12WK', '14WK', '4MO', '6MO', '9MO', '1YR', '15YR'};

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
data_save_dir = cat(2, base_dir, '/data/');

label_t1 = sequence_label{1};
label_t2star = sequence_label{2};
label_lge = sequence_label{3};

%% Main Body
for n = 4:4
    %n = 1:1
    name = Names{n};
    name_save_dir = cat(2, save_dir, name);
    if ~exist(name_save_dir, 'dir')
        mkdir(name_save_dir);
    end
    name_data_save_dir = cat(2, data_save_dir, name);
    if ~exist(name_data_save_dir, 'dir')
        mkdir(name_data_save_dir);
    end

    for tp = 8:8
        % tp = 9:9
        time_point = time_points{tp};
        tp_dir = cat(2, base_dir, '/ContourData/',  name, '/', name, '_', time_point,  '/');
        if ~exist(tp_dir, 'dir')
            disp(cat(2, 'No folder at: ', name, ' ', time_point));
        else
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
            
            % LGE
            myo_glob = glob(cat(2, tp_dir, label_lge, '/', anatomy_label{5}, '/*'));
            roi_glob = glob(cat(2, tp_dir, label_lge, '/',anatomy_label{3}, '/*'));
            remote_glob = glob(cat(2, tp_dir, label_lge, '/',anatomy_label{6}, '/*'));

            load(cat(2, tp_dir, label_lge, '/', label_lge, '_vol_img_3D.mat'));
            load(myo_glob{1});
            load(roi_glob{1});
            load(remote_glob{1});
            load(cat(2, tp_dir, label_lge, '/', label_lge, '_SliceLoc.mat'));
            [slc_array_t1, idx_reordered] = sort(slc_array);

            roi_in_myo_lge = mask_myocardium_3D .* freeROIMask_3D;
            remote_in_myo_lge = mask_myocardium_3D .* myoRefMask_3D;
            roi_lge = roi_in_myo_lge .* vol_img_3D;
            remote_lge = remote_in_myo_lge .* vol_img_3D;
            lge = vol_img_3D;
            myo_lge = mask_myocardium_3D;

            roi_in_myo_lge = roi_in_myo_lge(:,:,idx_reordered);
            remote_in_myo_lge = remote_in_myo_lge(:,:,idx_reordered);
            roi_lge = roi_lge(:,:,idx_reordered);
            remote_lge = remote_lge(:,:,idx_reordered);
            lge = lge(:,:,idx_reordered);
            myo_lge = myo_lge(:,:,idx_reordered);
%%
            slc = 3;
            figure();
            subplot(2,2,1);
            imagesc(lge(:,:,slc));axis image;
            subplot(2,2,2);
            imagesc(t1(:,:,slc));axis image;
            subplot(2,2,3);
            imagesc(ff(:,:,slc));axis image; caxis([0 40]);
            subplot(2,2,4);
            imagesc(r2star(:,:,slc));axis image; caxis([0 150]);

            se =  strel('disk',1);
            roi_lge = imerode(roi_in_myo_lge(:,:,slc), se);
            figure();
            subplot(1,2,1);
            imagesc(lge(:,:,slc).*roi_lge);axis image;
            subplot(1,2,2);
            imagesc(lge(:,:,slc).*roi_in_myo_lge(:,:,slc).*roi_lge);axis image;

            %%
            figure(); imagesc(t1(:,:,slc));axis image; caxis([0 2000]);
            colorbar;
            colormap(brewermap([],'*PuBuGn')); axis off;
%%
            figure(); imagesc(ff(:,:,slc));axis image; caxis([0 40]);
            colormap(brewermap([],'*RdYlBu')); axis off;
            figure(); imagesc(r2star(:,:,slc));axis image; caxis([0 200]);
            colormap(brewermap([],'*RdYlBu')); axis off;
%%
            slc = 3;
            figure(); imagesc(lge(:,:,slc));axis image;
            colorbar; axis off;
            colormap gray;
            
        end
    end
end


%% Merry scatter plot

ff_vec = [35.20333333; 22.5075; 0.768571429; 28.2675; 0.62]
r2star_vec = [32.62333333; 17.05; 30.29; 94.60125; 21.88]
figure();
h = scatter(ff_vec, r2star_vec, 72, 'LineWidth', 2);
set(h, 'FontSize', 12);