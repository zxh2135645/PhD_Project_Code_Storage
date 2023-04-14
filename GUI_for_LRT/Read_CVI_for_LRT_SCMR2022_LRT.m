clear all;
close all;
clc;
current_dir = pwd;
% Parse cvi from LRT images (DICOM_LGE_CVI)
%
addpath('../function/');
addpath('../T1NFF/');
base_dir = uigetdir;

name_glob = glob(cat(2, base_dir, '/*/*/*.cvi42wsx.xml'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'/');
    Names{i} = strings{end-1};
end

% sequence_label = {'PSIR', 'MAG'};
% sequence_label = {'DICOM_E', 'DICOM_L', 'PSIR', 'MAG'};
sequence_label = {'DICOM_LGE_CVI'};
time_label = {'D6', 'D8', 'WK8', 'WK8+2'};

names_to_rule_out = {};
RuleOutLabel = NameRuleOutFunc(Names, names_to_rule_out);
Names = Names(RuleOutLabel == 0);
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
% output_label = {'PSIR', 'MAG'};
output_label = {'LRT_LGE'};
%% Main Body
%for n = 1:length(Names)
for n = 19:20
    name = Names{n};
    % for t = 1:length(time_label)
    %cvi_glob = glob(cat(2, name_glob{n}, name, '_', time_label{t}, '/', name, '*.cvi42wsx.xml'));
    f_name = name_glob{n};
    if isempty(f_name)
        disp(cat(2, 'NO CVI XML: ', name));
    else
        cvi42wsx = char(f_name);
        con_cell = cell(0);
        for xml_ind = 1:size(cvi42wsx, 1)
            % As Yinyin reported, this one has two xml file because T1 and LGE are shown in different cvi42 directory
            % Thus, there are two different files
            con_cell{end+1} = CMR42ContourReader(cvi42wsx(xml_ind,:));
        end
        % Iterate through MAG, PSIR and LGE
        for con_idx = 1:length(con_cell)
            % A different label for Exvivo
            con = con_cell{con_idx};

            for ll = 1:length(sequence_label)
                label = sequence_label{ll};
                label_out = output_label{ll};
                % Check if the dstFolder
                dstFolder = GetFullPath(cat(2, f_name, '/../ContourData/', label_out, '/'));
                % dstFolder = cat(2, f_name, '/', '/ContourData/', label_out, '/');

                if strcmp(label, 'DICOM_E') || strcmp(label, 'DICOM_L')
                    dicom_glob = glob(cat(2, base_dir, '/', name, '/Data/LRT_recon_Seg15/', label, '/'));
                else
                    % dicom_glob = glob(cat(2, base_dir, '/', name, '/Data/OrigData/*', label, '_[0123456789]*/'));
                    dicom_glob = glob(GetFullPath(cat(2, f_name, '/../DICOM_LGE_CVI/')));
                end
                

                ReadCVI_Workflow_Cell_Func(con, dicom_glob, dstFolder, dicom_fields);
            end

        end
    end
end

%%
new_mask_myo = [];
new_vol_img = [];
new_freeROI = []
ct = 0;
for i = 1:84
    temp = mask_myocardium_3D(:,:,i);
    if any(temp(:))
        ct = ct + 1;
        new_mask_myo(:,:,ct) = mask_myocardium_3D(:,:,i);
        new_vol_img(:,:,ct) = vol_img_3D(:,:,i);
        new_freeROI(:,:,ct) = freeROIMask_3D(:,:,i);
    end
end

figure();
for i = 1:size(new_mask_myo,3)
    subplot(3,3,i);
    % imagesc(new_mask_myo(:,:,i).*new_vol_img(:,:,i)); axis image;
    imagesc(new_freeROI(:,:,i).*new_vol_img(:,:,i)); axis image;
end

figure();
for i = 1:size(new_mask_myo,3)
    subplot(3,3,i);
    imagesc(new_mask_myo(:,:,i).*new_vol_img(:,:,i)); axis image;
    % imagesc(new_freeROI(:,:,i).*new_vol_img(:,:,i)); axis image;
end
%% Plot MVO and Infarct size
wx = 64;
wy = 64;
% for n = 1:length(Names)
for n = 4:4 
    name = Names{n}
    %mvo_perc = zeros(length(sequence_label), 1);
    %infarct_perc = zeros(length(sequence_label), 1);
    f_name = name_glob{n};
    % dst_dir = cat(2, base_dir, '/', name,'/', name, '_', time_label{t});
    % dst_dir = cat(2, base_dir, '/', name,'/', name, '_', time_label{t});
    dst_dir = GetFullPath(cat(2, f_name, '/../'));
    if exist(cat(2, dst_dir, '/ContourData/'))
        
        save_dir = cat(2, dst_dir, '/Results/');
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end

        for i = 1:length(sequence_label)

            if i > 3 % DICOM_L or DICOM_E
                %             label = sequence_label{i};
                %             contour_dir = cat(2, base_dir, '/', name, '/ContourData/', label, '/' );
                %
                %             load(cat(2, contour_dir, label, '_vol_img_3D.mat'));
                %             load(cat(2, contour_dir, 'freeROI/freeROI.mat'));
                %             load(cat(2, contour_dir, 'MyoReference/myoRef.mat'));
                %             load(cat(2, contour_dir, 'Myocardium/mask_myocardium.mat'));
                %             load(cat(2, contour_dir, 'noReflowArea/noReflow.mat'));
                %
                %             slc_array = [];
                %             for ii = 1:size(mask_myocardium_3D, 3)
                %                 myo_temp = mask_myocardium_3D(:,:,ii);
                %                 if any(myo_temp(:))
                %                     slc_array = [slc_array, ii];
                %                 end
                %             end
                %
                %             mask_myocardium_tight = mask_myocardium_3D(:,:,slc_array);
                %             vol_img_tight = vol_img_3D(:,:,slc_array);
                %             noReflowMask_tight = noReflowMask_3D(:,:,slc_array);
                %             freeROIMask_tight = freeROIMask_3D(:,:,slc_array);
                %
                %             noReflowMask_tight = mask_myocardium_tight .* noReflowMask_tight;
                %             freeROIMask_tight = mask_myocardium_tight .* freeROIMask_tight;
                %
                %             mvo_perc(i) = numel(nonzeros(noReflowMask_tight)) ./ numel(nonzeros(mask_myocardium_tight));
                %             infarct_perc(i) = numel(nonzeros(freeROIMask_tight)) ./ numel(nonzeros(mask_myocardium_tight));
            else
                label = sequence_label{i};
                contour_dir = cat(2, dst_dir, '/ContourData/', label, '/');
                load(cat(2, contour_dir, label, '_vol_img_3D.mat'));
                load(cat(2, contour_dir, 'freeROI/freeROI.mat'));
                load(cat(2, contour_dir, 'MyoReference/myoRef.mat'));
                load(cat(2, contour_dir, 'Myocardium/mask_myocardium.mat'));
                load(cat(2, contour_dir, 'noReflowArea/noReflow.mat'));
                load(cat(2, contour_dir, label, '_SliceLoc.mat'));

                [slc_array_sorted, idx_sorted] = sort(slc_array);
                mask_myocardium_tight = [];
                vol_img_tight = [];
                noReflowMask_tight = [];
                freeROIMask_tight = [];
                ct = 0;

                for ii = 1:size(mask_myocardium_3D, 3)
                    temp = mask_myocardium_3D(:,:,ii);
                    if any(temp(:))
                        ct = ct + 1;
                        mask_myocardium_tight(:,:,ct) = mask_myocardium_3D(:,:,ii);
                        vol_img_tight(:,:,ct) = vol_img_3D(:,:,ii);
                        noReflowMask_tight(:,:,ct) = noReflowMask_3D(:,:,ii);
                        freeROIMask_tight(:,:,ct) = freeROIMask_3D(:,:,ii);
                    end
                end

                noReflowMask_tight = mask_myocardium_tight .* noReflowMask_tight;
                freeROIMask_tight = mask_myocardium_tight .* freeROIMask_tight;

                mvo_perc = zeros(size(noReflowMask_tight,3), 1);
                infarct_perc = zeros(size(noReflowMask_tight,3), 1);

                for slc = 1:size(mask_myocardium_tight, 3)
                    mvo_perc(slc) = numel(nonzeros(noReflowMask_tight(:,:,slc))) ./ numel(nonzeros(mask_myocardium_tight(:,:,slc)));
                    infarct_perc(slc) = numel(nonzeros(freeROIMask_tight(:,:,slc))) ./ numel(nonzeros(mask_myocardium_tight(:,:,slc)));
                end
                mvo_perc
                infarct_perc
                mvo_perc_total = numel(nonzeros(noReflowMask_tight)) ./ numel(nonzeros(mask_myocardium_tight))
                infarct_perc_total = numel(nonzeros(freeROIMask_tight)) ./ numel(nonzeros(mask_myocardium_tight))
                pause;
            end
        end
    end
end

%% Draw insertion point (Can be skipped if there is no data added)
for n = 1:length(Names)
% for n = 6:6
    name = Names{n}
    for t = 1:length(time_label)
        %mvo_perc = zeros(length(sequence_label), 1);
        %infarct_perc = zeros(length(sequence_label), 1);
        time = time_label{t}
        % dst_dir = cat(2, base_dir, '/', name,'/', name, '_', time_label{t});
        f_name = name_glob{n};
        dst_dir = GetFullPath(cat(2, f_name, '/../'));

        if exist(cat(2, dst_dir, '/ContourData/'))
            save_dir = cat(2, dst_dir, '/Results/');
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end

            for i = 1:length(sequence_label)
                label = sequence_label{i};
                contour_dir = cat(2, dst_dir, '/ContourData/', label, '/');
                load(cat(2, contour_dir, label, '_vol_img_3D.mat'));
                load(cat(2, contour_dir, label, '_SliceLoc.mat'));
                load(cat(2, contour_dir, 'Myocardium/mask_myocardium.mat'));

                [slc_array_sorted, idx_sorted] = sort(slc_array);
                mask_myocardium_tight = [];
                vol_img_tight = [];
                ct = 0;

                for ii = 1:size(mask_myocardium_3D, 3)
                    temp = mask_myocardium_3D(:,:,ii);
                    if any(temp(:))
                        ct = ct + 1;
                        mask_myocardium_tight(:,:,ct) = mask_myocardium_3D(:,:,ii);
                        vol_img_tight(:,:,ct) = vol_img_3D(:,:,ii);
                        % noReflowMask_tight(:,:,ct) = noReflowMask_3D(:,:,ii);
                        % freeROIMask_tight(:,:,ct) = freeROIMask_3D(:,:,ii);
                    end
                end

                ref_pts_save = cat(2, contour_dir, '/ref_pts.mat');
                if ~exist(ref_pts_save, 'file')
                    figure();

                    x = zeros(size(vol_img_tight, 3),1);
                    y = zeros(size(vol_img_tight, 3),1);
                    for slc = 1:size(vol_img_tight, 3)
                        imagesc(vol_img_tight(:,:,slc)); % caxis([0 100]);
                        [x(slc),y(slc)] = getpts(gca);
                    end
                    ref_pts.x = x;
                    ref_pts.y = y;
                    save(ref_pts_save, 'ref_pts');
                else
                    load(ref_pts_save);
                    x = ref_pts.x;
                    y = ref_pts.y;
                end
            end
        end
    end
end

%% 
%% Transmurality analysis
addpath('../AHA16Segment/');
% for n = 1:length(Names)
for n = 4:4
    name = Names{n}
    % for t = 1:length(time_label)
    %mvo_perc = zeros(length(sequence_label), 1);
    %infarct_perc = zeros(length(sequence_label), 1);
    %time = time_label{t}
    f_name = name_glob{n};
    dst_dir = GetFullPath(cat(2, f_name, '/../'));

    if exist(cat(2, dst_dir, '/ContourData/'))
        save_dir = cat(2, dst_dir, '/Results/');
        if ~exist(save_dir, 'dir')
            mkdir(save_dir);
        end

        for i = 1:length(sequence_label)

            label = sequence_label{i};
            contour_dir = cat(2, dst_dir, '/ContourData/', label, '/');
            load(cat(2, contour_dir, label, '_vol_img_3D.mat'));
            load(cat(2, contour_dir, label, '_SliceLoc.mat'));
            load(cat(2, contour_dir, 'Myocardium/mask_myocardium.mat'));
            load(cat(2, contour_dir, 'Heart/mask_heart.mat'));
            load(cat(2, contour_dir, 'freeROI/freeROI.mat'));

            [slc_array_sorted, idx_sorted] = sort(slc_array);
            mask_myocardium_tight = [];
            vol_img_tight = [];
            mask_heart_tight = [];
            freeROIMask_tight = [];
            ct = 0;

            for ii = 1:size(mask_myocardium_3D, 3)
                temp = mask_myocardium_3D(:,:,ii);
                if any(temp(:))
                    ct = ct + 1;
                    mask_myocardium_tight(:,:,ct) = mask_myocardium_3D(:,:,ii);
                    vol_img_tight(:,:,ct) = vol_img_3D(:,:,ii);
                    % noReflowMask_tight(:,:,ct) = noReflowMask_3D(:,:,ii);
                    % freeROIMask_tight(:,:,ct) = freeROIMask_3D(:,:,ii);
                    mask_heart_tight(:,:,ct) = mask_heart_3D(:,:,ii);
                    freeROIMask_tight(:,:,ct) = freeROIMask_3D(:,:,ii);
                end
            end

            freeROIMask_tight = mask_myocardium_tight .* freeROIMask_tight;

            ref_pts_save = cat(2, contour_dir, '/ref_pts.mat');

            load(ref_pts_save);
            x = ref_pts.x;
            y = ref_pts.y;

            Segn = 50;
            BaseGroove = zeros(size(vol_img_tight, 3), 1);
            transmurality_array = zeros(size(vol_img_tight, 3), 1);

            for slc = 1:size(vol_img_tight, 3)
                blood_pool_tight = (mask_myocardium_tight)>0;
                C = regionprops(blood_pool_tight(:,:,slc));
                x_centroid = C.Centroid(2);
                y_centroid = C.Centroid(1);
                BaseGroove(slc) = atan2(x(slc) - x_centroid, y(slc) - y_centroid) * 180 / pi;

                [Segmentpix_roi, stats, Mask_Segn] = AHASegmentation(freeROIMask_tight(:,:,slc), mask_myocardium_tight(:,:,slc)>0, Segn, BaseGroove(slc));
                [Segmentpix_myo, stats, Mask_Segn] = AHASegmentation(vol_img_tight(:,:,slc), mask_myocardium_tight(:,:,slc)>0, Segn, BaseGroove(slc));

                mi_numel = zeros(Segn, 1);
                myo_numel = zeros(Segn, 1);
                for seg = 1:Segn
                    mi_numel(seg) = sum(Segmentpix_roi{seg});
                    myo_numel(seg) = length(Segmentpix_myo{seg});
                end
                transmurality_temp = (mi_numel ./ myo_numel)*100;
                transmurality_temp(transmurality_temp < 1) = [];
                if isempty(transmurality_temp)
                    transmurality_array(slc) = 0;
                else
                    transmurality_array(slc) = mean(transmurality_temp);
                end
            end
            transmurality_array
            transmurality_array_total = mean(transmurality_array(transmurality_array>1))
            pause;
        end
        % end
    end
end