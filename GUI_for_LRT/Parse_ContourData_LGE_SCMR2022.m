clear all;
close all;
clc;
current_dir = pwd;
% Dog data configuration for Ting 06/18/2021
%
addpath('../function/');
addpath('../T1NFF/');
base_dir = uigetdir;

name_glob = glob(cat(2, base_dir, '/*'));
Names = cell(length(name_glob), 1);
for i = 1:length(name_glob)
    strings = strsplit(name_glob{i},'/');
    Names{i} = strings{end-1};
end

%sequence_label = {'MAG'};
sequence_label = {'LRT_LGE'};
% sequence_label = {'DICOM_E', 'DICOM_L', 'PSIR', 'MAG'};
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
% output_label = {'MAG'};
output_label = {'LRT_LGE'};
%% Plot MVO and Infarct size
wx = 64;
wy = 64;
for n = 11:length(Names)
% for n = 9:9
    name = Names{n}
    for t = 1:length(time_label)
        %mvo_perc = zeros(length(sequence_label), 1);
        %infarct_perc = zeros(length(sequence_label), 1);
        
        dst_dir = cat(2, base_dir, '/', name,'/', name, '_', time_label{t});

        if exist(cat(2, dst_dir, '/ContourData/'))
            time = time_label{t}
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
                    if ~strcmp(label, 'LRT_LGE')
                        mask_myocardium_tight = zeros(size(mask_myocardium_3D{1},1), size(mask_myocardium_3D{1},2), length(idx_sorted));
                        vol_img_tight = zeros(size(mask_myocardium_3D{1},1), size(mask_myocardium_3D{1},2), length(idx_sorted));
                        noReflowMask_tight = zeros(size(mask_myocardium_3D{1},1), size(mask_myocardium_3D{1},2), length(idx_sorted));
                        freeROIMask_tight = zeros(size(mask_myocardium_3D{1},1), size(mask_myocardium_3D{1},2), length(idx_sorted));
                    else
                        mask_myocardium_tight = zeros(size(mask_myocardium_3D,1), size(mask_myocardium_3D,2), length(idx_sorted));
                        vol_img_tight = zeros(size(mask_myocardium_3D,1), size(mask_myocardium_3D,2), length(idx_sorted));
                        noReflowMask_tight = zeros(size(mask_myocardium_3D,1), size(mask_myocardium_3D,2), length(idx_sorted));
                        freeROIMask_tight = zeros(size(mask_myocardium_3D,1), size(mask_myocardium_3D,2), length(idx_sorted));
                    end

                    if ~strcmp(label, 'LRT_LGE')
                        for ii = 1:length(mask_myocardium_3D)
                            mask_myocardium_tight(:,:,ii) = mask_myocardium_3D{idx_sorted(ii)};
                            vol_img_tight(:,:,ii) = vol_img_3D{idx_sorted(ii)};
                            noReflowMask_tight(:,:,ii) = noReflowMask_3D{idx_sorted(ii)};
                            freeROIMask_tight(:,:,ii) = freeROIMask_3D{idx_sorted(ii)};
                        end
                    else
                       for ii = 1:size(mask_myocardium_3D,3)
                            mask_myocardium_tight(:,:,ii) = mask_myocardium_3D(:,:,idx_sorted(ii));
                            vol_img_tight(:,:,ii) = vol_img_3D(:,:,idx_sorted(ii));
                            noReflowMask_tight(:,:,ii) = noReflowMask_3D(:,:,idx_sorted(ii));
                            freeROIMask_tight(:,:,ii) = freeROIMask_3D(:,:,idx_sorted(ii));
                        end
                    end

                    noReflowMask_tight = mask_myocardium_tight .* noReflowMask_tight;
                    freeROIMask_tight = mask_myocardium_tight .* freeROIMask_tight;
                   
                    mvo_perc = zeros(length(slc_array), 1);
                    infarct_perc = zeros(length(slc_array), 1);
                    
                    for slc = 1:size(mask_myocardium_tight, 3)
                        mvo_perc(slc) = numel(nonzeros(noReflowMask_tight(:,:,slc))) ./ numel(nonzeros(mask_myocardium_tight(:,:,slc)));
                        infarct_perc(slc) = numel(nonzeros(freeROIMask_tight(:,:,slc))) ./ numel(nonzeros(mask_myocardium_tight(:,:,slc)));
                    end
                    mvo_perc*100
                    infarct_perc*100
                    mvo_perc_total = numel(nonzeros(noReflowMask_tight)) ./ numel(nonzeros(mask_myocardium_tight))
                    infarct_perc_total = numel(nonzeros(freeROIMask_tight)) ./ numel(nonzeros(mask_myocardium_tight))
                    pause;
                end
            end
        end
    end
end
%%        
% subplot full image
figure('Position',[100, 100, 900, 1000]);
x = ceil(sqrt(length(slc_array)));
[ha, pos] = tight_subplot(x,x,[.01 0],[.01 .01],[.01 .01]);
for slc = 1:length(slc_array)
    axes(ha(slc));
    imagesc(vol_img_tight(:,:,slc)); axis image;
    colormap gray; axis off;
    title(cat(2, 'SA ', num2str(slc)));
end
saveas(gcf, cat(2, save_dir, label, '_orig.png'));


centroids = cell(length(slc_array), 1);
myo_crop = zeros(wx,wy,length(slc_array));
for slc = 1:length(slc_array)
    s = regionprops(mask_myocardium_tight(:,:,slc), 'centroid');
    centroids{slc} = round(s.Centroid);
    centroid = centroids{slc};
    myo_crop(:,:,slc) = imcrop(mask_myocardium_tight(:,:,slc), [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
end

% subplot cropped image
figure('Position',[100, 100, 900, 1000]);
x = ceil(sqrt(length(slc_array)));
[ha, pos] = tight_subplot(x,x,[.01 0],[.01 .01],[.01 .01]);
for slc = 1:length(slc_array)
    centroid = centroids{slc};
    axes(ha(slc));
    temp_crop = imcrop(vol_img_tight(:,:,slc), [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
    imagesc(temp_crop); axis image;
    colormap gray; axis off;
    title(cat(2, 'SA ', num2str(slc)));
end
saveas(gcf, cat(2, save_dir, label, '_cropped.png'));

% noReflowMask
figure('Position',[100, 100, 900, 1000]);
x = ceil(sqrt(length(slc_array)));
[ha, pos] = tight_subplot(x,x,[.01 0],[.01 .01],[.01 .01]);

for slc = 1:length(slc_array)
    centroid = centroids{slc};
    axes(ha(slc));
    temp_crop = imcrop(noReflowMask_tight(:,:,slc), [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
    freeROI_crop = imcrop(freeROIMask_tight(:,:,slc), [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
    comp_temp = temp_crop+myo_crop(:,:,slc)+freeROI_crop;
    %comp_temp(comp_temp == 0) = nan;
    imagesc(comp_temp); axis image;
    colormap jet; axis off;
    mvo = numel(nonzeros(noReflowMask_tight(:,:,slc))) ./ numel(nonzeros(mask_myocardium_tight(:,:,slc)));
    title(cat(2, 'MVO = ', num2str(mvo, '%.2f')));
end
saveas(gcf, cat(2, save_dir, label, '_masks.png'));

% Composite
figure('Position',[100, 100, 900, 1000]);
x = ceil(sqrt(length(slc_array)));
[ha, pos] = tight_subplot(x,x,[.01 0],[.01 .01],[.01 .01]);

for slc = 1:length(slc_array)
    centroid = centroids{slc};
    axes(ha(slc));
    temp_crop = imcrop(noReflowMask_tight(:,:,slc), [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
    freeROI_crop = imcrop(freeROIMask_tight(:,:,slc), [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
    img_crop = imcrop(vol_img_tight(:,:,slc), [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);

    img_rgb = zeros(wx, wy, 3);
    comp_temp = temp_crop+myo_crop(:,:,slc)+freeROI_crop;
    comp_temp(comp_temp == 0) = nan;
    comp = double2rgb(comp_temp, jet);
    img_gray = double2rgb(img_crop, gray);
    imagesc(comp+img_gray); axis image;
    axis off;
    mvo = numel(nonzeros(noReflowMask_tight(:,:,slc))) ./ numel(nonzeros(mask_myocardium_tight(:,:,slc))) * 100;
    infarct = numel(nonzeros(freeROIMask_tight(:,:,slc))) ./ numel(nonzeros(mask_myocardium_tight(:,:,slc))) * 100;
    title(cat(2, 'MVO = ', num2str(mvo, '%.1f'), '%, Infarct = ', num2str(infarct, '%.1f'), '%'));
end
saveas(gcf, cat(2, save_dir, label, '_composite.png'));

close all;

