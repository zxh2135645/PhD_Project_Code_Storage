clear all;
close all;

addpath('.\function\')
base_dir = 'C:\Users\xz100\Documents\Data\T2star_SimulationPhantom\';
vec = @(x) x(:);

load(cat(2, base_dir, '\Cylinder_Phantom.mat'));

%% See the noise level
detection_analysis = struct;
snr_array = [0.0500    0.1000    0.1500    0.2000    0.2500    0.3000    0.3500    0.4000    0.4500];
res_array = [0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4]; % in mm
res_through_array = [2, 4, 6, 8]; % in mm
transmural_array = [0.025, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8];

    %% Should iterate over different transmurality
    k = 1;
    transmural = transmural_array(k);
    VObj = Phantom_shape_cell{k};
    sz = size(VObj.t2star);
    
    Nz = sz(3);
    Nx = 1024;
    Ny = 1024;
    
    [X, Y] = meshgrid(1:sz(1));
    [Xq, Yq] = meshgrid(linspace(1, sz(1), Nx));
    
    t2star_crop_interp = interp2(X, Y, VObj.t2star(:,:,1), Xq, Yq, 'spline');
    
    hemo_mask_gt = repmat(t2star_crop_interp, [1, 1, sz(3)]);
    hemo_mask_gt(abs(hemo_mask_gt) >= 0.01375) = 0; % Blood Pool
    hemo_mask_gt(abs(hemo_mask_gt) < 0.0025 & abs(hemo_mask_gt) > 0) = 0;
    hemo_mask_gt(abs(hemo_mask_gt) < 0.01375 & abs(hemo_mask_gt) >= 0.0075) = 2; % Remote
    hemo_mask_gt(abs(hemo_mask_gt) < 0.0075 & abs(hemo_mask_gt) >= 0.0025) = 3;    % Hemo
    
    se1 = strel('disk', 2);
    hemo_mask_gt = imerode(hemo_mask_gt,se1);

    myo_mask_gt = hemo_mask_gt > 0;
    hemo_mask_gt_binary = hemo_mask_gt > 2;
    remote_mask_gt_binary = ~hemo_mask_gt_binary;
    air_mask_gt = hemo_mask_gt == 0;
    remote_mask_gt_binary_naned = remote_mask_gt_binary(myo_mask_gt);

    se2 = strel('disk', 30);
    remote_roi_binary = imerode(remote_mask_gt_binary.* myo_mask_gt, se2);
    
    %clear Phantom_shape_cell
    accuracy_array = zeros(length(res_array), length(res_through_array), length(snr_array));
    sensitivity_array = zeros(length(res_array), length(res_through_array), length(snr_array));
    dice_array = zeros(length(res_array), length(res_through_array), length(snr_array));
    specificity_array = zeros(length(res_array), length(res_through_array), length(snr_array));
    auc_array = zeros(length(res_array), length(res_through_array), length(snr_array));
    t2starnr_array = zeros(length(res_array), length(res_through_array), length(snr_array));

%%
for n = 1:length(snr_array)
    tic;
    load(cat(2, base_dir, 'T2starMap_Blocked_LinReg_Transmural', num2str(transmural_array(k)), '_NoiseLevel', num2str(n), '.mat'));

    % if n == 1
    %     figure();
    %     imagesc(abs(t2star_map(:,:,1,1,1))); axis off; axis equal; colormap gray; caxis([0 100]);
    %     roi = drawcircle('StripeColor','y');
    %     bw = createMask(roi);
    % end

% nte = 1;
% snr_array = zeros(length(res_array), length(res_through_array));
% for i = 1:length(res_array)
%     %for i = 8:8
%     for j = 1:length(res_through_array)
%         snr_array(i,j) = mean(nonzeros(bw .* abs(C(:,:,:,nte,i,j)))) / std(nonzeros(bw .* abs(C(:,:,:,nte,i,j))));
%     end
% end

    t2star_map(isnan(t2star_map)) = 0;
    t2star_map(t2star_map <= 0.198) = 0.198;
    t2star_map(t2star_map >= 100) = 100;
    
    
    for i = 1:length(res_array)
        %for i = 8:8
        for j = 1:length(res_through_array)
            t2starnr_array(i,j,n) = mean(nonzeros(remote_roi_binary .* abs(t2star_map(:,:,:,i,j)))) / std(nonzeros(remote_roi_binary .* abs(t2star_map(:,:,:,i,j))));
        end
    end

    %% Detection (Mean-2SD)
    thresh_array = zeros(length(res_array), length(res_through_array));
    hemo_mask = zeros(size(t2star_map));
    
    
    for j = 1:length(res_through_array)
        figure('Position', [100 100 1600 400]);
        for i = 1:length(res_array)
            thresh_array(i,j) = mean(nonzeros(remote_roi_binary .* abs(t2star_map(:,:,:,i,j)))) - 2*std(nonzeros(remote_roi_binary .* abs(t2star_map(:,:,:,i,j))));
            hemo_mask(:,:,:,i,j) = abs(t2star_map(:,:,:,i,j)) < thresh_array(i,j);

            subplot(2,length(res_array),i);
            imagesc(hemo_mask(:,:,1,i,j)); axis off; axis equal;
            subplot(2,length(res_array),i+length(res_array));
            imagesc(remote_roi_binary(:,:,1) .* abs(t2star_map(:,:,1,i,j))); axis off; axis equal;
        end
    end

    %%
    % convert for ROC analysis
    score_map = (abs(t2star_map) - min(vec(abs(t2star_map)))) / (max(vec(abs(t2star_map))) - min(vec(abs(t2star_map))));
    
    figure('Position', [100 100 1600 800]);
    for j = 1:length(res_through_array)
        %
        for i = 1:length(res_array)
            TP = sum(vec(hemo_mask(:,:,:,i,j) + hemo_mask_gt == 4));
            TN = sum(vec(hemo_mask(:,:,:,i,j) + hemo_mask_gt == 2));
            FP = sum(vec(2*hemo_mask(:,:,:,i,j) + hemo_mask_gt == 4));
            FN = sum(vec(2*hemo_mask(:,:,:,i,j) + hemo_mask_gt == 3));
            accuracy = (TP+TN) / (TP+TN+FP+FN);
            
            accuracy_array(i,j,n) = accuracy;
            sensitivity_array(i,j,n) = (TP) / (TP+FN);
            specificity_array(i,j,n) = (TN) / (TN+FP);
            dice_array(i,j,n) = (2*TP) / (2*TP+FP+FN);
    
            score_map_naned = score_map(:,:,:,i,j);
            score_map_naned(air_mask_gt == 1) = [];
    
    
            [X,Y,T,AUC,OPTROCPT] = perfcurve(remote_mask_gt_binary_naned, score_map_naned, 1);
            subplot(length(res_through_array),length(res_array),i+length(res_array)*(j-1));
            plot(X,Y); axis equal;
            auc_array(i,j,n) = AUC;
        end
    end
    toc;
end
%%
detection_analysis.accuracy_array = accuracy_array;
detection_analysis.sensitivity_array = sensitivity_array;
detection_analysis.specificity_array = specificity_array;
detection_analysis.dice_array = dice_array;
detection_analysis.auc_array = auc_array;
detection_analysis.t2starnr_array = t2starnr_array;

save_dir = cat(2, base_dir, 'Analysis\');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

save(cat(2, save_dir, 'Cylinder_Detection_Transmural', num2str(transmural_array(k)), '.mat'), 'detection_analysis', '-v7.3');


%% voxel size vs SNR
voxel_sz_array = sqrt((res_array.^2)' * res_through_array);
[voxel_sz_array_sorted, idx] = sort(voxel_sz_array(:));
t2starnr_array_reshape = reshape(t2starnr_array, [], size(t2starnr_array,3));
figure(); 
for i = 1:size(t2starnr_array,3)
    subplot(3,3,i);
    plot(voxel_sz_array_sorted, t2starnr_array_reshape(idx,i));
end

