clear all;
close all;

addpath('.\function\');
addpath('C:\Users\xz100\Documents\MATLAB\PhD_Project_Code_Storage-master\function\');
base_dir = 'C:\Users\xz100\Documents\Data\T2star_SimulationPhantom\';


vec = @(x) x(:);

load(cat(2, base_dir, '\Ellipsoid_Phantom_FurtherToEndo.mat'));
% remote_balanced = load(cat(2, base_dir, '\Ellipsoid_Phantom_BalanceRemote.mat'));
% remote_balanced = load(cat(2, base_dir, '\Ellipsoid_Phantom_ThirdMyo.mat'));
remote_balanced = load(cat(2, base_dir, '\Ellipsoid_Phantom_HalfMyo.mat'));
load('C:\Users\xz100\Documents\Data\T2star_SimulationPhantom\Analysis\Ellipsoid_RemoteROI.mat');
%Surrounding_Iron = load('C:\Users\xz100\Documents\Data\T2star_SimulationPhantom\Ellipsoid_Phantom_RemoteSurroundingIron+_2X.mat');
%% See the noise level
detection_analysis = struct;
% snr_array = T2star_Parameters.snr_array;
snr_array = [0.0250 0.0500    0.1000    0.1500    0.2000    0.2500    0.3000    0.3500    0.4000    0.4500 0.5000 1.0000, 2.0000];
res_array = [0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4]; % in mm
res_through_array = [2, 4, 6, 8]; % in mm
transmural_array = [0.025, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8];
TE_array = [2.55, 5.80, 9.90, 15.56, 21.22]';
select_indices = [1:13];
% select_indices = [2 4 7 13];
snr_array_selected = snr_array(select_indices);

%% See the acutally images
n = 3;
k = 1;
load(cat(2, 'C:\Users\xz100\Documents\Data\T2star_SimulationPhantom\Ellip_T2starMetrics_Blocked_LinReg_Transmural', num2str(transmural), '_ComprehensiveNoiseLevel_', num2str(n), '.mat'));

figure('Position', [100 100 1600 800]);
for j = 1:length(res_through_array)
    for i = 1:length(res_array)
        subplot(length(res_through_array), length(res_array), i+(j-1)*length(res_array))
        imagesc(t2star_metrics.t2star_map(:,:,1,i,j)); axis equal; axis off; clim([0 50]);
        %colormap(brewermap([],'*RdYlBu'));
        colormap(brewermap([],'RdBu'));
    end
end
%% Should iterate over different transmurality
for k = 1:1
%for k = 1:length(transmural_array)
    transmural = transmural_array(k);
    VObj = Phantom_shape_cell{k};
    sz = size(VObj.t2star);

    Nz = 48;
    Nx = 1024;
    Ny = 1024;

    [X, Y] = meshgrid(1:sz(1));
    [Xq, Yq] = meshgrid(linspace(1, sz(1), Nx));

    t2star_3d_interp = zeros([Nx, Ny, Nz]);

    for slc = 1:Nz
        t2star_3d_interp(:,:,slc) = interp2(X, Y, VObj.t2star(:,:,slc), Xq, Yq, 'spline');
    end

    hemo_mask_gt = t2star_3d_interp;
    hemo_mask_gt(abs(hemo_mask_gt) >= 0.01375) = 0; % Blood Pool
    hemo_mask_gt(abs(hemo_mask_gt) < 0.0025 & abs(hemo_mask_gt) > 0) = 0;
    hemo_mask_gt(abs(hemo_mask_gt) < 0.01375 & abs(hemo_mask_gt) >= 0.0075) = 2; % Remote
    hemo_mask_gt(abs(hemo_mask_gt) < 0.0075 & abs(hemo_mask_gt) >= 0.0025) = 3;    % Hemo

    se1 = strel('disk', 2);
    hemo_mask_gt = imerode(hemo_mask_gt,se1);

    myo_mask_gt_binary = hemo_mask_gt > 0;
    hemo_mask_gt_binary = hemo_mask_gt > 2;
    remote_mask_gt_binary = ~hemo_mask_gt_binary;
    air_mask_gt = hemo_mask_gt == 0;
    % remote_mask_gt_binary_naned = remote_mask_gt_binary(myo_mask_gt);

%% Try to balance the remote region
    %By eroding the edge doesn't work properly for remote ROIs
    %se2 = strel('disk', 30);
    %remote_roi_binary = imerode(remote_mask_gt_binary.* myo_mask_gt, se2);
    
    VObj = remote_balanced.Phantom_shape_cell{k};

    [X, Y] = meshgrid(1:sz(1));
    [Xq, Yq] = meshgrid(linspace(1, sz(1), Nx));

    t2star_3d_interp = zeros([Nx, Ny, Nz]);

    for slc = 1:Nz
        t2star_3d_interp(:,:,slc) = interp2(X, Y, VObj.t2star(:,:,slc), Xq, Yq, 'spline');
    end

    myo_mask_gt_balanced = t2star_3d_interp;
    myo_mask_gt_balanced(abs(myo_mask_gt_balanced) >= 0.01375) = 0; % Blood Pool
    myo_mask_gt_balanced(abs(myo_mask_gt_balanced) < 0.0025 & abs(myo_mask_gt_balanced) > 0) = 0;
    myo_mask_gt_balanced(abs(myo_mask_gt_balanced) < 0.01375 & abs(myo_mask_gt_balanced) >= 0.0075) = 2; % Remote
    myo_mask_gt_balanced(abs(myo_mask_gt_balanced) < 0.0075 & abs(myo_mask_gt_balanced) >= 0.0025) = 3;    % Hemo

    se1 = strel('disk', 2);
    myo_mask_gt_balanced = imerode(myo_mask_gt_balanced,se1);
    myo_mask_gt_balanced_binary = myo_mask_gt_balanced > 1;

    
    myo_mask_gt_balanced_binary = (hemo_mask_gt_binary + myo_mask_gt_balanced_binary) > 0;
    
    %% + hemo surrounding
    % VObj = Surrounding_Iron.Phantom_shape_cell{k};
    % [X, Y] = meshgrid(1:sz(1));
    % [Xq, Yq] = meshgrid(linspace(1, sz(1), Nx));
    % 
    % t2star_3d_interp = zeros([Nx, Ny, Nz]);
    % 
    % for slc = 1:Nz
    %     t2star_3d_interp(:,:,slc) = interp2(X, Y, VObj.t2star(:,:,slc), Xq, Yq, 'spline');
    % end
    % 
    % ironsur_mask_gt_balanced = t2star_3d_interp;
    % ironsur_mask_gt_balanced(abs(ironsur_mask_gt_balanced) >= 0.01375) = 0; % Blood Pool
    % ironsur_mask_gt_balanced(abs(ironsur_mask_gt_balanced) < 0.0025 & abs(ironsur_mask_gt_balanced) > 0) = 0;
    % ironsur_mask_gt_balanced(abs(ironsur_mask_gt_balanced) < 0.01375 & abs(ironsur_mask_gt_balanced) >= 0.0075) = 2; % Remote
    % ironsur_mask_gt_balanced(abs(ironsur_mask_gt_balanced) < 0.0075 & abs(ironsur_mask_gt_balanced) >= 0.0025) = 3;    % Hemo
    % 
    % se1 = strel('disk', 2);
    % ironsur_mask_gt_balanced = imerode(ironsur_mask_gt_balanced,se1);
    % ironsur_mask_gt_balanced_binary = ironsur_mask_gt_balanced > 2;

    %myo_mask_gt_balanced_binary = (ironsur_mask_gt_balanced_binary+myo_mask_gt_balanced_binary) > 0;
    %myo_mask_gt_balanced_binary = (ironsur_mask_gt_balanced_binary) > 0;
    %%
    layers = 1:Nz;
    %clear Phantom_shape_cell
    accuracy_array = zeros(length(res_array), length(res_through_array), length(snr_array_selected), length(layers));
    sensitivity_array = zeros(length(res_array), length(res_through_array), length(snr_array_selected), length(layers));
    dice_array = zeros(length(res_array), length(res_through_array), length(snr_array_selected), length(layers));
    specificity_array = zeros(length(res_array), length(res_through_array), length(snr_array_selected), length(layers));
    auc_array = zeros(length(res_array), length(res_through_array), length(snr_array_selected), length(layers));
    %t2starnr_array = zeros(length(res_array), length(res_through_array), length(snr_array), length(layers));


    %% Manual Mask

    % bw = zeros(size(t2star_map,[1 2 3]));
    % figure();
    % for slc = 1:size(t2star_map,3)
    %     imagesc(abs(t2star_map(:,:,slc,8,4))); axis off; axis equal; colormap gray; caxis([0 100]);
    %     %roi = drawcircle('StripeColor','y');
    %     roi = drawpolygon('StripeColor','y');
    %     bw(:,:,slc) = createMask(roi);
    % end

    %%
    for n = 1:length(snr_array_selected)
    %for n = 2:2
        %tic;
        select_idx = select_indices(n);

        load(cat(2, 'C:\Users\xz100\Documents\Data\T2star_SimulationPhantom\Ellip_T2starMetrics_Blocked_LinReg_Transmural', num2str(transmural), '_ComprehensiveNoiseLevel_', num2str(select_idx), '.mat'));

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

        % t2star_map(isnan(t2star_map)) = 0;
        % t2star_map(t2star_map <= 0.198) = 0.198;
        % t2star_map(t2star_map >= 100) = 100;
        t2star_map = t2star_metrics.t2star_map;

        % for i = 1:length(res_array)
        %     %for i = 8:8
        %     for j = 1:length(res_through_array)
        %         t2starnr_array(i,j,n) = mean(nonzeros(remote_roi_binary .* abs(t2star_map(:,:,:,i,j)))) / std(nonzeros(remote_roi_binary .* abs(t2star_map(:,:,:,i,j))));
        %     end
        % end

        %% Detection (Mean-2SD)
        % thresh_array = zeros(length(res_array), length(res_through_array));
        hemo_mask = zeros(size(t2star_map));


        for j = 1:length(res_through_array)
            %figure('Position', [100 100 1600 400]);
            for i = 1:length(res_array)
                %for m = 1:length(layers)
                % thresh_array(i,j) = mean(nonzeros(remote_roi_binary .* abs(t2star_map(:,:,:,i,j)))) - 2*std(nonzeros(remote_roi_binary .* abs(t2star_map(:,:,:,i,j))));
                thresh = mean(nonzeros(bw .* abs(t2star_map(:,:,:,i,j)))) - 2*std(nonzeros(bw .* abs(t2star_map(:,:,:,i,j))));
                hemo_mask(:,:,:,i,j) = abs(t2star_map(:,:,:,i,j)) < thresh;

                %subplot(2,length(res_array),i);
                %imagesc(hemo_mask(:,:,1,i,j)); axis off; axis equal;
                %subplot(2,length(res_array),i+length(res_array));
                %imagesc(remote_roi_binary(:,:,1) .* abs(t2star_map(:,:,1,i,j))); axis off; axis equal;
                %end
            end
        end

        %%
        % convert for ROC analysis
        layers = [1:Nz];
        %layers = [1 24 48];
        score_map = 1 - (abs(t2star_map) - min(vec(abs(t2star_map)))) / (max(vec(abs(t2star_map))) - min(vec(abs(t2star_map))));
        % T_trimmed = [0.673482106729676	0.673482106729676	0.673156034608656	0.673039135329634	0.672857429730983	0.671792127851913	0.657177012700051	0.656044628855927	0.654376778567421	0.653834385409557	0.635506489356952	0.635016784332838	0.633917834207751	0.632767949322981	0.616249472828360	0.613413539342221	0.605624430411888	0.605511662994926	0.603968271105997];
        hemo_mask_gt_balanced = 2*myo_mask_gt_balanced_binary+hemo_mask_gt;
        %figure('Position', [100 100 1600 800]);
        for j = 1:length(res_through_array)
            %
            for i = 1:length(res_array)

                for m = 1:length(layers)
                %for m = 1:1
                    layer = layers(m);
                    % TP = sum(vec(hemo_mask(:,:,layer,i,j) + hemo_mask_gt(:,:,layer) == 4));
                    % TN = sum(vec(hemo_mask(:,:,layer,i,j) + hemo_mask_gt(:,:,layer) == 2));
                    % FP = sum(vec(2*hemo_mask(:,:,layer,i,j) + hemo_mask_gt(:,:,layer) == 4));
                    % FN = sum(vec(2*hemo_mask(:,:,layer,i,j) + hemo_mask_gt(:,:,layer) == 3));

                    TP = sum(vec(hemo_mask(:,:,layer,i,j) + hemo_mask_gt_balanced(:,:,layer) == 6));
                    TN = sum(vec(hemo_mask(:,:,layer,i,j) + hemo_mask_gt_balanced(:,:,layer) == 4));
                    FP = sum(vec(2*hemo_mask(:,:,layer,i,j) + hemo_mask_gt_balanced(:,:,layer) == 6));
                    FN = sum(vec(2*hemo_mask(:,:,layer,i,j) + hemo_mask_gt_balanced(:,:,layer) == 5));

                    % TODO it's not balanced yet
                    accuracy = (TP+TN) / (TP+TN+FP+FN);

                    accuracy_array(i,j,n,layer) = accuracy;
                    sensitivity_array(i,j,n,layer) = (TP) / (TP+FN);
                    specificity_array(i,j,n,layer) = (TN) / (TN+FP);
                    dice_array(i,j,n,layer) = (2*TP) / (2*TP+FP+FN);

                    %score_map_naned = score_map(:,:,layer,i,j);
                    % score_map_naned(air_mask_gt(:,:,layer) == 1) = [];
                    %score_map_naned(myo_mask_gt_balanced_binary(:,:,layer) == 0) = [];

                    % score_map_naned = score_map(:,:,layer,i,j);
                    % score_map_naned((hemo_mask_gt_binary(:,:,layer) + myo_mask_gt_balanced_binary(:,:,layer)) == 0) = [];
                    
                    score_map_naned = score_map(:,:,layer,i,j);
                    %score_map_naned(myo_mask_gt_binary(:,:,layer) == 0) = [];
                    score_map_naned(myo_mask_gt_balanced_binary(:,:,layer) == 0) = [];

                    %remote_mask_gt_binary_layered = remote_mask_gt_binary(:,:,layer);
                    %remote_mask_gt_binary_naned = remote_mask_gt_binary_layered(myo_mask_gt_balanced_binary(:,:,layer)); % It means everything but hemo

                    %hemo_mask_gt_binary_naned = hemo_mask_gt_binary(:,:,layer);
                    %hemo_mask_gt_binary_naned(myo_mask_gt_balanced_binary(:,:,layer) == 0) = [];

                    % hemo_mask_gt_binary_naned = hemo_mask_gt_binary(:,:,layer);
                    % hemo_mask_gt_binary_naned((hemo_mask_gt_binary(:,:,layer) + myo_mask_gt_balanced_binary(:,:,layer)) == 0) = [];

                    hemo_mask_gt_binary_naned = hemo_mask_gt_binary(:,:,layer);
                    %hemo_mask_gt_binary_naned((myo_mask_gt_binary(:,:,layer)) == 0) = [];
                    hemo_mask_gt_binary_naned((myo_mask_gt_balanced_binary(:,:,layer)) == 0) = [];


                    [X,Y,T,AUC,OPTROCPT] = perfcurve(hemo_mask_gt_binary_naned, score_map_naned, 1);
                    %[X,Y,T,AUC,OPTROCPT] = perfcurve(hemo_mask_gt_binary_naned, score_map_naned, 1, 'TVals', T_trimmed);
                    %subplot(length(res_through_array),length(res_array),i+length(res_array)*(j-1));
                    %plot(X,Y); xlim([0 1]); ylim([0 1]); axis equal;
                    auc_array(i,j,n,layer) = AUC;
                end
            end
        end
        %toc;
    end
    %%
    detection_analysis.accuracy_array = accuracy_array;
    detection_analysis.sensitivity_array = sensitivity_array;
    detection_analysis.specificity_array = specificity_array;
    detection_analysis.dice_array = dice_array;
    detection_analysis.auc_array = auc_array;

    % detection_analysis.t2starnr_array = t2starnr_array;

    save_dir = cat(2, base_dir, 'Analysis\');
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    save(cat(2, save_dir, 'Ellipsoid_Detection_HalfMyo_Transmural', num2str(transmural_array(k)), '.mat'), 'detection_analysis', '-v7.3');
    %save(cat(2, save_dir, 'Ellipsoid_RemoteROI.mat'), 'bw', '-v7.3');
end

%% Check what's going on with detections
hemo_mask = zeros(size(t2star_map));
figure('Position', [100 100 1600 400]);

slc = 1;
for j = 1:length(res_through_array)
    for i = 1:length(res_array)
        %for m = 1:length(layers)
        % thresh_array(i,j) = mean(nonzeros(remote_roi_binary .* abs(t2star_map(:,:,:,i,j)))) - 2*std(nonzeros(remote_roi_binary .* abs(t2star_map(:,:,:,i,j))));
        thresh = mean(nonzeros(bw .* abs(t2star_map(:,:,:,i,j)))) - 2*std(nonzeros(bw .* abs(t2star_map(:,:,:,i,j))));
        hemo_mask(:,:,slc,i,j) = abs(t2star_map(:,:,slc,i,j)) < thresh;

        subplot(length(res_through_array),length(res_array),i+length(res_array)*(j-1));
        imagesc(hemo_mask(:,:,slc,i,j)+hemo_mask_gt_binary(:,:,slc)); axis off; axis equal;
        %imagesc(t2star_map(:,:,slc,i,j)); axis off; axis equal; clim([0 100]);
    end
end

%%
figure('Position', [100 100 400 400]);
slc = 1;
imagesc(myo_mask_gt_balanced_binary(:,:,slc) + hemo_mask_gt_binary(:,:,slc)); axis off; axis equal;

%% voxel size vs SNR
voxel_sz_array = sqrt((res_array.^2)' * res_through_array);
[voxel_sz_array_sorted, idx] = sort(voxel_sz_array(:));
t2starnr_array_reshape = reshape(squeeze(t2starnr_array(:,:,1,:)), [], size(t2starnr_array,4));
figure(); 
for i = 1:size(t2starnr_array,4)
    subplot(3,5,i);
    plot(voxel_sz_array_sorted, t2starnr_array_reshape(idx,i));
end

%% Voxel size vs T2wSNR
%% voxel size vs T2starwSNR
load(cat(2, 'C:\Users\xz100\Documents\Data\T2star_SimulationPhantom\Ellip_T2starWSNR_Metrics.mat'), 't2star_metrics', '-v7.3');
%%
% voxel_sz_array = sqrt((res_array.^2)' * res_through_array);
% [voxel_sz_array_sorted, idx] = sort(voxel_sz_array(:));
% 
% 
% for nte = 1:length(MxyTE_remote_echo)
%     t2starwnr_array_reshape = reshape(t2starwnr_array(:,:,nte,:), [], size(t2starwnr_array,4));
%     figure();
%     for i = 1:size(t2starwnr_array,4)
%         subplot(3,3,i);
%         plot(voxel_sz_array_sorted, t2starwnr_array_reshape(idx,i));
%     end
% end

%% Fit
% nte = 1;
% t2starwnr_array_reshape = reshape(t2starwnr_array(:,:,nte,:), [], size(t2starwnr_array,4));
% mdl_cell = cell(size(t2starwnr_array,4),1);
% slope_array = zeros(size(t2starwnr_array,4),1);
% figure();
% for i = 1:size(t2starwnr_array,4)
%     plot(voxel_sz_array_sorted, t2starwnr_array_reshape(idx,i));
%     hold on;
% 
%     mdl_cell{i} = fitlm(voxel_sz_array_sorted,t2starwnr_array_reshape(idx,i));
%     slope_array(i) = mdl_cell{i}.Coefficients.Estimate(2);
% end

%%  AUC
slc_start = 33;
slc_end = 48;
voxel_sz_array = (res_array.^2)' * res_through_array;
[voxel_sz_array_sorted, idx] = sort(voxel_sz_array(:));
auc_array_reshape = reshape(squeeze(auc_array(:,:,:,slc)), [], size(auc_array,3));
auc_array_reshape = reshape(mean(auc_array(:,:,:,slc_start:slc_end), 4), [], size(auc_array,3));

figure();
for i = 1:size(auc_array,3)
    subplot(3,3,i)
    plot(voxel_sz_array_sorted, auc_array_reshape(idx,i)); ylim([0.5 1])
end

%% AUC V2
inplane_res  = [0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4];
thrplane_res = [2, 4, 6, 8];   
vol_mat = ((inplane_res .* inplane_res)' * thrplane_res)';

n=2
slc_start = 33;
slc_end = 48;
AUC_exvivo_allavg = squeeze(mean(auc_array(:,:,n,slc_start:slc_end), 4));
figure('Position', [100 0 400 400]);

%Y = AUC_avg16_array(I);
%X = B;
%tbl = table(X, Y);
%modelfun = @(b,x) b(1) + b(2)*exp(b(3)*x);
%modelfun = @(b,x) b(1) + b(2)*x.^b(3);
%beta0 = [0 0 0];
%mdl = fitnlm(tbl,modelfun,beta0);

%ci = coefCI(mdl);
%b = mdl.Coefficients.Estimate;

%Y_pred = modelfun(b, X)
hold on;
%plot(X, Y_pred); %ylim([0.5 1])
%Y_lb = modelfun(ci(:,1), X);
%Y_ub = modelfun(ci(:,2), X);
%plot(X, Y_lb);
%plot(X, Y_ub);

% set(plotHandles_auc(:,1), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 12, ...
%     'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor' , [.75 .75 1]);
%set(plotHandles_auc(:,1), 'Visible','off');
plotHandles_auc(:,2) = plot(vol_mat(1,:), AUC_exvivo_allavg(:,1),'o');
set(plotHandles_auc(:,2), 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
plotHandles_auc(:,3) = plot(vol_mat(2,:), AUC_exvivo_allavg(:,2),'square');
set(plotHandles_auc(:,3), 'LineWidth', 1, 'Marker', 'square', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
plotHandles_auc(:,4) = plot(vol_mat(3,:), AUC_exvivo_allavg(:,3),'diamond');
set(plotHandles_auc(:,4), 'LineWidth', 1, 'Marker', 'diamond', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
plotHandles_auc(:,5) = plot(vol_mat(4,:), AUC_exvivo_allavg(:,4),'^');
set(plotHandles_auc(:,5), 'LineWidth', 1, 'Marker', '^', 'MarkerSize', 12, ...
    'MarkerEdgeColor', [0 0 0]/255, 'MarkerFaceColor' , [253,190,133]/255);
xlim([-1, 36]); 
ylim([0.5 1.0]);
%set(gca, 'XTick', [0.04, 0.10, 0.15, 0.20]);
set(gca, 'XTickLabels', []);
set(gca, 'XTick',[0 10 20 30 36]);
set(gca, 'YTick',[0.5 0.6 0.7 0.8 0.9 1]);
set(gca, 'YTickLabels', []);
set(gca,'LineWidth', 1.5,'TickLength',[0.02 0.02]);
set(gca,'TickDir','out', 'YGrid', 'on');
set(gca,'box','off');
%% Sensitivity
slc_start = 1;
slc_end = 48;
res_through_array = [2 4 6 8];
voxel_sz_array = (res_array.^2)' * res_through_array;
[voxel_sz_array_sorted, idx] = sort(voxel_sz_array(:));
sensitivity_array_reshape = reshape(squeeze(mean(sensitivity_array(:,:,:,slc_start:slc_end), 4)), [], size(sensitivity_array,3));


figure();
for i = 1:size(sensitivity_array,3)
    subplot(3,3,i)
    plot(voxel_sz_array_sorted, sensitivity_array_reshape(idx,i)); ylim([0 1])
end

%% Accuracy
slc = 1;
res_through_array = [2 4 6 8];
voxel_sz_array = (res_array.^2)' * res_through_array;
[voxel_sz_array_sorted, idx] = sort(voxel_sz_array(:));
accuracy_array_reshape = reshape(squeeze(accuracy_array(:,:,:,slc)), [], size(accuracy_array,3));
figure();
for i = 1:size(accuracy_array,3)
    subplot(3,3,i)
    plot(voxel_sz_array_sorted, accuracy_array_reshape(idx,i)); ylim([0 1])
end

%% Dice
slc = 1;
res_through_array = [2 4 6 8];
voxel_sz_array = (res_array.^2)' * res_through_array;
[voxel_sz_array_sorted, idx] = sort(voxel_sz_array(:));
dice_array_reshape = reshape(squeeze(dice_array(:,:,:,slc)), [], size(dice_array,3));
figure();
for i = 1:size(accuracy_array,3)
    subplot(3,3,i)
    plot(voxel_sz_array_sorted, dice_array_reshape(idx,i)); ylim([0 1])
end

%% Specificity
slc = 1;
res_through_array = [2 4 6 8];
voxel_sz_array = (res_array.^2)' * res_through_array;
[voxel_sz_array_sorted, idx] = sort(voxel_sz_array(:));
specificity_array_reshape = reshape(squeeze(specificity_array(:,:,:,slc)), [], size(specificity_array,3));
figure();
for i = 1:size(accuracy_array,3)
    subplot(3,3,i)
    plot(voxel_sz_array_sorted, specificity_array_reshape(idx,i)); ylim([0 1])
end

%% For display
transmural = transmural_array(k);
VObj = Phantom_shape_cell{k};
sz = size(VObj.t2star);

Nz = 72;
Nx = 1024;
Ny = 1024;

[X, Y] = meshgrid(1:sz(1));
[Xq, Yq] = meshgrid(linspace(1, sz(1), Nx));

t2star_3d_interp = zeros([Nx, Ny, Nz]);

for slc = 1:Nz
    t2star_3d_interp(:,:,slc) = interp2(X, Y, VObj.t2star(:,:,slc), Xq, Yq, 'spline');
end
%% SNR analysis
% k = 1;
% t2starwnr_array = zeros(length(res_array), length(res_through_array), length(TE_array), length(snr_array));
% t2starnr_array = zeros(length(res_array), length(res_through_array), length(snr_array));
% for n = 1:length(snr_array)
%     load(cat(2, base_dir, 'Ellip_T2starMetrics_Blocked_LinReg_Transmural', num2str(transmural_array(k)), '_NoiseLevel', num2str(n), '.mat'));
%     t2starwnr_array(:,:,:,n) = t2star_metrics.t2starwnr_array(:,:,:,n);
%     t2starnr_array(:,:,n) = t2star_metrics.t2starnr_array(:,:,n);
% end
% 
% % voxel size vs SNR
% voxel_sz_array = sqrt((res_array.^2)' * res_through_array);
% [voxel_sz_array_sorted, idx] = sort(voxel_sz_array(:));
% t2starnr_array_reshape = reshape(t2starnr_array, [], size(t2starnr_array,3));
% figure(); 
% for i = 1:size(t2starnr_array,3)
%     subplot(3,3,i);
%     plot(voxel_sz_array_sorted, t2starnr_array_reshape(idx,i));
% end
% 
% % voxel size vs T2wSNR
% voxel_sz_array = sqrt((res_array.^2)' * res_through_array);
% [voxel_sz_array_sorted, idx] = sort(voxel_sz_array(:));
% t2starwnr_array_squeeze = squeeze(t2starwnr_array(:,:,1,:));
% t2starwnr_array_reshape = reshape(t2starwnr_array_squeeze, [], size(t2starwnr_array,4));
% figure(); 
% for i = 1:size(t2starwnr_array,4)
%     subplot(4,8,i);
%     plot(voxel_sz_array_sorted, t2starwnr_array_reshape(idx,i));
% end


