clear all;
close all;

%% Load SimPhantom_04042021.mat
addpath('../function/');

base_dir = uigetdir;
f_to_read = cat(2, base_dir, '/SimPhantom_03162023.mat');
load(f_to_read);

res_array = SimPhantom_03162023.res_array;
sigma_array = SimPhantom_03162023.sigma_array;
width_max = SimPhantom_03162023.width_max;

C_t2star_fit_reshape = SimPhantom_03162023.tissue_canvas{1}(1).C_cell{6,1};
% Width x Height x Res x Sigma
figure();
imagesc(C_t2star_fit_reshape(:,:,1))
figure();
imagesc(C_t2star_fit_reshape(:,:,5))
figure();
imagesc(C_t2star_fit_reshape(:,:,10))
figure();
imagesc(C_t2star_fit_reshape(:,:,15))

%%
% Renew the hemo_mask and and myo_mask
dx = 0.1;
CNR_array = zeros(length(sigma_array), length(res_array), length(SimPhantom_03162023.tissue_canvas{1}));
SNR_array = zeros(length(sigma_array), length(res_array), length(SimPhantom_03162023.tissue_canvas{1}));
C_cell = SimPhantom_03162023.tissue_canvas{1}(1).C_cell;
Nx = size(C_cell{1,1}, 2);
Ny = size(C_cell{1,1}, 1);
myo_mask = zeros(Nx, Ny);

hemo_w = res_array(end)/dx;
hemo_ww = fix(hemo_w/2);
myo_mask(:, 1:(Nx/2-3*hemo_ww)) = ones(Ny, length(1:(Nx/2-3*hemo_ww)));
myo_mask(:, (Nx/2+3*hemo_ww)+1:end) = ones(Ny, length(1:(Nx/2-3*hemo_ww)));
t_hemo = zeros(Ny, Nx, length(width_max));
hemo_mask_cell = cell(length(res_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
hemo_mask_frame_cell = cell(length(res_array), length(width_max));

for iter = 1:length(SimPhantom_03162023.tissue_canvas)
    for w = 1:length(width_max)
        % res = res_array(i);
        t_width = width_max(w);
        C_cell = SimPhantom_03162023.tissue_canvas{iter}(w).C_cell;

        % Nx_hemo = res / dx;
        hemo_w = round(t_width/dx);
        hemo_ww = fix(hemo_w/2);

        if mod(hemo_w, 2) == 0
            t_hemo(:, (Nx/2-hemo_ww):(Nx/2+hemo_ww-1),w) = 1;
        else
            t_hemo(:, (Nx/2-hemo_ww):(Nx/2+hemo_ww),w) = 1;
        end

        for i = 1:size(C_cell, 1)
            res = res_array(i);
            blocs = Func_map_to_bloc_encoded_exaustive(dx, Nx, res, C_cell{i,1});
            slides = size(C_cell{i,1}, 3);
            hemo_mask_frame = zeros(Ny, Nx, slides);
            for sld = 1:slides
                blocs_overlap = t_hemo(:,:,w) .* blocs(:,:,sld);
                idx_blocs = nonzeros(unique(blocs_overlap));
                for k = 1:length(idx_blocs)
                    [row, col] = ind2sub([Ny Nx], find(blocs(:,:,sld) == idx_blocs(k)));
                    hemo_mask_frame(row, col, sld) = 1;
                end
            end
            hemo_mask_frame_cell{i,w} = hemo_mask_frame;
            hemo_mask = zeros(Ny, Nx, slides, length(sigma_array));
            for j = 1:size(C_cell, 2)
                for sld = 1:slides
                    thresh = mean(nonzeros(myo_mask .* C_cell{i,j}(:,:,sld))) - 2*std(nonzeros(myo_mask .* C_cell{i,j}(:,:,sld)));
                    hemo_mask(:,:,sld,j) = (C_cell{i,j}(:,:,sld) < thresh); %.* hemo_mask_frame(:,:,i,w);
                    %thresh_mat(:,:,i,j,w) = thresh;
                end
            end
            hemo_mask_cell{i,w,iter} = hemo_mask;
        end
    end
end

%%
hemo_mask_frame_w06 = hemo_mask_frame_cell{3,1};
hemo_mask_frame_w20 = hemo_mask_frame_cell{3,2};
figure(); 
for w = 1:size(hemo_mask_frame_w20, 3)
    subplot(8,2,2*(w-1)+1);
    imagesc(hemo_mask_frame_w06(:,:,w));
    axis off;
    subplot(8,2,2*(w-1)+2);
    imagesc(hemo_mask_frame_w20(:,:,w));
    axis off;
end

%% Sensitivity, Specificity, Accuracy, Precision
sens_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
spec_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
accu_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
prec_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
CNR_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
SNR_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
CNR_nom_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
auc_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
dice_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));

% VS truth in canvas
sens_cell2 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
spec_cell2 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
accu_cell2 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
prec_cell2 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
auc_cell2 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
dice_cell2 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));

accu_cell3 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
auc_cell3 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));

sens_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
spec_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
accu_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
prec_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
auc_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
dice_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
CNR_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
SNR_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
CNR_nom_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));

% VS truth in canvas
sens_matt2 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
spec_matt2 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
accu_matt2 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
prec_matt2 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
auc_matt2 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
dice_matt2 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));

accu_matt3 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));
auc_matt3 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03162023.tissue_canvas));

myo = myo_mask;
for iter = 1:length(SimPhantom_03162023.tissue_canvas)
    for w = 1:length(width_max)
        t_width = width_max(w);
        C_cell = SimPhantom_03162023.tissue_canvas{iter}(w).C_cell;
        truth_map2 = t_hemo(:,:,w);
        
        for i = 1:length(res_array)
            res = res_array(i);
            bloc_size = round(res/dx).^2;
            truth_map = hemo_mask_frame_cell{i,w}; % 3D truth map
            hemo_mask = hemo_mask_cell{i,w,iter};
            %  Ny x Nx x slides x sigma
            slides = size(hemo_mask, 3);
            for j = 1:length(sigma_array)
                sens_mat = zeros(slides, 1);
                spec_mat = zeros(slides, 1);
                accu_mat = zeros(slides, 1);
                prec_mat = zeros(slides, 1);
                CNR_mat = zeros(slides, 1);
                SNR_mat = zeros(slides, 1);
                CNR_nom_mat = zeros(slides, 1);

                % VS truth in canvas
                sens_mat2 = zeros(slides, 1);
                spec_mat2 = zeros(slides, 1);
                accu_mat2 = zeros(slides, 1);
                prec_mat2 = zeros(slides, 1);
                auc_mat2 = zeros(slides, 1);
                
                accu_mat3 = zeros(slides, 1);
                auc_mat3 = zeros(slides, 1);

                for sld = 1:slides
                    FP = find(truth_map(:,:,sld) - hemo_mask(:,:,sld,j) == -1);
                    FN = find(truth_map(:,:,sld) - hemo_mask(:,:,sld,j) == 1);
                    TP = find(2*truth_map(:,:,sld) - hemo_mask(:,:,sld,j) == 1);
                    TN = find(2*(~truth_map(:,:,sld)) - ~hemo_mask(:,:,sld,j) == 1);
                    
                    hemo_mask_temp = hemo_mask(:,:,sld,j);
                    truth_map_temp = truth_map(:,:,sld);
                    [X,Y,T,AUC,OPTROCPT] = perfcurve(truth_map_temp(:), hemo_mask_temp(:) , 1);

                    fp = length(FP)/bloc_size;
                    fn = length(FN)/bloc_size;
                    tp = length(TP)/bloc_size;
                    tn = length(TN)/bloc_size;

                    sensitivity = tp ./ (tp + fn);
                    specificity = tn ./ (tn + fp);
                    accuracy = (tp + tn) ./ (fp + fn + tp + tn);
                    precision = (tp) ./ (fp + tp);
                    dice = (2*tp) ./ (2*tp + fp + fn);

                    sens_mat(sld) = sensitivity;
                    spec_mat(sld) = specificity;
                    accu_mat(sld) = accuracy;
                    prec_mat(sld) = precision;
                    auc_mat(sld) = AUC;
                    dice_mat(sld) = dice;

                    FP = find(truth_map2 - hemo_mask(:,:,sld,j) == -1);
                    FN = find(truth_map2 - hemo_mask(:,:,sld,j) == 1);
                    TP = find(2*truth_map2 - hemo_mask(:,:,sld,j) == 1);
                    TN = find(2*(~truth_map2) - ~hemo_mask(:,:,sld,j) == 1);
                    
                    [X,Y,T,AUC,OPTROCPT] = perfcurve(truth_map2(:),hemo_mask_temp(:) ,1);
                    
                    fp = length(FP)/bloc_size;
                    fn = length(FN)/bloc_size;
                    tp = length(TP)/bloc_size;
                    tn = length(TN)/bloc_size;

                    sensitivity = tp ./ (tp + fn);
                    specificity = tn ./ (tn + fp);
                    accuracy = (tp + tn) ./ (fp + fn + tp + tn);
                    precision = (tp) ./ (fp + tp);
                    dice = (2*tp) ./ (2*tp + fp + fn);

                    sens_mat2(sld) = sensitivity;
                    spec_mat2(sld) = specificity;
                    accu_mat2(sld) = accuracy;
                    prec_mat2(sld) = precision;
                    auc_mat2(sld) = AUC;
                    dice_mat2(sld) = dice;
                    
                    % [row, col] = find(truth_map2 == 1);
                    FP = find((truth_map2 - hemo_mask(:,:,sld,j)).*truth_map(:,:,sld) == -1);
                    FN = find((truth_map2 - hemo_mask(:,:,sld,j)).*truth_map(:,:,sld) == 1);
                    TP = find((2*truth_map2 - hemo_mask(:,:,sld,j)).*truth_map(:,:,sld) == 1);
                    TN = find((2*(~truth_map2) - ~hemo_mask(:,:,sld,j)).*truth_map(:,:,sld) == 1);
                    
                    fp = length(FP)/bloc_size;
                    fn = length(FN)/bloc_size;
                    tp = length(TP)/bloc_size;
                    tn = length(TN)/bloc_size;
                    
                    hemo_mask_temp = hemo_mask(:,:,sld,j).*truth_map(:,:,sld);
                    
                    [X,Y,T,AUC,OPTROCPT] = perfcurve(truth_map2(:),hemo_mask_temp(:) ,1);
                    accuracy = (tp + tn) ./ (fp + fn + tp + tn);
                    accu_mat3(sld) = accuracy;
                    auc_mat3(sld) = AUC;

                    s_hemo = mean(nonzeros(hemo_mask(:,:,sld,j) .* truth_map(:,:,sld) .* C_cell{i,j}(:,:,sld)));
                    s_hemo_nom = mean(nonzeros(truth_map(:,:,sld) .* C_cell{i,j}(:,:,sld)));
                    s_myo = mean(nonzeros(myo .* C_cell{i,j}(:,:,sld)));
                    std_myo = std(nonzeros(myo .* C_cell{i,j}(:,:,sld)));
                    SNR_mat(sld) = s_myo ./ std_myo;
                    CNR_mat(sld) = (s_myo-s_hemo) ./ std_myo;
                    % A different way to do CNR
                    CNR_nom_mat(sld) = (s_myo-s_hemo_nom) ./ std_myo;
                end

                sens_cell{i,j,w,iter} = sens_mat;
                spec_cell{i,j,w,iter} = spec_mat;
                accu_cell{i,j,w,iter} = accu_mat;
                prec_cell{i,j,w,iter} = prec_mat;
                auc_cell{i,j,w,iter} = auc_mat;
                dice_cell{i,j,w,iter} = dice_mat;
                CNR_cell{i,j,w,iter} = CNR_mat;
                SNR_cell{i,j,w,iter} = SNR_mat;
                CNR_nom_cell{i,j,w,iter} = CNR_nom_mat;

                % VS truth in canvas
                sens_cell2{i,j,w,iter} = sens_mat2;
                spec_cell2{i,j,w,iter} = spec_mat2;
                accu_cell2{i,j,w,iter} = accu_mat2;
                prec_cell2{i,j,w,iter} = prec_mat2;
                auc_cell2{i,j,w,iter} = auc_mat2;
                dice_cell2{i,j,w,iter} = dice_mat2;
                
                accu_cell3{i,j,w,iter} = accu_mat3;
                auc_cell3{i,j,w,iter} = auc_mat3;

                sens_matt(i,j,w,iter) = mean(sens_mat);
                spec_matt(i,j,w,iter) = mean(spec_mat);
                accu_matt(i,j,w,iter) = mean(accu_mat);
                prec_matt(i,j,w,iter) = mean(prec_mat);
                auc_matt(i,j,w,iter) = mean(auc_mat);
                dice_matt(i,j,w,iter) = mean(dice_mat);

                CNR_matt(i,j,w,iter) = mean(CNR_mat);
                SNR_matt(i,j,w,iter) = mean(SNR_mat);
                CNR_nom_matt(i,j,w,iter) = mean(CNR_nom_mat);

                % VS truth in canvas
                sens_matt2(i,j,w,iter) = mean(sens_mat2);
                spec_matt2(i,j,w,iter) = mean(spec_mat2);
                accu_matt2(i,j,w,iter) = mean(accu_mat2);
                prec_matt2(i,j,w,iter) = mean(prec_mat2);
                auc_matt2(i,j,w,iter) = mean(auc_mat2);
                dice_matt2(i,j,w,iter) = mean(dice_mat2);
                
                accu_matt3(i,j,w,iter) = mean(accu_mat3);
                auc_matt3(i,j,w,iter) = mean(auc_mat3);
            end
        end
    end
end
%%
SNR_mat_06 = mean(SNR_matt(:,:,1,:),4);
SNR_mat_20 = mean(SNR_matt(:,:,2,:),4);
accu_matt2_06 = mean(accu_matt2(:,:,1,:),4);
accu_matt2_20 = mean(accu_matt2(:,:,2,:),4);
sens_matt2_06 = mean(sens_matt2(:,:,1,:),4);
sens_matt2_20 = mean(sens_matt2(:,:,2,:),4);
auc_matt2_06 = mean(auc_matt2(:,:,1,:),4);
auc_matt2_20 = mean(auc_matt2(:,:,2,:),4);
dice_matt2_06 = mean(dice_matt2(:,:,1,:),4);
dice_matt2_20 = mean(dice_matt2(:,:,2,:),4);
auc_matt3_06 = mean(auc_matt3(:,:,1,:),4);
auc_matt3_20 = mean(auc_matt3(:,:,2,:),4);
accu_matt3_06 = mean(accu_matt3(:,:,1,:),4);
accu_matt3_20 = mean(accu_matt3(:,:,2,:),4);
%%

figure();
imagesc(mean(SNR_matt(:,:,1,:),4)); axis image; axis off;
caxis([5 30]);
colormap(brewermap([],'*RdYlBu'));
%colormap(brewermap([],'*YlGnBu'));
colorbar;

[B,I] = sort(SNR_mat_06,1);
I_idx = zeros(size(I));
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        I_idx(I(i,j),j) = i;
    end
end
figure();
imagesc(I_idx);

[B,I] = sort(SNR_mat_20,1);
I_idx = zeros(size(I));
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        I_idx(I(i,j),j) = i;
    end
end
figure();
imagesc(I_idx);
%%
% This is the contrast in nominal hemo
CNR_mat_06 = mean(CNR_nom_matt(:,:,1,:),4);
CNR_mat_20 = mean(CNR_nom_matt(:,:,2,:),4);
figure();
imagesc(mean(CNR_nom_matt(:,:,2,:),4)); axis image; axis off;
caxis([3 20]);
colormap(brewermap([],'*RdYlBu'));
%colormap(brewermap([],'*YlGnBu'));
colorbar;

[B,I] = sort(CNR_mat_06,1);
I_idx = zeros(size(I));
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        I_idx(I(i,j),j) = i;
    end
end
figure();
imagesc(I_idx);

[B,I] = sort(CNR_mat_20,1);
I_idx = zeros(size(I));
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        I_idx(I(i,j),j) = i;
    end
end
figure();
imagesc(I_idx);

%%
figure();
imagesc(mean(sens_matt2(:,:,1,:),4)); axis image; axis off;
caxis([0 1]);
colormap(brewermap([],'*RdYlBu'));
%colormap(brewermap([],'*YlGnBu'));
colorbar;
%%
figure();
imagesc(mean(accu_matt2(:,:,2,:),4)); axis image; axis off;
caxis([0.8 1]);
colormap(brewermap([],'*RdYlBu'));
%colormap(brewermap([],'*YlGnBu'));
colorbar;

accu_matt2_06 = mean(accu_matt2(:,:,1,:),4);
accu_matt2_20 = mean(accu_matt2(:,:,2,:),4);

[B,I] = sort(accu_matt2_06,1);
I_idx = zeros(size(I));
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        I_idx(I(i,j),j) = i;
    end
end
figure();
imagesc(I_idx);

[B,I] = sort(accu_matt2_20,1);
I_idx = zeros(size(I));
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        I_idx(I(i,j),j) = i;
    end
end

figure();
imagesc(I_idx);
%% 
figure();
imagesc(mean(spec_matt2(:,:,2,:),4)); axis image; axis off;
caxis([0.8 1]);
colormap(brewermap([],'*RdYlBu'));
%colormap(brewermap([],'*YlGnBu'));
colorbar;

%% Load 03212023 SimPhantom
base_dir = uigetdir;
f_to_read = cat(2, base_dir, '/SimPhantom_03212023.mat');
load(f_to_read);

res_array = SimPhantom_03212023.res_array;
sigma_array = SimPhantom_03212023.sigma_array;
width_max = SimPhantom_03212023.width_max;

% 2,3,5,6 is not correct
C_t2star_fit_reshape = SimPhantom_03212023.tissue_canvas{1}(2).C_cell{1,1};
% Width x Height x Res x Sigma
figure();
imagesc(C_t2star_fit_reshape(:,:,1))
figure();
imagesc(C_t2star_fit_reshape(:,:,2))
figure();
imagesc(C_t2star_fit_reshape(:,:,3))
figure();
imagesc(C_t2star_fit_reshape(:,:,4))

%%
% Renew the hemo_mask and and myo_mask
dx = 0.1;
CNR_array = zeros(length(sigma_array), length(res_array), length(SimPhantom_03212023.tissue_canvas{1}));
SNR_array = zeros(length(sigma_array), length(res_array), length(SimPhantom_03212023.tissue_canvas{1}));
C_cell = SimPhantom_03212023.tissue_canvas{1}(1).C_cell;
Nx = size(C_cell{1,1}, 2);
Ny = size(C_cell{1,1}, 1);
myo_mask = zeros(Nx, Ny);

hemo_w = res_array(end)/dx;
hemo_ww = fix(hemo_w/2);
myo_mask(:, 1:(Nx/2-3*hemo_ww)) = ones(Ny, length(1:(Nx/2-3*hemo_ww)));
myo_mask(:, (Nx/2+3*hemo_ww)+1:end) = ones(Ny, length(1:(Nx/2-3*hemo_ww)));
t_hemo = zeros(Ny, Nx, length(width_max));
hemo_mask_cell = cell(length(res_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
hemo_mask_frame_cell = cell(length(res_array), length(width_max));

for iter = 1:length(SimPhantom_03212023.tissue_canvas)
    for w = 1:length(width_max)
        % res = res_array(i);
        t_width = width_max(w);
        C_cell = SimPhantom_03212023.tissue_canvas{iter}(w).C_cell;

        % Nx_hemo = res / dx;
        hemo_w = round(t_width/dx);
        hemo_ww = fix(hemo_w/2);

        if mod(hemo_w, 2) == 0
            t_hemo(:, (Nx/2-hemo_ww):(Nx/2+hemo_ww-1),w) = 1;
        else
            t_hemo(:, (Nx/2-hemo_ww):(Nx/2+hemo_ww),w) = 1;
        end

        for i = 1:size(C_cell, 1)
            res = res_array(i);
            blocs = Func_map_to_bloc_encoded_exaustive(dx, Nx, res, C_cell{i,1});
            slides = size(C_cell{i,1}, 3);
            hemo_mask_frame = zeros(Ny, Nx, slides);
            for sld = 1:slides
                blocs_overlap = t_hemo(:,:,w) .* blocs(:,:,sld);
                idx_blocs = nonzeros(unique(blocs_overlap));
                for k = 1:length(idx_blocs)
                    [row, col] = ind2sub([Ny Nx], find(blocs(:,:,sld) == idx_blocs(k)));
                    hemo_mask_frame(row, col, sld) = 1;
                end
            end
            hemo_mask_frame_cell{i,w} = hemo_mask_frame;
            hemo_mask = zeros(Ny, Nx, slides, length(sigma_array));
            for j = 1:size(C_cell, 2)
                for sld = 1:slides
                    thresh = mean(nonzeros(myo_mask .* C_cell{i,j}(:,:,sld))) - 2*std(nonzeros(myo_mask .* C_cell{i,j}(:,:,sld)));
                    hemo_mask(:,:,sld,j) = (C_cell{i,j}(:,:,sld) < thresh); %.* hemo_mask_frame(:,:,i,w);
                    %thresh_mat(:,:,i,j,w) = thresh;
                end
            end
            hemo_mask_cell{i,w,iter} = hemo_mask;
        end
    end
end

%%
hemo_mask_frame_w04 = hemo_mask_frame_cell{3,1};
hemo_mask_frame_w08 = hemo_mask_frame_cell{3,2};
hemo_mask_frame_w10 = hemo_mask_frame_cell{3,3};
hemo_mask_frame_w30 = hemo_mask_frame_cell{3,4};
hemo_mask_frame_w40 = hemo_mask_frame_cell{3,5};

figure(); 
for w = 1:size(hemo_mask_frame_w20, 3)
    subplot(8,2,2*(w-1)+1);
    imagesc(hemo_mask_frame_w04(:,:,w));
    axis off;
    subplot(8,2,2*(w-1)+2);
    imagesc(hemo_mask_frame_w30(:,:,w));
    axis off;
end
%% Sensitivity, Specificity, Accuracy, Precision
sens_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
spec_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
accu_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
prec_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
CNR_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
SNR_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
CNR_nom_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
auc_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
dice_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));

% VS truth in canvas
sens_cell2 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
spec_cell2 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
accu_cell2 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
prec_cell2 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
auc_cell2 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
dice_cell2 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));

accu_cell3 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
auc_cell3 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));

sens_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
spec_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
accu_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
prec_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
auc_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
dice_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
CNR_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
SNR_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
CNR_nom_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));

% VS truth in canvas
sens_matt2 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
spec_matt2 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
accu_matt2 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
prec_matt2 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
auc_matt2 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
dice_matt2 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));

accu_matt3 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));
auc_matt3 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_03212023.tissue_canvas));

myo = myo_mask;
for iter = 1:length(SimPhantom_03212023.tissue_canvas)
    for w = 1:length(width_max)
        t_width = width_max(w);
        C_cell = SimPhantom_03212023.tissue_canvas{iter}(w).C_cell;
        truth_map2 = t_hemo(:,:,w);
        
        for i = 1:length(res_array)
            res = res_array(i);
            bloc_size = round(res/dx).^2;
            truth_map = hemo_mask_frame_cell{i,w}; % 3D truth map
            hemo_mask = hemo_mask_cell{i,w,iter};
            %  Ny x Nx x slides x sigma
            slides = size(hemo_mask, 3);
            for j = 1:length(sigma_array)
                sens_mat = zeros(slides, 1);
                spec_mat = zeros(slides, 1);
                accu_mat = zeros(slides, 1);
                prec_mat = zeros(slides, 1);
                CNR_mat = zeros(slides, 1);
                SNR_mat = zeros(slides, 1);
                CNR_nom_mat = zeros(slides, 1);

                % VS truth in canvas
                sens_mat2 = zeros(slides, 1);
                spec_mat2 = zeros(slides, 1);
                accu_mat2 = zeros(slides, 1);
                prec_mat2 = zeros(slides, 1);
                auc_mat2 = zeros(slides, 1);
                
                accu_mat3 = zeros(slides, 1);
                auc_mat3 = zeros(slides, 1);

                for sld = 1:slides
                    FP = find(truth_map(:,:,sld) - hemo_mask(:,:,sld,j) == -1);
                    FN = find(truth_map(:,:,sld) - hemo_mask(:,:,sld,j) == 1);
                    TP = find(2*truth_map(:,:,sld) - hemo_mask(:,:,sld,j) == 1);
                    TN = find(2*(~truth_map(:,:,sld)) - ~hemo_mask(:,:,sld,j) == 1);
                    
                    hemo_mask_temp = hemo_mask(:,:,sld,j);
                    truth_map_temp = truth_map(:,:,sld);
                    [X,Y,T,AUC,OPTROCPT] = perfcurve(truth_map_temp(:), hemo_mask_temp(:) , 1);

                    fp = length(FP)/bloc_size;
                    fn = length(FN)/bloc_size;
                    tp = length(TP)/bloc_size;
                    tn = length(TN)/bloc_size;

                    sensitivity = tp ./ (tp + fn);
                    specificity = tn ./ (tn + fp);
                    accuracy = (tp + tn) ./ (fp + fn + tp + tn);
                    precision = (tp) ./ (fp + tp);
                    dice = (2*tp) ./ (2*tp + fp + fn);

                    sens_mat(sld) = sensitivity;
                    spec_mat(sld) = specificity;
                    accu_mat(sld) = accuracy;
                    prec_mat(sld) = precision;
                    auc_mat(sld) = AUC;
                    dice_mat(sld) = dice;

                    FP = find(truth_map2 - hemo_mask(:,:,sld,j) == -1);
                    FN = find(truth_map2 - hemo_mask(:,:,sld,j) == 1);
                    TP = find(2*truth_map2 - hemo_mask(:,:,sld,j) == 1);
                    TN = find(2*(~truth_map2) - ~hemo_mask(:,:,sld,j) == 1);
                    
                    [X,Y,T,AUC,OPTROCPT] = perfcurve(truth_map2(:),hemo_mask_temp(:) ,1);
                    
                    fp = length(FP)/bloc_size;
                    fn = length(FN)/bloc_size;
                    tp = length(TP)/bloc_size;
                    tn = length(TN)/bloc_size;

                    sensitivity = tp ./ (tp + fn);
                    specificity = tn ./ (tn + fp);
                    accuracy = (tp + tn) ./ (fp + fn + tp + tn);
                    precision = (tp) ./ (fp + tp);
                    dice = (2*tp) ./ (2*tp + fp + fn);

                    sens_mat2(sld) = sensitivity;
                    spec_mat2(sld) = specificity;
                    accu_mat2(sld) = accuracy;
                    prec_mat2(sld) = precision;
                    auc_mat2(sld) = AUC;
                    dice_mat2(sld) = dice;
                    
                    % [row, col] = find(truth_map2 == 1);
                    FP = find((truth_map2 - hemo_mask(:,:,sld,j)).*truth_map(:,:,sld) == -1);
                    FN = find((truth_map2 - hemo_mask(:,:,sld,j)).*truth_map(:,:,sld) == 1);
                    TP = find((2*truth_map2 - hemo_mask(:,:,sld,j)).*truth_map(:,:,sld) == 1);
                    TN = find((2*(~truth_map2) - ~hemo_mask(:,:,sld,j)).*truth_map(:,:,sld) == 1);
                    
                    fp = length(FP)/bloc_size;
                    fn = length(FN)/bloc_size;
                    tp = length(TP)/bloc_size;
                    tn = length(TN)/bloc_size;
                    
                    hemo_mask_temp = hemo_mask(:,:,sld,j).*truth_map(:,:,sld);
                    
                    [X,Y,T,AUC,OPTROCPT] = perfcurve(truth_map2(:),hemo_mask_temp(:) ,1);
                    accuracy = (tp + tn) ./ (fp + fn + tp + tn);
                    accu_mat3(sld) = accuracy;
                    auc_mat3(sld) = AUC;

                    s_hemo = mean(nonzeros(hemo_mask(:,:,sld,j) .* truth_map(:,:,sld) .* C_cell{i,j}(:,:,sld)));
                    s_hemo_nom = mean(nonzeros(truth_map(:,:,sld) .* C_cell{i,j}(:,:,sld)));
                    s_myo = mean(nonzeros(myo .* C_cell{i,j}(:,:,sld)));
                    std_myo = std(nonzeros(myo .* C_cell{i,j}(:,:,sld)));
                    SNR_mat(sld) = s_myo ./ std_myo;
                    CNR_mat(sld) = (s_myo-s_hemo) ./ std_myo;
                    % A different way to do CNR
                    CNR_nom_mat(sld) = (s_myo-s_hemo_nom) ./ std_myo;
                end

                sens_cell{i,j,w,iter} = sens_mat;
                spec_cell{i,j,w,iter} = spec_mat;
                accu_cell{i,j,w,iter} = accu_mat;
                prec_cell{i,j,w,iter} = prec_mat;
                auc_cell{i,j,w,iter} = auc_mat;
                dice_cell{i,j,w,iter} = dice_mat;
                CNR_cell{i,j,w,iter} = CNR_mat;
                SNR_cell{i,j,w,iter} = SNR_mat;
                CNR_nom_cell{i,j,w,iter} = CNR_nom_mat;

                % VS truth in canvas
                sens_cell2{i,j,w,iter} = sens_mat2;
                spec_cell2{i,j,w,iter} = spec_mat2;
                accu_cell2{i,j,w,iter} = accu_mat2;
                prec_cell2{i,j,w,iter} = prec_mat2;
                auc_cell2{i,j,w,iter} = auc_mat2;
                dice_cell2{i,j,w,iter} = dice_mat2;
                
                accu_cell3{i,j,w,iter} = accu_mat3;
                auc_cell3{i,j,w,iter} = auc_mat3;

                sens_matt(i,j,w,iter) = mean(sens_mat);
                spec_matt(i,j,w,iter) = mean(spec_mat);
                accu_matt(i,j,w,iter) = mean(accu_mat);
                prec_matt(i,j,w,iter) = mean(prec_mat);
                auc_matt(i,j,w,iter) = mean(auc_mat);
                dice_matt(i,j,w,iter) = mean(dice_mat);

                CNR_matt(i,j,w,iter) = mean(CNR_mat);
                SNR_matt(i,j,w,iter) = mean(SNR_mat);
                CNR_nom_matt(i,j,w,iter) = mean(CNR_nom_mat);

                % VS truth in canvas
                sens_matt2(i,j,w,iter) = mean(sens_mat2);
                spec_matt2(i,j,w,iter) = mean(spec_mat2);
                accu_matt2(i,j,w,iter) = mean(accu_mat2);
                prec_matt2(i,j,w,iter) = mean(prec_mat2);
                auc_matt2(i,j,w,iter) = mean(auc_mat2);
                dice_matt2(i,j,w,iter) = mean(dice_mat2);
                
                accu_matt3(i,j,w,iter) = mean(accu_mat3);
                auc_matt3(i,j,w,iter) = mean(auc_mat3);
            end
        end
    end
end
%%
SNR_mat_04 = mean(SNR_matt(:,:,1,:),4);
SNR_mat_08 = mean(SNR_matt(:,:,2,:),4);
SNR_mat_10 = mean(SNR_matt(:,:,3,:),4);
SNR_mat_30 = mean(SNR_matt(:,:,4,:),4);

accu_matt2_04 = mean(accu_matt2(:,:,1,:),4);
accu_matt2_08 = mean(accu_matt2(:,:,2,:),4);
accu_matt2_10 = mean(accu_matt2(:,:,3,:),4);
accu_matt2_30 = mean(accu_matt2(:,:,4,:),4);

sens_matt2_04 = mean(sens_matt2(:,:,1,:),4);
sens_matt2_08 = mean(sens_matt2(:,:,2,:),4);
sens_matt2_10 = mean(sens_matt2(:,:,3,:),4);
sens_matt2_30 = mean(sens_matt2(:,:,4,:),4);

auc_matt2_04 = mean(auc_matt2(:,:,1,:),4);
auc_matt2_08 = mean(auc_matt2(:,:,2,:),4);
auc_matt2_10 = mean(auc_matt2(:,:,3,:),4);
auc_matt2_30 = mean(auc_matt2(:,:,4,:),4);

dice_matt2_04 = mean(dice_matt2(:,:,1,:),4);
dice_matt2_08 = mean(dice_matt2(:,:,2,:),4);
dice_matt2_10 = mean(dice_matt2(:,:,3,:),4);
dice_matt2_30 = mean(dice_matt2(:,:,4,:),4);

auc_matt3_04 = mean(auc_matt3(:,:,1,:),4);
auc_matt3_08 = mean(auc_matt3(:,:,2,:),4);
auc_matt3_10 = mean(auc_matt3(:,:,3,:),4);
auc_matt3_30 = mean(auc_matt3(:,:,4,:),4);

accu_matt3_04 = mean(accu_matt3(:,:,1,:),4);
accu_matt3_08 = mean(accu_matt3(:,:,2,:),4);
accu_matt3_10 = mean(accu_matt3(:,:,3,:),4);
accu_matt3_30 = mean(accu_matt3(:,:,4,:),4);
%% CNR
% This is the contrast in nominal hemo
SNR_mat_04 = mean(SNR_matt(:,:,1,:),4);
SNR_mat_08 = mean(SNR_matt(:,:,2,:),4);
SNR_mat_10 = mean(SNR_matt(:,:,3,:),4);
SNR_mat_30 = mean(SNR_matt(:,:,4,:),4);

CNR_mat_04 = mean(CNR_nom_matt(:,:,1,:),4);
CNR_mat_08 = mean(CNR_nom_matt(:,:,2,:),4);
CNR_mat_10 = mean(CNR_nom_matt(:,:,3,:),4);
CNR_mat_30 = mean(CNR_nom_matt(:,:,4,:),4);
CNR_mat_40 = mean(CNR_nom_matt(:,:,5,:),4);

figure();
imagesc(mean(CNR_nom_matt(:,:,1,:),4)); axis image; axis off;
caxis([3 20]);
colormap(brewermap([],'*RdYlBu'));
%colormap(brewermap([],'*YlGnBu'));
colorbar;

[B,I] = sort(CNR_mat_04,1);
I_idx = zeros(size(I));
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        I_idx(I(i,j),j) = i;
    end
end
figure();
imagesc(I_idx);

[B,I] = sort(CNR_mat_08,1);
I_idx = zeros(size(I));
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        I_idx(I(i,j),j) = i;
    end
end
figure();
imagesc(I_idx);


[B,I] = sort(CNR_mat_10,1);
I_idx = zeros(size(I));
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        I_idx(I(i,j),j) = i;
    end
end
figure();
imagesc(I_idx);

[B,I] = sort(CNR_mat_30,1);
I_idx = zeros(size(I));
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        I_idx(I(i,j),j) = i;
    end
end
figure();
imagesc(I_idx);

%% Accuracy
figure();
imagesc(mean(accu_matt2(:,:,5,:),4)); axis image; axis off;
caxis([0.8 1]);
colormap(brewermap([],'*RdYlBu'));
%colormap(brewermap([],'*YlGnBu'));
colorbar;

accu_matt2_04 = mean(accu_matt2(:,:,1,:),4);
accu_matt2_08 = mean(accu_matt2(:,:,2,:),4);
accu_matt2_10 = mean(accu_matt2(:,:,3,:),4);
accu_matt2_30 = mean(accu_matt2(:,:,4,:),4);
accu_matt2_40 = mean(accu_matt2(:,:,5,:),4);

[B,I] = sort(accu_matt2_40,1);
I_idx = zeros(size(I));
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        I_idx(I(i,j),j) = i;
    end
end
figure();
imagesc(I_idx);

%% Sensitivity
figure();
imagesc(mean(sens_matt2(:,:,1,:),4)); axis image; axis off;
caxis([0.6 1]);
colormap(brewermap([],'*RdYlBu'));
%colormap(brewermap([],'*YlGnBu'));
colorbar;

sens_matt2_04 = mean(sens_matt2(:,:,1,:),4);
sens_matt2_08 = mean(sens_matt2(:,:,2,:),4);
sens_matt2_10 = mean(sens_matt2(:,:,3,:),4);
sens_matt2_30 = mean(sens_matt2(:,:,4,:),4);
sens_matt2_40 = mean(sens_matt2(:,:,5,:),4);

[B,I] = sort(sens_matt2_40,1);
I_idx = zeros(size(I));
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        I_idx(I(i,j),j) = i;
    end
end
figure();
imagesc(I_idx);
%% AUC2
figure();
imagesc(mean(auc_matt2(:,:,2,:),4)); axis image; axis off;
caxis([0.7 1]);
colormap(brewermap([],'*RdYlBu'));
%colormap(brewermap([],'*YlGnBu'));
colorbar;

auc_matt2_04 = mean(auc_matt2(:,:,1,:),4);
auc_matt2_08 = mean(auc_matt2(:,:,2,:),4);
auc_matt2_10 = mean(auc_matt2(:,:,3,:),4);
auc_matt2_30 = mean(auc_matt2(:,:,4,:),4);
auc_matt2_40 = mean(auc_matt2(:,:,5,:),4);

[B,I] = sort(auc_matt2_40,1);
I_idx = zeros(size(I));
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        I_idx(I(i,j),j) = i;
    end
end
figure();
imagesc(I_idx);

%% DICE2
figure();
imagesc(mean(dice_matt2(:,:,5,:),4)); axis image; axis off;
caxis([0 1]);
colormap(brewermap([],'*RdYlBu'));
%colormap(brewermap([],'*YlGnBu'));
colorbar;

dice_matt2_04 = mean(dice_matt2(:,:,1,:),4);
dice_matt2_08 = mean(dice_matt2(:,:,2,:),4);
dice_matt2_10 = mean(dice_matt2(:,:,3,:),4);
dice_matt2_30 = mean(dice_matt2(:,:,4,:),4);
dice_matt2_40 = mean(dice_matt2(:,:,5,:),4);

[B,I] = sort(dice_matt2_40,1);
I_idx = zeros(size(I));
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        I_idx(I(i,j),j) = i;
    end
end
figure();
imagesc(I_idx);
%% Putting them together
CNR_mat_compiled = double(CNR_mat_04>3) + double(CNR_mat_06>3) + double(CNR_mat_08>3) + double(CNR_mat_10>3) + double(CNR_mat_20>3) + double(CNR_mat_30>3);
CNR_mat_avg = (CNR_mat_04+CNR_mat_06+CNR_mat_08+CNR_mat_10+CNR_mat_20+CNR_mat_30)/6;

figure();
imagesc(CNR_mat_avg.');axis image; axis off;
caxis([3 20]);
colormap(brewermap([],'*RdYlBu'));
colorbar;

hold on;
x = 1:7;
y = 1:10;
[X,Y] = meshgrid(x,y);  
[C,h] = contour(X,Y,CNR_mat_compiled.', 'w', 'LineWidth', 1.5, 'ShowText', 'on')
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
clabel(C,h,'FontSize',15,'Color','white');
axis off;
%%
accu_matt2_avg = (accu_matt2_04+accu_matt2_06+accu_matt2_08+accu_matt2_10+accu_matt2_20+accu_matt2_30)/6;
figure();
imagesc(accu_matt2_avg.');axis image; axis off;
caxis([0.8 1]);
colormap(brewermap([],'*RdYlBu'));
colorbar;

[B,I] = sort(accu_matt2_04,1);
I_idx = zeros(size(I));
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        I_idx(I(i,j),j) = i;
    end
end
[idx] = find(I_idx == 7);
[row_04, col_04] = ind2sub(size(I), idx);

[B,I] = sort(accu_matt2_06,1);
I_idx = zeros(size(I));
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        I_idx(I(i,j),j) = i;
    end
end
[idx] = find(I_idx == 7);
[row_06, col_06] = ind2sub(size(I), idx)

[B,I] = sort(accu_matt2_08,1);
I_idx = zeros(size(I));
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        I_idx(I(i,j),j) = i;
    end
end
[idx] = find(I_idx == 7);
[row_08, col_08] = ind2sub(size(I), idx)

[B,I] = sort(accu_matt2_10,1);
I_idx = zeros(size(I));
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        I_idx(I(i,j),j) = i;
    end
end
[idx] = find(I_idx == 7);
[row_10, col_10] = ind2sub(size(I), idx)

[B,I] = sort(accu_matt2_20,1);
I_idx = zeros(size(I));
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        I_idx(I(i,j),j) = i;
    end
end
[idx] = find(I_idx == 7);
[row_20, col_20] = ind2sub(size(I), idx)

[B,I] = sort(accu_matt2_30,1);
I_idx = zeros(size(I));
for i = 1:size(I, 1)
    for j = 1:size(I, 2)
        I_idx(I(i,j),j) = i;
    end
end
[idx] = find(I_idx == 7);
[row_30, col_30] = ind2sub(size(I), idx)

%% 
col_04(end) = col_06(end);
row_04(end) = row_06(end);

figure(); plot(col_04,row_04, 'LineWidth', 2);
hold on;
plot(col_06+0.08,row_06+0.08, 'LineWidth', 2);
plot(col_08+0.16,row_08+0.16, 'LineWidth', 2);
plot(col_10+0.24,row_10+0.24, 'LineWidth', 2);
plot(col_20+0.32,row_20+0.32, 'LineWidth', 2);
plot(col_30+0.40,row_30+0.40, 'LineWidth', 2);
axis off;
hold off;
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
ylim([0 7]);
%%
SimPhantom_2023_analysis.SNR_mat_04 = SNR_mat_04;
SimPhantom_2023_analysis.SNR_mat_06 = SNR_mat_06;
SimPhantom_2023_analysis.SNR_mat_08 = SNR_mat_08;
SimPhantom_2023_analysis.SNR_mat_10 = SNR_mat_10;
SimPhantom_2023_analysis.SNR_mat_20 = SNR_mat_20;
SimPhantom_2023_analysis.SNR_mat_30 = SNR_mat_30;
SimPhantom_2023_analysis.res_array = res_array;
SimPhantom_2023_analysis.sigma_array = sigma_array;
SimPhantom_2023_analysis.width_max = [0.4, 0.6, 0.8, 1.0, 2.0, 3.0];
SimPhantom_2023_analysis.accu_matt2_04 = accu_matt2_04;
SimPhantom_2023_analysis.accu_matt2_06 = accu_matt2_06;
SimPhantom_2023_analysis.accu_matt2_08 = accu_matt2_08;
SimPhantom_2023_analysis.accu_matt2_10 = accu_matt2_10;
SimPhantom_2023_analysis.accu_matt2_20 = accu_matt2_20;
SimPhantom_2023_analysis.accu_matt2_30 = accu_matt2_30;

SimPhantom_2023_analysis.auc_matt2_04 = auc_matt2_04;
SimPhantom_2023_analysis.auc_matt2_06 = auc_matt2_06;
SimPhantom_2023_analysis.auc_matt2_08 = auc_matt2_08;
SimPhantom_2023_analysis.auc_matt2_10 = auc_matt2_10;
SimPhantom_2023_analysis.auc_matt2_20 = auc_matt2_20;
SimPhantom_2023_analysis.auc_matt2_30 = auc_matt2_30;

SimPhantom_2023_analysis.dice_matt2_04 = dice_matt2_04;
SimPhantom_2023_analysis.dice_matt2_06 = dice_matt2_06;
SimPhantom_2023_analysis.dice_matt2_08 = dice_matt2_08;
SimPhantom_2023_analysis.dice_matt2_10 = dice_matt2_10;
SimPhantom_2023_analysis.dice_matt2_20 = dice_matt2_20;
SimPhantom_2023_analysis.dice_matt2_30 = dice_matt2_30;

SimPhantom_2023_analysis.sens_matt2_04 = sens_matt2_04;
SimPhantom_2023_analysis.sens_matt2_06 = sens_matt2_06;
SimPhantom_2023_analysis.sens_matt2_08 = sens_matt2_08;
SimPhantom_2023_analysis.sens_matt2_10 = sens_matt2_10;
SimPhantom_2023_analysis.sens_matt2_20 = sens_matt2_20;
SimPhantom_2023_analysis.sens_matt2_30 = sens_matt2_30;

SimPhantom_2023_analysis.accu_matt3_04 = accu_matt3_04;
SimPhantom_2023_analysis.accu_matt3_06 = accu_matt3_06;
SimPhantom_2023_analysis.accu_matt3_08 = accu_matt3_08;
SimPhantom_2023_analysis.accu_matt3_10 = accu_matt3_10;
SimPhantom_2023_analysis.accu_matt3_20 = accu_matt3_20;
SimPhantom_2023_analysis.accu_matt3_30 = accu_matt3_30;

SimPhantom_2023_analysis.auc_matt3_04 = auc_matt3_04;
SimPhantom_2023_analysis.auc_matt3_06 = auc_matt3_06;
SimPhantom_2023_analysis.auc_matt3_08 = auc_matt3_08;
SimPhantom_2023_analysis.auc_matt3_10 = auc_matt3_10;
SimPhantom_2023_analysis.auc_matt3_20 = auc_matt3_20;
SimPhantom_2023_analysis.auc_matt3_30 = auc_matt3_30;

save_dir = cat(2, base_dir);
fname = 'SimPhantom_2023_analysis.mat';
save(cat(2, save_dir, '/', fname), 'SimPhantom_2023_analysis');

%% SNR vs Accuracy vs Width
SNR_mat_06 = SimPhantom_2023_analysis.SNR_mat_06;
res_array = SimPhantom_2023_analysis.res_array;
sigma_array = SimPhantom_2023_analysis.sigma_array;
width_max = SimPhantom_2023_analysis.width_max;
accu_matt2_04 = SimPhantom_2023_analysis.accu_matt2_04;
accu_matt2_06 = SimPhantom_2023_analysis.accu_matt2_06;
accu_matt2_08 = SimPhantom_2023_analysis.accu_matt2_08;
accu_matt2_10 = SimPhantom_2023_analysis.accu_matt2_10;
accu_matt2_20 = SimPhantom_2023_analysis.accu_matt2_20;
accu_matt2_30 = SimPhantom_2023_analysis.accu_matt2_30;

accu_matt2_04(6,5)
accu_matt2_06(6,5)
accu_matt2_08(6,5)
accu_matt2_10(6,5)
accu_matt2_20(6,5)
accu_matt2_30(6,5)
%%
noise_level = 10;
figure();
plot(SNR_mat_06(:,noise_level), accu_matt2_04(:,noise_level), 'LineWidth', 2);
hold on;
plot(SNR_mat_06(:,noise_level), accu_matt2_06(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_06(:,noise_level), accu_matt2_08(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_06(:,noise_level), accu_matt2_10(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_06(:,noise_level), accu_matt2_20(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_06(:,noise_level), accu_matt2_30(:,noise_level), 'LineWidth', 2);

noise_level = 5;
figure();
plot(SNR_mat_06(:,noise_level), accu_matt2_04(:,noise_level), 'LineWidth', 2);
hold on;
plot(SNR_mat_06(:,noise_level), accu_matt2_06(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_06(:,noise_level), accu_matt2_08(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_06(:,noise_level), accu_matt2_10(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_06(:,noise_level), accu_matt2_20(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_06(:,noise_level), accu_matt2_30(:,noise_level), 'LineWidth', 2);

noise_level = 1;
figure();
plot(SNR_mat_06(:,noise_level), accu_matt2_04(:,noise_level), 'LineWidth', 2);
hold on;
plot(SNR_mat_06(:,noise_level), accu_matt2_06(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_06(:,noise_level), accu_matt2_08(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_06(:,noise_level), accu_matt2_10(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_06(:,noise_level), accu_matt2_20(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_06(:,noise_level), accu_matt2_30(:,noise_level), 'LineWidth', 2);

%%
noise_level = 10;
figure('Position', [100 0 800 400]);
%figure();
plot(SNR_mat_04(:,noise_level), accu_matt2_04(:,noise_level), 'LineWidth', 2);
hold on;
plot(SNR_mat_06(:,noise_level), accu_matt2_06(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_08(:,noise_level), accu_matt2_08(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_10(:,noise_level), accu_matt2_10(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_20(:,noise_level), accu_matt2_20(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_30(:,noise_level), accu_matt2_30(:,noise_level), 'LineWidth', 2);

%figure();
noise_level = 5;
figure('Position', [100 0 800 400]);
plot(SNR_mat_04(:,noise_level), accu_matt2_04(:,noise_level), 'LineWidth', 2);
hold on;
plot(SNR_mat_06(:,noise_level), accu_matt2_06(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_08(:,noise_level), accu_matt2_08(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_10(:,noise_level), accu_matt2_10(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_20(:,noise_level), accu_matt2_20(:,noise_level), 'LineWidth', 2);
plot(SNR_mat_30(:,noise_level), accu_matt2_30(:,noise_level), 'LineWidth', 2);

%%
res_level = 6;
figure('Position', [100 0 800 400]);
%figure();
plot(SNR_mat_04(res_level,:), accu_matt2_04(res_level,:), 'LineWidth', 2);
hold on;
plot(SNR_mat_06(res_level,:), accu_matt2_06(res_level,:), 'LineWidth', 2);
plot(SNR_mat_08(res_level,:), accu_matt2_08(res_level,:), 'LineWidth', 2);
plot(SNR_mat_10(res_level,:), accu_matt2_10(res_level,:), 'LineWidth', 2);
plot(SNR_mat_20(res_level,:), accu_matt2_20(res_level,:), 'LineWidth', 2);
plot(SNR_mat_30(res_level,:), accu_matt2_30(res_level,:), 'LineWidth', 2);

res_level = 1;
figure('Position', [100 0 800 400]);
%figure();
plot(SNR_mat_04(res_level,:), accu_matt2_04(res_level,:), 'LineWidth', 2);
hold on;
plot(SNR_mat_06(res_level,:), accu_matt2_06(res_level,:), 'LineWidth', 2);
plot(SNR_mat_08(res_level,:), accu_matt2_08(res_level,:), 'LineWidth', 2);
plot(SNR_mat_10(res_level,:), accu_matt2_10(res_level,:), 'LineWidth', 2);
plot(SNR_mat_20(res_level,:), accu_matt2_20(res_level,:), 'LineWidth', 2);
plot(SNR_mat_30(res_level,:), accu_matt2_30(res_level,:), 'LineWidth', 2);