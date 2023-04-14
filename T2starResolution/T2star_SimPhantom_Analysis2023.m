close all;
clear all;

%% Load SimPhantom_04042021.mat
addpath('../function/');

base_dir = uigetdir;
f_to_read = cat(2, base_dir, '/SimPhantom_02282023.mat');
load(f_to_read);

res_array = SimPhantom_02282023.res_array;
sigma_array = SimPhantom_02282023.sigma_array;
width_max = SimPhantom_02282023.width_max;

C_t2star_fit_reshape = SimPhantom_02282023.tissue_canvas{1}(1).C_t2star_fit_reshape;
C_t2star_fit_reshape = SimPhantom_02282023.tissue_canvas{1}(3).C_t2star_fit_reshape;
% Width x Height x Res x Sigma
figure();
imagesc(C_t2star_fit_reshape(:,:,1,1))
figure();
imagesc(C_t2star_fit_reshape(:,:,7,1))
figure();
imagesc(C_t2star_fit_reshape(:,:,1,10))
figure();
imagesc(C_t2star_fit_reshape(:,:,7,10))

%%
% Renew the hemo_mask and and myo_mask
dx = 0.1;
CNR_array = zeros(length(sigma_array), length(res_array), length(SimPhantom_02282023.tissue_canvas{1}));
SNR_array = zeros(length(sigma_array), length(res_array), length(SimPhantom_02282023.tissue_canvas{1}));
C_t2star_fit_reshape = SimPhantom_02282023.tissue_canvas{1}(1).C_t2star_fit_reshape;
Nx = size(C_t2star_fit_reshape, 2);
Ny = size(C_t2star_fit_reshape, 1);
myo_mask = zeros(size(squeeze(C_t2star_fit_reshape(:,:,:,1))));

hemo_w = res_array(end)/dx;
hemo_ww = fix(hemo_w/2);
myo_mask(:, 1:(Nx/2-3*hemo_ww), :) = ones(Ny, length(1:(Nx/2-3*hemo_ww)), size(C_t2star_fit_reshape,3));
myo_mask(:, (Nx/2+3*hemo_ww)+1:end, :) = ones(Ny, length(1:(Nx/2-3*hemo_ww)), size(C_t2star_fit_reshape, 3));
t_hemo = zeros(Ny, Nx, length(width_max));
hemo_mask = zeros(Ny, Nx, length(res_array), length(sigma_array), length(width_max), length(SimPhantom_02282023.tissue_canvas));
%thresh_mat = zeros(Ny, Nx, length(res_array), length(sigma_array), length(width_max), length(SimPhantom_02282023.tissue_canvas));
hemo_mask_frame = zeros(Ny, Nx, length(res_array), length(width_max));

for iter = 1:length(SimPhantom_02282023.tissue_canvas)
    for w = 1:length(width_max)
        % res = res_array(i);
        t_width = width_max(w);
        C_t2star_fit_reshape = SimPhantom_02282023.tissue_canvas{iter}(w).C_t2star_fit_reshape;

        % Nx_hemo = res / dx;
        hemo_w = round(t_width/dx);
        hemo_ww = fix(hemo_w/2);

        if mod(hemo_w, 2) == 0
            t_hemo(:, (Nx/2-hemo_ww):(Nx/2+hemo_ww-1),w) = 1;
        else
            t_hemo(:, (Nx/2-hemo_ww):(Nx/2+hemo_ww),w) = 1;
        end

        for i = 1:length(res_array)
            res = res_array(i);
            blocs = Func_map_to_bloc_encoded(dx, Nx, res, C_t2star_fit_reshape(:,:,i,j));
            blocs_overlap = t_hemo(:,:,w) .* blocs;
            idx_blocs = nonzeros(unique(blocs_overlap));
            for k = 1:length(idx_blocs)
                [row, col] = ind2sub([Ny Nx], find(blocs == idx_blocs(k)));
                hemo_mask_frame(row, col, i, w) = 1;
            end

            for j = 1:length(sigma_array)
                thresh = mean(nonzeros(myo_mask(:,:,i) .* C_t2star_fit_reshape(:,:,i,j))) - 2*std(nonzeros(myo_mask(:,:,i) .* C_t2star_fit_reshape(:,:,i,j)));
                hemo_mask(:,:,i,j,w,iter) = (C_t2star_fit_reshape(:,:,i,j) < thresh); %.* hemo_mask_frame(:,:,i,w);
                %thresh_mat(:,:,i,j,w) = thresh;
            end
        end
    end
end
%%
figure(); 
for w = 1:length(width_max)
    subplot(9,2,2*(w-1)+1);
    imagesc(hemo_mask_frame(:,:,1,w));
    axis off;
    subplot(9,2,2*(w-1)+2);
    imagesc(t_hemo(:,:,w))
    axis off;
end
%%
figure(); 
for i = 1:length(res_array)
    subplot(7,2,2*(i-1)+1);
    imagesc(hemo_mask_frame(:,:,i,8));
    axis off;
    subplot(7,2,2*(i-1)+2);
    imagesc(t_hemo(:,:,8))
    axis off;
end
%% Sensitivity, Specificity, Accuracy, Precision
sens_mat = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_02282023.tissue_canvas));
spec_mat = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_02282023.tissue_canvas));
accu_mat = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_02282023.tissue_canvas));
prec_mat = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_02282023.tissue_canvas));
CNR_mat = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_02282023.tissue_canvas));
SNR_mat = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_02282023.tissue_canvas));
CNR_nom_mat = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_02282023.tissue_canvas));

% VS truth in canvas
sens_mat2 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_02282023.tissue_canvas));
spec_mat2 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_02282023.tissue_canvas));
accu_mat2 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_02282023.tissue_canvas));
prec_mat2 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_02282023.tissue_canvas));
for iter = 1:length(SimPhantom_02282023.tissue_canvas)

    for w = 1:length(width_max)
        t_width = width_max(w);
        C_t2star_fit_reshape = SimPhantom_02282023.tissue_canvas{iter}(w).C_t2star_fit_reshape;
        truth_map2 = t_hemo(:,:,w);
        
        for i = 1:length(res_array)
            res = res_array(i);
            bloc_size = round(res/dx).^2;
            myo = myo_mask(:,:,i);
            truth_map = hemo_mask_frame(:,:,i,w);

            for j = 1:length(sigma_array)
                
                FP = find(truth_map - hemo_mask(:,:,i,j,w,iter) == -1);
                FN = find(truth_map - hemo_mask(:,:,i,j,w,iter) == 1);
                TP = find(2*truth_map - hemo_mask(:,:,i,j,w,iter) == 1);
                TN = find(2*(~truth_map) - ~hemo_mask(:,:,i,j,w,iter) == 1);

                fp = length(FP)/bloc_size;
                fn = length(FN)/bloc_size;
                tp = length(TP)/bloc_size;
                tn = length(TN)/bloc_size;

                sensitivity = tp ./ (tp + fn);
                specificity = tn ./ (tn + fp);
                accuracy = (tp + tn) ./ (fp + fn + tp + tn);
                precision = (tp) ./ (fp + tp);

                sens_mat(i,j,w,iter) = sensitivity;
                spec_mat(i,j,w,iter) = specificity;
                accu_mat(i,j,w,iter) = accuracy;
                prec_mat(i,j,w,iter) = precision;

                FP = find(truth_map2 - hemo_mask(:,:,i,j,w,iter) == -1);
                FN = find(truth_map2 - hemo_mask(:,:,i,j,w,iter) == 1);
                TP = find(2*truth_map2 - hemo_mask(:,:,i,j,w,iter) == 1);
                TN = find(2*(~truth_map2) - ~hemo_mask(:,:,i,j,w,iter) == 1);

                fp = length(FP)/bloc_size;
                fn = length(FN)/bloc_size;
                tp = length(TP)/bloc_size;
                tn = length(TN)/bloc_size;

                sensitivity = tp ./ (tp + fn);
                specificity = tn ./ (tn + fp);
                accuracy = (tp + tn) ./ (fp + fn + tp + tn);
                precision = (tp) ./ (fp + tp);

                sens_mat2(i,j,w,iter) = sensitivity;
                spec_mat2(i,j,w,iter) = specificity;
                accu_mat2(i,j,w,iter) = accuracy;
                prec_mat2(i,j,w,iter) = precision;

                s_hemo = mean(nonzeros(hemo_mask(:,:,i,j,w,iter) .* truth_map .* C_t2star_fit_reshape(:,:,i,j)));
                s_hemo_nom = mean(nonzeros(truth_map .* C_t2star_fit_reshape(:,:,i,j)));
                s_myo = mean(nonzeros(myo .* C_t2star_fit_reshape(:,:,i,j)));
                std_myo = std(nonzeros(myo .* C_t2star_fit_reshape(:,:,i,j)));
                SNR_mat(i,j,w,iter) = s_myo ./ std_myo;
                CNR_mat(i,j,w,iter) = (s_myo-s_hemo) ./ std_myo;
                % A different way to do CNR
                CNR_nom_mat(i,j,w,iter) = (s_myo-s_hemo_nom) ./ std_myo;
            end
        end
    end
end

%% Save first and come back after SNR/CNR analysis on ex-vivo and in-vivo data
SimPhantom02282023_analysis = struct;
SimPhantom02282023_analysis.sens_mat = sens_mat;
SimPhantom02282023_analysis.spec_mat = spec_mat;
SimPhantom02282023_analysis.accu_mat = accu_mat;
SimPhantom02282023_analysis.prec_mat = prec_mat;
SimPhantom02282023_analysis.CNR_mat = CNR_mat;
SimPhantom02282023_analysis.SNR_mat = SNR_mat;
SimPhantom02282023_analysis.CNR_nom_mat = CNR_nom_mat;
SimPhantom02282023_analysis.sens_mat2 = sens_mat2;
SimPhantom02282023_analysis.spec_mat2 = spec_mat2;
SimPhantom02282023_analysis.accu_mat2 = accu_mat2;
SimPhantom02282023_analysis.prec_mat2 = prec_mat2;

save_dir = base_dir;
fname = 'SimPhantom02282023_analysis.mat';
save(cat(2, save_dir, fname), 'SimPhantom02282023_analysis');

%% Extaustive with sliding window
%%
% Renew the hemo_mask and and myo_mask
dx = 0.1;
CNR_array = zeros(length(sigma_array), length(res_array), length(SimPhantom_02282023.tissue_canvas{1}));
SNR_array = zeros(length(sigma_array), length(res_array), length(SimPhantom_02282023.tissue_canvas{1}));
C_t2star_fit_reshape = SimPhantom_02282023.tissue_canvas{1}(1).C_t2star_fit_reshape;
Nx = size(C_t2star_fit_reshape, 2);
Ny = size(C_t2star_fit_reshape, 1);
myo_mask = zeros(size(squeeze(C_t2star_fit_reshape(:,:,:,1))));

hemo_w = res_array(end)/dx;
hemo_ww = fix(hemo_w/2);
myo_mask(:, 1:(Nx/2-3*hemo_ww), :) = ones(Ny, length(1:(Nx/2-3*hemo_ww)), size(C_t2star_fit_reshape,3));
myo_mask(:, (Nx/2+3*hemo_ww)+1:end, :) = ones(Ny, length(1:(Nx/2-3*hemo_ww)), size(C_t2star_fit_reshape, 3));
t_hemo = zeros(Ny, Nx, length(width_max));
hemo_mask_cell = cell(length(res_array),1);
%thresh_mat = zeros(Ny, Nx, length(res_array), length(sigma_array), length(width_max), length(SimPhantom_02282023.tissue_canvas));
hemo_mask_frame_cell = cell(length(res_array),1);

for iter = 1:length(SimPhantom_02282023.tissue_canvas)
    for w = 1:length(width_max)
        % res = res_array(i);
        t_width = width_max(w);
        C_t2star_fit_reshape = SimPhantom_02282023.tissue_canvas{iter}(w).C_t2star_fit_reshape;

        % Nx_hemo = res / dx;
        hemo_w = round(t_width/dx);
        hemo_ww = fix(hemo_w/2);

        if mod(hemo_w, 2) == 0
            t_hemo(:, (Nx/2-hemo_ww):(Nx/2+hemo_ww-1),w) = 1;
        else
            t_hemo(:, (Nx/2-hemo_ww):(Nx/2+hemo_ww),w) = 1;
        end

        for i = 1:length(res_array)
            res = res_array(i);
            blocs = Func_map_to_bloc_encoded2(dx, Nx, res, C_t2star_fit_reshape(:,:,i,j));
            hemo_mask_frame = zeros(Ny, Nx, length(res_array), length(width_max), size(blocs,3));
            hemo_mask = zeros(Ny, Nx, length(res_array), length(sigma_array), length(width_max), length(SimPhantom_02282023.tissue_canvas), size(blocs,3));
            for xx = 1:size(blocs, 3)
                blocs_overlap = t_hemo(:,:,w) .* blocs(:,:,xx);
                idx_blocs = nonzeros(unique(blocs_overlap));
                for k = 1:length(idx_blocs)
                    [row, col] = ind2sub([Ny Nx], find(blocs == idx_blocs(k)));
                    hemo_mask_frame(row, col, i, w, xx) = 1;
                end

                for j = 1:length(sigma_array)
                    thresh = mean(nonzeros(myo_mask(:,:,i) .* C_t2star_fit_reshape(:,:,i,j))) - 2*std(nonzeros(myo_mask(:,:,i) .* C_t2star_fit_reshape(:,:,i,j)));
                    hemo_mask(:,:,i,j,w,iter,xx) = (C_t2star_fit_reshape(:,:,i,j) < thresh); %.* hemo_mask_frame(:,:,i,w);
                    %thresh_mat(:,:,i,j,w) = thresh;
                end
                hemo_mask_cell{xx} = hemo_mask;
                hemo_mask_frame_cell{xx} = hemo_mask_frame;
            end
        end
    end
end
%% Parse Phantom SimPhantom_03112023
addpath('../function/');
base_dir = uigetdir;
f_to_read = cat(2, base_dir, '/SimPhantom_03112023.mat');
load(f_to_read);

res_array = SimPhantom_03112023.res_array;
sigma_array = SimPhantom_03112023.sigma_array;
width_max = SimPhantom_03112023.width_max;

C_cell = SimPhantom_03112023.tissue_canvas{1}(1).C_cell;
% res_array x sigma_array

figure();
for i = 1:size(C_cell{1,1}, 3)
    subplot(2,2,i);
    imagesc(C_cell{1,1}(:,:,i)); caxis([0 40])
end 

figure();
for i = 1:size(C_cell{6,6}, 3)
    subplot(4,4,i);
    imagesc(C_cell{6,1}(:,:,i)); caxis([0 40])
end 
%% Plot
SNR_mat = SimPhantom02282023_analysis.SNR_mat;
SNR_mat_06 = mean(SNR_mat(:,:,4,:),4);
SNR_mat_20 = mean(SNR_mat(:,:,8,:),4);
figure();
imagesc(mean(SNR_mat(:,:,4,:),4)); axis image; axis off;
caxis([5 30]);
colormap(brewermap([],'*RdYlBu'));
%colormap(brewermap([],'*YlGnBu'));
colorbar;


[B,I] = sort(SNR_mat_06,1);
figure();
imagesc(I);

[B,I] = sort(SNR_mat_20,1);
figure();
imagesc(I);
%%
CNR_nom_mat = SimPhantom02282023_analysis.CNR_nom_mat;
CNR_mat_06 = mean(CNR_nom_mat(:,:,4,:),4);
CNR_mat_20 = mean(CNR_nom_mat(:,:,8,:),4);
figure();
imagesc(mean(CNR_nom_mat(:,:,4,:),4)); axis image; axis off;
caxis([3 20]);
colormap(brewermap([],'*RdYlBu'));
%colormap(brewermap([],'*YlGnBu'));
colorbar;

[B,I] = sort(CNR_mat_06,1);
figure();
imagesc(I);

[B,I] = sort(CNR_mat_20,1);
figure();
imagesc(I);
%%
figure();
imagesc(mean(sens_mat2(:,:,8,:),4)); axis image; axis off;
caxis([0 1]);
colormap(brewermap([],'*RdYlBu'));
%colormap(brewermap([],'*YlGnBu'));
colorbar;
%%
xslice = [];
yslice = [];
zslice = width_max;
figure();
[X, Y, Z] = meshgrid(sigma_array, res_array, width_max);
slice(X,Y,Z,double(CNR_nom_mat>3),xslice,yslice,zslice);
xlabel('Sigma'); ylabel('Resolution');zlabel('Hemo Width')
%caxis([3 20]);

%%
figure();
for i = 1:length(res_array)
    subplot(7,2,2*(i-1)+1);
    imagesc(hemo_mask(:,:,i) .* C_t2star_fit_reshape(:,:,i,10)); axis off;
    % axis equal;
    subplot(7,2,2*(i-1)+2);
    imagesc(C_t2star_fit_reshape(:,:,i,10)); axis off; % axis equal;
end