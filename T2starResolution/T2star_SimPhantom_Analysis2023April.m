clear all;
close all;

% Parse Phantom SimPhantom_03112023
addpath('../function/');
base_dir = uigetdir;
f_to_read = cat(2, base_dir, '/SimPhantom_04012023.mat');
load(f_to_read);
%%
res_array = SimPhantom_04012023.res_array;
sigma_array = SimPhantom_04012023.sigma_array;
width_max = SimPhantom_04012023.width_max;

C_cell = SimPhantom_04012023.tissue_canvas{1}(1).C_cell;
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

%%
% Renew the hemo_mask and and myo_mask
dx = 0.1;
CNR_array = zeros(length(sigma_array), length(res_array), length(SimPhantom_04012023.tissue_canvas{1}));
SNR_array = zeros(length(sigma_array), length(res_array), length(SimPhantom_04012023.tissue_canvas{1}));
C_cell = SimPhantom_04012023.tissue_canvas{1}(1).C_cell;
Nx = size(C_cell{1,1}, 2);
Ny = size(C_cell{1,1}, 1);
myo_mask = zeros(Nx, Ny);

hemo_w = res_array(end)/dx;
hemo_ww = fix(hemo_w/2);
myo_mask(:, 1:(Nx/2-3*hemo_ww-1)) = ones(Ny, length(1:(Nx/2-3*hemo_ww-1)));
myo_mask(:, (Nx/2+3*hemo_ww)+2:end) = ones(Ny, length(1:(Nx/2-3*hemo_ww-1)));
t_hemo = zeros(Ny, Nx, length(width_max));
hemo_mask_cell = cell(length(res_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
hemo_mask_frame_cell = cell(length(res_array), length(width_max));

for iter = 1:length(SimPhantom_04012023.tissue_canvas)
    for w = 1:length(width_max)
        % res = res_array(i);
        t_width = width_max(w);
        C_cell = SimPhantom_04012023.tissue_canvas{iter}(w).C_cell;

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
            hemo_mask_fwhm = zeros(Ny, Nx, slides, length(sigma_array));
            for j = 1:size(C_cell, 2)
                for sld = 1:slides
                    thresh = mean(nonzeros(myo_mask .* C_cell{i,j}(:,:,sld))) - 2*std(nonzeros(myo_mask .* C_cell{i,j}(:,:,sld)));
                    hemo_mask(:,:,sld,j) = (C_cell{i,j}(:,:,sld) < thresh); %.* hemo_mask_frame(:,:,i,w);
                    % thresh_mat(:,:,i,j,w) = thresh;
                    % thresh_fwhm = 0.5*max(max(C_cell{i,j}(:,:,sld)));
                    thresh_fwhm = min(min(C_cell{i,j}(:,:,sld))) + 0.5*(max(max(C_cell{i,j}(:,:,sld)))-min(min(C_cell{i,j}(:,:,sld))));
                    hemo_mask_fwhm(:,:,sld,j) = (C_cell{i,j}(:,:,sld) < thresh_fwhm); %.* hemo_mask_frame(:,:,i,w);
                end
            end
            hemo_mask_cell{i,w,iter} = hemo_mask;
            hemo_mask_fwhm_cell{i,w,iter} = hemo_mask_fwhm;
        end
    end
end

%%
hemo_mask_frame_w04 = hemo_mask_frame_cell{3,1};
hemo_mask_frame_w06 = hemo_mask_frame_cell{3,2};
hemo_mask_frame_w08 = hemo_mask_frame_cell{3,3};
hemo_mask_frame_w10 = hemo_mask_frame_cell{3,4};
hemo_mask_frame_w15 = hemo_mask_frame_cell{3,5};
hemo_mask_frame_w20 = hemo_mask_frame_cell{3,6};
hemo_mask_frame_w25 = hemo_mask_frame_cell{3,7};
hemo_mask_frame_w30 = hemo_mask_frame_cell{3,8};
hemo_mask_frame_w40 = hemo_mask_frame_cell{3,9};

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
sens_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
spec_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
accu_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
prec_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
CNR_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
T2star_diff_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));

SNR_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));

SD_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas)); % April 2024

CNR_nom_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
T2star_diff_nom_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));

auc_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
dice_cell = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));

% VS truth in canvas
sens_cell2 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
spec_cell2 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
accu_cell2 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
prec_cell2 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
auc_cell2 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
dice_cell2 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));

accu_cell3 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
auc_cell3 = cell(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));

sens_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
spec_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
accu_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
prec_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
auc_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
dice_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
CNR_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
T2star_diff_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas)); % April 2024

SNR_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));

SD_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas)); % April 2024

CNR_nom_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
T2star_diff_nom_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas)); % April 2024
T2star_div_nom_matt = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas)); % April 2024


% VS truth in canvas
sens_matt2 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
spec_matt2 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
accu_matt2 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
prec_matt2 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
auc_matt2 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
dice_matt2 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));

accu_matt3 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));
auc_matt3 = zeros(length(res_array), length(sigma_array), length(width_max), length(SimPhantom_04012023.tissue_canvas));

myo = myo_mask;
for iter = 1:length(SimPhantom_04012023.tissue_canvas)
    for w = 1:length(width_max)
        t_width = width_max(w);
        C_cell = SimPhantom_04012023.tissue_canvas{iter}(w).C_cell;
        truth_map2 = t_hemo(:,:,w);
        
        for i = 1:length(res_array)
        %for i = 1:1
            res = res_array(i);
            bloc_size = round(res/dx).^2;
            truth_map = hemo_mask_frame_cell{i,w}; % 3D truth map
            hemo_mask = hemo_mask_cell{i,w,iter};
            %hemo_mask = hemo_mask_fwhm_cell{i,w,iter};
            %  Ny x Nx x slides x sigma
            slides = size(hemo_mask, 3);
            for j = 1:length(sigma_array)
            %for j = 1:1
                sens_mat = zeros(slides, 1);
                spec_mat = zeros(slides, 1);
                accu_mat = zeros(slides, 1);
                prec_mat = zeros(slides, 1);
                CNR_mat = zeros(slides, 1);
                T2star_diff_mat = zeros(slides, 1); % April 2024
                SNR_mat = zeros(slides, 1);

                SD_mat = zeros(slides, 1); % April 2024
                
                CNR_nom_mat = zeros(slides, 1);

                T2star_diff_nom_mat = zeros(slides, 1); % April 2024
                T2star_div_nom_mat = zeros(slides, 1); % April 2024

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
                    SD_mat = std_myo; % April 2024
                    CNR_mat(sld) = (s_myo-s_hemo) ./ std_myo;

                    T2star_diff_mat(sld) = s_myo - s_hemo;
                    % A different way to do CNR
                    CNR_nom_mat(sld) = (s_myo-s_hemo_nom) ./ std_myo;
                    T2star_diff_nom_mat(sld) = s_myo - s_hemo_nom;
                    T2star_div_nom_mat(sld) = s_myo ./ s_hemo_nom;
                end

                sens_cell{i,j,w,iter} = sens_mat;
                spec_cell{i,j,w,iter} = spec_mat;
                accu_cell{i,j,w,iter} = accu_mat;
                prec_cell{i,j,w,iter} = prec_mat;
                auc_cell{i,j,w,iter} = auc_mat;
                dice_cell{i,j,w,iter} = dice_mat;
                CNR_cell{i,j,w,iter} = CNR_mat;
                SNR_cell{i,j,w,iter} = SNR_mat;
                SD_cell{i,j,w,iter} = SD_mat; % April 2024
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
                T2star_diff_matt(i,j,w,iter) = mean(T2star_diff_mat); % April 2024
                SNR_matt(i,j,w,iter) = mean(SNR_mat);
                SD_matt(i,j,w,iter) = mean(SD_mat); % April 2024
                CNR_nom_matt(i,j,w,iter) = mean(CNR_nom_mat);
                T2star_diff_nom_matt(i,j,w,iter) = mean(T2star_diff_nom_mat); % April 2024
                T2star_div_nom_matt(i,j,w,iter) = mean(T2star_div_nom_mat); % April 2024

                
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
SNR_mat_06 = mean(SNR_matt(:,:,2,:),4);
SNR_mat_08 = mean(SNR_matt(:,:,3,:),4);
SNR_mat_10 = mean(SNR_matt(:,:,4,:),4);
SNR_mat_15 = mean(SNR_matt(:,:,5,:),4);
SNR_mat_20 = mean(SNR_matt(:,:,6,:),4);
SNR_mat_25 = mean(SNR_matt(:,:,7,:),4);
SNR_mat_30 = mean(SNR_matt(:,:,8,:),4);
SNR_mat_40 = mean(SNR_matt(:,:,9,:),4);

% April 2024
SD_mat_04 = mean(SD_matt(:,:,1,:),4);
SD_mat_06 = mean(SD_matt(:,:,2,:),4);
SD_mat_08 = mean(SD_matt(:,:,3,:),4);
SD_mat_10 = mean(SD_matt(:,:,4,:),4);
SD_mat_15 = mean(SD_matt(:,:,5,:),4);
SD_mat_20 = mean(SD_matt(:,:,6,:),4);
SD_mat_25 = mean(SD_matt(:,:,7,:),4);
SD_mat_30 = mean(SD_matt(:,:,8,:),4);
SD_mat_40 = mean(SD_matt(:,:,9,:),4);

accu_matt2_04 = mean(accu_matt2(:,:,1,:),4);
accu_matt2_06 = mean(accu_matt2(:,:,2,:),4);
accu_matt2_08 = mean(accu_matt2(:,:,3,:),4);
accu_matt2_10 = mean(accu_matt2(:,:,4,:),4);
accu_matt2_15 = mean(accu_matt2(:,:,5,:),4);
accu_matt2_20 = mean(accu_matt2(:,:,6,:),4);
accu_matt2_25 = mean(accu_matt2(:,:,7,:),4);
accu_matt2_30 = mean(accu_matt2(:,:,8,:),4);
accu_matt2_40 = mean(accu_matt2(:,:,9,:),4);


sens_matt2_04 = mean(sens_matt2(:,:,1,:),4);
sens_matt2_06 = mean(sens_matt2(:,:,2,:),4);
sens_matt2_08 = mean(sens_matt2(:,:,3,:),4);
sens_matt2_10 = mean(sens_matt2(:,:,4,:),4);
sens_matt2_15 = mean(sens_matt2(:,:,5,:),4);
sens_matt2_20 = mean(sens_matt2(:,:,6,:),4);
sens_matt2_25 = mean(sens_matt2(:,:,7,:),4);
sens_matt2_30 = mean(sens_matt2(:,:,8,:),4);
sens_matt2_40 = mean(sens_matt2(:,:,9,:),4);

auc_matt2_04 = mean(auc_matt2(:,:,1,:),4);
auc_matt2_06 = mean(auc_matt2(:,:,2,:),4);
auc_matt2_08 = mean(auc_matt2(:,:,3,:),4);
auc_matt2_10 = mean(auc_matt2(:,:,4,:),4);
auc_matt2_15 = mean(auc_matt2(:,:,5,:),4);
auc_matt2_20 = mean(auc_matt2(:,:,6,:),4);
auc_matt2_25 = mean(auc_matt2(:,:,7,:),4);
auc_matt2_30 = mean(auc_matt2(:,:,8,:),4);
auc_matt2_40 = mean(auc_matt2(:,:,9,:),4);

dice_matt2_04 = mean(dice_matt2(:,:,1,:),4);
dice_matt2_06 = mean(dice_matt2(:,:,2,:),4);
dice_matt2_08 = mean(dice_matt2(:,:,3,:),4);
dice_matt2_10 = mean(dice_matt2(:,:,4,:),4);
dice_matt2_15 = mean(dice_matt2(:,:,5,:),4);
dice_matt2_20 = mean(dice_matt2(:,:,6,:),4);
dice_matt2_25 = mean(dice_matt2(:,:,7,:),4);
dice_matt2_30 = mean(dice_matt2(:,:,8,:),4);
dice_matt2_40 = mean(dice_matt2(:,:,9,:),4);

auc_matt3_04 = mean(auc_matt3(:,:,1,:),4);
auc_matt3_06 = mean(auc_matt3(:,:,2,:),4);
auc_matt3_08 = mean(auc_matt3(:,:,3,:),4);
auc_matt3_10 = mean(auc_matt3(:,:,4,:),4);
auc_matt3_15 = mean(auc_matt3(:,:,5,:),4);
auc_matt3_20 = mean(auc_matt3(:,:,6,:),4);
auc_matt3_25 = mean(auc_matt3(:,:,7,:),4);
auc_matt3_30 = mean(auc_matt3(:,:,8,:),4);
auc_matt3_40 = mean(auc_matt3(:,:,9,:),4);

accu_matt3_04 = mean(accu_matt3(:,:,1,:),4);
accu_matt3_06 = mean(accu_matt3(:,:,2,:),4);
accu_matt3_08 = mean(accu_matt3(:,:,3,:),4);
accu_matt3_10 = mean(accu_matt3(:,:,4,:),4);
accu_matt3_15 = mean(accu_matt3(:,:,5,:),4);
accu_matt3_20 = mean(accu_matt3(:,:,6,:),4);
accu_matt3_25 = mean(accu_matt3(:,:,7,:),4);
accu_matt3_30 = mean(accu_matt3(:,:,8,:),4);
accu_matt3_40 = mean(accu_matt3(:,:,9,:),4);

CNR_mat_04 = mean(CNR_nom_matt(:,:,1,:),4);
CNR_mat_06 = mean(CNR_nom_matt(:,:,2,:),4);
CNR_mat_08 = mean(CNR_nom_matt(:,:,3,:),4);
CNR_mat_10 = mean(CNR_nom_matt(:,:,4,:),4);
CNR_mat_15 = mean(CNR_nom_matt(:,:,5,:),4);
CNR_mat_20 = mean(CNR_nom_matt(:,:,6,:),4);
CNR_mat_25 = mean(CNR_nom_matt(:,:,7,:),4);
CNR_mat_30 = mean(CNR_nom_matt(:,:,8,:),4);
CNR_mat_40 = mean(CNR_nom_matt(:,:,9,:),4);

T2star_diff_nom_mat_04 = mean(T2star_diff_nom_matt(:,:,1,:),4);
T2star_diff_nom_mat_06 = mean(T2star_diff_nom_matt(:,:,2,:),4);
T2star_diff_nom_mat_08 = mean(T2star_diff_nom_matt(:,:,3,:),4);
T2star_diff_nom_mat_10 = mean(T2star_diff_nom_matt(:,:,4,:),4);
T2star_diff_nom_mat_15 = mean(T2star_diff_nom_matt(:,:,5,:),4);
T2star_diff_nom_mat_20 = mean(T2star_diff_nom_matt(:,:,6,:),4);
T2star_diff_nom_mat_25 = mean(T2star_diff_nom_matt(:,:,7,:),4);
T2star_diff_nom_mat_30 = mean(T2star_diff_nom_matt(:,:,8,:),4);
T2star_diff_nom_mat_40 = mean(T2star_diff_nom_matt(:,:,9,:),4);

T2star_div_nom_mat_04 = mean(T2star_div_nom_matt(:,:,1,:),4);
T2star_div_nom_mat_06 = mean(T2star_div_nom_matt(:,:,2,:),4);
T2star_div_nom_mat_08 = mean(T2star_div_nom_matt(:,:,3,:),4);
T2star_div_nom_mat_10 = mean(T2star_div_nom_matt(:,:,4,:),4);
T2star_div_nom_mat_15 = mean(T2star_div_nom_matt(:,:,5,:),4);
T2star_div_nom_mat_20 = mean(T2star_div_nom_matt(:,:,6,:),4);
T2star_div_nom_mat_25 = mean(T2star_div_nom_matt(:,:,7,:),4);
T2star_div_nom_mat_30 = mean(T2star_div_nom_matt(:,:,8,:),4);
T2star_div_nom_mat_40 = mean(T2star_div_nom_matt(:,:,9,:),4);

%% Putting them together
CNR_mat_compiled = double(CNR_mat_04>3) + double(CNR_mat_06>3) + double(CNR_mat_08>3) + double(CNR_mat_10>3)+ double(CNR_mat_15>3) + double(CNR_mat_20>3) + double(CNR_mat_25>3) + double(CNR_mat_30>3);
CNR_mat_avg = (CNR_mat_04+CNR_mat_06+CNR_mat_08+CNR_mat_10+CNR_mat_15+CNR_mat_20+CNR_mat_25+CNR_mat_30)/8;

figure();
imagesc(CNR_mat_avg.'); axis image; axis off;
caxis([3 20]);
colormap(brewermap([],'*RdYlBu'));
colorbar;

hold on;
x = 1:7;
y = 1:10;
[X,Y] = meshgrid(x,y);  
% [C,h] = contour(X,Y,CNR_mat_compiled.', [3 4 5 6 7 8], 'k', 'LineWidth', 1.5, 'ShowText', 'on')
[C,h] = contour(X,Y,CNR_mat_avg.', [3 4 5 6 7 8], 'k', 'LineWidth', 1.5, 'ShowText', 'on')

set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
clabel(C,h,'FontSize',15,'Color','black');
axis off;

%%
figure();
imagesc(CNR_mat_30.');axis image; axis off;
caxis([3 20]);
colormap(brewermap([],'*RdYlBu'));
colorbar;

figure();
imagesc(CNR_mat_04.');axis image; axis off;
caxis([3 20]);
colormap(brewermap([],'*RdYlBu'));
colorbar;
%% Just See T2* difference
T2star_diff_mat_compiled = double(T2star_diff_nom_mat_04>3) + double(T2star_diff_nom_mat_06>3) + double(T2star_diff_nom_mat_08>3) + double(T2star_diff_nom_mat_10>3)+ double(T2star_diff_nom_mat_15>3) + double(T2star_diff_nom_mat_20>3) + double(T2star_diff_nom_mat_25>3) + double(T2star_diff_nom_mat_30>3);
T2star_diff_mat_avg = (T2star_diff_nom_mat_04+T2star_diff_nom_mat_06+T2star_diff_nom_mat_08+T2star_diff_nom_mat_10+T2star_diff_nom_mat_15+T2star_diff_nom_mat_20+T2star_diff_nom_mat_25+T2star_diff_nom_mat_30)/8;

figure();
imagesc(T2star_diff_mat_avg.');axis image; axis off;
caxis([10 19]);
colormap(brewermap([],'*RdYlBu'));
colorbar;

hold on;
x = 1:7;
y = 1:10;
[X,Y] = meshgrid(x,y);  
[C,h] = contour(X,Y,T2star_diff_mat_avg.', 'k', 'LineWidth', 1.5, 'ShowText', 'on')
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
clabel(C,h,'FontSize',15,'Color','black');
axis off;

%% Hww about SD map
%SD_mat_compiled = double(SD_mat_04>3) + double(SD_mat_06>3) + double(SD_mat_08>3) + double(SD_mat_10>3)+ double(SD_mat_15>3) + double(SD_mat_20>3) + double(SD_mat_25>3) + double(SD_mat_30>3);
SD_mat_compiled = double(SD_mat_04>4) + double(SD_mat_06>4) + double(SD_mat_08>4) + double(SD_mat_10>4)+ double(SD_mat_15>4) + double(SD_mat_20>4) + double(SD_mat_25>4) + double(SD_mat_30>4);
SD_mat_avg = (SD_mat_04+SD_mat_06+SD_mat_08+SD_mat_10+SD_mat_15+SD_mat_20+SD_mat_25+SD_mat_30)/8;

figure();
imagesc(SD_mat_avg.');axis image; axis off;
%caxis([3 20]);
colormap(brewermap([],'*RdYlBu'));
colorbar;

hold on;
x = 1:7;
y = 1:10;
[X,Y] = meshgrid(x,y);  
% [C,h] = contour(X,Y,CNR_mat_compiled.', [3 4 5 6 7 8], 'k', 'LineWidth', 1.5, 'ShowText', 'on')
[C,h] = contour(X,Y,SD_mat_avg.', [2.5 5 10 15], 'k', 'LineWidth', 1.5, 'ShowText', 'on')
% [C,h] = contour(X,Y,SD_mat_compiled.', 'k', 'LineWidth', 1.5, 'ShowText', 'on')

set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
clabel(C,h,'FontSize',15,'Color','black');
axis off;

%% Just See T2* div
T2star_div_mat_compiled = double(T2star_div_nom_mat_04>3) + double(T2star_div_nom_mat_06>3) + double(T2star_div_nom_mat_08>3) + double(T2star_div_nom_mat_10>3)+ double(T2star_div_nom_mat_15>3) + double(T2star_div_nom_mat_20>3) + double(T2star_div_nom_mat_25>3) + double(T2star_div_nom_mat_30>3);
T2star_div_mat_avg = (T2star_div_nom_mat_04+T2star_div_nom_mat_06+T2star_div_nom_mat_08+T2star_div_nom_mat_10+T2star_div_nom_mat_15+T2star_div_nom_mat_20+T2star_div_nom_mat_25+T2star_div_nom_mat_25)/8;

figure();
imagesc(T2star_div_mat_avg.');axis image; axis off;
caxis([1 2]);
colormap(brewermap([],'*RdYlBu'));
colorbar;

hold on;
x = 1:7;
y = 1:10;
[X,Y] = meshgrid(x,y);  
[C,h] = contour(X,Y,T2star_div_mat_avg.',[1.5 1.8],  'k', 'LineWidth', 1.5, 'ShowText', 'on')
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
clabel(C,h,'FontSize',15,'Color','black');
axis off;
%%
SimPhantom_2023_analysis_April.SNR_mat_04 = SNR_mat_04;
SimPhantom_2023_analysis_April.SNR_mat_06 = SNR_mat_06;
SimPhantom_2023_analysis_April.SNR_mat_08 = SNR_mat_08;
SimPhantom_2023_analysis_April.SNR_mat_10 = SNR_mat_10;
SimPhantom_2023_analysis_April.SNR_mat_15 = SNR_mat_15;
SimPhantom_2023_analysis_April.SNR_mat_20 = SNR_mat_20;
SimPhantom_2023_analysis_April.SNR_mat_25 = SNR_mat_25;
SimPhantom_2023_analysis_April.SNR_mat_30 = SNR_mat_30;

SimPhantom_2023_analysis_April.SD_mat_04 = SD_mat_04;
SimPhantom_2023_analysis_April.SD_mat_06 = SD_mat_06;
SimPhantom_2023_analysis_April.SD_mat_08 = SD_mat_08;
SimPhantom_2023_analysis_April.SD_mat_10 = SD_mat_10;
SimPhantom_2023_analysis_April.SD_mat_15 = SD_mat_15;
SimPhantom_2023_analysis_April.SD_mat_20 = SD_mat_20;
SimPhantom_2023_analysis_April.SD_mat_25 = SD_mat_25;
SimPhantom_2023_analysis_April.SD_mat_30 = SD_mat_30;


SimPhantom_2023_analysis_April.res_array = res_array;
SimPhantom_2023_analysis_April.sigma_array = sigma_array;
SimPhantom_2023_analysis_April.width_max = [0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0];
SimPhantom_2023_analysis_April.accu_matt2_04 = accu_matt2_04;
SimPhantom_2023_analysis_April.accu_matt2_06 = accu_matt2_06;
SimPhantom_2023_analysis_April.accu_matt2_08 = accu_matt2_08;
SimPhantom_2023_analysis_April.accu_matt2_10 = accu_matt2_10;
SimPhantom_2023_analysis_April.accu_matt2_15 = accu_matt2_15;
SimPhantom_2023_analysis_April.accu_matt2_20 = accu_matt2_20;
SimPhantom_2023_analysis_April.accu_matt2_25 = accu_matt2_25;
SimPhantom_2023_analysis_April.accu_matt2_30 = accu_matt2_30;

SimPhantom_2023_analysis_April.auc_matt2_04 = auc_matt2_04;
SimPhantom_2023_analysis_April.auc_matt2_06 = auc_matt2_06;
SimPhantom_2023_analysis_April.auc_matt2_08 = auc_matt2_08;
SimPhantom_2023_analysis_April.auc_matt2_10 = auc_matt2_10;
SimPhantom_2023_analysis_April.auc_matt2_15 = auc_matt2_15;
SimPhantom_2023_analysis_April.auc_matt2_20 = auc_matt2_20;
SimPhantom_2023_analysis_April.auc_matt2_25 = auc_matt2_25;
SimPhantom_2023_analysis_April.auc_matt2_30 = auc_matt2_30;

SimPhantom_2023_analysis_April.dice_matt2_04 = dice_matt2_04;
SimPhantom_2023_analysis_April.dice_matt2_06 = dice_matt2_06;
SimPhantom_2023_analysis_April.dice_matt2_08 = dice_matt2_08;
SimPhantom_2023_analysis_April.dice_matt2_10 = dice_matt2_10;
SimPhantom_2023_analysis_April.dice_matt2_15 = dice_matt2_15;
SimPhantom_2023_analysis_April.dice_matt2_20 = dice_matt2_20;
SimPhantom_2023_analysis_April.dice_matt2_25 = dice_matt2_25;
SimPhantom_2023_analysis_April.dice_matt2_30 = dice_matt2_30;

SimPhantom_2023_analysis_April.sens_matt2_04 = sens_matt2_04;
SimPhantom_2023_analysis_April.sens_matt2_06 = sens_matt2_06;
SimPhantom_2023_analysis_April.sens_matt2_08 = sens_matt2_08;
SimPhantom_2023_analysis_April.sens_matt2_10 = sens_matt2_10;
SimPhantom_2023_analysis_April.sens_matt2_15 = sens_matt2_15;
SimPhantom_2023_analysis_April.sens_matt2_20 = sens_matt2_20;
SimPhantom_2023_analysis_April.sens_matt2_25 = sens_matt2_25;
SimPhantom_2023_analysis_April.sens_matt2_30 = sens_matt2_30;

SimPhantom_2023_analysis_April.accu_matt3_04 = accu_matt3_04;
SimPhantom_2023_analysis_April.accu_matt3_06 = accu_matt3_06;
SimPhantom_2023_analysis_April.accu_matt3_08 = accu_matt3_08;
SimPhantom_2023_analysis_April.accu_matt3_10 = accu_matt3_10;
SimPhantom_2023_analysis_April.accu_matt3_15 = accu_matt3_15;
SimPhantom_2023_analysis_April.accu_matt3_20 = accu_matt3_20;
SimPhantom_2023_analysis_April.accu_matt3_25 = accu_matt3_25;
SimPhantom_2023_analysis_April.accu_matt3_30 = accu_matt3_30;

SimPhantom_2023_analysis_April.auc_matt3_04 = auc_matt3_04;
SimPhantom_2023_analysis_April.auc_matt3_06 = auc_matt3_06;
SimPhantom_2023_analysis_April.auc_matt3_08 = auc_matt3_08;
SimPhantom_2023_analysis_April.auc_matt3_10 = auc_matt3_10;
SimPhantom_2023_analysis_April.auc_matt3_15 = auc_matt3_15;
SimPhantom_2023_analysis_April.auc_matt3_20 = auc_matt3_20;
SimPhantom_2023_analysis_April.auc_matt3_25 = auc_matt3_25;
SimPhantom_2023_analysis_April.auc_matt3_30 = auc_matt3_30;

save_dir = cat(2, base_dir);
%fname = 'SimPhantom_2023_analysis_April.mat';
%fname = 'SimPhantom_2023_analysis_April_fwhm_V2.mat';
fname = 'SimPhantom_2023_analysis_April_2024.mat';
save(cat(2, save_dir, '/', fname), 'SimPhantom_2023_analysis_April');
