close all;
clear all;

% AUC analysis across all cases

addpath('../function/');
%%%%%%%%%%%%%%%%%%%%%%%%%% input file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aha_analysis (mat file)
% from T2star_analysis_main.m
%%%%%%%%%%%%%%%%%%%%%%%%%% output file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SNR - struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 20P10
% The starting with 0014

labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};
label = labels{5};

proj_dir = GetFullPath(cat(2, pwd, '/../../T2star_Resolution_Project/'));
if ~exist(proj_dir, 'dir')
    mkdir(proj_dir)
end

save_dir = GetFullPath(cat(2, proj_dir, 'img/'));
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

data_dir = GetFullPath(cat(2, proj_dir, 'data/'));
if ~exist(data_dir, 'dir')
    mkdir(data_dir)
end

subject_name_cell = {'18P90', '18P93', '20P03_Exvivo5', '20P10_Exvivo7', '20P11_Exvivo6', '18P92', '18P94_Exvivo3', '18P95', '17P73', '20P48'};
avg_num_cell = {'Avg0016', 'Avg0001', 'Invivo'};
%avg_num = input('Please type average number here:  ');
%if isnumeric(avg_num)
%    avg_name = cat(2, 'Avg', num2str(avg_num, '%04.f'));
%else
%    avg_name = avg_num;
%end

auc_all16 = zeros(28, length(subject_name_cell));
auc_all01 = zeros(28, length(subject_name_cell));
auc_allvivo = zeros(20, length(subject_name_cell));

for i = 1:length(subject_name_cell)
    subject_name = subject_name_cell{i};
    subject_data_dir = GetFullPath(cat(2, data_dir, subject_name, '/'));
    for j = 1:length(avg_num_cell)
        avg_name = avg_num_cell{j};
        % subject_dir = GetFullPath(cat(2, save_dir, subject_name, '/'));
        
        if j == 1
            aha16 = load(cat(2, subject_data_dir, 'aha_analysis_', avg_name, '.mat'));
            auc_all16(:,i) = aha16.aha_analysis.auc_array_mi;
        elseif j == 2
            aha01 = load(cat(2, subject_data_dir, 'aha_analysis_', avg_name, '.mat'));
            auc_all01(:,i) = aha01.aha_analysis2.auc_array_mi;
        elseif j == 3
            aha_invivo = load(cat(2, subject_data_dir, 'aha_analysis_', avg_name, '.mat'));
            auc_allvivo(:,i) = aha_invivo.aha_analysis2.auc_array_mi;
        end
    end
end

auc_avg16 = mean(auc_all16, 2);
auc_avg01 = mean(auc_all01, 2);
auc_invivo = mean(auc_allvivo, 2);
auc_std_avg16 = std(auc_all16, 0, 2);
auc_std_avg01 = std(auc_all01, 0, 2);
auc_std_invivo = std(auc_allvivo, 0, 2);
%% Read SNR metrics
base_dir = uigetdir;
folder_glob = glob(cat(2, base_dir, '\*'));
[list_to_read, order_to_read] = NamePicker(folder_glob);
whatsinit = cell(length(list_to_read), 1);
slice_data = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit{i}, slice_data{i}] = dicom23D(f);
end

vx = zeros(length(slice_data), 1);
for i = 1:length(slice_data)
    x = slice_data{i}.PixelSpacing(1);
    y = slice_data{i}.PixelSpacing(2);
    z = slice_data{i}.SliceThickness;
    vx(i) = x*y*z;
end
%% Plot AUC
[vx_sorted, I] = sort(vx);
auc_avg16_reorder = auc_avg16(I);
auc_avg01_reorder = auc_avg01(I);
auc_invivo_expand = zeros(length(auc_avg16), 1);
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
for i = 1:length(auc_invivo)
    mask_idx = mask_idx_array(i);
    auc_invivo_expand(mask_idx) = auc_invivo(i);
end
auc_invivo_reorder = auc_invivo_expand(I);
auc_invivo_reorder(auc_invivo_reorder == 0) = nan;

res_array = {'03', '06', '08', '10', '13', '16', '21'};
slc_array = {'2', '4', '6', '8'};
lg = cell(length(whatsinit), 1);
for i = 1:length(whatsinit)
    idx = I(i);
    q = fix((idx-1)/7)+1;
    r = mod(idx, 7);
    if r == 0
        real_r = 7 - r;
    else
        real_r = r;
    end
    lg{i} = cat(2, slc_array{q}, '__', res_array{real_r});
    
end
idx_array = 1:length(auc_avg16_reorder);
figure('Position', [100 0 1600 1600]);
plot(auc_avg16_reorder, 'LineWidth', 2);
hold on;
plot(auc_avg01_reorder, 'LineWidth', 2);
plot(auc_invivo_reorder, 'LineWidth', 2);
scatter(idx_array, auc_avg16_reorder, 72, 'filled', 'MarkerFaceColor', [0, 0.4470, 0.7410]);
scatter(idx_array, auc_avg01_reorder, 72, 'filled', 'MarkerFaceColor', [0.8500, 0.3250, 0.0980]);
scatter(idx_array, auc_invivo_reorder, 72, 'filled', 'MarkerFaceColor', [0.9290, 0.6940, 0.1250]);
xticks(1:length(whatsinit));
xticklabels(lg);
legend({'Avg0016', 'Avg0001', 'Invivo'});
set(gca, 'FontSize', 12);
grid on;
%% AUC two-dimensional
lg = cell(length(whatsinit), 1);
for i = 1:length(whatsinit)
    q = fix((i-1)/7)+1;
    r = mod(i, 7);
    if r == 0
        real_r = 7 - r;
    else
        real_r = r;
    end
    lg{i} = cat(2, slc_array{q}, '__', res_array{real_r});
    
end
auc_invivo_expand(auc_invivo_expand == 0) = nan;
figure('Position', [100 0 1600 1600]);
plot(auc_avg16, 'LineWidth', 2);
hold on;
plot(auc_avg01, 'LineWidth', 2);
plot(auc_invivo_expand, 'LineWidth', 2);
xticks(1:length(whatsinit));
xticklabels(lg);
legend({'Avg0016', 'Avg0001', 'Invivo'});
set(gca, 'FontSize', 12); grid on;


%% Mean + SD plot 
d = length(whatsinit)/4;
auc_avg01_reshape = reshape(auc_avg01, d, length(whatsinit)/d).';
auc_avg16_reshape = reshape(auc_avg16, d, length(whatsinit)/d).';
d2 = length(auc_invivo)/4;
auc_invivo_reshape = reshape(auc_invivo, d2, length(auc_invivo)/d2).';

auc_std_avg01_reshape = reshape(auc_std_avg01, d, length(whatsinit)/d).';
auc_std_avg16_reshape = reshape(auc_std_avg16, d, length(whatsinit)/d).';
auc_std_invivo_reshape = reshape(auc_std_invivo, d2, length(auc_invivo)/d2).';

x = [0, d, 2*d, 3*d, 4*d] + [0 , 0.5, 0.5, 0.5, 1];
inplane_res = 1:d;
res = [inplane_res; inplane_res + d; inplane_res + 2*d; inplane_res + 3*d];
res_invivo = res(:,3:end);
figure('Position', [100 0 1600 1600]);
errorbar(auc_avg16, auc_std_avg16, 'LineStyle', 'none' );
hold on;
ylim_lb = min(ylim); ylim_ub = max(ylim);
patch([x(1) x(2) x(2) x(1)], [max(ylim) max(ylim) 0 0], [241 194 151]/255, 'FaceAlpha',.5)
patch([x(2) x(3) x(3) x(2)], [max(ylim) max(ylim) 0 0], [199 213 161]/255, 'FaceAlpha',.5)
patch([x(3) x(4) x(4) x(3)], [max(ylim) max(ylim) 0 0], [159 203 219]/255, 'FaceAlpha',.5)
patch([x(4) x(5) x(5) x(4)], [max(ylim) max(ylim) 0 0], [98 141 207]/255, 'FaceAlpha',.5)
e1 = errorbar(res', auc_avg16_reshape', auc_std_avg16_reshape', '-o', 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980]); 
e2 = errorbar(res', auc_avg01_reshape', auc_std_avg01_reshape', '-o', 'LineWidth', 2, 'Color', [0.9290, 0.6940, 0.1250]); 
e3 = errorbar(res_invivo', auc_invivo_reshape', auc_std_invivo_reshape', '-o', 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]); 

xticks([1.2 4 6.5 8.5 11 13.5 15.5 18 20.5 22.5 25 27.5]);
xticklabels({'0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1'})
xlim([0 x(5)]);ylim([ylim_lb, ylim_ub])

text(2,ylim_ub*0.92, 'Slice Thickness = 2 mm', 'FontSize', 16);
text(9,ylim_ub*0.92, 'Slice Thickness = 4 mm', 'FontSize', 16);
text(16,ylim_ub*0.92, 'Slice Thickness = 6 mm', 'FontSize', 16);
text(23,ylim_ub*0.92, 'Slice Thickness = 8 mm', 'FontSize', 16);
legend([e1(1), e2(1), e3(1)], {'Avg0016', 'Avg0001', 'Invivo'});
set(gca, 'FontSize', 16);
xlabel('Resolution (mm^2)', 'FontSize', 24); ylabel('T2^* (ms)', 'FontSize', 24);
hold off
%([0.3*0.3; 0.6*0.6; 0.8*0.8; 1.0*1.0; 1.3*1.3; 1.6*1.6; 2.1*2.1] * [2 4 6 8])'
%% Below are 1 dimensional voxel-wise plot, which is not necessary so far.
%% Draw contours @ epi, endo, MI, remote, fluid
img = whatsinit{1};
myo_coords_cell = cell(size(img, 3), 2);
subject_name = input('Please type subject name here:  ', 's');
subject_data_dir = GetFullPath(cat(2, data_dir, subject_name, '/'));
roi_save = cat(2, subject_data_dir, 'roi.mat');


if ~exist(roi_save, 'file')
for i = 1:(size(img, 3))
    disp('Epicardium: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 100])
    epi = drawpolygon(gca);
    epi_coords = epi.Position;
    
    disp('Endocardium: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 100])
    endo = drawpolygon(gca);
    endo_coords = endo.Position;
    
    disp('MI roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 100])
    mi = drawpolygon(gca);
    mi_coords = mi.Position;
    
    disp('remote roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 100])
    remote = drawpolygon(gca);
    remote_coords = remote.Position;
    
    disp('fluid roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 100])
    fluid = drawpolygon(gca);
    fluid_coords = fluid.Position;
    
    disp('Center line roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 100])
    center_line = drawpolygon(gca);
    center_coords = center_line.Position;
    
    disp('air roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 100])
    air = drawpolygon(gca);
    air_coords = air.Position;
    
    myo_coords_cell{i, 1} = epi.Position;
    myo_coords_cell{i, 2} = endo.Position;
    epi_mask = createMask(epi);
    endo_mask = createMask(endo);
    center_mask = createMask(center_line);
    
    close all;
end

roi.myo_coords_cell = myo_coords_cell;
roi.mi_coords = mi_coords;
roi.remote_coords = remote_coords;
roi.fluid_coords = fluid_coords;
roi.center_coords = center_coords;
roi.air_coords = air_coords;

save(roi_save, 'roi');

else
    load(roi_save);
    myo_coords_cell = roi.myo_coords_cell;
    mi_coords = roi.mi_coords;
    remote_coords = roi.remote_coords;
    fluid_coords = roi.fluid_coords;
    center_coords = roi.center_coords;
    air_coords = roi.air_coords;
end

% 28 different set of parameters
% Convert coords to masks for 28 images
img_size_truth = size(img);
mask_save = cat(2, subject_data_dir, 'mask.mat');

if ~exist(mask_save, 'file')
    figure();
    mask_struct = struct;
    for i = 1:length(whatsinit)
        img2 = whatsinit{i};
        img2_size = size(whatsinit{i});
        ratio = img_size_truth(1) ./ img2_size(1);
        imagesc(img2); caxis([0 100]);
        epi = drawpolygon(gca,'Position', [myo_coords_cell{1}(:,1)/ratio + (ratio-1)/ratio, myo_coords_cell{1}(:,2)/ratio + (ratio-1)/ratio]);
        endo = drawpolygon(gca,'Position', [myo_coords_cell{2}(:,1)/ratio + (ratio-1)/ratio, myo_coords_cell{2}(:,2)/ratio + (ratio-1)/ratio]);
        mi = drawpolygon(gca,'Position', [mi_coords(:,1)/ratio + (ratio-1)/ratio, mi_coords(:,2)/ratio + (ratio-1)/ratio]);
        remote = drawpolygon(gca,'Position', [remote_coords(:,1)/ratio + (ratio-1)/ratio, remote_coords(:,2)/ratio + (ratio-1)/ratio]);
        fluid = drawpolygon(gca,'Position', [fluid_coords(:,1)/ratio + (ratio-1)/ratio, fluid_coords(:,2)/ratio + (ratio-1)/ratio]);
        air = drawpolygon(gca,'Position', [air_coords(:,1)/ratio + (ratio-1)/ratio, air_coords(:,2)/ratio + (ratio-1)/ratio]);
        center_line = drawpolygon(gca,'Position', [center_coords(:,1)/ratio + (ratio-1)/ratio, center_coords(:,2)/ratio + (ratio-1)/ratio]);
        
        epi_mask = createMask(epi);
        endo_mask = createMask(endo);
        myo_mask = epi_mask - endo_mask;
        
        mi_mask = createMask(mi);
        remote_mask = createMask(remote);
        fluid_mask = createMask(fluid);
        air_mask = createMask(air);
        center_mask = createMask(center_line);
        
        mask_struct(i).myo_mask = myo_mask;
        mask_struct(i).mi_mask = mi_mask;
        mask_struct(i).remote_mask = remote_mask;
        mask_struct(i).fluid_mask = fluid_mask;
        mask_struct(i).air_mask = air_mask;
        
        mask_struct(i).epi_mask = epi_mask;
        mask_struct(i).endo_mask = endo_mask;
        
        myo_mask_endo = myo_mask .* center_mask;
        myo_mask_epi = myo_mask - myo_mask_endo;
        mask_struct(i).myo_mask_endo = myo_mask_endo;
        mask_struct(i).myo_mask_epi = myo_mask_epi;
    end
    
    save(mask_save, 'mask_struct');
else
    load(mask_save);
end

%% Mean + SD plot 
[vx_sorted, I] = sort(vx);
whatsinit_reordered = cell(length(list_to_read), 1);
for i = 1:length(slice_data)
    idx = I(i);
    whatsinit_reordered{i} = whatsinit{idx};
end

t2star_mean_array = zeros(1, length(whatsinit));
t2star_sd_array = zeros(1, length(whatsinit));
for i = 1:length(whatsinit)
    idx = I(i);
    t2star_mean_array(i) = mean(nonzeros(mask_struct(idx).remote_mask .* whatsinit{idx}));
    t2star_sd_array(i) = std(nonzeros(mask_struct(idx).remote_mask .* whatsinit{idx}));
end

figure('Position', [100 0 1600 1600]);
hold on;
ylim_lb = min(ylim); ylim_ub = max(ylim);
errorbar(t2star_mean_array, t2star_sd_array, '-o', 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980]); 
%xticks([1.2 4 6.5 8.5 11 13.5 15.5 18 20.5 22.5 25 27.5]);
%xticklabels({'0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1','0.3x0.3','---->','2.1x2.1'})
xlim([0 x(5)]);ylim([ylim_lb, ylim_ub])

set(gca, 'FontSize', 16);
xlabel('Resolution (mm^2)', 'FontSize', 24); ylabel('T2^* (ms)', 'FontSize', 24);
hold off

%% Mean + SD plot (boxplot)
clear t2star_remote_array
t2star_remote_array = [];
len_array = zeros(length(whatsinit), 1);
res_array = {'03', '06', '08', '10', '13', '16', '21'};
slc_array = {'2', '4', '6', '8'};
g = {};
for i = 1:length(whatsinit)
    idx = I(i);
    temp = nonzeros(mask_struct(idx).remote_mask .* whatsinit{idx});
    len_array(i) = length(temp);
    t2star_remote_array = [t2star_remote_array; temp];
    q = fix((idx-1)/7)+1;
    r = mod(idx, 7);
    if r == 0
        real_r = 7 - r;
    else
        real_r = r;
    end
    g_temp = repmat({cat(2, slc_array{q}, '_', res_array{real_r})}, len_array(i), 1);
    g = [g; g_temp];
end

y_lim = [min(t2star_remote_array) max(t2star_remote_array)];
figure('Position', [100 0 1600 1600]);
ylim_lb = min(ylim); ylim_ub = max(ylim);
h = boxplot(t2star_remote_array,g);
set(h,'LineWidth',2); grid on;
ylim([min(t2star_remote_array), max(t2star_remote_array)]);

%% AHA
addpath('../AHA16Segment/');
Segn = 50;
Groove = 0;
aha50 = struct;
for i = 1:length(whatsinit)
    idx = I(i);
    [Segmentpix, stats, Mask_Segn] = AHASegmentation(whatsinit{idx}, mask_struct(idx).myo_mask, Segn, Groove, mask_struct(idx).endo_mask);
    aha50(i).Segmentpix = Segmentpix;
    aha50(i).Mask_Segn = Mask_Segn;
end

%figure(); imagesc(aha50(1).Mask_Segn.* mask_struct(1).myo_mask);

% AHA analysis
perc_array = zeros(length(whatsinit),Segn);

for i = 1:length(whatsinit)
    idx = I(i);
    for j = 1:Segn
        thresh = mean(nonzeros(whatsinit{idx} .* mask_struct(idx).remote_mask)) - 2*std(nonzeros(whatsinit{idx} .* mask_struct(idx).remote_mask));
        perc_array(i,j) = sum(aha50(i).Segmentpix{j}<thresh) / length(aha50(i).Segmentpix{j});
    end
end

res = perc_array > 0.1;

%% 50 Segments in MI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aha_mi2 = struct;
for i = 1:length(whatsinit)
    idx = I(i);
    Mipix = cell(Segn, 1);
    Img = whatsinit{idx};
    for j = 1:Segn
        Mipix{j,1} = Img(aha50(i).Mask_Segn .* mask_struct(idx).mi_mask == j);
    end
    aha_mi2(i).Mipix = Mipix;
end

Mipix_flat = {};
for i = 1:length(whatsinit)
    idx = 1;
        for j = 1:Segn
            if ~isempty(aha_mi2(1).Mipix{j,1})
                Mipix_flat{idx} = aha_mi2(i).Mipix{j,1};
                idx = idx + 1;
            end
        end
    aha_mi2(i).Mipix_flat = Mipix_flat;
end

% AHA analysis
perc_array_mi2 = zeros(length(whatsinit),length(Mipix_flat));

for i = 1:length(whatsinit)
    idx = I(i);
    for j = 1:length(Mipix_flat)
        thresh = mean(nonzeros(whatsinit{idx} .* mask_struct(idx).remote_mask)) - 2*std(nonzeros(whatsinit{idx} .* mask_struct(idx).remote_mask));
        perc_array_mi2(i,j) = sum(aha_mi2(i).Mipix_flat{j}<thresh) / length(aha_mi2(i).Mipix_flat{j});
    end
end
res_mi2 = perc_array_mi2 > 0.1;
%% Confusion Matrix
figure('Position', [100 0 1600 1600]);
sens_mi2 = zeros(length(whatsinit),1);
spec_mi2 = zeros(length(whatsinit),1);
for i = 1:length(whatsinit)
    [cm, order] = confusionmat(res_mi2(1,:),res_mi2(i,:));
    subplot(4,7,i);confusionchart(cm, order); 
    set(gca, 'FontSize', 18);
    sens_mi2(i) = cm(2,2) / (cm(2,2) + cm(2,1));
    spec_mi2(i) = cm(1,1) / (cm(1,1) + cm(1,2));
end

% Plot sensitivity and specificity
figure();
plot(vx_sorted, sens_mi2, 'LineWidth', 2);
title('Sensitivity');
set(gca, 'FontSize', 18); 
hold on;
plot(vx_sorted, spec_mi2, 'LineWidth', 2);
title('Specificity');
set(gca, 'FontSize', 18); 
legend({'Sensitivity', 'Specificity'});

%aha_analysis.perc_array_mi = perc_array_mi2;
%aha_analysis.sens_mi = sens_mi2;
%aha_analysis.spec_mi = spec_mi2;
%aha_analysis.aha_mi = aha_mi2;

% ROC analysis %%
auc_array_mi2 = zeros(length(whatsinit),1);
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    [X,Y,T,AUC,OPTROCPT] = perfcurve(res_mi2(1,:),perc_array_mi2(i,:), 1);
    subplot(4,7,i);
    plot(X,Y, 'LineWidth', 2);
    xlabel('FPR');
    ylabel('TPR');
    text(0.5,0.5,num2str(round(AUC, 2)),'FontSize',24);
    %title('ROC for Classification')
    set(gca, 'FontSize', 16);
    auc_array_mi2(i) = AUC;
    
end

figure('Position', [100 0 1600 1600]);
lg = cell(length(whatsinit), 1);
for i = 1:length(whatsinit)
    idx = I(i);
    q = fix((idx-1)/7)+1;
    r = mod(idx, 7);
    if r == 0
        real_r = 7 - r;
    else
        real_r = r;
    end
    lg{i} = cat(2, slc_array{q}, '__', res_array{real_r});
    plot(auc_array_mi2, 'LineWidth', 2);
end
xticks(1:length(whatsinit));
xticklabels(lg);

%aha_analysis.auc_array_mi = auc_array_mi2;
%save(cat(2, subject_data_dir, 'aha_analysis_', avg_name, '.mat'), 'aha_analysis');
%res_mi2 = perc_array_mi2 > 0.1;