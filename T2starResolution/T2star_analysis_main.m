clear all;
close all;

% the main body for T2* resolution analysis
% integration of ../LineWidth_Analysis.m

addpath('../function/');

%%%%%%%%%%%%%%%%%%%%%%%%%% output file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% roi
% mask_struct
% aha_anlysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 20P10
% The starting with 0014
base_dir = uigetdir;
folder_glob = glob(cat(2, base_dir, '\*'));

labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};

label = labels{5};
idx_array = contains(folder_glob, label);

[list_to_read, order_to_read] = NamePicker(folder_glob);

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

subject_name = input('Please type subject name here:  ', 's');
subject_dir = GetFullPath(cat(2, save_dir, subject_name, '/'));
if ~exist(subject_dir, 'dir')
    mkdir(subject_dir)
end
subject_data_dir = GetFullPath(cat(2, data_dir, subject_name, '/'));
if ~exist(subject_data_dir, 'dir')
    mkdir(subject_data_dir)
end

%% Read T2* DICOM files
whatsinit = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    whatsinit{i} = dicom23D(f);
end

%% Display images
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    subplot(4,7,i);
    imagesc(whatsinit{i}); axis image;
    caxis([0 100])
end

%% Draw contours @ epi, endo, MI, remote, fluid
img = whatsinit{1};
myo_coords_cell = cell(size(img, 3), 2);
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
end



%% 28 different set of parameters
%% Convert coords to masks for 28 images
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
%% Analysis starts here
% Histogram analysis
% Myocardium
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    forhist = nonzeros(mask_struct(i).myo_mask .* whatsinit{i});
    forhist(forhist > 100) = 100;
    subplot(4,7,i)
    histogram(forhist,20);xlabel('T2* (ms)'); ylabel('Count');
    set(gca, 'FontSize', 16); % title('0.4x0.4 mm^2')
    xlim([0 100]);
end
%% Histogram MI
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    forhist = nonzeros(mask_struct(i).mi_mask .* whatsinit{i});
    forhist(forhist > 100) = 100;
    subplot(4,7,i)
    histogram(forhist, 20);xlabel('T2* (ms)'); ylabel('Count');
    set(gca, 'FontSize', 16); % title('0.4x0.4 mm^2')
    xlim([0 100]);
end
%% Histogram Remote
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    forhist = nonzeros(mask_struct(i).remote_mask .* whatsinit{i});
    forhist(forhist > 100) = 100;
    subplot(4,7,i)
    histogram(forhist);xlabel('T2* (ms)'); ylabel('Count');
    set(gca, 'FontSize', 16); % title('0.4x0.4 mm^2')
    xlim([0 50]);
end

%% AHA
addpath('../AHA16Segment/');
Segn = 50;
Groove = 0;
aha50 = struct;
for i = 1:length(whatsinit)
    [Segmentpix, stats, Mask_Segn] = AHASegmentation(whatsinit{i}, mask_struct(i).myo_mask, Segn, Groove, mask_struct(i).endo_mask);
    aha50(i).Segmentpix = Segmentpix;
    aha50(i).Mask_Segn = Mask_Segn;
end

%figure(); imagesc(aha50(1).Mask_Segn.* mask_struct(1).myo_mask);

% AHA analysis
perc_array = zeros(length(whatsinit),Segn);

for i = 1:length(whatsinit)
    for j = 1:Segn
        thresh = mean(nonzeros(whatsinit{i} .* mask_struct(i).remote_mask)) - 2*std(nonzeros(whatsinit{i} .* mask_struct(i).remote_mask));
        perc_array(i,j) = sum(aha50(i).Segmentpix{j}<thresh) / length(aha50(i).Segmentpix{j});
    end
end

res = perc_array > 0.1;

%% Confusion Matrix
figure('Position', [100 0 1600 1600]);
sens = zeros(length(whatsinit),1);
spec = zeros(length(whatsinit),1);
for i = 1:length(whatsinit)
    [cm, order] = confusionmat(res(1,:),res(i,:));
    subplot(4,7,i);confusionchart(cm, order); 
    set(gca, 'FontSize', 18);
    sens(i) = cm(2,2) / (cm(2,2) + cm(2,1));
    spec(i) = cm(1,1) / (cm(1,1) + cm(1,2));
end

%% Plot sensitivity and specificity
sens_reshape = reshape(sens, [7, 4])';
spec_reshape = reshape(spec, [7, 4])';
figure();
subplot(2,1,1);
imagesc(sens_reshape); axis image
title('Sensitivity');
set(gca, 'FontSize', 18); colorbar;
subplot(2,1,2);
imagesc(spec_reshape); axis image
title('Specificity');
set(gca, 'FontSize', 18); colorbar;

aha_analysis = struct;
aha_analysis.perc_array = perc_array;
aha_analysis.sens = sens;
aha_analysis.spec = spec;
aha_analysis.aha50 = aha50;

%% ROC analysis
auc_array = zeros(length(whatsinit),1);
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    [X,Y,T,AUC,OPTROCPT] = perfcurve(res(1,:),perc_array(i,:), 1);
    subplot(4,7,i);
    plot(X,Y, 'LineWidth', 2)
    xlabel('FPR')
    ylabel('TPR')
    text(0.5,0.5,num2str(round(AUC, 2)),'FontSize',24)
    % title('ROC for Classification')
    set(gca, 'FontSize', 16);
    auc_array(i) = AUC;
end

aha_analysis.auc_array = auc_array;
%% Plot AUC
auc_reshape = reshape(auc_array, [7, 4])';
figure();
imagesc(auc_reshape); axis image
title('Area Under Curve');
set(gca, 'FontSize', 18); colorbar;

%% Try 100 Segment (TODO)
% center_coords = DivideMyoInHalf(myo_coords_cell{1}, myo_coords_cell{2});
% Doesn't work if shape is not close to a concentric ring
%% 
figure(); imagesc(aha50(1).Mask_Segn.* mask_struct(1).myo_mask_endo);
figure(); imagesc(aha50(1).Mask_Segn.* mask_struct(1).myo_mask_epi);
% figure(); imagesc(aha50(28).Mask_Segn.* mask_struct(28).myo_mask);
aha100 = struct;
for i = 1:length(whatsinit)
    Segpix = cell(Segn, 2);
    Img = whatsinit{i};
    for j = 1:Segn
        Segpix{j,1} = Img(aha50(i).Mask_Segn .* mask_struct(i).myo_mask_endo == j);
        Segpix{j,2} = Img(aha50(i).Mask_Segn .* mask_struct(i).myo_mask_epi == j);
    end
    aha100(i).Segpix = Segpix;
end

% AHA analysis
perc_array_100 = zeros(length(whatsinit),Segn*2);

for k = 1:2
    for i = 1:length(whatsinit)
        for j = 1:Segn
            thresh = mean(nonzeros(whatsinit{i} .* mask_struct(i).remote_mask)) - 2*std(nonzeros(whatsinit{i} .* mask_struct(i).remote_mask));
            perc_array_100(i,j+(k-1)*Segn) = sum(aha100(i).Segpix{j,k}<thresh) / length(aha100(i).Segpix{j,k});
        end
    end
end

res_100 = perc_array_100 > 0.1;

%% Confusion Matrix
figure('Position', [100 0 1600 1600]);
sens_100 = zeros(length(whatsinit),1);
spec_100 = zeros(length(whatsinit),1);
for i = 1:length(whatsinit)
    [cm, order] = confusionmat(res_100(1,:),res_100(i,:));
    subplot(4,7,i);confusionchart(cm, order); 
    set(gca, 'FontSize', 18);
    sens_100(i) = cm(2,2) / (cm(2,2) + cm(2,1));
    spec_100(i) = cm(1,1) / (cm(1,1) + cm(1,2));
end

%% Plot sensitivity and specificity
sens_reshape = reshape(sens_100, [7, 4])';
spec_reshape = reshape(spec_100, [7, 4])';
figure();
subplot(2,1,1);
imagesc(sens_reshape); axis image
title('Sensitivity');
set(gca, 'FontSize', 18); colorbar;
subplot(2,1,2);
imagesc(spec_reshape); axis image
title('Specificity');
set(gca, 'FontSize', 18); colorbar;

aha_analysis.perc_array_100 = perc_array_100;
aha_analysis.sens_100 = sens_100;
aha_analysis.spec_100 = spec_100;
aha_analysis.aha100 = aha100;

%% ROC analysis
auc_array_100 = zeros(length(whatsinit),1);
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    [X,Y,T,AUC,OPTROCPT] = perfcurve(res_100(1,:),perc_array_100(i,:), 1);
    subplot(4,7,i);
    plot(X,Y, 'LineWidth', 2)
    xlabel('FPR')
    ylabel('TPR')
    text(0.5,0.5,num2str(round(AUC, 2)),'FontSize',24)
    % title('ROC for Classification')
    set(gca, 'FontSize', 16);
    auc_array_100(i) = AUC;
end

aha_analysis.auc_array_100 = auc_array_100;
%% Plot AUC
auc_reshape = reshape(auc_array_100, [7, 4])';
figure();
imagesc(auc_reshape); axis image
title('Area Under Curve');
set(gca, 'FontSize', 18); colorbar;

%% Only looking at MI region with 100 Segment 
%% This part is optional

aha_mi = struct;

% The AHASegmentation function doesn't work in this way

%Segn_mi = 10;
%for i = 1:length(whatsinit)
%    [Segmentpix, stats, Mask_Segn] = AHASegmentation(whatsinit{i}, mask_struct(i).mi_mask, Segn_mi, Groove, mask_struct(i).endo_mask);
%    aha_mi(i).Segmentpix = Segmentpix;
%    aha_mi(i).Mask_Segn = Mask_Segn;
%end
%figure(); imagesc(aha_mi(1).Mask_Segn)



for i = 1:length(whatsinit)
    Mipix = cell(Segn, 2);
    Img = whatsinit{i};
    for j = 1:Segn
        Mipix{j,1} = Img(aha50(i).Mask_Segn .* mask_struct(i).myo_mask_endo .* mask_struct(i).mi_mask  == j);
        Mipix{j,2} = Img(aha50(i).Mask_Segn .* mask_struct(i).myo_mask_epi .* mask_struct(i).mi_mask == j);
    end
    aha_mi(i).Mipix = Mipix;
end

% Only find non-empty segments for ground truth - 0.3x0.3x2
% But it appears that there will be a lot empty segments for it coarse
% resolution correspondence.
% For 20P10 Exvivo4, 37 Segments in ground truth but only 25 segments in
% 2.1x2.1x8.
Mipix_flat = {};
for i = 1:length(whatsinit)
    idx = 1;
    for k = 1:2
        for j = 1:Segn
            if ~isempty(aha_mi(1).Mipix{j,k})
                Mipix_flat{idx} = aha_mi(i).Mipix{j,k};
                idx = idx + 1;
            end
        end
    end
    aha_mi(i).Mipix_flat = Mipix_flat;
end


% AHA analysis
perc_array_mi = zeros(length(whatsinit),length(Mipix_flat));

for i = 1:length(whatsinit)
    for j = 1:length(Mipix_flat)
        thresh = mean(nonzeros(whatsinit{i} .* mask_struct(i).remote_mask)) - 2*std(nonzeros(whatsinit{i} .* mask_struct(i).remote_mask));
        perc_array_mi(i,j) = sum(aha_mi(i).Mipix_flat{j}<thresh) / length(aha_mi(i).Mipix_flat{j});
    end
end


res_mi = perc_array_mi > 0.1;
%%
figure(); imagesc(aha50(1).Mask_Segn.* mask_struct(1).mi_mask .* mask_struct(1).myo_mask_endo);
figure(); imagesc(aha50(1).Mask_Segn.* mask_struct(1).mi_mask .* mask_struct(1).myo_mask_epi);
%% Confusion Matrix

figure('Position', [100 0 1600 1600]);
sens_mi = zeros(length(whatsinit),1);
spec_mi = zeros(length(whatsinit),1);
for i = 1:length(whatsinit)
    [cm, order] = confusionmat(res_mi(1,:),res_mi(i,:));
    subplot(4,7,i);confusionchart(cm, order); 
    set(gca, 'FontSize', 18);
    sens_mi(i) = cm(2,2) / (cm(2,2) + cm(2,1));
    spec_mi(i) = cm(1,1) / (cm(1,1) + cm(1,2));
end

%% Plot sensitivity and specificity
sens_reshape = reshape(sens_mi, [7, 4])';
spec_reshape = reshape(spec_mi, [7, 4])';
figure();
subplot(2,1,1);
imagesc(sens_reshape); axis image
title('Sensitivity');
set(gca, 'FontSize', 18); colorbar;
subplot(2,1,2);
imagesc(spec_reshape); axis image
title('Specificity');
set(gca, 'FontSize', 18); colorbar;

aha_analysis.perc_array_mi = perc_array_mi;
aha_analysis.sens_mi = sens_mi;
aha_analysis.spec_mi = spec_mi;
aha_analysis.aha_mi = aha_mi;

%% ROC analysis
auc_array_mi = zeros(length(whatsinit),1);
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    [X,Y,T,AUC,OPTROCPT] = perfcurve(res_mi(1,:),perc_array_mi(i,:), 1);
    subplot(4,7,i);
    plot(X,Y, 'LineWidth', 2)
    xlabel('FPR')
    ylabel('TPR')
    text(0.5,0.5,num2str(round(AUC, 2)),'FontSize',24)
    % title('ROC for Classification')
    set(gca, 'FontSize', 16);
    auc_array_mi(i) = AUC;
end

aha_analysis.auc_array_mi = auc_array_mi;
%% Plot AUC
auc_reshape = reshape(auc_array_mi, [7, 4])';
figure();
imagesc(auc_reshape); axis image
title('Area Under Curve');
set(gca, 'FontSize', 18); colorbar;

%% Save the aha analysis
aha_analysis_save = cat(2, subject_data_dir, 'aha_analysis.mat');
save(aha_analysis_save, 'aha_analysis');