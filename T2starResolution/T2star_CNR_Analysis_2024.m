clear all;
close all;

% the main body for T2* SNR analysis
%%
addpath('../function/');
%%%%%%%%%%%%%%%%%%%%%%%%%% input file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% roi 
% mask_struct
% Both from T2star_analysis_main.m
%%%%%%%%%%%%%%%%%%%%%%%%%% output file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SNR - struct
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

avg_num = input('Please type average number here:  ');
if isnumeric(avg_num)
    avg_name = cat(2, 'Avg', num2str(avg_num, '%04.f'));
else
    avg_name = avg_num;
end
%% Read T2* mapping DICOM files
whatsinit = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    whatsinit{i} = dicom23D(f);
end

% Display images 
figure('Position', [100 0 1600 1600]);
row = 4;
col = length(whatsinit) / row;
for i = 1:length(whatsinit)
    subplot(row,col,i);
    imagesc(whatsinit{i}); axis image;
    caxis([0 100])
end


% Draw contours @ epi, endo, MI, remote, fluid
% Coords and Masks should be generated already in data folder
% But I'm still not assuming that
% images size is 3-dimensional
img = whatsinit{1};
myo_coords_cell = cell(size(img, 4), 2);
roi_save = cat(2, subject_data_dir, 'roi.mat');

if ~exist(roi_save, 'file')
for i = 1:(size(img, 4))
    disp('Epicardium: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,1,i)); axis image;
    %caxis([0 100])
    epi = drawpolygon(gca);
    epi_coords = epi.Position;
    
    disp('Endocardium: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,1,i)); axis image;
    %caxis([0 100])
    endo = drawpolygon(gca);
    endo_coords = endo.Position;
    
    disp('MI roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,1,i)); axis image;
    %caxis([0 100])
    mi = drawpolygon(gca);
    mi_coords = mi.Position;
    
    disp('remote roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,1,i)); axis image;
    %caxis([0 100])
    remote = drawpolygon(gca);
    remote_coords = remote.Position;
    
    disp('fluid roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,1,i)); axis image;
    %caxis([0 100])
    fluid = drawpolygon(gca);
    fluid_coords = fluid.Position;
    
    disp('Center line roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,1,i)); axis image;
    %caxis([0 100])
    center_line = drawpolygon(gca);
    center_coords = center_line.Position;
    
    disp('air roi: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,1,i)); axis image;
    %caxis([0 100])
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

% Convert coords to masks for 28 images
img_size_truth = size(img);
mask_save = cat(2, subject_data_dir, 'mask.mat');


if ~exist(mask_save, 'file')
    figure();
    mask_struct = struct;
    for i = 1:length(whatsinit)
        img2 = whatsinit{i};
        img2_size = size(img2);
        ratio = img_size_truth(1) ./ img2_size(1);
        imagesc(img2(:,:,1)); % caxis([0 100]);
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

%% Avg 16 only
img_size = size(whatsinit{1});
snr_remote = zeros(length(whatsinit), 1);
cnr_remote = zeros(length(whatsinit), 1);
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
mi_mask_truth = mask_struct(1).mi_mask;
remote_mask_truth = mask_struct(1).remote_mask;
img_size_truth = size(mask_struct(1).mi_mask);

row = 4;
col = length(whatsinit) / row;
figure();
for i = 1:length(whatsinit)
%for i = 28:28
    img = whatsinit{i};
    img2_size = size(mask_struct(i).mi_mask);
    ratio = img2_size(1) ./ img_size_truth(1);
    
    if i == 1
        thresh = mean(nonzeros(img .* mask_struct(1).remote_mask)) - 2*std(nonzeros(img .* mask_struct(1).remote_mask));
        hemo_mask_truth = img < thresh;
        hemo_mask_truth = hemo_mask_truth.*mask_struct(i).mi_mask;
    end

    hemo_mask_downsampled = imresize(hemo_mask_truth, ratio, 'bicubic') >= 0.25;

    subplot(row,col,i);
    imagesc(hemo_mask_downsampled);
    for j = 1:1
        
        if ~strcmp('Invivo', avg_name)
            idx_remote = find(mask_struct(i).remote_mask == 1);
            idx_hemo = find(hemo_mask_downsampled .* mask_struct(i).mi_mask  == 1);
            remote = mask_struct(i).remote_mask .* img(:,:,j);
            hemo = hemo_mask_downsampled .* img(:,:,j);
        else
            mask_idx = mask_idx_array(i);
            idx_remote = find(mask_struct(mask_idx).remote_mask == 1);
            idx_hemo = find(hemo_mask_downsampled .* mask_struct(mask_idx).mi_mask  == 1);
            remote = mask_struct(mask_idx).remote_mask .* img(:,:,j);
            hemo = hemo_mask_downsampled .* img(:,:,j);

        end
        
        snr_remote(i,j) = mean(remote(idx_remote)) / std(remote(idx_remote));
        cnr_remote(i,j) = (mean(remote(idx_remote)) - mean(hemo(idx_hemo))) / std(remote(idx_remote));

    end
end

%% Pull up images
len = 36;
img_size_truth = size(mask_struct(1).mi_mask);

if avg_num ~= 16
    idx_array = [4, 14, 20];
    mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
else
    idx_array = [1, 6, 20, 28];
    mask_idx_array = 1:length(whatsinit);
end

%figure();
base_x = size(mask_struct(7).myo_mask, 1);

se = strel('disk',1);

for i = 1:length(idx_array)
    %for i = 28:28
    idx = idx_array(i);
    img = whatsinit{idx};
    img = img .* mask_struct(mask_idx_array(idx)).epi_mask;


    %figure();
    %imagesc(img); caxis([0 100]); axis image; %colormap default; %colorbar;
    %axis off;
    %colormap(brewermap([],'RdBu'));
    %c = colorbar;
    %w = c.FontSize;
    %c.FontSize = 20;

    colormap_dir = cat(2, subject_dir, 'ColormapOverlaid_', avg_name, '/');
    if ~exist(colormap_dir, 'dir')
        mkdir(colormap_dir)
    end
    %saveas(gcf, cat(2, colormap_dir, num2str(idx), '.png'));


    stats = regionprops(mask_struct(mask_idx_array(idx)).myo_mask);
    centroid = round(stats.Centroid);

    ioi_x = size(mask_struct(mask_idx_array(idx)).myo_mask, 1);
    multiplier = ioi_x / base_x;
    len_mul = multiplier * len;
    img2_cropped = imcrop(img, [centroid(1) - len_mul/2, centroid(2) - len_mul/2, len_mul, len_mul]);
    figure();
    ax1 = axes;
    imagesc(img2_cropped); pbaspect([size(img2_cropped, 2), size(img2_cropped, 1) 1]);
    % imagesc(img2_cropped);
    caxis([0 100]); axis image;
    axis off;
    %colormap(ax1, 'gray');
    colormap(brewermap([],'RdBu'));
    set(ax1, 'xticklabel', []); set(ax1, 'yticklabel', []);
    


    img2_size = size(mask_struct(mask_idx_array(idx)).mi_mask);
    ratio = img2_size(1) ./ img_size_truth(1);

    if avg_num == 16
        if i == 1
            thresh = mean(nonzeros(img .* mask_struct(1).remote_mask)) - 2*std(nonzeros(img .* mask_struct(1).remote_mask));
            hemo_mask_truth = img < thresh;
            hemo_mask_truth = hemo_mask_truth.*mask_struct(i).mi_mask.*mask_struct(i).epi_mask;
            %hemo_mask_truth = imopen(hemo_mask_truth,se);
        end
    end

    hemo_mask_downsampled = imresize(hemo_mask_truth, ratio, 'bicubic') >= 0.25;

    if strcmp(subject_name, '18P95')
        [x,y] = find(hemo_mask_downsampled);
        hemo_mask_downsampled_v2 = zeros(size(hemo_mask_downsampled));
        for pts = 1:length(x)
            hemo_mask_downsampled_v2(x(pts),y(pts)+1) = 1;
        end
        hemo_mask_downsampled = hemo_mask_downsampled_v2;
    end

    hemo_mask_downsampled_cropped = imcrop(hemo_mask_downsampled, [centroid(1) - len_mul/2, centroid(2) - len_mul/2, len_mul, len_mul]);

    temp_double = double(hemo_mask_downsampled_cropped);
    temp_double(temp_double == 0) = nan;
    
    
    ax2 = axes;
    imagesc(ax2, temp_double .* hemo_mask_downsampled_cropped, 'AlphaData', 0.9*hemo_mask_downsampled_cropped);
    pbaspect([size(img2_cropped, 2), size(img2_cropped, 1) 1]); colormap(ax2, 'cool');
    caxis(ax1, [0 100]); caxis(ax2, [1 2]); linkprop([ax1 ax2], 'Position');
    ax2.Visible = 'off';


    %imagesc(hemo_mask_downsampled_cropped); axis image; axis off;

    saveas(gcf, cat(2, colormap_dir, num2str(idx), '_cropped.png'));

end

%% Pull up images Mean-2SD
len = 36;
img_size_truth = size(mask_struct(1).mi_mask);

if avg_num ~= 16
    idx_array = [4, 14, 20];
    mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
else
    idx_array = [1, 6, 20, 28];
    mask_idx_array = 1:length(whatsinit);
end

%figure();
base_x = size(mask_struct(7).myo_mask, 1);

se = strel('disk',1);

%for i = 1:length(idx_array)
for i = 1:1
    idx = idx_array(i);
    img = whatsinit{idx};
    img = img .* mask_struct(mask_idx_array(idx)).epi_mask;


    %figure();
    %imagesc(img); caxis([0 100]); axis image; %colormap default; %colorbar;
    %axis off;
    %colormap(brewermap([],'RdBu'));
    %c = colorbar;
    %w = c.FontSize;
    %c.FontSize = 20;

    colormap_dir = cat(2, subject_dir, 'ColormapOverlaid_Mean2SD_', avg_name, '/');
    if ~exist(colormap_dir, 'dir')
        mkdir(colormap_dir)
    end
    %saveas(gcf, cat(2, colormap_dir, num2str(idx), '.png'));


    stats = regionprops(mask_struct(mask_idx_array(idx)).myo_mask);
    centroid = round(stats.Centroid);

    ioi_x = size(mask_struct(mask_idx_array(idx)).myo_mask, 1);
    multiplier = ioi_x / base_x;
    len_mul = multiplier * len;
    img2_cropped = imcrop(img, [centroid(1) - len_mul/2, centroid(2) - len_mul/2, len_mul, len_mul]);
    figure();
    ax1 = axes;
    imagesc(img2_cropped); pbaspect([size(img2_cropped, 2), size(img2_cropped, 1) 1]);
    % imagesc(img2_cropped);
    caxis([0 100]); axis image;
    axis off;
    %colormap(ax1, 'gray');
    colormap(brewermap([],'RdBu'));
    set(ax1, 'xticklabel', []); set(ax1, 'yticklabel', []);
    

    img2_size = size(mask_struct(mask_idx_array(idx)).mi_mask);
    ratio = img2_size(1) ./ img_size_truth(1);

    thresh = mean(nonzeros(img .* mask_struct(mask_idx_array(idx)).remote_mask)) - 2*std(nonzeros(img .* mask_struct(mask_idx_array(idx)).remote_mask));
    hemo_mask = img < thresh;
    hemo_mask = hemo_mask.*mask_struct(mask_idx_array(idx)).mi_mask.*mask_struct(mask_idx_array(idx)).epi_mask;
    %hemo_mask_truth = imopen(hemo_mask_truth,se);

    hemo_mask_cropped = imcrop(hemo_mask, [centroid(1) - len_mul/2, centroid(2) - len_mul/2, len_mul, len_mul]);

    
    if strcmp(avg_num, 'Invivo') && i == 1 && strcmp(subject_name, '20P10_Exvivo7')
        hemo_mask_cropped(9,21) = 1;
        hemo_mask_cropped(9,23) = 1;
        hemo_mask_cropped(8,22) = 1;
    elseif strcmp(avg_num, 'Invivo') && i == 1 && strcmp(subject_name, '18P95')
        hemo_mask_cropped(16,43) = 1;
    end

    temp_double = double(hemo_mask_cropped);
    temp_double(temp_double == 0) = nan;
    
    ax2 = axes;
    imagesc(ax2, temp_double .* hemo_mask_cropped, 'AlphaData', 0.9*hemo_mask_cropped);
    pbaspect([size(img2_cropped, 2), size(img2_cropped, 1) 1]); colormap(ax2, 'cool');
    caxis(ax1, [0 100]); caxis(ax2, [0 1]); linkprop([ax1 ax2], 'Position');
    ax2.Visible = 'off';


    %imagesc(hemo_mask_downsampled_cropped); axis image; axis off;
    saveas(gcf, cat(2, colormap_dir, num2str(idx), '_cropped.png'));

end

%% CMR only
snr_remote = zeros(length(whatsinit), 1);
cnr_remote = zeros(length(whatsinit), 1);
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
img_size_truth = size(mask_struct(1).mi_mask);

row = 4;
col = length(whatsinit) / row;
figure();
for i = 1:length(whatsinit)
    mask_idx = mask_idx_array(i);
    img = whatsinit{i};
    img2_size = size(mask_struct(mask_idx).mi_mask);
    ratio = img2_size(1) ./ img_size_truth(1);
    
    hemo_mask_downsampled = imresize(hemo_mask_truth, ratio, 'bicubic') >= 0.25;

    subplot(row,col,i);
    imagesc(hemo_mask_downsampled);

    for j = 1:1
        
        if ~strcmp('Invivo', avg_name)
            idx_remote = find(mask_struct(i).remote_mask == 1);
            idx_hemo = find(hemo_mask_downsampled .* mask_struct(i).mi_mask  == 1);
            remote = mask_struct(i).remote_mask .* img(:,:,j);
            hemo = hemo_mask_downsampled .* img(:,:,j);
        else
            idx_remote = find(mask_struct(mask_idx).remote_mask == 1);
            idx_hemo = find(hemo_mask_downsampled .* mask_struct(mask_idx).mi_mask  == 1);
            remote = mask_struct(mask_idx).remote_mask .* img(:,:,j);
            hemo = hemo_mask_downsampled .* img(:,:,j);
        end
        
        snr_remote(i,j) = mean(remote(idx_remote)) / std(remote(idx_remote));
        cnr_remote(i,j) = (mean(remote(idx_remote)) - mean(hemo(idx_hemo))) / std(remote(idx_remote));
    end
end


%% CMR only (Trial)
snr_remote = zeros(length(whatsinit), 1);
cnr_remote = zeros(length(whatsinit), 1);
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
img_size_truth = size(mask_struct(1).mi_mask);

% 18P95
row = 4;
col = length(whatsinit) / row;
figure();
for i = 1:length(whatsinit)
    mask_idx = mask_idx_array(i);
    img = whatsinit{i};
    img2_size = size(mask_struct(mask_idx).mi_mask);
    ratio = img2_size(1) ./ img_size_truth(1);
    
    hemo_mask_downsampled = imresize(hemo_mask_truth, ratio, 'bicubic') >= 0.25;
    [x,y] = find(hemo_mask_downsampled);
    hemo_mask_downsampled_v2 = zeros(size(hemo_mask_downsampled));
    for pts = 1:length(x)
        hemo_mask_downsampled_v2(x(pts),y(pts)+1) = 1;
    end

    subplot(row,col,i);
    imagesc(hemo_mask_downsampled_v2);

    for j = 1:1
        
        if ~strcmp('Invivo', avg_name)
            idx_remote = find(mask_struct(i).remote_mask == 1);
            idx_fluid = find(mask_struct(i).fluid_mask == 1);
            idx_hemo = find(hemo_mask_downsampled_v2 .* mask_struct(i).mi_mask  == 1);
            remote = mask_struct(i).remote_mask .* img(:,:,j);
            fluid = mask_struct(i).fluid_mask .* img(:,:,j);
            hemo = hemo_mask_downsampled_v2 .* img(:,:,j);
        else
            idx_remote = find(mask_struct(mask_idx).remote_mask == 1);
            idx_fluid = find(mask_struct(mask_idx).fluid_mask == 1);
            idx_air = find(mask_struct(mask_idx).air_mask == 1);
            idx_hemo = find(hemo_mask_downsampled_v2 .* mask_struct(mask_idx).mi_mask  == 1);
            remote = mask_struct(mask_idx).remote_mask .* img(:,:,j);
            fluid = mask_struct(mask_idx).fluid_mask .* img(:,:,j);
            air = mask_struct(mask_idx).air_mask .* img(:,:,j);
            hemo = hemo_mask_downsampled_v2 .* img(:,:,j);
        end
        
        snr_remote(i,j) = mean(remote(idx_remote)) / std(remote(idx_remote));
        cnr_remote(i,j) = (mean(remote(idx_remote)) - mean(hemo(idx_hemo))) / std(remote(idx_remote));
    end
end

%% CMR only V2 add a constraint
snr_remote = zeros(length(whatsinit), 1);
cnr_remote = zeros(length(whatsinit), 1);
mask_idx_array = [3:7, [3:7]+7, [3:7]+7*2, [3:7]+7*3];
img_size_truth = size(mask_struct(1).mi_mask);

row = 4;
col = length(whatsinit) / row;
figure();
for i = 1:length(whatsinit)
    mask_idx = mask_idx_array(i);
    img = whatsinit{i};
    img2_size = size(mask_struct(mask_idx).mi_mask);
    ratio = img2_size(1) ./ img_size_truth(1);
    
    hemo_mask_downsampled = imresize(hemo_mask_truth, ratio, 'bicubic') >= 0.25;

    
    thresh = mean(nonzeros(img .* mask_struct(mask_idx).remote_mask)) - 2*std(nonzeros(img .* mask_struct(mask_idx).remote_mask));
    hemo_mask_var = img < thresh;
    hemo_mask_var = hemo_mask_var.*mask_struct(mask_idx).mi_mask;

    subplot(row,col,i);
    imagesc(hemo_mask_downsampled + hemo_mask_var);

    hemo_union = hemo_mask_downsampled + hemo_mask_var >=1;

    for j = 1:1
        
        if ~strcmp('Invivo', avg_name)
            idx_remote = find(mask_struct(i).remote_mask == 1);
            idx_hemo = find((hemo_mask_downsampled .* mask_struct(i).mi_mask)  == 1);
            remote = mask_struct(i).remote_mask .* img(:,:,j);
            hemo = hemo_mask_downsampled .* img(:,:,j);
        else
            idx_remote = find(mask_struct(mask_idx).remote_mask == 1);
            idx_mi = find(mask_struct(mask_idx).mi_mask == 1);
            idx_hemo = find((hemo_mask_downsampled .* mask_struct(mask_idx).mi_mask)  == 1);
            remote = mask_struct(mask_idx).remote_mask .* img(:,:,j);
            mi = mask_struct(mask_idx).mi_mask .* img(:,:,j);
            hemo = hemo_mask_downsampled .* img(:,:,j);
        end
        
        snr_remote(i,j) = mean(remote(idx_remote)) / std(remote(idx_remote));
        cnr_remote(i,j) = (mean(remote(idx_remote)) - mean(hemo(idx_hemo))) / std(remote(idx_remote));
    end
end

%% Try different remote

% Display images 
figure('Position', [100 0 1600 1600]);
row = 4;
col = length(whatsinit) / row;
for i = 1:length(whatsinit)
    img = whatsinit{i};
    mask_idx = mask_idx_array(i);
    subplot(row,col,i);
    imagesc(whatsinit{i} .* mask_struct(mask_idx).remote_mask); axis image;
    caxis([20 40])

    remote = mask_struct(mask_idx).remote_mask .* img(:,:,j);
    std(nonzeros(remote))
end

%%

cnr_cmr_mat = [0.132508958	0.520975116	0.516088398	0.237064541	0	0.319358585	0	0	0.615985053	0.36954666;...
    0.76840838	0.800871453	1.271419772	0.565877629	0.226765015	0	0	0	0.625588002	0;...
    1.187722043	1.011085438	2.766853618	0.672265462	0	0	0	0.661900967	0.322161146	0;...
    2.287996043	1.847311301	1.9775832	0.473635758	0	0	1.132102603	2.819262543	1.184752477	0;...
    0.553316745	1.397150032	0.54001827	0.26788407	0	0	0	0	0	0;...
    1.152418504	0.851482288	0.556985635	0	0	0	0	0	1.068317222	0;...
    1.718560005	1.114902159	1.339517594	0	0.073052846	0	0	0	1.757454556	0;...
    1.221563782	2.13563076	2.766117771	0.110779934	0	0	0.366923592	0.76015319	0.571349563	0;...
    3.241870364	2.32013459	1.919391628	0.510180264	0	0	0.234261832	2.050657215	1.981203602	0;...
    1.661259914	1.277793587	0	1.250965325	0	0	0	0	0.201171049	0;...
    1.196781083	0.971938572	0.405856572	0	0	0	0	0	0.949517703	0;...
    1.589127782	0.911707516	2.156875153	0	0.73123969	0	0	0	1.786179208	0;...
    1.454992124	2.404627089	1.966949415	0	0	0	0	0.096409021	0.877380321	0;...
    3.745534075	2.241826125	0.318086601	0	0	0	0	1.073900848	1.653468441	0;...
    0	1.630237675	0	0	0	0	0	0	0	0;...
    1.422096786	0.74011047	0	0	0.481042999	0	0	0	0.870633123	0;...
    0.856296103	0.363264181	0.927837575	0	1.017118451	0	1.161092641	0	1.37407163	0;...
    0.706880824	1.524345221	0.884054231	0	0	0	0	0	0.242653216	0;...
    3.05809111	1.113104099	0	0	0	0	0	0.524907632	1.12029768	0;...
    0.872228416	1.383399112	0	0	0	0	0	0	0	0];
labeled_mat = zeros(size(cnr_cmr_mat));

label1 = cnr_cmr_mat > 0;
labeled_mat = double(cnr_cmr_mat > 2) + double(cnr_cmr_mat > 1) + double(cnr_cmr_mat > 0);


mean(labeled_mat)

%%

cnr2_avg16 = [10	9	7	7	7	6	5	6	6	5	6	6	6	4	4	5	2	4	4	3	2	3	3	3	3	3	2	1];
cnr1_avg16 = [10	9	9	9	8	7	6	6	8	7	8	7	6	5	6	5	5	5	5	4	4	5	5	4	4	4	4	3];

cnr2_cmr = [0	0	1	2	0	0	0	3	2	0	0	2	1	1	0	1	0	0	0	0];
cnr1_cmr = [0	1	3	6	1	2	4	5	3	3	1	4	3	3	1	3	3	1	1	1];

voxel_size = [1.28	2	3.38	5.12	8.82	2.56	4	6.76	10.24	17.64	3.84	6	10.14	15.36	26.46	5.12	8	13.52	20.48	35.28];
[voxel_sorted, idx_vx] = sort(voxel_size);

cnr2_avg16_masked = cnr2_avg16(mask_idx_array);
figure('Position', [100 100 400 400]); 
plot(voxel_sorted, cnr2_avg16_masked(idx_vx), '.', 'MarkerSize', 42);
hold on;
plot(voxel_sorted, cnr2_cmr(idx_vx), '.', 'MarkerSize', 42);
ylim([0 10]);
grid on;

cnr1_avg16_masked = cnr1_avg16(mask_idx_array);
figure('Position', [100 100 400 400]); 
plot(voxel_sorted, cnr1_avg16_masked(idx_vx), '.', 'MarkerSize', 42);
hold on;
plot(voxel_sorted, cnr1_cmr(idx_vx), '.', 'MarkerSize', 42);
grid on;
ylim([0 10]);