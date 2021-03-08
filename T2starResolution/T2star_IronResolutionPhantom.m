close all;
clear all;
addpath('../function/');
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

% disp('Avg 0016 starts here: ');
% avg_num = input('Please type average number here:  ');
% avg_name = cat(2, 'Avg', num2str(avg_num, '%04.f'));
%% Read T2* DICOM files
whatsinit = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    whatsinit{i} = dicom23D(f);
end

% Display images
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    subplot(2,9,i);
    imagesc(whatsinit{i}.*0.1); axis image;
    caxis([0 100])
end
%% Draw ROIs in N vials (This can be skipped)
img = whatsinit{1}.*0.1;

roi_save = cat(2, subject_data_dir, 'roi.mat');
dim = input('Dimension of vials (1 or 2): ');

if ~exist(roi_save, 'file')
    if dim == 1
        roi_coords_cell = cell(N, 1);
        N = input('Number of vials: ');
        roi = cell(N, 1);
        for i = 1:length(roi)
            figure('Position', [100 0 1600 1600]); imagesc(img); axis image; caxis([0 100]); colorbar;
            temp = drawpolygon(gca); % You may need to use roipoly and its equivalent if MATLAB version is earlier than 2020a
            roi_coords_cell{i} = temp.Position;
            roi{i} = createMask(temp);
        end
    elseif dim == 2
        row = input('Number of rows: ');
        col = input('Number of cols: ');
        roi_coords_cell = cell(row, col, 1);
        N = row * col;
        roi = cell(row, col);
        for k = 1:row
            for j = 1:col
                figure('Position', [100 0 1600 1600]); imagesc(img); axis image; caxis([0 100]); colorbar;
                temp = drawpolygon(gca); % You may need to use roipoly and its equivalent if MATLAB version is earlier than 2020a
                roi_coords_cell{i, j} = temp.Position;
                roi{k, j} = createMask(temp);
            end
        end
    end
else
    load(roi_save);
    row = size(roi, 1);
    col = size(roi, 2);
    N = row * col;
end
%% convert to masks
img_size_truth = size(img);
mask_save = cat(2, subject_data_dir, 'mask.mat');
roi_cell = cell(N, length(whatsinit));

if ~exist(mask_save, 'file')
    figure('Position', [100 0 1600 1600]);
    mask_struct = struct;
    for i = 1:length(whatsinit)
        roi_mask = cell(N, 1);
        img2 = whatsinit{i} .* 0.1;
        img2_size = size(whatsinit{i});
        ratio = img_size_truth(1) ./ img2_size(1);
        imagesc(img2); caxis([0 100]);
        for n = 1:N
            temp = drawpolygon(gca,'Position', [roi_coords_cell{n}(:,1)/ratio + (ratio-1)/ratio, roi_coords_cell{n}(:,2)/ratio + (ratio-1)/ratio]);
            roi_mask{n} = createMask(temp);
        end
        
        mask_struct(i).roi_mask = roi_mask;
    end
    
    save(mask_save, 'mask_struct');
else
    load(mask_save);
end

%% Getting mean values of T2*
t2star_mean = zeros(N, length(whatsinit));
t2star_sd = zeros(N, length(whatsinit));

for i = 1:length(whatsinit)
    img2 = whatsinit{i} .* 0.1;
    for n = 1:N
        t2star_mean(n, i) = mean(nonzeros(img2 .* mask_struct(i).roi_mask{n}));
        t2star_sd(n, i) = std(nonzeros(img2 .* mask_struct(i).roi_mask{n}));
    end
end

d = length(whatsinit)/2;
t2star_mean_reshape = reshape(t2star_mean, N, d, length(whatsinit)/d);
t2star_sd_reshape = reshape(t2star_sd, N, d, length(whatsinit)/d);
x = [0, d, 2*d] + [0 , 0.5, 0.5];
inplane_res = 1:d;
res = [inplane_res; inplane_res + d];

figure('Position', [100 0 1600 1600]);
fe_conc = {'0', '10', '20', '30', '40', '50'};
for n = 1:N
    subplot(2,3,n);
    errorbar(t2star_mean(n, :), t2star_sd(n, :), 'LineStyle', 'none' );
    n_min = min(t2star_mean(n, :));
    ylim([n_min-10, n_min+35]);
    ylim_lb = min(ylim); ylim_ub = max(ylim);
    patch([x(1) x(2) x(2) x(1)], [ylim_ub ylim_ub 0 0], [241 194 151]/255, 'FaceAlpha',.5)
    patch([x(2) x(3) x(3) x(2)], [ylim_ub ylim_ub 0 0], [199 213 161]/255, 'FaceAlpha',.5)
    errorbar(res', squeeze(t2star_mean_reshape(n, :,:)), squeeze(t2star_sd_reshape(n, :,:)), 'LineWidth', 2);
    
    grid on;
    ylim([ylim_lb, ylim_ub]);
    %ylim([0, 150]);
    title(['Fe: ', fe_conc{n}, ' ug/L'])
end

%% Read T2*-weighted DICOM files
[list_to_read, order_to_read] = NamePicker(folder_glob);
whatsinit = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    whatsinit{i} = dicom23D(f);
end

% Display images
figure('Position', [100 0 1600 1600]);
for i = 1:length(whatsinit)
    subplot(2,9,i);
    imagesc(whatsinit{i}(:,:,1)); axis image;
    %caxis([0 100])
end


img_size_truth = size(img);
mask_save = cat(2, subject_data_dir, 'mask.mat');
roi_cell = cell(N, length(whatsinit));

if ~exist(mask_save, 'file')
    figure('Position', [100 0 1600 1600]);
    mask_struct = struct;
    for i = 1:length(whatsinit)
        roi_mask = cell(N, 1);
        img2 = whatsinit{i} .* 0.1;
        img2_size = size(whatsinit{i});
        ratio = img_size_truth(1) ./ img2_size(1);
        imagesc(img2); caxis([0 100]);
        for n = 1:N
            temp = drawpolygon(gca,'Position', [roi_coords_cell{n}(:,1)/ratio + (ratio-1)/ratio, roi_coords_cell{n}(:,2)/ratio + (ratio-1)/ratio]);
            roi_mask{n} = createMask(temp);
        end
        
        mask_struct(i).roi_mask = roi_mask;
    end
    
    save(mask_save, 'mask_struct');
else
    load(mask_save);
end

%% Getting mean values of T2*
t2star_w_mean = zeros(N, length(whatsinit), size(whatsinit{1}, 3));
t2star_w_sd = zeros(N, length(whatsinit), size(whatsinit{1}, 3));

for i = 1:length(whatsinit)
    img2 = whatsinit{i};
    for n = 1:N
        for t = 1:size(whatsinit{1}, 3)
            t2star_w_mean(n, i, t) = mean(nonzeros(img2(:,:,t) .* mask_struct(i).roi_mask{n}));
            t2star_w_sd(n, i, t) = std(nonzeros(img2(:,:,t) .* mask_struct(i).roi_mask{n}));
        end
    end
end

d = length(whatsinit)/2;
ax = cell(N, 1);
figure('Position', [100 0 1600 1600]);
fe_conc = {'0', '10', '20', '30', '40', '50'};
% Initialize video
subject_img_dir = GetFullPath(cat(2, save_dir, subject_name, '/'));

vid_save = cat(2, subject_img_dir, 'myVideoFile.mov');
myVideo = VideoWriter(vid_save); %open video file
myVideo.FrameRate = 2;  %can adjust this, 5 - 10 works well for me
open(myVideo)

for i = 1:length(whatsinit)/2
    for n = 1:N
        
        ax{n} = subplot(2,3,n);
        
        if i == 1
            hold on;
            grid on;
            ylim([0, 800]);
            title(['Fe: ', fe_conc{n}, ' ug/L']);
            pause(.5);
        end
        axes(ax{n})
        plot(squeeze(t2star_mean(n, i, :)), 'LineWidth', 2);
        
    end
    pause(.5);
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
end
close(myVideo)