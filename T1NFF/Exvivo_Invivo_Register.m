clear all;
close all;

%% Figure out if ex-vivo to in-vivo is  working
% Read ex-vivo data - T2* map
addpath('../function/');
dicom_dir = uigetdir;
folder_glob = glob(cat(2, dicom_dir, '\*'));

labels = {'T2Star'};
label = labels{1};

idx_array = contains(folder_glob, label);
[list_to_read, order_to_read] = NamePicker(folder_glob(idx_array));

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
    'TriggerTime',...
    'RepetitionTime',...
    };

whatsinit = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit{i} slice_data] = dicom23D(f, dicom_fields);
end
%% Draw Masks
img = whatsinit{1};
name = '18D16';
name = 'Sahara';
name = 'Mojave';
base_dir = uigetdir;
strings = strsplit(list_to_read{1}, '/');
strings = strsplit(strings{end-1}, '_');
num_label = strings{end};

mask_f = cat(2, base_dir, '/data/', name, '/Mask_Exvivo_', num_label, '.mat');
% fh = figure('Position', [100 100 300 400]);
% axis tight manual

mask_epi = zeros(size(img));
mask_endo = zeros(size(img));
mask_roi = zeros(size(img));

if ~exist(mask_f, 'file')
    figure();
    for i = 1:size(img, 3)
        imagesc(img(:,:,i)); axis image; axis off; colormap gray; caxis([0 50]); title(cat(2, 'ROI, slice ', num2str(i)));
        roi = drawpolygon;
        mask_roi(:,:,i) = createMask(roi);
    end
    
    save(mask_f, 'mask_roi');
else
    load(mask_f);
end

mask_f_lv = cat(2, base_dir, '/data/', name, '/Mask_Exvivo_', num_label, '_LV.mat');
if ~exist(mask_f_lv, 'file')
    figure();
    for i = 1:size(img, 3)
        imagesc(img(:,:,i)); axis image; axis off; colormap gray; caxis([0 50]); title(cat(2, 'ROI, slice ', num2str(i)));
        roi = drawpolygon;
        mask_lv(:,:,i) = createMask(roi);
    end
    
    save(mask_f_lv, 'mask_lv');
else
    load(mask_f_lv);
end


mask_f_blood = cat(2, base_dir, '/data/', name, '/Mask_Exvivo_', num_label, '_Blood.mat');
if ~exist(mask_f_blood, 'file')
    figure();
    for i = 1:size(img, 3)
        imagesc(img(:,:,i)); axis image; axis off; colormap gray; caxis([0 50]); title(cat(2, 'ROI, slice ', num2str(i)));
        roi = drawpolygon;
        mask_blood(:,:,i) = createMask(roi);
    end
    
    save(mask_f_blood, 'mask_blood');
else
    load(mask_f_blood);
end

mask_f_insert = cat(2, base_dir, '/data/', name, '/Mask_Exvivo_', num_label, '_InsertionPt.mat');
xs = zeros(1,size(img, 3));
ys = zeros(1,size(img, 3));
xi = zeros(1,size(img, 3));
yi = zeros(1,size(img, 3));

if ~exist(mask_f_insert, 'file')
    insert_pts = struct;
    figure();
    for i = 1:size(img, 3)
        imagesc(img(:,:,i)); axis image; axis off; colormap gray; caxis([0 50]); title(cat(2, 'Superior Insertion, slice ', num2str(i)));
        [xs(i),ys(i)] = getpts;
        imagesc(img(:,:,i)); axis image; axis off; colormap gray; caxis([0 50]); title(cat(2, 'Inferior Insertion, slice ', num2str(i)));
        [xi(i),yi(i)] = getpts;
    end
    insert_pts.xs = xs;
    insert_pts.ys = ys;
    insert_pts.xi = xi;
    insert_pts.yi = yi;
    
    save(mask_f_insert,'-struct', 'insert_pts');
else
    load(mask_f_insert);
end
%% Load invivo data (Skip for Cine)
dicom_dir = uigetdir;
folder_glob = glob(cat(2, dicom_dir, '\*'));

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
    'TriggerTime',...
    'RepetitionTime',...
    };

% img_invivo = zeros(size(V1, 1), size(V1, 2), 4);
% img_invivo(:,:,1) = squeeze(V(:,:,:,4));
% img_invivo(:,:,2) = squeeze(V3(:,:,:,4));
% img_invivo(:,:,3) = squeeze(V2(:,:,:,4));
% img_invivo(:,:,4) = squeeze(V4(:,:,:,4));

whatsinit = cell(length(folder_glob), 1);
for i = 1:length(folder_glob)
    f = folder_glob{i};
    [whatsinit{i} slice_data_invivo] = dicom23D(f, dicom_fields);
end

img_invivo = zeros(size(whatsinit{1},1), size(whatsinit{1},2), length(whatsinit));
for i = 1:length(whatsinit)
    temp = whatsinit{i};
    img_invivo(:,:,i) = temp(:,:,1);
end

%% Draw masks (Skip for Cine)
% name = '18D16';
strings = strsplit(folder_glob{1}, '/');
strings = strsplit(strings{end-1}, '_');
num_label_invivo = strings{end};

mask_roi_invivo = zeros(size(img_invivo));
mask_f_invivo = cat(2, base_dir, '/data/', name, '/Mask_Invivo_', num_label_invivo, '.mat');

if ~exist(mask_f_invivo, 'file')
    for i = 1:size(img_invivo, 3)
        imagesc(img_invivo(:,:,i)); axis image; axis off; colormap gray; title(cat(2, 'ROI, slice ', num2str(i)));
        roi = drawpolygon;
        mask_roi_invivo(:,:,i) = createMask(roi); % mask of
    end
    
    save(mask_f_invivo, 'mask_roi_invivo');
else
    load(mask_f_invivo);
end


figure();
mask_f_invivo_lv = cat(2, base_dir, '/data/', name, '/Mask_Invivo_', num_label_invivo, '_LV.mat');
if ~exist(mask_f_invivo_lv, 'file')
    for i = 1:size(img_invivo, 3)
        imagesc(img_invivo(:,:,i)); axis image; axis off; colormap gray; title(cat(2, 'LV, slice ', num2str(i)));
        roi = drawpolygon;
        mask_lv_invivo(:,:,i) = createMask(roi); % mask of
    end
    
    save(mask_f_invivo_lv, 'mask_lv_invivo');
else
    load(mask_f_invivo_lv);
end

figure();
mask_f_blood_invivo = cat(2, base_dir, '/data/', name, '/Mask_Invivo_', num_label_invivo, '_Blood.mat');
if ~exist(mask_f_blood_invivo, 'file')
    for i = 1:size(img_invivo, 3)
        imagesc(img_invivo(:,:,i)); axis image; axis off; colormap gray; title(cat(2, 'Blood, slice ', num2str(i)));
        roi = drawpolygon;
        mask_blood_invivo(:,:,i) = createMask(roi);
    end
    
    save(mask_f_blood_invivo, 'mask_blood_invivo');
else
    load(mask_f_blood_invivo);
end

mask_f_insert_invivo = cat(2, base_dir, '/data/', name, '/Mask_Invivo_', num_label_invivo, '_InsertionPt.mat');
xs_invivo = zeros(1,size(img_invivo, 3));
ys_invivo = zeros(1,size(img_invivo, 3));
xi_invivo = zeros(1,size(img_invivo, 3));
yi_invivo = zeros(1,size(img_invivo, 3));

if ~exist(mask_f_insert_invivo, 'file')
    insert_pts = struct;
    figure();
    for i = 1:size(img_invivo, 3)
        imagesc(img_invivo(:,:,i)); axis image; axis off; colormap gray; title(cat(2, 'Superior Insertion, slice ', num2str(i)));
        [xs_invivo(i),ys_invivo(i)] = getpts;
        imagesc(img_invivo(:,:,i)); axis image; axis off; colormap gray; title(cat(2, 'Inferior Insertion, slice ', num2str(i)));
        [xi_invivo(i),yi_invivo(i)] = getpts;
    end
    insert_pts.xs_invivo = xs_invivo;
    insert_pts.ys_invivo = ys_invivo;
    insert_pts.xi_invivo = xi_invivo;
    insert_pts.yi_invivo = yi_invivo;
    
    save(mask_f_insert_invivo,'-struct', 'insert_pts');
else
    load(mask_f_insert_invivo);
end

%% Exvivo display
mask_roi(mask_roi == 0) = nan;
mask_lv = double(mask_lv);
mask_lv(mask_lv == 0) = nan;
mask_roi_new = zeros(size(mask_roi));
mask_lv_new = zeros(size(mask_roi));
mask_lv_new = zeros(size(mask_roi));
mask_blood_new = zeros(size(mask_roi));
xs_new = zeros(size(xs));
ys_new = zeros(size(xs));
xi_new = zeros(size(xs));
yi_new = zeros(size(xs));

img_new = zeros(size(mask_roi));

for i = 1:size(mask_roi, 3)
    mask_roi_new(:,:,i) =  mask_roi(:,:,end-i+1);
    mask_lv_new(:,:,i) =  mask_lv(:,:,end-i+1);
    img_new(:,:,i) = img(:,:,end-i+1);
    mask_blood_new(:,:,i) = mask_blood(:,:,end-i+1);
    xs_new(i) = xs(end-i+1);
    ys_new(i) = ys(end-i+1);
    xi_new(i) = xi(end-i+1);
    yi_new(i) = yi(end-i+1);
end

D = double(squeeze(mask_roi_new.*img_new));

figure();
h = slice(D, [], [], 1:size(D,3)); caxis([0 50]);
set(h, 'EdgeColor','none', 'FaceColor','interp');
alpha(.3);
%light('position',[0,0,20],'style','local','color','w');lighting phong;

%% Invivo display (Skip for Cine)
sth = slice_data_invivo(1).SliceThickness;
mask_roi_invivo_interp = zeros(size(mask_roi_invivo,1), size(mask_roi_invivo,2), (size(mask_roi_invivo,3)-1)*sth+1);
img_invivo_interp = zeros(size(mask_roi_invivo,1), size(mask_roi_invivo,2), (size(mask_roi_invivo,3)-1)*sth+1);
img_invivo_gap = zeros(size(mask_roi_invivo,1), size(mask_roi_invivo,2), (size(mask_roi_invivo,3)-1)*sth+1);
mask_roi_invivo_gap = zeros(size(mask_roi_invivo,1), size(mask_roi_invivo,2), (size(mask_roi_invivo,3)-1)*sth+1);


for i = 1:size(img_invivo_interp,3)
    n = ceil(i/sth);
    mod_i = mod(i-1, sth);
   if  mod_i == 0
       mask_roi_invivo_interp(:,:,mod_i+sth*(n-1)+1) = nan;
       img_invivo_interp(:,:,mod_i+sth*(n-1)+1) = nan;
       
       mask_roi_invivo_gap(:,:,mod_i+sth*(n-1)+1) = mask_roi_invivo(:,:,n);
       img_invivo_gap(:,:,mod_i+sth*(n-1)+1) = img_invivo(:,:,n);
   else
       mask_roi_invivo_interp(:,:,mod_i+sth*(n-1)+1) = mask_roi_invivo(:,:,n) * (sth-mod_i)/sth + mask_roi_invivo(:,:,n+1) * mod_i/sth;
       img_invivo_interp(:,:,mod_i+sth*(n-1)+1) = img_invivo(:,:,n) * (sth-mod_i)/sth + img_invivo(:,:,n+1) * mod_i/sth;
   end
end

mask_roi_invivo(mask_roi_invivo == 0) = nan;
mask_lv_invivo = double(mask_lv_invivo);
mask_lv_invivo(mask_lv_invivo == 0) = nan;
mask_roi_new_invivo = zeros(size(mask_roi_invivo));
mask_lv_new_invivo = zeros(size(mask_lv_invivo));
mask_blood_new_invivo = zeros(size(mask_blood_invivo));
img_invivo_new = zeros(size(mask_roi_invivo));

mask_roi_invivo_interp(mask_roi_invivo_interp == 0) = nan;
mask_roi_new_invivo_interp = zeros(size(mask_roi_invivo_interp));
img_invivo_new_interp = zeros(size(mask_roi_invivo_interp));
mask_roi_new_invivo_draw = zeros(size(mask_roi_invivo_interp));


mask_roi_invivo_gap(mask_roi_invivo_gap == 0) = nan;
mask_roi_new_invivo_gap = zeros(size(mask_roi_invivo_gap));
img_invivo_new_gap = zeros(size(mask_roi_invivo_gap));

for i = 1:size(mask_roi_invivo, 3)
    mask_roi_new_invivo(:,:,i) =  mask_roi_invivo(:,:,end-i+1);
    mask_lv_new_invivo(:,:,i) =  mask_lv_invivo(:,:,end-i+1);
    mask_blood_new_invivo(:,:,i) = mask_blood_invivo(:,:,end-i+1);
    img_invivo_new(:,:,i) = img_invivo(:,:,end-i+1);
end

for i = 1:size(mask_roi_invivo_interp, 3)
    mask_roi_new_invivo_interp(:,:,i) =  mask_roi_invivo_interp(:,:,end-i+1);
    img_invivo_new_interp(:,:,i) = img_invivo_interp(:,:,end-i+1);
end

for i = 1:size(mask_roi_invivo_gap, 3)
    mask_roi_new_invivo_gap(:,:,i) =  mask_roi_invivo_gap(:,:,end-i+1);
    img_invivo_new_gap(:,:,i) = img_invivo_gap(:,:,end-i+1);
end

mask_f_roi_invivo_draw = cat(2, base_dir, '/data/', name, '/Mask_Invivo_', num_label_invivo, '_ForInterp.mat');

if ~exist(mask_f_roi_invivo_draw, 'file')
    figure();
    for i = 1:size(img_invivo_new_interp, 3)
        imagesc(img_invivo_new_interp(:,:,i)); axis image; axis off; colormap gray; title(cat(2, 'ROI, slice ', num2str(i)));
        roi = drawpolygon;
        mask_roi_new_invivo_draw(:,:,i) = createMask(roi);
    end
    
    save(mask_f_roi_invivo_draw, 'mask_roi_new_invivo_draw');
else
    load(mask_f_roi_invivo_draw);
end

mask_roi_new_invivo_draw(mask_roi_new_invivo_draw == 0) = nan;

D2 = double(squeeze(mask_roi_new_invivo_gap.*img_invivo_new_gap));
figure();
h1 = slice(D2, [], [], 1:size(D2,3)); %caxis([0 50]);
set(h1, 'EdgeColor','none', 'FaceColor','interp');
alpha(h1, .4);

%figure();
hold on;
D3 = double(squeeze(mask_roi_new_invivo_interp.*img_invivo_new_interp));
h2 = slice(D3, [], [], 1:size(D3,3)); %caxis([0 50]);
set(h2, 'EdgeColor','none', 'FaceColor','interp');
alpha(h2, .2);


D2 = double(squeeze(mask_roi_new_invivo_gap.*img_invivo_new_gap));
figure();
h1 = slice(D2, [], [], 1:size(D2,3)); %caxis([0 50]);
set(h1, 'EdgeColor','none', 'FaceColor','interp');
alpha(h1, .4);

%figure();
hold on;
D3 = double(squeeze(mask_roi_new_invivo_draw.*img_invivo_new_interp));
h2 = slice(D3, [], [], 1:size(D3,3)); %caxis([0 50]);
set(h2, 'EdgeColor','none', 'FaceColor','interp');
alpha(h2, .2);
%% Try to find a way to register between Invivo and Exvivo (Skip for Cine)
res = slice_data(1).PixelSpacing;
thickness = slice_data(1).SliceThickness;
res_invivo = slice_data_invivo(1).PixelSpacing;
thickness_invivo = slice_data_invivo(1).SliceThickness;

mm2_exvivo = res(1) * res(2);
mm2_invivo = res_invivo(1) * res_invivo(2);

nanIdx = isnan(mask_lv_new) | isnan(mask_roi_new); 
mask_lv_new(nanIdx) = 0;
mask_roi_new(isnan(mask_roi_new)) = 0;
nanIdx_invivo = isnan(mask_lv_new_invivo) | isnan(mask_roi_new_invivo); 
mask_lv_new_invivo(nanIdx_invivo) = 0;
mask_roi_new_invivo(isnan(mask_roi_new_invivo)) = 0;

mask_lv_union = mask_lv_new & mask_roi_new;
mask_lv_invivo_union = mask_lv_new_invivo & mask_roi_new_invivo;

mask_lv_union_reshape = reshape(mask_lv_union, [], size(mask_lv_union, 3));
mask_lv_invivo_union_reshape = reshape(mask_lv_invivo_union, [], size(mask_lv_invivo_union, 3));

mask_roi_new_reshape = reshape(mask_roi_new, [], size(mask_roi_new ,3)); % ex-vivo
mask_roi_new_invivo_reshape = reshape(mask_roi_new_invivo, [], size(mask_roi_new_invivo ,3)); % in-vivo

slice_area_exvivo = mm2_exvivo * sum(mask_roi_new_reshape, 'omitnan')/100;
slice_area_invivo = mm2_invivo * sum(mask_roi_new_invivo_reshape, 'omitnan')/100;
lv_area_exvivo = mm2_exvivo * sum(mask_lv_union_reshape, 'omitnan')/100;
lv_area_invivo = mm2_invivo * sum(mask_lv_invivo_union_reshape, 'omitnan')/100;

idx = lv_area_exvivo & slice_area_exvivo;
slice_area_exvivo_nonzero = slice_area_exvivo;
slice_area_exvivo_nonzero(slice_area_exvivo == 0) = nan;

interp_slice_invivo = ((length(slice_area_invivo)-1)*thickness_invivo+1);
slice_area_invivo_interp = interp1(1:thickness_invivo:interp_slice_invivo,  slice_area_invivo, 1:interp_slice_invivo);
figure();
plot(slice_area_exvivo_nonzero, 'LineWidth', 2); hold on;
plot(slice_area_invivo_interp, 'LineWidth', 2);


figure();
%plot(slice_area_exvivo); hold on; 
plot(diff(slice_area_exvivo)); ylim([-1 3]);
figure();
plot(diff(slice_area_invivo));

lv_perc_invivo = lv_area_invivo ./ slice_area_invivo;
lv_perc_invivo_interp = interp1(1:thickness_invivo:interp_slice_invivo, lv_perc_invivo, 1:interp_slice_invivo);
lv_perc_exvivo = lv_area_exvivo(idx) ./ slice_area_exvivo(idx);
figure();
plot(lv_perc_exvivo, 'LineWidth', 1.5); hold on; plot(lv_perc_invivo_interp, 'LineWidth', 1.5);
yline(max(lv_area_invivo ./ slice_area_invivo));

n = length(lv_perc_exvivo)-length(lv_perc_invivo_interp)+1;
RMS = zeros(1, n);

for i = 1:n
    RMS(i) = sum((lv_perc_exvivo(i:(length(lv_perc_invivo_interp)+i-1)) - lv_perc_invivo_interp).^2);
end

[a,b] = min(RMS);

x = b:(b+length(lv_perc_invivo_interp)-1);
figure();
plot(lv_perc_exvivo, 'LineWidth', 1.5); hold on; plot(x, lv_perc_invivo_interp, 'LineWidth', 1.5);
yline(max(lv_area_invivo ./ slice_area_invivo));
%% Find Superior and Inferior insertion point
%% Centerline (Skip for Cine)
[sep_ratio_invivo, centroids_invivo, mask_centerline_invivo] = ...
    Func_septum_ratio(mask_lv_new_invivo, mask_blood_new_invivo, xs_invivo, ys_invivo, xi_invivo, yi_invivo);

[sep_ratio_exvivo, centroids_exvivo, mask_centerline_exvivo] = ...
    Func_septum_ratio(mask_lv_new, mask_blood_new, xs, ys, xi, yi);

figure(); plot(sep_ratio_invivo, 'LineWidth', 1.5); hold on; plot(sep_ratio_exvivo, 'LineWidth', 1.5);

%% Load CINE images
addpath('../function/');
dicom_dir = uigetdir;
folder_glob = glob(cat(2, dicom_dir, '\*'));

labels = {'series'};
label = labels{1};

idx_array = contains(folder_glob, label);
% [list_to_read, order_to_read] = NamePicker(folder_glob(idx_array));
% [243,1,13,24,35,46,57,68,79,90,101] % 18D16
% [7, 9, 11, 13, 15, 17, 19, 21, 23] % Sahara
% [1, 2, 4, 5, 6] % Mojave
list_to_read = folder_glob([1, 2, 4, 5, 6]);
order_to_read = 1:length(list_to_read);

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
    'TriggerTime',...
    'RepetitionTime',...
    };

whatsinit = cell(length(list_to_read), 1);
slice_loc_cine = zeros(length(list_to_read), 1);;
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit{i} slice_data_cine] = dicom23D(f, dicom_fields);
    slice_loc_cine(i) = slice_data_cine(1).SliceLocation;
end

%% ES = 11
es = 11; % both 18D16 and Sahara
img_cine = zeros(size(whatsinit{1}, 1), size(whatsinit{1}, 2), length(whatsinit));
for i = 1:length(whatsinit)
   img_cine(:,:,i) = whatsinit{i}(:,:,es);
end

% figure();
% for i = 1:length(whatsinit)
%     subplot(3,4,i);
%     imagesc(img_cine(:,:,i));
% end

% name = '18D16';
strings = strsplit(folder_glob{1}, '/');
strings = strsplit(strings{end-1}, '_');
num_label_cine = strings{end};

mask_roi_cine = zeros(size(img_cine));
mask_f_cine = cat(2, base_dir, '/data/', name, '/Mask_Cine_', num_label_cine, '.mat');

if ~exist(mask_f_cine, 'file')
    figure();
    for i = 1:size(img_cine, 3)
        imagesc(img_cine(:,:,i)); axis image; axis off; colormap gray; title(cat(2, 'ROI, slice ', num2str(i)));
        roi = drawpolygon;
        mask_roi_cine(:,:,i) = createMask(roi); % mask of
    end
    
    save(mask_f_cine, 'mask_roi_cine');
else
    load(mask_f_cine);
end


mask_f_cine_lv = cat(2, base_dir, '/data/', name, '/Mask_Cine_', num_label_cine, '_LV.mat');
if ~exist(mask_f_cine_lv, 'file')
    figure();
    for i = 1:size(img_cine, 3)
        imagesc(img_cine(:,:,i)); axis image; axis off; colormap gray; title(cat(2, 'LV, slice ', num2str(i)));
        roi = drawpolygon;
        mask_lv_cine(:,:,i) = createMask(roi); % mask of
    end
    
    save(mask_f_cine_lv, 'mask_lv_cine');
else
    load(mask_f_cine_lv);
end


mask_f_blood_cine = cat(2, base_dir, '/data/', name, '/Mask_Cine_', num_label_cine, '_Blood.mat');
if ~exist(mask_f_blood_cine, 'file')
    figure();
    for i = 1:size(img_cine, 3)
        imagesc(img_cine(:,:,i)); axis image; axis off; colormap gray; title(cat(2, 'Blood, slice ', num2str(i)));
        roi = drawpolygon;
        mask_blood_cine(:,:,i) = createMask(roi);
    end
    
    save(mask_f_blood_cine, 'mask_blood_cine');
else
    load(mask_f_blood_cine);
end

mask_f_insert_cine = cat(2, base_dir, '/data/', name, '/Mask_Cine_', num_label_cine, '_InsertionPt.mat');
xs_cine = zeros(1,size(img_cine, 3));
ys_cine = zeros(1,size(img_cine, 3));
xi_cine = zeros(1,size(img_cine, 3));
yi_cine = zeros(1,size(img_cine, 3));

if ~exist(mask_f_insert_cine, 'file')
    insert_pts = struct;
    figure();
    for i = 1:size(img_cine, 3)
        imagesc(img_cine(:,:,i)); axis image; axis off; colormap gray; title(cat(2, 'Superior Insertion, slice ', num2str(i)));
        [xs_cine(i),ys_cine(i)] = getpts;
        imagesc(img_cine(:,:,i)); axis image; axis off; colormap gray; title(cat(2, 'Inferior Insertion, slice ', num2str(i)));
        [xi_cine(i),yi_cine(i)] = getpts;
    end
    insert_pts.xs_cine = xs_cine;
    insert_pts.ys_cine = ys_cine;
    insert_pts.xi_cine = xi_cine;
    insert_pts.yi_cine = yi_cine;
    
    save(mask_f_insert_cine,'-struct', 'insert_pts');
else
    load(mask_f_insert_cine);
end

%% CINE display
mask_roi_cine(mask_roi_cine == 0) = nan;
mask_lv_cine = double(mask_lv_cine);
mask_lv_cine(mask_lv_cine == 0) = nan;
mask_roi_new_cine = zeros(size(mask_roi_cine));
mask_lv_new_cine = zeros(size(mask_roi_cine));
mask_lv_new_cine = zeros(size(mask_roi_cine));
mask_blood_new_cine = zeros(size(mask_roi_cine));
xs_cine_new = zeros(size(xs_cine));
ys_cine_new = zeros(size(xs_cine));
xi_cine_new = zeros(size(xs_cine));
yi_cine_new = zeros(size(xs_cine));
slice_loc_cine_new = zeros(size(xs_cine,1), 1);

img_cine_new = zeros(size(mask_roi_cine));

for i = 1:size(mask_roi_cine, 3)
    mask_roi_new_cine(:,:,i) =  mask_roi_cine(:,:,end-i+1);
    mask_lv_new_cine(:,:,i) =  mask_lv_cine(:,:,end-i+1);
    img_cine_new(:,:,i) = img_cine(:,:,end-i+1);
    mask_blood_new_cine(:,:,i) = mask_blood_cine(:,:,end-i+1);
    xs_cine_new(i) = xs_cine(end-i+1);
    ys_cine_new(i) = ys_cine(end-i+1);
    xi_cine_new(i) = xi_cine(end-i+1);
    yi_cine_new(i) = yi_cine(end-i+1);
    slice_loc_cine_new(i) = slice_loc_cine(end-i+1);
end

D = double(squeeze(mask_roi_new_cine.*img_cine_new));

figure();
h = slice(D, [], [], 1:size(D,3)); % caxis([0 50]);
set(h, 'EdgeColor','none', 'FaceColor','interp');
alpha(.3);

%% Centerline and angles, find septum ratio
[sep_ratio_cine, centroids_cine, mask_centerline_cine, thera_i_array_cine, thera_s_array_cine] = ...
    Func_septum_ratio(mask_lv_new_cine, mask_blood_new_cine, xs_cine_new, ys_cine_new, xi_cine_new, yi_cine_new);

[sep_ratio_exvivo, centroids_exvivo, mask_centerline_exvivo, thera_i_array_exvivo, thera_s_array_exvivo] = ...
    Func_septum_ratio(mask_lv_new, mask_blood_new, xs_new, ys_new, xi_new, yi_new);

% sep_ratio_cine(sep_ratio_cine == 0) = 1;
% sep_ratio_cine(end-3) = []; % 18D16
% thera_s_array_cine(end-3) = [];
% thera_i_array_cine(end-3) = [];
% centroids_cine{end-3} = {};
% slice_loc_cine_new_new = slice_loc_cine_new;
% slice_loc_cine_new_new(end-3) = []; % for 18D16
% slice_loc_cine_new_new(thera_s_array_cine == 0) = [];
% thera_s_array_cine(thera_s_array_cine == 0) = [];
% thera_i_array_cine(thera_i_array_cine == 0) = [];
centroids_cine_new = centroids_cine(~cellfun(@isempty, centroids_cine));

slice_loc_cine_new_new = slice_loc_cine_new;
slice_loc_cine_new_new(thera_s_array_cine == 0) = [];

sep_ratio_exvivo_med = medfilt1(sep_ratio_exvivo, 3);
sep_ratio_exvivo_med(sep_ratio_exvivo == 0) = 0;
figure(); plot(sep_ratio_cine, 'LineWidth', 1.5); hold on; plot(sep_ratio_exvivo_med, 'LineWidth', 1.5);

%% Save as gif
mask_centerline_cine_nan = mask_centerline_cine;
mask_centerline_cine_nan(mask_centerline_cine_nan == 0) = nan; 
mask_centerline_cine_nan(mask_centerline_cine_nan >= 1) = 1;

cine_gif = cat(2, base_dir, '/data/', name, '/Cine_', num_label_cine, '.gif');

figure(101);
for slc = 1:size(img_cine_new,3)
    ax1 = axes;
    imagesc(img_cine_new(:,:,slc)); pbaspect([size(img_cine_new, 2), size(img_cine_new, 1) 1]);
    colormap(ax1, 'gray');
    set(ax1, 'xticklabel', []); set(ax1, 'yticklabel', []);
    ax2 = axes;
    
    imagesc(ax2, mask_centerline_cine_nan(:,:,slc) .* mask_centerline_cine(:,:,slc), 'AlphaData', mask_centerline_cine(:,:,slc));
    pbaspect([size(img_cine_new, 2), size(img_cine_new, 1) 1]); colormap(ax2, 'cool');
    
    %caxis(ax1, [0 100]); caxis(ax2, [-2 10]); 
    linkprop([ax1 ax2], 'Position');
    ax2.Visible = 'off';
    % pause(1);
    drawnow
    frame = getframe(101);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if slc == 1;
        imwrite(imind,cm,cine_gif,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,cine_gif,'gif','WriteMode','append');
    end
end


mask_centerline_exvivo_nan = mask_centerline_exvivo;
mask_centerline_exvivo_nan(mask_centerline_exvivo_nan == 0) = nan; 
mask_centerline_exvivo_nan(mask_centerline_exvivo_nan >= 1) = 1;

exvivo_gif = cat(2, base_dir, '/data/', name, '/Exvivo_', num_label, '.gif');

figure(102);
for slc = 1:size(img_new,3)
    ax1 = axes;
    imagesc(img_new(:,:,slc)); pbaspect([size(img_new, 2), size(img_new, 1) 1]);
    colormap(ax1, 'gray');
    set(ax1, 'xticklabel', []); set(ax1, 'yticklabel', []);
    ax2 = axes;
    
    imagesc(ax2, mask_centerline_exvivo_nan(:,:,slc) .* mask_centerline_exvivo(:,:,slc), 'AlphaData', mask_centerline_exvivo(:,:,slc));
    pbaspect([size(img_new, 2), size(img_new, 1) 1]); colormap(ax2, 'cool');
    
    caxis(ax1, [0 100]); %caxis(ax2, [-2 10]); 
    linkprop([ax1 ax2], 'Position');
    ax2.Visible = 'off';
    % pause(1);
    drawnow
    frame = getframe(102);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if slc == 1;
        imwrite(imind,cm,exvivo_gif,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,exvivo_gif,'gif','WriteMode','append');
    end
end

%% Match slices CINE vs gif
[exvivo_match, cine_match] = matching_slices(slice_data_cine, sep_ratio_cine, sep_ratio_exvivo_med)

img_for_rotate = img_new(:,:,exvivo_match);
mask_centerline_exvivo_nan_for_rotate = mask_centerline_exvivo_nan(:,:,exvivo_match);
mask_centerline_exvivo_for_rotate = mask_centerline_exvivo(:,:,exvivo_match);

% 18D16
% I2 = flipdim(img_for_rotate ,2);
% I3 = flipdim(mask_centerline_exvivo_nan_for_rotate,2);
% I4 = flipdim(mask_centerline_exvivo_for_rotate,2);

% Sahara
I2 = flipdim(img_for_rotate,1);
I3 = flipdim(mask_centerline_exvivo_nan_for_rotate,1);
I4 = flipdim(mask_centerline_exvivo_for_rotate,1);

% deg_rot = -(rad2deg(thera_i_array_exvivo(exvivo_match))+rad2deg(thera_s_array_cine)); % Because messed up with inferior and superior in 18D16 this line needs to be worked
deg_rot = -(rad2deg(thera_s_array_exvivo(exvivo_match(5)))+rad2deg(thera_s_array_cine(5)));

% img_for_rotate_new = cell(length(deg_rot), 1);
% I3_new = cell(length(deg_rot), 1);
% I4_new = cell(length(deg_rot), 1);
% 
% for i = 1:length(deg_rot)
%     img_for_rotate_new{i} = imrotate(I2(:,:,i), deg_rot(i));
%     I3_new{i} = imrotate(I3(:,:,i), deg_rot(i));
%     I4_new{i} = imrotate(I4(:,:,i), deg_rot(i));
% end

img_for_rotate_new = imrotate(I2, deg_rot);
I3_new = imrotate(I3, deg_rot);
I4_new = imrotate(I4, deg_rot);

exvivo_gif = cat(2, base_dir, '/data/', name, '/Exvivo_rotate', num_label_cine, '.gif');

figure(103);
for slc = 1:size(img_for_rotate_new, 3)
    ax1 = axes;
    imagesc(img_for_rotate_new(:,:,slc)); pbaspect([size(img_for_rotate_new(:,:,slc), 2), size(img_for_rotate_new(:,:,slc), 1) 1]);
    colormap(ax1, 'gray');
    set(ax1, 'xticklabel', []); set(ax1, 'yticklabel', []);
    ax2 = axes;
    
    imagesc(ax2, I3_new(:,:,slc) .* I4_new(:,:,slc), 'AlphaData', I4_new(:,:,slc));
    pbaspect([size(img_for_rotate_new(:,:,slc), 2), size(img_for_rotate_new(:,:,slc), 1) 1]); colormap(ax2, 'cool');
    
    caxis(ax1, [0 100]); %caxis(ax2, [-2 10]); 
    linkprop([ax1 ax2], 'Position');
    ax2.Visible = 'off';
    % pause(1);
    drawnow
    frame = getframe(103);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if slc == 1;
        imwrite(imind,cm,exvivo_gif,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,exvivo_gif,'gif','WriteMode','append');
    end
end

%% So what do we need for the invivo-exvivo matching
% insertion points
% centroids
% slice location of invivo and exvivo
slice_loc_cine_new_new = slice_loc_cine_new_new(cine_match);

reg_info = struct;
% reg_info.exvivo_SliceLoc = slice_loc_cine;
reg_info.cine_SliceLoc = slice_loc_cine_new_new;
reg_info.thera_i_array_cine = thera_i_array_cine;
reg_info.thera_s_array_cine = thera_s_array_cine;
reg_info.thera_i_array_exvivo = thera_i_array_exvivo;
reg_info.thera_s_array_exvivo = thera_s_array_exvivo;
reg_info.centroids_cine_new = centroids_cine_new;
reg_info.centroids_exvivo = centroids_exvivo;
reg_info.exvivo_match = exvivo_match;
reg_info_f = cat(2, base_dir, '/data/', name, '/RegInfo.mat');
save(reg_info_f, '-struct', 'reg_info');