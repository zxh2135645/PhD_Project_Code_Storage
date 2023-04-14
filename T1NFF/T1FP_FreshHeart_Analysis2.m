clear all;
close all;

addpath('../function/')
dicom_dir = uigetdir;
ff_glob = glob(GetFullPath(cat(2, dicom_dir, '/../Result/AllPhasemap_*_SinglePeak.mat')));
quant_map = load(GetFullPath(cat(2, dicom_dir, '/../Result/QuantitativeMap.mat')));
% ff_glob_chris = glob(GetFullPath(cat(2, dicom_dir, '/../Result/Jesse_IDEAL_4ave_exvivo_0902.mat')));
ff_map = load(ff_glob{1});

fat = ff_map.fat;
water = ff_map.water;
R2s = ff_map.R2s;
fat = mean(fat, 3);
water = mean(water, 3);
R2s = mean(R2s, 3);

ff = zeros(size(fat));
fat_flag = fat > water;
ff(fat_flag) = abs(fat(fat_flag)) ./ abs(fat(fat_flag) + water(fat_flag));
ff(~fat_flag) = 1 - (abs(water(~fat_flag))) ./ abs(fat(~fat_flag) + water(~fat_flag));

t1_map = quant_map.t1_map;
t2_map = quant_map.t2_map;
t2star_map = quant_map.t2star_map(:,:,4);

% ff_single_slice = mean(ff, 3);
R2s_flip = flip(imrotate(R2s,-90),2);
ff_flip = flip(imrotate(ff,-90),2);
figure(); imagesc(R2s_flip); caxis([0 100]);
figure(); imagesc(ff_flip); caxis([0 0.2]);


figure(); imagesc(t1_map); caxis([0 2000]);
figure(); imagesc(t2_map); caxis([0 100]);
figure(); imagesc(t2star_map); caxis([0,100]);
%% Need to skip right now, the matrix size is not matched
ff_map_ch = load(ff_glob_chris{1});
fat = ff_map_ch.wfat;
water = ff_map_ch.wwater;
%R2s = ff_map_ch.R2s;

ff = zeros(size(fat));
fat_flag = fat > water;
ff(fat_flag) = abs(fat(fat_flag)) ./ abs(fat(fat_flag) + water(fat_flag));
ff(~fat_flag) = 1 - (abs(water(~fat_flag))) ./ abs(fat(~fat_flag) + water(~fat_flag));
figure(); imagesc(ff(:,:,10)); caxis([0,0.2]);
% figure(); imagesc(R2s(:,:,5)); caxis([0 100]);
%% Draw MI and Remote
mask_f = GetFullPath(cat(2, dicom_dir, '/../mask_mi.mat'));
mi = zeros(size(R2s_flip,1), size(R2s_flip,2), size(R2s_flip,3));
remote = zeros(size(R2s_flip,1), size(R2s_flip,2), size(R2s_flip,3));
epi = zeros(size(R2s_flip,1), size(R2s_flip,2), size(R2s_flip,3));
endo = zeros(size(R2s_flip,1), size(R2s_flip,2), size(R2s_flip,3));

if ~exist(mask_f)
for i = 1:size(R2s_flip,3)
    figure();
    imagesc(R2s_flip); axis image; caxis([0 100]);
    roi = drawpolygon;
    mi(:,:,i) = createMask(roi);
    roi = drawpolygon;
    remote(:,:,i) = createMask(roi);
    roi = drawpolygon;
    epi(:,:,i) = createMask(roi);
    roi = drawpolygon;
    endo(:,:,i) = createMask(roi);
end
    mask.mi = mi;
    mask.remote = remote;
    mask.epi = epi;
    mask.endo = endo;
    save(mask_f, '-struct', 'mask');
else
    load(mask_f);
end

centroids = cell(size(epi,3) ,1);
figure();
for i = 1:length(centroids)
    imagesc(epi(:,:,i))
    s = regionprops(epi(:,:,i),'centroid');
    hold on;
    plot(s.Centroid(1), s.Centroid(2), 'r*');
    hold off;
    
    centroids{i} = round(s.Centroid);
end
%% Plot
%% Need to crop t1_map, R2s_flip and ff_flip
%
wx = 64;
wy = 64;
centroid = centroids{1};
t1_map_crop = imcrop(t1_map, [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
R2s_flip_crop = imcrop(R2s_flip, [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
ff_flip_crop = imcrop(ff_flip, [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
remote_crop = imcrop(remote, [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
mi_crop = imcrop(mi, [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);


ff_mi = nonzeros(ff_flip_crop .* mi_crop);
ff_mi(isnan(ff_mi)) = []; 
R2s_mi = nonzeros(R2s_flip_crop .* mi_crop);
t1_mi = nonzeros(t1_map_crop .* mi_crop);
t1_mi(isnan(t1_mi)) = [];
%t2_mi = nonzeros(t2_map .* mi);
%t2_mi(isnan(t2_mi)) = [];

ff_remote = nonzeros(ff_flip_crop .* remote_crop);
ff_remote(isnan(ff_remote)) = []; 
R2s_remote = nonzeros(R2s_flip_crop .* remote_crop);
t1_remote = nonzeros(t1_map_crop .* remote_crop);
t1_remote(isnan(t1_remote)) = [];
%t2_remote = nonzeros(t2_map .* remote);
%t2_remote(isnan(t2_remote)) = [];

figure('Position', [100 100 800 600]);
[ha0, pos0] = tight_subplot(1,2,[.01 -.1],[.01 .01],[.01 .01]);
[ha, pos] = tight_subplot(3,2,[.01 -.1],[.01 .01],[.01 .01]);

axes(ha0(1));
p1 = plot3(R2s_mi, ff_mi, t1_mi); 
p1.LineStyle = "none";
p1.Color = "red";
p1.Marker = "o";

hold on;
p2 = plot3(R2s_remote, ff_remote, t1_remote); 
p2.LineStyle = "none";
p2.Color = "blue";
p2.Marker = "o";
grid on;
xlabel('R2star (s^{-1})'); ylabel('FF'); zlabel('T1 (ms)');

axes(ha0(2));
axis off;
axes(ha(1));
axis off;
axes(ha(3));
axis off;
axes(ha(5));
axis off;
axes(ha(2));
axis off;
imagesc(t1_map_crop); caxis([0 2000]); axis image; axis off; colorbar;
pbaspect([size(t1_map, 2), size(t1_map, 1) 1]);
xticks([]); yticks([]);


axes(ha(4));
imagesc(R2s_flip_crop); caxis([0 100]); axis image; axis off; colorbar;
pbaspect([size(t1_map, 2), size(t1_map, 1) 1]);
xticks([]); yticks([]);

axes(ha(6));
imagesc(ff_flip_crop); caxis([0 0.2]); axis image; axis off; colorbar;
pbaspect([size(t1_map, 2), size(t1_map, 1) 1]);
xticks([]); yticks([]);
colormap gray;

% k=waitforbuttonpress;
% if k == 0
%     dt = findobj(p1,'Type','datatip');
% end

%% On click data tips
dt = findobj(p1,'Type','datatip');
dt_array = [];
for i = 1:length(dt)
    dt_array = [dt_array, dt(i).DataIndex];
end
dt_array = unique(dt_array);

nonzero_idx = find(R2s_flip_crop .* mi_crop ~= 0);
mask_binary = zeros(size(t1_map_crop));
mask_binary(nonzero_idx(dt_array)) = 1;

ax2 = axes('Position', ha(2).Position);
imagesc(ax2, mask_binary, 'AlphaData', mi_crop.*0.4);
pbaspect([size(t1_map_crop, 2), size(t1_map_crop, 1) 1]); colormap(ax2, 'cool');
linkprop([ha(2) ax2], 'Position');
ax2.Visible = 'off';

ax4 = axes('Position', ha(4).Position);
imagesc(ax4, mask_binary, 'AlphaData', mi_crop.*0.4);
pbaspect([size(t1_map_crop, 2), size(t1_map_crop, 1) 1]); colormap(ax4, 'cool');
linkprop([ha(4) ax4], 'Position');
ax4.Visible = 'off';

ax6 = axes('Position', ha(6).Position);
imagesc(ax6, mask_binary, 'AlphaData', mi_crop.*0.4);
pbaspect([size(t1_map_crop, 2), size(t1_map_crop, 1) 1]); colormap(ax6, 'cool');
linkprop([ha(6) ax6], 'Position');
ax6.Visible = 'off';
%% Are same t1 value pixels have the same recovery curve?
%% IR-T1
% t1_map_masked(t1_mi_idx1(t1_mi_idx2(t1_mi_idx3)))
dicom_glob = glob(cat(2, dicom_dir, '/IR_T1_TSE_TI*'));
dicom_glob_reorder = dicom_glob([7,3,5,6,1,2,4]);
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
    'EchoTime',...
    'EchoTrainLength',...
    'InversionTime'
    };

whatsinit = cell(length(dicom_glob_reorder), 1);
slice_data = cell(length(dicom_glob_reorder), 1);
IR_array = zeros(length(dicom_glob_reorder), 1);
for i = 1:length(dicom_glob_reorder)
    f = dicom_glob_reorder{i};
    [whatsinit{i} slice_data{i}] = dicom23D(f, dicom_fields);
    IR_array(i) = slice_data{i}.InversionTime;
end

img_t1_w = zeros(size(whatsinit{1},1), size(whatsinit{1},2), length(whatsinit));
for i = 1:length(whatsinit)
    img_t1_w(:,:,i) = whatsinit{i};  
end
%% Is it?
t1_narrow_mi = t1_mi > 1140 & t1_mi < 1160;
t1_narrow_remote = t1_remote > 1140 & t1_remote < 1160;

figure();
p = plot3(R2s_mi(~t1_narrow_mi), ff_mi(~t1_narrow_mi), t1_mi(~t1_narrow_mi)); 
p.LineStyle = "none";
p.Color = "red";
p.Marker = "o";
hold on;
p2 = plot3(R2s_mi(t1_narrow_mi), ff_mi(t1_narrow_mi), t1_mi(t1_narrow_mi)); 
p2.LineStyle = "none";
p2.Color = "green";
p2.Marker = "o";
grid on;
xlabel('R2star'); ylabel('FF'); zlabel('T1');

t1_map_masked = t1_map .* mi;
t1_mi_temp = nonzeros(t1_map .* mi);
t1_mi_idx1 = find(t1_map_masked ~= 0);
t1_mi_idx2 = find(~isnan(t1_mi_temp));
t1_mi_idx3 = find(t1_narrow_mi);

p3 = plot3(R2s_remote(~t1_narrow_remote), ff_remote(~t1_narrow_remote), t1_remote(~t1_narrow_remote)); 
p3.LineStyle = "none";
p3.Color = "blue";
p3.Marker = "x";

p4 = plot3(R2s_remote(t1_narrow_remote), ff_remote(t1_narrow_remote), t1_remote(t1_narrow_remote)); 
p4.LineStyle = "none";
p4.Color = "cyan";
p4.Marker = "x";
grid on;
xlabel('R2star'); ylabel('FF'); zlabel('T1');

t1_map_remote_masked = t1_map .* remote;
t1_remote_temp = nonzeros(t1_map .* remote);
t1_remote_idx1 = find(t1_map_remote_masked ~= 0);
t1_remote_idx2 = find(~isnan(t1_remote_temp));
t1_remote_idx3 = find(t1_narrow_remote);

img_t1_w_reshape = reshape(img_t1_w, [], length(whatsinit));
t1_w_array = img_t1_w_reshape(t1_mi_idx1(t1_mi_idx2(t1_mi_idx3)),:);
t1_w_array_remote = img_t1_w_reshape(t1_remote_idx1(t1_remote_idx2(t1_remote_idx3)),:);

figure(); 
pp1 = plot(IR_array, t1_w_array.', 'r'); hold on;
pp2 = plot(IR_array, t1_w_array_remote.', 'k');
% legend(pp1, {'MI'});legend(pp2, {'Remote'});

figure(); 
pp1 = plot(IR_array, (t1_w_array ./ max(t1_w_array(:))).', 'r'); hold on;
pp2 = plot(IR_array, (t1_w_array_remote ./ max(t1_w_array_remote(:))).', 'k');

ff_flip = flip(imrotate(ff,-90),2);
ff_masked = ff_flip .* mi;
ff_array = ff_masked(t1_mi_idx1(t1_mi_idx2(t1_mi_idx3)));
ff_remote_masked = ff_flip .* remote;
ff_remote_array = ff_remote_masked(t1_remote_idx1(t1_remote_idx2(t1_remote_idx3)));

R2s_flip = flip(imrotate(R2s,-90),2);
R2s_masked = R2s_flip .* mi;
R2s_array = R2s_masked(t1_mi_idx1(t1_mi_idx2(t1_mi_idx3)));
R2s_remote_masked = R2s_flip .* remote;
R2s_remote_array = R2s_remote_masked(t1_remote_idx1(t1_remote_idx2(t1_remote_idx3)));

t1_array = t1_map_masked(t1_mi_idx1(t1_mi_idx2(t1_mi_idx3)));
%% T1 MOLLI
dicom_glob = glob(cat(2, dicom_dir, '/T1MAP_SAX4_MOCO_T1_0018'));
dicom_glob_reorder = dicom_glob([1]);
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
    'EchoTime',...
    'EchoTrainLength',...
    'InversionTime'
    };

whatsinit = cell(length(dicom_glob_reorder), 1);
slice_data = cell(length(dicom_glob_reorder), 1);
%IR_array = zeros(length(dicom_glob_reorder), 1);
for i = 1:length(dicom_glob_reorder)
    f = dicom_glob_reorder{i};
    [whatsinit{i} slice_data{i}] = dicom23D(f, dicom_fields);
    %IR_array(i) = slice_data{i}.InversionTime;
end

figure(); imagesc(whatsinit{1}); axis image;
%% T1 MOLLI
img_t1_molli = whatsinit{1};
w_cut = 25;
h_cut = 24;
[h,w] = size(img_t1_molli);
img_t1_molli_crop = imcrop(img_t1_molli, [w_cut, h_cut, w-w_cut*2-1, h-h_cut*2-1]);
figure(); imagesc(img_t1_molli_crop);axis image;

img_t1_molli_resize = imresize(img_t1_molli_crop, [size(t1_map)]); 
figure(); imagesc(img_t1_molli_resize); axis image;

J = imtranslate(img_t1_molli_resize,[0, -4]);
figure(); imagesc(J .* mi); axis image; caxis([0 2000]);
figure(); imagesc(J .* remote); axis image; caxis([0 2000]);
%figure(); imagesc(t1_map .* mi); axis image; caxis([0 2000]);
%figure(); imagesc(t1_map); axis image; caxis([0 2000]);
t1_molli_mi = nonzeros(J .* mi);
t1_molli_remote = nonzeros(J .* remote);
figure(); plot(t1_mi, t1_molli_mi, 'o');xlim([800 2000]); ylim([800 2000]);
xlabel('T1SE');ylabel('T1MOLLI'); grid on;
hold on;
plot(t1_remote, t1_molli_remote, 'o');
%%
figure();
p = plot3(R2s_mi, ff_mi, t1_molli_mi); 
p.LineStyle = "none";
p.Color = "red";
p.Marker = "o";
hold on;
p1 = plot3(R2s_remote, ff_remote, t1_molli_remote); 
grid on;
p1.LineStyle = "none";
p1.Color = "blue";
p1.Marker = "x";

% %%
% function mybttnfcn(h,~)
% hf = get(h,'parent');
% b = get(hf,'selectiontype');
% xy = get(gca,'CurrentPoint');
% if strcmpi(b,'normal')
%     text(xy(1,1),xy(1,2),'Left click')
% elseif strcmpi(b,'alt')
%     text(xy(1,1),xy(1,2),'Right click')
% else
%     text(xy(1,1),xy(1,2),'Careful there, crazy man!')
% end
% end