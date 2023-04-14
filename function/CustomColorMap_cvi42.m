%clear all;
%close all;

function Cmap = CustomColorMap_cvi42()
Cmap_raw = [0 0 255; 0 252 255; 0 255 170; 0 255 139; 0 255 1; 14 255 0; 184 255 0; 212 255 0; 231 255 0; 254 255 0; 255 244 0; 255 152 0; 255 124 0; 255 51 0; 192 36 0; 192 2 0]./255;
       % 80        55++       55+        55-        48+      48-        40++       40+        40-        35+        35-         25+        25-        10+      10-       0

% rgbplot(Cmap);

vec =[80; 56; 55; 54; 48; 47; 41; 40; 39; 35; 34; 25; 24; 10; 9; 0]./80 * 100;
% [      100;       83;       68;       44;       30;       15;        0];
N = 256;
%vec_reverse = flip(vec,1);
%Cmap_raw_reverse = flip(Cmap_raw,1);
Cmap = interp1(vec, Cmap_raw, linspace(0, 100, N), 'pchip');
% Cmap = interp1(vec_reverse, Cmap_raw_reverse, linspace(100, 0, N), 'pchip');

% I = imread('pout.tif');
% figure();
% imagesc(I); colormap(Cmap);

end
%% Example
% 
% base_dir = uigetdir;
% folder_glob = glob(cat(2, base_dir, '\*'));
% 
% labels = {'T2STAR'};
% 
% label = labels{1};
% idx_array = contains(folder_glob, label);
% 
% [list_to_read, order_to_read] = NamePicker(folder_glob(idx_array));
% 
% whatsinit = cell(length(list_to_read), 1);
% slice_data = cell(length(list_to_read), 1);
% for i = 1:length(list_to_read)
%     f = list_to_read{order_to_read(i)};
%     [whatsinit{i}, slice_data{i}] = dicom23D(f);
% end
% 
% 
% figure(); imagesc(whatsinit{1}(:,:,2)*0.1); caxis([0 80]);
% figure(); imagesc(whatsinit{1}(:,:,2)*0.1); caxis([0 80]); colormap(Cmap);
