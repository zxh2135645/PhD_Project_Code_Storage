clear all;
close all;
% Get thickness of the hemorrhage

addpath('../function/')
labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};
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

subject_name_cell = {'18P90', '18P93', '20P40', '20P10_Exvivo7', '20P11_Exvivo6', '18P92', '18P94_Exvivo3', '18P95', '17P73', '20P48'};

% Load masks
mask_cell = cell(length(subject_name_cell), 1);
for i = 1:length(subject_name_cell)
    subject_data_dir = GetFullPath(cat(2, data_dir, subject_name_cell{i}, '/'));
    mask_cell{i} = load(cat(2, subject_data_dir, 'mask.mat'));
end
%% Read dicom files
subject_name_cell = {'18P90', '18P93', '20P40', '20P10_Exvivo7', '20P11_Exvivo6', '18P92', '18P94_Exvivo3', '18P95', '17P73', '20P48'};
avg_num_cell = {'Avg0016', 'Invivo'};
avg_name = avg_num_cell{1};
[70, 73, 69, 79, 71, 65, 94, 65, 79, 85]
whatsinit = cell(length(subject_name_cell), 1);
% for i = 1:length(subject_name_cell)
for i = 5:5
    subject_name = subject_name_cell{i};
    disp(subject_name);
    base_dir = uigetdir;
    folder_glob = glob(cat(2, base_dir, '\*'));
    [list_to_read, order_to_read] = NamePicker(folder_glob);
    f = list_to_read{order_to_read(1)};
    whatsinit{i} = dicom23D(f);
end

%% 
n = 5;
img = whatsinit{n};
myo_mask = mask_cell{n}.mask_struct(1).myo_mask;
mi_mask = mask_cell{n}.mask_struct(1).mi_mask;
remote_mask = mask_cell{n}.mask_struct(1).remote_mask;

thresh = mean(nonzeros(img.*remote_mask)) - 2*std(nonzeros(img.*remote_mask));
hemo_mask = (img < thresh).*myo_mask.*mi_mask;
se = strel('disk',1);
hemo_mask_eroded = imopen(hemo_mask,se);

out = bwskel(myo_mask>0);
stat = regionprops(myo_mask,'centroid');
x_centroid = stat.Centroid(1);
y_centroid = stat.Centroid(2);

figure();
imagesc(myo_mask + out.*mi_mask); axis image;
hold on; plot(x_centroid, y_centroid, 'rx');


centerline = out.*mi_mask;
indd = find(centerline);
[row, col] = ind2sub(size(centerline), indd);


y_coord = size(centerline, 1) - row + 1;
x_coord = col;
figure(); plot(x_coord, y_coord, 'o');

x_min = x_coord(1);
idx_x_min = find(x_coord == x_min);
idx_y_min = find(y_coord(idx_x_min) == min(y_coord(idx_x_min)));

y1 = y_coord(find(y_coord <= y_coord(idx_y_min)));
y2 = y_coord(find(y_coord > y_coord(idx_y_min)));
x1 = x_coord(find(y_coord <= y_coord(idx_y_min)));
x2 = x_coord(find(y_coord > y_coord(idx_y_min)));
hold on;
plot(x1, y1, 'rx');
%%
figure(); plot(x_coord, y_coord, 'o');
p1 = polyfit(x1,y1,4);
f1 = polyval(p1, x1);
p2 = polyfit(x2,y2,8);
f2 = polyval(p2, x2);
hold on;
plot(x1, f1, 'r-');
plot(x2, f2, 'b-');

k1 = 20;
dy = diff(f1)./diff(x1);
tang = (x1-x1(k1))*dy(k1)+f1(k1);
plot(x1, tang);
x_centroid = (max(x_coord) + min(x_coord))/2;
y_centroid = (max(y_coord) + min(y_coord))/2;
plot(x_centroid, y_centroid, 'rx');

z = x_centroid - x_coord + 1i*(y_centroid - y_coord);
theta = angle(z);

[theta_sorted idx_sorted] = sort(theta, 'descend');
x_coord_sorted = x_coord(idx_sorted);
y_coord_sorted = y_coord(idx_sorted);

%%
% figure(); plot(x_coord, y_coord, 'o');
% hold on;
% k = 1:5:length(x_coord_sorted);
% 
% for i = 1:length(k)
%     
%     plot([x_centroid, x_coord_sorted(k(i))], [y_centroid, y_coord_sorted(k(i))], 'k');
%     plot(x_coord_sorted(k(i)), y_coord_sorted(k(i)), 'ro');
%     pause(.5);
% end

%%
% Interpolate line
k = 1:5:length(x_coord_sorted);
figure();
ax1 = subplot(1,3,1);
imshow(myo_mask)
c_cell = cell(length(k), 1);
c_hemo_cell = cell(length(k), 1);
clen_prev = 0;
c_hemo_len_prev = 0;

for i = 1:length(k)
    slope = (y_centroid - y_coord_sorted(k(i))) ./ (x_centroid - x_coord_sorted(k(i)));
    b = -(slope * x_coord_sorted(k(i)) - y_coord_sorted(k(i)));
    x_interp = x_coord_sorted(k(i)) - (x_centroid - x_coord_sorted(k(i)));
    y_interp = x_interp*slope + b;

    x = [x_centroid, x_interp];
    y = size(centerline, 1) - [y_centroid, y_interp] + 1;
    c = improfile(myo_mask, x, y);
    % c_hemo = improfile(hemo_mask, x, y);
    c_hemo = improfile(hemo_mask_eroded, x, y);

    hold on;
    subplot(1,3,1)
    plot(x,y,'r');
    
    if i == 1
        ax2 = subplot(1,2,2);
    end

    subplot(1,2,2, ax2)
    plot(c(:,1,1),'r');
    %pause(.5);
    c_cell{i} = c;
    c_hemo_cell{i} = c_hemo;

    clen = length(c);
    if clen > clen_prev
        c_max = clen;
    end
    clen_hemo = length(c_hemo);
    if clen_hemo > c_hemo_len_prev
        c_max_hemo = clen_hemo;
    end
end

% zero-filling c_cell
c_mat = zeros(length(k), c_max);
c_hemo_mat = zeros(length(k), c_max_hemo);
for i = 1:length(k)
    l = length(c_cell{i});
    l_hemo = length(c_hemo_mat);

    c_mat(i,1:l) = c_cell{i};
    c_hemo_mat(i,1:l) = c_hemo_cell{i};
end



%
[X,Y] = meshgrid(1:size(c_mat,2),1:size(c_mat,1));
figure(); subplot(1,2,1); imagesc(c_mat);
subplot(1,2,2); mesh(X,Y,c_mat);
% 0.2717 % 18P93

myo_thick = sum(c_mat,2) .* 0.2604;
mean(nonzeros(myo_thick))

[X,Y] = meshgrid(1:size(c_hemo_mat,2),1:size(c_hemo_mat,1));
figure(); subplot(1,2,1); imagesc(c_hemo_mat);
subplot(1,2,2); mesh(X,Y,c_hemo_mat);

hemo_thick = sum(c_hemo_mat,2) .* 0.2604;
mean(nonzeros(hemo_thick))

%hold on;
%plot(medfilt1(x_coord,3), medfilt1(y_coord,3), 'x');
% dy = diff(y_coord)./diff(x_coord);
% t = 260:400;
% k = 100;
% tang = (t-t(k))*dy(k)+y_coord(k);
% hold on;
% plot(t, tang);
% hold off;
% xlim([0, size(centerline, 2)]); ylim([0, size(centerline, 1)]);