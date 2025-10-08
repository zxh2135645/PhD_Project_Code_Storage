clear all;
close all;
%% This script is used for LRT figure production 03/29/2025
% For producing representative images
% Mostly based off Lisbon Acute Data
%% Load Data
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params', 'Hidx', 'RR_int');
%% single slice - slice dimension
dispim = @(x)fftshift(x(:,:,19,:),1);
resp_phase = 1;
card_phase = 1;
temp = Gr\reshape(Phi(:,201,:,resp_phase,end), L, []);
temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
cw = 0.5*max(vec(abs(temp)));


ax1 = implay(abs(temp/cw));
% figure();
% ax2 = imagesc(abs(temp(:,:,1)/cw)); axis image; colormap gray;axis off;
%% Whole heart
for i = 1:Nz
    dispim = @(x)fftshift(x(:,:,i,:),1);
    temp = Gr\reshape(Phi(:,101,card_phase,resp_phase,end), L, []);
    temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
    if i == 1
        temp_wholeheart = temp;
    else
        temp_wholeheart(:,:,i) = temp;
    end
end
cw = 0.5*max(vec(abs(temp_wholeheart)));
temp = fftshift(temp_wholeheart,3);
%% Save representative images as gif
% Save representative images as gif
delay_time = 0.1;
save_path = cat(2, fid_path, 'representative_gif/');
if ~exist(save_path, 'dir')
    mkdir(save_path);
end

% Crop temp to Nx/2 x Ny/2 at the center before GIF creation
crop_Nx = floor(Nx/2);
crop_Ny = floor(Ny/2);
center_x = floor(Nx/2) + 1+20;
center_y = floor(Ny/2) + 1;
x_start = center_x - floor(crop_Nx/2);
y_start = center_y - floor(crop_Ny/2);

temp_cropped = temp(y_start:y_start+crop_Ny-1, x_start:x_start+crop_Nx-1, :);
cw = max(vec(abs(temp_cropped)));

temp_cropped = flip(imrotate(temp_cropped, 90, 'crop'), 2);
fh = figure('Position', [100 100 400 400]);
axis tight manual
for n = 1:size(temp_cropped , 3)
    imagesc(abs(temp_cropped(:,:,n)/cw)); axis image; colormap gray; axis off;
    set(gca,'LooseInset',get(gca,'TightInset'));
    hold on;
    % Overlay frame number at bottom left
    text(crop_Nx-10, crop_Ny-10, num2str(n), 'Color', 'w', 'FontSize', 18, 'FontWeight', 'bold', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    hold off;
    drawnow
    frame = getframe(fh);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    filename = cat(2, fid_file(1:17), ['_representative_#_Seg192_Cart_T1.gif']);
    if n == 1
        imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime',delay_time, 'Loopcount',inf);
    else
        imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime',delay_time, 'WriteMode','append');
    end
end
%%
% fh = figure('Position', [100 100 400 400]);
% axis tight manual
% for n = 1:size(temp_cropped , 3)
%     imagesc(abs(temp_cropped(:,:,n)/cw)); axis image; colormap gray; axis off;
%     set(gca,'LooseInset',get(gca,'TightInset'));
%     drawnow
%     frame = getframe(fh);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     filename = cat(2, fid_file(1:17), '_representative.gif');
%     if n == 1
%         imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime',0.1, 'Loopcount',inf);
%     else
%         imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime',0.1, 'WriteMode','append');
%     end
% end

%%
% cardiac phase and resp phase needs to be encoded
% Dave_D8      [12, 4]
% George_D6    [16, 4]
% George_WK8   [24, 1]
% Ginger_D8    [24, 4]
% Ginger_WK8   [7, 4]
% Carlos_D6    [20, 4]
% Paprika_D8   [1, 4]
% Paprika_WK8  [15, 4]
% Nutmeg_D6    [24, 1]
% Nutmeg_WK8   [18, 4]
% Cinnamon_D8  [6, 4]
% Cinnamon_WK8 [12, 4]
% Chili_D8     [11, 4]
% Sofia_D6     [1, 1]
% Sofia_WK8    [2, 1]
% Paris_D6     [24, 1]
% Paris_WK8    [1, 4]
% Lisbon_D6    [4, 1]
% Lisbon_WK8   [10, 1]
% Jesse_D8     [4, 1]
% Jesse_WK8    [14, 1]

% Donut_WK8     [ , 4]
% Latte_D7      [ , 4]
% Latte_WK8     [18, 1]
% Stilton_WK8   [ , 1]
% Twinkie_WK8   [ , 1]

%% single slice - slice dimension
addpath('../function/');
img_3d = zeros(Ny, Nx, Nz);

for slc = 1:14
    dispim = @(x)fftshift(x(:,:,slc,:),1);
    resp_phase = 1;
    card_phase = 4;
    temp = Gr\reshape(Phi(:,:,card_phase,resp_phase,1), L, []);
    temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
    cw = max(vec(abs(temp)));


    % ax1 = implay(abs(temp/cw));
    img_3d(:,:,slc) = abs(temp/cw);
end


img_3d_trunc = fftshift(img_3d,3);
img_3d_trunc = img_3d_trunc(:,:,2:Nz);
figure(); imshow3D(img_3d_trunc);

%% You can jump directly to Four-chamber View
%%
figure();
h = slice(img_3d_trunc, [], [], 1:13); 
set(h, 'EdgeColor','none', 'FaceColor','interp');
alpha(.3);

%% Heart rate
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'RR_int');
hr = 60/(RR_int/1000)
%% Cardiac Function
clear temp;
soi = [6, 5, 4, 3, 2, 1, 14, 13, 12];
%soi = [3];
%soi = 2;

temtemp_4D = zeros(Ny, Nx, size(Phi, 5), length(soi));
for i = 1:length(soi)
    slc = soi(i);
    dispim = @(x)fftshift(x(:,:,slc,:),1);
    
    temp = Gr\reshape(Phi(:,31,card_phase,resp_phase,:), L, []);
    temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
    
    temtemp_4D(:,:,:,i) = temp;
end

temp_4D = permute(temtemp_4D, [1 2 4 3]);

%% Decide ES and ED
mask = zeros(Ny, Nx, size(Phi, 3));
figure('Position', [100,100,1000,800]);
for i = 1:size(Phi, 3)
    
    imagesc(abs(temp_4D(:,:,3,i))); axis image; colormap gray;
    roi = drawpolygon;
    mask(:,:,i) = createMask(roi); % mask of
end

mask_array = sum(reshape(mask, [], size(Phi, 3)), 1);
[~, ED] = max(mask_array);
[~, ES] = min(mask_array);


%% ED & ES
% ED = 17;
% ES = 6;

ED = 6;
ES = 15;
mask_ed = zeros(Ny, Nx, length(soi));
mask_es = zeros(Ny, Nx, length(soi));
%%
figure('Position', [100,100,1000,800]);
for i = 1:length(soi)
    imagesc(abs(temp_4D(:,:,i,ED))); axis image; colormap gray;
    roi = drawpolygon;
    mask_ed(:,:,i) = createMask(roi);
end

figure('Position', [100,100,1000,800]);
for i = 1:length(soi)
    imagesc(abs(temp_4D(:,:,i,ES))); axis image; colormap gray;
    roi = drawpolygon;
    mask_es(:,:,i) = createMask(roi);
end


figure('Position', [100,100,1000,800]);
for i = 1:length(soi)
    imagesc(abs(temp_4D(:,:,1,ED))); axis image; colormap gray;
    roi = drawpolygon;
    mask_ed(:,:,i) = createMask(roi);
end
%% Calculate
dx = params.dReadoutFOV_mm ./ params.lBaseResolution;
dy = params.dPhaseFOV_mm ./ params.lBaseResolution;
dz = params.dThickness_mm ./ params.lPartitions;
voxel = dx * dy * dz;

mask_ed_array = sum(reshape(mask_ed, [], length(soi)), 1);
mask_es_array = sum(reshape(mask_es, [], length(soi)), 1);
EDV = sum(mask_ed_array) * voxel / 1000; % mL
ESV = sum(mask_es_array) * voxel / 1000; % mL
SV = EDV - ESV;
EF = SV / EDV;

save_path = cat(2, fid_path, 'cardiac_function/');
if ~exist(save_path, 'dir')
   mkdir(save_path); 
end

card_func = struct;
card_func.EDV = EDV;
card_func.ESV = ESV;
card_func.SV = SV;
card_func.EF = EF;
card_func.mask = mask;
card_func.mask_ed = mask_ed;
card_func.mask_es = mask_es;
card_func.ED = ED;
card_func.ES = ES;

filename = cat(2, fid_file(1:53), '_LRT_CardiacFunction', '.mat');
% save(cat(2, save_path, filename), '-struct', 'card_func');
% HR and CO
%% Four Chamber View
point = [Ny/2 Nx/2 13/2];

% [112 78] -> [82 108]
normal = [30 -30 0];
[B,x,y,z] = obliqueslice(img_3d_trunc,point,normal);

figure('Position', [100 100 800 200]);
surf(x,y,z,B,'EdgeColor','None','HandleVisibility','off');
grid on
view([-145 0])
colormap(gray)
xlabel('x-axis')
ylabel('y-axis');
zlabel('z-axis');
title('Slice in 3-D Coordinate Space')
% hold on
% plot3(point(1),point(2),point(3),'or','MarkerFaceColor','r');
% plot3(point(1)+[0 normal(1)],point(2)+[0 normal(2)],point(3)+[0 normal(3)], ...
%     '-b','MarkerFaceColor','b');
% hold off
%legend('Point in the volume','Normal vector')
zlim([-1 13]);

%% Two Chamber View
point = [Ny/2 Nx/2 13/2];

% [112 78] -> [82 108]
normal = [30 30 0];
[B,x,y,z] = obliqueslice(img_3d_trunc,point,normal);

figure('Position', [100 100 800 200]);
surf(x,y,z,B,'EdgeColor','None','HandleVisibility','off');
grid on
view([-45 0])
colormap(gray)
xlabel('x-axis')
ylabel('y-axis');
zlabel('z-axis');
title('Slice in 3-D Coordinate Space')
% hold on
% plot3(point(1),point(2),point(3),'or','MarkerFaceColor','r');
% plot3(point(1)+[0 normal(1)],point(2)+[0 normal(2)],point(3)+[0 normal(3)], ...
%     '-b','MarkerFaceColor','b');
% hold off
%legend('Point in the volume','Normal vector')
zlim([-1 13]);
%% Display LRT

HR = 60 / RR_int * 1000; % bpm

phase_weight = zeros(size(Phi, 3), 1);
for i = 1:size(Phi, 3)
    phase_weight(i) = sum(Hidx == i) ./ numel(Hidx);
end

delay_time = 1/(HR/60) * phase_weight;

save_path = cat(2, fid_path, 'cine_V2/');
if ~exist(save_path, 'dir')
   mkdir(save_path); 
end

wx = 80;
wy = 80;
strn = [2,3,4,5,6,7,8,9,10];

centroids = cell(length(soi) ,1);
figure();
for i = 1:length(soi)
    subplot(3,4,i)
    imagesc(mask_ed(:,:,i))
    s = regionprops(mask_ed(:,:,i),'centroid');
    hold on;
    plot(s.Centroid(1), s.Centroid(2), 'r*');
    hold off;
    
    centroids{i} = round(s.Centroid);
end

fh = figure('Position', [100 100 400 400]);
axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'testAnimated.gif';
for i = 1:length(soi)
    centroid = centroids{i};
    for n = 1:size(Phi, 3)
        % Draw plot for y = x.^n
        
        temp_crop = imcrop(abs(temp_4D(:,:,i,n)), [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
        temp_crop = flip(imrotate(temp_crop, 90),2);
        imagesc(abs(temp_crop)); axis image; colormap gray; axis off;
        set(gca,'LooseInset',get(gca,'TightInset'));
        drawnow
        % Capture the plot as an image
        frame = getframe(fh);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        filename = cat(2, fid_file(1:17), '_LRT_CINE_SA', num2str(strn(i)), '.gif');
        % Write to the GIF File
        if n == 1
            imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime',delay_time(n), 'Loopcount',inf);
        else
            imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime',delay_time(n), 'WriteMode','append');
        end
    end
end

%% Display LRT (T1 evolution)
HR = 60 / RR_int * 1000; % bpm

delay_time = repmat(1/(HR/6), [1, size(Phi,2)]);

save_path = cat(2, fid_path, 'T1_Evolution/');
if ~exist(save_path, 'dir')
   mkdir(save_path); 
end

wx = 80;
wy = 80;
strn = [2,3,4,5,6,7,8,9,10];

centroids = cell(length(soi) ,1);
figure();
for i = 1:length(soi)
    subplot(3,4,i)
    imagesc(mask_ed(:,:,i))
    s = regionprops(mask_ed(:,:,i),'centroid');
    hold on;
    plot(s.Centroid(1), s.Centroid(2), 'r*');
    hold off;
    
    centroids{i} = round(s.Centroid);
end

fh = figure('Position', [100 100 400 400]);
axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'testAnimated.gif';
for i = 1:length(soi)
    centroid = centroids{i};
    for n = 1:size(Phi, 2)
        % Draw plot for y = x.^n
        
        temp_crop = imcrop(abs(temp_4D(:,:,i,n)), [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
        temp_crop = flip(imrotate(temp_crop, 90),2);
        imagesc(abs(temp_crop)); axis image; colormap gray; axis off;
        set(gca,'LooseInset',get(gca,'TightInset'));
        drawnow
        % Capture the plot as an image
        frame = getframe(fh);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        filename = cat(2, fid_file(1:17), '_LRT_CINE_SA', num2str(strn(i)), '.gif');
        % Write to the GIF File
        if n == 1
            imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime',delay_time(n), 'Loopcount',inf);
        else
            imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime',delay_time(n), 'WriteMode','append');
        end
    end
end

%% Display LRT (Time evolution)
HR = 60 / RR_int * 1000; % bpm

delay_time = repmat(1/(HR/60), [1, size(Phi,2)]);

save_path = cat(2, fid_path, 'Time_Evolution/');
if ~exist(save_path, 'dir')
   mkdir(save_path); 
end

wx = 80;
wy = 80;
strn = [2,3,4,5,6,7,8,9,10];

centroids = cell(length(soi) ,1);
figure();
for i = 1:length(soi)
    subplot(3,4,i)
    imagesc(mask_ed(:,:,i))
    s = regionprops(mask_ed(:,:,i),'centroid');
    hold on;
    plot(s.Centroid(1), s.Centroid(2), 'r*');
    hold off;
    
    centroids{i} = round(s.Centroid);
end

fh = figure('Position', [100 100 400 400]);
axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'testAnimated.gif';
for i = 1:length(soi)
    centroid = centroids{i};
    for n = 1:size(Phi, 5)
        % Draw plot for y = x.^n
        
        temp_crop = imcrop(abs(temp_4D(:,:,i,n)), [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
        temp_crop = flip(imrotate(temp_crop, 90),2);
        imagesc(abs(temp_crop)); axis image; colormap gray; axis off;
        set(gca,'LooseInset',get(gca,'TightInset'));
        drawnow
        % Capture the plot as an image
        frame = getframe(fh);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        filename = cat(2, fid_file(1:17), '_LRT_CINE_SA', num2str(strn(i)), '.gif');
        % Write to the GIF File
        if n == 1
            imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime',delay_time(n), 'Loopcount',inf);
        else
            imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime',delay_time(n), 'WriteMode','append');
        end
    end
end

%% Display LRT (Respiratory Motion)
filename = 'FID17210_Lisbon_D_LRT_CardiacFunction.mat';
save_path = cat(2, fid_path, 'cardiac_function/');
load(cat(2, save_path, filename));

HR = 60 / RR_int * 1000; % bpm

soi = [6, 5, 4, 3, 2, 1, 14, 13, 12];
%soi = 2;

delay_time = 2*repmat(1/(HR/60), [1, size(Phi,2)]);

save_path = cat(2, fid_path, 'Resp_motion_V2/');
if ~exist(save_path, 'dir')
   mkdir(save_path); 
end

wx = 120;
wy = 120;
strn = [2,3,4,5,6,7,8,9,10];

centroids = cell(length(soi) ,1);
figure();
for i = 1:length(soi)
    subplot(3,4,i)
    imagesc(mask_ed(:,:,i))
    s = regionprops(mask_ed(:,:,i),'centroid');
    hold on;
    plot(s.Centroid(1), s.Centroid(2), 'r*');
    hold off;
    
    centroids{i} = round(s.Centroid);
end

fh = figure('Position', [100 100 400 400]);
axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'testAnimated.gif';
for i = 1:length(soi)
    centroid = centroids{i};
    for n = 1:size(Phi, 4)
        % Draw plot for y = x.^n
        
        temp_crop = imcrop(abs(temp_4D(:,:,i,n)), [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
        temp_crop = flip(imrotate(temp_crop, 90),2);
        imagesc(abs(temp_crop)); axis image; colormap gray; axis off;
        set(gca,'LooseInset',get(gca,'TightInset'));
        drawnow
        % Capture the plot as an image
        frame = getframe(fh);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        filename = cat(2, fid_file(1:17), '_LRT_CINE_SA', num2str(strn(i)), '.gif');
        % Write to the GIF File
        if n == 1
            imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime',delay_time(n), 'Loopcount',inf);
        else
            imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime',delay_time(n), 'WriteMode','append');
        end
    end
end

%% Save as tiffs
axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'testAnimated.gif';
save_path = cat(2, fid_path, 'Resp_motion_V2');
if ~exist(save_path, 'dir')
   mkdir(save_path); 
end

save_dir = cat(2, save_path, 'Tiffs/');
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

% soi = [3];
% soi = [6, 5, 4, 3, 2, 1, 14, 13, 12];

for i = 1:length(soi)
    for n = 1:size(Phi, 4)
        % Draw plot for y = x.^n

        temp_crop = flip(imrotate(temp_4D(:,:,i,n), 90),2);
        %imagesc(abs(temp_crop)); axis image; colormap gray; axis off;
        %set(gca,'LooseInset',get(gca,'TightInset'));
        cw = 0.8*max(vec(abs(temp_crop)));
        % Capture the plot as an image
        filename = cat(2, fid_file(1:19), '_LRT_CINE_SA', num2str(i), '_RespPhase', num2str(n), '_WholeHeart_DB.tif');
        % Write to the TIFF File
        imwrite(abs(temp_crop)/cw,cat(2,save_dir,filename),'tif');
    end
end
%% Full image LRT 
fh = figure('Position', [100 100 400 400]);
axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'testAnimated.gif';

soi = [3];

for i = 1:length(soi)
    for n = 1:size(Phi, 3)
        % Draw plot for y = x.^n
        
        temp_crop = flip(imrotate(temp_4D(:,:,i,n), 90),2);
        imagesc(abs(temp_crop)); axis image; colormap gray; axis off;
        set(gca,'LooseInset',get(gca,'TightInset'));
        drawnow
        % Capture the plot as an image
        frame = getframe(fh);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        filename = cat(2, fid_file(1:17), '_LRT_CINE_SA', num2str(strn(i)), '_WholeHeart.gif');
        % Write to the GIF File
        if n == 1
            imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime',delay_time(n), 'Loopcount',inf);
        else
            imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime',delay_time(n), 'WriteMode','append');
        end
    end
end
%% Display CINE from CMR
addpath('../function/');
dicom_dir = uigetdir;
folder_glob = glob(cat(2, dicom_dir, '\*'));

%labels = {'RETRO'};
labels = {'TRUEFISP'};
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
% [12, 13, 14, 15, 16, 17, 18, 19, 20]
whatsinit = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit{i} slice_data] = dicom23D(f, dicom_fields);
end

fh = figure('Position', [100 100 300 400]);
axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'testAnimated.gif';
RR = slice_data(end).TriggerTime + slice_data(end).RepetitionTime; % ms
time_resolution_cine = RR / length(slice_data) / 1000;
for i = 1:length(soi)
    for n = 1:size(whatsinit{1},3)
        % Draw plot for y = x.^n
        imagesc(abs(whatsinit{i}(:,:,n))); axis image; colormap gray; axis off;
        set(gca,'LooseInset',get(gca,'TightInset'));
        drawnow
        % Capture the plot as an image
        frame = getframe(fh);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        filename = cat(2, fid_file(1:17), '_CMR_CINE_SA', num2str(strn(i)), '_WholeHeart.gif');
        % Write to the GIF File
        if n == 1
            imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime', time_resolution_cine, 'Loopcount',inf);
        else
            imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime', time_resolution_cine, 'WriteMode','append');
        end
    end
end


%% Full image LRT  - Save each cardiac phase as a tif
% fh = figure('Position', [100 100 400 400]);
axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'testAnimated.gif';
save_path = cat(2, fid_path, 'cine/');
if ~exist(save_path, 'dir')
   mkdir(save_path); 
end

save_dir = cat(2, save_path, 'Tiffs/');
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

% soi = [3];
% soi = [6, 5, 4, 3, 2, 1, 14, 13, 12];

for i = 1:length(soi)
    for n = 1:size(Phi, 3)
        % Draw plot for y = x.^n

        temp_crop = flip(imrotate(temp_4D(:,:,i,n), 90),2);
        %imagesc(abs(temp_crop)); axis image; colormap gray; axis off;
        %set(gca,'LooseInset',get(gca,'TightInset'));
        cw = 0.6*max(vec(abs(temp_crop)));
        % Capture the plot as an image
        filename = cat(2, fid_file(1:19), '_LRT_CINE_SA', num2str(i), '_CardiacPhase', num2str(n), '_WholeHeart_DB.tif');
        % Write to the TIFF File
        imwrite(abs(temp_crop)/cw,cat(2,save_dir,filename),'tif');
    end
end

%%
NEco_old = params.NEco_old; % 6
% for i = 1:Nz
addpath('../function/');
echo_f_glob = glob(cat(2, fid_path, '*Echo???????.mat'));
N_seg = 15;
t2star_map = zeros(Ny, Nx, Nz, N_seg, sizes(5));

%%
% i = input(sprintf('Select Slice of Interest [%d]: ', 3));
%
NEco_old = params.NEco_old; % 6
% for i = 1:Nz
addpath('../function/');
echo_f_glob = glob(cat(2, fid_path, '*USR12_*Echo?????14_lambda0.5.mat'));
N_seg = 14;

for slc = 2:2
    dispim = @(x,st)fftshift(x(:,:,slc,:),1);
    mask_temp = mask_ed;
    if any(mask_temp(:))
        for nt = N_seg:N_seg
            for j = card_phase:card_phase
                for k = resp_phase:resp_phase
                    
                    for neco = 1:NEco_old
                        if neco == 1
                            temp_4D = zeros([192, 192, 1, NEco_old]);
                        end
                        load(echo_f_glob{neco}, 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'params', 'sizes');
                        temp = Gr\reshape(Phi(:,161,j,k,nt), L, []);
                        temp = reshape(reshape(dispim(reshape(U, Ny, Nx, Nz, [])),[],L) * temp, Ny, Nx, [], params.NEco);
                        %cw = 0.5*max(vec(abs(temp)));
                        %figure();
                        %implay(abs(temp)/cw); axis image; colormap gray;
                        temp_4D(:,:,:,neco) = temp;
                    end
                end
            end
        end
    end
end

%% Display LRT (TE evolution)
HR = 60 / RR_int * 1000; % bpm

delay_time = repmat(1/(HR/60), [1, size(Phi,2)]);

save_path = cat(2, fid_path, 'TE_Evolution/');
if ~exist(save_path, 'dir')
   mkdir(save_path); 
end

wx = 80;
wy = 80;
strn = [2,3,4,5,6,7,8,9,10];

centroids = cell(length(soi) ,1);
figure();
for i = 1:length(soi)
    subplot(3,4,i)
    imagesc(mask_ed(:,:,i))
    s = regionprops(mask_ed(:,:,i),'centroid');
    hold on;
    plot(s.Centroid(1), s.Centroid(2), 'r*');
    hold off;
    
    centroids{i} = round(s.Centroid);
end

strn = [3];
soi = 3;
fh = figure('Position', [100 100 400 400]);
axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'testAnimated.gif';
for i = 1:length(soi)
    centroid = centroids{i};
    for n = 1:NEco_old
        % Draw plot for y = x.^n
        
        temp_crop = imcrop(abs(temp_4D(:,:,i,n)), [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
        temp_crop = flip(imrotate(temp_crop, 90),2);
        imagesc(abs(temp_crop)); axis image; colormap gray; axis off;
        set(gca,'LooseInset',get(gca,'TightInset'));
        drawnow
        % Capture the plot as an image
        frame = getframe(fh);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        filename = cat(2, fid_file(1:17), '_LRT_CINE_SA', num2str(strn(i)), '.gif');
        % Write to the GIF File
        if n == 1
            imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime',delay_time(n), 'Loopcount',inf);
        else
            imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime',delay_time(n), 'WriteMode','append');
        end
    end
end


%% Display CINE from CMR (Cropped)
addpath('../function/');
dicom_dir = uigetdir;
folder_glob = glob(cat(2, dicom_dir, '\*'));

%labels = {'RETRO'};
labels = {'TRUEFISP'};
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
% [12, 13, 14, 15, 16, 17, 18, 19, 20]
whatsinit = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit{i} slice_data] = dicom23D(f, dicom_fields);
end

ED = 1;
Ny = size(whatsinit{1},1);
Nx = size(whatsinit{1},2);
mask_ed = zeros(Ny, Nx, length(soi));

figure('Position', [100,100,1000,800]);
for i = 1:length(soi)
    imagesc(whatsinit{1}(:,:,ED)); axis image; colormap gray;
    roi = drawpolygon;
    mask_ed(:,:,i) = createMask(roi);
end

%% 
wx = 100;
wy = 100;
strn = [2,3,4,5,6,7,8,9,10];

centroids = cell(length(soi) ,1);
figure();
for i = 1:length(soi)
    subplot(3,4,i)
    imagesc(mask_ed(:,:,i))
    s = regionprops(mask_ed(:,:,i),'centroid');
    hold on;
    plot(s.Centroid(1), s.Centroid(2), 'r*');
    hold off;
    
    centroids{i} = round(s.Centroid);
end


fh = figure('Position', [100 100 300 400]);
axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'testAnimated.gif';
RR = slice_data(end).TriggerTime + slice_data(end).RepetitionTime; % ms
time_resolution_cine = RR / length(slice_data) / 1000;
for i = 1:length(soi)
    for n = 1:size(whatsinit{1},3)
        % Draw plot for y = x.^n
        centroid = centroids{i};
        im_crop = imcrop(abs(whatsinit{i}(:,:,n)), [centroid(1)-wx/2, centroid(2)-wy/2, (wx-1), (wy-1)]);
        imagesc(im_crop); axis image; colormap gray; axis off;
        set(gca,'LooseInset',get(gca,'TightInset'));
        drawnow
        % Capture the plot as an image
        frame = getframe(fh);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        filename = cat(2, fid_file(1:17), '_CMR_CINE_SA', num2str(strn(i)), '_Cropped.gif');
        % Write to the GIF File
        if n == 1
            imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime', time_resolution_cine, 'Loopcount',inf);
        else
            imwrite(imind,cm,cat(2,save_path,filename),'gif', 'DelayTime', time_resolution_cine, 'WriteMode','append');
        end
    end
end