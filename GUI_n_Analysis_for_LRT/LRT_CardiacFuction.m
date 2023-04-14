clear all;
close all;

%% Load Data
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params', 'Hidx', 'RR_int');
%% single slice - slice dimension
dispim = @(x)fftshift(x(:,:,1,:),1);
resp_phase = 4;

temp = Gr\reshape(Phi(:,36,:,resp_phase,end), L, []);
temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
cw = max(vec(abs(temp)));


ax1 = implay(abs(temp/cw));
% figure();
% ax2 = imagesc(abs(temp(:,:,1)/cw)); axis image; colormap gray;axis off;
%% Heart rate
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'RR_int');
hr = 60/(RR_int/1000)
%% Cardiac Function
soi = [6, 5, 4, 3, 2, 1, 14, 13, 12];
temtemp_4D = zeros(Ny, Nx, size(Phi, 3), length(soi));
for i = 1:length(soi)
    slc = soi(i);
    dispim = @(x)fftshift(x(:,:,slc,:),1);
    
    temp = Gr\reshape(Phi(:,41,:,resp_phase,end), L, []);
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

ED = 24;
ES = 12;
mask_ed = zeros(Ny, Nx, length(soi));
mask_es = zeros(Ny, Nx, length(soi));

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
save(cat(2, save_path, filename), '-struct', 'card_func');
% HR and CO
%% Display LRT
HR = 60 / RR_int * 1000; % bpm

phase_weight = zeros(size(Phi, 3), 1);
for i = 1:size(Phi, 3)
    phase_weight(i) = sum(Hidx == i) ./ numel(Hidx);
end

delay_time = 1/(HR/60) * phase_weight;

save_path = cat(2, fid_path, 'cine/');
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

%% Full image LRT 
fh = figure('Position', [100 100 400 400]);
axis tight manual % this ensures that getframe() returns a consistent size
% filename = 'testAnimated.gif';
for i = 1:length(soi)
    centroid = centroids{i};
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
%% Display CINE
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
