clear all;
close all;
%% Load Data
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params', 'Hidx', 'RR_int');
%% single slice - slice dimension
% Sofia_D5
slc = 3;
dispim = @(x)fftshift(x(:,:,slc,:),1);
resp_phase = 1;
card_phase = 13;

% Lisbon_D6
slc = 2;
dispim = @(x)fftshift(x(:,:,slc,:),1);
resp_phase = 1;
card_phase = 15;

% Paris_D6
slc = 2;
dispim = @(x)fftshift(x(:,:,slc,:),1);
resp_phase = 1;
card_phase = 23;

% Jesse_D8
slc = 2;
dispim = @(x)fftshift(x(:,:,slc,:),1);
resp_phase = 1;
card_phase = 2;

% George_D6
slc = 2;
dispim = @(x)fftshift(x(:,:,slc,:),1);
resp_phase = 4;
card_phase = 16;

% Chili_D8
slc = 2;
dispim = @(x)fftshift(x(:,:,slc,:),1);
resp_phase = 4;
card_phase = 18;

% Nutmeg_D6
slc = 2;
dispim = @(x)fftshift(x(:,:,slc,:),1);
resp_phase = 1;
card_phase = 22;

% Ginger_D8
slc = 2;
dispim = @(x)fftshift(x(:,:,slc,:),1);
resp_phase = 4;
card_phase = 24;

% Dave_D8
slc = 3;
dispim = @(x)fftshift(x(:,:,slc,:),1);
resp_phase = 4;
card_phase = 12;

% Carlos_D6
slc = 3;
dispim = @(x)fftshift(x(:,:,slc,:),1);
resp_phase = 4;
card_phase = 1;

% Paprika_D8
slc = 3;
dispim = @(x)fftshift(x(:,:,slc,:),1);
resp_phase = 4;
card_phase = 1;

% Cinnamon_D8
slc = 3;
dispim = @(x)fftshift(x(:,:,slc,:),1);
resp_phase = 4;
card_phase = 6;

%%
num_dy = 15;
temp = Gr\reshape(Phi(:,36,card_phase,resp_phase,:), L, []);
temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
cw = max(vec(abs(temp)));

ax1 = implay(abs(temp/cw));
%%
temtemp_4D = zeros(Ny, Nx, size(Phi, 2), num_dy);
for i = 1:num_dy
    
    temp = Gr\reshape(Phi(:,:,card_phase,resp_phase,i), L, []);
    temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
    
    temtemp_4D(:,:,:,i) = temp;
end

temp_4D = temtemp_4D;
%% Remote ROI
cw = max(vec(abs(temp_4D)));
% Before that, need to match slices between LRT and conventional EGE 

figure('Position', [100,100,1000,800]);
for i = num_dy:num_dy
    imagesc(abs(temp_4D(:,:,36,i))); axis image; colormap gray;
    roi = drawpolygon;
    mask_remote = createMask(roi);
end

%% Find nulling point
num_t1 = 192;
null_array = zeros(num_t1, num_dy);
for i = 1:num_dy
    
    remote = temp_4D(:,:,:,i) .* mask_remote;
    
    for j = 1:num_t1
        null_array(j,i) = mean(abs(nonzeros(remote(:,:,j))));
    end
end

[M,I] = min(null_array);

figure();
for i = 1:num_dy
    ege = temp_4D(:,:,I(i),i);
    % cw = max(vec(abs(ege)));
    subplot(3,5,i);
    imagesc(abs(ege)); colormap gray;
end

%% Write EGE into DICOM
base_dir = uigetdir;
folder_glob = glob(cat(2, base_dir, '\*'));

labels = {'LRT'};

label = labels{1};
idx_array = contains(folder_glob, label);

[list_to_read, order_to_read] = NamePicker(folder_glob(idx_array));

whatsinit = cell(length(list_to_read), 1);
slice_data = cell(length(list_to_read), 1);
for i = 1:length(list_to_read)
    f = list_to_read{order_to_read(i)};
    [whatsinit{i}, slice_data{i}] = dicom23D(f);
end
%%
slc_array = [1 2 3 4 5 6 7 8 9 10 11 12 13 14];

% Specifically for Jesse_D8, George_D6, Clili_D8
I = repmat([36], [1,15]);
% For carlos_D6, George_D6 (slc=2)
I = repmat([31], [1,15]);
% Specifically for Ginger_D8, 
%I(end) = 31;
X = dicomread(slice_data{1}(1).Filename);
NEco_old = params.NEco_old;
len = length(slice_data{1})/NEco_old;

save_dir = GetFullPath(cat(2, fid_path, 'DICOM_EGE_CVI_V2/'));
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

num_dy_array = 1:num_dy;
te = 1;
for i = 1:num_dy
    if i > len
        te = 2;
    end
    % metadata = dicominfo(slice_data{1}(NEco_old*(i-1)+1).Filename);
    if i <= len
        metadata = dicominfo(slice_data{1}(NEco_old*(i-1)+te).Filename);
    else
        metadata = dicominfo(slice_data{1}(NEco_old*(i-len-1)+te).Filename);
    end
    %fname = GetFullPath(cat(2, save_dir, fid_file(1:18), 'Seg', num2str(num_seg), '_Slc', num2str(i), '_LGE', '.dcm'));
    fname = GetFullPath(cat(2, save_dir, fid_file(1:18), 'Seg', num2str(I(i)), '_Dynamic', num2str(i), '_EGE', '.dcm'));


    num_array_flip = flip(num_dy_array);
    slc = num_array_flip(i);
    temp = imrotate(temp_4D(:,:,I(i),i), 90);
    cw = max(vec(abs(temp_4D(:,:,I(i),i))));

    % temp = temp_4D(:,:,slc,te);
    temp = abs(flip(temp,2)./cw);
    temp = uint16(temp*4095);
    metadata.WindowCenter = 2048;
    metadata.WindowWidth = 4095;

    % dicomwrite(temp, fname, metadata, 'CreateMode', 'copy');
    metadata.SmallestImagePixelValue = min(temp(:));
    metadata.LargestImagePixelValue = max(temp(:));
    dicomwrite(temp, fname, metadata, 'CreateMode', 'copy');
end

