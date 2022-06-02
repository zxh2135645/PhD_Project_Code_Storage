clear all;
close all;
%% Read DICOM file
addpath('../function/');

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

%% Load LRT recon
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params');
%% Display image
slc = 3;
dispim = @(x)fftshift(x(:,:,slc,:),1);

num_seg = 21;
temp = Gr\reshape(Phi(:,num_seg,1,1,:), L, []);
temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);
cw = max(vec(abs(temp)));

figure();
temp = imrotate(temp, 90);
temp = abs(flip(temp(:,:,1),2)/cw); % First Section
temp = uint16(temp*4095);
ax2 = imagesc(temp); axis image; colormap gray;axis off;
%% 6TEs and then next sliceloc
X = dicomread(slice_data{1}(1).Filename);
NEco_old = params.NEco_old;
len = length(slice_data{1})/NEco_old;

% Enhanced DICOM
save_dir = cat(2, fid_path, 'DICOM_L/');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end


echo_f_glob = glob(cat(2, fid_path, '*Echo?_Seg??.mat'));

slc_array = [8 7 6 5 4 3 2 1 14 13 12 11 10 9];
% slc_array = [9 8 7 6 5 4 3 2 1 16 15 14 13 12 11 10];
n_seg = 15;

for te = 1:NEco_old
%for te = 1:1
    f = echo_f_glob{te};
    load(f, 'Gr', 'Phi', 'L', 'U', 'Nx', 'Ny', 'Nz', 'params', 'sizes');    
    for i = 1:len
        slc = slc_array(i);
        dispim = @(x)fftshift(x(:,:,slc,:),1);
        temp = Gr\reshape(Phi(:,num_seg,1,1,:), L, []);
        temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);  % temp is 192x192x15 (#ofSeg)
        
        metadata = dicominfo(slice_data{1}(NEco_old*(i-1)+te).Filename);
        fname = cat(2, save_dir, fid_file(1:18), 'Echo', num2str(te), '_Seg', num2str(num_seg), '_Slc', num2str(i), '_L', '.dcm');
        
        cw = max(vec(abs(temp(:,:,n_seg))));
        temp = imrotate(temp, 90);
        temp = abs(flip(temp(:,:,n_seg),2)/cw);
        temp = uint16(temp*4095);
        
        dicomwrite(temp, fname, metadata, 'CreateMode', 'copy');
    end
end

%% 02/24/2022 (single-slice)
slc_array = [8 7 6 5 4 3 2 1 14 13 12 11 10 9];
slc_array = [2 2 2 2 2 2 2 2 2 2 2 2 2 2];
num_seg_array = [21, 31, 41, 51];
n_seg = 15;
load(cat(2, fid_path, 'FID19652_Sofia_LRT_Mappings_Seg15.mat'));
NEco_old = params.NEco_old;
echo_f_glob = glob(cat(2, fid_path, '*Echo?_Seg??.mat'));

% |Echo1_Seg21|Echo1_Seg31|Echo1_Seg41|Echo1_Seg51|T1Map|T2*Map|
% Enhanced DICOM
save_dir = cat(2, fid_path, 'DICOM_CVI/');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

for te = 1:NEco_old
    %for te = 1:1
    if te < 5
       num_seg = num_seg_array(te); 
    end
    f = echo_f_glob{te};
    load(f, 'Gr', 'Phi', 'L', 'U', 'Nx', 'Ny', 'Nz', 'params', 'sizes');
    for i = 1:(n_seg-1) % should make n_seg = 14 so that n_slc equals n_seg
        slc = slc_array(i);
        dispim = @(x)fftshift(x(:,:,slc,:),1);
        temp = Gr\reshape(Phi(:,num_seg,1,1,:), L, []);
        temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);  % temp is 192x192x15 (#ofSeg)
        cw = max(vec(abs(temp(:,:,i))));
        
        metadata = dicominfo(slice_data{1}(NEco_old*(i-1)+te).Filename);
        if te < 5
            fname = cat(2, save_dir, fid_file(1:18), 'Echo', num2str(1), '_Seg', num2str(num_seg), '_Section', num2str(i), '_ForCVI', '.dcm');
            
            temp = imrotate(temp, 90);
            temp = abs(flip(temp(:,:,i),2)/cw);
            temp = uint16(temp*4095);
            metadata.WindowCenter = 2048;
            metadata.WindowWidth = 4095;
        elseif te == 5
            fname = cat(2, save_dir, fid_file(1:18), 'Echo', num2str(1), '_T1Map', '_Section', num2str(i), '_ForCVI', '.dcm');
            
            temp = uint16(squeeze(map_to_save.t1_map(:,:,slc,i)));
            
            metadata.WindowCenter = 500;
            metadata.WindowWidth = 1000;
        elseif te == 6
            fname = cat(2, save_dir, fid_file(1:18), '_T2starMap', '_Section', num2str(i), '_ForCVI', '.dcm');
            
            temp = uint16(squeeze(map_to_save.t2star(:,:,slc)));
            metadata.WindowCenter = 50;
            metadata.WindowWidth = 100;
        end
        metadata.SmallestImagePixelValue = min(temp(:));
        metadata.LargestImagePixelValue = max(temp(:));
        dicomwrite(temp, fname, metadata, 'CreateMode', 'copy');
    end
end

%% PSIR version
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params');
%%
NEco_old = params.NEco_old;
len = length(slice_data{1})/NEco_old;
slc_input = 3;
slc_array = repmat(slc_input, [1,len]);

Nseg = size(Phi,2);
Nsect = size(Phi,5);
save_dir = cat(2, fid_path, 'DICOM_PSIR_slc', num2str(slc_input), '/');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% odd_array = (1:Nseg/2) * 2 - 1;
% caxis_cell = {[0.3 0.9],[0.2 0.75]};
% Nsect_array = [1, 15];
num_seg_array = [16, 21, 26, 31, 36, 41];
% psir_mat = zeros(Nx, Ny, size(Phi, 2), size(Phi, 5));
for ii = 1:2
    f = echo_f_glob{ii};
    load(f, 'Gr', 'Phi', 'L', 'U', 'Nx', 'Ny', 'Nz', 'params', 'sizes');
for i = 1:1
% for i = len:len
    % nsec = Nsect_array(n);
    % figure('Position',[0, 100, 1800, 800]);
    % [ha, pos] = tight_subplot(1,3,[.01 -.1],[.01 .01],[.01 .01]);
    slc = slc_array(i);
    dispim = @(x) fftshift(x(:,:,slc,:), 1);
    temp = Gr\reshape(Phi(:,:,1,1,:), L, []);
    temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], Nsect);
    % cw = max(vec(abs(temp(:,:,i))));
    
    if ii == 1
        temp_sec1 = temp(:,:,:,i);
    elseif ii == 2
        temp_sec2 = temp(:,:,:,i);
    end
%end
%end
    temp_sec = temp_sec1;
    phase_diff = angle(squeeze(temp_sec1 ./ temp_sec1(:,:,end)));
    %[Iunwrap,seeds] = unwrap2mov(phase_diff(:,:,num_seg_array));
    
    % cos_phase_diff = cos(Iunwrap);
    cos_phase_diff(cos_phase_diff >= 0) = 1;
    cos_phase_diff(cos_phase_diff <  0) = -1;
    
    for te = 1:NEco_old
        num_seg = num_seg_array(te);
        % slc = soi(i);
        
        metadata = dicominfo(slice_data{1}(NEco_old*(i-1)+te).Filename);
        
        % temp = abs(flip(temp(:,:,i),2)/cw);
        % temp = uint16(temp*4095);
        
        phase_temp = angle(squeeze(temp_sec));
        phase_diff = phase_temp - phase_temp(:,:,end);
        
        figure();
        subplot(1,2,1);
        imagesc(angle(temp_sec(:,:,21))); colormap gray;
        subplot(1,2,2);
        imagesc(angle(temp_sec(:,:,end))); colormap gray; %caxis([-pi, pi]);
        
        figure(); imagesc(phase_diff(:,:,21)); colormap gray;
        %[Iunwrap,seeds] = unwrap2mov(phase_diff);

%         cos_phase_diff2 = cos(phase_diff2);
%         cos_phase_diff2(cos_phase_diff2 >= 0) = 1;
%         cos_phase_diff2(cos_phase_diff2 <  0) = -1;
        
        %             temp_sec = imrotate(squeeze(temp_sec), 90);
        %             temp_sec = flip(temp_sec,2);
        %             cos_phase_diff = imrotate(squeeze(cos_phase_diff), 90);
        %             cos_phase_diff = flip(cos_phase_diff, 2);
        %idx = find(soi == slc);
        %centroid = centroids{idx};
        
        %temp_crop = imcrop3(abs(temp(:,:,odd_array,nsec)), [centroid(1)-wx/2, centroid(2)-wy/2, 1, (wx-1), (wy-1), (Nseg/2-1)]);
        temp_psir = abs(temp_sec(:,:,te)) .* cos_phase_diff(:,:,te);
        psir = (temp_psir - min(temp_psir(:))) ./ (max(temp_psir(:)) - min(temp_psir(:)));
        % psir = squeeze(psir(:,:,num_seg));
        psir = imrotate(squeeze(psir), 90);
        psir = flip(psir, 2);
        
        %psir_crop = imcrop3(psir, [centroid(1)-wx/2, centroid(2)-wy/2, 1, (wx-1), (wy-1), (Nseg/2-1)]);
        % cos_phase_diff_crop = imcrop3(cos_phase_diff(:,:,odd_array,nsec), [centroid(1)-wx/2, centroid(2)-wy/2, 1, (wx-1), (wy-1), (Nseg/2-1)]);
        %phase_diff_crop = imcrop3(phase_diff(:,:,odd_array,nsec), [centroid(1)-wx/2, centroid(2)-wy/2, 1, (wx-1), (wy-1), (Nseg/2-1)]);
        %cw1 = 0.5*max(vec(abs(temp(:,:,odd_array,nsec))));
        %cw2 = 0.5*max(vec(phase_diff(:,:,odd_array,nsec)));
        %cw3 = 0.5*max(vec(psir));
        psir = uint16(psir*4095);
        
        metadata.WindowCenter = 2048;
        metadata.WindowWidth = 4095;
        metadata.SmallestImagePixelValue = min(psir(:));
        metadata.LargestImagePixelValue = max(psir(:));
        
        fname = cat(2, save_dir, fid_file(1:18), 'Echo', num2str(1), '_Seg', num2str(num_seg), '_Section', num2str(i), '_ForCVI', '.dcm');
        
        dicomwrite(psir, fname, metadata, 'CreateMode', 'copy');
        
        % psir_mat(:,:,:,i) = temp_sec .* cos_phase_diff;
    end
    
    close all;
end
end

%% save as mat for T1 fitting using psir
psir_mat = zeros(Nx, Ny, Nseg, Nsect);
slc_array = repmat(slc_input, [1,Nsect]);

for i = 1:Nsect
    slc = slc_array(i);
    dispim = @(x) fftshift(x(:,:,slc,:), 1);
    temp = Gr\reshape(Phi(:,:,1,1,:), L, []);
    temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], Nsect);
    
    temp_sec = temp(:,:,:,i);
    phase_temp = angle(squeeze(temp_sec));
    phase_diff = phase_temp - phase_temp(:,:,end);
    cos_phase_diff = cos(phase_diff);
    cos_phase_diff(cos_phase_diff >= 0) = 1;
    cos_phase_diff(cos_phase_diff <  0) = -1;
    
    temp_sec = imrotate(squeeze(temp_sec), 90);
    temp_sec = flip(temp_sec, 2);
    cos_phase_diff = imrotate(squeeze(cos_phase_diff), 90);
    cos_phase_diff = flip(cos_phase_diff, 2);
    
    psir_mat(:,:,:,i) = abs(temp_sec) .* cos_phase_diff;
end

save_dir = cat(2, fid_path, 'MAT_PSIR/');
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

fname = cat(2, save_dir, fid_file(1:18), 'Echo', num2str(1), '_ForT1Fitting', '_slc', num2str(slc_input), '.mat');
save(fname, 'psir_mat');
