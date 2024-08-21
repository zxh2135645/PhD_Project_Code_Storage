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
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params', 'Hidx', 'RR_int');
%% Display image (Write DICOM)
slc_array = [8 7 6 5 4 3 2 1 14 13 12 11 10 9];
slc_array = [1 2 3 4 5 6 7 8 9 10 11 12 13 14];
num_seg_array = [16, 21, 26, 31, 36, 41];
%num_seg_array = [141, 151, 161, 171, 181, 191];
temtemp_4D = zeros(Ny, Nx, length(num_seg_array), length(slc_array));

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

for i = 1:length(slc_array)
    slc = slc_array(i);
    dispim = @(x)fftshift(x(:,:,slc,:),1);
    for num = 1:length(num_seg_array)
        temp = Gr\reshape(Phi(:,num_seg_array(num),20,4,end), L, []);
        temp = reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp, Ny, Nx, [], params.NEco);

        temtemp_4D(:,:,num,i) = temp;
    end
end

temp_4D = permute(temtemp_4D, [1 2 4 3]);

% Write LGE into DICOM
X = dicomread(slice_data{1}(1).Filename);
NEco_old = params.NEco_old;
len = length(slice_data{1})/NEco_old;

save_dir = GetFullPath(cat(2, fid_path, 'DICOM_LGE_CVI/'));
%save_dir = GetFullPath(cat(2, fid_path, '../DICOM_T2star_CVI/'));
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

for te = 1:length(num_seg_array)
    num_seg = num_seg_array(te); 
    for i = 1:len
        %metadata = dicominfo(slice_data{1}(NEco_old*(i-1)+1).Filename);
        metadata = dicominfo(slice_data{1}(NEco_old*(i-1)+te).Filename);
        fname = GetFullPath(cat(2, save_dir, fid_file(1:18), 'Seg', num2str(num_seg), '_Slc', num2str(i), '_LGE', '.dcm'));
        %fname = GetFullPath(cat(2, save_dir, fid_file(1:18), 'Seg', num2str(num_seg), '_Slc', num2str(i), '_PostCon_T2star', '.dcm'));
        
        
        slc_array_flip = flip(slc_array);
        slc = slc_array_flip(i);
        temp = imrotate(temp_4D(:,:,slc,te), 90);
        cw = max(vec(abs(temp_4D(:,:,slc,te))));

        % temp = temp_4D(:,:,slc,te);
        temp = abs(flip(temp,2)./cw);
        temp = uint16(temp*4095);
        metadata.WindowCenter = 2048;
        metadata.WindowWidth = 4095;

        dicomwrite(temp, fname, metadata, 'CreateMode', 'copy');
        metadata.SmallestImagePixelValue = min(temp(:));
        metadata.LargestImagePixelValue = max(temp(:));
        dicomwrite(temp, fname, metadata, 'CreateMode', 'copy');
    end
end

%% Should be deprecated (Just for reference)
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