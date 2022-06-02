% Read mGRE DICOM file and massage it to be the input of FattyRiot

clear all;
close all;
clc;

% Need to run setup_FattyRiot.m in FattyRiot folder
%%
addpath('./function/');
addpath('./function/FattyRiot/');
base_dir = uigetdir;
folder_glob = glob(cat(2, base_dir, '\*'));

labels = {'MGRE'};
mGRE_dicom_fields = {...
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
    'EchoTime',...
    'MagneticFieldStrength',...
    };

for ll = 1:length(labels)
    label = labels{ll};
    idx_array = contains(folder_glob, label);
    
    if any(idx_array)
        num = numel(nonzeros(idx_array));
        ind_array = find(idx_array == 1);
        dst_files = cell(num, 1);
        for i = 1:num
            dst_files{i} = folder_glob{ind_array(i)};
        end
        
        dst_name = ExtractNames(dst_files);
        disp(dst_name);
        
        sel_array = input('Please add an array here:  ');
    end
    
    char_array = num2str(sel_array', '%04.f');
    
    ind_array2 = zeros(size(dst_name, 1), 1);
    for i = 1:size(char_array, 1)
        cha = char_array(i, :);
        ind_array2 = ind_array2 + contains(dst_name, cha);
    end
    
    ind_array3 = find(ind_array2 == 1);
    list_to_read = dst_files(ind_array3);
    
    name_to_compare = ExtractNames(list_to_read);
    
    order_to_read = zeros(length(list_to_read), 1);
    for i = 1:length(list_to_read)
        order_to_read(i) = find(contains(name_to_compare, char_array(i, :)) == 1);
    end
    
    save_dir = GetFullPath(cat(2, base_dir, '\..\img\'));
    if ~exist(save_dir, 'dir')
        mkdir(save_dir)
    end
    
    for i = 1:length(list_to_read)
        f = list_to_read{order_to_read(i)};
        f_vol = cat(2, f, 'VOLUME_IMAGE.mat');
        
        if exist(f_vol, 'file') == 2
            imData = load(f_vol);
            volume_image = imData.volume_image;
            slice_data = imData.slice_data;
            image_meta_data = imData.image_meta_data;
        else
            [volume_image, slice_data, image_meta_data] = dicom23D(f, mGRE_dicom_fields);
        end
        
        if strcmp(label, 'MGRE')
            volume_image = reshape(volume_image, size(volume_image, 1), size(volume_image, 2), 8, []); % 8 TE in mGRE
            volume_image = permute(volume_image, [1,2,4,3]);
            nx = size(volume_image, 1);
            ny = size(volume_image, 2);
            nz = size(volume_image, 3);
            nte = size(volume_image, 4);
            
            %mask = zeros(nx, ny, nz);
            %for slc = 1:nz
            ref_img = volume_image(:,:,40,1);
            ref_img = ref_img / max(ref_img(:));
            figure();
            mask = roipoly(ref_img);
            %end
        end
        
        te_array = zeros(1, 8);
        for num_te = 1:8
            te_array(num_te) = slice_data(num_te).EchoTime / 1000;
        end
        te_array(4) = te_array(1) * 4; %??
        te_array(5) = te_array(1) * 5; %??
        te_array(6) = te_array(1) * 6; %??
        te_array(7) = te_array(1) * 7; %??
        te_array(8) = te_array(1) * 8; %??
        
        mask_3d = repmat(mask, [1,1,nz]); 
        ncoils = 1;
        % To construct imDataParams for Fat-Water Separation
        % acquired images, array of size [nx, ny, nz, ncoils, nTE]
        clear imDataParams
        imDataParams.images = reshape(volume_image, nx, ny, nz, ncoils, nte);
        imDataParams.TE = te_array;
        imDataParams.FieldStrength = slice_data(1).MagneticFieldStrength;
        imDataParams.PrecessionIsClockwise = 1;
        imDataParams.mask = mask_3d;
        
        [FW, INFO] = FattyRiot(imDataParams);
    end
end

%% 

[nx ny nz ncoils nTE] = size(imDataParams.images);
F = FW(:,:,[1:nz]);
W = FW(:,:,[1:nz]+nz);
%%
figure();
subplot(2,1,1); imagesc(F(:,:,53)); axis image; colormap(gray); title('FAT'); %caxis([0 1]);
subplot(2,1,2); imagesc(W(:,:,53)); axis image; colormap(gray); title('WATER'); %caxis([0 1]);

%% 
FF1 = F ./ (F + W);
figure();
imagesc(FF1(:,:,60));

%%
FWs = FW - min(FW(:));
FWs = FWs / max(FWs(:));
Fs = FWs(:,:,[1:nz]);
Ws = FWs(:,:,[1:nz]+nz);
FF2 = Fs ./ (Fs + Ws);
figure();
imagesc(FF2(:,:,60));
