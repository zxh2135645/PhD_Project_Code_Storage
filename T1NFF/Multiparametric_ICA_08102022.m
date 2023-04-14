clear all;
close all;

dicom_dir = uigetdir();
t1_folder = 'T1MAP_SAX4_MOCO_0017';
t2_folder = 'T2MAP_SAX4_MOCO_0050';
t2star_folder = 'T2_MULTIECHO_2D_SAX4_0079';

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
    'EchoTime',...
    'MagneticFieldStrength',...
    'ImagingFrequency',...
    'ImageType','RescaleSlope','RescaleIntercept',...
    'EchoNumber',...
    'InversionTime',...
    };

[img_t1 slice_data_t1] = dicom23D(cat(2, dicom_dir, '/', t1_folder, '/'), dicom_fields);
[img_t2 slice_data_t2] = dicom23D(cat(2, dicom_dir, '/', t2_folder, '/'), dicom_fields);
[img_t2star slice_data_t2star] = dicom23D(cat(2, dicom_dir, '/', t2star_folder, '/'), dicom_fields);
[img_t2star_mag slice_data_t2star_mag] = dicom23D(cat(2, dicom_dir, '/', t2star_folder, '/'), dicom_fields);
[img_t2star_phase slice_data_t2star_phase] = dicom23D(cat(2, dicom_dir, '/', t2star_folder, '/'), dicom_fields);

img_t1_crop = img_t1(:,5:end-4,:);
img_t2_crop = img_t2(:,5:end-4,:);

%%
[B0file, B0path] = uigetfile('*.dat');
filename = [B0path, '/', B0file];
twix_obj_in = mapVBVD(filename, 'removeOS'); % return all image-data:
if (length(twix_obj_in)>1)% R.Y. avoid adj coil sensitivity
    for  k=1:length(twix_obj_in)
        if (~strcmp(twix_obj_in{k}.hdr.MeasYaps.tSequenceFileName,'%AdjustSeq%/AdjCoilSensSeq') )
            twix_obj=twix_obj_in{k};
        end
    end
else
    twix_obj=twix_obj_in;
end

rawdata = squeeze(twix_obj.image(''));

if twix_obj.image.NSeg == 1 && twix_obj.image.NAve == 1
    [NumRO, NumCh, NumPE, NumSlices, NumEchos] = size(rawdata);
elseif twix_obj.image.NSeg == 1 && twix_obj.image.NAve > 1
    [NumRO, NumCh, NumPE, NumSlices, NumAve, NumEchos] = size(rawdata);
elseif twix_obj.image.NAve == 1 && (twix_obj.image.NPar && twix_obj.image.NPar == 1)
    [NumRO, NumCh, NumPE, NumEchos, NSeg] = size(rawdata);
    NumSlices = 1;
end
clear twix_obj_in rawdata;

nameSession         = ['-' num2str(twix_obj.hdr.Config.MeasUID)];
PatientName         = twix_obj.hdr.Meas.tPatientName; % Need to work on it

xspace_glob = glob(cat(2, B0path, '/../xdata/', 'xspace', nameSession, '/*'));
xspace_cell = cell(length(xspace_glob), 1);

for eco = 1:length(xspace_glob)
    xspace_cell{eco} = load(xspace_glob{eco});
end

N = 1;
RawdataFT = zeros([size(xspace_cell{1}.data), 1, NumEchos]);

for eco = 1:NumEchos
    RawdataFT(:,:,:,eco) = xspace_cell{eco}.data;
end

RawdataFT = permute(RawdataFT, [2,3,4,1,5]);% NumRO, NumPE, NumSlices, NumCh, NumEcho
RawdataFT = flip(flip((RawdataFT),1),2);% NumRO, NumPE, NumSlices, NumCh, NumEcho
RawdataFT = permute(RawdataFT, [4,1,2,5,3]);

% XZ 06/08/2022
iField = permute(RawdataFT,[2 3 5 4 1])*100; % NumCh, NumRO, NumPE, NumEcho, NumSlices =>  NumRO, NumPE, NumSlices, NumEcho, NumCh


if size(iField,5)>1
    % combine multiple coils together, assuming the coil is the fifth dimension
    iField = sum(iField.*conj( repmat(iField(:,:,:,1,:),[1 1 1 size(iField,4) 1])),5);  %
    mag_iField = abs(iField);
    mag_iField_norm = mag_iField ./ max(mag_iField(:));

    Fieldmap_eddy  = sqrt(sqrt((iField(:,:,:,2)./iField(:,:,:,1))./(iField(:,:,:,3)./iField(:,:,:,2))));
    iField_uneddy = zeros(size(iField));
    iField_uneddy(:,:,:,1:2:end) = iField(:,:,:,1:2:end).*Fieldmap_eddy;
    iField_uneddy(:,:,:,2:2:end) = iField(:,:,:,2:2:end)./Fieldmap_eddy;
    iField = mag_iField_norm.*exp(1i*angle(iField_uneddy));
end

%figure(); imagesc(angle(iField(:,:,1,2)./iField(:,:,1,1))); caxis([-1 1]);
%figure(); imagesc(angle(iField(:,:,1,3)./iField(:,:,1,2))); caxis([-1 1]);

t2star_mag = squeeze(abs(iField));
t2star_phase = squeeze(angle(iField));

t2star_mag = squeeze(real(iField));
t2star_phase = squeeze(imag(iField));

%% Draw Masks
mask_f = GetFullPath(cat(2, dicom_dir, '/../mask_heart.mat'));
mask_myo = zeros(size(img_t1_crop,1), size(img_t1_crop,2));
if ~exist(mask_f)
    for slc = 1:1
        figure();
        imagesc(abs(img_t1_crop(:,:,slc)));axis image;
        heart = drawpolygon;
        mask_heart(:,:,slc) = createMask(heart);
    end
    mask_heart = mask_epi - mask_endo;
    save(mask_f, 'mask_heart');
else
    load(mask_f);
end


figure(); imagesc(img_t1_crop(:,:,1) .* mask_heart); axis image;
figure(); imagesc(img_t2_crop(:,:,1) .* mask_heart); axis image;
figure(); imagesc(img_t2star(:,:,1) .* mask_heart); axis image;
%%
img_t1_crop_masked = img_t1_crop .* mask_heart;
img_t2_crop_masked = img_t2_crop .* mask_heart;
img_t2star_masked = img_t2star .* mask_heart;
img_t2star_mag_masked = t2star_mag .* mask_heart;
img_t2star_phase_masked = t2star_phase .* mask_heart;


img_t1_crop_norm = img_t1_crop_masked ./ max(img_t1_crop_masked(:));
img_t2_crop_norm = img_t2_crop_masked ./ max(img_t2_crop_masked(:));
img_t2star_norm = img_t2star_masked ./ max(img_t2star_masked(:));
img_t2star_mag_norm = img_t2star_mag_masked ./ max(img_t2star_mag_masked(:));
img_t2star_phase_norm = (img_t2star_phase_masked - min(img_t2star_phase_masked(:))) ./ max(img_t2star_phase_masked(:));


sz = size(img_t2star);
sz = sz(1:2);

[x,y] = find(mask_heart>0);
%%
maskk = mask_heart;
img_t1_array = zeros(size(img_t1_crop_norm,3), sum(maskk(:)));
N = size(img_t1_crop_norm,3) + size(img_t2_crop_norm,3) + size(img_t2star_norm,3);
TE_array = zeros(N, sum(maskk(:)));
TI_array = zeros(N, sum(maskk(:)));
T2prep_dur = zeros(N, sum(maskk(:)));

for i = 1:size(img_t1_crop_norm,3)
   temp = img_t1_crop_norm(:,:,i);
   img_t1_array(i,:) = temp(maskk>0);
   TI_array(i,:) = repmat(slice_data_t1(i).InversionTime, 1, sum(maskk(:)));
   TE_array(i,:) = repmat(slice_data_t1(i).EchoTime, 1, sum(maskk(:)));
   T2prep_dur(i,:) = repmat(0, 1, sum(maskk(:)));
end

img_t2_array = zeros(size(img_t2_crop_norm,3), sum(maskk(:)));
tau_dur = [0, 30, 55];
for i = 1:size(img_t2_crop_norm,3)
   temp = img_t2_crop_norm(:,:,i);
   img_t2_array(i,:) = temp(maskk>0);
   TI_array(i+size(img_t1_crop_norm,3),:) = repmat(0, 1, sum(maskk(:)));
   TE_array(i+size(img_t1_crop_norm,3),:) = repmat(slice_data_t2(i).EchoTime, 1, sum(maskk(:)));
   T2prep_dur(i+size(img_t1_crop_norm,3),:) = repmat(tau_dur(i), 1, sum(maskk(:)));
end

img_t2star_mag_array = zeros(size(img_t2star_norm,3), sum(maskk(:)));
img_t2star_phase_array = zeros(size(img_t2star_norm,3), sum(maskk(:)));

for i = 1:size(img_t2star_norm,3)
    temp = img_t2star_mag_norm(:,:,i);
    img_t2star_mag_array(i,:) = temp(maskk>0);
    temp = img_t2star_phase_norm(:,:,i);
    img_t2star_phase_array(i,:) = temp(maskk>0);

    TI_array(i+size(img_t1_crop_norm,3)+size(img_t2_crop_norm,3),:) = repmat(0, 1, sum(maskk(:)));
    TE_array(i+size(img_t1_crop_norm,3)+size(img_t2_crop_norm,3),:) = repmat(slice_data_t2star(i).EchoTime, 1, sum(maskk(:)));
    T2prep_dur(i+size(img_t1_crop_norm,3)+size(img_t2_crop_norm,3),:) = repmat(0, 1, sum(maskk(:)));
end

X1 = [img_t1_array.', img_t2_array.', img_t2star_mag_array.', img_t2star_phase_array.']; % 19 features
X2 = [X1, TI_array.', TE_array.', T2prep_dur.']; % 76 features features

%% ica
q = 4;
% K1 = fft(X1);
% KK1 = real(K1);
% KK2 = imag(K1);
% KK = [KK1;KK2];
mdl = rica(X1, q);
Z1 = transform(mdl, X1);
feature_map = zeros(size(img_t1_crop,1), size(img_t1_crop,2), size(Z1, 2));
for i = 1:size(Z1, 1)
    for j = 1:size(Z1,2)
        % feature_map(x(i), y(i), j) = Z1(i,j) + 1i*Z1(i+size(Z1, 1)/2,j);
        feature_map(x(i), y(i), j) = Z1(i,j);
    end
end

% feature_map =ifft(feature_map);
figure();
for slc = 1:size(feature_map, 3)
    subplot(2,2,slc);
    imagesc(feature_map(:,:,slc)); axis image;
    title(['Feature ', num2str(slc)])
    axis off;
    colorbar;
end
%% ica 2
q = 4;
mdl = rica(X2, q);
 
Z2 = transform(mdl, X2);
feature_map = zeros(size(img_t1_crop,1), size(img_t1_crop,2), size(Z2, 2));
for i = 1:size(Z2, 1)
    for j = 1:size(Z2,2)
        feature_map(x(i), y(i), j) = Z2(i,j);
    end
end

figure();
for slc = 1:size(feature_map, 3)
    subplot(2,2,slc);
    imagesc(feature_map(:,:,slc)); axis image;
    title(['Feature ', num2str(slc)])
    axis off;
    colorbar;
end
