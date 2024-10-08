clear all;
close all;
%% Initialization & and image recon
addpath(genpath(pwd));
rmpath(genpath('./Example/'));
addpath(genpath('../reconstructionPipeline/'));
addpath('../function/');
FlagMask=1;
FlagSaveB0=1;
FlagLoadFolder=1;
IsManual_Mask=1;
Scan_Type='Cardiac';
%maskpercent=0.1;


SOS = @(x) sqrt(sum(abs(x).^2,3));

dicom_dir = uigetdir();
name_glob = glob(cat(2, dicom_dir, '/*'));
dicom_folder = '3D_IDEAL_Avg4_Bipolar';
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
    };

for nn = 1:length(name_glob)

strings = strsplit(name_glob{nn}, '/');
name = strings{end};

[img_mag slice_data_mag] = dicom23D(cat(2, name_glob{nn}, '/', dicom_folder, '/'), dicom_fields);
[img_phase slice_data_phase] = dicom23D(cat(2, name_glob{nn}, '/', dicom_folder, '_PHASE/'), dicom_fields);

matrix_size(1) = single(slice_data_mag(1).Width);
matrix_size(2) = single(slice_data_mag(1).Height);
matrix_size(3) = 1;

eco = slice_data_mag(1).EchoTime;
eco_array = [];
eco_array(1) = eco;
n = 2;

while eco ~= slice_data_mag(n).EchoTime
    eco_array(n) = slice_data_mag(n).EchoTime;
    n = n + 1;
end

num_eco = length(eco_array);

for i = 1:length(slice_data_mag)
    if mod(i, num_eco) == 0
        slice_data_mag(i).EchoNumber = num_eco;
        slice_data_phase(i).EchoNumber = num_eco;
    else
        slice_data_mag(i).EchoNumber = mod(i, num_eco);
        slice_data_phase(i).EchoNumber = mod(i, num_eco);
    end
end


voxel_size(1,1) = single(slice_data_mag(1).PixelSpacing(1));
voxel_size(2,1) = single(slice_data_mag(1).PixelSpacing(2));
voxel_size(3,1) = single(slice_data_mag(1).SliceThickness);
CF = slice_data_mag(1).ImagingFrequency *1e6;

minSlice = 1e10;
maxSlice = -1e10;
rctr = 0; ictr=0;
progress = '';
TE = eco_array;
NumEcho = num_eco;

p_ImagePositionPatient = zeros(3, matrix_size(3)*NumEcho);
m_ImagePositionPatient = zeros(3, matrix_size(3)*NumEcho);

for i = 1:length(slice_data_mag)
    if slice_data_mag(i).SliceLocation<minSlice
        minSlice = slice_data_mag(i).SliceLocation;
        minLoc = slice_data_mag(i).ImagePositionPatient;
    end
    if slice_data_mag(i).SliceLocation>maxSlice
        maxSlice = slice_data_mag(i).SliceLocation;
        maxLoc = slice_data_mag(i).ImagePositionPatient;
    end
    if (slice_data_phase(i).ImageType(18)=='P')||(slice_data_phase(i).ImageType(18)=='p')
        rctr = rctr + 1;
        for ii=1:length(progress); fprintf('\b'); end
        progress=sprintf('Reading file %d', rctr+ictr);
        fprintf(progress);
        p_ImagePositionPatient(:,rctr) = slice_data_phase(i).ImagePositionPatient;
        r_EchoNumber(rctr) = slice_data_phase(i).EchoNumber;
        ph = transpose(single(img_phase(:,:,rctr)));
        iPhase(:,:,rctr)  = (ph*slice_data_phase(i).RescaleSlope + slice_data_phase(i).RescaleIntercept)/single(max(ph(:)))*pi;%phase
    end
    if (slice_data_mag(i).ImageType(18)=='M')||(slice_data_mag(i).ImageType(18)=='m')
        ictr = ictr + 1;
        for ii=1:length(progress); fprintf('\b'); end
        progress=sprintf('Reading file %d', rctr+ictr);
        fprintf(progress);
        m_ImagePositionPatient(:,ictr) = slice_data_mag(i).ImagePositionPatient;
        i_EchoNumber(ictr) = slice_data_mag(i).EchoNumber;
        iMag(:,:,ictr)  = transpose(single(img_mag(:,:,ictr)));%magnitude
    end
end
fprintf('\n');

matrix_size(3) = round(norm(maxLoc - minLoc)/voxel_size(3))+1 ;

Affine2D = reshape(slice_data_mag(1).ImageOrientationPatient,[3 2]);
Affine3D = [Affine2D (maxLoc-minLoc)/( (matrix_size(3)-1)*voxel_size(3))];
B0_dir = Affine3D\[0 0 1]';

sz=size(p_ImagePositionPatient); sz(1)=1;
minLoc=repmat(minLoc, sz);

p_slice = int32(round(sqrt(sum((p_ImagePositionPatient-minLoc).^2,1))/voxel_size(3)) +1);
p_ind = sub2ind([matrix_size(3) NumEcho], p_slice(:), int32(r_EchoNumber(:)));
iPhase(:,:,p_ind) = iPhase;

m_slice = int32(round(sqrt(sum((m_ImagePositionPatient-minLoc).^2,1))/voxel_size(3)) +1);
m_ind = sub2ind([matrix_size(3) NumEcho], m_slice(:), int32(i_EchoNumber(:)));
iMag(:,:,m_ind) = iMag;

iField = reshape(iMag.*exp(-1i*iPhase), ...
    [matrix_size(1) matrix_size(2) matrix_size(3) NumEcho]);
% iField = permute(iField,[2 1 3 4 5]); %This is because the first dimension is row in DICOM but COLUMN in MATLAB
% iField(:,:,1:2:end,:) = -iField(:,:,1:2:end,:);
if length(TE)==1
    delta_TE = TE;
else
    delta_TE = TE(2) - TE(1);
end


iField_trunc = iField(:,:,49:52,:);
% figure(); imagesc(angle(iField_trunc(:,:,1,2)./iField_trunc(:,:,1,1))); caxis([-1 1]);
% figure(); imagesc(angle(iField_trunc(:,:,1,3)./iField_trunc(:,:,1,2))); caxis([-1 1]);
% figure(); imagesc(angle(iField_trunc(:,:,1,4)./iField_trunc(:,:,1,3))); caxis([-1 1]);
% figure(); imagesc(angle(iField_trunc(:,:,1,5)./iField_trunc(:,:,1,4))); caxis([-1 1]);
% figure(); imagesc(angle((iField_trunc(:,:,1,2)./iField_trunc(:,:,1,1))./(iField_trunc(:,:,1,3)./iField_trunc(:,:,1,2)))); caxis([-1 1]);

f_central = CF;
Mask_pre = ones((size(iField_trunc,1)-2), (size(iField_trunc, 2)-2), size(iField_trunc, 3));
Mask = zeros((size(iField_trunc,1)), (size(iField_trunc, 2)), size(iField_trunc, 3));
Mask(2:end-1, 2:end-1, :) = Mask_pre;

% correct for eddy current
mag_iField = abs(iField_trunc);
mag_iField_norm = mag_iField ./ max(mag_iField(:));

Fieldmap_eddy  = sqrt(sqrt((iField_trunc(:,:,:,2)./iField_trunc(:,:,:,1))./(iField_trunc(:,:,:,3)./iField_trunc(:,:,:,2))));
iField_uneddy = zeros(size(iField_trunc));
iField_uneddy(:,:,:,1:2:end) = iField_trunc(:,:,:,1:2:end).*Fieldmap_eddy;
iField_uneddy(:,:,:,2:2:end) = iField_trunc(:,:,:,2:2:end)./Fieldmap_eddy;
iField_trunc = mag_iField_norm.*exp(1i*angle(iField_uneddy));

figure(); imagesc(angle(iField_trunc(:,:,1,2)./iField_trunc(:,:,1,1))); caxis([-1 1]);
figure(); imagesc(angle(iField_trunc(:,:,1,3)./iField_trunc(:,:,1,2))); caxis([-1 1]);
figure(); imagesc(angle(iField_trunc(:,:,1,4)./iField_trunc(:,:,1,3))); caxis([-1 1]);
figure(); imagesc(angle(iField_trunc(:,:,1,5)./iField_trunc(:,:,1,4))); caxis([-1 1]);
figure(); imagesc(angle((iField_trunc(:,:,1,2)./iField_trunc(:,:,1,1))./(iField_trunc(:,:,1,3)./iField_trunc(:,:,1,2)))); caxis([-1 1]);


[water,fat,iFreq,unwph_uf,unwph,N_std] = ...
    spurs_gc(iField_trunc(:,:,:,:),TE*1e-3,f_central,voxel_size, Mask>0);

%% Save data
figure;
subplot(2,2,1);
imagesc(abs(fat(:,:,1)));colorbar;colormap jet;caxis([0,0.2]);
subplot(2,2,2);
imagesc(abs(water(:,:,1)));colorbar;colormap jet;caxis([0,0.5]);
subplot(2,2,3);
imagesc(abs(iFreq(:,:,1)));colorbar;colormap jet;caxis([-pi pi])
subplot(2,2,4);
imagesc(abs(unwph_uf(:,:,1)));colorbar;colormap jet;

AllOtherMaps = struct;
%AllOtherMaps.inputphase = inputphase;
dTE_2echo=diff(TE)/1e-3;
fmap_2echo = unwph/(2*pi*dTE_2echo(1)); % Hertz %using trimmed weighted mean
B0map = fmap_2echo/f_central; %ppm Central frequency

% AllOtherMaps.inputphase = inputphase;
AllOtherMaps.B0map = B0map;
AllOtherMaps.fmap_2echo = fmap_2echo;
AllOtherMaps.Voxelsize = voxel_size;
AllOtherMaps.water = water;
AllOtherMaps.fat = fat;
AllOtherMaps.unwph = unwph;

if FlagSaveB0
    resdir=GetFullPath([dicom_dir, '/..','/Result/']);
    if ~exist(resdir)
        mkdir(resdir);
    end
    strings = strsplit(dicom_dir, '/');
    subject_name = strings{end};
    save([resdir,'AllPhasemap_', subject_name, '.mat'], '-struct', 'AllOtherMaps','-v7.3');
end

end
