
%% Initialization & and image recon (Unfinished)
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

[B0file, B0path] = uigetfile('*.dat');
MRdat_path=B0path;
case_name = '';
directory_name = '';

clear ls
if FlagLoadFolder
    ls=dir(B0path);
else
    ls(3).name=B0file;
end

%% IDEAL FF map
gap = 3;
N = (length(ls)-gap);

for n = 1:N
    filename = [MRdat_path, case_name, directory_name, ls(n+gap).name];
    twix_obj_in = mapVBVD(filename, 'removeOS', 'ignoreSeg'); % return all image-data:
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

    [NumRO, NumCh, NumPE, NumSlices, NumEchos] = size(rawdata);

    % clear twix_obj_in rawdata;

    % temp = strsplit(B0file{1}, '_');
    nameSession         = ['-' num2str(twix_obj.hdr.Config.MeasUID)];
    PatientName         = twix_obj.hdr.Meas.tPatientName; % Need to work on it
    temp = strsplit(B0path, '/');
    strings = strsplit(temp{end-2},'_');
    PatientName_new = cat(2, strings{end-2}, '_', strings{end-1}, '_', strings{end});
    old_path = GetFullPath(cat(2, B0path, '../../FF_Data_FromIDEAL/', PatientName, '/'));
    new_path = GetFullPath(cat(2, B0path, '../../FF_Data_FromIDEAL/', PatientName_new, '/'));
    [status,message,messageId] = movefile( old_path, new_path);
    
    RawdataFT = rawdata;
    RawdataFT = permute(RawdataFT, [1,3,4,2,5]);% NumRO, NumPE, NumSlices, NumCh, NumEcho

    % full_kspace = complex(zeros(settings.parameters.nRO, settings.parameters.nPE, 16, 12, 8, 'single'));
    full_kspace = complex(zeros(128, 96, 16, 12, 8, 'single'));
    odd_idx = 1:2:size(full_kspace, 3);
    even_idx = 2:2:size(full_kspace, 3);

%     full_kspace(settings.parameters.uRO+1:end,:,odd_idx,:,:)  = RawdataFT(:,-settings.parameters.uPE:end-1,odd_idx,:,:);
%     full_kspace(1:end-settings.parameters.uRO,:,even_idx,:,:) = RawdataFT(:,-settings.parameters.uPE:end-1,even_idx,:,:);

%     full_kspace(22+1:end,:,odd_idx,:,:)  = RawdataFT(:,-settings.parameters.uPE:end-1,odd_idx,:,:);
%     full_kspace(1:end-settings.parameters.uRO,:,even_idx,:,:) = RawdataFT(:,-settings.parameters.uPE:end-1,even_idx,:,:);

    full_kspace(22+1:end,:,odd_idx,:,:)  = RawdataFT(:,2:end-1,odd_idx,:,:);
    full_kspace(1:end-22,:,even_idx,:,:) = RawdataFT(:,2:end-1,even_idx,:,:);
    RawdataFT_5D = zeros(size(full_kspace));

    for j = 1:size(full_kspace,4)
        for k = 1:size(full_kspace,5)
            for f = 1:size(full_kspace,3)
                % xspace(:,:,f) = ifftshift(ifft2(kspace(:,:,f)));
                xspace(:,:,f) = ifftshift(ifftshift(ifftshift(ifft(ifft(ifft(fftshift(fftshift(fftshift(full_kspace(:,:,f,j,k),1),2),3),[],1),[],2),[],3),1),2),3);
            end
            RawdataFT_5D(:,:,:,j,k) = flip(flip((xspace),1),2);
        end
    end
    

    AllPhasemap(n).Name = ls(n+2).name(1:end-4);
    AllPhasemap(n).compleximg = RawdataFT;% In object domain
    AllPhasemap(n).Freq = twix_obj.hdr.Dicom.lFrequency;%Convert unit into Hz by * Central Frequency
    for necho=1:NumEchos
        AllPhasemap(n).TE(necho)=twix_obj.hdr.MeasYaps.alTE{necho};%usec
    end

    AllPhasemap(n).RoFOV = twix_obj.hdr.Config.RoFOV;%usec
    AllPhasemap(n).PeFOV = twix_obj.hdr.Config.PeFOV;%usec
    AllPhasemap(n).Voxelsize=[AllPhasemap(n).RoFOV/twix_obj.hdr.Meas.NImageCols, AllPhasemap(n).PeFOV/twix_obj.hdr.Meas.NImageLins ];

    temp = double(sum(abs(AllPhasemap(n).compleximg(:,:,:,:,end)),4));
    temp = temp/(max(temp(:)));

    % I have on idea why this isn't working
    % maskpercent=multithresh(temp/(max(temp(:))),2); doesn't work on my laptop
    % AllPhasemap(n).Mask=temp>(max(temp(:))*min(maskpercent));% generate a rough mask
    % XZ 06/06/2022, use region growing instead
    % sz_tmp = size(temp);
    % A = regiongrowing(temp, sz_tmp(1)/2, 1);
    % B = regiongrowing(temp, sz_tmp(1)/2, sz_tmp(2));
    %
    % se = strel('disk', 2);
    % mask_rev = A | B;
    % for s = 1:size(mask_rev, 3)
    %     mask_rev(:,:,s)= bwareaopen(mask_rev(:,:,s),100,4); % Also bwareaopen isn't working too, XZ 06/06/2022
    %     %mask_rev(:,:,s)= imopen(mask_rev(:,:,s), se);
    % end

    % load masks from cvi contours
%     if n == 1
%         strings = strsplit(PatientName_new, '_');
%         name = strings{1};
%         [mask_file, mask_folder] = uigetfile('*.mat', 'Select a mask file', GetFullPath(cat(2, pwd, '/../../T1_Fat_Project/ContourData/', name, '/')));
%         load(cat(2, mask_folder, mask_file));
%     end
%     
%     figure(); imagesc(mask_heart_3D(:,:,n) .* temp);
%     AllPhasemap(n).Mask = mask_heart_3D(:,:,n);

    switch Scan_Type
        case 'ICD'
            se = strel('disk',4);
            AllPhasemap(n).Mask=imclose( AllPhasemap(n).Mask,se);
        case 'Cardiac'
            %         se = strel('disk',4);
            %         AllPhasemap(n).Mask=imclose( AllPhasemap(n).Mask,se);
        otherwise
    end

    if IsManual_Mask == 1
        AllPhasemap(n).ManualMask = 1;
    else
        AllPhasemap(n).ManualMask = 0;
    end

    f_central = AllPhasemap(n).Freq; % MHz
    
    %RawdataFT = permute(RawdataFT, [4,1,2,5,3]);
    %RawdataFT_5D(:,:,:,:,n) = RawdataFT;
    clear B0map temp
end

    %% !!!!! resolution need double check
    [Dicomfile,Dicompath] = uigetfile('*.ima', 'Selcet IMA', B0path);  %Charles 30Jan2019, speed up finding files
    %[Dicomfile,Dicompath] = uigetfile('*.dcm', 'Selcet DCM', B0path);
    Dicomhdr=dicominfo([Dicompath,Dicomfile]);
    voxel_size=[Dicomhdr.PixelSpacing;Dicomhdr.SliceThickness];

    %% Apply Mask
clear iField
n = 1;

% switch Scan_Type
%     case 'Cardiac'
%         Mask=AllPhasemap(n).Mask;
%     case 'ICD'
%         %do nothing
%         Mask=AllPhasemap(n).Mask;
%     case 'Brain'
%         Apply_Mask_Brain
%     case 'fMRI'
%         Apply_Mask_Brain
%         Apply_Mask_Transverse
%     otherwise
%         Apply_Mask
% end

% XZ 06/08/2022
iField = permute(RawdataFT_5D,[1 2 3 5 4])*100; % NumRO, NumPE, NumSlices, NumCh, NumEcho =>  NumRO, NumPE, NumSlices, NumEcho, NumCh

Mask_pre = ones((size(iField,1)-2), (size(iField, 2)-2), size(iField, 3));
Mask = zeros((size(iField,1)), (size(iField, 2)), size(iField, 3));
Mask(2:end-1, 2:end-1, :) = Mask_pre;

% iMag = iMag .* Mask;
% iMag = iMag/max(iMag(:));% Notmalization

% Spatial phase unwrapping (region-growing)
% inputphase = unwrapPhase(iMag, iFreq_raw, matrix_size);
%------------------------SPURS-----------------------------------------


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

% figure(); imagesc(angle(iField(:,:,1,2)./iField(:,:,1,1))); caxis([-1 1]);
% figure(); imagesc(angle(iField(:,:,1,3)./iField(:,:,1,2))); caxis([-1 1]);
% figure(); imagesc(angle((iField(:,:,1,2)./iField(:,:,1,1))./(iField(:,:,1,3)./iField(:,:,1,2)))); caxis([-1 1]);
% figure(); imagesc(angle(iField(:,:,1,4)./iField(:,:,1,3))); caxis([-1 1]);
% figure(); imagesc(angle(iField(:,:,1,5)./iField(:,:,1,4))); caxis([-1 1]);
% figure(); imagesc(angle(iField(:,:,1,8)./iField(:,:,1,7))); caxis([-1 1]);


% odd_idx = [1,3,5,7];
% iFreq = unwrapPhase(iMag, iFreq_raw, matrix_size);
SUBSAMPLE = 1;
dfat = -3.5e-6*f_central;
LABEL = 2;

[water,fat,iFreq,unwph_uf,unwph,N_std,R2s] = ...
    spurs_gc(iField(:,:,:,:),AllPhasemap(n).TE*1e-6,f_central,voxel_size, Mask>0,SUBSAMPLE,dfat,LABEL); %iFreq in rad
% TODO (06/16/2022)
% Phase unwrapping is not working properly
% figure(); imagesc(angle(iField(:,:,1,2)./iField(:,:,1,1))); caxis([-1 1]);
% figure(); imagesc(angle(iField(:,:,1,3)./iField(:,:,1,2))); caxis([-1 1]);
% figure(); imagesc(angle((iField(:,:,1,2)./iField(:,:,1,1))./(iField(:,:,1,3)./iField(:,:,1,2)))); caxis([-1 1]);

figure;
subplot(2,2,1);
imagesc(abs(fat(:,:,12)));colorbar;colormap jet;caxis([0,0.5]);
subplot(2,2,2);
imagesc(abs(water(:,:,12)));colorbar;colormap jet;caxis([0,0.5]);
subplot(2,2,3);
imagesc(abs(iFreq(:,:,12)));colorbar;colormap jet;caxis([-pi pi])
subplot(2,2,4);
imagesc(unwph_uf(:,:,12));colorbar;colormap jet;


% normalization
%unwph = unwph.*Mask;
%unwph = unwph/max(unwph(:));
%unwph_uf = unwph_uf.*Mask;
%unwph_uf = unwph_uf/max(unwph_uf(:));
%figure;imagesc(unwph(:,:,15));colorbar;colormap jet;caxis([-0.5,0.5]);
%figure;imagesc(unwph_uf(:,:,15));colorbar;colormap jet;caxis([-0.5,0.5])
%----------------------------------------------------------------------

AllOtherMaps = struct;
%AllOtherMaps.inputphase = inputphase;
dTE_2echo=diff(AllPhasemap(n).TE)/1e6;
fmap_2echo = unwph/(2*pi*dTE_2echo(1)); % Hertz %using trimmed weighted mean
B0map = fmap_2echo/f_central; %ppm Central frequency

% AllOtherMaps.inputphase = inputphase;
AllOtherMaps.B0map = B0map;
AllOtherMaps.fmap_2echo = fmap_2echo;
AllOtherMaps.Voxelsize = voxel_size;
AllOtherMaps.water = water;
AllOtherMaps.fat = fat;
AllOtherMaps.unwph = unwph;
AllOtherMaps.R2s = R2s;
AllOtherMaps.AllPhasemap = AllPhasemap;

%% save file
if FlagSaveB0
    resdir=GetFullPath([MRdat_path, '..','/Result/']);
    if ~exist(resdir)
        mkdir(resdir);
    end
    save([resdir,'AllPhasemap.mat'], '-struct', 'AllOtherMaps','-v7.3');
end


