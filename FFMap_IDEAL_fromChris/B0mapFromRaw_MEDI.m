%Calculate B0map from Raw data with the directories
% Randy Yang 5/10/2018
% Modified by Chris Huang 6/29/2019
%% Set up data path
clear all
clear AllPhasemap
%%
addpath(genpath(pwd));
FlagMask=1;
FlagSaveB0=1;
FlagLoadFolder=0;
IsManual_Mask=1;
Scan_Type='ICD';
%maskpercent=0.1;

SOS = @(x) sqrt(sum(abs(x).^2,3));

pathForSearching_B0=pwd;
pathForSearching_B0 = getFilePathFor_B0(); %Charles 30Jan2019, speed up finding files
[B0file,B0path] = uigetfile('*.dat', 'Select .dat', pathForSearching_B0);% Above, edited
MRdat_path=B0path;
case_name = '';
directory_name = '';

clear ls
if FlagLoadFolder
    ls=dir(B0path);
else
    ls(3).name=B0file;
end

%% XZ
[fid_file, fid_path] = uigetfile('*.dat');
MRdat_path=fid_path;
case_name = '';
directory_name = '';
B0file = fid_file;
B0path = fid_path;

clear ls
if FlagLoadFolder
    ls=dir(B0path);
else
    ls(3).name=B0file;
end
%addpath('U:\Randy\Unic\Invivo studies\Code\MEIDfunctions')
%% Load data and set up parameters

clear phase_unwrapped phase_unwrapped_mn_3d phase_unwrapped_wt_3d B0mapstd;
for n=1:(length(ls)-2)
    
    %%load Rawdata
    
    filename = [MRdat_path, case_name, directory_name, ls(n+2).name];
    %filename=[B0path,B0file];
    twix_obj_in = mapVBVD(filename); % return all image-data:
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
        rawdata = reshape(rawdata, [NumRO, NumCh, NumPE, NumSlices, NumEchos, NSeg]);
        Lin = twix_obj.image.Lin;
        Eco = twix_obj.image.Eco;
        Seg = twix_obj.image.Seg;
    end
    clear twix_obj_in
    %% R.Y. setup parameters
    % How to recon with iPAT  = 2 and NSeg ~= 1;
    if length(size(rawdata)) == 6
        Rawdatapermute = permute(rawdata,[1,3,4,2,5,6]);
    end
    Rawdatapermute=permute(rawdata,[1,3,4,2,5]); % NumRO, NumPE, NumSlices, NumCh, NumEcho
    RawdataSLft=(ifft(Rawdatapermute,[],3));
    RawdataFT=fftshift(ifft2(RawdataSLft));
    AllPhasemap(n).Name=B0file(1:end-4);
    %AllPhasemap(n).Shimch=filename(end-12:end-7);
    %AllPhasemap(n).Current=str2num(filename(end-5));
    
    AllPhasemap(n).compleximg=RawdataFT;% In object domain
    AllPhasemap(n).Freq=twix_obj.hdr.Dicom.lFrequency;%Convert unit into Hz by * Central Frequency
    for necho=1:NumEchos
    AllPhasemap(n).TE(necho)=twix_obj.hdr.MeasYaps.alTE{necho};%usec
    end
%     AllPhasemap(n).TE(2)=twix_obj.hdr.MeasYaps.alTE{2};%usec   
%     AllPhasemap(n).TE(3)=twix_obj.hdr.MeasYaps.alTE{3};%usec
%     AllPhasemap(n).TE(4)=twix_obj.hdr.MeasYaps.alTE{4};%usec
%     AllPhasemap(n).TE(5)=twix_obj.hdr.MeasYaps.alTE{5};%usec
%     AllPhasemap(n).TE(6)=twix_obj.hdr.MeasYaps.alTE{6};%usec
    AllPhasemap(n).RoFOV=twix_obj.hdr.Config.RoFOV;%usec
    AllPhasemap(n).PeFOV=twix_obj.hdr.Config.PeFOV;%usec
    
    %AllPhasemap(n).compleximg=permute(AllPhasemap(n).compleximg,[1 3 2 4 5]);
    [NumRO, NumPE, NumSlices, NumCh NumEchos] = size(AllPhasemap(n).compleximg);
    AllPhasemap(n).Voxelsize=[AllPhasemap(n).RoFOV/NumRO AllPhasemap(n).PeFOV/NumPE ];%usec
    %AllPhasemap(n).Phasediff=angle(AllPhasemap(n).compleximg(:,:,:,:,2)./AllPhasemap(n).compleximg(:,:,:,:,1));
    temp=double(sum(abs(AllPhasemap(n).compleximg(:,:,:,:,end)),4));
   % maskpercent=graythresh(temp/(max(temp(:))));
   maskpercent=multithresh(temp/(max(temp(:))),2);
    AllPhasemap(n).Mask=temp>(max(temp(:))*min(maskpercent));%generate a rough mask
    for s=1:size(AllPhasemap(n).Mask,3)
        AllPhasemap(n).Mask(:,:,s)= bwareaopen(AllPhasemap(n).Mask(:,:,s),100,4);
    end
      switch Scan_Type 
        case 'ICD'
      se = strel('disk',4);
     AllPhasemap(n).Mask=imclose( AllPhasemap(n).Mask,se);
      otherwise
      end
%    end
    if ~isfield(AllPhasemap(n),'ManualMask')
        AllPhasemap(n).ManualMask=0;
    end
    %AllPhasemap(n).Phasediff=AllPhasemap(n).Phasediff.*AllPhasemap(n).Mask;
    f_central = AllPhasemap(n).Freq; % MHz
    clear B0map temp
    
    
    
    
    %%  Creat B0map
    
    clear iField
    %[iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir]=Read_DICOM('DICOM_dir');
    iField=permute(RawdataFT,[1 2 3 5 4])*100;
    % Estimate the frequency offset in each of the voxel using a complex
    % fitting (even echo spacing)
    [iFreq_raw,N_std] = Fit_ppm_complex(iField);% Simple spicies fitting?

    % Compute magnitude image
    iMag = sqrt(sum(abs(iField(:,:,:,end,:)).^2,5)); 
    matrix_size=size(iMag);

    %% !!!!! resolution need double check
    %[Dicomfile,Dicompath] = uigetfile('*.ima', 'Selcet IMA', B0path);  %Charles 30Jan2019, speed up finding files
    [Dicomfile,Dicompath] = uigetfile('*.dcm', 'Selcet DCM', B0path);
    Dicomhdr=dicominfo([Dicompath,Dicomfile]);
    voxel_size=[Dicomhdr.PixelSpacing;Dicomhdr.SliceThickness];
    
    %% Apply Mask
    switch Scan_Type
        case 'Cardiac'
            Apply_Mask
        case 'ICD'
            %do nothing
            Mask=AllPhasemap(n).Mask;
        case 'Brain'
            Apply_Mask_Brain
        case 'fMRI'
            Apply_Mask_Brain
            Apply_Mask_Transverse
        otherwise
            Apply_Mask
    end

    iMag=iMag.*(Mask);
   iMag=iMag/max(iMag(:));% Notmalization

    % Spatial phase unwrapping (region-growing)
    inputphase = unwrapPhase(iMag, iFreq_raw, matrix_size);
    %------------------------SPURS-----------------------------------------

    if size(iField,5)>1
% combine multiple coils together, assuming the coil is the fifth dimension
    iField = sum(iField.*conj( repmat(iField(:,:,:,1,:),[1 1 1 size(iField,4) 1])),5);  % 
    iField = sqrt(abs(iField)).*exp(1i*angle(iField));
    end
    % iFreq = unwrapPhase(iMag, iFreq_raw, matrix_size);
[water,fat,iFreq,unwph_uf,unwph,N_std] = ...
    spurs_gc(iField,AllPhasemap(n).TE,f_central,voxel_size); %iFreq in rad
    % normalization
    %unwph = unwph.*Mask;
    %unwph = unwph/max(unwph(:));
    %unwph_uf = unwph_uf.*Mask;
    %unwph_uf = unwph_uf/max(unwph_uf(:));
    %figure;imagesc(unwph(:,:,15));colorbar;colormap jet;caxis([-0.5,0.5]);
    %figure;imagesc(unwph_uf(:,:,15));colorbar;colormap jet;caxis([-0.5,0.5])
    %----------------------------------------------------------------------
    %%
 switch Scan_Type
     case 'ICD'
%          for sl=1:size(iMag,3)
%          IM3D=iMag.*exp(i*iFreq_raw);
%          IM=IM3D(:,:,sl);
%          im_mask=Mask(:,:,sl);
%         GoldsteinUnwrap2D_r1;
%          inputphase(:,:,sl)=im_unwrapped;
%          end
    %   inputphase=  iFreq_raw;
     Apply_Mask
     otherwise
 end
    AllPhasemap(n).phase_unwrapped_3d=inputphase;
    dTE_2echo=diff(AllPhasemap(n).TE)/1e6;
    fmap_2echo = inputphase/(2*pi*dTE_2echo(1)); % Hertz %using trimmed weighted mean
    
    B0map = fmap_2echo/f_central; %ppm Central frequency
    AllPhasemap(n).B0map=B0map;%ppm
    AllPhasemap(n).fmap_2echo=fmap_2echo;
    AllPhasemap(n).Voxelsize=voxel_size;
    

    %%
    
    n
end

%%
clear UNIC_FreqMap UNIC_B0Map;

%parameters

UNIC_B0Map.Parameters.Rawdatapath=MRdat_path;
UNIC_B0Map.Parameters.Mask=AllPhasemap(1).Mask;
UNIC_B0Map.Parameters.CentralFeq=AllPhasemap(1).Freq;%Hz
UNIC_B0Map.Parameters.TEs=AllPhasemap(1).TE;%usec
UNIC_B0Map.Parameters.Voxelsize=voxel_size;
UNIC_B0Map.Parameters.CornerCoordinate=Dicomhdr.ImagePositionPatient;
UNIC_B0Map.Parameters.PhaseDir=Dicomhdr.InPlanePhaseEncodingDirection;
UNIC_B0Map.Parameters.MatrixSize=size(iMag);

UNIC_FreqMap.Parameters=UNIC_B0Map.Parameters;
%UNIC COIL(each channel minus baseline to only show the shim effect)
for n=1:length(AllPhasemap)
    Phasemap_coil=AllPhasemap(n).fmap_2echo;%Hz
    eval(['UNIC_FreqMap.Map=Phasemap_coil;']);
end
for n=1:length(AllPhasemap)
    Phasemap_coil=AllPhasemap(n).B0map;%ppm
    eval(['UNIC_B0Map.Map=Phasemap_coil;']);
end

%% save file
if FlagSaveB0
    resdir=[MRdat_path,AllPhasemap(n).Name,'Result\',];
    mkdir(resdir);
    %save([resdir,'AllPhasemap_b0.mat'],'AllPhasemap','-v7.3');
    
    %save([resdir,'UNIC_FreqMap',B0file(end-12:end-4),'.mat'],'UNIC_FreqMap');
    save([resdir,'UNIC_B0Map',B0file(end-12:end-4),'.mat'],'UNIC_B0Map');
end
