%Calculate B0map from Raw data with the directories
% Randy Yang 5/10/2018
%clear all
clear AllPhasemap
FlagMask=1;
FlagSaveB0=1
FlagLoadFolder=0;
IsManual_Mask=1;

maskpercent=0.1;
SOS = @(x) sqrt(sum(abs(x).^2,3));
[B0file,B0path] = uigetfile('*.dat');
MRdat_path=B0path;
case_name = '';
directory_name = '';

clear ls
if FlagLoadFolder
    ls=dir(B0path);
else
    ls(3).name=B0file;
end

%addpath('U:\Randy\Unic\Invivo studies\Code\MEIDfunctions')
%%

clear phase_unwrapped phase_unwrapped_mn_3d phase_unwrapped_wt_3d B0mapstd;
for n=1:(length(ls)-2)
    
    %%load Rawdata
    
    filename = [MRdat_path, case_name, directory_name, ls(n+2).name];
    %filename=[B0path,B0file];
    twix_obj_in = mapVBVD(filename);
    if (length(twix_obj_in)>1)% R.Y. avoid adj coil sensitivity
        for  k=1:length(twix_obj_in)
            if (~strcmp(twix_obj_in{k}.hdr.MeasYaps.tSequenceFileName,'%AdjustSeq%/AdjCoilSensSeq') )
                twix_obj=twix_obj_in{k}
            end
        end
    else
        twix_obj=twix_obj_in
    end
    
    rawdata = squeeze(twix_obj.image(''));
    [NumRO, NumCh, NumPE, NumSlices, NumEchos] = size(rawdata);
    clear twix_obj_in
    %% R.Y. setup parameters
    Rawdatapermute=permute(rawdata,[1,3,4,2,5]);
    RawdataSLft=(ifft(Rawdatapermute,[],3));
    RawdataFT=fftshift(ifft2(RawdataSLft));
    AllPhasemap(n).Name=B0file(1:end-4);
    %AllPhasemap(n).Shimch=filename(end-12:end-7);
    %AllPhasemap(n).Current=str2num(filename(end-5));
    
    AllPhasemap(n).compleximg=RawdataFT;
    AllPhasemap(n).Freq=twix_obj.hdr.Dicom.lFrequency;%Hz
    AllPhasemap(n).TE(1)=twix_obj.hdr.MeasYaps.alTE{1};%usec
    AllPhasemap(n).TE(2)=twix_obj.hdr.MeasYaps.alTE{2};%usec
    AllPhasemap(n).RoFOV=twix_obj.hdr.Config.RoFOV;%usec
    AllPhasemap(n).PeFOV=twix_obj.hdr.Config.PeFOV;%usec
    
    %AllPhasemap(n).compleximg=permute(AllPhasemap(n).compleximg,[1 3 2 4 5]);
    [NumRO, NumPE, NumSlices, NumCh NumEchos] = size(AllPhasemap(n).compleximg);
    AllPhasemap(n).Voxelsize=[AllPhasemap(n).RoFOV/NumRO AllPhasemap(n).PeFOV/NumPE ];%usec
    %AllPhasemap(n).Phasediff=angle(AllPhasemap(n).compleximg(:,:,:,:,2)./AllPhasemap(n).compleximg(:,:,:,:,1));
    temp=sum(abs(AllPhasemap(n).compleximg(:,:,:,:,1)),4);
    AllPhasemap(n).Mask=temp>(max(temp(:))*maskpercent);%generate a rough mask
    for s=1:size(AllPhasemap(n).Mask,3)
        AllPhasemap(n).Mask(:,:,s)= bwareaopen(AllPhasemap(n).Mask(:,:,s),100,4);
    end
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
    [iFreq_raw N_std] = Fit_ppm_complex(iField);
    % Compute magnitude image
    iMag = sqrt(sum(abs(iField(:,:,:,1,:)).^2,5));
    matrix_size=size(iMag);

    %% !!!!! resolution need double check
    [Dicomfile,Dicompath] = uigetfile('*.ima');
    Dicomhdr=dicominfo([Dicompath,Dicomfile]);
    voxel_size=[Dicomhdr.PixelSpacing;Dicomhdr.SliceThickness];
    
    %% Apply Mask
    switch Scan_Type
        case 'Cardiac'
            Apply_Mask
        case 'Brain'
            Apply_Mask_Brain
        otherwise
            Apply_Mask
    end
        iMag=iMag.*(Mask);

    % Spatial phase unwrapping (region-growing)
    inputphase = unwrapPhase(iMag, iFreq_raw, matrix_size);
    
    AllPhasemap(n).phase_unwrapped_3d=inputphase;
    dTE_2echo=diff(AllPhasemap(n).TE)/1e6;
    fmap_2echo = inputphase/(2*pi*dTE_2echo); % Hertz %using trimmed weighted mean
    
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
