%Calculate B0map from Raw data with the directories
% Randy Yang 7/1/2019 Cedars Sinai

%% clear needed parameters and functions
clear AllPhasemap
FlagMask=1;
FlagSaveB0=1;
FlagLoadFolder=0;
IsManual_Mask=1;
SOS = @(x) sqrt(sum(abs(x).^2,3));  %sum of squareroots, used for magnitude calculation for all 12 recieving coils

%Choose file path
%pathForSearching_B0=pwd; 
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

    %% !!!!! resolution need double check
    %Load Dicom image info
    %[Dicomfile,Dicompath] = uigetfile('*.ima', 'Selcet IMA', B0path);  %Charles 30Jan2019, speed up finding files 
    pathForResolutionmatch = extractBefore(pwd, 'Code'); %Leon 20Aug2019, update direct folder
    pathForResolutionmatch = char ( join([ pathForResolutionmatch, 'Data\dicom'], '') ); %Leon 20Aug2019, update direct folder
    [Dicomfile,Dicompath] = uigetfile('*.ima', 'Selcet IMA', pathForResolutionmatch);  %Leon 20Aug2019, Above, edited
    Dicomhdr=dicominfo([Dicompath,Dicomfile]);
    voxel_size=[Dicomhdr.PixelSpacing;Dicomhdr.SliceThickness]; %Load the voxel size of Dicom images
%% Read files

clear phase_unwrapped phase_unwrapped_mn_3d phase_unwrapped_wt_3d B0mapstd;
for n=1:(length(ls)-2)  % loop used for field map batch process

    
    %%load Rawdata
    
    filename = [MRdat_path, case_name, directory_name, ls(n+2).name];
    %filename=[B0path,B0file];
    twix_obj_in = mapVBVD(filename); %Reads Siemens raw .dat file from VB/VD MRI raw data
    if (length(twix_obj_in)>1)% R.Y. avoid adj coil sensitivity. to determin whether input is the result of pre-scan or actual scan
        for  k=1:length(twix_obj_in)
            if (~strcmp(twix_obj_in{k}.hdr.MeasYaps.tSequenceFileName,'%AdjustSeq%/AdjCoilSensSeq') )
                twix_obj=twix_obj_in{k}
            end
        end
    else
        twix_obj=twix_obj_in
    end
    
    rawdata = squeeze(twix_obj.image(''));  %Rawdata from kspace loaded (only two echos)
    [NumRO, NumCh, NumPE, NumSlices, NumEchos] = size(rawdata);  %extract the number of x,y,z,coil channels,echos
    clear twix_obj_in

    %% R.Y. setup parameters (Create AllPhasemap(n))
    Rawdatapermute=permute(rawdata,[1,3,4,2,5]);  %x,y,z,coil channel,echo
    RawdataSLft=(ifft(Rawdatapermute,[],3));  %Z dirention ifft
    RawdataFT=fftshift(ifft2(RawdataSLft));  %x,y dirention ifft
    AllPhasemap(n).Name=B0file(1:end-4);  %extract the file name before '.dat'
    AllPhasemap(n).compleximg=RawdataFT;
    AllPhasemap(n).Freq=twix_obj.hdr.Dicom.lFrequency;%Hz
    for necho=1:NumEchos
        AllPhasemap(n).TE(necho)=twix_obj.hdr.MeasYaps.alTE{necho};%usec
    end
    AllPhasemap(n).RoFOV=twix_obj.hdr.Config.RoFOV;%mm
    AllPhasemap(n).PeFOV=twix_obj.hdr.Config.PeFOV;%mm
    [NumRO, NumPE, NumSlices, NumCh NumEchos] = size(AllPhasemap(n).compleximg);
    AllPhasemap(n).Voxelsize=[AllPhasemap(n).RoFOV/NumRO AllPhasemap(n).PeFOV/NumPE];%mm
    temp=double(sum(abs(AllPhasemap(n).compleximg(:,:,:,:,end)),4));
     %maskpercent=multithresh(temp/(max(temp(:))),2);
     AllPhasemap(n).Mask=autoMask(temp,AllPhasemap(n).Voxelsize); %Otsu's method, mask out the air area

    if ~isfield(AllPhasemap(n),'ManualMask')
        AllPhasemap(n).ManualMask=0;  %whether ManualMask existed already
    end
    %AllPhasemap(n).Phasediff=AllPhasemap(n).Phasediff.*AllPhasemap(n).Mask;
    f_central = AllPhasemap(n).Freq; % MHz
    clear B0map temp

    %%  Creat B0map
    
    clear iField
    
    iField=permute(RawdataFT,[1 2 3 5 4])*100; %x,y,z,echo,coil channel
    %Universal B0Map from SPURS
    % calculate the nmap first and apply the shimming mask afterwards
    if size(iField,5)>1 %coil combination
% combine multiple coils together, assuming the coil is the fifth dimension
    iField = sum(iField.*conj( repmat(iField(:,:,:,size(iField,4)/2,:),[1 1 1 size(iField,4) 1])),5);  %replicate the iField in 'echo' direction,calcualte iField(e1)*iField(e1)_conj(phaseof result=0) and iField(e2)*iField(e1)_conj(phase of result=delta phase in e1 and e2) 
    iField = sqrt(abs(iField)).*exp(1i*angle(iField));
    end
    
    [unwph_uf unwph N_std] = spurs_gc_UNIC(iField,AllPhasemap(n).TE/1e6,AllPhasemap(n).Freq,voxel_size);
    inputphase=unwph_uf;% unwrap and fat corrected phase map, its the map of phase difference between echo1 and echo2, not field map nor phase image at certain echo
   
    
    %% Apply Mask
    switch Scan_Type
        case 'Cardiac'
            Apply_Mask_SPUR
        case 'ICD'
            Apply_Mask_SPUR
        case 'Brain'
            Apply_Mask_Brain
        case 'fMRI'
            Apply_Mask_Brain
            Apply_Mask_SPUR
        otherwise
            Apply_Mask_SPUR
    end
    iMag=abs(iField(:,:,:,1));
    iMag=iMag.*(Mask);
    iMag=iMag/max(iMag(:));

    AllPhasemap(n).phase_unwrapped_3d=inputphase;
    dTE_2echo=diff(AllPhasemap(n).TE)/1e6;
    fmap_2echo = inputphase/(2*pi*dTE_2echo(1)); % Hertz %using trimmed weighted mean &field map eqn
    
    B0map = fmap_2echo/f_central; %ppm Central frequency
    AllPhasemap(n).B0map=B0map;%ppm
    AllPhasemap(n).fmap_2echo=fmap_2echo;
    AllPhasemap(n).Voxelsize=voxel_size;
    
    
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
    save([resdir,'UNIC_B0Map',B0file(end-12:end-4),'.mat'],'UNIC_B0Map');
end
