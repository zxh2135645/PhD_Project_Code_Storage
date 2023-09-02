%Calculate B0map from Raw data with the directories
% Randy Yang 5/10/2018
% Modified by Chris Huang 6/29/2019
%% Set up data path
clear all
clear AllPhasemap
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
    [NumRO, NumCh, NumPE, NumSlices, NumEchos] = size(rawdata);
    clear twix_obj_in
    %% R.Y. setup parameters
    Rawdatapermute=permute(rawdata,[1,3,4,2,5]);
    % ifft
    for i = 1:NumEchos
        for j = 1:NumCh
            RawdataFT(:,:,:,j,i)=ifftshift(ifftn(squeeze(Rawdatapermute(:,:,:,j,i))));
        end
    end
    
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
    % 5D datat (x,y,z, NTe, Ncoil);
    iField=permute(RawdataFT,[1 2 3 5 4])*100;
    
    % Compute magnitude image
    iMag = sqrt(sum(abs(iField(:,:,:,end,:)).^2,5));
    matrix_size=size(iMag);
    
    voxel_size(1) = AllPhasemap.Voxelsize(1);
    voxel_size(2) = AllPhasemap.Voxelsize(1);
    voxel_size(3) = AllPhasemap.Voxelsize(2);
    Mask=autoMask(iMag,voxel_size);
    
    iMag=iMag.*(Mask);
    iMag=iMag/max(iMag(:));% Normalization
    
    %------------------------Coil Combination-----------------------------------------
    % 4D datat (x,y,z, NTe)
    if size(iField,5)>1
        % combine multiple coils together, assuming the coil is the fifth dimension
        iField_Coilcombined = sum(iField.*conj( repmat(iField(:,:,:,1,:),[1 1 1 size(iField,4) 1])),5);
        iField_Coilcombined = sqrt(abs(iField_Coilcombined)).*exp(1i*angle(iField_Coilcombined));
    end
    % Compute magnitude image
    iMag_Coilcombined = sqrt(sum(abs(iField_Coilcombined(:,:,:,end,:)).^2,5));
    matrix_size_Coilcombined=size(iMag_Coilcombined);
    
    Mask_Coilcombined=autoMask(iMag_Coilcombined,voxel_size);
    
    iMag_Coilcombined=iMag_Coilcombined.*(Mask_Coilcombined);
    iMag_Coilcombined=iMag_Coilcombined/max(iMag_Coilcombined(:));% Normalization
    
    
end
