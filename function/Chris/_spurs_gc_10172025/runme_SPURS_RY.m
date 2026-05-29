%clear all;clc;close all;
% STEP 1: Import data
%[iField,voxel_size,matrix_size,CF,delta_TE,TE,B0_dir]=...
%Read_GE_DICOM('AXL_QSM');
% load B0mapSpurs.mat
% delta_TE=B0mapSpurs.delta_TE;
% TE=B0mapSpurs.TE;
% CF=B0mapSpurs.CF;
% iField=B0mapSpurs.iField;
% voxel_size=B0mapSpurs.voxel_size;
% matrix_size=B0mapSpurs.matrix_size;
% Mask=B0mapSpurs.Mask;
Mask_exist=0;
Scan_Type='ICD';%'Cardiac'%'Brain';%'fMRI';%'Brain';%'Cardiac'
IS_Slice=0;



%% Load B0map data & Calculate B0map

%B0mapFromRaw_MEDI_fat
B0mapFromRaw_MEDI
MaskWhole=Mask;
%% 


if size(iField,5)>1
% combine multiple coils together, assuming the coil is the fifth dimension
    iField = sum(iField.*conj( repmat(iField(:,:,:,size(iField,4)/2,:),[1 1 1 size(iField,4) 1])),5);  
    iField = sqrt(abs(iField)).*exp(1i*angle(iField));
end
iField=iField/abs(max(iField(Mask)));
% %%%%% provide a Mask here if possible %%%%%%
% if (~exist('Mask','var'))                     
%     Mask = genMask(iField, voxel_size);
% end
% 
% %%%%% provide a noise_level here if possible %%%%%%
% if (~exist('noise_level','var'))
%     noise_level = calfieldnoise(iField, Mask)
% end
% 
% %%%%% normalize signal intensity by noise to get SNR %%%
% iField = iField/noise_level;
% 
% %%%% Generate the Magnitude image %%%%
% iMag = sqrt(sum(abs(iField).^2,4));

% STEP 2a: Field Map Estimation
%%%%%Estimate the frequency offset in each of the voxel using a 
%%%%%complex fitting %%%%
% if abs((TE(2)-TE(1))-(TE(3)-TE(2)))< 0.0002
%     [iFreq_raw N_std] = Fit_ppm_complex(iField);
% else
%     [iFreq_raw N_std] = Fit_ppm_complex_TE(iField);
% end
   
% STEP 2b: Spatial phase unwrapping %%%%
% iFreq = unwrapPhase(iMag, iFreq_raw, matrix_size);
[water fat iFreq unwph_uf unwph N_std] = spurs_gc_UNIC(iField(:,:,:,:,:),AllPhasemap(n).TE(:)/1e6,f_central,voxel_size);
%*** For high-resolution data with large matrix size, may try subsampling ***% 
% [water fat iFreq unwph_uf unwph N_std] = spurs_gc(iField,TE,CF,voxel_size,2);


% STEP 2c: Background Field Removal
%%%% Background field removal %%%%
%[RDF shim] = PDF(iFreq, N_std, Mask,matrix_size,voxel_size, B0_dir);


% STEP 3: Dipole Inversion
%save RDF.mat RDF iFreq iFreq_raw iMag N_std Mask matrix_size...
%     voxel_size delta_TE CF B0_dir;
%%%% run MEDI %%%%%
%QSM = MEDI_L1('lambda',1000);

