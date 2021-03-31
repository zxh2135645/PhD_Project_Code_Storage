[fid_file, fid_path] = uigetfile('*.dat');
dim = 2;

ScanTypes={'Cine', 'IR', 'T2IR', 'SR', 'T2prep'};
ScanType = ScanTypes{listdlg('PromptString','What type of scan?', ...
  'SelectionMode','single',...
  'ListString',ScanTypes,...
  'InitialValue',2,'ListSize',[160 80])};

file_string = [fid_path fid_file(1:(end-4)) '_' datestr(now,30)];
mkdir(file_string);
file_string = [file_string sep 'AC_recon'];
diary([file_string '.log']);

%%
fprintf('Loading raw data...');
addpath('/home/christodoula/Dropbox/Cedars/Research/Code/MMR/mapVBVD')
twix_obj_in = mapVBVD(strcat(fid_path,fid_file));
if (length(twix_obj_in)>1)% R.Y. avoid adj coil sensitivity
for n=1:length(twix_obj_in)
    if (~strcmp(twix_obj_in{n}.hdr.MeasYaps.tSequenceFileName,'%AdjustSeq%/AdjCoilSensSeq') )
    twix_obj=twix_obj_in{n}
    end
end
else
    twix_obj=twix_obj_in
end
clear twix_obj_in
% R.Y. 
kspace_data=twix_obj.image.unsorted(); %NCol, NCha, nsamples in acquisition order
kspace_data = permute(kspace_data,[3 1 2]);
kspace_data = reshape(kspace_data,size(kspace_data,1),size(kspace_data,2),[],size(kspace_data,3));

VerString=twix_obj.image.softwareVersion;
switch VerString
  case 'vb'
    params = twix_obj.hdr.Meas;
    params.alTR_seconds = params.alTR(1)*1e-6;
    params.alTE_seconds=params.alTE(1:2)*1e-6;
    params.dThickness_mm = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1,1}.dThickness;
    params.dPhaseFOV_mm = params.PhaseFoV;
    params.dReadoutFOV_mm = params.ReadFoV;
    params.lEchoSpacing=params.lEchoSpacing*1e-6;
    isVE=false;
  case 'vd' %Actually VE
    params = twix_obj.hdr.Config;
    params.alTR_seconds = params.TR*1e-6;
    params.alTE_seconds=twix_obj.hdr.MeasYaps.alTE{1}*1e-6;
    params.dThickness_mm = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1,1}.dThickness;
    params.adFlipAngleDegree=twix_obj.hdr.MeasYaps.adFlipAngleDegree{1};
    params.dPhaseFOV_mm = params.PhaseFOV;
    params.dReadoutFOV_mm = params.ReadoutFOV;
    params.lSegments=twix_obj.hdr.MeasYaps.sFastImaging.lSegments;
    params.lBaseResolution=params.BaseResolution;
    params.lEchoSpacing=params.alTR_seconds/params.lSegments
    isVE=true;
    fprintf('FIND ECHO SPACING IN HEADER!\n')
  otherwise
    fprintf('Version error!')
    keyboard;
end
fprintf('done.\n');

%%
switch ScanType
  case 'Cine'
    cutoff = 120*~isVE; %throw out approach to steady-state
    [~,coil]=max(sqrt(sum(sum(sum(abs(kspace_data).^2,1),2),3)));
    figure,plot(abs(kspace_data((1+~isVE):2:end,params.lBaseResolution+1,coil)));
    disp('Set cutoff with ''cutoff=first_segment*2''')
    keyboard
    cutoff_end = min(cutoff + ceil(total_time*params.lSegments/(2*params.alTR_seconds))*2,2*floor(size(kspace_data,1)/2));
  case {'IR','T2prep'}
    cutoff=params.lSegments*~isVE; %throw out first recovery period
    %RY 
    if start_time~=0
    cutoff=floor(start_time/(params.alTR_seconds))*(params.lSegments);
    end
    %RY end
    cutoff_end = min(cutoff + ceil(total_time/(params.alTR_seconds))*(params.lSegments),params.lSegments*floor(size(kspace_data,1)/params.lSegments));
  case 'T2IR'
    cutoff=5*params.lSegments*~isVE; %throw out first 5 recovery periods (unsolved for VE dummy hearbeats)
    cutoff_end = min(cutoff + ceil(total_time/params.alTR_seconds/5)*5*params.lSegments,5*params.lSegments*floor(size(kspace_data,1)/params.lSegments/5));
    if exist('rep','var')
      for j=1:rep-1
        cutoff=cutoff_end;
        cutoff_end = min(cutoff + ceil(total_time/params.alTR_seconds/5)*5*params.lSegments,5*params.lSegments*floor(size(kspace_data,1)/params.lSegments/5));
      end
      cutoff
      cutoff_end
    end
  case 'SR'
    cutoff=params.lSegments*~isVE; %throw out first recovery period
    cutoff_end = min(ceil(total_time/params.alTR_seconds)*params.lSegments,params.lSegments*floor(size(kspace_data,1)/params.lSegments));
    figure,imagesc(abs(kspace_data((params.lSegments-isVE):params.lSegments:end,:)));
    disp('Set cutoff_end with ''cutoff_end=last_segment*params.lSegments''')
    keyboard
end

if start_time~=0
    cutoff=floor(start_time/(params.alTR_seconds))*(params.lSegments);
end

if isVE
  nav_data = double(kspace_data((cutoff+1):2:cutoff_end,:,:,:));
  kspace_data = double(kspace_data((cutoff+2):2:cutoff_end,:,:,:));
else
  nav_data = double(kspace_data((cutoff+2):2:cutoff_end,:,:,:));
  kspace_data = double(kspace_data((cutoff+1):2:cutoff_end,:,:,:));
end
total_time = size(kspace_data,1)/params.lSegments*2*params.alTR_seconds
cutoffline=cutoff/2;
Ncoils=size(kspace_data,4);
Norig = params.lBaseResolution;
ovs = Norig; %FOV extension
N=Norig + ovs; %oversample
DC=Norig + 1;

Ny = N;
Nx = N;