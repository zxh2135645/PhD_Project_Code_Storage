[fid_file, fid_path] = uigetfile('*.dat');

ScanType = 'T2star';
Trajectory = 'Cartesian';
isradial = strcmp('Radial',Trajectory);
iscartesian = ~isradial;

file_string = [fid_path fid_file(1:(end-4)) '_' datestr(now,30)];
mkdir(file_string);
file_string = fullfile(file_string,'multitasking');
diary([file_string '.log']);

%%
fprintf('Loading raw data...');
twix_obj = mapVBVD(strcat(fid_path,fid_file));
if numel(twix_obj)>1
  twix_obj=twix_obj{end};
end
%%
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
    params.lEchoSpacing=params.lEchoSpacing*1e-6
    isVE=false;
  case 'vd' %Actually VE
    params = twix_obj.hdr.Config;
    params.alTR_seconds = params.TR(1)*1e-6;
    params.alTE_seconds=twix_obj.hdr.MeasYaps.alTE{1}*1e-6;
    params.dThickness_mm = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1,1}.dThickness;
    params.adFlipAngleDegree=twix_obj.hdr.MeasYaps.adFlipAngleDegree{1};
    try
      params.dPhaseFOV_mm = params.PhaseFOV;
    catch
      params.dPhaseFOV_mm = params.PhaseFoV; 
    end
    try
      params.dReadoutFOV_mm = params.ReadoutFOV;
    catch
      params.dReadoutFOV_mm = params.ReadFoV;
    end
    params.lSegments=twix_obj.hdr.MeasYaps.sFastImaging.lSegments;
    params.lBaseResolution=params.BaseResolution;
    [~,lES]=unix(sprintf('grep -a lEchoSpacing -A 2 ''%s'' | grep [0-9]',fullfile(fid_path,fid_file))); %manually extract it from header
    params.lEchoSpacing=str2double(lES)*1e-6
    if isnan(params.lEchoSpacing) %VIDA handling
      [~,lES]=unix(sprintf('grep -a lEchoSpacing ''%s'' | grep [0-9] | head -1 | awk -F{ ''{print $2}'' | awk ''{print $1}''',fullfile(fid_path,fid_file)))
      params.lEchoSpacing=str2double(lES)*1e-6
    end
    params.lPartitions=params.NPar
    isVE=true;
  otherwise
    fprintf('Version error!')
    keyboard;
end
fprintf('done.\n');

%%
switch ScanType
  case 'T2star'
      % TODO: This needs to be modified strictly
    cutoff = params.lSegments*params.NEco*2; %throw out approach to steady-state
    % cutoff = 0;
    fprintf('input cutoff_end_n');
%     keyboard;
    cutoff_end = size(kspace_data,1);% - params.lSegments*params.NEco*cutoff_end_n %-params.lSegments*params.NEco*120;
%     cutoff_end = 110592;
  case 'Cine'
    cutoff = 120; %throw out approach to steady-state
    [~,coil]=max(sqrt(sum(sum(sum(abs(kspace_data).^2,1),2),3)));
    figure,plot(abs(kspace_data((1+~isVE):2:end,params.lBaseResolution+1,coil)));
    disp('Set cutoff with ''cutoff=first_segment*2''')
    keyboard
    cutoff_end = min(cutoff + ceil(total_time*params.lSegments/(2*params.alTR_seconds))*2,2*floor(size(kspace_data,1)/2));
  case {'IR','T2prep'}
    cutoff=params.lSegments; %throw out first recovery period
    cutoff_end = min(cutoff + ceil(total_time/params.alTR_seconds)*params.lSegments,params.lSegments*floor(size(kspace_data,1)/params.lSegments));
  case 'T2IR'
    cutoff=5*params.lSegments; %throw out first 5 recovery periods
    cutoff_end = min(cutoff + ceil(total_time/params.alTR_seconds/5)*5*params.lSegments,5*params.lSegments*floor(size(kspace_data,1)/params.lSegments/5));
  case 'SR'
    cutoff=params.lSegments; %throw out first recovery period
    cutoff_end = min(ceil(total_time/params.alTR_seconds)*params.lSegments,params.lSegments*floor(size(kspace_data,1)/params.lSegments));
    figure,imagesc(abs(kspace_data((params.lSegments-isVE):params.lSegments:end,:)));
    disp('Set cutoff_end with ''cutoff_end=last_segment*params.lSegments''')
    keyboard
end
%%
if iscartesian
    pe_indices = [twix_obj.image.Lin(:), twix_obj.image.Par(:)];
    pe_indices = pe_indices((cutoff+1):cutoff_end,:);
    pe_indices(:,1) = (pe_indices(:,1) - pe_indices(2,1))/mode(diff(unique(pe_indices(:,1))));
    Nfft1 = max(abs(pe_indices(:,1)))*2 + 1;
    Nfft = max(params.lBaseResolution,Nfft1);
    DC_ky = max(abs(pe_indices(:,1))) + 1;
%     DC_ky = max(abs(pe_indices(:,1)));
    pe_indices(:,1) = pe_indices(:,1) + DC_ky;
    Ny = Nfft;
    Nx = params.lBaseResolution;
end

DC_kz=floor(params.lPartitions/2+1);
%%%%%%%%%%%%%%%%%%%%
% Hard-coded
SGblock = 2;
%%%%%%%%%%%%%%%%%%%%

% Interleave imaging and nav echos
sz=size(kspace_data);
kspace_data=reshape(kspace_data,params.NEco,SGblock,[]);
kspace_data=reshape(permute(kspace_data,[2 1 3]),sz);

sz=size(pe_indices);
pe_indices=reshape(pe_indices,params.NEco,SGblock,[]);
pe_indices=reshape(permute(pe_indices,[2 1 3]),sz);
% xz
figure(); scatter(pe_indices(1:SGblock:10000,1), pe_indices(1:SGblock:10000, 2));

cutoff_shift=1; %Segments may not be a multiple of SGblock

nav_data = double(kspace_data((cutoff+cutoff_shift):SGblock:cutoff_end,:,:,:));
kspace_data = double(kspace_data((cutoff+1):cutoff_end,:,:,:));
Nread=size(kspace_data,1);
total_time = Nread/params.lSegments*params.alTR_seconds/params.NEco
nav_indices = cutoff_shift:SGblock:Nread;
if ~iscartesian
  kspace_data(nav_indices,:,:,:)=[];
end

Ncoils=size(kspace_data,4); % 12
Norig = params.lBaseResolution;
% Norig = 192;
Nzorig = params.NImagePar;
ovs = Norig*~iscartesian; %FOV extension
sl_ovs=params.lPartitions-Nzorig; %Slice oversampling extension
N=Norig + ovs; %oversample

DC=Norig + 1;

if ~iscartesian
  Ny = N;
  Nx = N;
  Nz = params.lPartitions;
else
    Ny = params.NLin;
    Nx = params.lBaseResolution;
    Nz = params.lPartitions;
end

size(kspace_data)