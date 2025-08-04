function fitResult = runQSM(fitParams)


% This is an implementation of the Morphology Enabled Dipole Inversion (MEDI)
% method for reconstructing a Quantitative Susceptibility Map from MR data.
% The code is not fully optimized and is given for educational purpose.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To use this tool box, add MEDI_toolbox to your MATLAB Path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% USAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% EXAMPLE DATASETS %%%%%%%%%%%%%%%%%%%%%%
% Example datasets can be found MEDI_data 
% 01_Numerical_phantom contains the simulation in Neuroimage 2012;59(3):2560-8.
% 02_Wienieff_Liu contains a numerical brain
% 03_Invivo_GE contains a human brain dataset acquired from a GE scanner
% 04_Invivo_Siemens contains a human brain dataset acquired from a Siemens scanner
%%%%%%%%%%%%%%%%%%%%%%%%%% EXAMPLE DATASETS %%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% NAMING CONVENTION %%%%%%%%%%%%%%%%%%%%%%%%%
%
% high dimensional variable
%
% iField - 4, or 5 dimensional complex MRI dataset. 
%          the 4th dimension is echo 
%          the 5th dimension is channel
%
% 3D variables
%
% Mask - binary mask denoting the region of interest
% iMag - magnitude image, square root of squares of all echoes
% iFreq_raw - the raw field map, which may contain wrapping, 
%             unit in rad/echo
% N_std - estimated noise standard deviation on iFreq_raw
% iFreq - the unwrapped field map, aka total field
%         unit in rad/echo
% RDF - Relative Difference Field, aka local field
%       unit in rad/echo
% R2s - R2* map
% QSM - Quantitative Susceptibility Map, 
%       unit in parts per million, aka ppm
%
% vectors
%
% B0_dir - unit vector representing direction of B0 field 
% matrix_size - sizes ([x y z]) of the imaging volume
% voxel_size - size of the voxel
%              unit in mm
% TE - echo time, unit in sec
%
% scalars
%
% delta_TE - echo spacing, unit in sec
% CF - center frequency, unit in Hz
% B0_strength - magnetic field strength, unit in Tesla
%%%%%%%%%%%%%%%%%%%%%% NAMING CONVENTION %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% NECESSARY FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare multi-echo images
%% initialize fitting parameters

cphases = 1;
rphase = 1;

BETfrac = 0.35;
BETgrad = 0;

% SWI
power_weighting = 4;
threshold = 1;
filterSize = 65;

% QSM thresholds for tSWI calculation
chi_1 = 0;
chi_2 = 550;

% load fitting parameters, overwrite initial values
extractVarFromStruct(fitParams);

% recon parameters
L     = size(Gr,1);
Necho = size(Phi,6);
tempNecho = size(Phi,1)/L;

Phi = reshape(permute(Phi,[1 2 5 3 4 6]),L*tempNecho,Nseg*moduleLength,size(Phi,3),size(Phi,4),Necho);

[Ny,Nx,Nz,~,~] = size(U);
dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), floor(Nz/2)-floor(Nzdisp/2) + (1:Nzdisp), :, :);
Utemp = reshape(dispim(U),Nydisp,Nxdisp,Nzdisp,tempNecho,L);

% inline functions
vec = @(x) x(:);
row    = @(x) x(:).';
size2 = @(x) [size(x,1) size(x,2)];
size3 = @(x) [size(x,1) size(x,2) size(x,3)];

%% Choose image time points

rphase = max(min(rphase,size(Phi,4)),1);
cphases(cphases>size(Phi,3)) = [];
if isempty(cphases)
    cphases = 1;
end

% tIdx = Nseg:Nseg:Nseg*moduleLength;
tIdx = Nseg;

reconME = zeros(Nydisp,Nxdisp,Nzdisp,Necho);
for echo = 1:Necho
    if tempNecho == Necho
        PhiT1 = Gr\reshape(mean(Phi(echo:tempNecho:end,tIdx,cphases(1),rphase,echo),2),L,[]);
        temp = reshape(Utemp(:,:,:,echo,:),[],L);
    else
        PhiT1 = Gr\reshape(mean(Phi(:,params.linesPerShot:params.linesPerShot:end,cphases(1),rphase,echo),2),L,[]);
        temp = reshape(Utemp(:,:,:,1,:),[],L);
    end
    reconME(:,:,:,echo) = reshape(temp*PhiT1,Nydisp,Nxdisp,[],1);
end
reconME = reconME./max(abs(reconME(:)));

rawVoxelSpacing = params.rawVoxelSpacing;
[~,minSpacing] = min(rawVoxelSpacing);
temp = 1:3; temp(minSpacing) = [];
tempIdx1 = temp(1);
tempIdx2 = temp(2);
newRatio1 = rawVoxelSpacing(temp(1))/rawVoxelSpacing(minSpacing);
newRatio2 = rawVoxelSpacing(temp(2))/rawVoxelSpacing(minSpacing);
[~,Idx] = sort([[tempIdx1 tempIdx2 minSpacing]]);
newRatio = [newRatio1 newRatio2 1];
newRatio = newRatio(Idx);
newRatio(3) = 1;
size2 = @(x) [size(x,1) size(x,2)];
size3 = @(x) [size(x,1) size(x,2) size(x,3)];

reconME_re = imresize(real(reconME),newRatio(1:2).*size2(reconME));
reconME_im = imresize(imag(reconME),newRatio(1:2).*size2(reconME));
reconME = reconME_re + 1j*reconME_im;

QSMparams.iMag = sqrt(sum(abs(reconME).^2,4));

Affine3D = cat(2,vec(params.vecRow),vec(params.vecCol),vec(params.vecNorm));
B0_dir = Affine3D\[0 0 1]';

QSMparams.voxel_size = row(rawVoxelSpacing(:)./newRatio(:));
QSMparams.matrix_size = size3(reconME);
QSMparams.TEarray = TEarray;
QSMparams.delta_TE = TEarray(2) - TEarray(1);
QSMparams.CF = params.lResonanceFrequency;
QSMparams.B0_dir = B0_dir;

% Estimate the frequency offset in each of the voxel using a complex
% fitting (uneven echo spacing)
% Spatial phase unwrapping (region-growing)
fprintf('Phase unwrapping...\n');
QSMparams = QSMestimateTotalFieldfromME(reconME,QSMparams);

% Use FSL BET to extract brain mask
fprintf('Brain extraction (BET)...\n');
QSMparams.Mask = BET(QSMparams.iMag,QSMparams.matrix_size,QSMparams.voxel_size,BETfrac,BETgrad);

% Background field removal 
QSMparams = QSMbackgroundFieldRemoval(QSMparams);

totalField_ppm = QSMparams.totalField/(QSMparams.delta_TE*2*pi*QSMparams.CF)*1e6;
RDF_ppm = QSMparams.RDF/(QSMparams.delta_TE*2*pi*QSMparams.CF)*1e6;

figure;imagesc(imageOrientLPS(totalField_ppm(:,:,floor(Nz/2)+1),params,0),[-1 1]);axis equal tight;title('Total Field  (ppm)')
figure;imagesc(imageOrientLPS(RDF_ppm(:,:,floor(Nz/2)+1),params,0),[-0.15 0.15]);axis equal tight;title('Tissue Field  (ppm)')
drawnow;

% % R2* map needed for ventricular CSF mask
% R2s = arlo(TEarray, abs(reconME));
% fitResult.R2smap = R2s;

% Ventricular CSF mask for zero referencing 
temp = abs(reconME(:,:,:,1));
[idx,centroids] = kmeans(temp(:).^2,4);
idx = reshape(idx,size(temp));
[~,b] = sort(centroids);
temp1 = idx;temp1(idx~=b(2)) = 0;
QSMparams.Mask_CSF = temp1.*QSMparams.Mask;

QSM = QSMestimateChi(QSMparams);

figure;imagesc(imageOrientLPS(QSM(:,:,floor(Nz/2)+1),params,0),[-0.15 0.15]);axis equal tight;title('QSM')

fitResult.BETmask  = QSMparams.Mask;
fitResult.Mask_CSF = QSMparams.Mask_CSF;
fitResult.totalField_ppm = totalField_ppm;
fitResult.RDF_ppm = RDF_ppm;
fitResult.QSM  = QSM;

if flagDoSWI
    [pSWI,nSWI,tSWI,phase_hp] = SWI(QSMparams.iMag,angle(reconME(:,:,:,end)),threshold,power_weighting,QSM,chi_1,chi_2);
    fitResult.pSWI = pSWI;
    fitResult.nSWI = nSWI;
    fitResult.tSWI = tSWI;

    figure;imagesc(imageOrientLPS(tSWI(:,:,floor(Nz/2)+1),params,0),[-0.15 0.15]);axis equal tight;title('tSWI')
end

disp('Done.');




