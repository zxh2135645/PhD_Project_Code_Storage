function [reconME, reconME_iso, QSMparams] = genMEimagesfromFitParams(fitParams,cphase,rphase)

% load fitting parameters
extractVarFromStruct(fitParams);

if nargin < 3 || ~rphase>=1 || ~ rphase<=size(Phi,4)
    rphase = 1;
end
 
if nargin < 2 || ~cphase>=1 || ~ cphase<=size(Phi,3)
    cphase = 1;
end

% recon parameters
L     = size(Gr,1);
Necho = size(Phi,6);
tempNecho = size(Phi,1)/L;

Phi = reshape(permute(Phi,[1 2 5 3 4 6]),L*tempNecho,Nseg*moduleLength,size(Phi,3),size(Phi,4),Necho);

[Ny,Nx,Nz,~,~] = size(U);
dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), floor(Nz/2)-floor(Nzdisp/2) + (1:Nzdisp), :, :);
Utemp = reshape(dispim(U),Nydisp,Nxdisp,Nzdisp,tempNecho,L);

%% Choose image time points

% tIdx = Nseg:Nseg:Nseg*moduleLength;
tIdx = Nseg;

reconME = zeros(Nydisp,Nxdisp,Nzdisp,Necho);
for echo = 1:Necho
    if tempNecho == Necho
        PhiT1 = Gr\reshape(mean(Phi(echo:tempNecho:end,tIdx,cphase,rphase,echo),2),L,[]);
        temp = reshape(Utemp(:,:,:,echo,:),[],L);
    else
        PhiT1 = Gr\reshape(mean(Phi(:,params.linesPerShot:params.linesPerShot:end,cphase,rphase,echo),2),L,[]);
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
[~,Idx] = sort([tempIdx1 tempIdx2 minSpacing]);
newRatio = [newRatio1 newRatio2 1];
newRatio = newRatio(Idx);
newRatio(3) = 1;
size2 = @(x) [size(x,1) size(x,2)];

reconME_re  = imresize(real(reconME),newRatio(1:2).*size2(reconME));
reconME_im  = imresize(imag(reconME),newRatio(1:2).*size2(reconME));
reconME_iso = reconME_re + 1j*reconME_im;

Affine3D = cat(2,params.vecRow(:),params.vecCol(:),params.vecNorm(:));
B0_dir = Affine3D\[0 0 1]';

size3 = @(x) [size(x,1) size(x,2) size(x,3)];
QSMparams.voxel_size = rawVoxelSpacing(:).';
QSMparams.voxel_size_iso = rawVoxelSpacing(:).'./newRatio(:).';
QSMparams.matrix_size = size3(reconME);
QSMparams.matrix_size_iso = size3(reconME_iso);
QSMparams.TE = TEarray;
QSMparams.delta_TE = TEarray(2) - TEarray(1);
QSMparams.CF = params.lResonanceFrequency;
QSMparams.B0_dir = B0_dir;






