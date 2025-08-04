function QSM = QSMiSWIM(QSMparams)

Niter = 30;

extractVarFromStruct(QSMparams);

[Ny,Nx,Nz] = size(RDF);
[gridX,gridY,gridZ] = meshgrid(-Nx/2:Nx/2-1,-Ny/2:Ny/2-1,-Nz/2:Nz/2-1);
g_k = 1/3 - gridZ.^2./(gridX.^2+gridY.^2+gridZ.^2);

RDF = ifftshift(RDF);
g_k = ifftshift(g_k);

threshold_g = 0.1;
g_k(abs(g_k)<threshold_g) = threshold_g.*sign(g_k(abs(g_k)<threshold_g));
g_kinv = 1./g_k;
g_kinv(1,1,1) = 0;
g_kinv(~isfinite(g_kinv)) = 0;
maskGk = ones(size(g_kinv));
maskGk(abs(g_k)>threshold_g) = 0;

chi0 = ifftn(g_kinv.*fftn(-RDF.*1e6./(CF*2*pi*delta_TE)));

threshold01 = 0.07;
threshold02 = 0.15;
maskVessel = zeros(size(chi0));
maskVessel(chi0>threshold01) = 1;
maskVessel = imclose(maskVessel,strel('sphere',2));
maskVessel = wmedfilt3(maskVessel);
MIPchi0 = MIP(chi0,5);
maskMIP  = zeros(size(chi0));
maskMIP(MIPchi0>threshold02) = 1;
maskVessel = maskVessel.* maskMIP;

Chi = chi0;
for iter = 1:Niter
    tempkVessel = fftn(Chi.*maskVessel);
    tempkChi = fftn(Chi);
    tempkChi(maskGk>0) = tempkVessel(maskGk>0);
    Chi = ifftn(tempkChi);
end

QSM = fftshift(Chi);

function MIPimg = MIP(img,nslices)
if nargin < 2
    nslices = 5;
end
MIPimg = zeros(size(img));
Nz = size(img,3);
for slice = 1:Nz
    volume = img(:,:,max(1,slice-floor(nslices/2)):min(Nz,slice+floor(nslices/2)));
    MIPimg(:,:,slice) = max(volume,[],3);
end