function reconstructedImages = displayTensor_DCE(params,reconOptions,dataArray,temporalBasis,spatialCoeff,reconstructedImages,scale,tIdx,cIdx,rIdx,useL,slice,DCEIdx,dim,eIdx)

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(temporalBasis);
extractVarFromStruct(spatialCoeff);

Necho = size(Phi,6);
L = size(Phi,1)/Necho;

if nargin < 15 || isempty(eIdx) || eIdx > Necho || eIdx < 1
    eIdx = 1;
end

if nargin < 14 || dim > 3
    dim = 3;
end
permuteOrder = [1 2 3];
permuteOrder(dim) = [];
permuteOrder = [permuteOrder dim 4];
tempN  = [Ny Nx Nz];
Nzdisp = tempN(dim);

if nargin < 13 || DCEIdx > size(Phi,5) - 10 || DCEIdx < 11
    DCEIdx = size(Phi,5) - 10;
end
if size(Phi,5) < 21
    DCEIdx = size(Phi,5);
    DCEave = 1:size(Phi,5);
else
    DCEave = DCEIdx-10:DCEIdx+10;
end

if nargin < 12 || isempty(slice)
    if isfield(reconOptions,'dispSlice')
        slice = reconOptions.dispSlice;
    elseif Nzdisp > 3
        slice = floor(Nzdisp/2) + 1;
    else
        slice = 1:Nzdisp;
    end
end
if ~(1<=slice(1) && slice(end)<=Nzdisp)
    if Nz > 3
        slice = floor(Nzdisp/2) + 1;
    else
        slice = 1:Nzdisp;
    end
end

if nargin < 11 || isempty(useL)
    useL = 1:L;
end
for l = numel(useL):-1:1
    if useL(l) > size(Phi,1)
        useL(l) = [];
    end
end

if nargin < 10 || isempty(rIdx)
    rIdx = 1;
end
if ~(1<=rIdx && rIdx<=size(Phi,4))
    rIdx = 1;
end

if nargin < 9 || isempty(cIdx)
    cIdx = dataArray.diastoleIdx;
end
if ~(1<=cIdx && cIdx<=size(Phi,3))
    cIdx = 1;
end

if nargin < 8 || isempty(tIdx)
    tIdx = floor(1.5/lEchoSpacing);
end
if ~(1<=tIdx && tIdx<=size(Phi,2))
    tIdx = min(floor(1.5/lEchoSpacing),size(Phi,2));
end
if strcmp(ScanType,'Cine')
    tIdx = 1;
end

if nargin < 7
    scale = 1;
end

dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), :, :);

Utemp = reshape(U,Ny,Nx,Nz,[]);
if MBfactor == 1
    Utemp = fftshift(Utemp,3);
end
if params.isCartesian
    Utemp = fftshift(Utemp,1);
end
Utemp = permute(dispim(Utemp),permuteOrder);
Utemp = reshape(Utemp,Nydisp,Nxdisp,Nz,Necho,L);

tempNdisp = [Nydisp Nxdisp Nz L];
Nydisp = tempNdisp(permuteOrder(1)); Nxdisp = tempNdisp(permuteOrder(2));
Nzdisp = numel(slice);

Phi1 = Gr\mean(reshape(Phi(1:Necho:end,tIdx,cIdx,rIdx,DCEave,eIdx),L,[]),2);
reconVolume = reshape(Utemp(:,:,:,eIdx,:),[],L);
reconVolume = reshape(reconVolume(:,useL)*Phi1(useL,:),Nydisp,Nxdisp,[]);

reconEchos = zeros(Nydisp,Nxdisp,Nzdisp,Necho);
for n = 1:Necho
    Phi1 = Gr\reshape(Phi(n:Necho:end,tIdx,cIdx,rIdx,DCEIdx,n),L,[]);
    temp = reshape(Utemp(:,:,slice,n,:),[],L);
    reconEchos(:,:,:,n) = reshape(temp(:,useL)*Phi1(useL,:),Nydisp,Nxdisp,Nzdisp,[]);
end

Utemp = reshape(Utemp(:,:,slice,eIdx,:),[],L);
Phi   = Phi(eIdx:Necho:end,:,:,:,:,eIdx);

[~,T_T1,~,T_resp,T_DCE] = size(Phi);
PhiT1   = Gr\reshape(Phi(:,1:ceil(T_T1/10):end,cIdx,rIdx,DCEIdx),L,[]);
reconRecovery = reshape(Utemp(:,useL)*PhiT1(useL,:),Nydisp,Nxdisp,Nzdisp,[]);
reconRecovery = real(reconRecovery.*exp(-1j*angle(reconRecovery(:,:,:,end))));
reconRecovery = reconRecovery/prctile(abs(reconRecovery(:)),99.9);
% reconRecovery(reconRecovery<0) = 0;

PhiCard = Gr\reshape(Phi(:,tIdx,:,rIdx,DCEIdx),L,[]);
reconCard = reshape(Utemp(:,useL)*PhiCard(useL,:),Nydisp,Nxdisp,Nzdisp,[]);
reconCard = reconCard/prctile(abs(reconCard(:)),99.9);

PhiResp = Gr\reshape(Phi(:,tIdx,cIdx,:,DCEIdx),L,[]);
reconResp = reshape(Utemp(:,useL)*PhiResp(useL,:),Nydisp,Nxdisp,Nzdisp,[]);
reconResp = reconResp/prctile(abs(reconResp(:)),99.9);

PhiDCE  = Gr\reshape(Phi(:,tIdx,cIdx,rIdx,:),L,[]);
reconDCE = reshape(Utemp(:,useL)*PhiDCE(useL,:),Nydisp,Nxdisp,Nzdisp,[]);
reconDCE = reconDCE/prctile(abs(reconDCE(:)),99.9);

% % POCS
if strcmp(ScanType,'Cine')
    dim = 2;
    iterations = 4;
    reconRecovery = pocs(reconRecovery,dim,iterations);
    reconCard     = pocs(reconCard,dim,iterations);
    reconResp     = pocs(reconResp,dim,iterations);
    reconVolume = pocs(reconVolume,dim,iterations);
    reconEchos = pocs(reconEchos,dim,iterations);
end

                    
% Display
if flagCommandLine
    implayZoom(abs(reconRecovery(:,:,1,:))*scale,10,[],'Recovery Curve');
    if size(Phi,5) > 1
        implayZoom(abs(reconDCE(:,:,1,:))*scale,size(reconDCE,4)/5,[],'DCE');
    end
    if cbins > 1
        implayZoom(abs(reconCard(:,:,1,:))*scale,2,[],'Respiratory Bins');
    end
    if rbins > 1
        implayZoom(abs(reconResp(:,:,1,:))*scale,2,[],'Respiratory Bins');
    end
    if Nz > 1
        implayZoom(abs(reconVolume)*scale,ceil(Nz/8),[],'3D Volume');
    end
    if Necho > 1
        implayZoom(abs(reconEchos)*scale,2,[],'Multi-Echo');
    end
end

%% save gifs

fps = 20;
Nskip = floor(1/(params.lEchoSpacing*fps));

[~,filePrefix] = fileparts(params.filePath);    
try
    saveGif(scale*abs(reshape(reconRecovery(:,:,1,1:Nskip:end),Nydisp,Nxdisp,1,[])),filePath,[filePrefix '_reconRecovery.gif'],fps);
    saveGif(scale*abs(reconCard(:,:,1,:)),filePath,[filePrefix '_reconCard.gif'],min(16,size(reconCard,4)));
    saveGif(scale*abs(reconResp(:,:,1,:)),filePath,[filePrefix '_reconResp.gif'],2);
    if Nz > 1
        saveGif(scale*abs(reshape(reconVolume,Nydisp,Nxdisp,1,[])),filePath,[filePrefix '_reconVolume.gif'],max(2,Nz/8));
    end
    if Necho > 1
        saveGif(scale*abs(reconEchos(:,:,1,:)),filePath,[filePrefix '_reconEcho.gif'],2);
    end
    if size(Phi,5) > 1
        saveGif(scale*abs(reconDCE(:,:,1,:)),filePath,[filePrefix '_reconDCE.gif'],2);
    end
catch
    fprintf("Saving images to %s failed. Saving them in current folder.", filePath);
    saveGif(scale*abs(reshape(reconRecovery(:,:,1,1:Nskip:end),Nydisp,Nxdisp,1,[])),'.',[filePrefix '_reconRecovery.gif'],fps);
    saveGif(scale*abs(reconCard(:,:,1,:)),filePath,[filePrefix '_reconCard.gif'],min(16,size(reconCard,4)));
    saveGif(scale*abs(reconResp(:,:,1,:)),'.',[filePrefix '_reconResp.gif'],2);
    if Nz > 1
        saveGif(scale*abs(reshape(reconVolume,Nydisp,Nxdisp,1,[])),'.',[filePrefix '_reconVolume.gif'],max(2,Nz/8));
    end
    if Necho > 1
        saveGif(scale*abs(reconEchos(:,:,1,:)),'.',[filePrefix '_reconEcho.gif'],2);
    end
    if size(Phi,5) > 1
        saveGif(scale*abs(reconDCE(:,:,1,:)),'.',[filePrefix '_reconDCE.gif'],2);
    end
end

% Store recon images
reconstructedImages.reconRecovery = reconRecovery;
reconstructedImages.reconCard = reconCard;
reconstructedImages.reconResp = reconResp;
if Nz > 1
    reconstructedImages.reconVolume = reconVolume;
end
if Necho > 1
    reconstructedImages.reconEchos = reconEchos;
end
if size(Phi,5) > 1
    reconstructedImages.reconDCE = reconDCE;
end

% POCS
function img = pocs(img,dim,iterations)
N = size(img,dim);
tempk = fft(img,[],dim);
asymmask = true(1,N);
asymmask(N/2+1:floor(N*15/16)) = 0; 
centerwin = double(asymmask & flip(asymmask,dim));
centerwin(logical(centerwin)) = ifftshift(hamming(sum(centerwin),'periodic'));
imphase = sign(ifft(bsxfun(@times,tempk,centerwin),[],dim));
for j = 1:iterations %POCS iterations
    temp = abs(img).*imphase;
    temp = fft(temp,[],dim);
    temp(:,asymmask,:) = tempk(:,asymmask,:);
    img  = ifft(temp,[],dim);
end

