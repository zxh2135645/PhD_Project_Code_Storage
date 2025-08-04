function reconstructedImages = displayTensor(params,reconOptions,dataArray,temporalBasis,spatialCoeff,reconstructedImages,scale,tIdx,cIdx,rIdx,useL,slice,DCEIdx,eIdx,dim)

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(temporalBasis);
extractVarFromStruct(spatialCoeff);

L = size(Gr,1);
Necho = size(Phi,6);
tempNecho = size(Phi,1)/L;

if nargin < 15 || dim > 3
    dim = 3;
end
permuteOrder = [1 2 3];
permuteOrder(dim) = [];
permuteOrder = [permuteOrder dim 4];
tempN  = [Ny Nxdisp Nzdisp];
if dim ~= 3
    Nzmax = tempN(dim);
else
    Nzmax = Nzdisp;
end
rawVoxelSpacing = rawVoxelSpacing(permuteOrder(1:3));

if nargin < 14 || isempty(eIdx) || eIdx > Necho || eIdx < 1
    eIdx = 1;
end

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
    if Nzmax <= 3
        slice = 1:Nzmax;
    elseif isfield(reconOptions,'dispSlice')
        slice = reconOptions.dispSlice;
    else %if Nzmax > 3
        slice = floor(Nzmax/2) + 1;
    end
end
if ~(1<=slice(1) && slice(end)<=Nzmax)
    if Nz > 3
        slice = floor(Nzmax/2) + 1;
    else
        slice = 1:Nzmax;
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
elseif strcmp(ScanType,'CEST')
    tIdx = size(Phi,2) + 1 - linesPerShot;
end

if nargin < 7
    scale = 1;
end

dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), floor(Nz/2)-floor(Nzdisp/2) + (1:Nzdisp), :);

Utemp = reshape(U,Ny,Nx,Nz,[]);
if MBfactor == 1
    Utemp = fftshift(Utemp,3);
end
if params.isCartesian
    Utemp = fftshift(Utemp,1);
end
Utemp = permute(dispim(Utemp),permuteOrder);
[Nydisp,Nxdisp,Nzdisp,~] = size(Utemp);
Utemp = reshape(Utemp,Nydisp,Nxdisp,Nzdisp,tempNecho,L);

% Multi-echo image
reconEchos = zeros(Nydisp,Nxdisp,numel(slice),Necho);
if Necho > 1
    for n = 1:Necho
        if tempNecho == Necho
            Phi1 = Gr\reshape(Phi(n:tempNecho:end,tIdx,cIdx,rIdx,DCEIdx,n),L,[]);
            temp = reshape(Utemp(:,:,slice,n,:),[],L);
        else
            Phi1 = Gr\reshape(Phi(:,tIdx,cIdx,rIdx,DCEIdx,n),L,[]);
            temp = reshape(Utemp(:,:,slice,1,:),[],L);
        end
        reconEchos(:,:,:,n) = reshape(temp(:,useL)*Phi1(useL,:),Nydisp,Nxdisp,numel(slice));
    end
    reconEchos = reconEchos/prctile(abs(reconEchos(:)),99.9);
end

% keep only the selected echo
if tempNecho == Necho
    Phi = Phi(eIdx:Necho:end,:,:,:,:,eIdx);
else
    Phi = Phi(:,:,:,:,:,eIdx);
end

% 3D volume image
if Nzdisp > 1
    Phi1 = Gr\mean(reshape(Phi(:,tIdx,cIdx,rIdx,DCEave),L,[]),2);
    reconVolume = reshape(reshape(Utemp(:,:,:,eIdx,:),[],L)*Phi1,Nydisp,Nxdisp,Nzdisp,[]);
    reconVolume = real(reconVolume.*exp(-1j*angle(reconVolume(:,:,:,end))));
    %reconVolume = reshape(reconVolume(:,useL)*Phi1(useL,:),Nydisp,Nxdisp,[]);
    reconVolume = reconVolume/prctile(abs(reconVolume(:)),99.9);
end

% keep only the selected slice(s)
Utemp = reshape(Utemp(:,:,slice,eIdx,:),[],L);

% dynamic images
fps = 20;
Nskip = floor(1/(params.lEchoSpacing*fps));
Nskip = Nskip + mod(Nskip,2);

% Contrast change
PhiT1   = Gr\reshape(Phi(:,1:1:end,cIdx,rIdx,DCEave),L,[]);%Gr\reshape(Phi(:,:,cIdx,rIdx,end-7),L,[]);%
reconRecovery = reshape(Utemp(:,useL)*PhiT1(useL,:),Nydisp,Nxdisp,numel(slice),[]);
reconRecovery = real(reconRecovery.*exp(-1j*angle(reconRecovery(:,:,:,end))));
reconRecovery = reconRecovery/prctile(abs(reconRecovery(:)),99.9);
% reconRecovery(reconRecovery<0) = 0;

% cardiac phases
PhiCard = Gr\reshape(Phi(:,tIdx,:,rIdx,DCEIdx),L,[]);
reconCard = reshape(Utemp(:,useL)*PhiCard(useL,:),Nydisp,Nxdisp,numel(slice),[]);
reconCard = reconCard/prctile(abs(reconCard(:)),99.9);

% respiratory phases
PhiResp = Gr\reshape(Phi(:,tIdx,cIdx,:,DCEIdx),L,[]);
reconResp = reshape(Utemp(:,useL)*PhiResp(useL,:),Nydisp,Nxdisp,numel(slice),[]);
reconResp = reconResp/prctile(abs(reconResp(:)),99.9);

% Show coronal view for abdominal scans
% if Nz > 1
%     temp1 = reshape(temp(:,useL)*PhiResp(useL,:),1,Nxdisp,Nz,[]);
%     temp1 = flipud(permute(temp1,[3 2 1 4]));
%     temp1 = abs(temp1)/prctile(abs(temp1(:)),99.9);
%     temp1 = imresize(temp1,[floor(Nz*voxelSpacing(3)/voxelSpacing(2)) Nxdisp]);
%     implayZoom(0.8*temp1,2,[],'Respiratory Bins (coronal)');
%     [~,filePrefix] = fileparts(params.filePath);
%     saveGif(0.8*abs(temp1(:,:,1,:)),filePath,[filePrefix '_reconRespCoronal.gif'],2);
% end

% DCE
if strcmp(ScanType,'CEST')
    PhiDCE   = Gr\reshape(mean(Phi(:,size(Phi,2) + 1 - linesPerShot*(1:16),cIdx,rIdx,:),2),L,[]);
else
    PhiDCE   = Gr\reshape(Phi(:,tIdx,cIdx,rIdx,:),L,[]);
end
reconDCE = reshape(Utemp(:,useL)*PhiDCE(useL,:),Nydisp,Nxdisp,numel(slice),[]);
reconDCE = reconDCE/prctile(abs(reconDCE(:)),99.9);

% POCS
try
    tempflag = isAsymmEcho;
catch
    tempflag = strcmp(ScanType,'Cine');
end
if tempflag
    dim = 2;
    iterations = 4;
    reconRecovery = pocs(reconRecovery,dim,iterations);
    reconCard     = pocs(reconCard,dim,iterations);
    reconResp     = pocs(reconResp,dim,iterations);
    reconVolume   = pocs(reconVolume,dim,iterations);
    reconEchos    = pocs(reconEchos,dim,iterations);
end

% interpolate 2D image to isotropic pixels
[~,minSpacing] = min(rawVoxelSpacing);
temp = 1:3; temp(minSpacing) = [];
tempIdx1 = temp(1);
tempIdx2 = temp(2);
newRatio1 = rawVoxelSpacing(temp(1))/rawVoxelSpacing(minSpacing);
newRatio2 = rawVoxelSpacing(temp(2))/rawVoxelSpacing(minSpacing);
[~,Idx] = sort([tempIdx1 tempIdx2 minSpacing]);
newRatio = [newRatio1 newRatio2 1];
newRatio = newRatio(Idx(1:2));
size2 = @(x) [size(x,1) size(x,2)];
size3 = @(x) [size(x,1) size(x,2) size(x,3)];

dispSlice = floor(numel(slice)/2) + 1;

% display images
if flagCommandLine
    if size(Phi,2) > 1 && ~strcmp(ScanType,'CEST')
        implayZoom(imageOrientLPS(abs(reconRecovery(:,:,dispSlice,:))*scale,params),fps,[],'Recovery Curve');
    end
    if size(Phi,5) > 12
        if strcmp(ScanType,'CEST')
            implayZoom(imageOrientLPS(abs(reconDCE(:,:,dispSlice,:))*scale,params),size(reconDCE,4)/5,[],'Sat Freq Offset');            
        else
            implayZoom(imageOrientLPS(abs(reconDCE(:,:,dispSlice,:))*scale,params),size(reconDCE,4)/5,[],'DCE');
        end
    end
    if cbins > 1
        implayZoom(imageOrientLPS(abs(reconCard(:,:,dispSlice,:))*scale,params),20,[],'Cardiac Bins');
    end
    if rbins > 1
        implayZoom(imageOrientLPS(abs(reconResp(:,:,dispSlice,:))*scale,params),2,[],'Respiratory Bins');
    end
    if Nzdisp > 1
        implayZoom(imageOrientLPS(abs(reconVolume)*scale,params),ceil(Nz/8),[],'3D Volume');
    end
    if Necho > 1
        implayZoom(imageOrientLPS(abs(reconEchos(:,:,dispSlice,:))*scale,params),2,[],'Multi-Echo');
    end
end

%% save gifs

fps = 20;
Nskip = floor(1/(params.lEchoSpacing*fps));

try
    [~,filePrefix] = fileparts(params.filePath);
    if size(Phi,2) > 1
        saveGif(imageOrientLPS(scale*abs(reshape(reconRecovery(:,:,1,:),Nydisp,Nxdisp,1,[])),params),filePath,[filePrefix '_reconRecovery.gif'],fps);
    end
    if size(Phi,3) > 1
        saveGif(imageOrientLPS(scale*abs(reconCard(:,:,1,:)),params),filePath,[filePrefix '_reconCard.gif'],size(reconCard,4));
    end
    if size(Phi,4) > 1
        saveGif(imageOrientLPS(scale*abs(reconResp(:,:,1,:)),params),filePath,[filePrefix '_reconResp.gif'],2);
    end
    if size(Phi,5) > 12
        saveGif(imageOrientLPS(scale*abs(reconDCE(:,:,1,:)),params),filePath,[filePrefix '_reconDCE.gif'],2);
    end
    if Nz > 1
        saveGif(imageOrientLPS(scale*abs(reshape(reconVolume,Nydisp,Nxdisp,1,[])),params),filePath,[filePrefix '_reconVolume.gif'],max(2,Nz/8));
    end
    if Necho > 1
        saveGif(imageOrientLPS(scale*abs(reshape(reconEchos,Nydisp,Nxdisp,1,[])),params),filePath,[filePrefix '_reconEchos.gif'],2);
    end
catch
    fprintf("Saving images to %s failed. Saving them in current folder.\n", filePath);
    if size(Phi,2) > 1
%        saveGif(imageOrientLPS(scale*abs(reshape(reconRecovery(:,:,1,1:Nskip:end),params),Nydisp,Nxdisp,1,[])),'.','reconRecovery.gif',fps);
    end
    if size(Phi,3) > 1
        saveGif(imageOrientLPS(scale*abs(reconCard(:,:,1,:)),params),'.','reconCard.gif',size(reconCard,4));
    end
    if size(Phi,4) > 1
        saveGif(imageOrientLPS(scale*abs(reconResp(:,:,1,:)),params),'.','reconResp.gif',2);
    end
    if size(Phi,5) > 12
        saveGif(imageOrientLPS(scale*abs(reconDCE(:,:,1,:)),params),'.',[filePrefix '_reconDCE.gif'],2);
    end
    if Nz > 1
        saveGif(imageOrientLPS(scale*abs(reshape(reconVolume,Nydisp,Nxdisp,1,[])),params),'.','reconVolume.gif',max(2,Nz/8));
    end
    if Necho > 1
        saveGif(imageOrientLPS(scale*abs(reshape(reconEchos,Nydisp,Nxdisp,1,[])),params),'.','reconEchos.gif',2);
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
    reconstructedImages.reconEchos  = reconEchos;
end
if size(Phi,5) > 12
    if strcmp(ScanType,'CEST')
        reconstructedImages.reconCEST = reconDCE;            
    else
        reconstructedImages.reconDCE = reconDCE;            
    end
end
    

% POCS
function img = pocs(img,dim,iterations)
N = size(img,dim);
tempk = ft1d(img,dim);
centerwinwidth = floor(N/16) * 2;
asymmask = zeros(N,1);
asymmask(N/2-centerwinwidth/2+1:end) = 1;
newdim = [1 2 3]; newdim(1) = dim; newdim(dim) = 1;
asymmask = permute(asymmask,newdim);
asymmask = bsxfun(@times,asymmask,ones(size(tempk)));
asymmask = logical(asymmask);
centerwin = zeros(N,1);
centerwin(floor(N/2-centerwinwidth/2)+(1:centerwinwidth)) = hamming(centerwinwidth); 
centerwin = permute(centerwin,newdim);
imphase = sign(ift1d(bsxfun(@times,tempk,centerwin),dim));
for j = 1:iterations    % POCS iterations
    temp = abs(img).*imphase;
    temp = ft1d(temp,dim);
    temp(asymmask) = tempk(asymmask);
    img  = ift1d(temp,dim);
end

