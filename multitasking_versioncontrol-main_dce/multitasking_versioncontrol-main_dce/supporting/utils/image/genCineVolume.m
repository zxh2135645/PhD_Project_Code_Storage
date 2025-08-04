function reconCard = genCineVolume(params,reconOptions,dataArray,temporalBasis,spatialCoeff,tIdx,rIdx,useL)

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(temporalBasis);
extractVarFromStruct(spatialCoeff);

L = size(Phi,1);

if nargin < 8 || isempty(useL)
    useL = 1:L;
end
for l = numel(useL):-1:1
    if useL(l) > L
        useL(l) = [];
    end
end

if nargin < 7 || isempty(rIdx)
    rIdx = 1;
end
if ~(1<=rIdx && rIdx<=size(Phi,4))
    rIdx = 1;
end

if nargin < 6 || isempty(tIdx)
    tIdx = floor(1.5/lEchoSpacing);
end
if ~(1<=tIdx && tIdx<=size(Phi,2))
    tIdx = floor(size(Phi,2)/2) + 1;
end
if strcmp(ScanType,'Cine')
    tIdx = 1;
end

dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), floor(Nz/2)-floor(Nzorig/2) + (1:Nzorig), :);

Utemp = reshape(U,Ny,Nx,Nz,[]);
if MBfactor == 1
    Utemp = fftshift(Utemp,3);
end
if params.isCartesian
    Utemp = fftshift(Utemp,1);
end

Utemp = reshape(dispim(Utemp),[],L);

PhiCard   = Gr\reshape(Phi(:,tIdx,:,rIdx,end),L,[]);
reconCard = reshape(Utemp(:,useL)*PhiCard(useL,:),Nydisp,Nxdisp,Nzorig,[]);
  
% % POCS
if strcmp(ScanType,'Cine')
    dim = 2;
    iterations = 4;
    reconCard = pocs(reconCard,dim,iterations);
end

% Intensity normalization
reconCard     = reconCard/prctile(abs(reconCard(:)),99);

% Display
if flagCommandLine
    implayZoom(abs(reconCard(:,:,floor(Nzorig/2)+1,:)),size(reconCard,4));
end


%% POCS
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

