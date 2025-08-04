function reconstructedImages = genImageTensor(params,reconOptions,dataArray,temporalBasis,spatialCoeff,reconstructedImages,dim,tIdx,cardiacIdx,respIdx,flagBB,flagFullSegment)
 
if nargin < 12
    flagFullSegment = false;
end

if nargin < 11
    flagBB = false;
end

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(temporalBasis);
extractVarFromStruct(spatialCoeff);
    
dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), :, :);

L   =  size(Phi,1);
useL = 1:L;

Utemp = reshape(U,Ny,Nx,Nz,[]);
if MBfactor == 1
    Utemp = fftshift(Utemp,3);
end
if isCartesian
    Utemp = fftshift(Utemp,1);
end
Utemp = reshape(dispim(Utemp),[],L);

if flagFullSegment
    ds = 1;
else
    fps = 20;
    ds = round(1/(fps*params.lEchoSpacing));    %downsampling to get to 20fps
end

if flagBB
    if ~isfield(dataArray,'bIdx')
        dataArray.bIdx = getBlackBloodIdx(params,reconOptions,temporalBasis,spatialCoeff);
        tIdx = abs(dataArray.bIdx);
    else
        tIdx = abs(dataArray.bIdx);
    end
end
if dim == 4
    PhiDisplay = Gr\reshape(Phi(:,tIdx,cardiacIdx,:,end),L,[]);
    recon = zeros(Nydisp,Nxdisp,Nz,size(Phi,4));
    for n = 1:size(Phi,4)
        PhiEnd = Gr\reshape(Phi(:,end,cardiacIdx,n,end),L,[]);
        temp  = reshape(Utemp(:,useL)*PhiEnd(useL,:),Nydisp,Nxdisp,Nz);
        recon(:,:,:,n) = reshape(Utemp(:,useL)*PhiDisplay(useL,n),Nydisp,Nxdisp,Nz,[]).*exp(-1j*angle(temp));
    end
elseif dim == 3
    PhiDisplay = Gr\reshape(Phi(:,tIdx,:,respIdx,end),L,[]);
    recon = zeros(Nydisp,Nxdisp,Nz,size(Phi,3));
    for n = 1:size(Phi,3)
        PhiEnd = Gr\reshape(Phi(:,end,n,respIdx,end),L,[]);
        temp  = reshape(Utemp(:,useL)*PhiEnd(useL,:),Nydisp,Nxdisp,Nz);
        recon(:,:,:,n) = reshape(Utemp(:,useL)*PhiDisplay(useL,n),Nydisp,Nxdisp,Nz,[]).*exp(-1j*angle(temp));
    end
elseif dim == 2
    PhiDisplay = Gr\reshape(Phi(:,1:ds:end,cardiacIdx,respIdx,:),L,[]);
    PhiEnd = Gr\reshape(Phi(:,end,cardiacIdx,respIdx,end),L,[]);
    temp  = reshape(Utemp(:,useL)*PhiEnd(useL,:),Nydisp,Nxdisp,Nz);
    recon = reshape(Utemp(:,useL)*PhiDisplay(useL,:),Nydisp,Nxdisp,Nz,[]).*exp(-1j*angle(temp));
end
if flagBB 
    if sign(dataArray.bIdx) == -1 && dim ~= 2
        recon = -recon;
    end
    recon = real(recon);  
elseif flagFullSegment
    recon = real(recon);
else
    recon = abs(recon);  
end
recon = recon/max(recon(:));

% Store recon images
reconstructedImages.reconTensor = recon;


