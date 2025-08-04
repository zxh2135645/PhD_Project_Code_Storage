function reconstructedImages = displayTensorBlackBlood(params,reconOptions,dataArray,temporalBasis,spatialCoeff,reconstructedImages,scale,cIdx,rIdx,slice)

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(temporalBasis);
extractVarFromStruct(spatialCoeff);

if nargin < 10 
    slice = 1;
end
if ~(1<=slice && slice<=Nz)
    slice = 1;
end

if nargin < 9 || isempty(rIdx)
    rIdx = 1;
end
if ~(1<=rIdx && rIdx<=size(Phi,4))
    rIdx = 1;
end

if nargin < 8 || isempty(cIdx)
    if ~isempty(dataArray)
        cIdx = dataArray.diastoleIdx;
    else
        cIdx = 1;
    end
end
if ~(1<=cIdx && cIdx<=size(Phi,3))
    cIdx = 1;
end

if nargin < 7
    scale = 1;
end

    
dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), slice, :);

L   =  size(Phi,1);
useL = 1:L;
tic;
fprintf('Phase correction for %d cardiac phases',size(Phi,3));
Utemp = reshape(U,Ny,Nx,Nz,[]);
if MBfactor == 1
    Utemp = fftshift(Utemp,3);
end
if params.isCartesian
    Utemp = fftshift(Utemp,1);
end
recon = reshape(dispim(Utemp),[],L);

reconRecovery = zeros(Nydisp,Nxdisp,1,linesPerShot);
%PCmaps = zeros(Nydisp,Nxdisp,Nz,linesPerShot,size(Phi,3));
for n = 1 %:size(Phi,3)
    fprintf(', %d',n);
    PhiT1 = Gr\reshape(Phi(:,1:linesPerShot,n,rIdx,:),L,[]);
    temp = reshape(recon(:,useL)*PhiT1(useL,:),Nydisp,Nxdisp,1,[]);
    reconRecovery(:,:,:,:,n) = real(temp.*exp(-1j*angle(temp(:,:,:,end))));
    %PCmaps(:,:,:,:,n) = reshape(reshape(reconRecovery(:,:,:,:,n),[],linesPerShot)*pinv(curvePhi(1:linesPerShot,:).'),Nydisp,Nxdisp,Nz,[]);
end
reconRecovery = reconRecovery/max(abs(reconRecovery(:)));
reconRecovery(reconRecovery<0) = 0;
if flagCommandLine
    implayZoom(abs(reconRecovery(:,:,1,:))*3,20);
end

fprintf('\n');
toc;
tIdx = floor(1.5/lEchoSpacing);
h = figure;imagesc(abs(reconRecovery(:,:,1,tIdx,1)));axis equal tight;colormap('gray');title('Select blood pool')
roiBlood = imellipse;
maskBlood = createMask(roiBlood);
close(h)
[~,tIdx] = find(sum(sum(reconRecovery(:,:,1,:,1).*maskBlood,1),2) > 0, 1);

if isempty(tIdx)
    tIdx = 1;
end
fprintf('Black blood tIdx = %d (%f sec after inversion)\n',tIdx,tIdx*lEchoSpacing);

PhiEnd = Gr\reshape(Phi(:,end,:,rIdx,end),L,[]);
temp = reshape(recon(:,useL)*PhiEnd(useL,:),Nydisp,Nxdisp,1,[]);
PhiCard = Gr\reshape(Phi(:,tIdx,:,rIdx,:),L,[]);
reconCardBlackBlood = reshape(recon(:,useL)*PhiCard(useL,:),Nydisp,Nxdisp,1,[]);
reconCardBlackBlood = real(reconCardBlackBlood.*exp(-1j*angle(temp)));

PhiEnd = Gr\reshape(Phi(:,end,cIdx,:,end),L,[]);
temp = reshape(recon(:,useL)*PhiEnd(useL,:),Nydisp,Nxdisp,1,[]);
PhiResp = Gr\reshape(Phi(:,tIdx,cIdx,:,end),L,[]);
reconRespBlackBlood = reshape(recon(:,useL)*PhiResp(useL,:),Nydisp,Nxdisp,1,[]);
reconRespBlackBlood = real(reconRespBlackBlood.*exp(-1j*angle(temp)));

if sum(reconCardBlackBlood(:)) > 0
    reconCardBlackBlood(reconCardBlackBlood<0) = 0;
    reconRespBlackBlood(reconRespBlackBlood<0) = 0;
else
    reconCardBlackBlood(reconCardBlackBlood>0) = 0;
    reconRespBlackBlood(reconRespBlackBlood>0) = 0;
end

reconCardBlackBlood  = reconCardBlackBlood/prctile(abs(reconCardBlackBlood(:)),99);
reconRespBlackBlood  = reconRespBlackBlood/prctile(abs(reconRespBlackBlood(:)),99);

if flagCommandLine
    implayZoom(abs(reconCardBlackBlood(:,:,1,:))*scale,20);
    implayZoom(abs(reconRespBlackBlood(:,:,1,:))*scale,2);
end

%% save gifs

fps = 20;
Nskip = floor(1/(params.lEchoSpacing*fps));

[~,filePrefix] = fileparts(params.filePath);
try
    saveGif(abs(reconCardBlackBlood(:,:,1,:)),filePath,[filePrefix '_reconCard.gif'],min(16,size(reconCardBlackBlood,4)));
    saveGif(abs(reconRespBlackBlood(:,:,1,:)),filePath,[filePrefix '_reconResp.gif'],2);
catch
    fprintf("Saving images to %s failed. Saving them in current folder.", filePath);
    saveGif(abs(reconCardBlackBlood(:,:,1,:)),'.',[filePrefix '_reconCard.gif'],min(16,size(reconCardBlackBlood,4)));
    saveGif(abs(reconRespBlackBlood(:,:,1,:)),'.',[filePrefix '_reconResp.gif'],2);
end

% imwrite(uint8(255*abs(reshape(reconCardBlackBlood(:,:,slice,:),Nydisp,Nxdisp,1,[]))),[filePath '_reconCard.gif'],'loopcount',inf','delaytime',1/16)
% imwrite(uint8(255*abs(reshape(reconRespBlackBlood(:,:,slice,:),Nydisp,Nxdisp,1,[]))),[filePath '_reconResp.gif'],'loopcount',inf','delaytime',1/2)

% Store recon images
reconstructedImages.tIdxDarkBlood = tIdx;
reconstructedImages.reconCard = reconCardBlackBlood;
reconstructedImages.reconResp = reconRespBlackBlood;
% reconstructedImages.PCmaps = PCmaps;

