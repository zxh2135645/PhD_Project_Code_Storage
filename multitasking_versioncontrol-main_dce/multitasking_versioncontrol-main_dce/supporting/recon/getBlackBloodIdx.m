function bIdx = getBlackBloodIdx(params,reconOptions,temporalBasis,spatialCoeff)
 
extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(temporalBasis);
extractVarFromStruct(spatialCoeff);
    
dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), dispSlice, :);

L   =  size(Phi,1);
useL = 1:L;

Utemp = reshape(U,Ny,Nx,Nz,[]);
if MBfactor == 1
    Utemp = fftshift(Utemp,3);
end
if params.isCartesian
    Utemp = fftshift(Utemp,1);
end
recon = reshape(dispim(Utemp),[],L);

reconRecovery = zeros(Nydisp,Nxdisp,1,linesPerShot);
for n = 1%:size(Phi,3)
    PhiT1 = Gr\reshape(Phi(:,1:linesPerShot,n,1,end),L,[]);
    temp = reshape(recon(:,useL)*PhiT1(useL,:),Nydisp,Nxdisp,1,[]);
    tempend = exp(-1j*angle(temp(:,:,:,end)));
    reconRecovery(:,:,:,:,n) = temp.*tempend;
end
fprintf('\n');
toc;
tIdx = floor(1.5/lEchoSpacing);
h = figure;imagesc(abs(reconRecovery(:,:,1,tIdx,1)));axis equal tight;colormap('gray');title('Draw blood pool ROI')
roiBlood = imellipse;
maskBlood = createMask(roiBlood);
close(h)
[~,bIdx] = find(sum(sum(reconRecovery(:,:,1,:,1).*maskBlood,1),2) < 0, 1, 'last');

fprintf('Black blood tIdx = %d (%f sec after inversion)\n',tIdx,tIdx*lEchoSpacing);

reconCardBlackBlood = reconRecovery(:,:,:,bIdx,:);
if sum(reconCardBlackBlood(:)) > 0
    bIdx = bIdx-10;
    if bIdx < 1
        bIdx = 1;
    end
else
    bIdx = -(bIdx+10);
    if abs(bIdx) > linesPerShot
        bIdx = -linesPerShot;
    end
end


