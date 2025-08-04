function fitResult = fit_CESTLorentzian(fitParams)

% load fitting parameters, overwrite initial values
extractVarFromStruct(fitParams);

% recon parameters
Norig = floor(min(Norig,Ny)/2);
Necho = size(Phi,6);
L     = size(Gr,1);
tempNecho = size(Phi,1)/L;

% Use only one echo for fitting
echo = 1;

Phi   = Phi(min(echo,tempNecho):tempNecho:end,:,:,:,:,echo);
% if size(Phi,2) ~= Nseg*moduleLength
%     Phi = reshape(permute(Phi,[1 2 5 3 4]),L,Nseg*moduleLength,size(Phi,3),[]);
% end
PhiCEST = mean(Phi(:,Nseg*(CESTNumRep-16:CESTNumRep-1)+1,:,:,:),2);

for n = numel(fitSlice):-1:1
    if fitSlice(n) > Nzdisp || fitSlice(n) < 1
        fitSlice(n) = [];
    end
end

[Ny,Nx,Nz,~,~] = size(U);
dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), floor(Nz/2)-floor(Nzdisp/2) + (1:Nzdisp), :, :);
Utemp = reshape(dispim(U),Nydisp,Nxdisp,Nzdisp,tempNecho,L);
Utemp = permute(Utemp(:,:,fitSlice,min(echo,tempNecho),:),[1,2,3,5,4]);
[Nydisp,Nxdisp,Nzdisp,~] = size(Utemp);

% inline functions
row    = @(x) x(:).';

%% 

opts = [];
opts.MaxFunEvals = 1000;
opts.Display = 'off';

x0  = ParamsCEST.x0;
xlb = ParamsCEST.xlb;
xub = ParamsCEST.xub;

% Initialize fitResult
fitResult.imgCEST      = zeros(Nydisp,Nxdisp,Nzdisp,numel(cphases));
fitResult.imgCEST_corr = zeros(Nydisp,Nxdisp,Nzdisp,numel(cphases));
fitResult.B0map        = zeros(Nydisp,Nxdisp,Nzdisp,numel(cphases));
fitResult.MTRasym      = zeros(Nydisp,Nxdisp,Nzdisp,numel(cphases));
fitResult.MTRrex       = zeros(Nydisp,Nxdisp,Nzdisp,numel(cphases));
fitResult.Resmap       = zeros(Nydisp,Nxdisp,Nzdisp,numel(cphases));
fitResult.mask         = zeros(Nydisp,Nxdisp,Nzdisp,numel(cphases));


%% Do fitting
% Fit for each cardiac phase
for cphase = 1:numel(cphases)
    fprintf('Resp phase: %d, Cardiac phase: %d / %d, ', rphase, cphase, numel(cphases));
    
    % generate mask
    volume = abs(Utemp(:,:,:,1));
    volume(volume>prctile(volume(:),80)) = prctile(volume(:),80);

    if flagLargeROI > 4 || flagLargeROI < 0
        flagLargeROI = 0;
    end
    if flagLargeROI == 4
        size3 = size(volume);size3(end+1:3) = 1;
        maskBrain = BET(volume,size3,params.rawVoxelSpacing,BETfrac,BETgrad);
        im_mask = imfill(maskBrain,4);
    elseif flagLargeROI == 3
        im_mask = ones(size(recon));
    else
        if flagLargeROI == 0
            Nbins = 4;
        elseif flagLargeROI == 1
            Nbins = 6;
        elseif flagLargeROI == 2
            Nbins = 8;
        end
        roi = ones(size(volume));
        im_mask = roi.*volume;
        im_mask = genMask3D(im_mask,Nbins);
    end
    if abs(params.zDir) == 1
        masktemp = ones(Nydisp+2,Nxdisp,Nzdisp);
        masktemp(2:end-1,:,:) = im_mask;
        masktemp = imfill(masktemp,4);
        masktemp = masktemp(2:end-1,:,:);
    elseif abs(params.zDir) == 2
        masktemp = ones(Nydisp,Nxdisp+2,Nzdisp);
        masktemp(:,2:end-1,:) = im_mask;
        masktemp = imfill(masktemp,4);
        masktemp = masktemp(:,2:end-1,:);
    elseif abs(params.zDir) == 3
        masktemp = ones(Nydisp,Nxdisp,Nzdisp+2);
        masktemp(:,:,2:end-1) = im_mask;
        masktemp = ipermute(imfill(permute(masktemp,abs([params.zDir,params.xDir,params.yDir])),4),abs([params.zDir,params.xDir,params.yDir]));
        masktemp = masktemp(:,:,2:end-1);
    end
    fitResult.mask(:,:,:,cphase) = logical(masktemp);

    % loop through slices
    Phitemp = Gr\reshape(PhiCEST(:,:,cphases(cphase),rphase,:),L,[]);   
    for sl = 1:Nzdisp
        fprintf('slice: %d / %d, ', sl, Nzdisp);

        reconsl = reshape(reshape(Utemp(:,:,sl,:),[],L)*Phitemp,Nydisp,Nxdisp,[]);
        reconsl = abs(reconsl)./abs(reconsl(:,:,end));
        
        tempmask = hann(Nydisp)*hann(Nxdisp)';
        tempmask(tempmask < 0.4) = 0;
        im_mask = fitResult.mask(:,:,sl,cphase).*tempmask;
        im_mask = logical(im_mask);
        im_mask_temp = imageOrientLPS(im_mask,params);
        figure(11),imshow(im_mask_temp),title('image mask'),drawnow;

        reconsl_corr = reconsl;
        [reconsl_corr(:,:,1:end-1),B0map,zSpectFitRef] = CESTdoB0Correction2Pools(reconsl(:,:,1:end-1),CESTSatFreqOffsetppmList(1:end-1),im_mask,x0,xlb,xub);
    
        fitResult.B0map(:,:,sl) = B0map;
        fitResult.imgCEST(:,:,sl,:) = reconsl;
        fitResult.imgCEST_corr(:,:,sl,:) = reconsl_corr;

        % Lorentzian fitting
        MTR_LD = zSpectFitRef - reconsl_corr(:,:,1:end-1);
        tempZ = reshape(MTR_LD,[],size(MTR_LD,3));
%         tempZ = reshape(reconsl_corr(:,:,1:end-1),[],size(reconsl_corr(:,:,1:end-1),3));
        tempfit = zeros(size(tempZ,1),numel(x0)+3);
        tempZ = tempZ(im_mask(:),1:end-1);
        fitmat = zeros(size(tempZ,1),size(tempfit,2));
        offsets = CESTfreqOffsetppmList(1:end-1);
        parfor j = 1:size(tempZ,1)
            zSpect = tempZ(j,:);
            cost = @(x) (lineshape_Lorentzian(x,zSpect,offsets)-zSpect);
            [fitx, res] = lsqnonlin(cost, x0, xlb, xub, opts);
            [~,MTRasym,MTRrex,~] = lineshape_Lorentzian(fitx,zSpect,offsets);
            fitmat(j,:) = cat(2,row(fitx),res,MTRasym(1),MTRrex(1));
        end
        tempfit(im_mask(:),:) = fitmat;
        tempfit = reshape(tempfit,[size(im_mask) numel(x0)+4]);

        % MTRasym in %
        fitResult.MTRasym(:,:,sl,cphase) = tempfit(:,:,end-1)*100;

        % MTRrex in %
        fitResult.MTRrex(:,:,sl,cphase) = tempfit(:,:,end)*100;

        % display maps
        B0map = imageOrientLPS(abs(fitResult.MTRasym(:,:,sl,cphase)),params);
        figure(21),imagesc(B0map),axis equal tight,colormap(parula);colorbar;title('B0 map (ppm)');
       
        MTRasym = imageOrientLPS(abs((fitResult.B1map(:,:,sl,cphase))),params);
        titleString = sprintf('MTRasym %s [%%]',CESTMetabolite);
        figure(23),imagesc(MTRasym,[0 100]),axis equal tight,colormap(parula);colorbar;title(titleString);

        MTRrex = imageOrientLPS(abs((fitResult.MTRrex(:,:,sl,cphase))),params);
        titleString = sprintf('MTRrex %s [%%]',CESTMetabolite);
        figure(24),imagesc(MTRrex),axis equal tight,colormap(parula);colorbar;title(titleString);

        drawnow;
        fprintf(' done.\n');
    end
end


