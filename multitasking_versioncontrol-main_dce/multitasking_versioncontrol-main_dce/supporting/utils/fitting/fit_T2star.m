function fitResult = fit_T2star(fitParams)

%% initialize fitting parameters
% slice parameters
fitSlice = 1:fitParams.Nz;
B1_coeff = ones(1,96);

% BET config
BETfrac = 0.35;
BETgrad = 0;

% load fitting parameters, overwrite initial values
extractVarFromStruct(fitParams);

% recon parameters
Norig = floor(min(Norig,Ny)/2);
Necho = size(Phi,6);
L     = size(Gr,1);
tempNecho = size(Phi,1)/L;

Phi = reshape(permute(Phi,[1 2 5 3 4 6]),L*tempNecho,Nseg*moduleLength,size(Phi,3),size(Phi,4),Necho);

for n = numel(fitSlice):-1:1
    if fitSlice(n) > Nzdisp || fitSlice(n) < 1
        fitSlice(n) = [];
    end
end

[Ny,Nx,Nz,~,~] = size(U);
dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), floor(Nz/2)-floor(Nzdisp/2) + (1:Nzdisp), :, :);
Utemp = reshape(dispim(U),Nydisp,Nxdisp,Nzdisp,tempNecho,L);
Utemp = Utemp(:,:,fitSlice,:,:);

% interpolate 2D image to isotropic pixels
rawVoxelSpacing = params.rawVoxelSpacing;
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

% slice profile
if Nz > 1 && numel(B1_coeff) ~= Nz
    halfSlice = numel(B1_coeff)/Nz/2;
    B1_coeff  = interp1((1:numel(B1_coeff))-0.5,B1_coeff,linspace(halfSlice,numel(B1_coeff)-halfSlice,Nz),'pchip');
elseif Nz == 1
    B1_coeff = 1;
end

% inline functions
row    = @(x) x(:).';

%% Choose image time points

totaln    = length(floor(3*Nseg/4):Nseg);
Nsegnew   = ones(1,totaln*moduleLength);
fitw_temp = ones(1,totaln*moduleLength);
for ii = 1: moduleLength
    temp = floor(3*Nseg/4)+Nseg*(ii-1):Nseg*ii;
    Nsegnew(totaln*(ii-1)+(1:totaln)) = temp(1:totaln);
end

if exist('fitw_full','var')
    if numel(fitw_full) == Nseg*moduleLength
        fitw = fitw_full(Nsegnew).*fitw_temp;
    elseif numel(fitw_full) == Nseg
        fitw = fitw_full(imstart:imskip:Nseg).*fitw_temp;
    end
else
    fitw = fitw_temp;
end

%% Setup fitparams
opts = [];
opts.MaxFunEvals = 1000;
opts.Display = 'off';

rphase = max(min(rphase,size(Phi,4)),1);
cphases(cphases>size(Phi,3)) = [];
if isempty(cphases)
    cphases = 1;
end
% store fitParams in fitResult
fitParams.rphase      = rphase;
fitParams.cphases     = cphases;
fitParams.fitw        = fitw;
fitParams.U           = Utemp;
fitParams.Phi         = Phi;
fitParams.fitSlice    = fitSlice;
fitResult.fitParams   = fitParams;

% Initialize fitResult
fitResult.T2starmap    = zeros(Nydisp,Nxdisp,numel(fitSlice),numel(cphases));
fitResult.ResT2starmap = zeros(Nydisp,Nxdisp,numel(fitSlice),numel(cphases));
fitResult.mask         = zeros(Nydisp,Nxdisp,numel(fitSlice),numel(cphases));
fitResult.T2map        = zeros(Nydisp,Nxdisp,numel(fitSlice),numel(cphases));

normTEarray = TEarray - TEarray(1);
normT2prepDuration = TEs - TEs(1);
deltaTE = TEarray(1);

%% Do fitting
% Fit for each cardiac phase
for cphase = 1:numel(cphases)
    fprintf('Resp phase: %d, Cardiac phase: %d / %d, ', rphase, cphase, numel(cphases));
    
    if flagLargeROI == 4 && ~(exist('mask','var') && numel(mask) == Nydisp*Nxdisp*numel(fitSlice))
        size3 = size(Utemp(:,:,:,1,1)); size3(end+1:3) = 1;
        maskBrain = BET(abs(Utemp(:,:,:,1,1)),size3,params.rawVoxelSpacing,BETfrac,BETgrad);
    end

    % generate images
    for sl = 1:numel(fitSlice)
        reconsl = zeros(Nydisp,Nxdisp,numel(Nsegnew),Necho);
        fprintf('slice: %d / %d, ', sl, numel(fitSlice));
        for echo = 1:Necho
%             PhiT1 = Gr\reshape(Phi(echo:Necho:end,Nsegnew,cphases(cphase),rphase,echo),L,[]);
            if tempNecho == Necho
                PhiT2star = Gr\reshape(Phi(echo:tempNecho:end,Nsegnew,cphases(cphase),rphase,echo),L,[]);
                temp = reshape(Utemp(:,:,sl,echo,:),[],L);
            else
                PhiT2star = Gr\reshape(Phi(:,Nsegnew,cphases(cphase),rphase,echo),L,[]);
                temp = reshape(Utemp(:,:,sl,1,:),[],L);
            end
            reconsl(:,:,:,echo) = reshape(temp*PhiT2star,Nydisp,Nxdisp,[],1);
        end
        
        if tempNecho == Necho
            PhiT2 = Gr\reshape(Phi(echo:tempNecho:end,1:Nseg:end,cphases(cphase),rphase,2),L,[]);
            temp = reshape(Utemp(:,:,sl,2,:),[],L);
        else
            PhiT2 = Gr\reshape(Phi(:,1:Nseg:end,cphases(cphase),rphase,2),L,[]);
            temp = reshape(Utemp(:,:,sl,1,:),[],L);
        end
        reconslT2 =  reshape(temp*PhiT2,Nydisp,Nxdisp,[],1);

        % image mask for current slice
        if exist('mask','var') && numel(mask) == Nydisp*Nxdisp*numel(fitSlice)
            im_mask = mask(:,:,sl,cphase);
            im_mask = logical(im_mask);
        else
            if flagLargeROI < 3
                if flagLargeROI == 0
                    Nbins = 4;
                elseif flagLargeROI == 1
                    Nbins = 6;
                elseif flagLargeROI == 2
                    Nbins = 8;
                end
                roi = ones(size(reconsl(:,:,1)));
                im_mask = roi.*abs(mean(reconsl(:,:,totaln:totaln:end),3));
                im_mask(im_mask>prctile(im_mask(:),80)) = prctile(im_mask(:),80);
                im_mask = genMask2D(im_mask,Nbins);
            elseif flagLargeROI == 4
                im_mask = maskBrain(:,:,sl);
            else
                im_mask = ones(size(reconsl(:,:,1)));
            end
            im_mask = logical(im_mask);
            im_mask_temp = imageOrientLPS(im_mask,params);
            figure(11),imshow(im_mask_temp),title('image mask'),drawnow;
            fitResult.mask(:,:,sl,cphase) = im_mask;
        end

        % mask images
        recontemp = reshape(abs(reconsl),size(reconsl,1)*size(reconsl,2),[]);
        recontemp = recontemp(im_mask(:),:);
        recontemp = reshape(recontemp,size(recontemp,1),[],Necho);
        
        recontempT2 = reshape(abs(reconslT2),size(reconslT2,1)*size(reconslT2,2),[]);
        recontempT2 = recontempT2(im_mask(:),:);

        T2star = zeros(size(recontemp,1),1);
        rss    = zeros(size(recontemp,1),1);
        T2     = zeros(size(recontemp,1),1);

        fprintf('fitting %d voxels... ', size(recontemp,1));
        parfor j = 1:size(recontemp,1)
            curves = squeeze(double(recontemp(j,:,:)));
            curves = curves.*fitw.';
            normcurves = curves(:,1)\curves;
            T2star(j)  = -normTEarray/log(normcurves);
            temp       = curves(:,1)*exp(-normTEarray/T2star(j));
            rss(j)     = sum((temp(:) - curves(:)).^2);

            curves = double(recontempT2(j,:));
            normcurves = curves/curves(1);
            T2(j)  = -(normT2prepDuration)/log(normcurves);
        end
        
        T2star(~isfinite(T2star)) = 0;
        rss(~isfinite(rss)) = 0;
        T2(~isfinite(T2)) = 0;

        tempfit = zeros(size(reconsl,1)*size(reconsl,2), 3);
        tempfit(im_mask(:),1) = T2star;
        tempfit(im_mask(:),2) = rss;
        tempfit(im_mask(:),3) = T2;
        tempfit = reshape(tempfit,size(reconsl,1),size(reconsl,2),[]);
        
        % create maps
        fitResult.T2starmap(:,:,sl,cphase) = tempfit(:,:,1);
        fitResult.RSSmap(:,:,sl,cphase)    = tempfit(:,:,2);
%         fitResult.T2map(:,:,sl,cphase)     = tempfit(:,:,3);

        % display maps
        maxRSSdisplay = prctile(rss,97);     
        if maxRSSdisplay == 0 && max(rss) > 0
            maxRSSdisplay = max(rss);
        elseif maxRSSdisplay == 0 || ~isfinite(maxRSSdisplay)
            maxRSSdisplay = 1;
        end

        figure(31),imagesc(wmedfilt2(imageOrientLPS(fitResult.T2starmap(:,:,sl,cphase)*1000,params)),[0 150]),axis equal tight,colormap(warmmetal);colorbar;title('T2star');
        figure(32),imagesc(wmedfilt2(imageOrientLPS(fitResult.RSSmap(:,:,sl,cphase),params)),[0 maxRSSdisplay]),axis equal tight,colormap(jet);colorbar;title('RSS T2star');
%         figure(33),imagesc(wmedfilt2(imageOrientLPS(fitResult.T2map(:,:,sl,cphase),params)),[0 0.3]),axis equal tight,colormap(warmmetal);colorbar;title('T2');

        drawnow;
        fprintf(' done.\n');
    end    
end


