function fitResult = fit_T1rhoT2IRVFA_dic(fitParams)

%% initialize fitting parameters
% Choose whether use BT2 in fitting
% 0: fit all; 1: constant BT2; 2: BT2 = BIR, 3: constant BT2&Beta, 4: BT2 = BIR, constant Beta
ChooseConst = 4;    

% initial values
initT1   = 1;
initT2   = 50e-3;
initBIR  = 1;
initBT2  = 1;
initGlobalBeta = 0.6;

% lower & upper bounds
minT1 = 100e-3;
maxT1 = 3;
minT2 = 20e-3;
maxT2 = 300e-3;
minBIR = 0.5;
maxBIR = 1.5;
minBeta = 0.3;
maxBeta = 1;
minBT2 = 0.5;
maxBT2 = 1.5;

% slice parameters
fitSlice = 1:fitParams.Nz;
B1_coeff = ones(1,96);

% load fitting parameters, overwrite initial values
extractVarFromStruct(fitParams);

% recon parameters
Norig = floor(min(Norig,Ny)/2);
L     = size(Phi,1);

% cardiac/respiratory phases to fit
if ~isfield(fitParams,'cphases') || max(fitParams.cphases) > size(Phi,3)
    cphases = 1:size(Phi,3);
end
if ~isfield(fitParams,'rphase') || max(fitParams.rphase) > size(Phi,4)
    rphase = 1;
end

% slice profile
halfSlice = numel(B1_coeff)/Nz/2;
B1_coeff  = interp1((1:numel(B1_coeff))-0.5,B1_coeff,linspace(halfSlice,numel(B1_coeff)-halfSlice,Nz),'pchip');

% inline functions
dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), :, :);


%% load dictionary

[curves_compress,basis_compress,fitParams] = genGREdictionary_T1rhoT2IRVFA(fitParams);

norm_curves = vecnorm(curves_compress,2);
curves_compress = curves_compress./norm_curves;

%% Do fitting

% Initialize fitResult
fitResult.T1map  = zeros(Nydisp,Nxdisp,numel(fitSlice),numel(cphases));
fitResult.T2map  = zeros(Nydisp,Nxdisp,numel(fitSlice),numel(cphases));
fitResult.B1map  = zeros(Nydisp,Nxdisp,numel(fitSlice),numel(cphases));
fitResult.BIRmap = zeros(Nydisp,Nxdisp,numel(fitSlice),numel(cphases));

% store fitParams in fitResult
fitResult.fitParams = fitParams;
fitResult.fitParams.ChooseConst = ChooseConst;
fitResult.fitParams.fitSlice    = fitSlice;
fitResult.fitParams.initT1      = initT1;
fitResult.fitParams.initT2      = initT2;
fitResult.fitParams.initBIR     = initBIR;
fitResult.fitParams.initBT2     = initBT2;
fitResult.fitParams.initGlobalBeta = initGlobalBeta;
fitResult.fitParams.B1_coeff       = B1_coeff;

% Fit for each cardiac phase
for cphase = 1:numel(cphases)
    fprintf('Cardiac Motion: %d / %d\n', cphases(cphase), numel(cphases));
    
    % generate images
    Phitemp   = Gr\reshape(Phi(:,end,cphases(cphase),rphase,:),L,[]);
    recontemp = reshape(reshape(U,[],L)*Phitemp,Ny,Nx,Nz,[]);
    PhiT1 = Gr\reshape(Phi(:,Nsegnew,cphases(cphase),rphase,:),L,[]);
    recon = reshape(reshape(U,[],L)*PhiT1,Ny,Nx,Nz,[]);
    recon = real(recon.*exp(-1j*angle(recontemp)));
    
    for sl = 1:numel(fitSlice)
        fprintf('slice: %d / %d', sl, numel(fitSlice));
        reconsl = recon(:,:,fitSlice(sl),:);
        if Ny ~= Nynew || Nx ~= Nxnew
            reconsl = imresize(reconsl,[Nynew Nxnew]);
        end
        reconsl = squeeze(dispim(reconsl));

        % image mask for current slice
        roi = ones(size(reconsl(:,:,1)));
        im_mask = roi.*abs(reconsl(:,:,end));
        [im_mask,centroids] = kmeans(im_mask(:),2);
        im_mask = reshape(im_mask == 1+(centroids(2)>centroids(1)),size(reconsl,1),size(reconsl,2));
        im_mask = im_mask/max(abs(im_mask(:)));
        im_mask = im_mask > .4*graythresh(im_mask);
        if flagLargeROI == 1
            se0 = strel('disk',Norig/8);     
            im_mask(3*Norig/8+(2:Norig/4),3*Norig/8+(2:Norig/4)) = im_mask(3*Norig/8+(2:Norig/4),3*Norig/8+(2:Norig/4)) + double(se0.Neighborhood);
        elseif flagLargeROI == 2
            se1 = strel('disk',Norig/4);     
            im_mask(1*Norig/4+(2:Norig/2),1*Norig/4+(2:Norig/2)) = im_mask(1*Norig/4+(2:Norig/2),1*Norig/4+(2:Norig/2)) + double(se1.Neighborhood);
        end
        im_mask = logical(im_mask);
        figure(11),imshow(im_mask),title('image mask'),drawnow;
        fitResult.mask(:,:,sl,cphase) = im_mask;

        % mask images
        recontemp = reshape(reconsl,size(reconsl,1)*size(reconsl,2),[]);
        recontemp = recontemp(im_mask(:),:);
        recontemp = recontemp*basis_compress;
        %recontemp = realify(recontemp,'rows');
        
        % use a patch based method
        patch = 500;
        fitmat = zeros(size(recontemp,1),4);
        tic

        for j = 1:ceil(size(recontemp,1)/patch)

            if j*patch > size(recontemp,1)
                curve = (double(recontemp((j-1)*patch+1:end,:)));
            else
                curve = (double(recontemp((j-1)*patch+1:j*patch,:)));
            end
            curve = (curve./curve(:,end));
            curve = curve./vecnorm(curve,2);
            
            res = abs((curves_compress*curve'));
            for patch_ind = 1:size(res,2)
                [~,ind] = max(res(:,patch_ind));
                [I1,I2,I3,I4] = ind2sub([numel(fitParams.T1s),numel(fitParams.T2s),numel(fitParams.Betas),numel(fitParams.BIRs)],ind);
                fitmat((j-1)*patch+patch_ind,:) = [fitParams.T1s(I1),fitParams.T2s(I2),fitParams.Betas(I3),fitParams.BIRs(I4)];
            end
        end
        toc
 
        tempfit = zeros(size(reconsl,1)*size(reconsl,2), 4);
        tempfit(im_mask(:),:) = fitmat;
        tempfit = reshape(tempfit,size(reconsl,1),size(reconsl,2),[]);
        
        % create maps
        fitResult.T1map(:,:,sl,cphase)  = tempfit(:,:,1);   %reshape(tempfit(:,1),size(reconsl,1),size(reconsl,2));
        fitResult.T2map(:,:,sl,cphase)  = tempfit(:,:,2);   %reshape(tempfit(:,2),size(reconsl,1),size(reconsl,2));
        fitResult.B1map(:,:,sl,cphase)  = tempfit(:,:,3); 
        fitResult.BIRmap(:,:,sl,cphase) = tempfit(:,:,4); 
        
        % display maps
        figure(21),imagesc(wmedfilt2(fitResult.T1map(:,:,sl,cphase))*1000,[0 3000]),axis image,colormap(T1colormap);title('T1');
        figure(22),imagesc(wmedfilt2(fitResult.T2map(:,:,sl,cphase))*1000,[0 100]), axis image,colormap(T2colormap);title('T2');
        figure(23),imagesc(wmedfilt2(fitResult.BIRmap(:,:,sl,cphase)),[0 1.5]), axis image,colormap(parula);title('B');
        figure(25),imagesc(wmedfilt2(fitResult.B1map(:,:,sl,cphase)),[0 0.8]),axis image,colormap(viridis);title('beta');
        drawnow;
    end
end

end

function n = vecnorm(x,dim)
% calculate l2 norm for 2D array
if dim == 1 
    for i=1:size(x,2)
        n(i) = norm(x(:,i));
    end
else
    for i =1:size(x,1)
        n(i) = norm(x(i,:));
    end
    n = n(:);
end
        
end
