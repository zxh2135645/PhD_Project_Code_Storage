function fitResult = fit_T1_1FA(fitParams)

% Choose if fit for flip angle
ChooseConst = 1;

% initial values
initT1   = 1;
initBIR  = 1;
initGlobalBeta = 0.6;

% lower & upper bounds
minT1   = 100e-3;
maxT1   = 3;
minBIR  = 0.5;
maxBIR  = 2;
minBeta = 0.3;
maxBeta = 1;

% slice parameters
fitSlice = 1:fitParams.Nz;
B1_coeff = ones(1,96);

% BET config
BETfrac = 0.35;
BETgrad = 0;

% load fitting parameters, overwrite initial values
extractVarFromStruct(fitParams);

% use 1st echo
Necho = size(Phi,6);
echo  = min(2,Necho);
Phi   = Phi(echo:Necho:end,:,:,:,:,echo);
L     = size(Phi,1);

for n = numel(fitSlice):-1:1
    if fitSlice(n) > Nzdisp || fitSlice(n) < 1
        fitSlice(n) = [];
    end
end

[Ny,Nx,Nz,~,~] = size(U);
dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), floor(Nz/2)-floor(Nzdisp/2) + (1:Nzdisp), :, :);
Utemp = reshape(dispim(U),Nydisp,Nxdisp,Nzdisp,Necho,L);
Utemp = permute(Utemp(:,:,fitSlice,echo,:),[1,2,3,5,4]);

% slice profile
if Nz > 1
    halfSlice = numel(B1_coeff)/Nz/2;
    B1_coeff  = interp1((1:numel(B1_coeff))-0.5,B1_coeff,linspace(halfSlice,numel(B1_coeff)-halfSlice,Nz),'pchip');
    B1_coeff  = B1_coeff(floor(Nz/2)-floor(Nzdisp/2) + (1:Nzdisp));
    B1_coeff  = B1_coeff(fitSlice);
else
    B1_coeff = ones(1,Nzdisp);
end

[Nydisp,Nxdisp,Nzdisp,~] = size(Utemp);

% inline functions
row   = @(x) x(:).';
ppinv = @(x,y)(y*x')/norm(x)^2; %fast right-sided pseudoinverse function (for later)


%% Choose image time points

totaln = length(imstart:imskip:Nseg);
Nsegnew = [];
for ii = 1: moduleLength
    temp = imstart+Nseg*(ii-1):imskip:Nseg*ii;
    temp = temp(1:totaln);
    Nsegnew = [Nsegnew temp];
end

fitw = reshape(fitw_full(Nsegnew),1,[]);

%% Setup functions

ns = imstart:imskip:Nseg;
ns = ns(1:totaln);

e1  = @(T1) exp(-TR /T1);
Mss = @(e1,alpha)(1-e1)/(1-cos(alpha)*e1);

Sint = @(A,e1,B,Balpha) A * Mss(e1,Balpha) * (1 - (B+1)*(e1*cos(Balpha)).^(ns-1)) * sin(Balpha);
S    = @(A,T1,B,Beta)   row(Sint(A,e1(T1),B,Beta*alphaArray(1)));


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
fitParams.ns          = ns;
fitParams.fitw        = fitw;
fitParams.ChooseConst = ChooseConst;
fitParams.fitSlice    = fitSlice;
fitParams.initT1      = initT1;
fitParams.initBIR     = initBIR;
fitParams.B1_coeff    = B1_coeff;
fitParams.initGlobalBeta = initGlobalBeta;
fitParams.opts = opts;
fitParams.U    = Utemp;
fitParams.Phi  = Phi;
fitResult.fitParams = fitParams;

%alpha0 = initGlobalBeta * alphaArray(1);

fitResult.Avpmap = zeros(Nydisp,Nxdisp,Nzdisp,numel(cphases));
fitResult.T1map  = zeros(Nydisp,Nxdisp,Nzdisp,numel(cphases));
fitResult.BIRmap = zeros(Nydisp,Nxdisp,Nzdisp,numel(cphases));
fitResult.B1map  = zeros(Nydisp,Nxdisp,Nzdisp,numel(cphases));
fitResult.Resmap = zeros(Nydisp,Nxdisp,Nzdisp,numel(cphases));
fitResult.mask   = zeros(Nydisp,Nxdisp,Nzdisp,numel(cphases));

%% Fit for each cardiac phase

for cphase = 1:numel(cphases)
    fprintf('Resp phase: %d, Cardiac phase: %d / %d, ', rphase, cphase, numel(cphases));

    % generate mask
    volume = abs(Utemp(:,:,fitSlice,1));
    volume(volume>prctile(volume(:),80)) = prctile(volume(:),80);

    if flagLargeROI > 4 || flagLargeROI < 0
        flagLargeROI = 0;
    end
    if flagLargeROI == 4
        size3 = size(volume);size3(end+1:3) = 1;
        maskBrain = BET(volume,size3,params.rawVoxelSpacing,BETfrac,BETgrad);
        im_mask = imfill(maskBrain,4);
    elseif flagLargeROI == 3
        im_mask = ones(size(volume));
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
    PhiT1 = Gr\reshape(Phi(:,Nsegnew,cphases(cphase),rphase,:),L,[]);
    for sl = 1:Nzdisp
        fprintf('slice: %d / %d\n', sl, Nzdisp);     
        
        im_mask = logical(fitResult.mask(:,:,sl,cphase));
        im_mask_temp = imageOrientLPS(im_mask,params);
        figure(11),imshow(im_mask_temp),title('image mask'),drawnow;
        
        reconsl = reshape(reshape(Utemp(:,:,sl,:),[],L)*PhiT1,Nydisp,Nxdisp,[]);
        
        recontemp = reshape(reconsl,size(reconsl,1)*size(reconsl,2),[]);
        recontemp = recontemp(im_mask(:),:);
        recontemp = realify(recontemp,'rows');
           
        initBeta = initGlobalBeta * B1_coeff(sl);
        %alpha = B1_coeff(fitSlice(sl)) * alpha0 * ones(size(recontemp,1),1);   

        if ChooseConst == 0
            xlb = [minT1, minBIR, minBeta];
            xub = [maxT1, maxBIR, maxBeta];       
        else
            xlb = double([minT1, minBIR, initBeta]);
            xub = double([maxT1, maxBIR, initBeta]);
        end
        
        fitmat = zeros(size(recontemp,1),numel(xlb)+2);
        fprintf('fitting %d voxels... ', size(recontemp,1));
        parfor j = 1:size(recontemp,1)
            curve = double(recontemp(j,:));
            normcurve = curve(end);
            curve = curve./normcurve;
            Avp = @(T1,BIR,Beta) ppinv(S(1,T1,BIR,Beta).*fitw,curve.*fitw);      % parameterize solution to A as function of R1,B
            [tempfit, res] = lsqnonlin(@(x)abs(S(Avp(x(1),x(2),x(3)),x(1),x(2),x(3))-curve).*fitw,...
                                       [initT1,median([minBIR, abs(curve(1)), maxBIR]),initBeta],... 
                                       xlb, xub, opts);               
            tempfit = [Avp(tempfit(1),tempfit(2),tempfit(3))*normcurve tempfit res];
            fitmat(j,:) = tempfit;
        end
        tempfit = zeros(size(reconsl,1)*size(reconsl,2),numel(xlb)+2);
        tempfit(:,2) = 32767;
        tempfit(im_mask(:),:) = fitmat;

        fitResult.Avpmap(:,:,sl,cphase) = reshape(tempfit(:,1),size(reconsl,1),size(reconsl,2));
        fitResult.T1map (:,:,sl,cphase) = reshape(tempfit(:,2),size(reconsl,1),size(reconsl,2));
        fitResult.BIRmap(:,:,sl,cphase) = reshape(tempfit(:,3),size(reconsl,1),size(reconsl,2));
        fitResult.B1map (:,:,sl,cphase) = reshape(tempfit(:,4),size(reconsl,1),size(reconsl,2));
        fitResult.Resmap(:,:,sl,cphase) = reshape(tempfit(:,end),size(reconsl,1),size(reconsl,2));
        
        
        figure(21),imagesc((fitResult.T1map(:,:,sl,cphase))*1000,[0 2000]),axis image,colormap(warmmetal);title('T1');
        figure(22),imagesc((fitResult.BIRmap(:,:,sl,cphase)),[0 1.5]),axis image,colormap(parula(256));title('B');
        drawnow;
        fprintf(' done.\n');
    end
end
