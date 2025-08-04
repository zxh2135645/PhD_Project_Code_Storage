function fitResult = fit_T1rhoT2IRVFA(fitParams)

%% initialize fitting parameters
% Choose whether use BT2 in fitting
% 0: fit all; 1: constant BT2; 2: BT2 = BIR, 3: constant BT2&Beta, 4: BT2 = BIR, constant Beta
ChooseConst = 4;    

% initial values
initT1    = 1;
initT2    = 50e-3;
initT1rho = 50e-3;
initBIR   = 1;
initBT2   = 1;
initGlobalBeta = 0.6;

% lower & upper bounds
minT1    = 100e-3;
maxT1    = 3;
minT2    = 30e-3;
maxT2    = 300e-3;
minT1rho = 30e-3;
maxT1rho = 300e-3;
minBIR   = 0.5;
maxBIR   = 1;
minBeta  = 0.3;
maxBeta  = 1;
minBT2   = 0.5;
maxBT2   = 1;

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

% Use only one echo for T1/T2 fitting
echo = 1;

Phi   = Phi(min(echo,tempNecho):tempNecho:end,:,:,:,:,echo);
if size(Phi,2) ~= Nseg*moduleLength
    Phi = reshape(permute(Phi,[1 2 5 3 4]),L,Nseg*moduleLength,size(Phi,3),[]);
end

for n = numel(fitSlice):-1:1
    if fitSlice(n) > Nzdisp || fitSlice(n) < 1
        fitSlice(n) = [];
    end
end

[Ny,Nx,Nz,~,~] = size(U);
dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), floor(Nz/2)-floor(Nzdisp/2) + (1:Nzdisp), :, :);
Utemp = reshape(dispim(U),Nydisp,Nxdisp,Nzdisp,tempNecho,L);
Utemp = permute(Utemp(:,:,fitSlice,min(echo,tempNecho),:),[1,2,3,5,4]);

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

% check ChooseConst
if ChooseConst > 5
    ChooseConst = 5;
end
if (numFA > 1) && (ChooseConst > 2) && ~exist('B1map','var')
    ChooseConst = ChooseConst - 2;
end

% No T2 prep
if numT2prep == 0
    minT2 = initT2;
    maxT2 = initT2;
end
% No T1rho prep
if numT1rhoPrep == 0
    minT1rho = initT1rho;
    maxT1rho = initT1rho;
end
% constant BT2
if (ChooseConst == 1) || (ChooseConst == 3)
    minBT2 = initBT2;
    maxBT2 = initBT2;
end
% BT2 = BIR
if (ChooseConst == 2) || (ChooseConst == 4)
    initBT2 = initBIR;
    minBT2  = minBIR;
    maxBT2  = maxBIR;
end

% inline functions
row    = @(x) x(:).';

%% Choose image time points

totaln    = length(imstart:imskip:Nseg);
Nsegnew   = ones(1,totaln*moduleLength);
fitw_temp = ones(1,totaln*moduleLength);
for ii = 1: moduleLength
    temp = imstart+Nseg*(ii-1):imskip:Nseg*ii;
    Nsegnew(totaln*(ii-1)+(1:totaln)) = temp(1:totaln);
%     if (TEs(ii) == 0 && TSLs(ii) == 0 && moduleLength > 1)
%         fitw_temp(totaln*(ii-1)+(1:totaln)) = 0;
%     end
end

if exist('fitw_full','var') && numel(fitw_full) == Nseg*moduleLength
    fitw = fitw_full(Nsegnew).*fitw_temp;
elseif exist('fitw_full','var') && numel(fitw_full) == Nseg
    fitw = fitw_full(mod(Nsegnew-1,Nseg)+1).*fitw_temp;
else
    fitw = fitw_temp;
end

ns = imstart:imskip:Nseg;
ns = ns(1:totaln);

%% Signal equation

e1    = @(T1) exp(-TR /T1);
e2    = @(T2) exp(-TEs/T2);
e1rho = @(T1rho)exp(-TSLs/T1rho);
Mss      = @(e1,alphas) (1-e1) ./ (1-cos(alphas)*e1);
step     = @(e1,alphas) (bsxfun(@power, e1*cos(alphas).', (ns-1))).';
sin_step = @(alphas)    sin(alphas);

Sint = @(A,e1,e2,e1rho,BalphaArray,BIR,BT2,Eff) A .* Mss(e1,BalphaArray) .* (1 + (step(e1,BalphaArray)) .* ((cos(BIR*pi).*IRs + (invSign.*sin(BT2*pi/2)^2).*(T2s+T1rhos).*e2.*e1rho + cos(BT2*pi/2)^2.*(T2s+T1rhos)).*Eff-1)) .* sin_step(BalphaArray);
%Sint = @(A,e1,e2,e1rho,BalphaArray,BIR,BT2,Eff) A .* Mss(e1,BalphaArray) .* (1 + (step(e1,BalphaArray)) .* ((repmat([cos(BIR*pi)*ones(1,numIRns) (-sin(BT2*pi/2)^2)*ones(1,numT2prep+numT1rhoPrep)],1,rep_VE).*e2.*e1rho + repmat([zeros(1,numIRns) cos(BT2*pi/2)^2*ones(1,numT2prep+numT1rhoPrep)],1,rep_VE)).*Eff-1)) .* sin_step(BalphaArray);

if Nz > 1 && MBfactor == 1
    S = @(A,T1,T2,T1rho,Beta,BIR,BT2) row(Sint(  A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray, BIR, BT2, Eff_cellfunc(T1,T2,T1rho,Beta*alphaArray,BIR,BT2,TR,TEs,TSLs,Nseg,invSign)));
else
    S = @(A,T1,T2,T1rho,Beta,BIR,BT2) row(Sint(  A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.9387, BIR, BT2, Eff_cellfunc(T1,T2,T1rho,Beta*alphaArray*.9387,BIR,BT2,TR,TEs,TSLs,Nseg,invSign))...
                                        + Sint(2*A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.5049, BIR, BT2, Eff_cellfunc(T1,T2,T1rho,Beta*alphaArray*.5049,BIR,BT2,TR,TEs,TSLs,Nseg,invSign))...
                                        + Sint(2*A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.0525, BIR, BT2, Eff_cellfunc(T1,T2,T1rho,Beta*alphaArray*.0525,BIR,BT2,TR,TEs,TSLs,Nseg,invSign)));
end

ppinv = @(x,y)(y*x')/norm(x)^2;     % fast right-sided pseudoinverse function (for later)

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
fitParams.U           = Utemp;
fitParams.Phi         = Phi;
fitParams.ChooseConst = ChooseConst;
fitParams.fitSlice    = fitSlice;
fitParams.initT1      = initT1;
fitParams.initBIR     = initBIR;
fitParams.B1_coeff    = B1_coeff;
fitParams.initGlobalBeta = initGlobalBeta;
fitParams.opts = opts;
fitResult.fitParams = fitParams;

% Initialize fitResult
fitResult.Avpmap   = zeros(Nydisp,Nxdisp,Nzdisp,numel(cphases));
fitResult.T1map    = zeros(Nydisp,Nxdisp,Nzdisp,numel(cphases));
fitResult.T2map    = zeros(Nydisp,Nxdisp,Nzdisp,numel(cphases));
fitResult.T1rhomap = zeros(Nydisp,Nxdisp,Nzdisp,numel(cphases));
fitResult.B1map    = zeros(Nydisp,Nxdisp,Nzdisp,numel(cphases));
fitResult.Resmap   = zeros(Nydisp,Nxdisp,Nzdisp,numel(cphases));
fitResult.BIRmap   = zeros(Nydisp,Nxdisp,Nzdisp,numel(cphases));
fitResult.BT2map   = zeros(Nydisp,Nxdisp,Nzdisp,numel(cphases));
fitResult.mask     = zeros(Nydisp,Nxdisp,Nzdisp,numel(cphases));


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
        fprintf('slice: %d / %d, ', sl, Nzdisp);

        im_mask = logical(fitResult.mask(:,:,sl,cphase));
        im_mask_temp = imageOrientLPS(im_mask,params);
        figure(11),imshow(im_mask_temp),title('image mask'),drawnow;
        
        reconsl = reshape(reshape(Utemp(:,:,sl,:),[],L)*PhiT1,Nydisp,Nxdisp,[]);
        
        % mask images
        recontemp = reshape(reconsl,size(reconsl,1)*size(reconsl,2),[]);
        recontemp = recontemp(im_mask(:),:);
        recontemp = realify(recontemp,'rows');
           
        % Setup fitting parameters
        initBeta = initGlobalBeta*B1_coeff(sl);
        if exist('B1map','var') && size(B1map,3) == Nzdisp && size(B1map,4) == numel(cphases)
            betatemp = reshape(B1map(:,:,sl,cphase),size(reconsl,1)*size(reconsl,2),1);
            betatemp = betatemp(im_mask(:));
        else
            betatemp = initBeta*ones(size(recontemp,1),1);
        end
        betatemp(~isfinite(betatemp)) = mean(betatemp(isfinite(betatemp)));
        betatemp(betatemp <= 0) = 0.001;
        
        % constant beta
        if (ChooseConst == 3) || (ChooseConst == 4)
            minBeta  = initBeta;
            maxBeta  = initBeta;
        end
        
        if ChooseConst == 0
            fitmat = zeros(size(recontemp,1),8);
        else
            fitmat = zeros(size(recontemp,1),7);
        end
        
        % voxel-wise fitting
%         fitmat = [];
        fprintf('fitting %d voxels... ', size(recontemp,1));
        parfor j = 1:size(recontemp,1)
            curve = double(recontemp(j,:));
            normcurve = curve(end);
            curve = curve/normcurve;
            betaj = betatemp(j);
            
            if ChooseConst == 0
                x0  = [initT1, initT2, initT1rho, initBeta, initBIR, initBT2];  % initial values
                xlb = [ minT1,  minT2,  minT1rho,  minBeta,  minBIR,  minBT2];  % lower bounds
                xub = [ maxT1,  maxT2,  maxT1rho,  maxBeta,  maxBIR,  maxBT2];  % upper bounds
            elseif (ChooseConst == 1) || (ChooseConst == 2) || (ChooseConst == 6)               
                x0  = [initT1, initT2, initT1rho, initBeta, initBIR];           % 1: constant BT2
                xlb = [ minT1,  minT2,  minT1rho,  minBeta,  minBIR];           % 2: BT2 = BIR
                xub = [ maxT1,  maxT2,  maxT1rho,  maxBeta,  maxBIR];  
            elseif (ChooseConst == 3) || (ChooseConst == 4)    
                x0  = [initT1, initT2, initT1rho,    betaj, initBIR];           % 3: constant BT2 & beta 
                xlb = [ minT1,  minT2,  minT1rho,    betaj,  minBIR];           % 4: BT2 = BIR, constant beta
                xub = [ maxT1,  maxT2,  maxT1rho,    betaj,  maxBIR];   
            elseif (ChooseConst == 5)   
                x0  = [initT1, initT2, initT1rho,  initBIR];           % 5: BT2 = BIR = beta
                xlb = [ minT1,  minT2,  minT1rho,   minBIR];          
                xub = [ maxT1,  maxT2,  maxT1rho,   maxBIR];   
            end
        
            Avp  = @(T1,T2,T1rho,Beta,BIR,BT2) ppinv(S(1,T1,T2,T1rho,Beta,BIR,BT2).*fitw,curve.*fitw);
            if ChooseConst == 0         % fit all parameters
                cost = @(x)abs(S(Avp(x(1),x(2),x(3),x(4),x(5),x(6)),x(1),x(2),x(3),x(4),x(5),x(6))-curve).*fitw;
                [tempfit, res] = lsqnonlin(cost, x0, xlb, xub, opts);
                tempfit = [Avp(tempfit(1),tempfit(2),tempfit(3),tempfit(4),tempfit(5),tempfit(6))*normcurve, tempfit];
                tempfit(numel(x0)+2) = sqrt(res)*abs(normcurve); 
            elseif ChooseConst == 1 || ChooseConst == 3     % constant BT2
                cost = @(x)abs(S(Avp(x(1),x(2),x(3),x(4),x(5),initBT2),x(1),x(2),x(3),x(4),x(5),initBT2)-curve).*fitw;
                [tempfit, res] = lsqnonlin(cost, x0, xlb, xub, opts);
                tempfit = [Avp(tempfit(1),tempfit(2),tempfit(3),tempfit(4),tempfit(5),initBT2)*normcurve, tempfit];
                tempfit(numel(x0)+2) = sqrt(res)*abs(normcurve); 
            elseif ChooseConst == 2 || ChooseConst == 4     % BT2 = BIR
                cost = @(x)abs(S(Avp(x(1),x(2),x(3),x(4),x(5),x(5)),x(1),x(2),x(3),x(4),x(5),x(5))-curve).*fitw;
                [tempfit, res] = lsqnonlin(cost, x0, xlb, xub, opts);
                tempfit = [Avp(tempfit(1),tempfit(2),tempfit(3),tempfit(4),tempfit(5),tempfit(5))*normcurve, tempfit];
                tempfit(numel(x0)+2) = sqrt(res)*abs(normcurve); 
            elseif ChooseConst == 5
                cost = @(x)abs(S(Avp(x(1),x(2),x(3),x(4),x(4),x(4)),x(1),x(2),x(3),x(4),x(4),x(4))-curve).*fitw;
                [tempfit, res] = lsqnonlin(cost, x0, xlb, xub, opts);
                tempfit = [Avp(tempfit(1),tempfit(2),tempfit(3),tempfit(4),tempfit(4),tempfit(4))*normcurve, tempfit tempfit(4)];
                tempfit(numel(x0)+3) = sqrt(res)*abs(normcurve); 
%             elseif ChooseConst == 3     % constant BT2 & beta
%                 cost = @(x)abs(S(Avp(x(1),x(2),x(3),x(4),x(5), initBT2),x(1),x(2),x(3),x(4),x(5),initBT2)-curve).*fitw;
%                 [tempfit, res] = lsqnonlin(cost, x0, xlb, xub, opts);
%                 tempfit = [Avp(tempfit(1),tempfit(2),tempfit(3),tempfit(4),tempfit(5),initBT2)*normcurve, tempfit];
%                 tempfit(numel(x0)+2) = sqrt(res)*abs(normcurve);
%             elseif ChooseConst == 4     % BT2 = BIR, constant beta
%                 cost = @(x)abs(S(Avp(x(1),x(2),x(3),x(4),x(5),x(5)),x(1),x(2),x(3),x(4),x(5),x(5))-curve).*fitw;
%                 [tempfit, res] = lsqnonlin(cost, x0, xlb, xub, opts);
%                 tempfit = [Avp(tempfit(1),tempfit(2),tempfit(3),tempfit(4),tempfit(5),tempfit(5))*normcurve, tempfit];
%                 tempfit(numel(x0)+2) = sqrt(res)*abs(normcurve);
            end
            fitmat(j,:) = tempfit;
        end
        
        tempfit = zeros(size(reconsl,1)*size(reconsl,2), size(fitmat,2));
        tempfit(:,2) = 32767;
        tempfit(:,3) = 32767;
        tempfit(:,4) = 32767;
        tempfit(im_mask(:),:) = fitmat;
        tempfit = reshape(tempfit,size(reconsl,1),size(reconsl,2),[]);
        
        % create maps
        fitResult.Avpmap(:,:,sl,cphase)   = tempfit(:,:,1);
        fitResult.T1map(:,:,sl,cphase)    = tempfit(:,:,2); 
        fitResult.T2map(:,:,sl,cphase)    = tempfit(:,:,3); 
        fitResult.T1rhomap(:,:,sl,cphase) = tempfit(:,:,4); 
        fitResult.B1map(:,:,sl,cphase)    = tempfit(:,:,5);
        fitResult.BIRmap(:,:,sl,cphase)   = tempfit(:,:,6); 
        fitResult.Resmap(:,:,sl,cphase)   = tempfit(:,:,end);     
        
        if ChooseConst == 0
            fitResult.BT2map(:,:,sl,cphase) = tempfit(:,:,7);
            BT2map = imageOrientLPS(abs((fitResult.BT2map(:,:,sl,cphase))),params);
            figure(25),imagesc(BT2map,[0.5 maxBT2]),axis equal tight,colormap(parula);colorbar;title('T2prep Inversion Efficiency');
        elseif (ChooseConst == 1) || (ChooseConst == 3)
            temp = zeros(size(im_mask));
            temp(im_mask) = initBT2;
            fitResult.BT2map(:,:,sl,cphase) = temp;
        elseif (ChooseConst == 2) || (ChooseConst == 4) || (ChooseConst == 5)
            fitResult.BT2map(:,:,sl,cphase) = fitResult.BIRmap(:,:,sl,cphase);
        end
        
        % display maps
        T1map = imageOrientLPS(abs((fitResult.T1map(:,:,sl,cphase))*1000),params);
        figure(21),imagesc(T1map,[0 maxT1*1000]),axis equal tight,colormap(warmmetal);colorbar;title('T1 (msec)');
        
        if numT2prep > 0
            T2map = imageOrientLPS(abs((fitResult.T2map(:,:,sl,cphase))*1000),params);
            figure(22),imagesc(T2map,[0 maxT2*1000]),axis equal tight,colormap(warmmetal);colorbar;title('T2 (msec)');
        end
        
        if numT1rhoPrep > 0
            T1rhomap = imageOrientLPS(abs((fitResult.T1rhomap(:,:,sl,cphase))*1000),params);
            figure(26),imagesc(T1rhomap,[0 maxT1rho*1000]),axis equal tight,colormap(warmmetal);colorbar;title('T1rho (msec)');
        end

        B1map = imageOrientLPS(abs((fitResult.B1map(:,:,sl,cphase))),params);
        figure(23),imagesc(B1map,[0 maxBeta]),axis equal tight,colormap(viridis);colorbar;title('FA Scaling');

        BIRmap = imageOrientLPS(abs((fitResult.BIRmap(:,:,sl,cphase))),params);
        figure(24),imagesc(BIRmap,[min(minBIR,0.5) maxBIR]),axis equal tight,colormap(parula);colorbar;title('Inversion Efficiency');

        drawnow;
        fprintf(' done.\n');
    end
end


