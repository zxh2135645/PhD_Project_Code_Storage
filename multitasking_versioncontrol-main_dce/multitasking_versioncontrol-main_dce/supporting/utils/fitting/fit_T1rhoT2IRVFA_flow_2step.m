function fitResult = fit_T1rhoT2IRVFA_flow_2step(fitParams)

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

% load fitting parameters, overwrite initial values
extractVarFromStruct(fitParams);

% recon parameters
Norig = floor(min(Norig,Ny)/2);
Necho = size(Phi,6);
Phi   = Phi(1:Necho:end,:,:,:,:,1);
L     = size(Phi,1);

% slice profile
if Nz > 1 && numel(B1_coeff) ~= Nz
    halfSlice = numel(B1_coeff)/Nz/2;
    B1_coeff  = interp1((1:numel(B1_coeff))-0.5,B1_coeff,linspace(halfSlice,numel(B1_coeff)-halfSlice,Nz),'pchip');
elseif Nz == 1
    B1_coeff = 1;
end

% check ChooseConst
if ChooseConst > 4
    ChooseConst = 4;
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

% interpolate image    
[Nynew,Nxnew,~] = size(U);
Nydisp = floor(interpFactor(1)*Nydisp);
Nxdisp = floor(interpFactor(2)*Nxdisp);
dispim = @(x) x(floor(Nynew/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nxnew/2)-floor(Nxdisp/2) + (1:Nxdisp), fitSlice, :, :);

Utemp = dispim(U(:,:,:,1,:));

% inline functions
row    = @(x) x(:).';

%% Choose image time points

totaln    = length(imstart:imskip:Nseg);
Nsegnew   = ones(1,totaln*moduleLength);
fitw_temp = ones(1,totaln*moduleLength);
minFA = min(alphaArray);
for ii = 1: moduleLength
    temp = imstart+Nseg*(ii-1):imskip:Nseg*ii;
    Nsegnew(totaln*(ii-1)+(1:totaln)) = temp(1:totaln);
%     if (TEs(ii) == 0 && TSLs(ii) == 0 && moduleLength > 1)
%         fitw_temp(totaln*(ii-1)+(1:totaln)) = 0;
%     end
    if alphaArray(ii) > minFA
        fitw_temp(totaln*(ii-1)+(1:totaln)) = 0;
    end
end

fitw  = fitw_full(Nsegnew);
fitw2 = fitw_full(Nsegnew).*fitw_temp;
%fitw(:) = 1;

ns = imstart:imskip:Nseg;
ns = ns(1:totaln);

%% Signal equation

e1    = @(T1) exp(-TR /T1);
e2    = @(T2) exp(-TEs/T2);
e1rho = @(T1rho)exp(-TSLs/T1rho);
Mss      = @(e1,alphas) (1-e1) ./ (1-cos(alphas)*e1);
step     = @(e1,alphas) (bsxfun(@power, e1*cos(alphas).', (ns-1))).';
sin_step = @(alphas)    sin(alphas);

Sint = @(A,e1,e2,e1rho,BalphaArray,BIR,BT2,Eff) A .* Mss(e1,BalphaArray) .* (1 + (step(e1,BalphaArray)) .* ((BIR*cos(pi)*IRs + (invSign.*BT2.*sin(pi/2)^2).*(T2s+T1rhos).*e2.*e1rho).*Eff-1)) .* sin_step(BalphaArray);

%Sint = @(A,e1,e2,e1rho,BalphaArray,BIR,BT2,Eff) A .* Mss(e1,BalphaArray) .* (1 + (step(e1,BalphaArray)) .* ((repmat([BIR*cos(pi)*ones(1,numIRns) (-BT2*sin(pi/2)^2)*ones(1,numT2prep+numT1rhoPrep)],1,rep_VE).*e2.*e1rho).*Eff-1)) .* sin_step(BalphaArray);

if Nz > 1 && MBfactor == 1
    S = @(A,T1,T2,T1rho,Beta,BIR,BT2) row(Sint(  A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray, BIR, BT2, Eff_cellfunc_flow(T1,T2,T1rho,Beta*alphaArray,BIR,BT2,TR,TEs,TSLs,Nseg,invSign)));
else
    S = @(A,T1,T2,T1rho,Beta,BIR,BT2) row(Sint(  A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.9387, BIR, BT2, Eff_cellfunc_flow(T1,T2,T1rho,Beta*alphaArray*.9387,BIR,BT2,TR,TEs,TSLs,Nseg,invSign))...
                                        + Sint(2*A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.5049, BIR, BT2, Eff_cellfunc_flow(T1,T2,T1rho,Beta*alphaArray*.5049,BIR,BT2,TR,TEs,TSLs,Nseg,invSign))...
                                        + Sint(2*A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.0525, BIR, BT2, Eff_cellfunc_flow(T1,T2,T1rho,Beta*alphaArray*.0525,BIR,BT2,TR,TEs,TSLs,Nseg,invSign)));
end

ppinv = @(x,y)(y*x')/norm(x)^2;     % fast right-sided pseudoinverse function (for later)


%% Setup fitparams
opts = [];
opts.MaxFunEvals = 1000;
opts.Display = 'off';
        
% store fitParams in fitResult
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
fitResult.Avpmap   = zeros(Nydisp,Nxdisp,numel(fitSlice),numel(cphases));
fitResult.T1map    = zeros(Nydisp,Nxdisp,numel(fitSlice),numel(cphases));
fitResult.T2map    = zeros(Nydisp,Nxdisp,numel(fitSlice),numel(cphases));
fitResult.T1rhomap = zeros(Nydisp,Nxdisp,numel(fitSlice),numel(cphases));
fitResult.B1map    = zeros(Nydisp,Nxdisp,numel(fitSlice),numel(cphases));
fitResult.Resmap   = zeros(Nydisp,Nxdisp,numel(fitSlice),numel(cphases));
fitResult.BIRmap   = zeros(Nydisp,Nxdisp,numel(fitSlice),numel(cphases));
fitResult.BT2map   = zeros(Nydisp,Nxdisp,numel(fitSlice),numel(cphases));
fitResult.mask     = zeros(Nydisp,Nxdisp,numel(fitSlice),numel(cphases));


%% Do fitting

% Fit for each cardiac phase
for cphase = 1:numel(cphases)
    fprintf('Cardiac phase: %d / %d\n', cphase, numel(cphases));
    
    % generate images
    PhiT1 = Gr\reshape(Phi(:,Nsegnew,cphases(cphase),rphase,:),L,[]);
    
    for sl = 1:numel(fitSlice)
        fprintf('slice: %d / %d, ', sl, numel(fitSlice));
        reconsl = reshape(reshape(Utemp(:,:,sl,:),[],L)*PhiT1,Nydisp,Nxdisp,[]);
        
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
%             n = 100;
%             se1 = strel('disk',n);     
%             im_mask(1*(Nydisp-n*2)/2+(2:n*2),1*(Nxdisp-n*2)/2+(2:n*2)) = im_mask(1*(Nydisp-n*2)/2+(2:n*2),1*(Nxdisp-n*2)/2+(2:n*2)) + double(se1.Neighborhood);
        end
        im_mask = logical(im_mask);
        figure(11),imshow(im_mask),title('image mask'),drawnow;
        fitResult.mask(:,:,sl,cphase) = im_mask;
        
        % mask images
        recontemp = reshape(reconsl,size(reconsl,1)*size(reconsl,2),[]);
        recontemp = recontemp(im_mask(:),:);
        recontemp = realify(recontemp,'rows');
           
        % Setup fitting parameters
        initBeta = initGlobalBeta*B1_coeff(fitSlice(sl));
        
        % constant beta
        if (ChooseConst == 3) || (ChooseConst == 4)
            minBeta  = initBeta;
            maxBeta  = initBeta;
        end
        
        if ChooseConst == 0
            x0  = [initT1, initT2, initT1rho, initBeta, initBIR, initBT2];  % initial values
            xlb = [ minT1,  minT2,  minT1rho,  minBeta,  minBIR,  minBT2];  % lower bounds
            xub = [ maxT1,  maxT2,  maxT1rho,  maxBeta,  maxBIR,  maxBT2];  % upper bounds
            fitmat = zeros(size(recontemp,1),8);
        else
            x0  = [initT1, initT2, initT1rho, initBeta, initBIR];           % 1: constant BT2
            xlb = [ minT1,  minT2,  minT1rho,  minBeta,  minBIR];           % 2: BT2 = BIR
            xub = [ maxT1,  maxT2,  maxT1rho,  maxBeta,  maxBIR];  
            fitmat = zeros(size(recontemp,1),7);
        end
        
        % voxel-wise fitting
        fitmat = [];
        fprintf('fitting %d voxels... ', size(recontemp,1));
        parfor j = 1:size(recontemp,1)
            curve = double(recontemp(j,:));
            normcurve = curve(end);
            curve = curve/normcurve;

            Avp  = @(T1,T2,T1rho,Beta,BIR,BT2) ppinv(S(1,T1,T2,T1rho,Beta,BIR,BT2).*fitw,curve.*fitw);
            if ChooseConst == 0         % fit all parameters
                cost = @(x)abs(S(Avp(x(1),x(2),x(3),x(4),x(5),x(6)),x(1),x(2),x(3),x(4),x(5),x(6))-curve).*fitw;
                [tempfit, ~] = lsqnonlin(cost, x0, xlb, xub, opts);
                xinit2 = x0; xinit2(4) = tempfit(4);
                xlb2 = xlb;  xlb2(4)   = tempfit(4);
                xub2 = xub;  xub2(4)   = tempfit(4);
                Avp  = @(T1,T2,T1rho,Beta,BIR,BT2) ppinv(S(1,T1,T2,T1rho,Beta,BIR,BT2).*fitw2,curve.*fitw2);
                cost = @(x)abs(S(Avp(x(1),x(2),x(3),x(4),x(5),x(6)),x(1),x(2),x(3),x(4),x(5),x(6))-curve).*fitw2;
                [tempfit, res] = lsqnonlin(cost, xinit2, xlb2, xub2, opts);
                tempfit = [Avp(tempfit(1),tempfit(2),tempfit(3),tempfit(4),tempfit(5),tempfit(6))*normcurve, tempfit];
                tempfit(numel(x0)+2) = sqrt(res)*abs(normcurve); 
            elseif ChooseConst == 1 || ChooseConst == 3     % constant BT2
                cost = @(x)abs(S(Avp(x(1),x(2),x(3),x(4),x(5),initBT2),x(1),x(2),x(3),x(4),x(5),initBT2)-curve).*fitw;
                [tempfit, ~] = lsqnonlin(cost, x0, xlb, xub, opts);
                xinit2 = x0; xinit2(4) = tempfit(4);
                xlb2 = xlb;  xlb2(4)   = tempfit(4);
                xub2 = xub;  xub2(4)   = tempfit(4);
                Avp  = @(T1,T2,T1rho,Beta,BIR,BT2) ppinv(S(1,T1,T2,T1rho,Beta,BIR,BT2).*fitw2,curve.*fitw2);
                cost = @(x)abs(S(Avp(x(1),x(2),x(3),x(4),x(5),initBT2),x(1),x(2),x(3),x(4),x(5),initBT2)-curve).*fitw2;
                [tempfit, res] = lsqnonlin(cost, xinit2, xlb2, xub2, opts);
                tempfit = [Avp(tempfit(1),tempfit(2),tempfit(3),tempfit(4),tempfit(5),initBT2)*normcurve, tempfit];
                tempfit(numel(x0)+2) = sqrt(res)*abs(normcurve); 
            elseif ChooseConst == 2 || ChooseConst == 4     % BT2 = BIR
                cost = @(x)abs(S(Avp(x(1),x(2),x(3),x(4),x(5),x(5)),x(1),x(2),x(3),x(4),x(5),x(5))-curve).*fitw;
                [tempfit, ~] = lsqnonlin(cost, x0, xlb, xub, opts);
                xinit2 = x0; xinit2(4) = tempfit(4);
                xlb2 = xlb;  xlb2(4)   = tempfit(4);
                xub2 = xub;  xub2(4)   = tempfit(4);
                Avp  = @(T1,T2,T1rho,Beta,BIR,BT2) ppinv(S(1,T1,T2,T1rho,Beta,BIR,BT2).*fitw2,curve.*fitw2);
                cost = @(x)abs(S(Avp(x(1),x(2),x(3),x(4),x(5),x(5)),x(1),x(2),x(3),x(4),x(5),x(5))-curve).*fitw2;
                [tempfit, res] = lsqnonlin(cost, xinit2, xlb2, xub2, opts);
                tempfit = [Avp(tempfit(1),tempfit(2),tempfit(3),tempfit(4),tempfit(5),tempfit(5))*normcurve, tempfit];
                tempfit(numel(x0)+2) = sqrt(res)*abs(normcurve); 
            end
            fitmat(j,:) = tempfit;
        end
        
        tempfit = zeros(size(reconsl,1)*size(reconsl,2), size(fitmat,2)+2);
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
            figure(25),imagesc(fitResult.BT2map(:,:,sl,cphase),[0.5 2]),axis image,colormap(parula);title('BT2');
        elseif (ChooseConst == 1) || (ChooseConst == 3)
            temp = zeros(size(im_mask));
            temp(im_mask) = initBT2;
            fitResult.BT2map(:,:,sl,cphase) = temp;
        elseif (ChooseConst == 2) || (ChooseConst == 4)
            fitResult.BT2map(:,:,sl,cphase) = fitResult.BIRmap(:,:,sl,cphase);
        end
        
        % display maps
        figure(21),imagesc(wmedfilt2(fitResult.T1map(:,:,sl,cphase))*1000,[0 3000]),axis image,colormap(warmmetal);title('T1');
        if numT2prep > 0
            figure(22),imagesc(wmedfilt2(fitResult.T2map(:,:,sl,cphase))*1000,[0 250]),axis image,colormap(warmmetal);title('T2');
        end
        if numT1rhoPrep > 0
            figure(26),imagesc(wmedfilt2(fitResult.T1rhomap(:,:,sl,cphase))*1000,[0 250]),axis image,colormap(warmmetal);title('T1rho');
        end
        figure(23),imagesc(wmedfilt2(fitResult.B1map(:,:,sl,cphase)),[0 1]),axis image,colormap(viridis);title('Beta');
        figure(24),imagesc(wmedfilt2(fitResult.BIRmap(:,:,sl,cphase)),[0.5 2]),axis image,colormap(parula);title('BIR');
        
        drawnow;
        fprintf(' done.\n');
    end
end


