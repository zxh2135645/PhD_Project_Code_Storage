function fitResult = fit_T1_VFA(fitParams)

% Choose if fit for flip angle
ChooseConst = 1;

% initial values
initT1   = 1;
initBIR  = 1;
initGlobalBeta = 0.6;

% lower & upper bounds
minT1 = 100e-3;
maxT1 = 3;
minBIR = 0.5;
maxBIR = 2;
minBeta = 0.3;
maxBeta = 1;

% slice parameters
fitSlice = 1:fitParams.Nz;
B1_coeff = ones(1,96);

% load fitting parameters, overwrite initial values
extractVarFromStruct(fitParams);

% recon parameters
Necho = size(Phi,6);
Phi   = Phi(1:Necho:end,:,:,:,:,1);
L     = size(Phi,1);

for n = numel(fitSlice):-1:1
    if fitSlice(n) > Nzdisp || fitSlice(n) < 1
        fitSlice(n) = [];
    end
end

[Ny,Nx,Nz,~,~] = size(U);
dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), floor(Nz/2)-floor(Nzdisp/2) + (1:Nzdisp), :, :);
Utemp = reshape(dispim(U),Nydisp,Nxdisp,Nzdisp,Necho,L);
Utemp = permute(Utemp(:,:,fitSlice,1,:),[1,2,3,5,4]);

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
row    = @(x) x(:).';

%% Choose image time points

totaln = length(imstart:imskip:Nseg);
Nsegnew = zeros(1,totaln*moduleLength);
for ii = 1: moduleLength
    temp = (ii-1)*Nseg + (imstart:imskip:Nseg);
    Nsegnew((ii-1)*totaln+(1:totaln)) = temp;
end

fitw = reshape(fitw_full(Nsegnew),1,[]);

%% Setup functions

ns   = imstart:imskip:Nseg;
ns   = ns(1:totaln);

e1  = @(T1) exp(-TR /T1);
Mss = @(e1,alpha)(1-e1)./(1-cos(alpha).*e1);
Mss_scale = @(e1, alpha, alpha_prev) (1-cos(alpha_prev)*e1) / (1-cos(alpha)*e1);

step     = @(e1,alphas) (bsxfun(@power, e1*cos(alphas).', (ns-1))).';
sin_step = @(alphas)    sin(alphas);

Sint = @(A,e1,BalphaArray,BIR,Eff) A .* Mss(e1,BalphaArray) .* (1 + (step(e1,BalphaArray)) .* (-BIR.*Eff-1)) .* sin_step(BalphaArray);

if Nz > 1 && MBfactor == 1
    S = @(A,T1,Beta,BIR) row(Sint(  A, e1(T1), Beta*alphaArray, BIR, Eff_cellfunc_flow(T1,1,1,Beta*alphaArray,BIR,1,TR,TEs,TSLs,Nseg)));
else
    S = @(A,T1,Beta,BIR) row(Sint(  A, e1(T1), Beta*alphaArray, BIR, Eff_cellfunc_flow(T1,1,1,Beta*alphaArray,BIR,1,TR,TEs,TSLs,Nseg)));
%    S = @(A,T1,Beta,BIR) row(Sint(  A, e1(T1), Beta*alphaArray*.9387, BIR, Eff_cellfunc_flow(T1,1,1,Beta*alphaArray*.9387,BIR,1,TR,TEs,TSLs,Nseg))...
%                           + Sint(2*A, e1(T1), Beta*alphaArray*.5049, BIR, Eff_cellfunc_flow(T1,1,1,Beta*alphaArray*.5049,BIR,1,TR,TEs,TSLs,Nseg))...
%                           + Sint(2*A, e1(T1), Beta*alphaArray*.0525, BIR, Eff_cellfunc_flow(T1,1,1,Beta*alphaArray*.0525,BIR,1,TR,TEs,TSLs,Nseg)));
end

ppinv = @(x,y)(y*x')/norm(x)^2; %fast right-sided pseudoinverse function (for later)


%%

% syms Eff_last R1 BIR Balpha
% Eff_prev = Eff_last;
% 
% for shot = 1:moduleLength
%     shot_prev = mod(shot-2, moduleLength) + 1;
%     if mod(shot-1,Ncontrast) < numIRns   
%         % IR
%         Eff_sym(shot) = Mss_scale(exp(-TR*R1), Balpha*alphaArray(shot), Balpha*alphaArray(shot_prev)) * (1 + ((exp(-TR*R1)*cos(Balpha*alphaArray(shot)))^Nseg).'*(BIR*cos(pi)*Eff_prev-1)); 
%     end
% end
% 
% Eff_prev     = Eff_sym(shot);
% Eff_last     = solve(Eff_last == Eff_sym(end), Eff_last);
% Eff_last     = simplifyFraction(Eff_last,'Expand',true);
% Eff_sym(end) = Eff_last;
% 
% evalstring = 'Eff = @(R1,Balpha,BIR) cat(2';
% for shot = 1:moduleLength
%     evalstring = strcat(evalstring,sprintf(',%s',char(subs(Eff_sym(shot)))));
% end
% evalstring = strcat(evalstring,');');
% eval(evalstring);
% 
% Mss      = @(e1,alphas) (1-e1) ./ (1-cos(alphas)*e1);
% step     = @(e1,alphas) bsxfun(@power, e1*cos(alphas).', (ns-1));
% sin_step = @(alphas)    sin(alphas);
% 
% Sint = @(A,e1,e2,Balpha,BIR,Eff) A .* Mss(e1,Balpha*alphaArray) .* (1 + (step(e1,Balpha*alphaArray).') .* ((repmat(BIR*cos(pi)*ones(1,numIRns),1,rep_VE).*e2).*Eff-1)) .* sin_step(Balpha*alphaArray);
% S1   = @(A,R1,R2,Balpha,BIR)     vec(Sint(A,e1(R1),e2(R2), Balpha, BIR, circshift(Eff(R1,Balpha,BIR), [0, 1])));
% 

%% Setup fitparams
%alpha0 = initGlobalBeta * alphaArray(1);

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
fitResult.fitParams = fitParams;

% Initialize fitResult
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
        
        opts = [];
        opts.MaxFunEvals = 1000;
        opts.Display = 'off';
        
        initBeta = initGlobalBeta * B1_coeff(sl);
        if exist('B1map','var') && size(B1map,3) == Nzdisp && size(B1map,4) == numel(cphases)
            betatemp = reshape(B1map(:,:,sl,cphase),size(reconsl,1)*size(reconsl,2),1);
            betatemp = betatemp(im_mask(:));
        else
            betatemp = initBeta*ones(size(recontemp,1),1);
        end
        betatemp(~isfinite(betatemp)) = mean(betatemp(isfinite(betatemp)));
        betatemp(betatemp <= 0) = 0.001;
        
        fitmat = [];
        fprintf('fitting %d voxels... ', size(recontemp,1));
        parfor j = 1:size(recontemp,1)
            curve = double(recontemp(j,:));
            normcurve = curve(end);
            curve = curve./normcurve;
            betaj = betatemp(j);
            
            xinit = [initT1,initBIR,betaj];
            if ChooseConst == 0
                xlb = [minT1, minBIR, minBeta];
                xub = [maxT1, maxBIR, maxBeta];             
            else
                xlb = [minT1, minBIR, betaj];
                xub = [maxT1, maxBIR, betaj]; 
            end
            
            Avp = @(T1,BIR,Beta) ppinv(S(1,T1,BIR,Beta).*fitw,curve.*fitw);      % parameterize solution to A as function of R1,B
            [tempfit, res] = lsqnonlin(@(x)abs(S(Avp(x(1),x(2),x(3)),x(1),x(2),x(3))-curve).*fitw,...
                                           xinit, xlb, xub, opts); 
            tempfit = [Avp(tempfit(1),tempfit(2),tempfit(3))*normcurve tempfit res];
            fitmat(j,:) = tempfit;
        end
        tempfit = zeros(size(reconsl,1)*size(reconsl,2),5);
        tempfit(:,2) = 32767;
        tempfit(im_mask(:),:) = fitmat;

        fitResult.Avpmap(:,:,sl,cphase) = reshape(tempfit(:,1),size(reconsl,1),size(reconsl,2));
        fitResult.T1map (:,:,sl,cphase) = reshape(tempfit(:,2),size(reconsl,1),size(reconsl,2));
        fitResult.BIRmap(:,:,sl,cphase) = reshape(tempfit(:,3),size(reconsl,1),size(reconsl,2));
        fitResult.B1map(:,:,sl,cphase)  = reshape(tempfit(:,4),size(reconsl,1),size(reconsl,2));
        fitResult.Resmap(:,:,sl,cphase) = reshape(tempfit(:,end),size(reconsl,1),size(reconsl,2));   
        
        figure(21),imagesc(wmedfilt2(fitResult.T1map(:,:,sl,cphase))*1000,[0 3000]),axis image,colormap(warmmetal);title('T1');
        figure(22),imagesc(wmedfilt2(fitResult.BIRmap(:,:,sl,cphase)),[0 1.5]),axis image,colormap(parula(256));title('BIR');
        figure(23),imagesc(wmedfilt2(fitResult.B1map(:,:,sl,cphase)),[0 1.2]),axis image,colormap(parula(256));title('Beta');
        drawnow;
        fprintf(' done.\n');
    end
end

