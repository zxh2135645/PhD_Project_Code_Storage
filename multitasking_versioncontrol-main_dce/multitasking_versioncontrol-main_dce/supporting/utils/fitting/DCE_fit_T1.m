function fitResult = DCE_fit_T1(fitParams)

% inline functions
vec   = @(x) x(:); 
row   = @(x) x(:).';
ppinv = @(x,y)(y*x')/norm(x)^2; %fast right-sided pseudoinverse function (for later)
crop = @(x,st)x(1:st.Nd(1),1:st.Nd(2),:,:,:);
    
% Choose if fit for flip angle
% 1: flip angle is constant
% 0: fit for flip angle
ChooseConst = 1;

% initial values
initT1   = 1;
initBIR  = 0;
initGlobalBeta = 1;

% lower & upper bounds
minT1 = 100e-3;
maxT1 = 3;
minBIR = -0.1;
maxBIR = 0.1;
minBeta = 0.3;
maxBeta = 1;

% slice parameters
fitSlice = 1:fitParams.Nz;
B1_coeff = ones(1,96);

% DCE fitting parameters
preCutoff   = 51;
contrastInj = 35;

% load fitting parameters, overwrite initial values
extractVarFromStruct(fitParams);

% recon parameters
Necho = size(Phi,6);
Phi = Phi(1:Necho:end,:,:,:,:,1);
oriLength = size(Phi,5); 
Phi = Phi(:,:,:,:,[preCutoff:end]);
L = size(Phi,1);

[Ny,Nx,Nz,~,~] = size(U);
dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), floor(Nz/2)-floor(Nzdisp/2) + (1:Nzdisp), :, :);
for n = numel(fitSlice):-1:1
    if fitSlice(n) > Nzdisp || fitSlice(n) < 1
        fitSlice(n) = [];
    end
end

Utemp = reshape(dispim(U),Nydisp,Nxdisp,Nzdisp,Necho,L);
Utemp = permute(Utemp(:,:,fitSlice,1,:),[1,2,3,5,4]);

if ~exist('mask','var') || size(mask,1)~= Nydisp || size(mask,2)~= Nxdisp || size(mask,3)~= Nz
    mask = zeros(Nydisp,Nxdisp,numel(fitSlice));
    for n = 1:numel(fitSlice)
        [im_mask,centroids] = kmeans(vec(abs(Utemp(:,:,n,1))),4);
        [~,air] = min(centroids);
        im_mask = reshape(im_mask,Nydisp,Nxdisp);
        im_mask = (im_mask~=air);
        mask(:,:,n) = im_mask;
    end
    mask = logical(mask);
else
    mask = mask(:,:,fitSlice);
end

% slice profile
if Nz > 1
    halfSlice = numel(B1_coeff)/Nz/2;
    B1_coeff  = interp1((1:numel(B1_coeff))-0.5,B1_coeff,linspace(halfSlice,numel(B1_coeff)-halfSlice,Nz),'pchip');
else
    B1_coeff = 1;
end

    
%% Choose image time points

totaln = length(imstart:imskip:Nseg);
Nsegnew = [];
for ii = 1: moduleLength
    temp = imstart+Nseg*(ii-1):imskip:Nseg*ii;
    temp = temp(1:totaln);
    Nsegnew = [Nsegnew temp];
end

fitw = fitw_full(1,Nsegnew,preCutoff:end);
fitw(:,:,end+1) = fitw(:,:,end);
fitw(:,:,end+1) = fitw(:,:,end);

%% Signal model
ns = vec(1:Nseg);

e1  = @(T1) exp(-TR./T1);
Mss = @(e1,alpha) (1-e1)./(1-cos(alpha).*e1);

Sint = @(A,e1,B,Balpha) vec(A .* Mss(e1,Balpha) .* (1-(B+1)*(e1*cos(Balpha)).^(vec(ns)-1)) .* sin(Balpha));
S    = @(A,T1,B,Beta)   row(Sint(A,e1(T1),B,Beta*alphaArray(1)));


%% Setup fitparams
opts = [];
opts.MaxFunEvals = 1000;
opts.Display = 'off';
        
fitParams.ns = ns;
fitParams.fitw = fitw;
fitParams.ChooseConst = ChooseConst;
fitParams.fitSlice    = fitSlice;
fitParams.initT1      = initT1;
fitParams.initBIR     = initBIR;
fitParams.B1_coeff    = B1_coeff;
fitParams.initGlobalBeta = initGlobalBeta;
fitParams.opts = opts;

fitResult.fitParams = fitParams;

%% initial fitting


% Set initial fitting DCE indices
wallClock = cat(2, 1:10:contrastInj, (contrastInj+1):(contrastInj+60), (contrastInj+61):10:size(Phi,5));
wallClock(wallClock > size(Phi,5)) = [];
wallClock_count = numel(wallClock);

% Generate images
PhiDCE  = Gr\realify(reshape(Phi(:,Nsegnew,1,rphase,wallClock),L,[]),'rows');   
Utemp   = reshape(Utemp,[],L);
curves  = Utemp(vec(mask),:)*PhiDCE;
curves  = realify(curves,'rows');
num_all = size(curves,1);

% Fitting parameters
if ChooseConst == 1
    minBeta = initGlobalBeta;
    maxBeta = initGlobalBeta;
end
xinit = double([initT1*ones(1,wallClock_count), initBIR, initGlobalBeta]);
xlb = double([minT1*ones(1,wallClock_count), minBIR, minBeta]);
xub = double([maxT1*ones(1,wallClock_count), maxBIR, maxBeta]);

% fitting weight
fitw_wallClock = row(fitw(1,:,wallClock));

% Do fitting
tic; fprintf('Initial coarse fitting for %d voxels, %d DCE time points... ', num_all, wallClock_count);
fits_init = zeros(num_all,wallClock_count+2);
parfor j = 1:num_all
    curve = double(curves(j,:));
    normcurve = curve(end);
    curve = curve/normcurve;

    Avp = @(T1,BIR,Beta) ppinv(S(1,T1,BIR,Beta).*fitw_wallClock,curve.*fitw_wallClock);      % parameterize solution to A as function of R1,B
    [tempfit, ~] = lsqnonlin(@(x) abs(S(Avp(x(1:wallClock_count),x(wallClock_count+1),x(wallClock_count+2)),x(1:wallClock_count),x(wallClock_count+1),x(wallClock_count+2))-curve).*fitw_wallClock,...
                                    xinit, xlb, xub, opts);               

    fits_init(j,:) = tempfit;
end
fprintf('Done. '); toc;

%% Second fitting use previous T1 result as initial value
% New DCE indices
wallClocknew = cat(2, 1:5:contrastInj, (contrastInj+1):(contrastInj+60), (contrastInj+61):3:size(Phi,5));
wallClocknew(wallClocknew > size(Phi,5)) = [];
wallClock_count = numel(wallClocknew);

% Generate images
Phitemp = zeros(L,size(Phi,2),size(Phi,3),size(Phi,4),numel(wallClocknew));
for i = 1:wallClock_count-1
    Phitemp(:,:,:,:,i) = mean(Phi(:,:,:,:,wallClocknew(i):wallClocknew(i+1)-1),5);
end
Phitemp(:,:,:,:,end) = mean(Phi(:,:,:,:,wallClocknew(i+1):end),5);
PhiDCE  = Gr\reshape(Phitemp(:,:,1,rphase,:),L,[]);
curves  = Utemp(vec(mask),:)*PhiDCE;
curves  = realify(curves,'rows');
num_all = size(curves,1);

% Interpolate initial T1 values
X  = vec(wallClock);
Xq = vec(wallClocknew);
V  = reshape(fits_init(:,1:end-2),num_all,[]).';
Vq = interp1(X, V, Xq, 'linear', 'extrap').';
xinit = cat(2, Vq, fits_init(:,(end-1):end)); % As initial value for next fitting
xlb = double([minT1*ones(1,wallClock_count), minBIR, minBeta]);
xub = double([maxT1*ones(1,wallClock_count), maxBIR, maxBeta]);

% fitting weight
fitw_wallClock = ones(1,Nseg,wallClock_count);
for i = 1:wallClock_count-1
    fitw_wallClock(1,:,i) = mean(fitw(:,:,wallClocknew(i):wallClocknew(i+1)-1),3);
end
fitw_wallClock(:,:,end) = mean(fitw(:,:,wallClocknew(i+1):end),3);
fitw_wallClock = row(fitw_wallClock);

% Do fitting 
tic; fprintf('Second Fitting for %d voxels, %d DCE time points... ', num_all, wallClock_count);
fits = zeros(num_all,wallClock_count+4);
parfor j = 1:num_all
    curve = double(curves(j,:));
    normcurve = curve(end);
    curve = curve/normcurve;
    
    Avp = @(T1,BIR,Beta) ppinv(S(1,T1,BIR,Beta).*fitw_wallClock,curve.*fitw_wallClock);      % parameterize solution to A as function of R1,B
    [tempfit, res] = lsqnonlin(@(x) abs(S(Avp(x(1:wallClock_count),x(wallClock_count+1),x(wallClock_count+2)),x(1:wallClock_count),x(wallClock_count+1),x(wallClock_count+2))-curve).*fitw_wallClock,...
                                    xinit(j,:), xlb, xub, opts);       

    fits(j,:) = [Avp(tempfit(1:wallClock_count),tempfit(wallClock_count+1),tempfit(wallClock_count+2))*normcurve tempfit res];
end
toc; fprintf('Done.\n');

% Interpolate into full DCE temporal resolution (= 2xTR)
X  = vec(wallClocknew);
Xq = vec(1:size(Phi,5));
V  = reshape(fits(:,2:end-3),num_all,[]).';
Vq = interp1(X, V, Xq, 'linear', 'extrap').';
fits_new = zeros(num_all,size(Phi,5)+4);
fits_new(:,1) = fits(:,1); 
fits_new(:,2:end-3) = reshape(Vq,num_all,[]);
fits_new(:,end-2:end) = fits(:,end-2:end);
fits = fits_new;


%% fitResult
fitResult.mask = mask;

% T1 field is a 2D matrix. The 1st dimension corresponds to the voxels in
% 'mask', the 2nd dimension is the T1 value as a function of time
fitResult.T1 = reshape(Vq, num_all,[]);

% The other fields are vectors correspond to each voxel in 'mask'. 
fitResult.Avp  = fits(:,1);     % intensity scalar
fitResult.BIR  = fits(:,end-2); % Saturation efficiency (assumed perfect saturation)
fitResult.Beta = fits(:,end-1);   % flip angle scalar
fitResult.res  = fits(:,end);


