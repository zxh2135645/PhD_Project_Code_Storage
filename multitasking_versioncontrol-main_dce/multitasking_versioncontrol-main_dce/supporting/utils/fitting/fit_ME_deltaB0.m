function fitResult = fit_ME_deltaB0(fitParams)
% inline functions
vec   = @(x) x(:); 

% load fitting parameters, overwrite initial values
extractVarFromStruct(fitParams);

% recon parameters
L     = size(Gr,1);
Necho = size(Phi,6);
tempNecho = size(Phi,1)/L;

Phi = reshape(permute(Phi,[1 2 5 3 4 6]),L*tempNecho,Nseg*moduleLength,size(Phi,3),size(Phi,4),Necho);

for n = numel(fitSlice):-1:1
    if fitSlice(n) > Nzdisp || fitSlice(n) < 1
        fitSlice(n) = [];
    end
end

Nzdisp = numel(fitSlice);
[Ny,Nx,Nz,~,~] = size(U);
dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), fitSlice, :, :);
Utemp = reshape(dispim(U),Nydisp,Nxdisp,Nzdisp,tempNecho,L);


% inline functions
vec = @(x) x(:);
row    = @(x) x(:).';
size2 = @(x) [size(x,1) size(x,2)];
size3 = @(x) [size(x,1) size(x,2) size(x,3)];

%% Choose image time points

rphase = max(min(rphase,size(Phi,4)),1);
cphases(cphases>size(Phi,3)) = [];
if isempty(cphases)
    cphases = 1;
end

% tIdx = Nseg:Nseg:Nseg*moduleLength;
tIdx = Nseg;

reconME = zeros(Nydisp,Nxdisp,Nzdisp,Necho);
for echo = 1:Necho
    if tempNecho == Necho
        PhiME = Gr\reshape(mean(Phi(echo:tempNecho:end,tIdx,cphases(1),rphase,echo),2),L,[]);
        temp = reshape(Utemp(:,:,:,echo,:),[],L);
    else
        PhiME = Gr\reshape(mean(Phi(:,tIdx,cphases(1),rphase,echo),2),L,[]);
        temp = reshape(Utemp(:,:,:,1,:),[],L);
    end
    reconME(:,:,:,echo) = reshape(temp*PhiME,Nydisp,Nxdisp,[],1);
end
reconME = reconME./max(abs(reconME(:)));

if ~exist('mask','var') || size(mask,1)~= Nydisp || size(mask,2)~= Nxdisp || size(mask,3)~= numel(fitSlice)
    mask = genMask2D(reconME(:,:,:,1));
    mask = logical(mask);
else
    mask = logical(mask);
end


%% Set fitting parameters
% Image parameters
imDataParams.FieldStrength = params.MagneticFieldStrength;
imDataParams.PrecessionIsClockwise = 0;  % check this!!!!
imDataParams.TE = TEarray;

% General fitting parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
algoParams.species(2).relAmps   = [0.087 0.693 0.128 0.004 0.039 0.048];

% Algorithm-specific parameters
algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
algoParams.range_r2star = [1 150]; % Range of R2* values, was [0,150]
algoParams.NUM_R2STARS = 76; % Numbre of R2* values for quantization, was 65
algoParams.range_fm = [-450 450]; % Range of field map values
algoParams.NUM_FMS = 301; % Number of field map values to discretize
algoParams.NUM_ITERS = 40; % Number of graph cut iterations
algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
algoParams.lambda = 0.05; % Regularization parameter
algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
algoParams.TRY_PERIODIC_RESIDUAL = 0;
THRESHOLD = 0.01;

%% Recon -- graph cut 
%  Hernando D, Kellman P, Haldar JP, Liang ZP. Robust water/fat separation in the presence of large 
%  field inhomogeneities using a graph cut algorithm. Magn Reson Med. 2010 Jan;63(1):79-90.)

fitResult.mask     = mask;
fitResult.water    = zeros(Nydisp,Nxdisp,numel(fitSlice));
fitResult.fat      = zeros(Nydisp,Nxdisp,numel(fitSlice));
fitResult.delta_B0 = zeros(Nydisp,Nxdisp,numel(fitSlice));
fitResult.ff       = zeros(Nydisp,Nxdisp,numel(fitSlice));
fitResult.T2star   = zeros(Nydisp,Nxdisp,numel(fitSlice));

tic; fprintf('Calculating delta B0: cphase %d, rphase %d', cphases(1), rphase);
for n = 1:numel(fitSlice)
    fprintf(' - slice %d/%d... \n',n,numel(fitSlice));
    imDataParams.images = reshape(double(reconME(:,:,n,:)),Nydisp,Nxdisp,1,1,[]);
  
    outParams = fw_i2cm1i_3pluspoint_hernando_graphcut( imDataParams, algoParams );

    %% Recon -- mixed fit for phase error correction
    % Initialize mixed fitting to graph cut solution
    algoParams.fieldmap  = outParams.fieldmap;
    algoParams.r2starmap = outParams.r2starmap;
    algoParams.NUM_MAGN = 1;
    algoParams.THRESHOLD = 0.04;
    algoParams.range_r2star = [1 150];

    % Do mixed fitting
    % Hernando D, Hines CDG, Yu H, Reeder SB. Addressing phase errors in fat-water imaging 
    % using a mixed magnitude/complex fitting method. Magn Reson Med; 2011.)
    % outParamsMixed = fw_i2xm1c_3pluspoint_hernando_mixedfit( imDataParams, algoParams );

    %% Calculate fat fraction
    water = outParams.species(1).amps;
    fat   = outParams.species(2).amps;
    
    denom = (abs(fat + water));
    denom2 = denom;
    denom2(denom==0) = 1; % To avoid divide-by-zero issues
    % FF = F/(F+W), F>W
    %      1- W/(F+W). F<W
    ff_bias_corr = 100*abs(fat)./denom2;
    ff_bias_corr(abs(fat)<abs(water)) = 100*(1-abs(water(abs(fat)<abs(water)))./denom2(abs(fat)<abs(water)));
    FW_ratio_bias_corr = zeros(size(ff_bias_corr));
    FW_ratio_bias_corr(abs(fat)<abs(water)) = 100*(denom2(abs(fat)<abs(water))./abs(water(abs(fat)<abs(water)))-1);
    FW_ratio_bias_corr(abs(fat)>=abs(water)) = 100*abs(fat(abs(fat)>=abs(water)))./(denom2(abs(fat)>=abs(water))-abs(fat(abs(fat)>=abs(water))));

    denom = (abs(fat) + abs(water));
    denom2 = denom;
    denom2(denom==0) = 1; % To avoid divide-by-zero issues
    ff = 100*abs(fat)./denom2;
    FW_ratio = 100*abs(fat)./abs(water);
    
    delta_B0 = outParams.fieldmap;    
    
    R2s = outParams.r2starmap;
    T2star = 1./R2s;
    T2star(~isfinite(T2star)) = 0;
    
    % display result
    B0map = imageOrientLPS(wmedfilt2(mask(:,:,n).*delta_B0),params);
    figure(53);imagesc(B0map,[-150 150]); colormap('parula'); colorbar; axis equal off;title('delta B_0 (Hz)');

    % store result
    fitResult.water(:,:,n)    = water;
    fitResult.fat(:,:,n)      = fat;
    fitResult.delta_B0(:,:,n) = delta_B0;
    fitResult.ff(:,:,n)       = ff;
%     fitResult.R2star(:,:,n)   = R2s;
    fitResult.T2star(:,:,n)   = T2star; % Saturation efficiency (assumed perfect saturation)

end
toc;


