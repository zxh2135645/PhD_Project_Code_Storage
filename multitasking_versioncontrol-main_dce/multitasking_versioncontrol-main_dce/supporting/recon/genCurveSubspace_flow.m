function [curveU, curveS, fitParams] = genCurveSubspace_flow(params, reconOptions, fitParams)

% License check
if ~checkLicense(reconOptions)
    dlg = errordlg('Multitasking license check failed');
    waitfor(dlg);
    return;
end

if nargin < 2
    flagUseBT2 = false;
else
    flagUseBT2 = reconOptions.flagUseBT2;
end

tic;

if ~exist('isT2IR','var')
    isT2IR = 1;
end
if ~exist('isT1rhoIR','var')
    isT1rhoIR = 1;
end

% Initialize fitting parameter
minT1 = 250e-3;
maxT1 = 5;
minT2 = 20e-3;
maxT2 = 200e-3;

% Load fitting parameters, overwrite initial values
extractVarFromStruct(fitParams);

% Scan parameters
Nseg = params.linesPerShot;
ns   = 1:Nseg;

moduleLength = params.moduleLength;
numIRns      = params.numIRns;
numT2prep    = params.numT2prep;
numT1rhoPrep = params.numT1rhoPrep;

TR = params.lEchoSpacing;

IRs    = [];
T2s    = [];
T1rhos = [];
TEs  = [];
TSLs = [];
if numIRns > 0
    IRs    = [IRs ones(1,numIRns)];
    T2s    = [T2s zeros(1,numIRns)];
    T1rhos = [T1rhos zeros(1,numIRns)];
    TEs  = [TEs zeros(1,numIRns)];
    TSLs = [TSLs zeros(1,numIRns)];
end
if numT2prep > 0
    IRs    = [IRs zeros(1,numT2prep)];
    T2s    = [T2s ones(1,numT2prep)]*isT2IR;
    T1rhos = [T1rhos zeros(1,numT2prep)];
    TEs  = [TEs  params.T2prepDuration];
    TSLs = [TSLs params.T2prepDuration*0];
end
if numT1rhoPrep > 0
    IRs    = [IRs zeros(1,numT1rhoPrep)];
    T2s    = [T2s ones(1,numT1rhoPrep)]*isT1rhoIR;
    T1rhos = [T1rhos zeros(1,numT1rhoPrep)];
    TEs  = [TEs  params.T1rhoDuration*0];
    TSLs = [TSLs params.T1rhoDuration];
end
Ncontrast = length(TEs);
rep_VE    = moduleLength/Ncontrast;
IRs       = repmat(IRs,    1, rep_VE);
T2s       = repmat(T2s,    1, rep_VE);
T1rhos    = repmat(T1rhos, 1, rep_VE);
TEs       = repmat(TEs,  1, rep_VE) *1e-3;
TSLs      = repmat(TSLs, 1, rep_VE) *1e-3;

rep_VFA    = moduleLength/numFA;
alphaArray = repmat(flipAngleArray, 1, rep_VFA)*pi/180;

invSign = -(T2s.*isT2IR + T1rhos.*isT1rhoIR);

% Update fitParams
fitParams.rep_VFA = rep_VFA;
fitParams.rep_VE  = rep_VE;
fitParams.Ncontrast = Ncontrast;
fitParams.IRs  = IRs;
fitParams.T2s  = T2s;
fitParams.T1rhos = T1rhos;
fitParams.TEs  = TEs;
fitParams.TSLs = TSLs;
fitParams.alphaArray = alphaArray;
fitParams.invSign = invSign;

row = @(x) x(:).';

%% Signal equation

e1  = @(T1) exp(-TR /T1);
e2  = @(T2) exp(-TEs/T2);
e1rho = @(T1rho)exp(-TSLs/T1rho);

Mss      = @(e1,alphas) (1-e1) ./ (1-cos(alphas)*e1);
step     = @(e1,alphas) (bsxfun(@power, e1*cos(alphas).', (ns-1))).';
sin_step = @(alphas)    sin(alphas);

Sint = @(A,e1,e2,e1rho,BalphaArray,BIR,BT2,Eff) A .* Mss(e1,BalphaArray) .* (1 + (step(e1,BalphaArray)) .* ((BIR*cos(pi)*IRs + (invSign.*BT2.*sin(pi/2)^2).*(T2s+T1rhos).*e2.*e1rho).*Eff-1)) .* sin_step(BalphaArray);

% if Nz > 1 && MBfactor == 1
    S = @(A,T1,T2,T1rho,Beta,BIR,BT2) row(Sint(  A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray, BIR, BT2, Eff_cellfunc_flow(T1,T2,T1rho,Beta*alphaArray,BIR,BT2,TR,TEs,TSLs,Nseg,invSign)));
% else
%     S = @(A,T1,T2,T1rho,Beta,BIR,BT2) row(Sint(  A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.9387, BIR, BT2, Eff_cellfunc_flow(T1,T2,T1rho,Beta*alphaArray*.9387,BIR,BT2,TR,TEs,TSLs,Nseg,invSign))...
%                                         + Sint(2*A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.5049, BIR, BT2, Eff_cellfunc_flow(T1,T2,T1rho,Beta*alphaArray*.5049,BIR,BT2,TR,TEs,TSLs,Nseg,invSign))...
%                                         + Sint(2*A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.0525, BIR, BT2, Eff_cellfunc_flow(T1,T2,T1rho,Beta*alphaArray*.0525,BIR,BT2,TR,TEs,TSLs,Nseg,invSign)));
% end


%% curves of different params

T1s = logspace(log10(minT1),log10(maxT1),21);
T1w = ones(size(T1s));

% T1s = 0.1:0.04:2.6;
% T1w = [488,731,1409,2513,4471,7823,5186,1125,843,684,670,744,678,794,838,896,1408,2038,1877,1507,1383,1247,1281,1352,1307,1399,1430,1532,1596,1549,1353,1292,1178,1160,1285,1283,1046,951,795,804,691,664,696,672,647,733,828,1258,1417,934,496,344,226,214,184,186,159,141,118,90,104,82,78,59,56,56,44,48,48,52,45,40,30,38,25,27,26,13,30,22,15,23,26,16,20,17,15,16,28,17,17,17,24,24,13,18,8,20];

% T1s = 0.5:0.04:2.6;
% T1w = [670,744,678,794,838,896,1408,2038,1877,1507,1383,1247,1281,1352,1307,1399,1430,1532,1596,1549,1353,1292,1178,1160,1285,1283,1046,951,795,804,691,664,696,672,647,733,828,1258,1417,934,496,344,226,214,184,186,159,141,118,90,104,82,78];

if numT2prep > 0 
    T2s = logspace(log10(minT2),log10(maxT2),21);
else
    T2s = 1;
end

if numT1rhoPrep > 0
    T1rhos = logspace(log10(minT2),log10(maxT2),21);
else
    T1rhos = 1;
end

Betas = linspace(0.2, 1.0, 9);
if strcmp(ScanType, 'SR')
    BIRs = linspace(0.3, 0.6, 4);
else
    BIRs = linspace(0.6, 1, 5);
end
BT2s  = linspace(0.6, 1, 5);
Bflows = 1;

fprintf('Bloch simulation: ');
strProgress = sprintf('0/%d ',numel(T1s)); fprintf(strProgress);
if flagUseBT2
    curves = zeros(Nseg*moduleLength, numel(T1s), numel(T2s),numel(T1rhos), numel(Betas), numel(BIRs), numel(BT2s),numel(Bflows));
    for j = 1:numel(T1s)
        prevLength = numel(strProgress);
        strProgress = sprintf('%d/%d ',j,numel(T1s));
        fprintf([repmat('\b',1,prevLength) '%s'],strProgress);
        for k = 1:numel(T2s)
            for l = 1:numel(T1rhos)
                for m = 1:numel(Betas)
                    for n = 1:numel(BIRs)
                        for o = 1:numel(BT2s)
                            for p = 1:numel(Bflows)
                                curves(:,j,k,l,m,n,o,p) = S(T1w(j),T1s(j),T2s(k),T1rhos(l),Betas(m),BIRs(n),BT2s(o));
                                %curves(:,j,k,l,n,m,p) = S(1,T1s(j),T2s(m),Balphas(k),BIRs(l),BT2s(n),Bflows(p));
                            end
                        end
                    end
                end
            end
        end  
    end
else
    curves = zeros(Nseg*moduleLength, numel(T1s), numel(T2s), numel(T1rhos), numel(Betas), numel(BIRs), numel(Bflows));
    for j = 1:numel(T1s)
        prevLength = numel(strProgress);
        strProgress = sprintf('%d/%d ',j,numel(T1s));
        fprintf([repmat('\b',1,prevLength) '%s'],strProgress);
        for k = 1:numel(T2s)
            for l = 1:numel(T1rhos)
                for m = 1:numel(Betas)
                    for n = 1:numel(BIRs)
                        for p = 1:numel(Bflows)
                            curves(:,j,k,l,m,n,p) = S(T1w(j),T1s(j),T2s(k),T1rhos(l),Betas(m),BIRs(n),BIRs(n));
                            %curves(:,j,k,l,m,p) = S(1,T1s(j),T2s(m),Balphas(k),BIRs(l),BIRs(l),Bflows(p));
                        end
                    end
                end
            end
        end
    end
end

fprintf('done. Do SVD... ');

[~,curveS,curveU] = svd(curves(:,:)','econ');
curveS = diag(curveS);

toc;


%% License check
function isValid = checkLicense(reconOptions)

% get machine ID
[~,strAdd] = genID;

md = java.security.MessageDigest.getInstance('MD5');    
try
    % load license file
    licenseID = loadLicenseFile(reconOptions);
    [hashes,dateNr] = getHash(licenseID);
    
    % check hash
    for n = 1:length(strAdd)
        ID   = double(strAdd{n})*sum(dateNr);
        hash = dec2hex(uint8(double(md.digest(ID))+128));
        hash = hash(:).';
        isValid = contains(hashes,hash);
        if isValid; break; end
    end
    if ~isValid
        fprintf(2,'License check error: invalid hash.\n');
        disp(['Current date: ' date]);
        disp('Current system info:');
        for n = 1:length(strAdd)
            disp(['    ' strAdd{n}]);
        end
    end
catch errormsg
    fprintf(2,'License check error: %s\n', errormsg.message);
    isValid = false;
end
