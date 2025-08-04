function [curveU, curveS, fitParams] = genCurveSubspace_T1_1FAn(params, fitParams, reconOptions)

% License check
% if ~checkLicense(reconOptions)
%     dlg = errordlg('Multitasking license check failed');
%     waitfor(dlg);
%     return;
% end

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
fitParams.reconOptions = reconOptions;
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

%% Signal equation

ns  = 1:Nseg;
e1  = @(T1) exp(-TR /T1);
Mss = @(e1,alpha)(1-e1) / (1-cos(alpha)*e1);

e1crusher = @(T1) exp(-0/T1);

Sint = @(A,e1,B,Balpha,e1crusher)A * Mss(e1,Balpha) * (1 - e1crusher*(B+1)*(e1*cos(Balpha)).^(ns-1)) * sin(Balpha);
S    = @(A,T1,B,Beta) Sint(A,e1(T1),B,Beta*alphaArray(1),e1crusher(T1));
% S = @(A,T1,alpha,B)Sint(A,e(T1),alpha,B)+Sint(A,e(T1),alpha/2,B); %slice prof (unnecessary! already linear combination of multiple flip angles)


%% curves of different params

T1s = logspace(log10(minT1),log10(maxT1),101);%logspace(log10(minT1),log10(maxT1),21);
T1w = ones(size(T1s));

% T1s = 0.1:0.04:4;
% T1w = [488,731,1409,2513,4471,7823,5186,1125,843,684,670,744,678,794,838,896,1408,2038,1877,1507,1383,1247,1281,1352,1307,1399,1430,1532,1596,1549,1353,1292,1178,1160,1285,1283,1046,951,795,804,691,664,696,672,647,733,828,1258,1417,934,496,344,226,214,184,186,159,141,118,90,104,82,78,59,56,56,44,48,48,52,45,40,30,38,25,27,26,13,30,22,15,23,26,16,20,17,15,16,28,17,17,17,24,24,13,18,8,20];
% 
% T1s = 0.5:0.04:2.6;
% T1w = [670,744,678,794,838,896,1408,2038,1877,1507,1383,1247,1281,1352,1307,1399,1430,1532,1596,1549,1353,1292,1178,1160,1285,1283,1046,951,795,804,691,664,696,672,647,733,828,1258,1417,934,496,344,226,214,184,186,159,141,118,90,104,82,78];
 
switch(params.ScanType)
    case {'IR','T2IR'}
        Bs  = linspace(0.5, 1.5, 21); % linspace(0.5, 1.5, 11);
        %alphas = (.5*pi/180:.5*pi/180:(alphaArray*1.5));
    case 'SR'
        Bs  = linspace(-0.25, 0.1, 7);
end
Betas = linspace(0.2, 1.5, 14);

fprintf('Bloch simulation: ');
strProgress = sprintf('0/%d ',numel(T1s)); fprintf(strProgress);
curves = zeros(Nseg, numel(T1s), numel(Bs), numel(Betas));
for j = 1:numel(T1s)
    prevLength = numel(strProgress);
    strProgress = sprintf('%d/%d ',j,numel(T1s));
    fprintf([repmat('\b',1,prevLength) '%s'],strProgress);
    for k = 1:numel(Bs)
        for l = 1:numel(Betas)
            curves(:,j,k,l) = S(T1w(j),T1s(j),Bs(k),Betas(l));
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