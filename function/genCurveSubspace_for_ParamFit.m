function [curves, T1s] = genCurveSubspace_for_ParamFit(params, fitParams, reconOptions)

% License check
% if ~checkLicense(reconOptions)
%     dlg = errordlg('Multitasking license check failed');
%     waitfor(dlg);
%     return;
% end

if nargin < 4
    flagUseBT2 = false;
else
    flagUseBT2 = reconOptions.flagUseBT2;
end

vec = @(x) x(:);
row = @(x) x(:).';

tic;

if ~exist('isT2IR','var')
    isT2IR = 1;
end
if ~exist('isT1rhoIR','var')
    isT1rhoIR = 1;
end

if ~isfield(params,'isPrep2Seg')
    params.isPrep2Seg = 0;
end

% Initialize fitting parameter
minT1 = 250e-3;
maxT1 = 5;
minT2 = 20e-3;
maxT2 = 200e-3;
minT1rho = 100e-3;
maxT1rho = 3;
minD = 0;
maxD = 2e-3;

% Load fitting parameters, overwrite initial values
fitParams.reconOptions = reconOptions;
extractVarFromStruct(fitParams);

minD = minD(1);
maxD = maxD(1);

% Scan parameters
Nseg = params.linesPerShot;
ns   = 1:Nseg;

moduleLength = params.moduleLength;
numSRns      = params.numSRns;
numIRns      = params.numIRns;
numT2prep    = params.numT2prep;
numT1rhoPrep = params.numT1rhoPrep;
numDiffPrep  = params.numDiffPrep;
numDiffPrepDirs = params.numDiffPrepDirs;

Necho = params.Necho;
isVTR = params.isVTR;
SGBlock = params.SGBlock;

TR = params.lEchoSpacing;
if Necho > 1 && isVTR && SGBlock == 2
    TRnav = params.lEchoSpacing - (params.alTE_seconds(end)-params.alTE_seconds(1));
    TR2 = [TR TRnav];
    if mod(Nseg,2) > 0
        fprintf(2,'Warning: VTR and odd linesPerShot.')
    end
elseif Necho > 1 && isVTR && SGBlock > 2
    fprintf(2,'Warning: VTR and SGBlock > 2.')
    TR = (params.lEchoSpacing*SGBlock - (params.alTE_seconds(end)-params.alTE_seconds(1)))/SGBlock;
end

SRs    = [];
IRs    = [];
T2s    = [];
T1rhos = [];
Diffs  = [];
TEs  = [];
TSLs = [];
invSign = [];
if numSRns > 0
    SRs    = [SRs ones(1,numSRns)];  
    IRs    = [IRs zeros(1,numSRns)];
    T2s    = [T2s zeros(1,numSRns)];
    T1rhos = [T1rhos zeros(1,numSRns)];
    Diffs  = [Diffs zeros(1,numSRns)];
    TEs    = [TEs zeros(1,numSRns)];
    TSLs   = [TSLs zeros(1,numSRns)];
    invSign = [invSign 0.5*ones(1,numSRns)];
end
if numIRns > 0
    SRs    = [SRs zeros(1,numIRns)];  
    IRs    = [IRs ones(1,numIRns)];
    T2s    = [T2s zeros(1,numIRns)];
    T1rhos = [T1rhos zeros(1,numIRns)];
    Diffs  = [Diffs zeros(1,numIRns)];
    TEs    = [TEs zeros(1,numIRns)];
    TSLs   = [TSLs zeros(1,numIRns)];
    invSign = [invSign -ones(1,numIRns)];
end
if numT2prep > 0
    SRs    = [SRs zeros(1,numT2prep)];  
    IRs    = [IRs zeros(1,numT2prep)];
    T2s    = [T2s ones(1,numT2prep)];
    T1rhos = [T1rhos zeros(1,numT2prep)];
    Diffs  = [Diffs zeros(1,numT2prep)];
    TEs    = [TEs  params.T2prepDuration];
    TSLs   = [TSLs params.T2prepDuration*0];
    invSign = [invSign -ones(1,numT2prep)*(params.isT2IR-0.5)*2];
end
if numT1rhoPrep > 0
    SRs    = [SRs zeros(1,numT1rhoPrep)];  
    IRs    = [IRs zeros(1,numT1rhoPrep)];
    T2s    = [T2s zeros(1,numT1rhoPrep)];
    T1rhos = [T1rhos ones(1,numT1rhoPrep)];
    Diffs  = [Diffs zeros(1,numT1rhoPrep)];
    TEs    = [TEs  params.T1rhoDuration*0];
    TSLs   = [TSLs params.T1rhoDuration];
    invSign = [invSign -ones(1,numT1rhoPrep)*(params.isT1rhoIR-0.5)*2];
end
if numDiffPrep > 0
    SRs    = [SRs zeros(1,numDiffPrep)];  
    IRs    = [IRs zeros(1,numDiffPrep)];
    T2s    = [T2s ones(1,numDiffPrep)];
    T1rhos = [T1rhos zeros(1,numDiffPrep)];
    Diffs  = [Diffs ones(1,numDiffPrep)];
    diffGradNorm = sqrt(sum(abs(reshape(params.DiffPrepGradTable,3,[])).^2,1));
    DiffGrads  = [zeros(3,numel(TEs)) reshape(params.DiffPrepGradTable,3,[])./diffGradNorm];
    DiffBvalue = [zeros(1,numel(TEs)) reshape(params.DiffPrepBvalueTable,1,[])];
    DiffGrads(~isfinite(DiffGrads)) = 0;
    TEs  = [TEs  ones(1,numDiffPrep)*params.DiffPrepTE];
    TSLs = [TSLs zeros(1,numDiffPrep)];
    invSign = [invSign -ones(1,numDiffPrep)*(params.isDiffPrepIR-0.5)*2];
end

Ncontrast = length(TEs);
rep_VE    = moduleLength/Ncontrast;
SRs       = repmat(SRs,    1, rep_VE);
IRs       = repmat(IRs,    1, rep_VE);
T2s       = repmat(T2s,    1, rep_VE);
T1rhos    = repmat(T1rhos, 1, rep_VE);
Diffs     = repmat(Diffs,  1, rep_VE);
TEs       = repmat(TEs,  1, rep_VE) *1e-3;
TSLs      = repmat(TSLs, 1, rep_VE) *1e-3;
invSign   = repmat(invSign, 1, rep_VE);
if numDiffPrep > 0
    DiffGrads  = repmat(DiffGrads, 1, rep_VE);
    DiffBvalue = repmat(DiffBvalue, 1, rep_VE);
end
rep_VFA    = moduleLength/numFA;
alphaArray = repmat(flipAngleArray, 1, rep_VFA)*pi/180;

Preps = (SRs + IRs*1 + T2s*2 + T1rhos*3 + Diffs*2).*invSign;
if params.isPrep2Seg
    Preps = [Preps;ones(size(Preps))*params.BlankPrepDurationMs*1e-3];
end

% Update fitParams
fitParams.rep_VFA = rep_VFA;
fitParams.rep_VE  = rep_VE;
fitParams.Ncontrast = Ncontrast;
fitParams.Preps = Preps;
fitParams.SRs  = SRs;
fitParams.IRs  = IRs;
fitParams.T2s  = T2s;
fitParams.T1rhos = T1rhos;
fitParams.invSign = invSign;
fitParams.TEs  = TEs;
fitParams.TSLs = TSLs;
fitParams.alphaArray = alphaArray;
if numDiffPrep > 0
    fitParams.Diffs = Diffs;
    fitParams.DiffGrads = DiffGrads;
    fitParams.DiffBvalue = DiffBvalue;
    fitParams.numDiffPrepDirs = numDiffPrepDirs;
end

if params.isPrep2Seg
    SRs    = row([SRs;zeros(size(SRs))]);
    IRs    = row([IRs;ones(size(IRs))]);
    T2s    = row([T2s;zeros(size(T2s))]);
    T1rhos = row([T1rhos;zeros(size(T1rhos))]);
    Diffs  = row([Diffs;zeros(size(Diffs))]);
    TEs    = row([TEs;zeros(size(TEs))]);
    TSLs   = row([TSLs;zeros(size(TSLs))]);
    alphaArray = row([alphaArray,alphaArray]);
    invSign    = row([invSign;zeros(size(invSign))]);
end

%% Signal equation
 
if Necho > 1 && isVTR && SGBlock == 2
    e1    = @(T1) exp(-TR /T1);
    e1nav = @(T1) exp(-TRnav/T1);
    e2    = @(T2) exp(-TEs/T2);
    e1rho = @(T1rho)exp(-TSLs/T1rho);
    
    Mss1  = @(T1,alphas) ((1-e1nav(T1))+(1-e1(T1))*e1nav(T1)*cos(alphas))./ (1-e1(T1)*e1nav(T1)*cos(alphas).^2);
    Mss2  = @(T1,alphas) ((1-e1(T1))+(1-e1nav(T1))*e1(T1)*cos(alphas))./ (1-e1(T1)*e1nav(T1)*cos(alphas).^2);

    cropNseg = @(x) x(1:Nseg,:);
    ns    = 1:ceil(Nseg/2);
    step  = @(T1,alphas) (bsxfun(@power, e1(T1)*e1nav(T1)*cos(alphas).^2, vec(ns-1)));
    
    sin_step = @(alphas) sin(alphas);
    
    Sint1 = @(A,T1,e2,e1rho,BalphaArray,BIR,BT2,Eff)      A .* Mss1(T1,BalphaArray) .* (1 + step(T1,BalphaArray) .* ((cos(invSign.*BIR*pi).*(IRs+SRs) + (invSign.*sin(BT2*pi/2).^2).*(T2s+T1rhos).*e2.*e1rho + cos(BT2*pi/2).^2.*(T2s+T1rhos)).*Eff-1)) .* sin_step(BalphaArray);
    Sint2 = @(A,T1,e2,e1rho,BalphaArray,BIR,BT2,Eff)  (   A .* Mss1(T1,BalphaArray) .*e1(T1) .* cos(BalphaArray) .* (step(T1,BalphaArray) .* ((cos(invSign.*BIR*pi).*(IRs+SRs) + (invSign.*sin(BT2*pi/2).^2).*(T2s+T1rhos).*e2.*e1rho + cos(BT2*pi/2).^2.*(T2s+T1rhos)).*Eff-1)) ...
                                                        + A .* Mss2(T1,BalphaArray) ) .* sin_step(BalphaArray);
                                                                        
    Sint = @(A,T1,e2,e1rho,BalphaArray,BIR,BT2,Eff)   row(cropNseg(reshape([row(Sint1(A,T1,e2,e1rho,BalphaArray,BIR,BT2,Eff));row(Sint2(A,T1,e2,e1rho,BalphaArray,BIR,BT2,Eff))],ceil(Nseg/2)*2,[])));
    
    if Nz > 1 && MBfactor == 1
        S = @(A,T1,T2,T1rho,Beta,BIR,BT2) row( Sint(  A, T1, e2(T2), e1rho(T1rho), Beta*alphaArray, BIR, BT2, Eff_cellfunc_VTR(T1,T2,T1rho,Beta*alphaArray,BIR,BT2,TR2,TEs,TSLs,Nseg,Preps)));
    elseif Nz > 1 && MBfactor > 1   % 2D-SMS
        S = @(A,T1,T2,T1rho,Beta,BIR,BT2) row( Sint(  A, T1, e2(T2), e1rho(T1rho), Beta*alphaArray*.9683, BIR, BT2, Eff_cellfunc_VTR(T1,T2,T1rho,Beta*alphaArray*.9683,BIR,BT2,TR2,TEs,TSLs,Nseg,Preps))...
                                             + Sint(  A, T1, e2(T2), e1rho(T1rho), Beta*alphaArray*.7831, BIR, BT2, Eff_cellfunc_VTR(T1,T2,T1rho,Beta*alphaArray*.7831,BIR,BT2,TR2,TEs,TSLs,Nseg,Preps))...
                                             + Sint(  A, T1, e2(T2), e1rho(T1rho), Beta*alphaArray*.4996, BIR, BT2, Eff_cellfunc_VTR(T1,T2,T1rho,Beta*alphaArray*.4996,BIR,BT2,TR2,TEs,TSLs,Nseg,Preps))...
                                             + Sint(  A, T1, e2(T2), e1rho(T1rho), Beta*alphaArray*.2354, BIR, BT2, Eff_cellfunc_VTR(T1,T2,T1rho,Beta*alphaArray*.2354,BIR,BT2,TR2,TEs,TSLs,Nseg,Preps))...
                                             + Sint(  A, T1, e2(T2), e1rho(T1rho), Beta*alphaArray*.0751, BIR, BT2, Eff_cellfunc_VTR(T1,T2,T1rho,Beta*alphaArray*.0751,BIR,BT2,TR2,TEs,TSLs,Nseg,Preps))); 
    else
        S = @(A,T1,T2,T1rho,Beta,BIR,BT2) row( Sint(  A, T1, e2(T2), e1rho(T1rho), Beta*alphaArray*.9387, BIR, BT2, Eff_cellfunc_VTR(T1,T2,T1rho,Beta*alphaArray*.9387,BIR,BT2,TR2,TEs,TSLs,Nseg,Preps))...
                                             + Sint(2*A, T1, e2(T2), e1rho(T1rho), Beta*alphaArray*.5049, BIR, BT2, Eff_cellfunc_VTR(T1,T2,T1rho,Beta*alphaArray*.5049,BIR,BT2,TR2,TEs,TSLs,Nseg,Preps))...
                                             + Sint(2*A, T1, e2(T2), e1rho(T1rho), Beta*alphaArray*.0525, BIR, BT2, Eff_cellfunc_VTR(T1,T2,T1rho,Beta*alphaArray*.0525,BIR,BT2,TR2,TEs,TSLs,Nseg,Preps)));
    end
else
    e1  = @(T1) exp(-TR /T1);
    e2  = @(T2) exp(-TEs/T2);
    e1rho = @(T1rho)exp(-TSLs/T1rho);
    
    Mss      = @(e1,alphas) (1-e1) ./ (1-cos(alphas)*e1);
    step     = @(e1,alphas) (bsxfun(@power, e1*cos(alphas).', (ns-1))).';
    sin_step = @(alphas)    sin(alphas);

    if params.isPrep2Seg
        Sint = @(A,e1,e2,e1rho,BalphaArray,BIR,BT2,Eff) A .* Mss(e1,BalphaArray) .* (1 + step(e1,BalphaArray) .* ((cos(invSign.*BIR*pi).*(IRs+SRs) + (invSign.*sin(BT2*pi/2).^2).*(T2s+T1rhos).*e2.*e1rho + cos(BT2*pi/2).^2.*(T2s+T1rhos)).*Eff-1)) .* sin_step(BalphaArray);
        
        if Nz > 1 && MBfactor == 1
            S = @(A,T1,T2,T1rho,Beta,BIR,BT2) row( Sint(  A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray, BIR, BT2, Eff_cellfunc_2seg(T1,T2,T1rho,Beta*alphaArray,BIR,BT2,TR,TEs,TSLs,Nseg,Preps)));
        elseif Nz > 1 && MBfactor > 1   % 2D-SMS
            S = @(A,T1,T2,T1rho,Beta,BIR,BT2) row( Sint(  A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.9683, BIR, BT2, Eff_cellfunc_2seg(T1,T2,T1rho,Beta*alphaArray*.9683,BIR,BT2,TR,TEs,TSLs,Nseg,Preps))...
                                                 + Sint(  A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.7831, BIR, BT2, Eff_cellfunc_2seg(T1,T2,T1rho,Beta*alphaArray*.7831,BIR,BT2,TR,TEs,TSLs,Nseg,Preps))...
                                                 + Sint(  A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.4996, BIR, BT2, Eff_cellfunc_2seg(T1,T2,T1rho,Beta*alphaArray*.4996,BIR,BT2,TR,TEs,TSLs,Nseg,Preps))...
                                                 + Sint(  A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.2354, BIR, BT2, Eff_cellfunc_2seg(T1,T2,T1rho,Beta*alphaArray*.2354,BIR,BT2,TR,TEs,TSLs,Nseg,Preps))...
                                                 + Sint(  A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.0751, BIR, BT2, Eff_cellfunc_2seg(T1,T2,T1rho,Beta*alphaArray*.0751,BIR,BT2,TR,TEs,TSLs,Nseg,Preps))); 
        else
            S = @(A,T1,T2,T1rho,Beta,BIR,BT2) row( Sint(  A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.9387, BIR, BT2, Eff_cellfunc_2seg(T1,T2,T1rho,Beta*alphaArray*.9387,BIR,BT2,TR,TEs,TSLs,Nseg,Preps))...
                                                 + Sint(2*A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.5049, BIR, BT2, Eff_cellfunc_2seg(T1,T2,T1rho,Beta*alphaArray*.5049,BIR,BT2,TR,TEs,TSLs,Nseg,Preps))...
                                                 + Sint(2*A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.0525, BIR, BT2, Eff_cellfunc_2seg(T1,T2,T1rho,Beta*alphaArray*.0525,BIR,BT2,TR,TEs,TSLs,Nseg,Preps)));
        end
    else
        Sint = @(A,e1,e2,e1rho,BalphaArray,BIR,BT2,Eff) A .* Mss(e1,BalphaArray) .* (1 + step(e1,BalphaArray) .* ((cos(invSign.*BIR*pi).*(IRs+SRs) + (invSign.*sin(BT2*pi/2).^2).*(T2s+T1rhos).*e2.*e1rho + cos(BT2*pi/2).^2.*(T2s+T1rhos)).*Eff-1)) .* sin_step(BalphaArray);
        
        %Sint = @(A,e1,e2,e1rho,BalphaArray,BIR,BT2,Eff) A .* Mss(e1,BalphaArray) .* (1 + step(e1,BalphaArray) .* ((repmat([cos(BIR*pi)*ones(1,numIRns) (((-1)^isT2IR)*sin(BT2*pi/2)^2)*ones(1,numT2prep) (((-1)^isT1rhoIR)*sin(BT2*pi/2)^2)*ones(1,numT1rhoPrep)],1,rep_VE).*e2.*e1rho + repmat([zeros(1,numIRns) cos(BT2*pi/2)^2*ones(1,numT2prep+numT1rhoPrep)],1,rep_VE)).*Eff-1)) .* sin_step(BalphaArray);
        
        if Nz > 1 && MBfactor == 1
            S = @(A,T1,T2,T1rho,Beta,BIR,BT2) row( Sint(  A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray, BIR, BT2, Eff_cellfunc(T1,T2,T1rho,Beta*alphaArray,BIR,BT2,TR,TEs,TSLs,Nseg,Preps)));
        elseif Nz > 1 && MBfactor > 1   % 2D-SMS
            S = @(A,T1,T2,T1rho,Beta,BIR,BT2) row( Sint(  A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.9683, BIR, BT2, Eff_cellfunc(T1,T2,T1rho,Beta*alphaArray*.9683,BIR,BT2,TR,TEs,TSLs,Nseg,Preps))...
                                                 + Sint(  A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.7831, BIR, BT2, Eff_cellfunc(T1,T2,T1rho,Beta*alphaArray*.7831,BIR,BT2,TR,TEs,TSLs,Nseg,Preps))...
                                                 + Sint(  A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.4996, BIR, BT2, Eff_cellfunc(T1,T2,T1rho,Beta*alphaArray*.4996,BIR,BT2,TR,TEs,TSLs,Nseg,Preps))...
                                                 + Sint(  A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.2354, BIR, BT2, Eff_cellfunc(T1,T2,T1rho,Beta*alphaArray*.2354,BIR,BT2,TR,TEs,TSLs,Nseg,Preps))...
                                                 + Sint(  A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.0751, BIR, BT2, Eff_cellfunc(T1,T2,T1rho,Beta*alphaArray*.0751,BIR,BT2,TR,TEs,TSLs,Nseg,Preps))); 
        else
            S = @(A,T1,T2,T1rho,Beta,BIR,BT2) row( Sint(  A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.9387, BIR, BT2, Eff_cellfunc(T1,T2,T1rho,Beta*alphaArray*.9387,BIR,BT2,TR,TEs,TSLs,Nseg,Preps))...
                                                 + Sint(2*A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.5049, BIR, BT2, Eff_cellfunc(T1,T2,T1rho,Beta*alphaArray*.5049,BIR,BT2,TR,TEs,TSLs,Nseg,Preps))...
                                                 + Sint(2*A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.0525, BIR, BT2, Eff_cellfunc(T1,T2,T1rho,Beta*alphaArray*.0525,BIR,BT2,TR,TEs,TSLs,Nseg,Preps)));
        end
    end
end
%% curves of different params

T1s = logspace(log10(minT1),log10(maxT1),101);
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
    T1rhos = logspace(log10(minT1rho),log10(maxT1rho),21);
else
    T1rhos = 1;
end
if numDiffPrep > 0
    DiffCoeff = logspace(log10(0.0001),log10(maxD),5);
    T2s = logspace(log10(minT2),log10(maxT2),11);
end

Betas = linspace(0.2, 1.0, 5);
if strcmp(ScanType, 'SR')
    BIRs = linspace(0.3, 0.6, 4);
else
    BIRs = linspace(0.6, 1, 3);
end
BT2s  = linspace(0.6, 1, 3);
Bflows = 1;

fprintf('Bloch simulation: ');
strProgress = sprintf('0/%d ',numel(T1s)); fprintf(strProgress);
if flagUseBT2
    curves = zeros(Nseg*numel(IRs), numel(T1s), numel(T2s),numel(T1rhos), numel(Betas), numel(BIRs), numel(BT2s),numel(Bflows));
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
    if numDiffPrep > 0
        curves = zeros(Nseg*moduleLength, numel(T1s), numel(T2s), numel(DiffCoeff), numel(Betas), numel(BIRs), numel(Bflows));
    else
        curves = zeros(Nseg*moduleLength, numel(T1s), numel(T2s), numel(T1rhos), numel(Betas), numel(BIRs), numel(Bflows));
    end
    for j = 1:numel(T1s)
        prevLength = numel(strProgress);
        strProgress = sprintf('%d/%d ',j,numel(T1s));
        fprintf([repmat('\b',1,prevLength) '%s'],strProgress);
        for k = 1:numel(T2s)
            for l = 1:numel(T1rhos)
                for m = 1:numel(Betas)
                    for n = 1:numel(BIRs)
                        for p = 1:numel(Bflows)
                            if numDiffPrep > 0
                                for q = 1:numel(DiffCoeff)^(numDiffPrepDirs)
                                    DArray = zeros(1,numDiffPrepDirs);
                                    for dir = 1:numDiffPrepDirs
                                        DArray(dir) = DiffCoeff(ceil((mod((q-1),numel(DiffCoeff)^dir) + 1)/numel(DiffCoeff)^(dir-1)));
                                    end
                                    curves(:,j,k,q,m,n,p) = S(T1w(j),T1s(j),T2s(k),BIRs(n),BIRs(n),Betas(m),DArray);
                                end
                            else
                                curves(:,j,k,l,m,n,p) = S(T1w(j),T1s(j),T2s(k),T1rhos(l),Betas(m),BIRs(n),BIRs(n));
                                %curves(:,j,k,l,n,m,p) = S(1,T1s(j),T2s(m),Balphas(k),BIRs(l),BT2s(n),Bflows(p));
                            end
                        end
                    end
                end
            end
        end
    end
end

% fprintf('done. Do SVD... ');
% 
% [~,curveS,curveU] = svd(curves(:,:)','econ');
% curveS = diag(curveS);

toc;

%% Build D vector for 3-scan trace
function vecD = genDvec(DiffBvalue,DArray)

vecD = zeros(1,numel(DiffBvalue));
[b,~,idxc] = unique(DiffBvalue);
for n = 1:numel(b)
    if b(n) > 0
        repD = numel(find(idxc==n))/numel(DArray);
        vecD(idxc==n) = repmat(DArray,[1,repD]);
    end
end


% %% License check
% function isValid = checkLicense(reconOptions)
% 
% % get machine ID
% [~,strAdd] = genID;
% 
% md = java.security.MessageDigest.getInstance('MD5');    
% try
%     % load license file
%     licenseID = loadLicenseFile(reconOptions);
%     [hashes,dateNr] = getHash(licenseID);
% 
%     % check hash
%     for n = 1:length(strAdd)
%         ID   = double(strAdd{n})*sum(dateNr);
%         hash = dec2hex(uint8(double(md.digest(ID))+128));
%         hash = hash(:).';
%         isValid = contains(hashes,hash);
%         if isValid; break; end
%     end
%     if ~isValid
%         fprintf(2,'License check error: invalid hash.\n');
%         disp(['Current date: ' date]);
%         disp('Current system info:');
%         for n = 1:length(strAdd)
%             disp(['    ' strAdd{n}]);
%         end
%     end
% catch errormsg
%     fprintf(2,'License check error: %s\n', errormsg.message);
%     isValid = false;
% end
