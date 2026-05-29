function dict = genT1Dictionary_VTR(params, fitParams, options)
% genT1Dictionary_VTR  T1 dictionary generation for VTR multitasking fitting.
%
% This function replaces the simple single-TR recovery model used in
% LRT_T1_recovery.m with the VTR Bloch model from gen_bloch_subspace_VTR.m.
% The returned dictionary is arranged as [time points x T1 candidates], which
% matches T1Mapping_DictFit_Func2_Anzhen.m.
%
% Example:
%   dict = genT1Dictionary_VTR(FitParams.params, FitParams);
%   Mz_dict_norm_abs_truc = dict.Mz_dict_norm_abs_truc;
%   T1s = dict.T1s;

if nargin < 2 || isempty(fitParams)
    fitParams = struct();
end
if nargin < 3 || isempty(options)
    options = struct();
end

row = @(x) x(:).';
vec = @(x) x(:);

params = setDefault(params, 'isVTR', true);
params = setDefault(params, 'isPrep2Seg', 0);
params = setDefault(params, 'moduleLength', 1);
params = setDefault(params, 'numSRns', 0);
params = setDefault(params, 'numIRns', 1);
params = setDefault(params, 'numT2prep', 0);
params = setDefault(params, 'numT1rhoPrep', 0);
params = setDefault(params, 'numDiffPrep', 0);
params = setDefault(params, 'numDiffPrepDirs', 0);
params = setDefault(params, 'isT2IR', 1);
params = setDefault(params, 'isT1rhoIR', 1);
params = setDefault(params, 'isDiffPrepIR', 1);

Nseg = getFirstField(params, {'linesPerShot', 'Nseg'}, []);
if isempty(Nseg)
    error('genT1Dictionary_VTR:MissingNseg', ...
        'params.linesPerShot or params.Nseg is required.');
end

Necho = getFirstField(params, {'Necho', 'NEco', 'NEco_old', 'NEcoOld'}, 1);
SGBlock = getFirstField(params, {'SGBlock', 'SGblock'}, 2);
TR = getFirstField(params, {'lEchoSpacing', 'TR'}, []);
if isempty(TR)
    error('genT1Dictionary_VTR:MissingTR', ...
        'params.lEchoSpacing or params.TR is required.');
end

flipAngleArray = getFirstField(params, {'adFlipAngleDegree', 'flipAngleArray', 'flipAngleDegree'}, 5);
flipAngleArray = row(flipAngleArray);
numFA = getOption(options, fitParams, 'numFA', numel(flipAngleArray));
if numFA == 1 && numel(flipAngleArray) > 1
    numFA = numel(flipAngleArray);
end

MBfactor = getOption(options, fitParams, 'MBfactor', 1);
Nz = getOption(options, fitParams, 'Nz', 1);
flagUseBT2 = getOption(options, fitParams, 'flagUseBT2', true);
ScanType = getOption(options, fitParams, 'ScanType', 'IR');
cutoff = getOption(options, fitParams, 'cutoff', 10);

minT1 = getOption(options, fitParams, 'minT1', 100e-3);
maxT1 = getOption(options, fitParams, 'maxT1', 2);
minT2 = getOption(options, fitParams, 'minT2', 20e-3);
maxT2 = getOption(options, fitParams, 'maxT2', 200e-3);
minT1rho = getOption(options, fitParams, 'minT1rho', 100e-3);
maxT1rho = getOption(options, fitParams, 'maxT1rho', 3);

T1s = getOption(options, fitParams, 'T1s', logspace(log10(minT1), log10(maxT1), 101));
T1w = getOption(options, fitParams, 'T1w', ones(size(T1s)));
Betas = getOption(options, fitParams, 'Betas', linspace(0.2, 1.0, 5));
if strcmpi(ScanType, 'SR')
    BIRs = getOption(options, fitParams, 'BIRs', linspace(0.3, 0.6, 4));
else
    BIRs = getOption(options, fitParams, 'BIRs', linspace(0.6, 1, 3));
end
BT2s = getOption(options, fitParams, 'BT2s', linspace(0.6, 1, 3));
Bflows = 1;

TR2 = TR;
if Necho >= 1 && params.isVTR && SGBlock == 2
    TRnav = getOption(options, fitParams, 'TRnav', []);
    if isempty(TRnav)
        alTE_seconds = getFirstField(params, {'alTE_seconds'}, []);
        if isempty(alTE_seconds)
            warning('genT1Dictionary_VTR:MissingTRnav', ...
                ['params.alTE_seconds was not found. Using TRnav = TR; ', ...
                'pass options.TRnav if the navigator TR is known.']);
            TRnav = TR;
        else
            TRnav = alTE_seconds(1)*2;
        end
    end
    TR2 = [TR TRnav];
    if mod(Nseg, 2) > 0
        warning('genT1Dictionary_VTR:OddNseg', ...
            'VTR with odd linesPerShot may leave an unmatched TR pair.');
    end
elseif Necho > 1 && params.isVTR && SGBlock > 2
    alTE_seconds = getFirstField(params, {'alTE_seconds'}, []);
    if ~isempty(alTE_seconds)
        TR = (TR*SGBlock - (alTE_seconds(end) - alTE_seconds(1))) / SGBlock;
    else
        warning('genT1Dictionary_VTR:MissingEchoTimes', ...
            'SGBlock > 2 VTR requested, but params.alTE_seconds is missing.');
    end
end

[SRflags, IRflags, T2flags, T1rhoFlags, TEs, TSLs, invSign, Preps, alphaArray] = ...
    buildPrepSchedule(params, flipAngleArray, numFA, row);

moduleLength = numel(TEs);
if params.isPrep2Seg
    SRflags = row([SRflags; zeros(size(SRflags))]);
    IRflags = row([IRflags; ones(size(IRflags))]);
    T2flags = row([T2flags; zeros(size(T2flags))]);
    T1rhoFlags = row([T1rhoFlags; zeros(size(T1rhoFlags))]);
    TEs = row([TEs; zeros(size(TEs))]);
    TSLs = row([TSLs; zeros(size(TSLs))]);
    alphaArray = row([alphaArray; alphaArray]);
    invSign = row([invSign; zeros(size(invSign))]);
end
moduleLength = numel(TEs);

if params.numT2prep > 0
    T2Candidates = logspace(log10(minT2), log10(maxT2), 21);
else
    T2Candidates = 1;
end
if params.numT1rhoPrep > 0
    T1rhoCandidates = logspace(log10(minT1rho), log10(maxT1rho), 21);
else
    T1rhoCandidates = 1;
end

S = makeSignalFunction(params, Necho, SGBlock, Nz, MBfactor, Nseg, TR, TR2, ...
    TEs, TSLs, Preps, alphaArray, SRflags, IRflags, T2flags, T1rhoFlags, ...
    invSign, row, vec);

fprintf('VTR T1 dictionary simulation: ');
strProgress = sprintf('0/%d ', numel(T1s));
fprintf(strProgress);

if flagUseBT2
    curves = zeros(Nseg*moduleLength, numel(T1s), numel(T2Candidates), ...
        numel(T1rhoCandidates), numel(Betas), numel(BIRs), numel(BT2s), numel(Bflows));
    for j = 1:numel(T1s)
        strProgress = printProgress(strProgress, j, numel(T1s));
        for k = 1:numel(T2Candidates)
            for l = 1:numel(T1rhoCandidates)
                for m = 1:numel(Betas)
                    for n = 1:numel(BIRs)
                        for o = 1:numel(BT2s)
                            curves(:, j, k, l, m, n, o, 1) = S( ...
                                T1w(j), T1s(j), T2Candidates(k), ...
                                T1rhoCandidates(l), Betas(m), BIRs(n), BT2s(o));
                        end
                    end
                end
            end
        end
    end
else
    curves = zeros(Nseg*moduleLength, numel(T1s), numel(T2Candidates), ...
        numel(T1rhoCandidates), numel(Betas), numel(BIRs), numel(Bflows));
    for j = 1:numel(T1s)
        strProgress = printProgress(strProgress, j, numel(T1s));
        for k = 1:numel(T2Candidates)
            for l = 1:numel(T1rhoCandidates)
                for m = 1:numel(Betas)
                    for n = 1:numel(BIRs)
                        curves(:, j, k, l, m, n, 1) = S( ...
                            T1w(j), T1s(j), T2Candidates(k), ...
                            T1rhoCandidates(l), Betas(m), BIRs(n), BIRs(n));
                    end
                end
            end
        end
    end
end
fprintf('done.\n');

betaIdx = findClosestIndex(Betas, getOption(options, fitParams, 'dictBeta', 1));
BIRIdx = findClosestIndex(BIRs, getOption(options, fitParams, 'dictBIR', max(BIRs)));
BT2Idx = findClosestIndex(BT2s, getOption(options, fitParams, 'dictBT2', max(BT2s)));

if flagUseBT2
    Mz_dict = squeeze(curves(:, :, 1, 1, betaIdx, BIRIdx, BT2Idx, 1));
else
    Mz_dict = squeeze(curves(:, :, 1, 1, betaIdx, BIRIdx, 1));
end

Mz_dict_norm_abs = abs(Mz_dict) ./ max(abs(Mz_dict), [], 1);
Mz_dict_norm_abs(~isfinite(Mz_dict_norm_abs)) = 0;
if cutoff > 0
    Mz_dict_norm_abs_truc = Mz_dict_norm_abs((cutoff + 1):end, :);
else
    Mz_dict_norm_abs_truc = Mz_dict_norm_abs;
end

dict = struct();
dict.curves = curves;
dict.T1s = T1s;
dict.T2s = T2Candidates;
dict.T1rhos = T1rhoCandidates;
dict.Betas = Betas;
dict.BIRs = BIRs;
dict.BT2s = BT2s;
dict.cutoff = cutoff;
dict.Mz_dict = Mz_dict;
dict.Mz_dict_norm_abs = Mz_dict_norm_abs;
dict.Mz_dict_norm_abs_truc = Mz_dict_norm_abs_truc;
dict.selectedIndex = struct('Beta', betaIdx, 'BIR', BIRIdx, 'BT2', BT2Idx);
dict.TR = TR;
dict.TR2 = TR2;
dict.Preps = Preps;
dict.alphaArray = alphaArray;

end

function S = makeSignalFunction(params, Necho, SGBlock, Nz, MBfactor, Nseg, TR, TR2, ...
    TEs, TSLs, Preps, alphaArray, SRflags, IRflags, T2flags, T1rhoFlags, ...
    invSign, row, vec)

if Necho >= 1 && params.isVTR && SGBlock == 2
    TRnav = TR2(2);
    e1 = @(T1) exp(-TR/T1);
    e1nav = @(T1) exp(-TRnav/T1);
    e2 = @(T2) exp(-TEs/T2);
    e1rho = @(T1rho) exp(-TSLs/T1rho);

    Mss1 = @(T1, alphas) ((1 - e1nav(T1)) + (1 - e1(T1))*e1nav(T1)*cos(alphas)) ./ ...
        (1 - e1(T1)*e1nav(T1)*cos(alphas).^2);
    Mss2 = @(T1, alphas) ((1 - e1(T1)) + (1 - e1nav(T1))*e1(T1)*cos(alphas)) ./ ...
        (1 - e1(T1)*e1nav(T1)*cos(alphas).^2);

    ns = 1:ceil(Nseg/2);
    cropNseg = @(x) x(1:Nseg, :);
    step = @(T1, alphas) bsxfun(@power, e1(T1)*e1nav(T1)*cos(alphas).^2, vec(ns - 1));
    sin_step = @(alphas) sin(alphas);

    prepTerm = @(e2v, e1rhov, BIR, BT2, Eff) ...
        (cos(invSign.*BIR*pi).*(IRflags + SRflags) + ...
        (invSign.*sin(BT2*pi/2).^2).*(T2flags + T1rhoFlags).*e2v.*e1rhov + ...
        cos(BT2*pi/2).^2.*(T2flags + T1rhoFlags)).*Eff - 1;

    Sint1 = @(A, T1, e2v, e1rhov, BalphaArray, BIR, BT2, Eff) ...
        A.*Mss1(T1, BalphaArray).*(1 + step(T1, BalphaArray).*prepTerm(e2v, e1rhov, BIR, BT2, Eff)).*sin_step(BalphaArray);
    Sint2 = @(A, T1, e2v, e1rhov, BalphaArray, BIR, BT2, Eff) ...
        (A.*Mss1(T1, BalphaArray).*e1(T1).*cos(BalphaArray).* ...
        (step(T1, BalphaArray).*prepTerm(e2v, e1rhov, BIR, BT2, Eff)) + ...
        A.*Mss2(T1, BalphaArray)).*sin_step(BalphaArray);
    Sint = @(A, T1, e2v, e1rhov, BalphaArray, BIR, BT2, Eff) ...
        row(cropNseg(reshape([row(Sint1(A, T1, e2v, e1rhov, BalphaArray, BIR, BT2, Eff)); ...
        row(Sint2(A, T1, e2v, e1rhov, BalphaArray, BIR, BT2, Eff))], ceil(Nseg/2)*2, [])));

    if Nz > 1 && MBfactor == 1
        S = @(A, T1, T2, T1rho, Beta, BIR, BT2) row(Sint(A, T1, e2(T2), e1rho(T1rho), ...
            Beta*alphaArray, BIR, BT2, Eff_cellfunc_VTR_local(T1, T2, T1rho, Beta*alphaArray, BIR, BT2, TR2, TEs, TSLs, Nseg, Preps)));
    elseif Nz > 1 && MBfactor > 1
        smsScales = [0.9683 0.7831 0.4996 0.2354 0.0751];
        S = @(A, T1, T2, T1rho, Beta, BIR, BT2) sumSMS(Sint, A, T1, T2, T1rho, Beta, BIR, BT2, ...
            smsScales, e2, e1rho, alphaArray, TR2, TEs, TSLs, Nseg, Preps, row);
    else
        smsScales = [0.9387 0.5049 0.0525];
        smsWeights = [1 2 2];
        S = @(A, T1, T2, T1rho, Beta, BIR, BT2) sumSMS(Sint, A, T1, T2, T1rho, Beta, BIR, BT2, ...
            smsScales, e2, e1rho, alphaArray, TR2, TEs, TSLs, Nseg, Preps, row, smsWeights);
    end
else
    error('genT1Dictionary_VTR:UnsupportedModel', ...
        'This generator is for the isVTR, multi-echo, SGBlock == 2 signal model.');
end

end

function y = sumSMS(Sint, A, T1, T2, T1rho, Beta, BIR, BT2, scales, e2, e1rho, ...
    alphaArray, TR2, TEs, TSLs, Nseg, Preps, row, weights)
if nargin < 19
    weights = ones(size(scales));
end
y = 0;
for ii = 1:numel(scales)
    BalphaArray = Beta*alphaArray*scales(ii);
    Eff = Eff_cellfunc_VTR_local(T1, T2, T1rho, BalphaArray, BIR, BT2, TR2, TEs, TSLs, Nseg, Preps);
    y = y + Sint(weights(ii)*A, T1, e2(T2), e1rho(T1rho), BalphaArray, BIR, BT2, Eff);
end
y = row(y);
end

function [SRflags, IRflags, T2flags, T1rhoFlags, TEs, TSLs, invSign, Preps, alphaArray] = ...
    buildPrepSchedule(params, flipAngleArray, numFA, row)

SRflags = [];
IRflags = [];
T2flags = [];
T1rhoFlags = [];
DiffFlags = [];
TEs = [];
TSLs = [];
invSign = [];

if params.numSRns > 0
    SRflags = [SRflags ones(1, params.numSRns)];
    IRflags = [IRflags zeros(1, params.numSRns)];
    T2flags = [T2flags zeros(1, params.numSRns)];
    T1rhoFlags = [T1rhoFlags zeros(1, params.numSRns)];
    DiffFlags = [DiffFlags zeros(1, params.numSRns)];
    TEs = [TEs zeros(1, params.numSRns)];
    TSLs = [TSLs zeros(1, params.numSRns)];
    invSign = [invSign 0.5*ones(1, params.numSRns)];
end
if params.numIRns > 0
    SRflags = [SRflags zeros(1, params.numIRns)];
    IRflags = [IRflags ones(1, params.numIRns)];
    T2flags = [T2flags zeros(1, params.numIRns)];
    T1rhoFlags = [T1rhoFlags zeros(1, params.numIRns)];
    DiffFlags = [DiffFlags zeros(1, params.numIRns)];
    TEs = [TEs zeros(1, params.numIRns)];
    TSLs = [TSLs zeros(1, params.numIRns)];
    invSign = [invSign -ones(1, params.numIRns)];
end
if params.numT2prep > 0
    T2prepDuration = getFirstField(params, {'T2prepDuration'}, 0);
    SRflags = [SRflags zeros(1, params.numT2prep)];
    IRflags = [IRflags zeros(1, params.numT2prep)];
    T2flags = [T2flags ones(1, params.numT2prep)];
    T1rhoFlags = [T1rhoFlags zeros(1, params.numT2prep)];
    DiffFlags = [DiffFlags zeros(1, params.numT2prep)];
    TEs = [TEs T2prepDuration];
    TSLs = [TSLs T2prepDuration*0];
    invSign = [invSign -ones(1, params.numT2prep)*(params.isT2IR - 0.5)*2];
end
if params.numT1rhoPrep > 0
    T1rhoDuration = getFirstField(params, {'T1rhoDuration'}, 0);
    SRflags = [SRflags zeros(1, params.numT1rhoPrep)];
    IRflags = [IRflags zeros(1, params.numT1rhoPrep)];
    T2flags = [T2flags zeros(1, params.numT1rhoPrep)];
    T1rhoFlags = [T1rhoFlags ones(1, params.numT1rhoPrep)];
    DiffFlags = [DiffFlags zeros(1, params.numT1rhoPrep)];
    TEs = [TEs T1rhoDuration*0];
    TSLs = [TSLs T1rhoDuration];
    invSign = [invSign -ones(1, params.numT1rhoPrep)*(params.isT1rhoIR - 0.5)*2];
end

Ncontrast = numel(TEs);
repVE = params.moduleLength/Ncontrast;
if abs(repVE - round(repVE)) > eps
    error('genT1Dictionary_VTR:BadModuleLength', ...
        'params.moduleLength must be an integer multiple of the prep contrast count.');
end
repVE = round(repVE);

SRflags = repmat(SRflags, 1, repVE);
IRflags = repmat(IRflags, 1, repVE);
T2flags = repmat(T2flags, 1, repVE);
T1rhoFlags = repmat(T1rhoFlags, 1, repVE);
DiffFlags = repmat(DiffFlags, 1, repVE);
TEs = repmat(TEs, 1, repVE)*1e-3;
TSLs = repmat(TSLs, 1, repVE)*1e-3;
invSign = repmat(invSign, 1, repVE);

repVFA = params.moduleLength/numFA;
if abs(repVFA - round(repVFA)) > eps
    error('genT1Dictionary_VTR:BadFlipAngleCount', ...
        'params.moduleLength must be an integer multiple of numFA.');
end
alphaArray = repmat(flipAngleArray, 1, round(repVFA))*pi/180;

Preps = (SRflags + IRflags + T2flags*2 + T1rhoFlags*3 + DiffFlags*2).*invSign;
if params.isPrep2Seg
    BlankPrepDurationMs = getFirstField(params, {'BlankPrepDurationMs'}, 0);
    Preps = [Preps; ones(size(Preps))*BlankPrepDurationMs*1e-3];
end
Preps = row(Preps);

end

function Eff_result = Eff_cellfunc_VTR_local(T1, T2, T1rho, alphaArray, BIR, BT2, TR, TEs, TSLs, Nseg, Preps)

invSign = Preps./abs(Preps);
invSign(Preps == 0.5) = 0.5;

moduleLength = numel(TEs);
Eff_result = zeros(1, moduleLength);
Eff_result(1) = solve_Eff_VTR_local(T1, T2, T1rho, alphaArray, BIR, BT2, TR, TEs, TSLs, Nseg, Preps);

if numel(TR) > 1
    TRnav = TR(2);
    TR = TR(1);
else
    TRnav = TR;
end

e1 = @(T1v) exp(-TR/T1v);
e1nav = @(T1v) exp(-TRnav/T1v);
e2 = @(T2v, shot) exp(-TEs(shot - 1)/T2v);
e1rho = @(T1rhov, shot) exp(-TSLs(shot - 1)/T1rhov);

Mss1 = @(e1v, e1navv, alpha, alphaPrev) ((1 - e1navv) + (1 - e1v)*e1navv*cos(alpha))* ...
    (1 - e1v*e1navv*cos(alphaPrev).^2) / ...
    ((1 - e1v*e1navv*cos(alpha).^2)*((1 - e1navv) + (1 - e1v)*e1navv*cos(alphaPrev)));
stepEnd = @(T1v, alpha) bsxfun(@power, e1(T1v)*e1nav(T1v)*cos(alpha).^2, ceil((Nseg + 1)/2));

for shot = 2:moduleLength
    FA1 = alphaArray(mod(shot - 2, moduleLength) + 1);
    FA2 = alphaArray(mod(shot - 3, moduleLength) + 1);

    if abs(Preps(shot - 1)) <= 1
        Eff_result(shot) = Mss1(e1(T1), e1nav(T1), FA1, FA2)* ...
            (1 + stepEnd(T1, FA1)*(cos(invSign(shot - 1)*BIR*pi)*Eff_result(shot - 1) - 1));
    else
        Eff_result(shot) = Mss1(e1(T1), e1nav(T1), FA1, FA2)* ...
            (1 + stepEnd(T1, FA1)*((cos(BT2*pi/2)^2 + ...
            invSign(shot - 1)*sin(BT2*pi/2)^2*e2(T2, shot)*e1rho(T1rho, shot))*Eff_result(shot - 1) - 1));
    end
end

end

function res = solve_Eff_VTR_local(T1, T2, T1rho, alphaArray, BIR, BT2, TR, TEs, TSLs, Nseg, Preps)

invSign = Preps./abs(Preps);
invSign(Preps == 0.5) = 0.5;

if numel(TR) > 1
    TRnav = TR(2);
    TR = TR(1);
else
    TRnav = TR;
end

e1 = @(T1v) exp(-TR/T1v);
e1nav = @(T1v) exp(-TRnav/T1v);
Mss1 = @(e1v, e1navv, alpha, alphaPrev) ((1 - e1navv) + (1 - e1v)*e1navv*cos(alpha))* ...
    (1 - e1v*e1navv*cos(alphaPrev).^2) / ...
    ((1 - e1v*e1navv*cos(alpha).^2)*((1 - e1navv) + (1 - e1v)*e1navv*cos(alphaPrev)));
stepEnd = @(T1v, alpha) bsxfun(@power, e1(T1v)*e1nav(T1v)*cos(alpha).^2, ceil((Nseg + 1)/2));

moduleLength = numel(TEs);
a = 1;
b = 0;
for shot = 1:moduleLength
    FA1 = alphaArray(mod(shot - 1, moduleLength) + 1);
    FA2 = alphaArray(mod(shot - 2, moduleLength) + 1);

    if abs(Preps(shot)) <= 1
        c = Mss1(e1(T1), e1nav(T1), FA1, FA2)*stepEnd(T1, FA1)*cos(invSign(shot)*BIR*pi);
    else
        c = Mss1(e1(T1), e1nav(T1), FA1, FA2)*stepEnd(T1, FA1)* ...
            (cos(BT2*pi/2)^2 + invSign(shot)*sin(BT2*pi/2)^2*exp(-TEs(shot)/T2)*exp(-TSLs(shot)/T1rho));
    end
    d = Mss1(e1(T1), e1nav(T1), FA1, FA2)*(1 + stepEnd(T1, FA1)*(-1));
    a = c*a;
    b = c*b + d;
end
res = -b/(a - 1);

end

function value = getOption(options, fitParams, fieldName, defaultValue)
if isstruct(options) && isfield(options, fieldName) && ~isempty(options.(fieldName))
    value = options.(fieldName);
elseif isstruct(fitParams) && isfield(fitParams, fieldName) && ~isempty(fitParams.(fieldName))
    value = fitParams.(fieldName);
elseif isstruct(fitParams) && isfield(fitParams, 'reconOptions') && ...
        isfield(fitParams.reconOptions, fieldName) && ~isempty(fitParams.reconOptions.(fieldName))
    value = fitParams.reconOptions.(fieldName);
else
    value = defaultValue;
end
end

function value = getFirstField(s, fields, defaultValue)
value = defaultValue;
for ii = 1:numel(fields)
    if isstruct(s) && isfield(s, fields{ii}) && ~isempty(s.(fields{ii}))
        value = s.(fields{ii});
        return;
    end
end
end

function s = setDefault(s, fieldName, defaultValue)
if ~isfield(s, fieldName) || isempty(s.(fieldName))
    s.(fieldName) = defaultValue;
end
end

function idx = findClosestIndex(values, target)
[~, idx] = min(abs(values - target));
end

function strProgress = printProgress(strProgress, current, total)
prevLength = numel(strProgress);
strProgress = sprintf('%d/%d ', current, total);
fprintf([repmat('\b', 1, prevLength) '%s'], strProgress);
end
