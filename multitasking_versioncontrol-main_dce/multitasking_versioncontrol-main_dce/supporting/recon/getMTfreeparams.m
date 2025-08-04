function paramsMTfree = getMTfreeparams(twixObj, params)

if nargin > 1
    paramsMTfree = params;
end

row = @(x) x(:)';

paramsMTfree.SoftwareVersions = twixObj.hdr.Dicom.SoftwareVersions;
temp = regexp(twixObj.hdr.Dicom.SoftwareVersions,' ','split');
VerString = temp{end};

% initialize paramsMTfree
paramsMTfree.SGBlock  = 2;
paramsMTfree.navShift = 1;
paramsMTfree.linesPerShot = twixObj.hdr.Meas.lSegments;

paramsMTfree.MBfactor = 1;
paramsMTfree.Ntrajs   = 0;
paramsMTfree.trajInc  = 0;
try
    paramsMTfree.linOrder_raw = twixObj.image.iceParam(11,:);
catch
    paramsMTfree.linOrder_raw = twixObj.image.iceParam(1,:)*0;
end
paramsMTfree.asymm_echo = 0.5;

paramsMTfree.moduleLength = 1;
paramsMTfree.numFA        = 1;
paramsMTfree.flipAngleArray = twixObj.hdr.MeasYaps.adFlipAngleDegree{1};
paramsMTfree.numIRns      = 0;
paramsMTfree.numSRns      = 0;
paramsMTfree.numIRsel     = 0;
paramsMTfree.numSRsel     = 0;
paramsMTfree.numT2prep    = 0;
paramsMTfree.T2prepDuration = [];
paramsMTfree.isT2IR       = 0;
paramsMTfree.numT1rhoPrep = 0;
paramsMTfree.T1rhoDuration = [];
paramsMTfree.isT1rhoIR     = 0;

paramsMTfree.isTinyGoldenAngle = 0;
paramsMTfree.isNoNav           = 0;
paramsMTfree.isVTR             = 0;

try
    % parse tFree as a json string
    tFree = twixObj.hdr.Meas.tFree;
    tFree(end+1) = '}';
    temp = jsondecode(tFree);
    vars = fieldnames(temp);
    for n = 1:numel(vars)
        paramsMTfree.(vars{n}) = temp.(vars{n});
    end

    if isfield(paramsMTfree,'numIR')
        paramsMTfree.numIRns = paramsMTfree.numIR;
        paramsMTfree = rmfield(paramsMTfree,'numIR');
    end
    if isfield(paramsMTfree,'numSR')
        paramsMTfree.numSRns = paramsMTfree.numSR;
        paramsMTfree = rmfield(paramsMTfree,'numSR');
    end
    if ~isfield(paramsMTfree,'numIRsel')
        paramsMTfree.numIRsel = 0;
    end
    if isfield(paramsMTfree,'CESTSatFreqShiftppmList')
        paramsMTfree.CESTSatFreqOffsetppmList = paramsMTfree.CESTSatFreqShiftppmList;
        paramsMTfree = rmfield(paramsMTfree,'CESTSatFreqShiftppmList');
    end
    for n = 1:paramsMTfree.numFA
        paramsMTfree.flipAngleArray(n) = twixObj.hdr.MeasYaps.adFlipAngleDegree{n};
    end

    paramsMTfree.T2prepDuration = row(paramsMTfree.T2prepDuration(1:paramsMTfree.numT2prep));
    paramsMTfree.T1rhoDuration  = row(paramsMTfree.T1rhoDuration(1:paramsMTfree.numT1rhoPrep));  
    paramsMTfree.moduleLength   = paramsMTfree.numSRns + paramsMTfree.numIRns + paramsMTfree.numT2prep + paramsMTfree.numT1rhoPrep;

    if paramsMTfree.numFA > 0 && paramsMTfree.moduleLength > 0
        paramsMTfree.moduleLength = lcm(paramsMTfree.numFA, paramsMTfree.moduleLength);
    elseif paramsMTfree.numFA > 0
        paramsMTfree.moduleLength = paramsMTfree.numFA;
    end

    if isfield(paramsMTfree,'navShift')
        paramsMTfree.navShift = mod(paramsMTfree.navShift,paramsMTfree.SGBlock) + 1;
    end
    
catch
    % tFree is not a json string, extract MT parameters from sWipMemBlock
    try
        numIRns_idx   = 12;
        numSR_idx     = numIRns_idx + 1;
        numIRsel_idx  = 39;
        numT2prep_idx = 18;
        numT1rho_idx  = 48;

        if strcmp(twixObj.image.softwareVersion,'vb')
            sWipMemBlock = twixObj.hdr.MeasYaps.sWiPMemBlock;
        else
            sWipMemBlock = twixObj.hdr.MeasYaps.sWipMemBlock;
        end
        Nalfree = numel(sWipMemBlock.alFree);

        % SG block length (# of nav lines + # of imaging lines)
        if Nalfree >= 26 && ~isempty(sWipMemBlock.alFree{26})
            paramsMTfree.SGBlock  = sWipMemBlock.alFree{26};
        else
            paramsMTfree.SGBlock  = sWipMemBlock.alFree{3} - 2 + 1;     % SGBlock + ImaLinesPerSGBlock - (abs(minVal)+1) + 1
        end

        % Lines per shot / segments
        if Nalfree >= 32 && ~isempty(sWipMemBlock.alFree{32})
            paramsMTfree.linesPerShot = sWipMemBlock.alFree{32};
        else
            paramsMTfree.linesPerShot = sWipMemBlock.alFree{7} - 9;     % abs(minVal)+1 = 9
        end

        % Number of flip angles and flip angle array
        if Nalfree >= 10 && ~isempty(sWipMemBlock.alFree{10})
            paramsMTfree.numFA = sWipMemBlock.alFree{10} - 2;
        end
        paramsMTfree.numFA = min(paramsMTfree.numFA,numel(twixObj.hdr.MeasYaps.adFlipAngleDegree));
        for n = 1:paramsMTfree.numFA
            paramsMTfree.flipAngleArray(n) = twixObj.hdr.MeasYaps.adFlipAngleDegree{n};
        end
    catch
        paramsMTfree.SGBlock  = 2;
        paramsMTfree.navShift = 1;
        paramsMTfree.linesPerShot = twixObj.hdr.Meas.lSegments;
        paramsMTfree.numFA = 1;
        paramsMTfree.flipAngleArray = twixObj.hdr.MeasYaps.adFlipAngleDegree{1};
    end
    try
        % Gap time for MATCH protocol
        if Nalfree >= 8 && ~isempty(sWipMemBlock.alFree{8})
            paramsMTfree.MATCH_GapTime = sWipMemBlock.alFree{8} - 1;
        else
            paramsMTfree.MATCH_GapTime = 0;
        end

        % MB factor (only relevant in BEAT_MT)
        if Nalfree >= 36 && ~isempty(sWipMemBlock.alFree{36})
            paramsMTfree.MBfactor = sWipMemBlock.alFree{36};
        else
            paramsMTfree.MBfactor = 1;
        end

        % Radial trajectory 
        if Nalfree >= 34 && ~isempty(sWipMemBlock.alFree{34})
            paramsMTfree.Ntrajs  = sWipMemBlock.alFree{33};
            paramsMTfree.trajInc = sWipMemBlock.alFree{34};
        else
            paramsMTfree.Ntrajs  = 0;
            paramsMTfree.trajInc = 0;
        end
    catch
        paramsMTfree.MBfactor = 1;
        paramsMTfree.Ntrajs  = 0;
        paramsMTfree.trajInc = 0;
    end
    try
        paramsMTfree.isTinyGoldenAngle = (sWipMemBlock.alFree{62} > 1);
        paramsMTfree.isNoNav           = (sWipMemBlock.alFree{62} == 2);
        paramsMTfree.isVTR             = sWipMemBlock.alFree{63} - 1;
    catch
        paramsMTfree.isTinyGoldenAngle = 0;
        paramsMTfree.isNoNav           = 0;
        paramsMTfree.isVTR             = 0;
    end
    try
        paramsMTfree.moduleLength = 0;
        % Number of non-selective inversion modules
        if Nalfree >= numIRns_idx && ~isempty(sWipMemBlock.alFree{numIRns_idx})
            paramsMTfree.numIRns = sWipMemBlock.alFree{numIRns_idx} - 1;
            paramsMTfree.moduleLength = paramsMTfree.moduleLength + paramsMTfree.numIRns;
        else
            fprintf('Reading alFree: IRns module not specified.\n');
            paramsMTfree.numIRns = 0;
        end
        % Number of slice-selective inversion modules
        if Nalfree >= numIRsel_idx && ~isempty(sWipMemBlock.alFree{numIRsel_idx})
            paramsMTfree.numIRsel = sWipMemBlock.alFree{numIRsel_idx} - 1;
            paramsMTfree.moduleLength = paramsMTfree.moduleLength + paramsMTfree.numIRsel;
        else
            fprintf('Reading alFree: IRsel module not specified.\n');
            paramsMTfree.numIRsel = 0;
        end
        % Number of saturation modules
        if Nalfree >= numSR_idx && ~isempty(sWipMemBlock.alFree{numSR_idx})
            paramsMTfree.numSRns = sWipMemBlock.alFree{numSR_idx} - 1;
            paramsMTfree.moduleLength = paramsMTfree.moduleLength + paramsMTfree.numSRns;
        else
            fprintf('Reading alFree: SR module not specified.\n');
            paramsMTfree.numSRns = 0;
        end
        % Number of T2/T2IR preps and T2 prep duration array
        if Nalfree >= numT2prep_idx + 6 && ~isempty(sWipMemBlock.alFree{numT2prep_idx})
            minT2prep = 16;
            if Nalfree >= 64 && strcmp(VerString,'E11')
                if sWipMemBlock.alFree{64} > 2208
                    minT2prep = 1;
                end     
            end
            paramsMTfree.numT2prep = sWipMemBlock.alFree{numT2prep_idx} - 1;
            paramsMTfree.T2prepDuration = [];
            for n = 1:paramsMTfree.numT2prep
                paramsMTfree.T2prepDuration(n) = sWipMemBlock.alFree{numT2prep_idx+n} - minT2prep;
            end
            paramsMTfree.isT2IR = sign(sWipMemBlock.alFree{numT2prep_idx-2} - 1);
            paramsMTfree.moduleLength = paramsMTfree.moduleLength + paramsMTfree.numT2prep;
        else
            fprintf('Reading alFree: T2-IR prep module not specified.\n');
            paramsMTfree.numT2prep = 0;
            paramsMTfree.T2prepDuration = [];
            paramsMTfree.isT2IR = 0;
        end
        % Number of T1rho preps and T1rho prep duration array
        if Nalfree >= numT1rho_idx + 6 && ~isempty(sWipMemBlock.alFree{numT1rho_idx})
            paramsMTfree.numT1rhoPrep = sWipMemBlock.alFree{numT1rho_idx} - 1;
            paramsMTfree.T1rhoDuration = [];
            for n = 1:paramsMTfree.numT1rhoPrep
                paramsMTfree.T1rhoDuration(n) = sWipMemBlock.alFree{numT1rho_idx+n} - 7;
            end
            paramsMTfree.isT1rhoIR = sign(sWipMemBlock.alFree{numT1rho_idx-5} - 1);
            paramsMTfree.moduleLength = paramsMTfree.moduleLength + paramsMTfree.numT1rhoPrep;
        else
            fprintf('Reading alFree: T1rho prep module not specified.\n');
            paramsMTfree.numT1rhoPrep = 0;
            paramsMTfree.T1rhoDuration = [];
            paramsMTfree.isT1rhoIR = 0;
        end
        if (paramsMTfree.numFA > 0) && (paramsMTfree.moduleLength > 0)
            paramsMTfree.moduleLength = lcm(paramsMTfree.numFA, paramsMTfree.moduleLength);
        else
            paramsMTfree.moduleLength = 1;
        end

%         % Asymmetric echo length (0.5~1.0)
%         if isfield(sWipMemBlock,'adFree') && numel(sWipMemBlock.adFree) >= 12
%             paramsMTfree.asymm_echo = sWipMemBlock.adFree{12} - 0.001;
%         else
%             fprintf('Reading adFree: No Assymmetric echo info.\n');
%             paramsMTfree.asymm_echo = 0.5;
%         end
    catch errormsg
        fprintf(' Automatic extraction of MT parmeters failed.\n');
        fprintf('%s\n', errormsg.message);

        paramsMTfree.moduleLength = 1;
        paramsMTfree.numIRns = 0;
        paramsMTfree.numSRns = 0;
        paramsMTfree.numIRsel = 0;
        paramsMTfree.numSRsel = 0;
        paramsMTfree.numT2prep = 0;
        paramsMTfree.T2prepDuration = [];
        paramsMTfree.numT1rhoPrep = 0;
        paramsMTfree.T1rhoDuration = [];
        paramsMTfree.asymm_echo = 0.5;
    end

end

% for spiralVIBE_MT
if strcmp(params.Trajectory,'Spiral') && numel(sWipMemBlock.alFree)>38
    Nalfree = numel(sWipMemBlock.alFree);
    for n = 1:Nalfree
        if isempty(sWipMemBlock.alFree{n})
            sWipMemBlock.alFree{n} = 0;
        end
    end
    paramsMTfree.SGBlock = sWipMemBlock.alFree{39} + 1;
    paramsMTfree.linesPerShot = sWipMemBlock.alFree{40};
    paramsMTfree.lEchoSpacing = paramsMTfree.alTR_seconds;
    paramsMTfree.alTR_seconds = paramsMTfree.lEchoSpacing*paramsMTfree.linesPerShot;
    paramsMTfree.numFA  = numel(twixObj.hdr.MeasYaps.adFlipAngleDegree);
    for n = 1:paramsMTfree.numFA 
        paramsMTfree.flipAngleArray(n) = [twixObj.hdr.MeasYaps.adFlipAngleDegree{n}];
    end
    if sum(diff(paramsMTfree.flipAngleArray)) == 0
        paramsMTfree.numFA = 1;
        paramsMTfree.flipAngleArray = twixObj.hdr.MeasYaps.adFlipAngleDegree{1};
    end
    paramsMTfree.moduleLength = paramsMTfree.numFA;

    if numel(sWipMemBlock.alFree) > 42
        paramsMTfree.numIRns  = (twixObj.hdr.Meas.ucInversion == 2) * sWipMemBlock.alFree{43};
        paramsMTfree.numIRsel = (twixObj.hdr.Meas.ucInversion == 1) * sWipMemBlock.alFree{43};
    else
        paramsMTfree.numIRns  = (twixObj.hdr.Meas.ucInversion == 2);
        paramsMTfree.numIRsel = (twixObj.hdr.Meas.ucInversion == 1);
    end
    if numel(sWipMemBlock.alFree) > 47
        paramsMTfree.numT2prep = sWipMemBlock.alFree{44};
        paramsMTfree.T2prepDuration = [];
        for n = 1:paramsMTfree.numT2prep
            paramsMTfree.T2prepDuration(n) = sWipMemBlock.alFree{44+n};
        end
    else
        paramsMTfree.numT2prep = 0;
        paramsMTfree.T2prepDuration = [];
    end
    if numel(sWipMemBlock.alFree) > 48
        paramsMTfree.isT2IR = sign(sWipMemBlock.alFree{49});
    else
        paramsMTfree.isT2IR = 0;
    end
    try
        paramsMTfree.moduleLength = lcm((paramsMTfree.numIRns+paramsMTfree.numIRsel+paramsMTfree.numT2prep), paramsMTfree.moduleLength);
    catch
    end
elseif strcmp(params.Trajectory,'Spiral')
    paramsMTfree.lEchoSpacing = paramsMTfree.alTR_seconds;
    paramsMTfree.moduleLength = 1;
    paramsMTfree.SGBlock = 1;
    paramsMTfree.linesPerShot = paramsMTfree.lSegments;
    paramsMTfree.numIRns  = (twixObj.hdr.Meas.ucInversion == 2);
    paramsMTfree.numIRsel = (twixObj.hdr.Meas.ucInversion == 1);
end

% for hard-coded sequence in VB19B -->
if strcmp(VerString,'B19')
    paramsMTfree.numSRns = 0;
    paramsMTfree.numIRns = 1;
    paramsMTfree.numT2prep = 4;
    paramsMTfree.T2prepDuration = [30 40 50 60];
    if numel(twixObj.hdr.MeasYaps.adFlipAngleDegree) > 1 && (twixObj.hdr.MeasYaps.adFlipAngleDegree{1} ~= twixObj.hdr.MeasYaps.adFlipAngleDegree{2})
        paramsMTfree.numFA = 2;
    else
        paramsMTfree.numFA = 1;
    end
    for n = 1:paramsMTfree.numFA
        paramsMTfree.flipAngleArray(n) = twixObj.hdr.MeasYaps.adFlipAngleDegree{n};
    end
    paramsMTfree.moduleLength = lcm(paramsMTfree.numFA, 5);
    paramsMTfree.MBfactor = 1;
    paramsMTfree.isT2IR = 1;
    paramsMTfree.isT1rhoIR = 0;
end
% <-- for hard-coded sequence in VB19B

if ~isfield(paramsMTfree,'isCEST')
    paramsMTfree.isCEST = 0;
end
if paramsMTfree.isCEST
    paramsMTfree.ScanType = 'CEST';
elseif paramsMTfree.numSRns == 0 && paramsMTfree.numIRns == 0 && paramsMTfree.numT2prep == 0 && paramsMTfree.numT1rhoPrep == 0
    paramsMTfree.ScanType = 'Cine';
elseif paramsMTfree.numIRns > 0 && paramsMTfree.numT2prep == 0 && paramsMTfree.numT1rhoPrep == 0 && paramsMTfree.numFA == 1
    paramsMTfree.ScanType = 'IR';
elseif paramsMTfree.numIRns > 0 && paramsMTfree.numT2prep == 0 && paramsMTfree.numT1rhoPrep == 0 && paramsMTfree.numFA > 1
    paramsMTfree.ScanType = 'IR_VFA';
elseif paramsMTfree.numT2prep > 0 && paramsMTfree.numT1rhoPrep == 0 && paramsMTfree.numFA == 1
    paramsMTfree.ScanType = 'T2IR';
elseif paramsMTfree.numT2prep > 0 && paramsMTfree.numT1rhoPrep == 0 && paramsMTfree.numFA > 1
    paramsMTfree.ScanType = 'T2IR_VFA';
elseif paramsMTfree.numT2prep == 0 && paramsMTfree.numT1rhoPrep > 0 && paramsMTfree.numFA == 1
    paramsMTfree.ScanType = 'T1rho';
elseif paramsMTfree.numT2prep == 0 && paramsMTfree.numT1rhoPrep > 0 && paramsMTfree.numFA > 1
    paramsMTfree.ScanType = 'T1rho_VFA';
elseif paramsMTfree.numT2prep > 0 && paramsMTfree.numT1rhoPrep > 0 && paramsMTfree.numFA == 1
    paramsMTfree.ScanType = 'T1rhoT2IR';
elseif paramsMTfree.numT2prep > 0 && paramsMTfree.numT1rhoPrep > 0 && paramsMTfree.numFA > 1
    paramsMTfree.ScanType = 'T1rhoT2IR_VFA';
elseif paramsMTfree.numSRns > 0
    paramsMTfree.ScanType = 'SR';
end 
