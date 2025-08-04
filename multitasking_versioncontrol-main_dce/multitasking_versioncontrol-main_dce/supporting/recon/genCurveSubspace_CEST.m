function [curvePhi, curvePhi_binning_hybrid, temporalBasis, fitParams, curvePhi_alt] = genCurveSubspace_CEST(params, fitParams, reconOptions, dataArray, temporalBasis)

% License check
if ~checkLicense(reconOptions)
    dlg = errordlg('Multitasking license check failed');
    waitfor(dlg);
    return;
end

%%
tic;

fitParams.params = params;
fitParams.reconOptions = reconOptions;
extractVarFromStruct(fitParams);
extractVarFromStruct(params);
extractVarFromStruct(reconOptions);

vec = @(x) x(:);

alphaArray = flipAngleArray(1)*pi/180;
fitParams.alphaArray = alphaArray;

CESTfreqStart = find((CESTSatFreqOffsetppmList==CESTSatFreqOffsetppmList(1))==0,1);
CESTfreqEnd   = find((CESTSatFreqOffsetppmList==CESTSatFreqOffsetppmList(1))==0,1,'last')+1;

curveMaskNav = zeros(size(dataArray.navData,1),1);
curveMaskNav(CESTNumRep*linesPerShot/SGBlock*(CESTfreqStart-1)+1:CESTNumRep*linesPerShot/SGBlock*CESTfreqEnd) = 1;
curveMaskFull = vec(repmat(curveMaskNav',[SGBlock,1]));

temporalBasis.CESTfreqStart = CESTfreqStart;
temporalBasis.CESTfreqEnd   = CESTfreqEnd;
temporalBasis.curveMaskNav  = curveMaskNav;
temporalBasis.curveMaskFull = curveMaskFull;


%% Simulation for T1 relaxation

% Set up CEST pools and saturation pulse
% switch CESTMetabolite
%     case 'APT'
%         metaboliteFreqOffset = 3.5;
%     case 'Creatine'
%         metaboliteFreqOffset = 1.8;
%     case 'Glucose'
%         metaboliteFreqOffset = 0.8;
%     case 'GAG'
%         metaboliteFreqOffset = 1.0;
% end

ParamsCEST   = SetParamsCESTnPools(CESTnPools,CESTMetabolite,params);
SatPulseCEST = SetCESTSatPulse(CESTSatFA);

fitParams.ParamsCEST = ParamsCEST;
fitParams.SatPulseCEST = SatPulseCEST;

% initialize simulation space and signal vector
R1A_list = 1./logspace(log10(.3),log10(3),21);
NR1A = length(R1A_list);

Rkex_list = 20:20:500;    % exchange rate
NRkex = length(Rkex_list);

signalRO = zeros(linesPerShot,CESTNumRep,CESTnumFreqOffsets,NR1A,NRkex);  

% start simulation
tic;
fprintf('Bloch-McConnell simulation: ');
strProgress = sprintf('0/%d ',NR1A); fprintf(strProgress);
for iR1A = 1:NR1A
    prevLength = numel(strProgress);
    strProgress = sprintf('%d/%d... ',iR1A,NR1A);
    fprintf([repmat('\b',1,prevLength) '%s'],strProgress);

    % Bloch simulation of GRE sequence 
    Atr = zeros(3,3,ParamsCEST.nPools);
    Btr = zeros(3,ParamsCEST.nPools);
    for n = 1:ParamsCEST.nPools
        [Atr(:,:,n),Btr(:,n)] = freeprecess(lEchoSpacing,1/ParamsCEST.R1(n),1/ParamsCEST.R2(n),ParamsCEST.PoolOffsets(n)/2/pi);
    end
    spoilerG = [0 0 0;0 0 0;0 0 1];

    for iRkex = 1:NRkex
        ParamsCEST.R1(:) = R1A_list(iR1A); 
        ParamsCEST.RateNA(2) = Rkex_list(iRkex);
        ParamsCEST.RateAN = ParamsCEST.RateNA.*ParamsCEST.fConc; ParamsCEST.RateAN(1) = 0;
        ParamsCEST.RateNA(1) = sum(ParamsCEST.RateAN(2:end));

        M = ParamsCEST.M0;

        % loop through frequency offsets
        for iOffset = 1:CESTnumFreqOffsets
            SatPulseCEST.offsetppm = CESTSatFreqOffsetppmList(iOffset);
            SatPulseRotA = CESTsimuBlochMcConnell_pulse(ParamsCEST,SatPulseCEST);

            % loop through of the sat pulse repetition
            for iCESTrep = 1:CESTNumRep                
                % CEST saturation and spoiler gradients
                M = SatPulseRotA*M;
                M(1:ParamsCEST.nPools*2) = 0;
                
                % loop through readouts
                Mtemp = reshape(M(1:end-1),ParamsCEST.nPools,3)';
                for iRO = 1:linesPerShot
                    signalRO(iRO,iCESTrep,iOffset,iR1A,iRkex) = Mtemp(3,1)*sin(alphaArray(1));
                    for n = 1:ParamsCEST.nPools
                        Mtemp(:,n) = spoilerG*(Atr(:,:,n)*(yrot(alphaArray(1))*Mtemp(:,n)) + Btr(:,n)*ParamsCEST.fConc(n));
                    end
                end
                M(1:end-1) = reshape(Mtemp',[],1);
            end            
        end
    end
end
toc;

signalRO = reshape(signalRO,[],NR1A,NRkex);
figure, plot(signalRO(1:SGBlock:end,round(end/2),round(end/2)),'.-'), title('Simulated signal evolution');

fitParams.CESTdictionary.signal = signalRO;
fitParams.CESTdictionary.R1A  = R1A_list;
fitParams.CESTdictionary.Rkex = Rkex_list;

% save('ssCEST_dictionary.mat','signal_readout');
% save(fullfile(file_string,'ssCEST_dictionary.mat'),'signalRO');


%% Generate curvePhi

% T1 basis (simulation) from masked curve
temp_simu = reshape(signalRO.*curveMaskFull,linesPerShot*CESTNumRep,[]);
if SGBlock >= linesPerShot
    temp_simu = interp1(1:SGBlock:size(temp_simu,1),temp_simu(1:SGBlock:end,:),1:size(temp_simu,1),'pchip','extrap');
end
[curvePhi_simu,~,~] = svd(temp_simu(:,:),'econ'); 

% -- No mask for binning 
% [curvePhi_binning,~,~] = svd(signalRO(:,:),'econ'); 
% Data-driven z-basis
temp = reshape(dataArray.navData,CESTNumRep*linesPerShot/SGBlock,CESTnumFreqOffsets,[]);
[curvez_dd,dd,~] = svd(reshape(permute(temp,[2 1 3]),CESTnumFreqOffsets,[]),'econ'); 

bothPhi_hybrid = kron(curvez_dd(:,1:10),curvePhi_simu(:,1:8)); % prev: 5 & 3
CUk_hybrid = pinv(bothPhi_hybrid(1:SGBlock:end,:))*(dataArray.navData(:,:));

[C1_hybrid,C2_hybrid,~] = svd(CUk_hybrid,'econ');
C_hybrid = C1_hybrid*C2_hybrid; 
curvePhi_binning_hybrid = bothPhi_hybrid*C_hybrid;

% -- Masked curves to reduce rank
% Data-driven z-basis
temp = reshape(dataArray.navData.*curveMaskNav,CESTNumRep*linesPerShot/SGBlock,CESTnumFreqOffsets,[]);
[curvez_dd,~,~] = svd(reshape(permute(temp,[2 1 3]),CESTnumFreqOffsets,[]),'econ'); 

% Hybrid curvePhi
bothPhi_hybrid = kron(curvez_dd(:,1:8),curvePhi_simu(:,1:5)); % prev: 5 & 3
CUk_hybrid = pinv(bothPhi_hybrid(1:SGBlock:end,:))*(dataArray.navData(:,:).*curveMaskNav);
[C1_hybrid,C2_hybrid,~] = svd(CUk_hybrid,'econ');
C_hybrid = C1_hybrid*C2_hybrid; 
curvePhi_hybrid = bothPhi_hybrid*C_hybrid;

% figure, plot(realify(bothPhi_hybrid(:,1)),'.-')
% figure, plot(realify(curvePhi_binning_hybrid(:,1)),'.-')

if flagDataDriven
    curvePhi = curvePhi_simu;
    curvePhi_alt = curvePhi_hybrid;
else
    curvePhi = curvePhi_hybrid;     % reshape(curvePhi_binning_hybrid,linesPerShot*CESTNumRep*Noffsets,[]);
    curvePhi_alt = curvePhi_simu;
end
% figure, plot(realify(curvePhi(:,1:3)),'.-');


%% Function simulates free precession and decay
function [Afp,Bfp] = freeprecess(T,T1,T2,df)
%	over a time interval T, given relaxation times T1 and T2
%	and off-resonance df.  Times in s, off-resonance in Hz.

phi = 2*pi*df*T;	% Resonant precession, radians.
E1 = exp(-T/T1);	
E2 = exp(-T/T2);

Afp = [E2 0 0;0 E2 0;0 0 E1]*zrot(phi);
Bfp = [0 0 1-E1]';


%% Function simulates RF rotation about y axis
function Ry = yrot(phi)
Ry = [ cos(phi) 0   sin(phi);
       0        1   0;
      -sin(phi) 0   cos(phi)];

%% Function simulates RF rotation about z axis
function Rz = zrot(phi)

Rz = [ cos(phi)    -sin(phi)    0;
       sin(phi)     cos(phi)    0; 
       0            0           1];

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