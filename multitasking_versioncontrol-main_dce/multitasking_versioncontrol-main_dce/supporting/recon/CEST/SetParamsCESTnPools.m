function ParamsCEST = SetParamsCESTnPools(nPools,CESTMetabolite,params)
% exchange using least-squares regression
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Syntax
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% Inputs
%
%       magField, Magnetic Field [Tesla]
%       satTime, Saturation Time [Seconds]
%       satPower, Saturation Power [Micro Tesla]
%       nPools, Number of Pools
%
%       x-Data, Saturation Offsets [ppm]                           
%       y-Data, Z-magnetization (Mz/Mo) [%]
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Adapted from the CESTMRI.mat tool by Authors:
% Jaden Johnston                    Julio Cardenas-Rodriguez
% University of Arizona             University of Arizona
% jadenjohnston@email.arizona.edu   cardenaj@email.arizona.edu
%                       
% www.cardenaslab.org/resources
% v1.0 05/27/2016

% Code
%% Allocate variables

gyroRatH = 42.576;                        % Gyromagnetic Ratio [MHz/T]
concH = 111;                              % Molarity [M]

% water pool
T1 = [1];                                % Seconds
T2 = [0.07];                             % Seconds
fConc = [1];                               % Fractional concentration
ExchangeRates = [0];                       % Exchange rates from pools to water [Hz]
PoolOffsets = [0];                         % ppm
x0  = [0.67,  0,  1.8];
xlb = [0.02, -1,  0.3];
xub = [1,     1, 10];

% MT pool
T1 = [T1 1];                             % Seconds
T2 = [T2 4.5e-5];                              % Seconds
fConc = [fConc 15/concH];                        
ExchangeRates= [ExchangeRates 30];                 % Hz
PoolOffsets= [PoolOffsets -1];                     % ppm
x0  = [x0,  0.15, -1,   40];
xlb = [xlb, 0,    -2.5, 30];
xub = [xub, 0.5,   0,   60];

if nPools > 5
    fprintf(2,'%d-pool model not (yet) supported in this recon package.',nPools);
    nPools = 5;
elseif nPools < 2
    fprintf(2,'%d-pool model not (yet) supported in this recon package.',nPools);
    nPools = 2;
end

if nPools >= 3 
    % CEST metabolite pool
    if strcmp(CESTMetabolite,'APT')
        T1 = [T1 1];                             % Seconds
        T2 = [T2 0.1];                             % Seconds
        fConc = [fConc 72e-3/concH];                       
        ExchangeRates= [ExchangeRates 100];                    % Hz
        PoolOffsets= [PoolOffsets 3.5];                     % ppm
        x0  = [x0,  0.05, 3.5,  4.5];
        xlb = [xlb, 0,    3.2,  0.4];
        xub = [xub, 0.2,  3.8,  6.0];
    elseif strcmp(CESTMetabolite,'Creatine')
        T1 = [T1 1];                             % Seconds
        T2 = [T2 0.1];                             % Seconds
        fConc = [fConc 72e-3/concH];                       
        ExchangeRates= [ExchangeRates 100];                 % Hz
        PoolOffsets= [PoolOffsets 1.8];                     % ppm
        x0  = [x0,  0.05, 1.8,  4.5];
        xlb = [xlb, 0,    1.5,  0.4];
        xub = [xub, 0.2,  2.1,  6.0];
    elseif strcmp(CESTMetabolite,'Glucose')
        T1 = [T1 1];                             % Seconds
        T2 = [T2 0.1];                             % Seconds
        fConc = [fConc 72e-3/concH];                       
        ExchangeRates= [ExchangeRates 100];                    % Hz
        PoolOffsets= [PoolOffsets 0.8];                     % ppm
        x0  = [x0,  0.05, 0.8,  4.5];
        xlb = [xlb, 0,    0.5,  0.4];
        xub = [xub, 0.2,  1.1,  6.0];
    elseif strcmp(CESTMetabolite,'GAG')
        T1 = [T1 1];                             % Seconds
        T2 = [T2 0.1];                             % Seconds
        fConc = [fConc 72e-3/concH];                       
        ExchangeRates= [ExchangeRates 100];                    % Hz
        PoolOffsets= [PoolOffsets 1.0];                     % ppm
    else
        T1 = [T1 1];                             % Seconds
        T2 = [T2 0.1];                             % Seconds
        fConc = [fConc 72e-3/concH];                       
        ExchangeRates= [ExchangeRates 100];                    % Hz
        PoolOffsets= [PoolOffsets 3.5];                     % ppm    
    end
end

if nPools >= 4      % add NOE pool
    T1 = [T1 1];                             % Seconds
    T2 = [T2 0.005];                             % Seconds
    fConc = [fConc 500e-3/concH];                        
    ExchangeRates= [ExchangeRates 16];                      % Hz
    PoolOffsets= [PoolOffsets -3.5];                     % ppm
    x0  = [x0,  0.1, -3.5,  4.5];
    xlb = [xlb, 0,   -4.0,  1.0];
    xub = [xub, 0.2, -3.0,  6.0];
end

if nPools >= 5      % add amine pool
    T1 = [T1 1];                               % Seconds
    T2 = [T2 0.1];                             % Seconds
    fConc = [fConc 20e-3/concH];                        
    ExchangeRates= [ExchangeRates 5000];              % Hz
    PoolOffsets= [PoolOffsets 2];                     % ppm
    x0  = [x0,  0.05, 2.0,  1.0];
    xlb = [xlb, 0,    1.8,  0.4];
    xub = [xub, 0.1,  2.2,  6.0];
end

M0 = zeros(nPools*3+1,1);
M0(nPools*2+1:nPools*3) = fConc;
M0(end) = 1;

ParamsCEST.M0 = M0;
ParamsCEST.R1 = 1./ T1(:);                           % R1 Relaxation Rates [Hz]
ParamsCEST.R2 = 1./ T2(:);                           % R2 Relaxation Rates [Hz]
ParamsCEST.fConc = fConc(:);                         % Fractional concentration of the solute protons. Water has ratio of 1
ParamsCEST.RateNA = ExchangeRates(:);            
ParamsCEST.RateAN = ExchangeRates(:) .* fConc(:);    
ParamsCEST.RateNA(1) = sum(ExchangeRates(2:end).*fConc(2:end));                           % Sum of exchange rates from water to pools [Hz]
ParamsCEST.PoolOffsetsppm = PoolOffsets(:);
ParamsCEST.PoolOffsets = PoolOffsets(:).* 1e-6 .* params.lResonanceFrequency .* 2.* pi;             % Larmor frequencies [rad/s] of all the pool
ParamsCEST.lResonanceFrequency = params.lResonanceFrequency;
ParamsCEST.gyroRatH = gyroRatH;
ParamsCEST.nPools = nPools;
ParamsCEST.x0  = x0;
ParamsCEST.xlb = xlb;
ParamsCEST.xub = xub;

