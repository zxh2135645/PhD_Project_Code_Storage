function [SatPulseRotA, Mz] = CESTsimuBlochMcConnell_pulse(ParamsCEST,SatPulseCEST)
% Bloch-McConnell simulation of the propagation matrix for the pool & saturation pulse 
% defined in ParamsCEST & SatPulseCEST
%
% Input:
%
%   ParamsCEST has the following field:
%   - R1                    T1 Relaxation Rates [Hz]
%   - R2                    T2 Relaxation Rates [Hz]
%   - fConc                 Fractional concentration of the solute protons
%   - RateNA                Exchange rates from pools to water [Hz]
%   - RateAN                Exchange rates from water to pools [Hz]
%   - PoolOffsets           Larmor frequencies [rad/s]
%   - lResonanceFrequency   Larmor frequency [Hz]
%
%   SatPulseCEST has the following field:
%   - B1train               abs B1 value waveform of the pulse [T]
%   - dt 			        time step [s]
%   - offsetppm             frequency offset of the pulse [ppm]
%
%
% Output:
%
%   SatPulseRotA            propogation matrix of the whole pulse
%   Mz                      water z component after saturation
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Reference
% Murase K, and Tanki Nobuyoshi. Magnetic Resonance Imaging 29 (2011) 126; 131
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Code modified from CESTMRI.mat toolbox
% https://github.com/JCardenasRdz/CESTMRI.mat
%
% Authors:
% Jaden Johnston                    Julio Cardenas-Rodriguez
% University of Arizona             University of Arizona
% jadenjohnston@email.arizona.edu   cardenaj@email.arizona.edu
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

extractVarFromStruct(ParamsCEST);
extractVarFromStruct(SatPulseCEST);

SatPulseOffset = offsetppm*1e-6 * lResonanceFrequency * 2*pi;   % Offset frequency of the pulse [rad/s]

% initialize magnetization and propogation matrix
M = M0;
SatPulseRotA = eye(numel(M0));

% start simulation
for tidx = 1:numel(B1train)

    w1 = B1train(tidx) * gyroRatH*1e6 * 2*pi;  % Nutation rate of RF irradiation [rad/s]

    % Allocate Propagation Matrix
    %{ 
        Matrix A is in the form (see Fig 4 in Murase):
             A1     A5     0     0
            -A5     A1    A4     0
              0    -A4    A2    A3
              0      0     0     0
    %}
        
    CA = cell(4);
    
    A1 = -diag(R2 + RateNA);
         A1(1,2:end) = RateNA(2:end);
         A1(2:end,1) = RateAN(2:end);
    A2 = -diag(R1 + RateNA);
         A2(1,2:end) = RateNA(2:end);
         A2(2:end,1) = RateAN(2:end);
    A3 = R1.* fConc;
    A4 = w1*eye(nPools);
    
    CA{1,1} = A1;
    CA{2,2} = A1;
    CA{3,3} = A2;
    CA{3,4} = A3;
    CA{2,3} = A4;
    CA{3,2} = -A4;
    
    CA{1,3} = zeros(nPools);
    CA{1,4} = zeros(nPools,1);
    CA{2,4} = zeros(nPools,1);
    CA{3,1} = zeros(nPools);
    CA{4,1} = zeros(1,nPools);
    CA{4,2} = zeros(1,nPools);
    CA{4,3} = zeros(1,nPools);
    CA{4,4} = 0;
        
    A5 = diag(PoolOffsets - SatPulseOffset);
    
    CA{1,2} = A5;
    CA{2,1} = -A5;
    
    A = cell2mat(CA);

    SatPulseRotA = expm(A*dt) * SatPulseRotA;
    M = expm(A*dt) * M;
end

Mz = M(nPools*2+1);
