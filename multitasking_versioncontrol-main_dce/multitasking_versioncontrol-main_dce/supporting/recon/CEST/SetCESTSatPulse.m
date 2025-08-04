function SatPulseCEST = SetCESTSatPulse(SatFA,B1FileName)

if nargin < 2
    B1FileName = 'B1_td30_fa100.mat';
end

filePath = fullfile(fileparts(mfilename('fullpath')),B1FileName);

% Assign parameters for CEST saturation pulse
load(filePath);				% 30ms / 100deg

% CEST pulse
SatPulseCEST.B1train = B1_td30_unit_100 * SatFA/B1_FA;     % abs B1 value [T]
SatPulseCEST.B1train = [SatPulseCEST.B1train(:); zeros(floor(B1_spoilerDuration/B1_dt),1)];
SatPulseCEST.dt = B1_dt;					% time step [s]
% SatPulseCEST.duration = SatPulseCEST.dt*numel(SatPulseCEST.B1train);	% duration of a single pulse [s]
SatPulseCEST.offsetppm = 0;
