function [mx0, my0, mz0] = AdiabaticPulse(adiabatic, T1, T2, df, M0)

%%% User inputs:
mu = adiabatic.mu;   % Phase modulation parameter [dimensionless]
beta1 = adiabatic.beta1;   % Frequency modulation parameter [rad/s]
pulseWidth = adiabatic.pulseWidth;   % RF pulse duration [ms] % According to siemens 3T
A0 = adiabatic.A0;   % Peak B1 amplitude [Gauss].

%%%%%%

nSamples = 1024;        % number of samples in the RF pulse
dt = pulseWidth/nSamples/1000;  % time step, [seconds] % 2e-5
tim_sech = linspace(-pulseWidth/2,pulseWidth/2,nSamples)./1000';  
% time scale to calculate the RF waveforms in seconds. 

% Amplitude modulation function B1(t):
B1 = A0.* sech(beta1.*tim_sech);

% Carrier frequency modulation function w(t):
w = -mu.*beta1.*tanh(beta1.*tim_sech)./(2*pi);
% The 2*PI scaling factor at the end converts the unit from rad/s to Hz

% Phase modulation function phi(t):
phi = mu .* log(sech(beta1.*tim_sech));

% Put together complex RF pulse waveform:
rf_pulse = B1 .* exp(1i.*phi);

% Generate a time scale for the Bloch simulation:
tim_bloch = [0:(nSamples-1)]*dt;

%% The Bloch simulator requires a gradient input. For our simulation,
% gradient will be zero, as we are simulating a non-selective RF pulse.

f_max = 2000;   % off-resonance frequency range [Hz]
freq_range  = linspace(-f_max,f_max,1000);   % off-resonance frequency range [Hz]

mx0 = repmat([M0(1)], [1, length(freq_range)]);
my0 = repmat([M0(2)], [1, length(freq_range)]);
mz0 = repmat([M0(3)], [1, length(freq_range)]);

grad_pulse = zeros(1,length(rf_pulse));  
mod = 2; 
% 0 -- /* Transient state AND end time points. */
% 1 -- /* Steady state? AND end time points. */
% 2 -- /* Transient state AND record all time points. */
% 3 -- /* Steady state AND record all time points. */

[mx,my,mz] = bloch(rf_pulse,  grad_pulse,  dt,  ...
    T1/1000,  T2/1000,  freq_range,  0,  mod,  mx0,  my0,  mz0);

[tt, idx] = min(abs(freq_range - df));
mx0 = mx(:,idx); my0 = my(:,idx); mz0 = mz(:,idx);
end