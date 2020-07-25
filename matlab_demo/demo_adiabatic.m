% This script will allow you to simulate the hyperbolic secant adiabatic 
% inversion pulse and view the RF waveforms and resultant magnetization 
% profiles. You should explore how the pulse waveform and magnetization 
% profiles are altered by changing the parameters.
%
% Refer to the Handbook MRI Pulse Sequences by Bernstein et al., page 194,
% for the equations for the hyperbolic secant pulse.
%
% To use this script, make sure you have the Bloch simulation code by Dr. Brian Hargreaves. 
% These files have already been provided to you, but you can also download them from:
% http://www-mrsrl.stanford.edu/~brian/blochsim/
% 
% The compiled mex files for Mac OS X and Linux are provided,
% but if you are using Windows, you need to compile your own mex
% file from Matlab by running 
% mex bloch.c 
% in the Matlab command prompt. % This should produce a new mex file such as bloch.mexw64
%
% To get the most out of this simulation, read the instructions in readme_first.m and go through the exercises. 
%%
clear all; close all;

%%% User inputs:
mu = 5;   % Phase modulation parameter [dimensionless]
beta1 = 672;   % Frequency modulation parameter [rad/s]
pulseWidth = 10.24;   % RF pulse duration [ms]
A0 = 0.12;   % Peak B1 amplitude [Gauss].

%%%%%%

nSamples = 512;        % number of samples in the RF pulse
dt = pulseWidth/nSamples/1000;  % time step, [seconds]
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
T1_value   = 10000;      % [ms]
T2_value   = 10000;      % [ms]

f_max = 4000;   % off-resonance frequency range [Hz]
freq_range  = linspace(-f_max,f_max,1000);   % off-resonance frequency range [Hz]

grad_pulse = zeros(1,length(rf_pulse));  
mod = 0;
[mx1,my1,mz1] = bloch(rf_pulse,  grad_pulse,  dt,  ...
    T1_value/1000,  T2_value/1000,  freq_range,  0,  mod,  0,  0,  0);


%% Create line figures:

LineWidthVal = 1.5;

% Plotting the RF amplitude modulation function.
% rf_pulse is multiplied by 100 to express the units in microTesla
figure(1);
subplot(211);
plot(tim_bloch.*1000, abs(rf_pulse).*100,'k','LineWidth',LineWidthVal);
grid on;
title('RF Amplitude Modulation Function'); xlabel('Time (ms)'); ylabel('B1 (\muT)');
v = axis; axis([v(1) v(2) v(3) v(4)+1]);

% Plotting the RF frequency modulation function.
subplot(212);
plot(tim_bloch.*1000,w,'k','LineWidth',LineWidthVal); %title('RF Frequency Modulation'); 
grid on;
title('RF Frequency Modulation Function'); xlabel('Time (ms)'); ylabel('RF Carrier Frequency (Hz)')

% Plotting the longitudinal magnetization to see the inversion profile as a
% function of resonance frequency
figure(2);
plot(freq_range, mz1,'k','LineWidth',LineWidthVal);
title('Inversion Profile'); xlabel('Frequency (Hz)'); ylabel('M_z'); grid on;
v = axis; axis([v(1) v(2) -1.05 1.05]);


