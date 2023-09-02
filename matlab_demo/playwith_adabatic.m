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
beta1 = 750;   % Frequency modulation parameter [rad/s]
pulseWidth = 10.24*2;   % RF pulse duration [ms] % According to siemens 3T
A0 = 0.12;   % Peak B1 amplitude [Gauss].

%%%%%%

nSamples = 512;        % number of samples in the RF pulse
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
T1_myo   = 1000;      % [ms]
T2_myo   = 45;      % [ms]
T1_fat15t   = 260;      % [ms]
T2_fat15t   = 60;      % [ms]
T1_fat3t   = 400;      % [ms]
T2_fat3t   = 100;      % [ms]

f_max = 2000;   % off-resonance frequency range [Hz]
freq_range  = linspace(-f_max,f_max,1000);   % off-resonance frequency range [Hz]

grad_pulse = zeros(1,length(rf_pulse));  
mod = 2; 
% 0 -- /* Transient state AND end time points. */
% 1 -- /* Steady state? AND end time points. */
% 2 -- /* Transient state AND record all time points. */
% 3 -- /* Steady state AND record all time points. */

[mx1,my1,mz1] = bloch(rf_pulse,  grad_pulse,  dt,  ...
    T1_myo/1000,  T2_myo/1000,  freq_range,  0,  mod,  0,  0,  0);
[mx2,my2,mz2] = bloch(rf_pulse,  grad_pulse,  dt,  ...
    T1_fat15t/1000,  T2_fat15t/1000,  freq_range,  0,  mod,  0,  0,  0);
[mx3,my3,mz3] = bloch(rf_pulse,  grad_pulse,  dt,  ...
    T1_fat3t/1000,  T2_fat3t/1000,  freq_range,  0,  mod,  0,  0,  0);

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
set(gca, 'FontSize', 18);

% Plotting the RF frequency modulation function.
subplot(212);
plot(tim_bloch.*1000,w,'k','LineWidth',LineWidthVal); %title('RF Frequency Modulation'); 
grid on;
title('RF Frequency Modulation Function'); xlabel('Time (ms)'); ylabel('RF Carrier Frequency (Hz)')
set(gca, 'FontSize', 18);
%%
[tt1, idx] = min(abs(freq_range - 0));
[tt2, idx_fat15t] = min(abs(freq_range - 210));
[tt3, idx_fat3t] = min(abs(freq_range - 420));
% Plotting the longitudinal magnetization to see the inversion profile as a
% function of resonance frequency
figure();
plot(freq_range, mz1(end, :),'LineWidth',LineWidthVal);
title('Inversion Profile'); xlabel('Frequency (Hz)'); ylabel('M_z'); grid on;
v = axis; axis([v(1) v(2) -1.05 1.05]);
hold on;


plot(freq_range, mz2(end, :),'LineWidth',LineWidthVal);
plot(freq_range, mz3(end, :),'LineWidth',LineWidthVal);

plot([0 0], [-1.05 1.05], 'k', 'LineWidth', 0.5)
plot([210 210], [-1.05 1.05], 'k','LineWidth', 0.5)
plot([420 420], [-1.05 1.05], 'k', 'LineWidth', 0.5)

plot(freq_range(idx), mz1(end,idx), 'r*', 'MarkerSize', 12);
plot(freq_range(idx_fat15t), mz2(end,idx_fat15t),  'r*', 'MarkerSize', 12);
plot(freq_range(idx_fat3t), mz3(end,idx_fat3t),  'r*', 'MarkerSize', 12);

legend({'T1=1000, T2=45', 'T1=260, T2=60', 'T1=400, T2=100', '0 Hz', '210 Hz', '420 Hz'});
set(gca, 'FontSize', 18);
%%
loops = numel(tim_bloch);
figure();
plot3(mx1(:, idx), my1(:, idx), mz1(:, idx), '.-'); axis equal;
grid on;
hold on;
for i = 1:loops
    hline = plot3([0 mx1(i, idx)], [0 my1(i, idx)], [0 mz1(i, idx)], 'r');
    xlabel('Mx'); ylabel('My');zlabel('Mz')
    xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);
    pause(.02)
    if i ~= loops
        set(hline,'visible','off');
    end
end
%% Fat 3T
loops = numel(tim_bloch);
figure();
plot3(mx3(:, idx_fat3t), my3(:, idx_fat3t), mz3(:, idx_fat3t), '.-'); axis equal;
grid on;
hold on;
for i = 1:loops
    hline = plot3([0 mx3(i, idx_fat3t)], [0 my3(i, idx_fat3t)], [0 mz3(i, idx_fat3t)], 'r');
    xlabel('Mx'); ylabel('My');zlabel('Mz')
    xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);
    pause(.02)
    if i ~= loops
        set(hline,'visible','off');
    end
end
%% Fat 1.5T
loops = numel(tim_bloch);
figure();
plot3(mx2(:, idx_fat15t), my2(:, idx_fat15t), mz2(:, idx_fat15t), '.-'); axis equal;
grid on;
hold on;
for i = 1:loops
    hline = plot3([0 mx2(i, idx_fat15t)], [0 my2(i, idx_fat15t)], [0 mz2(i, idx_fat15t)], 'r');
    xlabel('Mx'); ylabel('My');zlabel('Mz')
    xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);
    pause(.02)
    if i ~= loops
        set(hline,'visible','off');
    end
end

%% Try the function
%%% User inputs:
adiabatic.mu = 5;   % Phase modulation parameter [dimensionless]
adiabatic.beta1 = 750;   % Frequency modulation parameter [rad/s]
adiabatic.pulseWidth = 10.24*2;   % RF pulse duration [ms] % According to siemens 3T
adiabatic.A0 = 0.12;   % Peak B1 amplitude [Gauss].

T1 = 1000;
T2 = 45;
df = 0;

[mx0, my0, mz0] = AdiabaticPulse(adiabatic, T1, T2, df);
ad_dt = adiabatic.pulseWidth / length(mx0);
