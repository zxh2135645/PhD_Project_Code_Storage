clear all; clc; close all;

%% Design of Windowed Sinc RF Pulses
tbw = 4;
samples = 512;

rf = wsinc(tbw, samples);

%% Plot RF Amplitude
tbw = 4;
samples = 512;

rf = (pi/2)*wsinc(tbw,samples);

pulseduration = 1; % in ms
rfs = rfscaleg(rf, pulseduration); % Scaled to Gauss

dt = pulseduration/samples;
t = (1:length(rfs))*dt;  % in msec

b1    = [rfs zeros(1,samples/2)];                 % in Gauss
g     = [ones(1,samples) -ones(1,samples/2)];     % in G/cm
t_all = (1:length(g))*dt;  % in msec

figure(1); 
subplot(311); plot(t_all,b1); grid on; 
xlabel('Time (msec)'); ylabel('B_1 (Gauss)');

subplot(312); plot(t_all,g); grid on; 
xlabel('Time (msec)'); ylabel('G_z (Gauss/cm)');

%% Simulate Slice Profile using Bloch Simulation
x = (-4:.1:4);          % in cm
f = 0;                  % in Hz
dt = pulseduration/samples/1e3;
t = (1:length(b1))*dt;  % in usec

% Bloch Simulation
[mx,my,mz] = bloch(b1(:),g(:),t(:),1,.2,f(:),x(:),0); 

% Transverse Magnetization
mxy = mx+1i*my;

figure(1); subplot(313); plot(x,abs(mxy)); 
grid on; xlabel('z position (cm)');

%% Simulate Slice Profile using Small Tip Approximation



