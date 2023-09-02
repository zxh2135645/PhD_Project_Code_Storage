clear all; clc; close all;

%% Design of Windowed Sinc RF Pulses
tbw = 4;
samples = 512;

rf = wsinc(tbw, samples);

figure(1); 
plot(rf); grid on; 

%% Plot RF Amplitude
tbw = 4;
samples = 512;

rf = (pi/2)*wsinc(tbw,samples);

pulseduration = 1; % in ms
rfs = rfscaleg(rf, pulseduration); % Scaled to Gauss

dt = pulseduration/samples;
t = (1:length(rfs))*dt;  % in msec

figure(2); 
plot(t,rfs); grid on; 
xlabel('Time (msec)'); ylabel('B_1 (Gauss)');

%% Simulate Slice Profile
tbw = 4;
samples = 512;

rf = (pi/2)*wsinc(tbw,samples);
pulseduration = 1; % ms

rfs = rfscaleg(rf, pulseduration);          % Scaled to Gauss
b1  = [rfs zeros(1,samples/2)];              % in Gauss
g   = [ones(1,samples) -ones(1,samples/2)];   % in G/cm

x = (-4:.1:4);          % in cm
f = (-250:5:250);       % in Hz
dt = pulseduration/samples/1e3;
t = (1:length(b1))*dt;  % in usec

% Bloch Simulation
[mx,my,mz] = bloch(b1(:),g(:),t(:),1,.2,f(:),x(:),0); 

% Transverse Magnetization
mxy = mx+1i*my;

figure(3); imshow(abs(mxy),[]);
xlabel('Off-resonance'); ylabel('z position (cm)');

figure(4); plot(x,abs(mxy(:,round(end/2)))); 
grid on; xlabel('z position (cm)');



