% 2D Spiral RF Pulse
% Resolution = 1cm
% FOV = 16cm
clear all; clc; close all;

%% Design of 2D Spiral Gradients
gamma = 2*pi*4.257;  % krad/G
samples = 1024;
kmax = 0.5;     % cycles/cm
N = 8;

% Constant angular design
t = (1:samples)/samples;
ka = kmax * (1-t) .* exp(1i*2*pi*N*(1-t));

% Constant slew-rate design
ks = csg(ka,4,15);

figure(1); plot(ka);
grid on; xlabel('k_x (cycles/cm)'); ylabel('k_y (cycles/cm)');

%% Design of RF Pulse
pulseduration = 5; %2.985;      % in msec
dt = pulseduration/samples; 
% Convert from k-space to gradients
g  = ktog(ks,dt);
t = (1:length(g))*dt;       % in msec

rf = abs(g);
rf = rf/sum(rf);
b1 = rfscaleg(rf, pulseduration);          % Scaled to Gauss

figure(2); plot(t,b1);
grid on; xlabel('Time (msec)'); ylabel('B_1 (Gauss)');

%% 1D Bloch Simulation of 2D Pulse
x = (-16:.1:16);        % in cm
dt = pulseduration/samples/1e3; 
t = (1:length(g))*dt;   % in usec

% Bloch Simulation
[mx, my, mz] = bloch(b1(:),[real(g(:)) imag(g(:))],t(:),1,.2,0,x(:),0); 

% Transverse Magnetization
mxy = mx+1i*my;

figure(3); plot(x,abs(mxy));
grid on; xlabel('x position (cm)');

%% 2D Bloch Simulation of 2D Pulse
x = (-16:.1:16); % in cm
y = (-16:.1:16); % in cm
[yall, xall] = meshgrid(y,x);

dt = pulseduration/samples/1e3; 
t = (1:length(g))*dt;   % in usec

% Bloch Simulation
[mx,my,mz] = bloch(b1(:),[real(g(:)) imag(g(:))],t(:),1,.2,0,[xall(:) yall(:)],0); 

% Transverse Magnetization
mxy = mx+1i*my;
mxy = reshape(mxy,length(x),length(y));

figure(4); imshow(abs(mxy),[]);
xlabel('y position(cm)'); ylabel('x position (cm)');

