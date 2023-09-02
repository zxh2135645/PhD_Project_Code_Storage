clear all; clc; close all;
addpath('D:\src\M219\')
gamma=42.58e6;% Hz/T
Gz=2.8*10e-4;% T/cm
delta_z=1.2;% cm
delta_omega=gamma*Gz*delta_z;% Hz


%% Design of Windowed Sinc RF Pulses
% tbw = 4*10e3;
% samples = 512;
% 
% rf = wsinc(tbw, samples);

%% Plot RF Amplitude
tbw = 10;
beta = 672;%rad/s
u = 5;
% samples = 512;

% rf = (pi/2)*wsinc(tbw,samples);
% %delta_omega=
% pulseduration = 10e3*tbw/delta_omega; % in ms
% rfs = rfscaleg(rf, pulseduration); % Scaled to Gauss
% 
% dt = pulseduration/samples;
% t = (1:length(rfs))*dt;  % in msec
% 
% b1    = [rfs zeros(1,samples/2)];                 % in Gauss
% g     = 4*[ones(1,samples) -ones(1,samples/2)];     % in G/cm
% t_all = (1:length(g))*dt;  % in msec
% 
% figure(1); 
% subplot(311); plot(t_all,b1); grid on; 
% xlabel('Time (msec)'); ylabel('B_1 (Gauss)');
% 
% subplot(312); plot(t_all,g); grid on; 
% xlabel('Time (msec)'); ylabel('G_z (Gauss/cm)');
% 

pulseduration  = 10.240; %ms
dt = 0.002;
t  = (0:dt:10.24)./1000;%s
B1 = ((100.*sech(beta.*(t-0.00512))).^(1+1j*u))./100;%Gauss
G  = zeros(1,length(t));
G1  = .3.*ones(1,length(t));
f  = (0:50:2000);
x = (-2:.1:2);
%x = ones(1,length(t));
y = ones(1,length(t));


mx0 = zeros(1,length(x));
my0 = zeros(1,length(x));
mz0 = 1.*ones(1,length(x));



figure(1); 
subplot(411); plot(t.*1000,abs(B1)); grid on; 
xlabel('Time (msec)'); ylabel('B_1 (Gauss)');

omega = -4.*beta.*tanh(beta.*(t-0.00512));

figure(1); subplot(412); plot(t.*1000,omega); grid on; 
xlabel('Time (msec)'); ylabel('G_z (Gauss/cm)');


%

[mx,my,mz] = bloch_adiabatic(B1(:),G1(:),dt/1000,1.535,.058,f,x(:),0,omega); 


mxy = mx+1i*my;

figure(1); subplot(411); plot(x,abs(mxy)); 
grid on; xlabel('z position (cm) mxy');

figure(1); subplot(412); plot(x,mz); 
grid on; xlabel('z position (cm)');

figure(1); subplot(413); plot(f,mz(26,:));
%%
[mx,my,mz] = bloch_adiabatic(B1(:),G1(:),dt/1000,1.535,.058,f(:),x(:),0,-omega,mx0(:),my0(:),mz); 




mz((length(mz)+1)/2)

mxy = mx+1i*my;

figure(1); subplot(413); plot(x,abs(mxy)); 
grid on; xlabel('z position (cm)');

figure(1); subplot(414); plot(x,mz); 
grid on; xlabel('z position (cm)');

% %% Simulate Slice Profile using Bloch Simulation
% % x = (-4:.1:4);          % in cm
% % 
% % f = 0;                  % in Hz
% % dt = pulseduration/samples/1e3;
% % t = (1:length(b1))*dt;  % in usec
% % 
% % Bloch Simulation
% % [mx,my,mz] = bloch(b1(:),g(:),t(:),1,.2,f(:),x(:),0); 
% % 
% % Transverse Magnetization
% % mxy = mx+1i*my;
% % 
% % figure(1); subplot(313); plot(x,abs(mxy)); 
% % grid on; xlabel('z position (cm)');
% % 
% % % Simulate Slice Profile using Small Tip Approximation

%%
% theta = atan(abs(B1)./(omega./26752));
% plot(t,theta);

% 
% T1 = 1227;
% T = (1:1:2000);
% Mz = .8+.2*(1-exp(-T./T1));
% plot(T, Mz)
% 
% 
% %% 
% a = 0.8517;
% T1 = 1224;
% mz_0 = 1;
% for i = 1:10
%     mz_1(i) = mz_0 * a  +(1- mz_0*a)*(1 - exp(-1000/T1));
%     mz_0 = mz_1(i);
% end
% 
% %
% T2star = 28;
% mz_2 = mz_1.*exp(-10/T2star)
% 
% %%
% T2star = 28;
% mz_2 = exp(-10/T2star)

























