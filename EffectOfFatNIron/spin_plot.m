clear all;
close all;

addpath('D:\Data\Exvivo_Phantom\lib\')
addpath('D:\Data\Exvivo_Phantom\EPGX-src')
addpath('D:\Data\Exvivo_Phantom');
addpath('D:\src\BlochSimDemo');
%%

% npulse = 60; % Single-shot 
% A final half-alpha 'restore pulse' to return the magnetization into Mz
% t_delay = TI_array(1); % TI = 650 ms
% flip = 180; % prep flip angle = 180 degree
alpha = 35;
T1 = 1000;
T2 = 45;
T1_fat = 260;
T2_fat = 60;
df = 0; %Hz
f_vec = [-600:10:600];

M0 = [0 0 1]';
FA = alpha;
NTR = 100;
TR = 2.65;
dt = 0.1;
PLOT_EACHSPIN = 1;
PLOT_SS = 0;
[MxyTE, BSIM, plotTEidx, bSSFPcat, Msimf] = bSSFP_engine(M0, T1, T2, df, TR, FA, NTR, dt, f_vec, PLOT_EACHSPIN, PLOT_SS);
%% 
t = [double((1:1:(BSIM.NstepTot+bSSFPcat.NstepTot)))*dt];
idx = find(f_vec == df);
idx = find(f_vec == 10);
figure();
plot(t, Msimf(1,:, idx), t, Msimf(2, :, idx), t, Msimf(3, :, idx));
legend('Mx', 'My', 'Mz'); ylim([-1 1]); grid on;
xlabel('ms'); ylabel('Normalized Magnetization');

figure();
plot3(Msimf(1,:, idx), Msimf(2,:, idx), Msimf(3,:, idx), '.-'); axis equal;
grid on;
hold on;
for i = 1:numel(t)
    hline = plot3([0 Msimf(1,i, idx)], [0 Msimf(2,i, idx)], [0 Msimf(3,i, idx)], 'r');
    xlabel('Mx'); ylabel('My');zlabel('Mz')
    xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);
    pause(.02)
    set(hline,'visible','off');
end