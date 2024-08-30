close all;
clear all;
%%
dx = 0.04; % in mm
dz = 0.5; % mm
res_array = [0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0, 2.4]; % in mm
res_through_array = [2, 4, 6, 8]; % in mm
transmural_array = [0.025, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6];

MxyTE_remote_echo = T2star_Parameters.MxyTE_remote_echo;
MxyTE_hemo_echo = T2star_Parameters.MxyTE_hemo_echo;
MxyTE_lvblood_echo = T2star_Parameters.MxyTE_lvblood_echo;
snr_array = T2star_Parameters.snr_array;

Nx = 1024;
Ny = 1024;
% noise_map = T2star_Parameters.noise_map(:,:,:,1);
% clear T2star_Parameters

TE_array = [2.55, 5.80, 9.90, 15.56, 21.22]';
A = TE_array .^ [0, 1];

Nz = 48; % Effective slides location


%%
k = 1;
n = 2;

transmural = transmural_array(k);

load(cat(2, 'C:\Users\xz100\Documents\Data\T2star_SimulationPhantom\Ellip_T2starMetrics_Blocked_LinReg_Transmural', num2str(transmural), '_ComprehensiveNoiseLevel_', num2str(n), '.mat'));
t2star_map = t2star_metrics.t2star_map;

%%
figure('Position', [100 100 400 400]); imagesc(t2star_map(:,:,48,1,1));
axis image; axis off;
colormap(brewermap([],'RdBu')); clim([0 100]);