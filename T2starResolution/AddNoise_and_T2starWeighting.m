clear all;
close all;

%% GRE simulation (Remote)
T1_remote = 1125; % Liliana Tribuna et al. 2021
T2_remote = 40.1; % Brianna et al. 2022

T1_hemo = 1125; % Liliana Tribuna et al. 2021
T2_hemo = 9.4; % Brianna et al. 2022

TE_array = [2.56, 5.80, 9.90, 15.56, 21.22];
MxyTE_remote_cell = cell(1,length(TE_array));
MxyTE_hemo_cell = cell(1,length(TE_array));
f_vec = [0];
M0 = [0 0 1]';
FA = 5;
NTR = 192;
TR = 80;
dt = 0.02;
PLOT_EACHSPIN = 0;
PLOT_SS = 0;

for i = 1:length(TE_array)
    TE = TE_array(i);
    [MxyTE, BSIM, plotTEidx, SPGRcat, Msimf] = SPGR_engine(M0, T1_remote, T2_remote, TR, TE, FA, NTR, dt, f_vec, PLOT_EACHSPIN, PLOT_SS);
    MxyTE_remote_cell{i} = MxyTE;
    [MxyTE, BSIM, plotTEidx, SPGRcat, Msimf] = SPGR_engine(M0, T1_hemo, T2_hemo, TR, TE, FA, NTR, dt, f_vec, PLOT_EACHSPIN, PLOT_SS);
    MxyTE_hemo_cell{i} = MxyTE;
end


figure();
MxyTE_remote_echo = zeros(1,length(TE_array));
MxyTE_hemo_echo = zeros(1,length(TE_array));

for i = 1:length(TE_array)
    MxyTE_remote_echo(i) = MxyTE_remote_cell{i}(NTR);
    MxyTE_hemo_echo(i) = MxyTE_hemo_cell{i}(NTR);
end

plot(TE_array, abs(MxyTE_remote_echo), 'LineWidth', 2);
hold on;
plot(TE_array, abs(MxyTE_hemo_echo), 'LineWidth', 2);
xlim([0 30]);


snr_array = 5:5:45;
sigma_array = abs(MxyTE_remote_echo(1)) ./ snr_array;


T2star_Parameters = struct;
T2star_Parameters.MxyTE_remote_echo = MxyTE_remote_echo;
T2star_Parameters.MxyTE_hemo_echo = MxyTE_hemo_echo;
T2star_Parameters.snr_array = snr_array;

Nx = 1024;
Ny = 1024;
Nz = 24;
noise_map = zeros(Nx, Ny, Nz, length(sigma_array));
for n = 1:length(sigma_array)
    noise_map(:,:,:,n) = randn(Nx, Ny, Nz) * sigma_array(n) + 1i * randn(Nx, Ny, Nz) * sigma_array(n);
end

T2star_Parameters.noise_map = noise_map;
save('/Users/jameszhang/Documents/MATLAB/T2star_Resolution_Project/Simulation_Results/3D_MagPurtabation/T2star_Parameters_1024x1024x24.mat', 'T2star_Parameters', '-v7.3');
%% Plot noise map
figure();
noise_map = randn(Nx, Ny) * sigma_array(1) + 1i * randn(Nx, Ny) * sigma_array(1);
imagesc(abs(noise_map)); colormap jet;

%% Fit for T2*
noise_map = randn(1, length(TE_array)) * sigma_array(1) + 1i * randn(1, length(TE_array)) * sigma_array(1);

MxyTE_remote_echo_noised = MxyTE_remote_echo + noise_map;
MxyTE_hemo_echo_noised = MxyTE_hemo_echo + noise_map;

plot(TE_array, abs(MxyTE_remote_echo+noise_map), 'LineWidth', 2);
xlim([0 30]);
f = fit(TE_array',abs(MxyTE_remote_echo_noised)','exp1'); % a*exp(b*x)
T2star = -1/f.b
hold on;
plot(f);

plot(TE_array, abs(MxyTE_hemo_echo+noise_map), 'LineWidth', 2);
xlim([0 30]);
f = fit(TE_array',abs(MxyTE_hemo_echo_noised)','exp1'); % a*exp(b*x)
T2star = -1/f.b
hold on;
plot(f);
