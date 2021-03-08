clear all;
close all;
%% IR-bSSFP simulation
addpath(genpath('lib'));
addpath(genpath('EPGX-src'));
num_rampup = 10;
PhaseEnc = 60;
npulse = PhaseEnc + num_rampup; % Single-shot 
% A final half-alpha 'restore pulse' to return the magnetization into Mz
t_delay = 100:50:1200; % TI = 650 ms
flip = 180; % prep flip angle = 180 degree
alpha_array = 15:5:90;
T1_bound = 1000;
T2_bound = 8.5e-3;
T1 = 1133.2;
T2 = 35;
TR = 3.0;
restore_pulse = 0;

prep = struct;
prep.flip = d2r(flip);

Mz0 = 1;
HR = 60;
window = round(60*1000 / HR);
acq_win = TR*(PhaseEnc+num_rampup+restore_pulse);

gam = 267.5221 *1e-3; % rad /ms /uT
MT_para_remote = struct;
MT_para_remote.T1x = [T1 T1];
MT_para_remote.T2x = [T2, 8.1e-3];
% MT_para_remote.b1sqrdtau_array = b1sqrdtau_array;
MT_para_remote.F = 0.097;
MT_para_remote.Kf = 5.2e-3;
MT_para_remote.trf = 0.600; % ms

MT_prep = struct;
MT_prep.flip = d2r(flip);

% Assuming pulse duration is 20 ms
trf_prep = 20.00;
alpha_inv = 180;
M0_remote = [0 0 1-MT_para_remote.F MT_para_remote.F]';

T1_mi = 1384.7;
T2_mi = 34;
MT_para_mi = struct;
MT_para_mi.T1x = [T1_mi T1_mi];
MT_para_mi.T2x = [T2_mi, 8.1e-3];
% MT_para_remote.b1sqrdtau_array = b1sqrdtau_array;
MT_para_mi.F = 0.02;
MT_para_mi.Kf = 2e-3;
MT_para_mi.trf = 0.600; % ms
M0_mi = [0 0 1-MT_para_mi.F MT_para_mi.F]';

Mxy = zeros(length(t_delay), length(alpha_array));
Mxymt = zeros(length(t_delay), length(alpha_array));
Mxymi = zeros(length(t_delay), length(alpha_array));

for i = 1:length(t_delay)
    prep.t_delay = t_delay(i);
    trigger = 500;
    TD = trigger; % ms
    MT_prep.t_delay = t_delay(i);
    MT_prep.B1SqrdTau = 2^2 * (d2r(alpha_inv)./(trf_prep.*gam)).^2.*trf_prep;
    for j = 1:length(alpha_array)
        alpha = alpha_array(j);
        [t, Mz_total, mxys_array, mzs_array, Mxy_readout, t_readout] = ...
            seq_block_noMT(TD, npulse, T1, T2, alpha, TR, prep, num_rampup, Mz0, restore_pulse);
        [~, Mzmt_total, mxysmt_array, mzsmt_array, Mzmt_bound_total, Mxymt_readout, tmt_readout] = ...
            seq_block_MT(TD, npulse, alpha, TR, MT_para_remote, MT_prep, M0_remote, num_rampup, restore_pulse);
        [~, Mzmi_total, mxysmi_array, mzsmi_array, Mzmi_bound_total, Mxymi_readout, tmi_readout] = ...
            seq_block_MT(TD, npulse, alpha, TR, MT_para_mi, MT_prep, M0_mi, num_rampup, restore_pulse);
        Mxy(i,j) = Mxy_readout;
        Mxymt(i,j) = Mxymt_readout;
        Mxymi(i,j) = Mxymi_readout;
    end
end
%% 
figure();
imagesc(Mxy);xlabel('FA (Degree)'); ylabel('TI (ms)')
figure();
imagesc(Mxymt);xlabel('FA (Degree)'); ylabel('TI (ms)')
figure();
imagesc(Mxymi);xlabel('FA (Degree)'); ylabel('TI (ms)')
figure();
imagesc(Mxymi-Mxymt);xlabel('FA (Degree)'); ylabel('TI (ms)')
yticklabels({'150', '250', '350', '450', '550', '650', '750', '850', '950', '1050', '1150'})
xticklabels({'20', '30', '40', '50', '60', '70', '80', '90'})
Mxy_diff = Mxymi-Mxymt;
max(Mxy_diff(:))
[x,y] = find(max(Mxy_diff(:)) == Mxy_diff);

