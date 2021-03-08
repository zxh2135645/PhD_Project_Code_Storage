clear all;
close all;
%% IR-bSSFP simulation
addpath(genpath('lib'));
addpath(genpath('EPGX-src'));
num_rampup = 1;
PhaseEnc = 300;
npulse = PhaseEnc + num_rampup; % Single-shot 
% A final half-alpha 'restore pulse' to return the magnetization into Mz
flip = 180; % prep flip angle = 180 degree
alpha = 45;
T1 = 400;
T2 = 70;
TR = 6.86; % On resonance: 2.2, 4.5, 6.8
          % Off resonance: 3.4, 5.7, 8.0
restore_pulse = 0;

prep = struct;
prep.flip = d2r(flip);

Mz0 = 1;

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
MT_para_mi.F = 0.2;
MT_para_mi.Kf = 20e-3;
MT_para_mi.trf = 0.600; % ms
M0_mi = [0 0 1-MT_para_mi.F MT_para_mi.F]';

phi = RF_phase_cycle(npulse,'balanced');

% Initialization
fa_flow = d2r(alpha)*ones(npulse,1);
fa_flow(1) = d2r(alpha/2);
initM = [0 0 Mz0]';

[~,fn0, Zn0] = EPG_GRE(fa_flow,phi,TR,T1,T2,'kmax',inf, 'initM', initM);

mxys = size(fn0,1)*ifftshift(ifft(ifftshift(fn0,1),[],1),1);
mzs = size(Zn0,1)*ifftshift(ifft(ifftshift(Zn0,1),[],1),1);
mxys_array = abs(mxys(floor(size(mxys, 1)/2),:));
if mod(npulse, 2) == 0
    mzs_array = real(mzs(floor(size(mzs, 1)/2+1),:));
else
    mzs_array = real(mzs(floor(size(mzs, 1)/2),:));
end


freq_period1 = linspace(-1/(2*TR*10^-3),1/(2*TR*10^-3),size(mxys,1));
freq_period2 = linspace(1/(2*TR*10^-3),3/(2*TR*10^-3),size(mxys,1));
freq_period3 = linspace(3/(2*TR*10^-3),5/(2*TR*10^-3),size(mxys,1));
freq_period4 = linspace(5/(2*TR*10^-3),7/(2*TR*10^-3),size(mxys,1));
freq_period = [freq_period1 freq_period2 freq_period3 freq_period4];
Mxy_all = repmat(squeeze(abs(mxys(:,end))),4,1);
figure();
plot(freq_period,Mxy_all,'linewidth',2);
grid on;
xlabel('Off-Resonance (Hz)');ylabel('M_{xy}');
hold on;
vline(440,'r','Fat')

[a, b] = min(abs(freq_period - 440));
M_fast = Mxy_all(b);
% TR = 4.4 ms
TR = 4.38;
[~,fn0, Zn0] = EPG_GRE(fa_flow,phi,TR,T1,T2,'kmax',inf, 'initM', initM);

mxys = size(fn0,1)*ifftshift(ifft(ifftshift(fn0,1),[],1),1);
mzs = size(Zn0,1)*ifftshift(ifft(ifftshift(Zn0,1),[],1),1);
mxys_array = abs(mxys(floor(size(mxys, 1)/2),:));
if mod(npulse, 2) == 0
    mzs_array = real(mzs(floor(size(mzs, 1)/2+1),:));
else
    mzs_array = real(mzs(floor(size(mzs, 1)/2),:));
end


freq_period1 = linspace(-1/(2*TR*10^-3),1/(2*TR*10^-3),size(mxys,1));
freq_period2 = linspace(1/(2*TR*10^-3),3/(2*TR*10^-3),size(mxys,1));
freq_period3 = linspace(3/(2*TR*10^-3),5/(2*TR*10^-3),size(mxys,1));
freq_period4 = linspace(5/(2*TR*10^-3),7/(2*TR*10^-3),size(mxys,1));
freq_period = [freq_period1 freq_period2 freq_period3 freq_period4];
Mxy_all = repmat(squeeze(abs(mxys(:,end))),4,1);
figure();
plot(freq_period,Mxy_all,'linewidth',2);
grid on;
xlabel('Off-Resonance (Hz)');ylabel('M_{xy}');
hold on;
vline(440,'r','Fat')

[a, b] = min(abs(freq_period - 440));
M_lowsar = Mxy_all(b);

MTR = (M_lowsar - M_fast) / M_lowsar

%% With MT

% Initialization
TR = 3.2;
fa_flow = d2r(alpha)*ones(npulse,1);
fa_flow(1) = d2r(alpha/2);
initM = [0 0 1-MT_para_remote.F MT_para_remote.F]';
gam = 267.5221 *1e-3; % rad /ms /uT
trf = MT_para_remote.trf; % ms
b1 = d2r(alpha)./(trf.*gam); % uT
b1sqrdtau = 2^2 * b1.^2.*trf;
b1sqrdtau0 = 2^2 * (d2r(alpha/2)./(trf.*gam)).^2.*trf;
b1sqrdtau_array = [b1sqrdtau0'; b1sqrdtau * ones(npulse-1,1)];

[ff,G] = SuperLorentzian(MT_para_remote.T2x(2)*1e-3);
idx = find(ff == min(abs(ff)));
MT_para_remote.G = G(idx);

[~,fnmt0,Znmt0] = EPGX_GRE_MT(fa_flow,phi,b1sqrdtau_array,...
    TR,MT_para_remote.T1x,MT_para_remote.T2x(1),MT_para_remote.F,MT_para_remote.Kf,MT_para_remote.G,...
    'kmax',inf, 'initMmt', initM);

mxymt = size(fnmt0,1)*ifftshift(ifft(ifftshift(fnmt0,1),[],1),1);
mzmt = size(Znmt0,1)*ifftshift(ifft(ifftshift(Znmt0,1),[],1),1);
mxymt_array = abs(mxymt(floor(size(mxymt, 1)/2),:));
if mod(npulse, 2) == 0
    mzmt_array = real(mzmt(floor(size(mzmt, 1)/2+1),:));
else
    mzmt_array = real(mzmt(floor(size(mzmt, 1)/2),:));
end


freq_period1 = linspace(-1/(2*TR*10^-3),1/(2*TR*10^-3),size(mxymt,1));
freq_period2 = linspace(1/(2*TR*10^-3),3/(2*TR*10^-3),size(mxymt,1));
freq_period = [freq_period1 freq_period2];
Mxymt_all = repmat(squeeze(abs(mxymt(:,end))),2,1);
figure();
plot(freq_period,Mxymt_all,'linewidth',2);
grid on;
xlabel('Off-Resonance (Hz)');ylabel('M_{xy}');
hold on;
vline(440,'r','Fat')

[a, b] = min(abs(freq_period - 440));
M_fast = Mxymt_all(b);

% TR = 4.4 ms
TR = 4.4;
fa_flow = d2r(alpha)*ones(npulse,1);
fa_flow(1) = d2r(alpha/2);
initM = [0 0 1-MT_para_remote.F MT_para_remote.F]';
gam = 267.5221 *1e-3; % rad /ms /uT
trf = MT_para_remote.trf; % ms
b1 = d2r(alpha)./(trf.*gam); % uT
b1sqrdtau = 2^2 * b1.^2.*trf;
b1sqrdtau0 = 2^2 * (d2r(alpha/2)./(trf.*gam)).^2.*trf;
b1sqrdtau_array = [b1sqrdtau0'; b1sqrdtau * ones(npulse-1,1)];

[ff,G] = SuperLorentzian(MT_para_remote.T2x(2)*1e-3);
idx = find(ff == min(abs(ff)));
MT_para_remote.G = G(idx);

[~,fnmt0,Znmt0] = EPGX_GRE_MT(fa_flow,phi,b1sqrdtau_array,...
    TR,MT_para_remote.T1x,MT_para_remote.T2x(1),MT_para_remote.F,MT_para_remote.Kf,MT_para_remote.G,...
    'kmax',inf, 'initMmt', initM);

mxymt = size(fnmt0,1)*ifftshift(ifft(ifftshift(fnmt0,1),[],1),1);
mzmt = size(Znmt0,1)*ifftshift(ifft(ifftshift(Znmt0,1),[],1),1);
mxymt_array = abs(mxymt(floor(size(mxymt, 1)/2),:));
if mod(npulse, 2) == 0
    mzmt_array = real(mzmt(floor(size(mzmt, 1)/2+1),:));
else
    mzmt_array = real(mzmt(floor(size(mzmt, 1)/2),:));
end


freq_period1 = linspace(-1/(2*TR*10^-3),1/(2*TR*10^-3),size(mxymt,1));
freq_period2 = linspace(1/(2*TR*10^-3),3/(2*TR*10^-3),size(mxymt,1));
freq_period3 = linspace(3/(2*TR*10^-3),5/(2*TR*10^-3),size(mxymt,1));
freq_period = [freq_period1 freq_period2 freq_period3];
Mxymt_all = repmat(squeeze(abs(mxymt(:,end))),3,1);
figure();
plot(freq_period,Mxymt_all,'linewidth',2);
grid on;
xlabel('Off-Resonance (Hz)');ylabel('M_{xy}');
hold on;
vline(440,'r','Fat')

[a, b] = min(abs(freq_period - 440));
M_lowsar = Mxy_all(b);

MTR = (M_lowsar - M_fast) / M_lowsar