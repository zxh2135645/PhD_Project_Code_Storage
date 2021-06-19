clear all;
close all;
% From seq_block_noMT3
TD = 561.8;
npulse = 64;
T1 = 1000;
T2 = 45;
alpha = 35;
TR = 2.4;
prep = struct;
prep.flip = 3.1416;
prep.t_delay = 102;
M0 = [0;0;1];
df = 0;
RAMP_DOWN = 1;
adiabatic = struct;
adiabatic.mu = 5;
adiabatic.beta1 = 750;
adiabatic.pulseWidth = 20.4; % 20.48 For temporal resolution 0.1 ms
adiabatic.A0 = 0.12;
ddt = 0.1;

% Input for bSSFP_engine
M0 = [0; 0; 1];
T1 = 1000;
T2 = 100;
df = 0;
TR = 2.4;
FA = 35;
NTR = 200;
ddt = 0.1;
f_vec = [-1000:10:1000];
PLOT_EACHSPIN = 0;
PLOT_SS = 1;
RAMP_DOWN = 1;
dt = ddt;
%% From seq_block_MT
TD = 561.8;
npulse = 60;
alpha = 35;
TR = 2.4;
MT_para = struct;
MT_para.T1x = [1133.2, 1133.2];
MT_para.T2x = [35, 0.0081];
MT_para.F = 0.097;
MT_para.Kf = 0.0052;
MT_para.trf = 0.6;

MT_prep = struct;
MT_prep.flip = 3.1416;
MT_prep.t_delay = 102;
MT_prep.B1SqrdTau = 27.581;
M0 = [0;0;0.903;0.097];
num_rampup = 5;
RAMP_DOWN = 1;
ddt = 0.1;


phi = RF_phase_cycle(npulse,'balanced');
fa_flow = d2r(alpha)*ones(npulse,1);
rampup = linspace(0, alpha, num_rampup+1);
rampup = rampup(end-num_rampup+1:end);
fa_flow(1:num_rampup) = d2r(rampup);
t_delay = MT_prep.t_delay;

% Initiation
npulse = npulse + num_rampup + RAMP_DOWN;
dt = ddt;
t_delay = MT_prep.t_delay;
phi = RF_phase_cycle(npulse,'balanced');

% linear ramp-up
fa_flow = d2r(alpha)*ones(npulse,1);
rampup = linspace(0, alpha, num_rampup+1);
rampup = rampup(end-num_rampup+1:end);
fa_flow(1:num_rampup) = d2r(rampup);

% MT Energy deposition
gam = 267.5221 *1e-3; % rad /ms /uT
trf = MT_para.trf; % ms
b1 = d2r(alpha)./(trf.*gam); % uT
b1sqrdtau = b1.^2.*trf; 
b1sqrdtau0 = (d2r(rampup)./(trf.*gam)).^2.*trf;

if RAMP_DOWN == 1
    fa_flow(end) = d2r(alpha/2);
    b1sqrdtau_end = (d2r(alpha/2)./(trf.*gam)).^2.*trf;
    b1sqrdtau_array = [b1sqrdtau0'; b1sqrdtau * ones(npulse-1-num_rampup,1); b1sqrdtau_end];
    % fprintf('%d\n', length(b1sqrdtau_array))
else
    b1sqrdtau_array = [b1sqrdtau0'; b1sqrdtau * ones(npulse-num_rampup,1)];
end

% Super Lorentzian absorption lineshape
[ff,G] = SuperLorentzian(MT_para.T2x(2)*1e-3);
idx = find(ff == min(abs(ff)));
MT_para.G = G(idx);

% For better accuracy
dt = ddt;
total_time = TD + t_delay + TR * npulse ;
one_rep_t = 0:dt:(total_time-dt);
one_rep_Mzmt = ones(1, length(one_rep_t));
one_rep_Mxymt = ones(1, length(one_rep_t));

% 1) Add recovery here
recovery_timepoint = round(TD/dt);
fa_flow0 = zeros(1, recovery_timepoint);
phi0 = zeros(1, recovery_timepoint);
b1sqrdtau_array0 = zeros(1, recovery_timepoint);

[~,fnmt0,Znmt0] = EPGX_GRE_MT(fa_flow0,phi0,b1sqrdtau_array0,...
    dt,MT_para.T1x,MT_para.T2x(1),MT_para.F,MT_para.Kf,MT_para.G,...
    'kmax',inf, 'initMmt', M0);
mxymt0 = size(fnmt0,1)*ifftshift(ifft(ifftshift(fnmt0,1),[],1),1);
mzmt0 = size(Znmt0,1)*ifftshift(ifft(ifftshift(Znmt0,1),[],1),1);

mxymt_array0 = mxymt0(floor(size(mxymt0, 1)/2),:);
if mod(recovery_timepoint, 2) == 0
    mzmt_array0 = real(mzmt0(floor(size(mzmt0, 1)/2+1),:));
else
    mzmt_array0 = real(mzmt0(floor(size(mzmt0, 1)/2),:));
end
    

% MT_prep 
initM1 = [0 0 mzmt_array0(recovery_timepoint) mzmt_array0(end)]';
% Sparsify fa_flow, phi, b1sqrdtau_array and TR
fa_flow_reshape = reshape([fa_flow'; zeros(round(TR/dt)-1, size(fa_flow',2))], [], 1);
phi_reshape = reshape([phi; zeros(round(TR/dt)-1, size(phi,2))], [], 1);
b1sqrdtau_array_reshape = reshape([b1sqrdtau_array'; zeros(round(TR/dt)-1, size(b1sqrdtau_array',2))], [], 1);
% IR-bSSFP simulation
[~,fnmt,Znmt] = EPGX_GRE_MT(fa_flow_reshape,phi_reshape, b1sqrdtau_array_reshape,...
    dt,MT_para.T1x,MT_para.T2x(1),MT_para.F,MT_para.Kf,MT_para.G,...
    'kmax',inf, 'prep', MT_prep, 'initMmt', initM1);
mxymt = size(fnmt,1)*ifftshift(ifft(ifftshift(fnmt,1),[],1),1);
mzmt = size(Znmt,1)*ifftshift(ifft(ifftshift(Znmt,1),[],1),1);

mxymt_array = mxymt(floor(size(mxymt, 1)/2),:);
mzmt_array = real(mzmt(floor(size(mzmt, 1)/2+1),:)); % size is 1x(2*npulse)
% Thus, first half is free pool Mz
% Second half is bound pool Mz
% Need to use mzmt_array(end) as input for next EPG simulation

TI_timepoint = length(one_rep_Mzmt(end-npulse*(TR/dt)-round(t_delay/dt)+1:end-npulse*(TR/dt)));
one_rep_Mzmt_bound = MT_para.F * ones(1, length(one_rep_t));

% 
if TI_timepoint ~= 0
    ti_flow = zeros(TI_timepoint, 1);
    phi_ti = zeros(TI_timepoint, 1);
    b1sqrdtau_array = zeros(TI_timepoint,1);
    MT_prep2 = MT_prep;
    MT_prep2.t_delay = 0;
    initM2 = [0 0 mzmt_array0(recovery_timepoint) mzmt_array0(end)]';
    [~,fnmt_ti,Znmt_ti] = EPGX_GRE_MT(ti_flow,phi_ti,b1sqrdtau_array,...
        dt,MT_para.T1x,MT_para.T2x(1),MT_para.F,MT_para.Kf,MT_para.G,...
        'kmax',inf, 'prep', MT_prep2, 'initMmt', initM2);
    mxymt_ti = size(fnmt_ti,1)*ifftshift(ifft(ifftshift(fnmt_ti,1),[],1),1);
    mzmt_ti = size(Znmt_ti,1)*ifftshift(ifft(ifftshift(Znmt_ti,1),[],1),1);
    
    mxymt_array_ti = mxymt_ti(floor(size(mxymt_ti, 1)/2),:);
    if mod(TI_timepoint, 2) == 0
        mzmt_array_ti = real(mzmt_ti(floor(size(mzmt_ti, 1)/2+1),:));
    else
        mzmt_array_ti = real(mzmt_ti(floor(size(mzmt_ti, 1)/2),:));
    end
    % To sum up
    one_rep_Mzmt(end-npulse*round(TR/dt)-TI_timepoint+1:end-npulse*round(TR/dt)) = mzmt_array_ti(1:TI_timepoint);
    one_rep_Mxymt(end-npulse*round(TR/dt)-TI_timepoint+1:end-npulse*round(TR/dt)) = mxymt_array_ti(1:TI_timepoint);
    % Sum up bound pool
    
    one_rep_Mzmt_bound(end-npulse*round(TR/dt)-TI_timepoint+1:end-npulse*round(TR/dt)) = mzmt_array_ti(TI_timepoint+1:end);
    % one_rep_Mxymt(end-npulse*round(TR/dt)-TI_timepoint+1:end-npulse*round(TR/dt)) = mxymt_array_ti(TI_timepoint+1:end);

end

% Sum up first shot
one_rep_Mzmt(1:recovery_timepoint) = mzmt_array0(1:recovery_timepoint);
one_rep_Mzmt(end-npulse*round(TR/dt)+1:end) = mzmt_array(1:npulse*round(TR/dt));
one_rep_Mxymt(1:recovery_timepoint) = mxymt_array0(1:recovery_timepoint);
one_rep_Mxymt(end-npulse*round(TR/dt)+1:end) = mxymt_array(1:npulse*round(TR/dt));

one_rep_Mzmt_bound(1:recovery_timepoint) = mzmt_array0(1+recovery_timepoint:end);
one_rep_Mzmt_bound(end-npulse*round(TR/dt)+1:end) = mzmt_array(1+npulse*round(TR/dt):end);

% combine two shots
t = one_rep_t;
Mzmt_total = one_rep_Mzmt;
Mzmt_bound_total = one_rep_Mzmt_bound;

TEidx = round(((num_rampup + round((npulse - num_rampup - RAMP_DOWN) / 2)) * TR + TR/2) / dt);
Mxy_readout = mxymt_array(TEidx);
corespond_t = t(recovery_timepoint+TI_timepoint+TEidx); % Readout t

Mxymt_total = one_rep_Mxymt;

%% Dictionary generation starts here
clear all;
%% Fig. A
addpath('../lib_EPGX/')
addpath('../EPGX-src/')
addpath('../T1NFF/')
addpath('../BlochSimDemo/');
addpath('../M219/');
addpath('../MT/');
addpath('../EffectOfFatNIron/');
TI_array = [102, 182, 935, 1010, 1762, 1840, 2587, 3410];
figure();
b1 = 750;
TR = 2.4;
PhaseEnc = 60;
num_rampup = 5;
restore_pulse = 1;
HR = 60*1000/(((TI_array(8) - TI_array(7)) + (TI_array(7) - TI_array(6)) + (TI_array(6) - TI_array(4)) + ...
    (TI_array(4) - TI_array(2)) + (TI_array(5) - TI_array(3)) + (TI_array(3) - TI_array(1))) / 6);
window = round(60*1000 / HR);
acq_win = TR*(PhaseEnc+num_rampup+restore_pulse);
trigger = window-TI_array(1)-acq_win;
line([0,10],[0,0],'Color', 'black')
hold on;
plot([trigger trigger]/1000, [-b1 b1], 'LineWidth', 1.5)

% 3 bSSFP readout
rectangle('Position', [(window-acq_win)/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])
rectangle('Position', [(trigger+TI_array(3))/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])
rectangle('Position', [(trigger+TI_array(5))/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])

% 3 recovery time
% 5 bSSFP readout
trigger2 = 7 * window - TI_array(2) - acq_win;
plot([trigger2 trigger2]/1000, [-b1 b1], 'LineWidth', 1.5)
rectangle('Position', [(7*window-acq_win)/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])
rectangle('Position', [(trigger2+TI_array(4))/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])
rectangle('Position', [(trigger2+TI_array(6))/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])
rectangle('Position', [(trigger2+TI_array(7))/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])
rectangle('Position', [(trigger2+TI_array(8))/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])

xlabel('Time (s)')
ylabel('B_1 (Hz)')

half_readout = zeros(8, 1);
center_k = TR * (num_rampup + PhaseEnc/2) / 1000;
half_readout(1) = (window-acq_win)/1000 + center_k;
half_readout(2) = (trigger+TI_array(3))/1000 + center_k;
half_readout(3) = (trigger+TI_array(5))/1000 + center_k;
half_readout(4) = (7*window-acq_win)/1000 + center_k;
half_readout(5) = (trigger2+TI_array(4))/1000 + center_k;
half_readout(6) = (trigger2+TI_array(6))/1000 + center_k;
half_readout(7) = (trigger2+TI_array(7))/1000 + center_k;
half_readout(8) = (trigger2+TI_array(8))/1000 + center_k;
%% Initiate parameters
num_rampup = 5;
TI_array = [102, 935, 1762, 182, 1010, 1840, 2587, 3410];
npulse = 60; % Single-shot 
% A final half-alpha 'restore pulse' to return the magnetization into Mz
t_delay = TI_array(1); % TI = 650 ms
flip = 180; % prep flip angle = 180 degree
alpha = 35;
T1_bound = 1000;
T2_bound = 8.1e-3;
T1 = 1133.2;
T2 = 45;
TR = 2.4;

prep = struct;
prep.flip = d2r(flip);
prep.t_delay = t_delay;
Mz0 = 1;
TD = trigger; % ms
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
MT_prep.t_delay = t_delay;
%% Dictionary generation needs to run on workstation
% Assuming pulse duration is 20 ms
addpath('../T1NFF/');
trf_prep = 20.00;
alpha_inv = 180;
MT_prep.B1SqrdTau = (d2r(alpha_inv)./(trf_prep.*gam)).^2.*trf_prep; 
%ddt = 0.6;
M0_remote = [0 0 1-MT_para_remote.F MT_para_remote.F]';

%[t_total, Mzmt_total_total, Mxymt_total_total, t_readout_mt, Mxy_readout_mt] = seq_T1MOLLI_MT2(TI_array, TD, npulse,...
%    alpha, TR, MT_para_remote, MT_prep, num_rampup, M0_remote, restore_pulse, trigger, trigger2, ddt);

F_array = 0:0.002:0.1;
Kf_array = 0:0.6:10.2;
%for f = 1:length(F_array)
for f = 1:1
    %for k = 1:length(Kf_array)
    for k = 1:1
        MT_para_remote = struct;
        MT_para_remote.T1x = [T1 T1];
        MT_para_remote.T2x = [T2, T2_bound];
        % MT_para_remote.b1sqrdtau_array = b1sqrdtau_array;
        MT_para_remote.F = F_array(f);
        MT_para_remote.Kf = Kf_array(k);
        MT_para_remote.trf = 0.600; % ms

        MT_prep = struct;
        MT_prep.flip = d2r(flip);
        MT_prep.t_delay = t_delay;
        % Assuming pulse duration is 20 ms
        trf_prep = 20.00;
        alpha_inv = 180;
        MT_prep.B1SqrdTau = (d2r(alpha_inv)./(trf_prep.*gam)).^2.*trf_prep; 
        ddt = 0.4;

        % Looping F and Kf (20x20)
        M0_remote = [0 0 1-MT_para_remote.F MT_para_remote.F]';
        
        disp(['F: ', num2str(F_array(f)), '  Kf: ', num2str(Kf_array(k))]);
        tic;
        [t_total, Mzmt_total_total, Mxymt_total_total, t_readout_mt, Mxy_readout_mt] = seq_T1MOLLI_MT2(TI_array, TD, npulse,...
            alpha, TR, MT_para_remote, MT_prep, num_rampup, M0_remote, restore_pulse, trigger, trigger2, ddt);
        toc;
        
        if f == 1 && k == 1
            Mzmt_dict = zeros(length(F_array), length(Kf_array), length(Mzmt_total_total));
            Mxymt_dict = zeros(length(F_array), length(Kf_array), length(Mxymt_total_total));
        end
        
        Mzmt_dict(f, k, :) = Mzmt_total_total;
        Mxymt_dict(f, k, :) = Mxymt_total_total;
        
    end
end

%% Save as dictionary
Dict.Mzmt_dict = Mzmt_dict;
Dict.Mxymt_dict = Mxymt_dict;
Dict.ddt = ddt;
save_dir = uigetdir;
f_save = 'MT_MOLLI_Dict.mat';
save(cat(2, save_dir, '/', f_save), '-struct', 'Dict');
%% Should create another script for parsing the dictionary
% We can refer to the script below 
% 05/01/2021

% load Dict
% ddt = 0.4;

ddt = 0.2;
f = 1;
k = 1;


TI_array = [102, 182, 935, 1010, 1762, 1840, 2587, 3410];
TR = 2.4;
PhaseEnc = 60;
num_rampup = 5;
restore_pulse = 1;
HR = 60*1000/(((TI_array(8) - TI_array(7)) + (TI_array(7) - TI_array(6)) + (TI_array(6) - TI_array(4)) + ...
    (TI_array(4) - TI_array(2)) + (TI_array(5) - TI_array(3)) + (TI_array(3) - TI_array(1))) / 6);
window = round(60*1000 / HR);
acq_win = TR*(PhaseEnc+num_rampup+restore_pulse);
trigger = window-TI_array(1)-acq_win;
trigger2 = 7 * window - TI_array(2) - acq_win;

TI_array = [102, 935, 1762, 182, 1010, 1840, 2587, 3410];
TD = trigger;
npulse = 60; % Single-shot 
t_delay = TI_array(1); % TI = 650 ms
flip = 180; % prep flip angle = 180 degree
alpha = 35;
% T1_bound = 1000;
T2_bound = 8.1e-3;
T1 = 1133.2;
T2 = 45;
%TR = 2.4;

%prep = struct;
%prep.flip = d2r(flip);
%prep.t_delay = t_delay;
Mz0 = 1;
gam = 267.5221 *1e-3; % rad /ms /uT

MT_para_remote = struct;
MT_para_remote.T1x = [T1 T1];
MT_para_remote.T2x = [T2, T2_bound];
% MT_para_remote.b1sqrdtau_array = b1sqrdtau_array;
MT_para_remote.F = F_array(f);
MT_para_remote.Kf = Kf_array(k);
MT_para_remote.trf = 0.600; % ms

MT_prep = struct;
MT_prep.flip = d2r(flip);
MT_prep.t_delay = t_delay;
% Assuming pulse duration is 20 ms
trf_prep = 20.00;
alpha_inv = 180;
MT_prep.B1SqrdTau = (d2r(alpha_inv)./(trf_prep.*gam)).^2.*trf_prep;

M0_remote = [0 0 1-MT_para_remote.F MT_para_remote.F]';

% t_readout_mt is useful here
[t_total, Mzmt_total_total, Mxymt_total_total, t_readout_mt, Mxy_readout_mt] = seq_T1MOLLI_MT2(TI_array, TD, npulse,...
    alpha, TR, MT_para_remote, MT_prep, num_rampup, M0_remote, restore_pulse, trigger, trigger2, ddt);

%% Try to plot
load('../../T1_Fat_Project/Results/MT_MOLLI_Dict.mat');

F_array = 0:0.002:0.1;
Kf_array = 0:0.6:10.2;

f = 50;
k = 18;
Mzmt_total_total = squeeze(Mzmt_dict(f,k,:));
Mxymt_total_total = squeeze(Mxymt_dict(f,k,:));
%%
figure(); plot(Mzmt_total_total);
hold on;
plot(abs(Mxymt_total_total));
readout_array = round(t_readout_mt/ddt);

for i = 1:length(readout_array)
    xline(readout_array(i));
end

%% To do/ TO fit
addpath('../EffectOfFatNIron/');
%Mzmt_readout = Mzmt_dict(:,:,readout_array);
%Mxymt_readout = Mxymt_dict(:,:,readout_array);
Mzmt_readout = Mzmt_total_total(readout_array);
Mxymt_readout = Mxymt_total_total(readout_array);

%f = 1;
%k = 1;
native_t1_mat = zeros(length(F_array), length(Kf_array));
%for f = 1:length(F_array)
%    for k = 1:length(Kf_array)
        %f = F_array(i);
        %k = Kf_array(j);
        %figure();
        %plot(squeeze(Mzmt_readout));
        Mzmt_readout_reordered = squeeze(MOLLI_readout_reorder(Mzmt_readout(:)));
        Mxymt_readout_reordered = squeeze(MOLLI_readout_reorder(abs(Mxymt_readout(:))));
        %figure();
        %plot(Mzmt_readout_reordered);
        %hold on;
        %plot(Mxymt_readout_reordered);
        
        idx = find(Mxymt_readout_reordered == min(Mxymt_readout_reordered));
        sign_array = ones(1, length(Mxymt_readout_reordered));
        for i = 1:length(Mxymt_readout_reordered)
            if i < idx
                sign_array(i) = -1;
            end
        end
        Mxymt_readout_reordered_sign = Mxymt_readout_reordered .* sign_array';
        TI_array_sorted = sort(TI_array);
        g = fittype('a-b*exp(-c*x)');
        f0 = fit(TI_array_sorted',Mxymt_readout_reordered_sign,g,'StartPoint',[.0;.0; 0.001]);
        coef = coeffvalues(f0);
        % native_t1_mat(f,k) = 1/coef(3) * (coef(2) / coef(1) - 1);
        native_t1 = 1/coef(3) * (coef(2) / coef(1) - 1)
        
        xx = linspace(1,3500,100);
        figure();
        plot(TI_array_sorted',Mxymt_readout_reordered_sign,'ro',xx,f0(xx),'b-', 'LineWidth', 1.5);
 %   end
%end
%% 
figure();
imagesc(native_t1_mat(:,2:end));
%% This part isn't working for Linux
%% Because bloch is compiled for OS/Windows
% A final half-alpha 'restore pulse' to return the magnetization into Mz
T1_fat3t = 400;
T2_fat3t = 100;
T1_fat15t = 260;
T2_fat15t = 60;
df = 0; %Hz
df_fat3t = 420;
df_fat15t = 210;

M0 = [0 0 1]';
RAMP_DOWN = 1;
npulse = 60 + num_rampup + RAMP_DOWN;

%%% User inputs for adiabatic pulse:
adiabatic.mu = 5;   % Phase modulation parameter [dimensionless]
adiabatic.beta1 = 750;   % Frequency modulation parameter [rad/s]
% adiabatic.pulseWidth = 10*2;   % For temporal resolution 0.1ms
adiabatic.pulseWidth = trf_prep;  % RF pulse duration [ms] % According to siemens 3T
adiabatic.A0 = 0.12; 

[t_total, M_total_total_fat3t, t_readout, Mxy_readout_fat3t] = seq_T1MOLLI_noMT_bloch3(TI_array, TD, npulse, T1_fat3t, T2_fat3t, alpha, TR, prep, M0, trigger, trigger2, df_fat3t, adiabatic, RAMP_DOWN, ddt);
[t_total, M_total_total, t_readout, Mxy_readout] = seq_T1MOLLI_noMT_bloch3(TI_array, TD, npulse, T1, T2, alpha, TR, prep, M0, trigger, trigger2, df, adiabatic, RAMP_DOWN, ddt);
%% Plot
figure();
plot(t_total/1000, M_total_total(3, :), 'LineWidth', 2);
hold on;
plot(t_total/1000, Mzmt_total_total', 'LineWidth', 2);
plot(t_total/1000, M_total_total_fat3t(3, :), '-.', 'LineWidth', 2);

xlabel('Time (s)'); ylabel('M_z/M_0');
grid on;
xlim([0 10]);
plot([half_readout half_readout], [-1 1], 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1)
hold off;
set(gca,'fontsize', 18)
legend({'Myo noMT', 'Myo MT', 'Fat'}, 'Location', 'SouthEast');
xlabel('TI (s)'); ylabel('Signal')


TI_array_sorted = sort(TI_array) + TR * ((npulse-num_rampup-RAMP_DOWN) / 2 + num_rampup);

Mxy_readout_array = MOLLI_readout_reorder(Mxy_readout);
Mxy_readout_array_PhaseEnc = Mxy_PhaseEnc(Mxy_readout_array);

Mxy_readout_array_fat3t = MOLLI_readout_reorder(Mxy_readout_fat3t);
Mxy_readout_array_PhaseEnc_fat3t = Mxy_PhaseEnc(Mxy_readout_array_fat3t);

Mxy_readout_array_mt = MOLLI_readout_reorder(Mxy_readout_mt*1i);
Mxy_readout_array_PhaseEnc_mt = Mxy_PhaseEnc(Mxy_readout_array_mt);

figure();
plot(TI_array_sorted/1000, Mxy_readout_array_PhaseEnc, 'o-', 'LineWidth', 1.5)
hold on;
plot(TI_array_sorted/1000, Mxy_readout_array_PhaseEnc_mt, 'o-', 'LineWidth', 1.5)
plot(TI_array_sorted/1000, Mxy_readout_array_PhaseEnc_fat3t, '*-', 'LineWidth', 1.5)
legend({'Myocardium', 'Myo MT', 'Fat 3T'}, 'Location', 'SouthEast');
xlabel('TI (s)'); ylabel('Signal')
grid on;

%% Kellman figure 2 (change between 1.5T and 3T here)
FF = [0:1:8]/8;
Mxy_readout_mt_new = - Mxy_readout_mt;
figure();
Mxy_readout_array_PhaseEnc = zeros(numel(FF), numel(Mxy_readout_mt_new));
for i = 1:numel(FF)
    ff = FF(i);
    Mxy_readout_comp = (1-ff) * Mxy_readout_mt_new + ff * Mxy_readout_fat3t;
    Mxy_readout_array = MOLLI_readout_reorder(Mxy_readout_comp);
    Mxy_readout_array_PhaseEnc(i, :) = Mxy_PhaseEnc(Mxy_readout_array);
    subplot(3,3,i);
    plot(TI_array_sorted/1000, Mxy_readout_array_PhaseEnc(i, :), 'r-', 'LineWidth', 2)
    hold on;
    plot(TI_array_sorted/1000, Mxy_readout_array_PhaseEnc(i, :), 'o', 'LineWidth', 2)
    ylim([-0.4 0.3])
    grid on;
    xlabel('Inversion Time (s)'); ylabel('Mxy');
    title(cat(2, 'FF = ', num2str(ff)));
    set(gca,'fontsize', 18) 
end

%% Exp fitting (kellman figure 2)
x=TI_array_sorted';
native_t1_noMT_array = zeros(1, numel(FF));
figure();
for i = 1:numel(FF)
y=Mxy_readout_array_PhaseEnc(i, :).';
g = fittype('a-b*exp(-c*x)');
f0 = fit(x,y,g,'StartPoint',[.0;.0; 0.001]);
xx = linspace(1,3500,100);

subplot(3,3,i)
plot(x,y,'ro',xx,f0(xx),'b-', 'LineWidth', 1.5);
grid on;
ylim([-0.4 0.3]);
xlabel('Inversion Time (s)'); ylabel('Mxy');
title(cat(2, 'FF = ', num2str(FF(i))));
set(gca,'fontsize', 18)
coef = coeffvalues(f0);
native_t1_noMT_array(i) = 1/coef(3) * (coef(2) / coef(1) - 1);

txt = cat(2, 'T1 = ', num2str(round(native_t1_noMT_array(i))), ' ms');
if i == 4 || i == 5
    x_text = 3000;
    y_text = 0.1;
else
    x_text = 3000;
    y_text = 0;
end
text(x_text, y_text, txt,'HorizontalAlignment','right', 'FontSize', 16)
end