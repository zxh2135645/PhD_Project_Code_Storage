clear all;
close all;
%% Fig. A
addpath('D:\Data\Exvivo_Phantom\lib\')
addpath('D:\Data\Exvivo_Phantom\EPGX-src')
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
%% Fig. B (without MT)
TI_array = [102, 935, 1762, 182, 1010, 1840, 2587, 3410];
npulse = 60 + num_rampup; % Single-shot 
% A final half-alpha 'restore pulse' to return the magnetization into Mz
t_delay = TI_array(1); % TI = 650 ms
flip = 180; % prep flip angle = 180 degree
alpha = 35;
T1_bound = 1000;
T2_bound = 8.5e-3;
T1 = 1133.2;
T2 = 35;

prep = struct;
prep.flip = d2r(flip);
prep.t_delay = t_delay;
Mz0 = 1;
TD = trigger; % ms

[t_total, Mz_total_total, t_readout, Mxy_readout] = seq_T1MOLLI_noMT(TI_array, TD, npulse, T1, T2,...
    alpha, TR, prep, num_rampup, Mz0, restore_pulse, trigger, trigger2);
[tr_total, Mzr_total_total, tr_readout, Mxyr_readout] = seq_T1MOLLI_noMT(TI_array, TD, npulse, T1_bound, T2_bound,...
    alpha, TR, prep, num_rampup, Mz0, restore_pulse, trigger, trigger2);
%% Plot for Fig. B
figure();
plot(t_total/1000, Mz_total_total, 'LineWidth', 2)
hold on;
plot(tr_total/1000, Mzr_total_total, '-.', 'LineWidth', 2)

legend({'Remote Free (no MT)', 'Remote Bound (no MT)'});
xlabel('Time (s)'); ylabel('M_z/M_0')
grid on;
xlim([0 10])
%% Fig. C (with MT)
npulse = 60 + num_rampup;
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
% Assuming pulse duration is 20 ms
trf_prep = 20.00;
alpha_inv = 180;
MT_prep.B1SqrdTau = 2^2 * (d2r(alpha_inv)./(trf_prep.*gam)).^2.*trf_prep; 

M0_remote = [0 0 1-MT_para_remote.F MT_para_remote.F]';

[t_total, Mzmt_total_total, t_readout_mt, Mxy_readout_mt] = seq_T1MOLLI_MT(TI_array, TD, npulse,...
    alpha, TR, MT_para_remote, MT_prep, num_rampup, M0_remote, restore_pulse, trigger, trigger2);

%% Plot for Fig. C
figure();
plot(t_total/1000, Mz_total_total, 'LineWidth', 2)
hold on;
plot(t_total/1000, Mzmt_total_total, '-.', 'LineWidth', 2)

legend({'Free pool no MT', 'Free pool with MT'});
xlabel('Time (s)'); ylabel('M_z/M_0')
grid on;

%% Fig. D
TI_array = sort(TI_array)+ TR * (npulse - num_rampup) / 2; % + TR * (npulse - num_rampup) / 2
Mxy_readout_array_mt = Mxy_readout_mt;
Mxy_readout_array_mt(1) = -Mxy_readout_mt(1);
Mxy_readout_array_mt(2) = -Mxy_readout_mt(4);
Mxy_readout_array_mt(3) = Mxy_readout_mt(2);
Mxy_readout_array_mt(4) = Mxy_readout_mt(5);
Mxy_readout_array_mt(5) = Mxy_readout_mt(3);
Mxy_readout_array = Mxy_readout;
Mxy_readout_array(1) = -Mxy_readout(1);
Mxy_readout_array(2) = -Mxy_readout(4);
Mxy_readout_array(3) = Mxy_readout(2);
Mxy_readout_array(4) = Mxy_readout(5);
Mxy_readout_array(5) = Mxy_readout(3);
figure();
plot(TI_array/1000, Mxy_readout_array, 'o-', 'LineWidth', 1.5)
hold on;
plot(TI_array/1000, Mxy_readout_array_mt, 'x-', 'LineWidth', 1.5)

legend({'Singal no MT', 'Signal with MT'}, 'Location', 'SouthEast');
xlabel('TI (s)'); ylabel('Signal')
grid on;

%% Put it into one Figure
figure('Renderer', 'painters', 'Position', [200 200 1350 900]);
subplot(2,2,1)
line([0,10],[0,0],'Color', 'black')
hold on;
plot([trigger trigger]/1000, [-b1 b1], 'LineWidth', 1.5)


% 3 bSSFP readout
rectangle('Position', [(window-acq_win)/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])
rectangle('Position', [(trigger+TI_array(3))/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])
rectangle('Position', [(trigger+TI_array(5))/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])

% 3 recovery time
% 5 bSSFP readout
plot([trigger2 trigger2]/1000, [-b1 b1], 'LineWidth', 1.5)
rectangle('Position', [(7*window-acq_win)/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])
rectangle('Position', [(trigger2+TI_array(4))/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])
rectangle('Position', [(trigger2+TI_array(6))/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])
rectangle('Position', [(trigger2+TI_array(7))/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])
rectangle('Position', [(trigger2+TI_array(8))/1000 -b1 acq_win/1000 b1*2], 'FaceColor',[0.2 0.2 0.2])
xlabel('Time (s)')

ylabel('B_1 (Hz)')
set(gca,'FontSize',16)

subplot(2,2,2)
plot(t_total/1000, Mz_total_total, 'LineWidth', 2)
hold on;
plot(tr_total/1000, Mzr_total_total, '-.', 'LineWidth', 2)
legend({'Remote Free (no MT)', 'Remote Bound (no MT)'}, 'Location', 'SouthEast');
xlabel('Time (s)'); ylabel('M_z/M_0')
grid on;
xlim([0 10])
set(gca,'FontSize',16)

subplot(2,2,3)
plot(t_total/1000, Mz_total_total, 'LineWidth', 2)
hold on;
plot(t_total/1000, Mzmt_total_total, '-.', 'LineWidth', 2)
legend({'Free pool no MT', 'Free pool with MT'}, 'Location', 'SouthEast');
xlabel('Time (s)'); ylabel('M_z/M_0')
grid on;
set(gca,'FontSize',16)

subplot(2,2,4)
plot(TI_array/1000, Mxy_readout_array, 'o-', 'LineWidth', 1.5)
hold on;
plot(TI_array/1000, Mxy_readout_array_mt, 'x-', 'LineWidth', 1.5)
legend({'Singal no MT', 'Signal with MT'}, 'Location', 'SouthEast');
xlabel('TI (s)'); ylabel('Signal')
grid on;
set(gca,'FontSize',16)

%% Exp fitting
x=TI_array';
y=Mxy_readout_array.';
g = fittype('a-b*exp(-c*x)');
f0 = fit(x,y,g,'StartPoint',[.5;.5; 0.001]);
xx = linspace(1,3500,100);
figure()
plot(x,y,'o',xx,f0(xx),'b-', 'LineWidth', 1.5);

coef = coeffvalues(f0);
native_t1_noMT = 1/coef(3) * (coef(2) / coef(1) - 1)

x=TI_array';
y=Mxy_readout_array_mt.';
g = fittype('a-b*exp(-c*x)');
f0 = fit(x,y,g,'StartPoint',[.5;.5; 0.001]);
xx = linspace(1,3500,100);
hold on;
plot(x,y,'x',xx,f0(xx),'r-', 'LineWidth', 1.5);
grid on;
legend({'Singal no MT','', 'Signal with MT',''}, 'Location', 'SouthEast');coef = coeffvalues(f0);
native_t1_MT = 1/coef(3) * (coef(2) / coef(1) - 1)
xlabel('TI (s)'); ylabel('Signal')

%% Fig. C (with MT)
TI_array = [102, 182, 935, 1010, 1762, 1840, 2587, 3410];
trigger = window-TI_array(1)-acq_win;
trigger2 = 7 * window - TI_array(2) - acq_win;
TI_array = [102, 935, 1762, 182, 1010, 1840, 2587, 3410];

npulse = 60 + num_rampup;
gam = 267.5221 *1e-3; % rad /ms /uT

MT_para_remote = struct;
MT_para_remote.T1x = [T1 T1];
MT_para_remote.T2x = [T2, 8.1e-3];
% MT_para_remote.b1sqrdtau_array = b1sqrdtau_array;
MT_para_remote.F = 0.097;
MT_para_remote.Kf = 5.2e-3;
MT_para_remote.trf = 0.600; % ms

T1_mi = 1384.7;
T2_mi = 34;
MT_para_mi = struct;
MT_para_mi.T1x = [T1_mi T1_mi];
MT_para_mi.T2x = [T2_mi, 8.1e-3];
% MT_para_remote.b1sqrdtau_array = b1sqrdtau_array;
MT_para_mi.F = 0.2;
MT_para_mi.Kf = 20e-3;
MT_para_mi.trf = 0.600; % ms

MT_prep = struct;
MT_prep.flip = d2r(flip);
MT_prep.t_delay = t_delay;
% Assuming pulse duration is 20 ms
trf_prep = 20.00;
alpha_inv = 180;
MT_prep.B1SqrdTau = 2^2 * (d2r(alpha_inv)./(trf_prep.*gam)).^2.*trf_prep; 

M0_remote = [0 0 1-MT_para_remote.F MT_para_remote.F]';
M0_mi = [0 0 1-MT_para_mi.F MT_para_mi.F]';
[t_total, Mzmt_total_total, t_readout_mt, Mxy_readout_mt] = seq_T1MOLLI_MT(TI_array, TD, npulse,...
    alpha, TR, MT_para_remote, MT_prep, num_rampup, M0_remote, restore_pulse, trigger, trigger2);
[t_total, Mzmi_total_total, t_readout_mi, Mxy_readout_mi] = seq_T1MOLLI_MT(TI_array, TD, npulse,...
    alpha, TR, MT_para_mi, MT_prep, num_rampup, M0_mi, restore_pulse, trigger, trigger2);
%% Plot for Fig. C
figure();
plot(t_total/1000, Mzmt_total_total, 'LineWidth', 2)
hold on;
plot(t_total/1000, Mzmi_total_total, '-.', 'LineWidth', 2)

legend({'Remote with MT', 'MI with MT'});
xlabel('Time (s)'); ylabel('M_z/M_0')
grid on;
%% Fig. D
TI_array = sort(TI_array)+ TR * (npulse - num_rampup) / 2; % + TR * (npulse - num_rampup) / 2
Mxy_readout_array_mt = Mxy_readout_mt;
Mxy_readout_array_mt(1) = -Mxy_readout_mt(1);
Mxy_readout_array_mt(2) = -Mxy_readout_mt(4);
Mxy_readout_array_mt(3) = Mxy_readout_mt(2);
Mxy_readout_array_mt(4) = Mxy_readout_mt(5);
Mxy_readout_array_mt(5) = Mxy_readout_mt(3);
Mxy_readout_array = Mxy_readout_mi;
Mxy_readout_array(1) = -Mxy_readout_mi(1);
Mxy_readout_array(2) = -Mxy_readout_mi(4);
Mxy_readout_array(3) = Mxy_readout_mi(2);
Mxy_readout_array(4) = Mxy_readout_mi(5);
Mxy_readout_array(5) = Mxy_readout_mi(3);
figure();
plot(TI_array/1000, Mxy_readout_array_mt, 'x-', 'LineWidth', 1.5)
hold on;
plot(TI_array/1000, Mxy_readout_array, 'o-', 'LineWidth', 1.5)
legend({'MI with MT', 'Remote with MT'}, 'Location', 'SouthEast');
xlabel('TI (s)'); ylabel('Signal')
grid on;
%% Exp fitting
x=TI_array';
y=Mxy_readout_array.';
g = fittype('a-b*exp(-c*x)');
f0 = fit(x,y,g,'StartPoint',[.5;.5; 0.001]);
xx = linspace(1,3500,100);
figure()
plot(x,y,'o',xx,f0(xx),'b-', 'LineWidth', 1.5);

coef = coeffvalues(f0);
native_t1_mi_MT = 1/coef(3) * (coef(2) / coef(1) - 1)

x=TI_array';
y=Mxy_readout_array_mt.';
g = fittype('a-b*exp(-c*x)');
f0 = fit(x,y,g,'StartPoint',[.5;.5; 0.001]);
xx = linspace(1,3500,100);
hold on;
plot(x,y,'x',xx,f0(xx),'r-', 'LineWidth', 1.5);
grid on;
legend({'MI with MT','', 'Remote with MT',''}, 'Location', 'SouthEast');coef = coeffvalues(f0);
native_t1_remote_MT = 1/coef(3) * (coef(2) / coef(1) - 1)
xlabel('TI (s)'); ylabel('Signal')