clear all;
close all;
% This script is kept as a demo for kellman2015 fig2
% But its simulation engine is not up-to-date 
%%
addpath('../lib_EPGX/')
addpath('../EPGX-src/')
%addpath('D:\Data\Exvivo_Phantom');

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

%% Fig. B Remote (without MT)
addpath('../BlochSimDemo/');
TI_array = [102, 935, 1762, 182, 1010, 1840, 2587, 3410];
npulse = 60 + num_rampup; % Single-shot 
% A final half-alpha 'restore pulse' to return the magnetization into Mz
t_delay = TI_array(1); % TI = 650 ms
flip = 180; % prep flip angle = 180 degree
alpha = 35;
T1 = 1133.2;
T2 = 35;
%T1 = 260;
%T2 = 60;
df = 0; %Hz

prep = struct;
prep.flip = d2r(flip);
prep.t_delay = t_delay;
Mz0 = 1;
TD = trigger; % ms

%[t_total, Mz_total_total, t_readout, Mxy_readout] = seq_T1MOLLI_noMT(TI_array, TD, npulse, T1, T2,...
%    alpha, TR, prep, num_rampup, Mz0, restore_pulse, trigger, trigger2);

[t_total, Mz_total_total, t_readout, Mxy_readout] = seq_T1MOLLI_noMT_bloch(TI_array, TD, npulse, T1, T2, alpha, TR, prep, Mz0, trigger, trigger2, df);

%[tr_total, Mzr_total_total, tr_readout, Mxyr_readout] = seq_T1MOLLI_noMT(TI_array, TD, npulse, T1_bound, T2_bound,...
%    alpha, TR, prep, num_rampup, Mz0, restore_pulse, trigger, trigger2);
% Plot for Fig. B
figure();
plot(t_total/1000, Mz_total_total, 'LineWidth', 2)
%hold on;
%plot(tr_total/1000, Mzr_total_total, '-.', 'LineWidth', 2)

legend({'Remote (no MT)'});
xlabel('Time (s)'); ylabel('M_z/M_0')
grid on;
xlim([0 10])

TI_array_sorted = sort(TI_array)+ TR * (npulse - num_rampup) / 2; % + TR * (npulse - num_rampup) / 2

abs_Mxy_readout = abs(Mxy_readout);
Mxy_readout_array = abs_Mxy_readout;
Mxy_readout_array(1) = -abs_Mxy_readout(1);
Mxy_readout_array(2) = -abs_Mxy_readout(4);
Mxy_readout_array(3) = abs_Mxy_readout(2);
Mxy_readout_array(4) = abs_Mxy_readout(5);
Mxy_readout_array(5) = abs_Mxy_readout(3);
figure();
plot(TI_array_sorted/1000, Mxy_readout_array, 'o-', 'LineWidth', 1.5)

legend({'Singal no MT'}, 'Location', 'SouthEast');
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
legend({'Remote Free (no MT)'}, 'Location', 'SouthEast');
xlabel('Time (s)'); ylabel('M_z/M_0')
grid on;
xlim([0 10])
set(gca,'FontSize',16)

subplot(2,2,3)
plot(t_total/1000, Mz_total_total, 'LineWidth', 2)
legend({'Free pool no MT'}, 'Location', 'SouthEast');
xlabel('Time (s)'); ylabel('M_z/M_0')
grid on;
set(gca,'FontSize',16)

subplot(2,2,4)
plot(TI_array_sorted/1000, Mxy_readout_array, 'o-', 'LineWidth', 1.5)
legend({'Singal no MT'}, 'Location', 'SouthEast');
xlabel('TI (s)'); ylabel('Signal')
grid on;
set(gca,'FontSize',16)

%% Exp fitting
x=TI_array_sorted';
y=Mxy_readout_array.';
g = fittype('a-b*exp(-c*x)');
f0 = fit(x,y,g,'StartPoint',[.5;.5; 0.001]);
xx = linspace(1,3500,100);
figure()
plot(x,y,'o',xx,f0(xx),'b-', 'LineWidth', 1.5);

coef = coeffvalues(f0);
native_t1_noMT = 1/coef(3) * (coef(2) / coef(1) - 1)

%% 
M0 = [0 0 1]';
FA = alpha;
NTR = npulse - num_rampup;
NTR = 60;
TR = 2.4;
dt = 0.1;
PLOT_EACHSPIN = 0;
PLOT_SS = 0;
T1 = 400;
T2 = 100;
f_vec = [-600:10:600];

[MxyTE, BSIM, plotTEidx, bSSFPcat, Msimf] = bSSFP_engine(M0, T1, T2, df, TR, FA, NTR, dt, f_vec, PLOT_EACHSPIN, PLOT_SS);
% plot signal at TE, over all TRs, for specific df
[Y, fidx1] = min( abs(f_vec-0) ); % locate the pass-band idx
[Y, fidx2] = min( abs(f_vec-1000/2/BSIM.TR) ); % locate the stop-band idx
fidx3 = round( mean([fidx1, fidx2]) );

figure; 
subplot(2,1,1);
plot( 1:numel(plotTEidx),abs(MxyTE(:,fidx1)), 1:numel(plotTEidx),abs(MxyTE(:,fidx3)), 1:numel(plotTEidx),abs(MxyTE(:,fidx2)) );
hold on; line(1+[bSSFPcat.NcatTR bSSFPcat.NcatTR], [0 1], 'Color','k');
xlim([1 numel(plotTEidx)]); legend(num2str(f_vec([fidx1 fidx3 fidx2]).')); ylim([0 1]);
xlabel('TR #'); ylabel('Mxy (normalized)');
title('evolution of pass-band and stop-band over time');

subplot(2,1,2);
plot( 1:numel(plotTEidx),angle(MxyTE(:,fidx1)), 1:numel(plotTEidx),angle(MxyTE(:,fidx3)), 1:numel(plotTEidx),angle(MxyTE(:,fidx2)) );
hold on; line(1+[bSSFPcat.NcatTR bSSFPcat.NcatTR], [-pi pi], 'Color','k');
xlim([1 numel(plotTEidx)]); legend(num2str(f_vec([fidx1 fidx3 fidx2]).')); ylim([-pi pi]);
xlabel('TR #'); ylabel('Mxy (phase)');
title('evolution of pass-band and stop-band over time');

%
% extract Mz
Mz = squeeze( Msimf(3, :, :));
MzTE = Mz(int32(plotTEidx), :);

figure; 
subplot(2,1,1);
plot( 1:numel(plotTEidx),MzTE(:,fidx1), 1:numel(plotTEidx),MzTE(:,fidx3), 1:numel(plotTEidx),MzTE(:,fidx2) );
hold on; line(1+[bSSFPcat.NcatTR bSSFPcat.NcatTR], [-1 1], 'Color','k');
xlim([1 numel(plotTEidx)]); legend(num2str(f_vec([fidx1 fidx3 fidx2]).')); ylim([-1 1]);
xlabel('TR #'); ylabel('Mz (normalized)');
title('evolution of pass-band and stop-band over time');

% frequency profile
