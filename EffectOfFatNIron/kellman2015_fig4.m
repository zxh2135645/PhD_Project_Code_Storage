close all;
clear all;

% 1. To plot fat fraction vs estimated T1 
% 2. Save simulated results for estimated T1 mapping - ff_map_t1

addpath('../lib_EPGX/');
addpath('../EPGX-src/');
addpath('../BlochSimDemo/');
addpath('../M219/');

%% T1 mapping MOLLI
TI_array = [102, 182, 935, 1010, 1762, 1840, 2587, 3410];
figure('Position', [100 100 400 400]);
b1 = 750;
TR = 2.65;
PhaseEnc = 58;
num_rampup = 5;

HR = 60*1000/(((TI_array(8) - TI_array(7)) + (TI_array(7) - TI_array(6)) + (TI_array(6) - TI_array(4)) + ...
    (TI_array(4) - TI_array(2)) + (TI_array(5) - TI_array(3)) + (TI_array(3) - TI_array(1))) / 6);
window = round(60*1000 / HR);
acq_win = TR*(PhaseEnc+num_rampup);
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
%% main body simulation
TI_array = [102, 935, 1762, 182, 1010, 1840, 2587, 3410];
num_rampup = 5;
% npulse = 58 + num_rampup; % Single-shot 
% A final half-alpha 'restore pulse' to return the magnetization into Mz
t_delay = TI_array(1); % TI = 650 ms
flip = 180; % prep flip angle = 180 degree
alpha = 35;
T1 = 1000;
T2 = 45;
T1_fat3t = 400;
T2_fat3t = 100;
T1_fat15t = 260;
T2_fat15t = 60;
df = 0; %Hz
df_fat3t = 420;
df_fat15t = 210;


prep = struct;
prep.flip = d2r(flip);
prep.t_delay = t_delay;
M0 = [0 0 1]';
TD = trigger; % ms

RAMP_DOWN = 1;
npulse = 58 + num_rampup + RAMP_DOWN;

%%% User inputs for adiabatic pulse:
adiabatic.mu = 5;   % Phase modulation parameter [dimensionless]
adiabatic.beta1 = 750;   % Frequency modulation parameter [rad/s]
adiabatic.pulseWidth = 10.24*2;   % RF pulse duration [ms] % According to siemens 3T
adiabatic.A0 = 0.12; 

[t_total, M_total_total, t_readout, Mxy_readout] = seq_T1MOLLI_noMT_bloch2(TI_array, TD, npulse, T1, T2, alpha, TR, prep, M0, trigger, trigger2, df, adiabatic, RAMP_DOWN);
[t_total, M_total_total_fat3t, t_readout, Mxy_readout_fat3t] = seq_T1MOLLI_noMT_bloch2(TI_array, TD, npulse, T1_fat3t, T2_fat3t, alpha, TR, prep, M0, trigger, trigger2, df_fat3t, adiabatic, RAMP_DOWN);
[t_total, M_total_total_fat15t, t_readout, Mxy_readout_fat15t] = seq_T1MOLLI_noMT_bloch2(TI_array, TD, npulse, T1_fat15t, T2_fat15t, alpha, TR, prep, M0, trigger, trigger2, df_fat15t, adiabatic, RAMP_DOWN);

%%
Mz_total_total = M_total_total(3, :);
Mz_total_total_fat3t = M_total_total_fat3t(3, :);
Mz_total_total_fat15t = M_total_total_fat15t(3, :);
% Plot for Fig. B
figure();
plot(t_total/1000, Mz_total_total, 'LineWidth', 2.5)
hold on;
%plot(t_total/1000, Mz_total_total_fat3t, '-.', 'LineWidth', 2.5)
plot(t_total/1000, Mz_total_total_fat15t, '-.', 'LineWidth', 2.5)

% legend({'Remote (no MT)'});
xlabel('Time (s)'); ylabel('M_z/M_0');
grid on;
xlim([0 10]);

TI_array_sorted = sort(TI_array) + TR * ((npulse-num_rampup-RAMP_DOWN) / 2 + num_rampup);
%TI_array_sorted/1000 + center_k
hold on;
plot([half_readout half_readout], [-1 1], 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1)
hold off;
set(gca,'fontsize', 18)

Mxy_readout_array = MOLLI_readout_reorder(Mxy_readout);
Mxy_readout_array_PhaseEnc = Mxy_PhaseEnc(Mxy_readout_array);

Mxy_readout_array_fat15t = MOLLI_readout_reorder(Mxy_readout_fat15t);
Mxy_readout_array_PhaseEnc_fat15t = Mxy_PhaseEnc(Mxy_readout_array_fat15t);

figure();
plot(TI_array_sorted/1000, Mxy_readout_array_PhaseEnc, 'o-', 'LineWidth', 1.5)
hold on;
plot(TI_array_sorted/1000, Mxy_readout_array_PhaseEnc_fat15t, '*-', 'LineWidth', 1.5)
legend({'Myocardium', 'Fat 1.5T'}, 'Location', 'SouthEast');
xlabel('TI (s)'); ylabel('Signal')
grid on;

%% Kellman figure 4 (computing)
FF = 0:0.01:1;
Mxy_readout_array_PhaseEnc = zeros(numel(FF), numel(Mxy_readout));
x=TI_array_sorted';
native_t1_noMT_array = zeros(1, numel(FF));

for i = 1:numel(FF)
    ff = FF(i);
    Mxy_readout_comp = (1-ff) * Mxy_readout + ff * Mxy_readout_fat3t;
    Mxy_readout_array = MOLLI_readout_reorder(Mxy_readout_comp);
    Mxy_readout_array_PhaseEnc(i, :) = Mxy_PhaseEnc(Mxy_readout_array);
    
    % Fitting intialization
    y=Mxy_readout_array_PhaseEnc(i, :).';
    g = fittype('a-b*exp(-c*x)');
    f0 = fit(x,y,g,'StartPoint',[.0;.0; 0.001]);
    
    coef = coeffvalues(f0);
    native_t1_noMT_array(i) = 1/coef(3) * (coef(2) / coef(1) - 1);
end

% Plot figure 4
figure();
subplot(1,2,1);
plot(FF*100, native_t1_noMT_array, '-', 'LineWidth', 2)
grid on;
xlabel('Fat Fraction (%)'); ylabel('Estiamted T1 (ms)');
title('FF vs Estimated T1');
set(gca,'fontsize', 18)
subplot(1,2,2);
plot(FF*100, native_t1_noMT_array, '-', 'LineWidth', 2)
ylim([0 2500])
grid on;
xlabel('Fat Fraction (%)'); ylabel('Estiamted T1 (ms)');
title('FF vs Estimated T1');
set(gca,'fontsize', 18)

%% Save simulated results
addpath('../function/')
ff_map_t1 = struct;
ff_map_t1.FF = FF;
ff_map_t1.native_t1_noMT_array = native_t1_noMT_array;

save_dir_L1 = GetFullPath(cat(2, pwd, '/../../T1_Fat_Project/'));
if ~exist(save_dir_L1, 'dir')
   mkdir(save_dir_L1); 
end

save_dir_L2 = GetFullPath(cat(2, save_dir_L1, 'simulation/'));
if ~exist(save_dir_L2, 'dir')
   mkdir(save_dir_L2); 
end
save(cat(2, save_dir_L2, 'kellman2015_fig4.mat'), 'ff_map_t1');