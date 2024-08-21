close all;
clear all;

% Simulation of T1 MOLLI with different iron concentration (off-resonance effect)
addpath('../lib_EPGX/');
addpath('../EPGX-src/');
addpath('../BlochSimDemo/');
addpath('../M219/');

%% T1 mapping MOLLI
TI_array = [102, 182, 935, 1010, 1762, 1840, 2587, 3410];
b1 = 750;
TR = 2.65;
PhaseEnc = 58;
num_rampup = 5;

[trigger, trigger2] = MOLLI_SchemeFunc(TI_array, b1, TR, PhaseEnc, num_rampup);

% main body simulation
% Parameter initialization
% A final half-alpha 'restore pulse' to return the magnetization into Mz
TI_array = [102, 935, 1762, 182, 1010, 1840, 2587, 3410];
t_delay = TI_array(1); % TI = 650 ms
flip = 180; % prep flip angle = 180 degree
alpha = 35;
T1 = 1000;
T2 = 45;
T1_fat3t = 400;
T2_fat3t = 100;
T1_fat15t = 260;
T2_fat15t = 60;

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

df_array = 0:10:500; %Hz
M_total_total_cell =  cell(length(df_array), 1);
Mxy_readout_cell = cell(length(df_array), 1);

for i = 1:numel(df_array)
    df = df_array(i);
    [t_total, M_total_total, t_readout, Mxy_readout] = seq_T1MOLLI_noMT_bloch2(TI_array, TD, npulse, T1, T2, alpha, TR, prep, M0, trigger, trigger2, df, adiabatic, RAMP_DOWN);
    M_total_total_cell{i} = M_total_total;
    Mxy_readout_cell{i} = Mxy_readout;
end
%%
num = 41;
M_total_total = M_total_total_cell{num};
Mxy_readout = Mxy_readout_cell{num};

Mz_total_total = M_total_total(3, :);
% Plot for Fig. B

figure();
plot(t_total/1000, Mz_total_total, 'LineWidth', 2.5)
hold on;

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

figure();
plot(TI_array_sorted/1000, Mxy_readout_array_PhaseEnc, 'o-', 'LineWidth', 1.5)
hold on;
legend({'Myocardium'}, 'Location', 'SouthEast');
xlabel('TI (s)'); ylabel('Signal')
grid on;

%% Fitting
Mxy_readout_array_PhaseEnc = zeros(numel(df_array), numel(Mxy_readout_cell{1}));
x=TI_array_sorted';
native_t1_noMT_array = zeros(1, numel(df_array));

for i = 1:numel(df_array)
    
    Mxy_readout_array = MOLLI_readout_reorder(Mxy_readout_cell{i});
    Mxy_readout_array_PhaseEnc(i, :) = Mxy_PhaseEnc(Mxy_readout_array);
    
    % Fitting intialization
    y=Mxy_readout_array_PhaseEnc(i, :).';
    g = fittype('a-b*exp(-c*x)');
    f0 = fit(x,y,g,'StartPoint',[.0;.0; 0.001]);
    
    coef = coeffvalues(f0);
    native_t1_noMT_array(i) = 1/coef(3) * (coef(2) / coef(1) - 1);
end

%%
figure();
plot(df_array, native_t1_noMT_array, '-', 'LineWidth', 2)
grid on;
xlabel('Frequency Difference'); ylabel('Estiamted T1 (ms)');
title('DF vs Estimated T1');
set(gca,'fontsize', 18); ylim([0 2000])
