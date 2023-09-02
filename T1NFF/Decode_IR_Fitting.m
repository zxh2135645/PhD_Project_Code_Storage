% This script is trying to decode what is Inversion recovery algorithm from qMRLab 
% doing. Getting information about how to fit the data when we only know
% the absolute intensity.

% Previous node: T1Fitting.m
%%%%%%%%%%%%%%%%%%%%%%%%%% input  file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ir_weighted_metrics.mat
%%%%%%%%%%%%%%%%%%%%%%%%%% output file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

addpath('../function/');

labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT', 'IR'};
label = labels{11};


proj_dir = GetFullPath(cat(2, pwd, '/../../T1_Fat_Project/'));
if ~exist(proj_dir, 'dir')
    mkdir(proj_dir)
end

data_dir = GetFullPath(cat(2, proj_dir, 'data/'));
if ~exist(data_dir, 'dir')
    mkdir(data_dir)
end

subject_name = input('Please type subject name here:  ', 's');
subject_data_dir = GetFullPath(cat(2, data_dir, subject_name, '/'));
if ~exist(subject_data_dir, 'dir')
    mkdir(subject_data_dir)
end

%% Load IR weighted metrics
load(cat(2, subject_data_dir, 'ir_weighted_metrics.mat'));

% RD-NLS-PR (Reduced-Dimension Non-Linear Least Squares with Polarity Restoration)
mean_t1_irse_w = ir_weighted_metrics.mean_t1_irse_w;
mean_t1_irsewe_w = ir_weighted_metrics.mean_t1_irsewe_w;
mean_t1_molli_w = ir_weighted_metrics.mean_t1_molli_w;

T_IR = ir_weighted_metrics.IR_array;
T_IR_molli = ir_weighted_metrics.IR_array_molli;
tp = length(T_IR);
tp_molli = length(T_IR_molli);
method = 'Magnitude';

mean_t1_psirse_w = zeros(size(mean_t1_irse_w));
mean_t1_psmolli_w = zeros(size(mean_t1_molli_w));
mean_t1_psirsewe_w = zeros(size(mean_t1_irse_w));
dim = input('Dimension of vials (1 or 2): ');

if dim == 1
    for i = 1:size(mean_t1_irse_w,1)
        data = mean_t1_irse_w(i,:);
        data_molli = mean_t1_molli_w(i,:);
        data_we = mean_t1_irsewe_w(i,:);
        [T1,b,a,res,idx]=fitT1_IR(data,T_IR,method);
        [T1_molli,b,a,res,idx_molli]=fitT1_IR(data_molli,T_IR_molli,method);
        [T1_we,b,a,res,idx_we]=fitT1_IR(data_we,T_IR,method);
        %     -T1: T1 value
        %     -rb: b parameter
        %     -ra: a prameter
        %     -res: residual of the fit
        %     -idx: index of last polarity restored datapoint (only used for magnitude data)
        sign_array = cat(2, -1 * ones(1,idx), ones(1, tp-idx));
        sign_array_molli = cat(2, -1 * ones(1,idx_molli), ones(1, tp_molli-idx_molli));
        sign_array_we = cat(2, -1 * ones(1,idx_we), ones(1, tp-idx_we));
        mean_t1_psirse_w(i,:) = sign_array .* data;
        mean_t1_psmolli_w(i,:) = sign_array_molli .* data_molli;
        mean_t1_psirsewe_w(i,:) = sign_array_we .* data_we;
    end
elseif dim == 2
    for i = 1:size(mean_t1_irse_w,1)
        for j = 1:size(mean_t1_irse_w, 2)
            data = mean_t1_irse_w(i,j,:);
            data_molli = mean_t1_molli_w(i,j,:);
            data_we = mean_t1_irsewe_w(i,j,:);
            [T1,b,a,res,idx]=fitT1_IR(data,T_IR,method);
            [T1_molli,b,a,res,idx_molli]=fitT1_IR(data_molli,T_IR_molli,method);
            [T1_we,b,a,res,idx_we]=fitT1_IR(data_we,T_IR,method);

            sign_array = cat(2, -1 * ones(1,idx), ones(1, tp-idx));
            sign_array_molli = cat(2, -1 * ones(1,idx_molli), ones(1, tp_molli-idx_molli));
            sign_array_we = cat(2, -1 * ones(1,idx_we), ones(1, tp-idx_we));
            mean_t1_psirse_w(i,j,:) = sign_array .* squeeze(data)';
            mean_t1_psmolli_w(i,j,:) = sign_array_molli .* squeeze(data_molli)';
            mean_t1_psirsewe_w(i,j,:) = sign_array_we .* squeeze(data_we)';
        end
    end
    
end
%% Plot
ff = {'2.5', '5', '10', '20', '30', '40', '100'};
if dim == 1
    figure();
    for i = 1:size(mean_t1_psirse_w, 1)
        subplot(3,3,i);
        plot(T_IR_molli, mean_t1_psmolli_w(i, :)', 'LineWidth', 2);ylim([-1000 1000]);
        hold on;
        yyaxis right;
        plot(T_IR, mean_t1_psirse_w(i, :)', 'LineWidth', 2);ylim([-2000 2000]);
        %plot(T_IR, mean_t1_psirsewe_w(i,:)', 'LineWidth', 2);
        title(cat(2, 'Fat Fraction: ', ff{i}, '%'));
        legend({'MOLLI', 'IRSE'});
    end
elseif dim == 2
    for i = 1:size(mean_t1_psirse_w, 1)
        figure();
        for j = 1:size(mean_t1_psirse_w, 2)
            subplot(3,3,j);
            plot(T_IR_molli, squeeze(mean_t1_psmolli_w(i,j, :)), 'LineWidth', 2);ylim([-1000 1000]);
            hold on;
            yyaxis right;
            plot(T_IR, squeeze(mean_t1_psirse_w(i,j, :)), 'LineWidth', 2);ylim([-2000 4000]);
            plot(T_IR, squeeze(mean_t1_psirsewe_w(i,j,:)), 'LineWidth', 2);
            title(cat(2, 'Fat Fraction: ', ff{j}, '%'));
            legend({'MOLLI', 'IRSE'}, 'Location', 'SouthEast');
        end
    end
end
%% Compare MOLLI weighted image to T1 simulation [5(3)3] (OPTIONAL)
addpath('../EffectOfFatNIron/');
TI_array = T_IR_molli;
b1 = 750;
TR = 2.4;
PhaseEnc = 58;
num_rampup = 5;
[trigger, trigger2, half_readout] = MOLLI533_SchemeFunc(TI_array, b1, TR, PhaseEnc, num_rampup);

TI_array = [T_IR_molli(1), T_IR_molli(3), T_IR_molli(5), T_IR_molli(7), T_IR_molli(8), T_IR_molli(2), T_IR_molli(4), T_IR_molli(6)];
t_delay = TI_array(1); % TI = 650 ms
flip = 180; % prep flip angle = 180 degree
alpha = 35;
T1 = 2543;
T2 = 85;
T1_fat3t = 320;
T2_fat3t = 100;
df = 0; %Hz
df_fat3t = 420;


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

[t_total, M_total_total, t_readout, Mxy_readout] = seq_T1MOLLI_noMT_bloch533(TI_array, TD, npulse, T1, T2, alpha, TR, prep, M0, trigger, trigger2, df, adiabatic, RAMP_DOWN);
[t_total, M_total_total_fat3t, t_readout, Mxy_readout_fat3t] = seq_T1MOLLI_noMT_bloch533(TI_array, TD, npulse, T1_fat3t, T2_fat3t, alpha, TR, prep, M0, trigger, trigger2, df_fat3t, adiabatic, RAMP_DOWN);

%%
Mz_total_total = M_total_total(3, :);
Mz_total_total_fat3t = M_total_total_fat3t(3, :);
% Plot for Fig. B
figure();
plot(t_total/1000, Mz_total_total, 'LineWidth', 2.5)
hold on;
%plot(t_total/1000, Mz_total_total_fat3t, '-.', 'LineWidth', 2.5)
plot(t_total/1000, Mz_total_total_fat3t, '-.', 'LineWidth', 2.5)

% legend({'Remote (no MT)'});
xlabel('Time (s)'); ylabel('M_z/M_0');
grid on;
xlim([0 12]);

TI_array_sorted = sort(TI_array); %TR * ((npulse-num_rampup-RAMP_DOWN) / 2 + num_rampup);
%TI_array_sorted/1000 + center_k
hold on;
plot([half_readout half_readout], [-1 1], 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1)
hold off;
set(gca,'fontsize', 18)

Mxy_readout_array = MOLLI533_readout_reorder(Mxy_readout);
Mxy_readout_array_PhaseEnc = Mxy_PhaseEnc(Mxy_readout_array);

Mxy_readout_array_fat3t = MOLLI533_readout_reorder(Mxy_readout_fat3t);
Mxy_readout_array_PhaseEnc_fat3t = Mxy_PhaseEnc(Mxy_readout_array_fat3t);

figure();
plot(TI_array_sorted/1000, Mxy_readout_array_PhaseEnc, 'o-', 'LineWidth', 1.5)
hold on;
plot(TI_array_sorted/1000, Mxy_readout_array_PhaseEnc_fat3t, '*-', 'LineWidth', 1.5)
legend({'Myocardium', 'Fat 3T'}, 'Location', 'SouthEast');
xlabel('TI (s)'); ylabel('Signal')
grid on;

%% Kellman figure 2 (3T here)
FF = [0, 0.025, 0.05, 0.10, 0.20, 0.30, 0.40, 1];
figure();
Mxy_readout_array_PhaseEnc = zeros(numel(FF), numel(Mxy_readout));
for i = 1:numel(FF)
    ff = FF(i);
    Mxy_readout_comp = (1-ff) * Mxy_readout + ff * Mxy_readout_fat3t;
    Mxy_readout_array = MOLLI533_readout_reorder(Mxy_readout_comp);
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

%% Plot data and simulation together
Mxy_readout_array_PhaseEnc_norm = Mxy_readout_array_PhaseEnc ./ max(abs(Mxy_readout_array_PhaseEnc(:)));
Mxy_readout_array_PhaseEnc_norm(end,:) = -Mxy_readout_array_PhaseEnc_norm(end,:);
mean_t1_psmolli_w_norm = mean_t1_psmolli_w ./ max(mean_t1_molli_w(:));
figure();
for i = 1:numel(FF)
    ff = FF(i);
    subplot(3,3,i);
    h1=plot(TI_array_sorted/1000, Mxy_readout_array_PhaseEnc_norm(i, :), 'r-', 'LineWidth', 2)
    hold on;
    plot(TI_array_sorted/1000, Mxy_readout_array_PhaseEnc_norm(i, :), 'o', 'LineWidth', 2);
    ylim([-1 1]);
    
    h2=plot(TI_array_sorted/1000, mean_t1_psmolli_w_norm(i, :), 'b-', 'LineWidth', 2);
    plot(TI_array_sorted/1000, mean_t1_psmolli_w_norm(i, :), 'o', 'LineWidth', 2);
    grid on;
    legend([h1, h2], {'Simulation', 'Data'}, 'Location', 'SouthEast');
    title(cat(2, 'FF = ', num2str(ff)));
end

%% Exp fitting (kellman figure 2) Old way of fitiing
x=TI_array_sorted';
native_t1_noMT_array = zeros(1, numel(FF));
figure();
for i = 1:numel(FF)
y=Mxy_readout_array_PhaseEnc(i, :).';
g = fittype('a-b*exp(-c*x)');
f0 = fit(x,y,g,'StartPoint',[.0;.0; 0.001]);
xx = linspace(1,5000,100);

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
%% New way of fitting
T1_molli = zeros(1, size(mean_t1_irse_w,1));
T1_sim_molli = zeros(1, size(mean_t1_irse_w,1));

for i = 1:size(mean_t1_irse_w,1)
    data_molli = mean_t1_molli_w(i,:);
    sim_molli = Mxy_readout_array_PhaseEnc(i,:);
    [T1_molli(i),b,a,res,idx]=fitT1_IR(data_molli,T_IR_molli,method);
    T1_molli(i) = T1_molli(i) .* (-b/a - 1);
    [T1_sim_molli(i),b,a,res,idx]=fitT1_IR(sim_molli,T_IR_molli,method);
    T1_sim_molli(i) = T1_sim_molli(i) .* (-b/a - 1);
end

figure();
plot(FF, T1_molli, 'LineWidth', 2); hold on;
plot(FF, T1_sim_molli, 'LineWidth', 2); ylim([0 5000]);