clear all;
close all;

addpath('../lib_EPGX/')
addpath('../EPGX-src/')
%addpath('D:\Data\Exvivo_Phantom');
addpath('../BlochSimDemo/');
% Resolved issues from kellman2015.m
% Fitting more accurately and resolved second breath-hold opposite sign
% problem

% using different set of parameters (same as kellman's)
%%

% npulse = 60; % Single-shot 
% A final half-alpha 'restore pulse' to return the magnetization into Mz
% t_delay = TI_array(1); % TI = 650 ms
% flip = 180; % prep flip angle = 180 degree
alpha = 35;
T1 = 1000;
T2 = 45;
T1_fat = 260;
T2_fat = 60;
df = 0; %Hz
f_vec = [-600:10:600];

M0 = [0 0 1]';
FA = alpha;
NTR = 30; % NTR is not the definition (It's "Half" of repetition)
TR = 2.4;
dt = 0.1;
PLOT_EACHSPIN = 0;
PLOT_SS = 0;
[MxyTE, BSIM, plotTEidx, bSSFPcat, Msimf, bSSFPend] = bSSFP_engine(M0, T1, T2, df, TR, FA, NTR, dt, f_vec, PLOT_EACHSPIN, PLOT_SS);

t = [double((1:1:(BSIM.NstepTot+bSSFPcat.NstepTot+bSSFPend.NstepTot)))*dt];
idx = find(f_vec == df);
%idx = find(f_vec == 10);
figure('Position', [100 100 400 400]);
plot(t, Msimf(1,:, idx), t, Msimf(2, :, idx), t, Msimf(3, :, idx));
legend('Mx', 'My', 'Mz'); ylim([-1 1]); grid on;
xlabel('ms'); ylabel('Normalized Magnetization');

% figure();
% plot3(Msimf(1,:, idx), Msimf(2,:, idx), Msimf(3,:, idx), '.-'); axis equal;
% grid on;
% hold on;
% for i = 1:numel(t)
%     hline = plot3([0 Msimf(1,i, idx)], [0 Msimf(2,i, idx)], [0 Msimf(3,i, idx)], 'r');
%     xlabel('Mx'); ylabel('My');zlabel('Mz')
%     xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);
%     pause(.02)
%     if i ~= numel(t)
%         set(hline,'visible','off');
%     end
% end

f_vec_3t = [-600:10:600] + 420;
[MxyTE_fat3t, BSIM_fat3t, plotTEidx_fat3t, bSSFPcat_fat3t, Msimf_fat3t, bSSFPend_fat3t] = bSSFP_engine(M0, T1_fat, T2_fat, df, TR, FA, NTR, dt, f_vec_3t, PLOT_EACHSPIN, PLOT_SS);

f_vec_15t = [-600:10:600] + 210;
[MxyTE_fat15t, BSIM_fat15t, plotTEidx_fat15t, bSSFPcat_fat15t, Msimf_fat15t, bSSFPend_fat15t] = bSSFP_engine(M0, T1_fat, T2_fat, df, TR, FA, NTR, dt, f_vec_15t, PLOT_EACHSPIN, PLOT_SS);
fig = figure('Position', [100 100 400 400]);

%pos_idx = real(MxyTElast) >= 0;
%neg_idx = real(MxyTElast) < 0;

MxyTElast = squeeze( MxyTE(end, :) );
MxyTElast_fat3t = squeeze( MxyTE_fat3t(end, :) );
MxyTElast_fat15t = squeeze( MxyTE_fat15t(end, :) );

pos_idx = imag(MxyTElast) >= 0;
neg_idx = imag(MxyTElast) < 0;
pos_idx_fat3t = imag(MxyTElast_fat3t) >= 0;
neg_idx_fat3t = imag(MxyTElast_fat3t) < 0;
pos_idx_fat15t = imag(MxyTElast_fat15t) >= 0;
neg_idx_fat15t = imag(MxyTElast_fat15t) < 0;
MxyTElast_phaseenc = -abs(MxyTElast .* pos_idx ) + abs(MxyTElast .* neg_idx);
MxyTElast_phaseenc_fat3t = -abs(MxyTElast_fat3t .* pos_idx_fat3t ) + abs(MxyTElast_fat3t .* neg_idx_fat3t);
MxyTElast_phaseenc_fat15t = -abs(MxyTElast_fat15t .* pos_idx_fat15t ) + abs(MxyTElast_fat15t .* neg_idx_fat15t);
% MxyTElast_phaseenc = -abs(MxyTElast .* pos_idx) + abs(MxyTElast .* neg_idx);

idx = find(f_vec == 0);
idx_fat3t = find(f_vec_3t == 420);
idx_fat15t = find(f_vec_15t == 210);
h1 = plot(f_vec, MxyTElast_phaseenc);
hold on;
h2 = plot(f_vec, MxyTElast_phaseenc_fat3t);
%h3 = plot(f_vec, MxyTElast_phaseenc_fat15t);

plot(0, MxyTElast_phaseenc(idx), 'r*', 'linewidth', 1.5, 'MarkerSize', 16);
plot(0, MxyTElast_phaseenc_fat3t(idx_fat3t), 'r*', 'linewidth', 1.5, 'MarkerSize', 16);
%plot(0, MxyTElast_phaseenc_fat15t(idx_fat15t), 'r*');
%title('Balanced SSFP frequency profile: phased-enc');
%xlabel('Hz'); ylabel('Normalized |Mxy|');
grid on;
set([h1], 'linewidth', 2.5)
%set([h2,h3], 'linewidth', 2.5)
set([h2], 'linewidth', 2.5)
%legend({'Myocardium', 'Fat 3T', 'Fat 1.5T'})
set(gca, 'FontSize', 18);

%% T1 mapping MOLLI
TI_array = [102, 182, 935, 1010, 1762, 1840, 2587, 3410];
figure();
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

%% Fig. 2
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
Mz0 = 1;
TD = trigger; % ms

RAMP_DOWN = 1;
npulse = 58 + num_rampup + RAMP_DOWN;
%[t_total, Mz_total_total, t_readout, Mxy_readout] = seq_T1MOLLI_noMT(TI_array, TD, npulse, T1, T2,...
%    alpha, TR, prep, num_rampup, Mz0, restore_pulse, trigger, trigger2);

[t_total, Mz_total_total, t_readout, Mxy_readout] = seq_T1MOLLI_noMT_bloch(TI_array, TD, npulse, T1, T2, alpha, TR, prep, Mz0, trigger, trigger2, df, RAMP_DOWN);
[t_total, Mz_total_total_fat3t, t_readout, Mxy_readout_fat3t] = seq_T1MOLLI_noMT_bloch(TI_array, TD, npulse, T1_fat3t, T2_fat3t, alpha, TR, prep, Mz0, trigger, trigger2, df_fat3t, RAMP_DOWN);
[t_total, Mz_total_total_fat15t, t_readout, Mxy_readout_fat15t] = seq_T1MOLLI_noMT_bloch(TI_array, TD, npulse, T1_fat15t, T2_fat15t, alpha, TR, prep, Mz0, trigger, trigger2, df_fat15t, RAMP_DOWN);

%[tr_total, Mzr_total_total, tr_readout, Mxyr_readout] = seq_T1MOLLI_noMT(TI_array, TD, npulse, T1_bound, T2_bound,...
%    alpha, TR, prep, num_rampup, Mz0, restore_pulse, trigger, trigger2);
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


%% Kellman figure 2
FF = [0:1:8]/8;
figure();
Mxy_readout_array_PhaseEnc = zeros(numel(FF), numel(Mxy_readout));
for i = 1:numel(FF)
    ff = FF(i);
    Mxy_readout_comp = (1-ff) * Mxy_readout + ff * Mxy_readout_fat15t;
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
figure('Position', [100 100 800 800]);
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

%% Plot of myocardium vs fat @ 3T
%% For Fat 3T

figure('Position', [100 100 800 800]);
plot(t_total/1000, Mz_total_total, 'LineWidth', 2.5)
hold on;
%plot(t_total/1000, Mz_total_total_fat3t, '-.', 'LineWidth', 2.5)
plot(t_total/1000, Mz_total_total_fat3t, '-.', 'LineWidth', 2.5)

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

Mxy_readout_array_fat3t = MOLLI_readout_reorder(Mxy_readout_fat3t);
Mxy_readout_array_PhaseEnc_fat3t = Mxy_PhaseEnc(Mxy_readout_array_fat3t);

figure();
plot(TI_array_sorted/1000, Mxy_readout_array_PhaseEnc, 'o-', 'LineWidth', 1.5)
hold on;
plot(TI_array_sorted/1000, Mxy_readout_array_PhaseEnc_fat3t, '*-', 'LineWidth', 1.5)
legend({'Myocardium', 'Fat 1.5T'}, 'Location', 'SouthEast');
xlabel('TI (s)'); ylabel('Signal')
grid on;

%% Kellman figure 2 @ fat 3T
FF = [0:1:10]/10;
figure();
Mxy_readout_array_PhaseEnc = zeros(numel(FF), numel(Mxy_readout));
for i = 1:numel(FF)
    ff = FF(i);
    Mxy_readout_comp = (1-ff) * Mxy_readout + ff * Mxy_readout_fat3t;
    Mxy_readout_array = MOLLI_readout_reorder(Mxy_readout_comp);
    Mxy_readout_array_PhaseEnc(i, :) = Mxy_PhaseEnc(Mxy_readout_array);
    subplot(6,2,i);
    plot(TI_array_sorted/1000, Mxy_readout_array_PhaseEnc(i, :), 'r-', 'LineWidth', 2)
    hold on;
    plot(TI_array_sorted/1000, Mxy_readout_array_PhaseEnc(i, :), 'o', 'LineWidth', 2)
    ylim([-0.4 0.3])
    grid on;
    xlabel('Inversion Time (s)'); ylabel('Mxy');
    title(cat(2, 'FF = ', num2str(ff)));
    set(gca,'fontsize', 18)
    
end

%% Exp fitting (kellman figure 2) @ fat 3T
x=TI_array_sorted';
native_t1_noMT_array = zeros(1, numel(FF));
figure('Position', [100 100 1200 400]);
for i = 1:numel(FF)
y=Mxy_readout_array_PhaseEnc(i, :).';
g = fittype('a-b*exp(-c*x)');
f0 = fit(x,y,g,'StartPoint',[.0;.0; 0.001]);
xx = linspace(1,3500,100);

subplot(2,6,i)
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

%% Exp fitting (kellman figure 2) @ fat 3T
x=TI_array_sorted';
native_t1_noMT_array = zeros(1, numel(FF));
figure('Position', [100 100 1200 400]);
for i = 1:numel(FF)
y=Mxy_readout_array_PhaseEnc(i, :).';
g = fittype('a-b*exp(-c*x)');
f0 = fit(x,y,g,'StartPoint',[.0;.0; 0.001]);
xx = linspace(1,3500,100);

subplot(2,6,i)
plot(x,y,'ro', xx,f0(xx),'b-', 'LineWidth', 2);
grid on;
ylim([-0.4 0.3]);
%xlabel('Inversion Time (s)'); ylabel('Mxy');
%title(cat(2, 'FF = ', num2str(FF(i))));
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
%text(x_text, y_text, txt,'HorizontalAlignment','right', 'FontSize', 16)
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})
end

%% fig 4.
figure('Position', [100 100 400 400]);
plot(FF, native_t1_noMT_array, 'LineWidth', 2); ylim([0 2500]);
grid on;
%hold on;
%plot(FF, native_t1_noMT_array(1) * (1-FF) + native_t1_noMT_array(end) * FF, 'LineWidth', 2)
