clear all;
close all;

addpath('D:\Data\Exvivo_Phantom\lib\')
addpath('D:\Data\Exvivo_Phantom\EPGX-src')
addpath('D:\Data\Exvivo_Phantom');
addpath('D:\src\BlochSimDemo');
addpath('D:\src\M219')
%% T1 mapping MOLLI
TI_array = [102, 182, 935, 1010, 1762, 1840, 2587, 3410];
figure();
b1 = 750;
TR = 2.4;
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
%% Kellman figure 2
FF = [0:1:8]/8;
figure();
Mxy_readout_array_PhaseEnc = zeros(numel(FF), numel(Mxy_readout));
for i = 1:numel(FF)
    ff = FF(i);
    Mxy_readout_comp = (1-ff) * Mxy_readout + ff * Mxy_readout_fat3t;
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

%% OPTIONAL FROM BELOW
% OPTIONAL FROM BELOW
% WRITE MOVIES
%%
Mxy = M_total_total_fat15t(1,:) + 1i * (M_total_total_fat15t(2,:));
Mxy_PhaseEnc_array = Mxy_PhaseEnc(Mxy);
figure();
curve = animatedline('Color', 'k');
curve2 = animatedline('Color', 'r');
set(gca, 'XLim', [0 10], 'YLim', [-1 1]);
%legend({'Mz', 'Mxy'})
grid on;
for i = 1:length(t_total)
    addpoints(curve, t_total(i)/1000, Mz_total_total(i));
    addpoints(curve2, t_total(i)/1000, Mxy_PhaseEnc_array(i));
    drawnow limitrate
    pause(0.002)
end

%% Mxy
Mxy = M_total_total(1,:) + 1i * (M_total_total(2,:));
Mxy_PhaseEnc_array = Mxy_PhaseEnc(Mxy);
figure();
curve = animatedline('Color', 'k');
set(gca, 'XLim', [0 10], 'YLim', [-1 1]);
grid on;
for i = 1:length(t_total)
    addpoints(curve, t_total(i)/1000, Mxy_PhaseEnc_array(i));
    drawnow limitrate
    pause(0.0025)
end

%%
figure();

for i = 1:length(t_total)
    plot(t_total(1:i)/1000, Mz_total_total(1:i), 'b-');
    set(gca, 'XLim', [0 10], 'YLim', [-1 1.05]);
    hold on;
    plot(t_total(1:i)/1000, Mz_total_total_fat15t(1:i), 'r--');
    hold off;
    legend({'Myocardium', 'Fat 1.5T'}, 'Location', 'SouthEast');
    grid on;
    %hold on;
    % drawnow limitrate
    % pause(0.002)
    movieVector(i) = getframe;
end

 myWriter = VideoWriter('curve');
 myWriter.FrameRate = 200;
 
 % Open the videowriter object, write the movie and close the file
 open(myWriter);
 writeVideo(myWriter, movieVector);
 close(myWriter);
 
%%
figure();
set(gca, 'XLim', [0 10], 'YLim', [-1 1]);
grid on;
hold on;
for i = 1:length(t_total)
    plot(t_total(i)/1000, Mz_total_total(i), 'o');
    %hold on;
    % drawnow limitrate
    % pause(0.002)
    movieVector(i) = getframe;
end

 myWriter = VideoWriter('curve');
 myWriter.FrameRate = 200;
 
 % Open the videowriter object, write the movie and close the file
 open(myWriter);
 writeVideo(myWriter, movieVector);
 close(myWriter);