clear all;
close all;

addpath('D:\Data\Exvivo_Phantom\lib\')
addpath('D:\Data\Exvivo_Phantom\EPGX-src')
addpath('D:\Data\Exvivo_Phantom');
addpath('D:\src\BlochSimDemo');
% This code is pretty old?

%%
TI_array = [102, 182, 935, 1010, 1762, 1840, 2587, 3410];
PhaseEnc = 60;
num_rampup = 5;

HR = 60*1000/(((TI_array(8) - TI_array(7)) + (TI_array(7) - TI_array(6)) + (TI_array(6) - TI_array(4)) + ...
    (TI_array(4) - TI_array(2)) + (TI_array(5) - TI_array(3)) + (TI_array(3) - TI_array(1))) / 6);
window = round(60*1000 / HR);

TI_array = [102, 935, 1762, 182, 1010, 1840, 2587, 3410];
npulse = 60 + num_rampup; % Single-shot 
% A final half-alpha 'restore pulse' to return the magnetization into Mz
t_delay = TI_array(1); % TI = 650 ms
flip = 180; % prep flip angle = 180 degree
alpha = 35;
T1 = 1133.2;
T2 = 45;
T1_fat3t = 400;
T2_fat3t = 100;
df = 0; %Hz
df_fat3t = 420;

prep = struct;
prep.flip = d2r(flip);
prep.t_delay = t_delay;
Mz0 = 1;


TR_array = [2:0.1:10];
alpha_array = [10:5:90];
% for i = 1:length(TR_array)
Mxy_readout = zeros(numel(TR_array), numel(alpha_array), numel(TI_array));
Mxy_readout_fat3t = zeros(numel(TR_array), numel(alpha_array), numel(TI_array));
Mxy_readout_array_PhaseEnc = zeros(numel(TR_array), numel(alpha_array), numel(TI_array));
Mxy_readout_array_PhaseEnc_fat3t = zeros(numel(TR_array), numel(alpha_array), numel(TI_array));

for i = 1:length(TR_array)
    TR = TR_array(i);
    acq_win = TR*(PhaseEnc+num_rampup);
    trigger = window-TI_array(1)-acq_win;
    trigger2 = 7 * window - TI_array(2) - acq_win;
    TD = trigger; % ms
    
    for j = 1:length(alpha_array)
        alpha = alpha_array(j);
        [t_total, Mz_total_total, t_readout, Mxy_readout(i,j,:)] = seq_T1MOLLI_noMT_bloch(TI_array, TD, npulse, T1, T2, alpha, TR, prep, Mz0, trigger, trigger2, df);
        [t_total, Mz_total_total_fat3t, t_readout, Mxy_readout_fat3t(i,j,:)] = seq_T1MOLLI_noMT_bloch(TI_array, TD, npulse, T1_fat3t, T2_fat3t, alpha, TR, prep, Mz0, trigger, trigger2, df_fat3t);
        
        Mxy_readout_array = MOLLI_readout_reorder(Mxy_readout(i,j,:));
        Mxy_readout_array_PhaseEnc(i,j,:) = Mxy_PhaseEnc(Mxy_readout_array);
        
        Mxy_readout_array_fat3t = MOLLI_readout_reorder(Mxy_readout_fat3t(i,j,:));
        Mxy_readout_array_PhaseEnc_fat3t(i,j,:) = Mxy_PhaseEnc(Mxy_readout_array_fat3t);
    end
end

%%
figure();
idx = 2;
jdx = 6;
TR = TR_array(idx);
TI_array_sorted = sort(TI_array) + TR * ((npulse-num_rampup) / 2 + num_rampup);
plot(TI_array_sorted/1000, squeeze(Mxy_readout_array_PhaseEnc(idx,jdx, :)), 'o-', 'LineWidth', 1.5)
hold on;
plot(TI_array_sorted/1000, squeeze(Mxy_readout_array_PhaseEnc_fat3t(idx,jdx, :)), '*-', 'LineWidth', 1.5)
legend({'Myocardium', 'Fat 3T'}, 'Location', 'SouthEast');
xlabel('TI (s)'); ylabel('Signal')
grid on;

%% Exp fitting
TI_array_sorted = zeros(numel(TR_array), numel(TI_array));
for i = 1:numel(TR_array)
    TR = TR_array(i);
    TI_array_sorted(i,:) = sort(TI_array) + TR * ((npulse-num_rampup) / 2 + num_rampup);
end


g = fittype('a-b*exp(-c*x)');

native_t1_array = zeros(numel(TR_array), numel(alpha_array));
native_t1_fat_array = zeros(numel(TR_array), numel(alpha_array));

for i = 1:length(TR_array)
    x=squeeze(TI_array_sorted(i,:)).';
    for j = 1:length(alpha_array)
        y=squeeze(Mxy_readout_array_PhaseEnc(i,j,:));
        f0 = fit(x,y,g,'StartPoint',[.5;.5; 0.001]);
        coef = coeffvalues(f0);
        native_t1_array(i, j) = 1/coef(3) * (coef(2) / coef(1) - 1);
        
        y_fat = squeeze(Mxy_readout_array_PhaseEnc_fat3t(i,j, :));
        f0_fat = fit(x,y_fat,g,'StartPoint',[.5;0; 0.001]);
        coef = coeffvalues(f0_fat);
        native_t1_fat_array(i, j) = 1/coef(3) * (coef(2) / coef(1) - 1);
    end
end
%%
figure();
plot(TR_array, native_t1_array(:,6), 'LineWidth', 2);
hold on;
plot(TR_array, native_t1_fat_array(:,6), 'LineWidth', 2);
legend({'Myocardium', 'Fat 3T'}, 'Location', 'SouthEast');
xlabel('TR (ms)'); ylabel('T1')
grid on;

%% myocardium
figure()
[X,Y] = meshgrid(alpha_array, TR_array);
imagesc(alpha_array, TR_array,  native_t1_array);
set(gca,'XTick',10:10:90);
set(gca,'YTick',2:1:10);
xlabel('Flip Angle (degree)'); ylabel('TR (ms)');
hold on;

[C,h] = contour(X,Y,native_t1_array, 'k', 'ShowText','on'); 
clabel(C,h,'FontSize',16)
set(gca, 'FontSize', 18)

%% fat 3t
figure()
[X,Y] = meshgrid(alpha_array, TR_array);
imagesc(alpha_array, TR_array,  native_t1_fat_array);
set(gca,'XTick',10:10:90);
set(gca,'YTick',2:1:10);
xlabel('Flip Angle (degree)'); ylabel('TR (ms)');
hold on;

[C,h] = contour(X,Y,native_t1_fat_array, 'k', 'ShowText','on'); 
clabel(C,h,'FontSize',16)
set(gca, 'FontSize', 18)

%% Mixture
%%%% TODO

Mxy_readout_array_comp_PhaseEnc = zeros(numel(TR_array), numel(alpha_array), numel(TI_array));
native_t1_array_comp = zeros(numel(TR_array), numel(alpha_array), numel(FF));
g = fittype('a-b*exp(-c*x)');
FF = [0:1:8]/8;

for f = 1:numel(FF)
    ff = FF(f);
    Mxy_readout_comp = ff * Mxy_readout_fat3t + (1-ff) * Mxy_readout;
    for i = 1:length(TR_array)
        x=squeeze(TI_array_sorted(i,:)).';
        for j = 1:length(alpha_array)
            Mxy_readout_array = MOLLI_readout_reorder(Mxy_readout_comp(i,j,:));
            Mxy_readout_array_comp_PhaseEnc(i,j,:) = Mxy_PhaseEnc(Mxy_readout_array);
            y=squeeze(Mxy_readout_array_comp_PhaseEnc(i,j,:));
            try
                f0 = fit(x,y,g,'StartPoint',[.5; 0; 0.001]);
                coef = coeffvalues(f0);
                native_t1_array_comp(i, j, f) = 1/coef(3) * (coef(2) / coef(1) - 1);
            catch
                disp('Fitting condition not suffice, set to NaN');
                native_t1_array_comp(i, j, f) = NaN;
            end
            
        end
    end
end
%%
figure();
[X,Y] = meshgrid(alpha_array, TR_array);


for i = 1:numel(FF)
    ff = FF(i);
    subplot(3,3,i)
    imagesc(alpha_array, TR_array,  native_t1_array_comp(:,:,i));
    set(gca,'XTick',10:10:90);
    set(gca,'YTick',2:1:10);
    xlabel('Flip Angle (degree)'); ylabel('TR (ms)');
    xlabel('Flip Angle (degree)'); ylabel('TR (ms)');
    hold on;
    [C,h] = contour(X,Y, native_t1_array_comp(:,:,i), 'k', 'ShowText','on');
    clabel(C,h,'FontSize',16)
    set(gca, 'FontSize', 18)
    title(cat(2, 'FF = ', num2str(ff)));
    colorbar;
    caxis([1000 1500]);
end