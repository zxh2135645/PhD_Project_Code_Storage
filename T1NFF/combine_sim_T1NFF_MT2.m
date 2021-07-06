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
T1_array = 1000:20:1500;
for f = 1:length(F_array)
%for f = 1:1
    for k = 1:length(T1_array)
    %for k = 1:1
        T1 = T1_array(k);
        MT_para_remote = struct;
        MT_para_remote.T1x = [T1 T1];
        MT_para_remote.T2x = [T2, T2_bound];
        % MT_para_remote.b1sqrdtau_array = b1sqrdtau_array;
        MT_para_remote.F = F_array(f);
        MT_para_remote.Kf = 5.2e-3; % ms-1
        MT_para_remote.trf = 0.600; % ms

        MT_prep = struct;
        MT_prep.flip = d2r(flip);
        MT_prep.t_delay = t_delay;
        % Assuming pulse duration is 20 ms
        trf_prep = 20.00;
        alpha_inv = 180;
        MT_prep.B1SqrdTau = (d2r(alpha_inv)./(trf_prep.*gam)).^2.*trf_prep; 
        ddt = 0.2;

        % Looping F and Kf (20x20)
        M0_remote = [0 0 1-MT_para_remote.F MT_para_remote.F]';
        
        disp(['F: ', num2str(F_array(f)), '  T1: ', num2str(T1_array(k))]);
        tic;
        [t_total, Mzmt_total_total, Mxymt_total_total, t_readout_mt, Mxy_readout_mt] = seq_T1MOLLI_MT2(TI_array, TD, npulse,...
            alpha, TR, MT_para_remote, MT_prep, num_rampup, M0_remote, restore_pulse, trigger, trigger2, ddt);
        toc;
        
        if f == 1 && k == 1
            Mzmt_dict = zeros(length(F_array), length(T1_array), length(Mzmt_total_total));
            Mxymt_dict = zeros(length(F_array), length(T1_array), length(Mxymt_total_total));
        end
        
        Mzmt_dict(f, k, :) = Mzmt_total_total;
        Mxymt_dict(f, k, :) = Mxymt_total_total;
        
    end
end

%% Save as dictionary
Dict.Mzmt_dict = Mzmt_dict;
Dict.Mxymt_dict = Mxymt_dict;
Dict.ddt = ddt;
Dict.Mxy_readout_mt = Mxy_readout_mt;
Dict.t_readout_mt = t_readout_mt;
save_dir = uigetdir;
f_save = 'MT_MOLLI_Dict2.mat';
save(cat(2, save_dir, '/', f_save), '-struct', 'Dict');

%% Plot 1
f = 1;
k = 10;
Mzmt_total_total = squeeze(Mzmt_dict(f,k,:));
Mxymt_total_total = squeeze(Mxymt_dict(f,k,:));

figure();
subplot(2,1,1);
plot(Mzmt_total_total);
subplot(2,1,2);
plot(abs(Mxymt_total_total));
% Should discard Kf and introduce FF
%% Plot 2
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

f = 1;
k = 1;
Mzmt_total_total = squeeze(Mzmt_dict(f,k,:));
Mxymt_total_total = squeeze(Mxymt_dict(f,k,:));
Mzmt_readout = Mzmt_total_total(readout_array);
Mxymt_readout = Mxymt_total_total(readout_array);

%f = 1;
%k = 1;
native_t1_mat = zeros(length(F_array), length(T1_array));
for f = 1:length(F_array)
    for k = 1:length(T1_array)
        %f = F_array(i);
        %k = T1_array(j);
        %figure();
        %plot(squeeze(Mzmt_readout));
        Mzmt_total_total = squeeze(Mzmt_dict(f,k,:));
        Mxymt_total_total = squeeze(Mxymt_dict(f,k,:));
        Mzmt_readout = Mzmt_total_total(readout_array);
        Mxymt_readout = Mxymt_total_total(readout_array);
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
        TI_array_sorted = sort(TI_array)+npulse*TR/2;
        g = fittype('a-b*exp(-c*x)');
        f0 = fit(TI_array_sorted',Mxymt_readout_reordered_sign,g,'StartPoint',[.0;.0; 0.001]);
        coef = coeffvalues(f0);
        native_t1_mat(f,k) = 1/coef(3) * (coef(2) / coef(1) - 1);
         %native_t1 = 1/coef(3) * (coef(2) / coef(1) - 1)
        
        %xx = linspace(1,3500,100);
        %figure();
        %plot(TI_array_sorted',Mxymt_readout_reordered_sign,'ro',xx,f0(xx),'b-', 'LineWidth', 1.5);
    end
end
%% Plot heat map of T1 with only MT simulation
figure();
imagesc(native_t1_mat);
xt = get(gca, 'XTick');
xtlbl = round(linspace(1000+500/26*xt(1), 1000+500/26*xt(end), numel(xt)));
set(gca, 'XTick', xt, 'XTickLabel', xtlbl);
xlabel('T1 (ms)');
yt = get(gca, 'YTick');
ytlbl = round(linspace(0+0.10/51*yt(1), 0+0.10/51*yt(end), numel(yt)), 3);
set(gca, 'YTick', yt, 'YTickLabel', ytlbl);
ylabel('FF');