clear all;
close all;
% Xinheng Zhang 
% 04/30/2024
addpath('../EPGX-src/')
addpath('../BlochSimDemo/');
%% GRE simulation
T1 = 1000;
T2 = 10;
TE_array = [2.56, 5.80, 9.90, 15.56, 21.22];
MxyTE_cell = cell(1,length(TE_array));

for i = 1:length(TE_array)
    
    f_vec = [-300:5:300];
    M0 = [0 0 1]';
    FA = 5;
    NTR = 192;
    TR = 80;
    TE = TE_array(i);
    dt = 0.02;
    PLOT_EACHSPIN = 0;
    PLOT_SS = 0;
    [MxyTE, BSIM, plotTEidx, SPGRcat, Msimf] = SPGR_engine(M0, T1, T2, TR, TE, FA, NTR, dt, f_vec, PLOT_EACHSPIN, PLOT_SS);
    MxyTE_cell{i} = MxyTE;
end
%% 
t = [double((1:1:(BSIM.NstepTot+SPGRcat.NstepTot)))*dt];
df = 0;
idx = find(f_vec == df);
figure();
plot(t, Msimf(1,:, idx), t, Msimf(2, :, idx), t, Msimf(3, :, idx));
legend('Mx', 'My', 'Mz'); ylim([-1 1]); grid on;
xlabel('ms'); ylabel('Normalized Magnetization');

%% Show Off-resonance vs Phase
% df = 300;
% idx = find(f_vec == df);
% figure();
% plot(abs(MxyTE(NTR,:))); 


figure();
plot(angle(MxyTE(NTR,:))); 

%% Uniform distribution
df1 = 0;
df2_array = 0:5:50;
T2star_array = zeros(1, length(df2_array));
idx1 = find(f_vec == df1);

figure();
for df2_idx = 1:length(df2_array)

    df2 = df2_array(df2_idx);
    idx2 = find(f_vec == df2);

    MxyTE_ensemble = zeros(1,length(TE_array));

    for i = 1:length(TE_array)
        MxyTE_ensemble(i) = sum(MxyTE_cell{i}(NTR,idx1:idx2)) / (idx2 - idx1 + 1);
    end

    hold on;
    plot(TE_array, abs(MxyTE_ensemble), 'LineWidth', 2);
    xlim([0 30]);

    f = fit(TE_array',abs(MxyTE_ensemble)','exp1'); % a*exp(b*x)
    T2star_array(df2_idx) = -1/f.b;
end

%% Normal Distribution
df1 = 0;
df2_array = 0:5:50;
T2star_array = zeros(1, length(df2_array));
idx1 = find(f_vec == df1);

% Normal distribution pdf
pd = makedist('Normal');

figure();
for df2_idx = 1:length(df2_array)

    df2 = df2_array(df2_idx);
    idx2 = find(f_vec == df2);
    
    x = df1:5:df2;
    x_normalized = x/5 - mean(x/5);
    pdf_normal = pdf(pd,x_normalized);

    if df2_idx ~= 1
        pdf_normal_to_1 = 1/sum(pdf_normal) .* pdf_normal;
    else
        pdf_normal_to_1 = 1;
    end

    MxyTE_ensemble = zeros(1,length(TE_array));

    for i = 1:length(TE_array)
        MxyTE_ensemble(i) = MxyTE_cell{i}(NTR,idx1:idx2) * pdf_normal_to_1.';
    end

    hold on;
    plot(TE_array, abs(MxyTE_ensemble), 'LineWidth', 2);
    xlim([0 30]);

    f = fit(TE_array',abs(MxyTE_ensemble)','exp1'); % a*exp(b*x)
    T2star_array(df2_idx) = -1/f.b;
end
%%

% Values from Brinna et al 2022
R2star_remote_B = 25;
R2star_hemo_B = 106;

T2star_remote_B = 1000/R2star_remote_B; % 
T2star_hemo_B = 1000/R2star_hemo_B; % 

lambda_deltaB0 = 1./T2star_array - 1/T2;
chi = 110 * 10^-9;
gamma = 2.67 * 10^8; %rad/T/sec
f_gamma = 42.58 * 3 * 10^6; % 42.58 MHz/T * 3;

figure();
plot(df2_array, lambda_deltaB0, 'LineWidth', 2);



(R2star_hemo_B - 1/T2)/f_gamma * 10^9 % = 830 ppb / 0.830 ppm
(R2star_remote_B - 1/T2)/f_gamma * 10^9 % = 195 ppb / 0.195 ppm

mdl = fitlm(df2_array, lambda_deltaB0);
hold on;
plot(mdl)

%% The steady-state signal
rho = 1;
alpha = FA *(pi/180);
Sig = rho * sin(alpha) * (1 - exp(-TR/T1)) / (1 - cos(alpha)*exp(-TR/T1)) * exp(-TE/T2)
% fits well with perfect spoiling

%% Mz from echo 1 to echo 192
figure();
plot3(Msimf(1,:, idx), Msimf(2,:, idx), Msimf(3,:, idx), '.-'); axis equal;
grid on;
hold on;
for i = 1:numel(t)
    hline = plot3([0 Msimf(1,i, idx)], [0 Msimf(2,i, idx)], [0 Msimf(3,i, idx)], 'r');
    xlabel('Mx'); ylabel('My');zlabel('Mz')
    xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);
    pause(.02)
    set(hline,'visible','off');
end