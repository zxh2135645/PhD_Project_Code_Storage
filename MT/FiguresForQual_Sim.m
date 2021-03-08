%% Simulation %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Theoretical magnetization for different BW
clear all; close all;
addpath(genpath('lib'));
addpath(genpath('EPGX-src'));

TR = 3.18:0.1:4.38;
T1_remote = 1250;
T2 = 45; 
alpha = 45;
E1 = exp(-TR / T1_remote);
E2 = exp(-TR / T2);
rho = 1;
M_remote = rho * sin(d2r(alpha)) * (1 - E1) .* sqrt(E2) ./ (1 - (E1 - E2) .* cos(d2r(alpha)) - E1 .* E2);

T1_MI = 1470;
T2 = 45;
E1 = exp(-TR / T1_MI);
E2 = exp(-TR / T2);
rho = 1;

M_mi = rho * sin(d2r(alpha)) * (1 - E1) .* sqrt(E2) ./ (1 - (E1 - E2) .* cos(d2r(alpha)) - E1 .* E2);

%% MT parameters initiation
%%% Relaxation parameters: exchange
T1x = [T1_remote T1_remote]; % ms
T2x = [T2 8.1e-3]; % ms
T2_MT = T2;
kx = 5.2e-3; % ms
fx = 0.097;

T1y = [T1_MI T1_MI]; % ms
T2y = [T2 12.9e-3];
ky = 4.1e-3;
fy = 0.07;

%%% RF saturation factor for MT
[ff,G] = SuperLorentzian(T2x(2)*1e-3);% us super-Lorentzian absorption lineshape
% [ff_gauss,G_gauss] = GaussianLineShape(T2x(2)*1e-3);

gam = 267.5221 *1e-3; %< rad /ms /uT
trf = linspace(0.800, 2.000, length(TR)); % ms
b1 = d2r(alpha)./(trf.*gam); % uT
b1sqrdtau = 2^2*b1.^2.*trf;
b1sqrdtau2 = zeros(1, length(b1sqrdtau));
idx = find(ff == min(abs(ff)));
G = G(idx);


% Gloor analytical solution
% analytic form
mxyGloor = zeros(length(TR), 1);
mxyGloor_noMT = zeros(length(TR), 1);
mxyGloor_MI = zeros(length(TR), 1);
for i = 1:length(TR)
    tr = TR(i);
    b1sqrdtau_iter = b1sqrdtau(i);
    mxyGloor(i) = ssSSFP_Gloor(d2r(alpha),b1sqrdtau_iter, tr, T1x,T2_MT,fx,kx,G);
    mxyGloor_MI(i) = ssSSFP_Gloor(d2r(alpha),b1sqrdtau_iter, tr, T1y,T2_MT,fy,ky,G);
    mxyGloor_noMT(i) = ssSSFP_Gloor(d2r(alpha),b1sqrdtau2(i), tr, T1x,T2_MT,fx,kx,G);
end


%% EPG simulations
% Either EPG simulation and Bloch-McConnell simulation does not take RF
% duration (Power deposition) into account
npulse=500;
dw=0;
mxymt = zeros(2*npulse-1, npulse, length(TR));
ssmt = zeros(length(TR), 1);

for i = 1:length(TR)
    tr = TR(i);
    b1sqrdtau_iter = b1sqrdtau(i);
    phi = RF_phase_cycle(npulse,'balanced'); % alpha, -alpha, alpha, -alpha ...
    % EPGX-BM
    % There is no input to specify RF duration
    [~,fnmt] = EPGX_GRE_MT(d2r(alpha)*ones(npulse,1),phi,b1sqrdtau_iter*ones(npulse,1),...
    tr,T1x,T2_MT,fx,kx,G,'kmax',inf);
    mxymt(:,:,i) = size(fnmt,1)*ifftshift(ifft(ifftshift(fnmt,1),[],1),1);
    % mxyx(:,:,:,i) = ifftshift(fft(ifftshift(fnx,1),[],1),1);
    % mxyx1(:,:,i) = sum(mxyx(:,:,:,i),3); % sum compartments
    ssmt(i) = ssSSFP_MT(d2r(alpha),b1sqrdtau_iter,tr,dw,T1x,T2_MT,fx,kx,G); % Bloch-McConnell simulation to steady-state
end

mxymt_array = zeros(length(TR), 1);
for i = 1:length(TR)
    mxymt_array(i) = abs(mxymt(floor(size(mxymt, 1)/2),end,i));
end

%% Plots
figure();
plot(TR, M_remote, 'LineWidth', 2);
grid on;
xlabel('TR (ms)'); ylabel('Magnetization (M_{xy})');
hold on;
% plot(TR, M_mi, 'LineWidth', 2);
% title('Theoretical Magnetization for Remote and MI');
plot(TR, mxyGloor, 'LineWidth', 2)
% plot(TR, mxyGloor_MI, 'LineWidth', 2)
plot(TR, mxymt_array, '--', 'LineWidth', 2)
% plot(TR, abs(ssmt), 'LineWidth', 2)
ylim([0 0.08])
legend({'Remote', 'Remote(MT) analytic', 'Remote(MT) EPG-X'}, 'Location', 'SouthEast');
% legend({'Remote', 'MI', 'Remote(MT) analytic', 'MI(MT) analytic'});
set(gca,'FontSize',16)
x0=10;
y0=10;
width=550;
height=300;
set(gcf,'position',[x0,y0,width,height])
disp(cat(2, "MTR of remote is:", num2str((mxyGloor(end) - mxyGloor(1)) / (mxyGloor(end)))))
disp(cat(2, "MTR of MI is:", num2str((mxyGloor_MI(end) - mxyGloor_MI(1)) / (mxyGloor_MI(end)))))

%% Range of some parameters for MI simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% exchange rate varies
% Theorectical steady-state (no MT)
T1_MI = 1470;
T2 = 45;
alpha = 45;
tr = [3.18, 4.38]; % ms
E1 = exp(-tr / T1_MI);
E2 = exp(-tr / T2);
rho = 1;

M_mi = rho * sin(d2r(alpha)) * (1 - E1) .* sqrt(E2) ./ (1 - (E1 - E2) .* cos(d2r(alpha)) - E1 .* E2);

gam = 267.5221 *1e-3;
T1y = [T1_MI T1_MI]; % ms
T2y = [T2 8.1e-3];
T2_MT = T2;
ky = (0.1:0.1:10)*10^(-3);
fy = 0.097;

trf = [0.800, 2.000]; % ms
b1 = d2r(alpha)./(trf.*gam); % uT
b1sqrdtau_iter = 2^2*b1.^2.*trf;
% Get the lineshape
[ff,G] = SuperLorentzian(T2y(2)*1e-3);
idx = find(ff == min(abs(ff)));
G = G(idx);

mxyGloor_MI_fast = zeros(length(ky), 1);
mxyGloor_MI_lowsar = zeros(length(ky), 1);
for i = 1:length(ky)
    ky_iter = ky(i);
    mxyGloor_MI_fast(i) = ssSSFP_Gloor(d2r(alpha), b1sqrdtau_iter(1), tr(1), T1y, T2_MT, fy, ky_iter, G);
    mxyGloor_MI_lowsar(i) = ssSSFP_Gloor(d2r(alpha), b1sqrdtau_iter(2), tr(2), T1y, T2_MT, fy, ky_iter, G);
end

figure();
M_mi_array = linspace(M_mi(1), M_mi(2), length(ky));
plot(ky, M_mi_array, 'LineWidth', 2)
hold on;
plot(ky, mxyGloor_MI_fast, 'LineWidth', 2)
plot(ky, mxyGloor_MI_lowsar, 'LineWidth', 2)
grid on;
ylabel('Magnetization');
yyaxis right
mtr = (mxyGloor_MI_lowsar - mxyGloor_MI_fast) ./ mxyGloor_MI_lowsar;
plot(ky, mtr, 'LineWidth', 2)
legend({'MI no MT', 'MI Fast', 'MI LowSAR', 'MI MTR'});
xlabel('K_f s^{-1}');ylabel('MTR');
%% macro-molecule pool size varies
% Theorectical steady-state (no MT)
T1_MI = 1470;
T2 = 45;
alpha = 45;
tr = [3.18, 4.38]; % ms
E1 = exp(-tr / T1_MI);
E2 = exp(-tr / T2);
rho = 1;

M_mi = rho * sin(d2r(alpha)) * (1 - E1) .* sqrt(E2) ./ (1 - (E1 - E2) .* cos(d2r(alpha)) - E1 .* E2);

gam = 267.5221 *1e-3;
T1y = [T1_MI 1000]; % ms
T2y = [T2 8.1e-3];
T2_MT = T2;
ky = 5.2e-3;
fy = 0.01:0.01:1;

trf = [0.800, 2.000]; % ms
b1 = d2r(alpha)./(trf.*gam); % uT
b1sqrdtau_iter = 2^2*b1.^2.*trf;
% Get the lineshape
[ff,G] = SuperLorentzian(T2y(2)*1e-3);
idx = find(ff == min(abs(ff)));
G = G(idx);

mxyGloor_MI_fast = zeros(length(fy), 1);
mxyGloor_MI_lowsar = zeros(length(fy), 1);
for i = 1:length(fy)
    fy_iter = fy(i);
    mxyGloor_MI_fast(i) = ssSSFP_Gloor(d2r(alpha), b1sqrdtau_iter(1), tr(1), T1y, T2_MT, fy_iter, ky, G);
    mxyGloor_MI_lowsar(i) = ssSSFP_Gloor(d2r(alpha), b1sqrdtau_iter(2), tr(2), T1y, T2_MT, fy_iter, ky, G);
end

figure();
M_mi_array = linspace(M_mi(1), M_mi(2), length(fy));
plot(fy, M_mi_array, 'LineWidth', 2)
hold on;
plot(fy, mxyGloor_MI_fast, 'LineWidth', 2)
plot(fy, mxyGloor_MI_lowsar, 'LineWidth', 2)
grid on;
ylabel('Magnetization');
yyaxis right
mtr = (mxyGloor_MI_lowsar - mxyGloor_MI_fast) ./ mxyGloor_MI_lowsar;
plot(fy, mtr, 'LineWidth', 2)
legend({'MI no MT', 'MI Fast', 'MI LowSAR', 'MI MTR'});
xlabel('F');ylabel('MTR');
%% T2r varies
% Theorectical steady-state (no MT)
T1_MI = 1470;
T2 = 45;
alpha = 45;
tr = [3.18, 4.38]; % ms
E1 = exp(-tr / T1_MI);
E2 = exp(-tr / T2);
rho = 1;

M_mi = rho * sin(d2r(alpha)) * (1 - E1) .* sqrt(E2) ./ (1 - (E1 - E2) .* cos(d2r(alpha)) - E1 .* E2);

gam = 267.5221 *1e-3;
T1y = [T1_MI 1000]; % ms
T2y = [T2 8.1e-3];
T2r = (1:0.1:15) * 10^(-3);
T2_MT = T2;
ky = 5.2e-3;
fy = 0.097;

trf = [0.800, 2.000]; % ms
b1 = d2r(alpha)./(trf.*gam); % uT
b1sqrdtau_iter = 2^2*b1.^2.*trf;
% Get the lineshape
G_array = zeros(length(T2r), 1);
for i = 1:length(T2r)
    [ff,G] = SuperLorentzian(T2r(i)*1e-3);
    idx = find(ff == min(abs(ff)));
    G_array(i) = G(idx);
end

mxyGloor_MI_fast = zeros(length(T2r), 1);
mxyGloor_MI_lowsar = zeros(length(T2r), 1);
for i = 1:length(T2r)
    G = G_array(i);
    mxyGloor_MI_fast(i) = ssSSFP_Gloor(d2r(alpha), b1sqrdtau_iter(1), tr(1), T1y, T2_MT, fy, ky, G);
    mxyGloor_MI_lowsar(i) = ssSSFP_Gloor(d2r(alpha), b1sqrdtau_iter(2), tr(2), T1y, T2_MT, fy, ky, G);
end

figure();
M_mi_array = linspace(M_mi(1), M_mi(2), length(T2r));
plot(T2r, M_mi_array, 'LineWidth', 2)
hold on;
plot(T2r, mxyGloor_MI_fast, 'LineWidth', 2)
plot(T2r, mxyGloor_MI_lowsar, 'LineWidth', 2)
grid on;
ylabel('Magnetization');
yyaxis right
mtr = (mxyGloor_MI_lowsar - mxyGloor_MI_fast) ./ mxyGloor_MI_lowsar;
plot(T2r, mtr, 'LineWidth', 2)
legend({'MI no MT', 'MI Fast', 'MI LowSAR', 'MI MTR'});
xlabel('T2_r s');ylabel('MTR');
%% T1r varies
% Theorectical steady-state (no MT)
clear G;
T1_MI = 1470;
T2 = 45;
alpha = 45;
tr = [3.18, 4.38]; % ms
E1 = exp(-tr / T1_MI);
E2 = exp(-tr / T2);
rho = 1;

M_mi = rho * sin(d2r(alpha)) * (1 - E1) .* sqrt(E2) ./ (1 - (E1 - E2) .* cos(d2r(alpha)) - E1 .* E2);

gam = 267.5221 *1e-3; % rad /ms /uT
T1r = 500:10:1500;
T1y = zeros(length(T1r),2);
for i = 1:length(T1r)
    T1y(i,:) = [T1_MI, T1r(i)];
end
% ms
T2y = [T2 8.1e-3];
T2_MT = T2;
ky = 5.2e-3;
fy = 0.097;

trf = [0.800, 2.000]; % ms
b1 = d2r(alpha)./(trf.*gam); % uT
b1sqrdtau_iter = 2^2*b1.^2.*trf;
% Get the lineshape
[ff,G] = SuperLorentzian(T2y(2)*1e-3);
idx = find(ff == min(abs(ff)));
G = G(idx);

mxyGloor_MI_fast = zeros(length(T1r), 1);
mxyGloor_MI_lowsar = zeros(length(T1r), 1);
for i = 1:length(T1r)
    T1y_iter = T1y(i,:);
    mxyGloor_MI_fast(i) = ssSSFP_Gloor(d2r(alpha), b1sqrdtau_iter(1), tr(1), T1y_iter, T2_MT, fy, ky, G);
    mxyGloor_MI_lowsar(i) = ssSSFP_Gloor(d2r(alpha), b1sqrdtau_iter(2), tr(2), T1y_iter, T2_MT, fy, ky, G);
end
figure();
M_mi_array = linspace(M_mi(1), M_mi(2), length(T1r));
plot(T1r, M_mi_array, 'LineWidth', 2)
hold on;
plot(T1r, mxyGloor_MI_fast, 'LineWidth', 2)
plot(T1r, mxyGloor_MI_lowsar, 'LineWidth', 2)
grid on;
ylabel('Magnetization');
yyaxis right
mtr = (mxyGloor_MI_lowsar - mxyGloor_MI_fast) ./ mxyGloor_MI_lowsar;
plot(T1r, mtr, 'LineWidth', 2)
legend({'MI no MT', 'MI Fast', 'MI LowSAR', 'MI MTR'});
xlabel('T1_r s');ylabel('MTR');
%% Kf and F varies
% Theorectical steady-state (no MT)
T1_MI = 1470;
T2 = 45;
alpha = 45;
tr = [3.18, 4.38]; % ms
E1 = exp(-tr / T1_MI);
E2 = exp(-tr / T2);
rho = 1;

M_mi = rho * sin(d2r(alpha)) * (1 - E1) .* sqrt(E2) ./ (1 - (E1 - E2) .* cos(d2r(alpha)) - E1 .* E2);

gam = 267.5221 *1e-3;
T1y = [T1_MI T1_MI]; % ms
T2y = [T2 8.1e-3];
T2_MT = T2;
ky = (0.1:0.1:100)*10^(-3);
fy = 0.01:0.01:1;

trf = [0.800, 2.000]; % ms
b1 = d2r(alpha)./(trf.*gam); % uT
b1sqrdtau_iter = 2^2*b1.^2.*trf;
% Get the lineshape
[ff,G] = SuperLorentzian(T2y(2)*1e-3);
idx = find(ff == min(abs(ff)));
G = G(idx);

mxyGloor_MI_fast = zeros(length(ky), length(fy));
mxyGloor_MI_lowsar = zeros(length(ky), length(fy));
for i = 1:length(ky)
    ky_iter = ky(i);
    for j = 1:length(fy)
        fy_iter = fy(j);
        mxyGloor_MI_fast(i, j) = ssSSFP_Gloor(d2r(alpha), b1sqrdtau_iter(1), tr(1), T1y, T2_MT, fy_iter, ky_iter, G);
        mxyGloor_MI_lowsar(i, j) = ssSSFP_Gloor(d2r(alpha), b1sqrdtau_iter(2), tr(2), T1y, T2_MT, fy_iter, ky_iter, G);
    end
end


mtr = (mxyGloor_MI_lowsar - mxyGloor_MI_fast) ./ mxyGloor_MI_lowsar;
figure('Renderer', 'painters', 'Position', [10 10 400 300])
imagesc(mtr);colormap jet;
c = colorbar;
x1 = get(gca,'position');
x = get(c,'Position');
x(3) = 0.05;
set(c,'Position',x)
set(gca,'position',x1)
%c.Label.String = 'MTR';
hold on;
[C, h] = contour(mtr, 'ShowText','on');
h.LineColor = [0.2 0.2 0.2];
h.LineWidth = 1;
clabel(C,h,'FontSize',14,'Color',[0.2 0.2 0.2])
xlabel('F (%)'); ylabel('K_f (s^{-1})');
set(gca,'XTickLabel',[20 40 60 80 100]);
set(gca,'YTickLabel',[20 40 60 80 100]);
set(gca, 'FontSize', 14)
save_dir = 'C:\Users\ZhangX1\Dropbox\James\PhD_Qualifying\img\';
print(cat(2, save_dir, 'fig_MTR.png') ,'-dpng', '-r300')

%% Simulation of the sequence for methods session
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read approx T1 values of the infarct
base_dir = 'D:\Data\London_Canine\';
read_dir = cat(2, base_dir, 'MERRY_FEPPA_MERRY_26JUL19\HEART_CEDARS_20190726_131018_242000\');
save_dir = cat(2, base_dir, 'MERRY_FEPPA_MERRY_26JUL19\Figures\');
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

glob_str = 'T2MAP_*T2*';
[img_flow, hdr_flow] = ImgPick(read_dir, glob_str);
%% Display image and Draw ROI
for_mask = img_flow; % at the LowSAR
rescaled_for_mask = for_mask / (max(for_mask(:)) - min(for_mask(:)));
figure();
roi = roipoly(rescaled_for_mask);
%% T1 value
% roi_mi = roi;
mean(nonzeros(roi.*img_flow))
% T1 value is 1384.7 ms
% T1 value of remote is 1133.2 ms
% T2 value of mi is 34 ms
% T2 vslue of remote is 35 ms
%% Starting from here %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% simulation without MT
clear all; close all;
addpath(genpath('lib'));
addpath(genpath('EPGX-src'));
TD = 2212.5; % ms
num_rampup = 10;
npulse = 54 + num_rampup; % 10 linear ramp-up pulses
t_delay = 650; % TI = 650 ms
flip = 180; % prep flip angle = 180 degree
T1 = 1133.2;
T2 = 45;
alpha = 90;
% alpha, -alpha, alpha, -alpha ...
phi = RF_phase_cycle(npulse,'balanced');
TR = 3;
prep = struct;
prep.flip = d2r(flip);
prep.t_delay = t_delay;
fa_flow = d2r(alpha)*ones(npulse,1);

rampup = linspace(0, alpha, num_rampup+1);
rampup = rampup(end-num_rampup+1:end);
fa_flow(1:num_rampup) = d2r(rampup);
[~,fn0, Zn0] = EPG_GRE(fa_flow,phi,TR,T1,T2,'kmax',inf, 'prep', prep);
rho = 1;
anal = rho * (1 - 2 * exp(-t_delay / T1));

mxys = size(fn0,1)*ifftshift(ifft(ifftshift(fn0,1),[],1),1);
mzs = size(Zn0,1)*ifftshift(ifft(ifftshift(Zn0,1),[],1),1);
mxys_array = abs(mxys(floor(size(mxys, 1)/2),:));
mzs_array = real(mzs(floor(size(mzs, 1)/2+1),:));
figure();
subplot(2,2,1)
plot(mxys_array, 'LineWidth', 2)
title('1st shot Mxy')
xlabel('Number of TR'); ylabel('M_{xy}')
grid on;
subplot(2,2,2)
plot(mzs_array, 'LineWidth', 2)
title('1st shot Mz')
xlabel('Number of TR'); ylabel('M_{z}')
grid on;

% second shot
Mz = mzs_array(end);
Mz2 = 1 - (1 - Mz).*exp(-TD/T1);
initM = [0 0 Mz2]';
[~,fn02, Zn02] = EPG_GRE(fa_flow,phi,TR,T1,T2,'kmax',inf, 'prep', prep, 'initM', initM);
mxys2 = size(fn02,1)*ifftshift(ifft(ifftshift(fn02,1),[],1),1);
mzs2 = size(Zn02,1)*ifftshift(ifft(ifftshift(Zn02,1),[],1),1);
mxys_array2 = abs(mxys2(floor(size(mxys2, 1)/2),:));
mzs_array2 = real(mzs2(floor(size(mzs2, 1)/2+1),:));
subplot(2,2,3)
plot(mxys_array2, 'LineWidth', 2)
title('2nd shot Mxy')
xlabel('Number of TR'); ylabel('M_{xy}')
grid on;
subplot(2,2,4)
plot(mzs_array2, 'LineWidth', 2)
title('2nd shot Mz')
xlabel('Number of TR'); ylabel('M_{z}')
grid on;

% For better visualization
dt = TR;
total_time = TD + TR * npulse + t_delay;
one_rep_t = 0:dt:total_time;
one_rep_Mz = ones(1, length(one_rep_t));

one_rep_Mz(end-npulse+1:end) = mzs_array;


TI_timepoint = length(one_rep_Mz(end-npulse-round(t_delay/dt)+1:end-npulse));
E1 = exp(-dt / T1);
Mz_TI_array = zeros(1, TI_timepoint);
for i = 1:TI_timepoint
    if i == 1
        Mz_TI_array(i) = -1;
    else
        Mz_TI_array(i) = 1 - E1 + Mz_TI_array(i-1) * E1;
    end
end

% Linear Curves
% Mz_TI_array = linspace(-1, one_rep_Mz(end-npulse+1),TI_timepoint);

% Add inversion pulse
one_rep_Mz(end-npulse-round(t_delay/dt)+1:end-npulse) = Mz_TI_array;

figure();
plot(one_rep_t, one_rep_Mz, 'LineWidth', 2)
grid on;

% Second rep
two_rep_t = 0:dt:total_time;
two_rep_Mz = zeros(1, length(two_rep_t));
% two_rep_Mz(end-npulse+1:end) = mzs_array2;

% Trigger Delay
TD_timepoint = length(two_rep_Mz(1:end-npulse-round(t_delay/dt)));
Mz_TD_array2 = zeros(1, TD_timepoint);
for i = 1:TD_timepoint
    if i == 1
        Mz_TD_array2(i) = 1 - E1 + one_rep_Mz(end) * E1;
    else
        Mz_TD_array2(i) = 1 - E1 + Mz_TD_array2(i-1) * E1;
    end
end

% Inversion time
TI_timepoint = length(two_rep_Mz(end-npulse-round(t_delay/dt)+1:end-npulse));
Mz_TI_array2 = zeros(1, TI_timepoint);
for i = 1:TI_timepoint
    if i == 1
        Mz_TI_array2(i) = 1 - E1 + (-Mz_TD_array2(end)) * E1;
    else
        Mz_TI_array2(i) = 1 - E1 + Mz_TI_array2(i-1) * E1;
    end
end

two_rep_Mz = [Mz_TD_array2, Mz_TI_array2, mzs_array2];
figure();
plot(two_rep_t, two_rep_Mz, 'LineWidth', 2)
grid on;

t = [one_rep_t, one_rep_t(end)+two_rep_t];
Mz_total = [one_rep_Mz, two_rep_Mz];
figure();
plot(t, Mz_total, 'LineWidth', 2)
grid on;
%% simulation with MT
MT_para_remote = struct;
MT_para_remote.T1x = [T1 T1];
MT_para_remote.T2x = [T2, 8.1e-3];
gam = 267.5221 *1e-3; % rad /ms /uT
trf = 0.600; % ms
b1 = d2r(alpha)./(trf.*gam); % uT
b1sqrdtau = 2^2 * b1.^2.*trf;
b1sqrdtau0 = 2^2 * (d2r(rampup)./(trf.*gam)).^2.*trf;
b1sqrdtau_array = [b1sqrdtau0'; b1sqrdtau * ones(npulse-num_rampup,1)];
MT_para_remote.b1sqrdtau_array = b1sqrdtau_array;
MT_para_remote.F = 0.097;
MT_para_remote.Kf = 5.2e-3;
[ff,G] = SuperLorentzian(MT_para_remote.T2x(2)*1e-3);
idx = find(ff == min(abs(ff)));
MT_para_remote.G = G(idx);
MT_prep = struct;
MT_prep.flip = d2r(flip);
MT_prep.t_delay = t_delay;
% Assuming pulse duration is 20 ms
MT_prep.B1SqrdTau = 2^2 * (d2r(alpha)./(20.00.*gam)).^2.*20.00; 
[~,fnmt,Znmt] = EPGX_GRE_MT(fa_flow,phi,MT_para_remote.b1sqrdtau_array,...
    TR,MT_para_remote.T1x,MT_para_remote.T2x(1),MT_para_remote.F,MT_para_remote.Kf,MT_para_remote.G,...
    'kmax',inf, 'prep', MT_prep);

mxymt = size(fnmt,1)*ifftshift(ifft(ifftshift(fnmt,1),[],1),1);
mzmt = size(Znmt,1)*ifftshift(ifft(ifftshift(Znmt,1),[],1),1);

mxymt_array = abs(mxymt(floor(size(mxymt, 1)/2),:));
mzmt_array = real(mzmt(floor(size(mzmt, 1)/2+1),:));
figure();
subplot(2,2,1)
plot(mxymt_array, 'LineWidth', 2)
title('1st shot Mxy')
xlabel('Number of TR'); ylabel('M_{xy}')
grid on;
subplot(2,2,2)
plot(mzmt_array(1:npulse), 'LineWidth', 2)
title('1st shot Mz')
xlabel('Number of TR'); ylabel('M_{z}')
grid on;


% second shot
Mzmt = mzmt_array(npulse);
Mzmt2 = 1 - (1 - Mzmt).*exp(-TD/T1); % This is not accurate
initMmt = [0 0 Mzmt2]';
[~,fnmt2,Znmt2] = EPGX_GRE_MT(fa_flow,phi,MT_para_remote.b1sqrdtau_array,...
    TR,MT_para_remote.T1x,MT_para_remote.T2x(1),MT_para_remote.F,MT_para_remote.Kf,MT_para_remote.G,...
    'kmax',inf, 'prep', MT_prep, 'initMmt', initMmt);

mxymt2 = size(fnmt2,1)*ifftshift(ifft(ifftshift(fnmt2,1),[],1),1);
mzmt2 = size(Znmt2,1)*ifftshift(ifft(ifftshift(Znmt2,1),[],1),1);

mxymt_array2 = abs(mxymt2(floor(size(mxymt2, 1)/2),:));
mzmt_array2 = real(mzmt2(floor(size(mzmt2, 1)/2+1),:));
subplot(2,2,3)
plot(mxymt_array2, 'LineWidth', 2)
title('2nd shot Mxy')
xlabel('Number of TR'); ylabel('M_{xy}')
grid on;
subplot(2,2,4)
plot(mzmt_array2(1:npulse), 'LineWidth', 2)
title('2nd shot Mz')
xlabel('Number of TR'); ylabel('M_{z}')
grid on;

% For better visualization
ti_flow = zeros(TI_timepoint, 1);
phi_ti = zeros(TI_timepoint, 1);
b1sqrdtau_array = zeros(TI_timepoint,1);
MT_prep2 = MT_prep;
MT_prep2.t_delay = 0;
[~,~,Znmt_ti] = EPGX_GRE_MT(ti_flow,phi_ti,b1sqrdtau_array,...
    TR,MT_para_remote.T1x,MT_para_remote.T2x(1),MT_para_remote.F,MT_para_remote.Kf,MT_para_remote.G,...
    'kmax',inf, 'prep', MT_prep2);
mzmt_ti = size(Znmt_ti,1)*ifftshift(ifft(ifftshift(Znmt_ti,1),[],1),1);
mzmt_array_ti = real(mzmt_ti(floor(size(mzmt_ti, 1)/2),:));
figure();
plot(mzmt_array_ti(1:TI_timepoint), 'LineWidth', 2)
grid on;

one_rep_Mzmt = ones(1, length(one_rep_t));
one_rep_Mzmt(end-npulse+1:end) = mzmt_array(1:npulse);
one_rep_Mzmt(end-npulse-TI_timepoint+1:end-npulse) = mzmt_array_ti(1:TI_timepoint);
figure();
plot(one_rep_t, one_rep_Mzmt, 'LineWidth', 2)
grid on;

% For rep two
two_rep_Mz = zeros(1, length(two_rep_t));
% For better visualization
td_flow = zeros(TD_timepoint, 1);
phi_td = zeros(TD_timepoint, 1);
b1sqrdtau_array = zeros(TD_timepoint,1);
initMmt2 = [0 0 one_rep_Mzmt(end)]';
[~,~,Znmt_td] = EPGX_GRE_MT(td_flow,phi_td,b1sqrdtau_array,TR,...
    MT_para_remote.T1x,MT_para_remote.T2x(1),MT_para_remote.F,MT_para_remote.Kf,MT_para_remote.G,...
    'kmax',inf, 'initMmt', initMmt2);
mzmt_td = size(Znmt_td,1)*ifftshift(ifft(ifftshift(Znmt_td,1),[],1),1);
mzmt_array_td = real(mzmt_td(floor(size(mzmt_td, 1)/2+1),:));
figure();
plot(mzmt_array_td(1:TD_timepoint), 'LineWidth', 2)
grid on;

% The second TI
ti_flow2 = zeros(TI_timepoint, 1);
phi_ti2 = zeros(TI_timepoint, 1);
b1sqrdtau_array2 = zeros(TI_timepoint,1);
MT_prep3 = MT_prep;
MT_prep3.t_delay = 0;
initMmt3 = [0 0 mzmt_array_td(TD_timepoint)]';
[~,~,Znmt_ti2] = EPGX_GRE_MT(ti_flow2,phi_ti2,b1sqrdtau_array2,TR,...
    MT_para_remote.T1x,MT_para_remote.T2x(1),MT_para_remote.F,MT_para_remote.Kf,MT_para_remote.G,...
    'kmax',inf, 'prep', MT_prep3, 'initMmt', initMmt3);
mzmt_ti2 = size(Znmt_ti2,1)*ifftshift(ifft(ifftshift(Znmt_ti2,1),[],1),1);
mzmt_array_ti2 = real(mzmt_ti2(floor(size(mzmt_ti2, 1)/2),:));
figure();
plot(mzmt_array_ti2(1:TI_timepoint), 'LineWidth', 2)
grid on;

two_rep_Mzmt = [mzmt_array_td(1:TD_timepoint), mzmt_array_ti2(1:TI_timepoint), mzmt_array2(1:npulse)];
figure();
plot(two_rep_t, two_rep_Mzmt, 'Linewidth', 2)
grid on;

Mzmt_total = [one_rep_Mzmt, two_rep_Mzmt];
figure();
plot(t, Mz_total, '--', 'LineWidth', 2)
xlabel('Time (ms)'); ylabel('M_z')
hold on;
plot(t, Mzmt_total, 'LineWidth', 2)
grid on;
legend({'No MT', 'With MT'});

%% Fig readout magnetization
figure();
subplot(2,2,1)
plot(mxys_array, 'LineWidth', 2)
hold on;
plot(mxymt_array, 'LineWidth', 2)
title('1st shot Mxy')
xlabel('Number of TR'); ylabel('M_{xy}')
legend({'w/o MT', 'with MT'});
grid on;

subplot(2,2,2)
plot(mzs_array, 'LineWidth', 2)
hold on;
plot(mzmt_array(1:npulse), 'LineWidth', 2)
title('1st shot Mz')
xlabel('Number of TR'); ylabel('M_{z}')
legend({'w/o MT', 'with MT'});
grid on;

subplot(2,2,3)
plot(mxys_array2, 'LineWidth', 2)
hold on;
plot(mxymt_array2, 'LineWidth', 2)
title('2nd shot Mxy')
xlabel('Number of TR'); ylabel('M_{xy}')
legend({'w/o MT', 'with MT'});
grid on;

subplot(2,2,4)
plot(mzs_array2, 'LineWidth', 2)
hold on;
plot(mzmt_array2(1:npulse), 'LineWidth', 2)
title('2nd shot Mz')
xlabel('Number of TR'); ylabel('M_{z}')
legend({'w/o MT', 'with MT'});
grid on;

%% For MI %%
%% Actual starting point is here 10/01/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Made into functions
%% (without MT)
num_rampup = 10;
TD = 2212.5; % ms
npulse = 54 + num_rampup;
t_delay = 650; % TI = 650 ms
flip = 180; % prep flip angle = 180 degree
alpha = 90;
TR = 3;
T1_mi = 1384.7;
T2_mi = 34;
T1 = 1133.2;
T2 = 35;

prep = struct;
prep.flip = d2r(flip);
prep.t_delay = t_delay;


[t, Mz_total, mxys_array, mxys_array2, mzs_array, mzs_array2] = ...
    seq_plot_noMT(TD, npulse, T1, T2, alpha, TR, prep);
[~, Mzmi_total, mxymis_array, mxymis_array2, mzmis_array, mzmis_array2] = ...
    seq_plot_noMT(TD, npulse, T1_mi, T2_mi, alpha, TR, prep);

%% (with MT)

b1 = d2r(alpha)./(trf.*gam); % uT
b1sqrdtau = 2^2 * b1.^2.*trf;
b1sqrdtau0 = 2^2 * (d2r(rampup)./(trf.*gam)).^2.*trf;
b1sqrdtau_array = [b1sqrdtau0'; b1sqrdtau * ones(npulse-1,1)];

MT_para_remote = struct;
MT_para_remote.T1x = [T1 T1];
MT_para_remote.T2x = [T2, 8.1e-3];
MT_para_remote.b1sqrdtau_array = b1sqrdtau_array;
MT_para_remote.F = 0.097;
MT_para_remote.Kf = 5.2e-3;
MT_para_remote.trf = 0.600; % ms

MT_para_mi = struct;
MT_para_mi.T1x = [T1_mi T1_mi];
MT_para_mi.T2x = [T2_mi, 8.1e-3];
MT_para_mi.b1sqrdtau_array = b1sqrdtau_array;
MT_para_mi.F = 0.02;
MT_para_mi.Kf = 5.2e-3;
MT_para_mi.trf = 0.600; % ms

MT_prep = struct;
MT_prep.flip = d2r(flip);
MT_prep.t_delay = t_delay;
% Assuming pulse duration is 20 ms
trf_prep = 20.00;
MT_prep.B1SqrdTau = 2^2 * (d2r(alpha)./(trf_prep.*gam)).^2.*trf_prep; 

[~, Mzmt_total, mxymt_array, mxymt_array2, mzmt_array, mzmt_array2, Mzmt_bound_total] = ...
    seq_plot_MT(TD, npulse, alpha, TR, MT_para_remote, MT_prep);
[~, Mzmimt_total, mxymimt_array, mxymimt_array2, mzmimt_array, mzmimt_array2, Mzmimt_bound_total] = ...
    seq_plot_MT(TD, npulse, alpha, TR, MT_para_mi, MT_prep);
%% Plots
figure();
plot(t,Mzmt_total, 'LineWidth', 2)
hold on;
plot(t, Mzmimt_total, 'LineWidth', 2)
plot(t, Mz_total, '--', 'LineWidth', 2)
plot(t, Mzmi_total, '--', 'LineWidth', 2)
legend({'Remote (MT)', 'MI (MT)', 'Remote (no MT)', 'MI (no MT)'});
xlabel('Time (ms)'); ylabel('M_z')
grid on;

%% Fig readout magnetization
figure();
subplot(2,2,1)
plot(mxymt_array, 'LineWidth', 2)
hold on;
plot(mxymimt_array, 'LineWidth', 2)
plot(mxys_array, '--', 'LineWidth', 2)
plot(mxymis_array, '--', 'LineWidth', 2)
title('1st shot Mxy')
xlabel('Number of TR'); ylabel('M_{xy}')
% legend({'w/o MT', 'with MT'});
grid on;

subplot(2,2,2)
plot(mzmt_array(1:npulse), 'LineWidth', 2)
hold on;
plot(mzmimt_array(1:npulse), 'LineWidth', 2)
plot(mzs_array(1:npulse), '--', 'LineWidth', 2)
plot(mzmis_array(1:npulse), '--', 'LineWidth', 2)
title('1st shot Mz')
xlabel('Number of TR'); ylabel('M_{xy}')
% legend({'w/o MT', 'with MT'});
grid on;

subplot(2,2,3)
plot(mxymt_array2, 'LineWidth', 2)
hold on;
plot(mxymimt_array2, 'LineWidth', 2)
plot(mxys_array2, '--', 'LineWidth', 2)
plot(mxymis_array2, '--', 'LineWidth', 2)
title('2nd shot Mxy')
xlabel('Number of TR'); ylabel('M_{xy}')
% legend({'w/o MT', 'with MT'});
grid on;

subplot(2,2,4)
plot(mzmt_array2(1:npulse), 'LineWidth', 2)
hold on;
plot(mzmimt_array2(1:npulse), 'LineWidth', 2)
plot(mzs_array2(1:npulse), '--', 'LineWidth', 2)
plot(mzmis_array2(1:npulse), '--', 'LineWidth', 2)
title('2nd shot Mz')
xlabel('Number of TR'); ylabel('M_{xy}')
% legend({'w/o MT', 'with MT'});
grid on;

%% Show free pool vs bound pool
%% Plots
figure();
plot(t, Mzmimt_total, 'LineWidth', 2)
hold on;
plot(t, Mzmimt_bound_total, 'LineWidth', 2)
legend({'Free pool (MT)', 'Bound pool (MT) F = 0.02'});
xlabel('Time (ms)'); ylabel('M_z')
grid on;

%% Create a seq block 10/01/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (without MT)
num_rampup = 10;
TD = 2212.5; % ms
npulse = 60 + num_rampup; % Single-shot
t_delay = 650; % TI = 650 ms
flip = 180; % prep flip angle = 180 degree
alpha = 90;
TR = 3;
T1_mi = 1384.7;
T2_mi = 34;
T1 = 1133.2;
T2 = 35;

prep = struct;
prep.flip = d2r(flip);
prep.t_delay = t_delay;

Mz0 = 0;

[t, Mz_total, mxys_array, mzs_array] = ...
    seq_block_noMT(TD, npulse, T1, T2, alpha, TR, prep, Mz0);
[~, Mzmi_total, mxymi_array, mzmi_array] = ...
    seq_block_noMT(TD, npulse, T1_mi, T2_mi, alpha, TR, prep, Mz0);

%% Plots
figure();
plot(t, Mz_total, 'LineWidth', 2)
hold on;
plot(t, Mzmi_total, 'LineWidth', 2)

legend({'Remote (no MT)', 'MI (no MT)'});
xlabel('Time (ms)'); ylabel('M_z')
grid on;

%% (with MT)
trf = 0.600; % ms
gam = 267.5221 *1e-3; % rad /ms /uT
b1 = d2r(alpha)./(trf.*gam); % uT
b1sqrdtau = 2^2 * b1.^2.*trf;
b1sqrdtau0 = 2^2 * (d2r(rampup)./(trf.*gam)).^2.*trf;
b1sqrdtau_array = [b1sqrdtau0'; b1sqrdtau * ones(npulse-1,1)];

MT_para_remote = struct;
MT_para_remote.T1x = [T1 T1];
MT_para_remote.T2x = [T2, 8.1e-3];
MT_para_remote.b1sqrdtau_array = b1sqrdtau_array;
MT_para_remote.F = 0.097;
MT_para_remote.Kf = 5.2e-3;
MT_para_remote.trf = 0.600; % ms

MT_para_mi = struct;
MT_para_mi.T1x = [T1_mi T1_mi];
MT_para_mi.T2x = [T2_mi, 8.1e-3];
MT_para_mi.b1sqrdtau_array = b1sqrdtau_array;
MT_para_mi.F = 0.02;
MT_para_mi.Kf = 5.2e-3;
MT_para_mi.trf = 0.600; % ms

MT_prep = struct;
MT_prep.flip = d2r(flip);
MT_prep.t_delay = t_delay;
% Assuming pulse duration is 20 ms
trf_prep = 20.00;
MT_prep.B1SqrdTau = 2^2 * (d2r(alpha)./(trf_prep.*gam)).^2.*trf_prep; 

M0_remote = [0 0 1-MT_para_remote.F MT_para_remote.F]';
M0_mi = [0 0 1-MT_para_mi.F MT_para_mi.F]';

[~, Mzmt_total, mxymt_array, mzmt_array, Mzmt_bound_total] = ...
    seq_block_MT(TD, npulse, alpha, TR, MT_para_remote, MT_prep, M0_remote);
[~, Mzmimt_total, mxymimt_array, mzmimt_array, Mzmimt_bound_total] = ...
    seq_block_MT(TD, npulse, alpha, TR, MT_para_mi, MT_prep, M0_mi);

%% Plots
figure();
plot(t, Mzmt_total, 'LineWidth', 2)
hold on;
plot(t, Mzmimt_total, 'LineWidth', 2)

legend({'Remote (MT)', 'MI (MT)'});
xlabel('Time (ms)'); ylabel('M_z')
grid on;

%% EPG will not be needed for prelim simulation
npulse=200;
dw=0;
mxymt = zeros(2*npulse-1, npulse, length(ky));
ssmt = zeros(length(ky), 1);

for i = 1:length(ky)
    ky_iter = ky(i);
    phi = RF_phase_cycle(npulse,'balanced'); % alpha, -alpha, alpha, -alpha ...
    % EPGX-BM
    % There is no input to specify RF duration
    [~,fnmt] = EPGX_GRE_MT(d2r(alpha)*ones(npulse,1),phi,b1sqrdtau_iter*ones(npulse,1),...
    tr,T1y,T2_MT,fy,ky_iter,G,'kmax',inf);
    mxymt(:,:,i) = size(fnmt,1)*ifftshift(ifft(ifftshift(fnmt,1),[],1),1);
end
mxymt_array = zeros(length(ky), 1);
for i = 1:length(ky)
    mxymt_array(i) = abs(mxymt(floor(size(mxymt, 1)/2),end,i));
end
figure();
plot(ky, mxymt_array, 'LineWidth', 2)