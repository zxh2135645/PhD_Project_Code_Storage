%% Reproduce Gloor 2008
clear all;
close all;
addpath(genpath('lib'));
addpath(genpath('EPGX-src'));
%% General Parameters
TR = 2.92;
T1 = 585;
T2 = 81; 
alpha = 0:5:120;
E1 = exp(-TR / T1);
E2 = exp(-TR / T2);
rho = 1;
M = rho * sin(d2r(alpha)) * (1 - E1) .* sqrt(E2) ./ (1 - (E1 - E2) .* cos(d2r(alpha)) - E1 .* E2);
%% MT parameters initiation
%%% Relaxation parameters: exchange
T1x = [T1 1000]; % ms
T2x = [T2 12e-3]; % ms
T2_MT = T2;
kx = 4.45e-3; % ms
fx = 0.157;
%%% RF saturation factor for MT
[ff,G] = SuperLorentzian(T2x(2)*1e-3);% us super-Lorentzian absorption lineshape
% [ff_gauss,G_gauss] = GaussianLineShape(T2x(2)*1e-3);

gam = 267.5221 *1e-3; %< rad /ms /uT
trf = 0.230; % ms
b1 = d2r(alpha)./(trf.*gam); % uT
b1sqrdtau = b1.^2.*trf;
b1sqrdtau2 = 0;
% b1sqrdtau = zeros(1, length(b1sqrdtau));
% figure();
% plot(ff, G, 'LineWidth', 2)
% title('Absorption Line Shape');
% xlabel('Offset (Hz)'); ylabel('us')
% grid on;
% hold on;
% plot(ff_gauss, G_gauss, 'LineWidth', 2);
% legend({'Super-Lorenztian', 'Gaussian'});
idx = find(ff == min(abs(ff)));
G = G(idx);


% Gloor analytical solution
% analytic form
mxyGloor = zeros(length(alpha), 1);
mxyGloor_noMT = zeros(length(alpha), 1);
for i = 1:length(alpha)
    al = alpha(i);
    mxyGloor(i) = ssSSFP_Gloor(d2r(al),b1sqrdtau,TR,T1x,T2_MT,fx,kx,G);
    mxyGloor_noMT(i) = ssSSFP_Gloor(d2r(al),b1sqrdtau2, TR, T1x,T2_MT,fx,kx,G);
end
%% Plots
figure();
plot(alpha, M, 'LineWidth', 2);
grid on;
xlabel('alpha (degree)'); ylabel('Magnetization (M_{xy})');
hold on;
title('Theoretical Magnetization From Myocardium and Water');
plot(alpha, mxyGloor, 'LineWidth', 2)
plot(alpha, mxyGloor_noMT, 'LineWidth', 2)
MTR = (M' - mxyGloor) ./ M';
yyaxis right;
plot(alpha, MTR, 'LineWidth', 2)
legend({'M', 'MT', 'no MT', 'MTR'});
%% Fig 2C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% General Parameters
TR = 2.72:0.1:4.72;
T1 = 585;
T2 = 81; 
alpha = 35;
E1 = exp(-TR / T1);
E2 = exp(-TR / T2);
rho = 1;
M = rho * sin(d2r(alpha)) * (1 - E1) .* sqrt(E2) ./ (1 - (E1 - E2) .* cos(d2r(alpha)) - E1 .* E2);

%% MT parameters initiation
%%% Relaxation parameters: exchange
T1x = [T1 1000]; % ms
T2x = [T2 12e-3]; % ms
T2_MT = T2;
kx = 4.45e-3; % ms
fx = 0.157;
%%% RF saturation factor for MT
[ff,G] = SuperLorentzian(T2x(2)*1e-3);% us super-Lorentzian absorption lineshape
% [ff_gauss,G_gauss] = GaussianLineShape(T2x(2)*1e-3);

gam = 267.5221 *1e-3; %< rad /ms /uT
trf = linspace(0.030, 2.030, length(TR)); % ms
b1 = d2r(alpha)./(trf.*gam); % uT
b1sqrdtau = b1.^2.*trf;
b1sqrdtau2 = zeros(length(TR), 1);
% b1sqrdtau = zeros(1, length(b1sqrdtau));
% figure();
% plot(ff, G, 'LineWidth', 2)
% title('Absorption Line Shape');
% xlabel('Offset (Hz)'); ylabel('us')
% grid on;
% hold on;
% plot(ff_gauss, G_gauss, 'LineWidth', 2);
% legend({'Super-Lorenztian', 'Gaussian'});
idx = find(ff == min(abs(ff)));
G = G(idx);


% Gloor analytical solution
% analytic form
mxyGloor = zeros(length(TR), 1);
mxyGloor_noMT = zeros(length(TR), 1);
for i = 1:length(TR)
    tr = TR(i);
    b1sqrdtau_iter = b1sqrdtau(i);
    b1sqrdtau_iter2 = b1sqrdtau2(i);
    mxyGloor(i) = ssSSFP_Gloor(d2r(alpha),b1sqrdtau_iter,tr,T1x,T2_MT,fx,kx,G);
    mxyGloor_noMT(i) = ssSSFP_Gloor(d2r(alpha),b1sqrdtau_iter2, tr, T1x,T2_MT,fx,kx,G);
end

%% Idealized Upper Bound
T1_id = 1/(1/T1 + kx);
rho_id = 1/T1 * T1_id;
E1_id = exp(-TR / T1_id);
M_up = rho_id * sin(d2r(alpha)) * (1 - E1_id) .* sqrt(E2) ./ (1 - (E1_id - E2) .* cos(d2r(alpha)) - E1_id .* E2);

%% Plots
figure();
plot(TR, M, 'LineWidth', 2);
grid on;
xlabel('TR (ms)'); ylabel('Magnetization (M_{xy})');
hold on;
title('Theoretical Magnetization From Myocardium and Water');
plot(TR, mxyGloor, 'LineWidth', 2)
plot(TR, mxyGloor_noMT, 'LineWidth', 2)
plot(TR, M_up, 'LineWidth', 2)
MTR = (M' - mxyGloor) ./ M';
yyaxis right;
plot(TR, MTR, 'LineWidth', 2)
legend({'M', 'MT', 'no MT', 'LowerBound','MTR'});

%% Comments
% The simulation seems matched for most of the points, but there is
% inconsistency when RF duration -> 0 or flip angle close to 0 (more
% discrepancy making making MTR curve not matching)