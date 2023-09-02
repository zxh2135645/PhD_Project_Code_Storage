%% Simulation of two-pool 4% Agar
addpath(genpath('lib'));
addpath(genpath('EPGX-src'));

%%% Sequences
% TR = hdr_flow_fa_lowBW{1}{1}.EchoTime * 2; % ms
% TRmin = hdr_flow_fa_highBW{1}{1}.EchoTime * 2; % ms
TR = 9.5;
TRmin = 4;
alpha = 60;

agar_t1 = 1000;
agar_t2 = 43;
%%% Relaxation parameters: single pool, copy MT model
T1=agar_t1;
T2=agar_t2;

%%% Relaxation parameters: MT
% T1_MT = [agar_t1 agar_t1]; % ms
% f_MT = 0.011; % Unitless
% k_MT = 1.8e-3; % ms^-1
% T2_MT = agar_t2; % ms
T1_MT = [agar_t1 agar_t1]; % ms
%f_MT = 0.07; % Unitless
%k_MT = 4.3e-3; % ms^-1
f_MT = 0.011; % Unitless
k_MT = 1.8e-3; % ms^-1
T2_MT = agar_t2; % ms

%%% Relaxation parameters: exchange
T1x = [agar_t1 agar_t1]; % ms
T2x = [agar_t2 12.9e-3]; % ms
%T2x = [agar_t2 10]
kx = 1.8e-3; % ms
%fx = 0.2;
fx = 0.011;

%%% RF saturation factor for MT
[ff,G] = SuperLorentzian(T2x(2)*1e-3);% us super-Lorentzian absorption lineshape
[ff_gauss,G_gauss] = GaussianLineShape(T2x(2)*1e-3);
b1 = 13; % uT
gam = 267.5221 *1e-3; %< rad /ms /uT
trf = d2r(alpha)/(gam*b1);% ms
b1sqrdtau = b1^2*trf;

figure();
plot(ff, G, 'LineWidth', 2)
title('Absorption Line Shape');
xlabel('Offset (Hz)'); ylabel('us')
grid on;
hold on;
plot(ff_gauss, G_gauss, 'LineWidth', 2);
legend({'Super-Lorenztian', 'Gaussian'});

print -dpng -r300 bin/LineShapes.png
idx = find(ff == min(abs(ff)));
G = G_gauss(idx); % G(0) = ? us
%% bSSFP with same parameters as above

%%% Do four cases: EPG, EPG-X(MT), EPG-X(BM) with delta=0 and 0.1ppm
delta = 1e-6 * 3 * 42.6e3;

% EPG simulations
npulse=floor(5*T1/TR);
%npulse = 779;
phi = RF_phase_cycle(npulse,'balanced'); % alpha, -alpha, alpha, -alpha ...

% single pool
[~,fn0] = EPG_GRE(d2r(alpha)*ones(npulse,1),phi,TR,T1,T2,'kmax',inf);
mxys = size(fn0,1)*ifftshift(ifft(ifftshift(fn0,1),[],1),1);

% MT
[~,fnmt] = EPGX_GRE_MT(d2r(alpha)*ones(npulse,1),phi,b1sqrdtau*ones(npulse,1),...
    TR,T1_MT,T2_MT,f_MT,k_MT,G,'kmax',inf);
mxymt = size(fnmt,1)*ifftshift(ifft(ifftshift(fnmt,1),[],1),1);

[~,fnmt2] = EPGX_GRE_MT(d2r(alpha)*ones(npulse,1),phi,b1sqrdtau*ones(npulse,1),...
    TRmin,T1_MT,T2_MT,f_MT,k_MT,G,'kmax',inf);
mxymt2 = size(fnmt2,1)*ifftshift(ifft(ifftshift(fnmt2,1),[],1),1);

% BM
[~,fnx] = EPGX_GRE_BM(d2r(alpha)*ones(npulse,1),phi,TR,T1x,T2x,fx,kx,'kmax',inf);
mxyx = ifftshift(fft(ifftshift(fnx,1),[],1),1);
mxyx1 = sum(mxyx,3); % sum compartments

% delta = 0.1ppm
[~,fnx] = EPGX_GRE_BM(d2r(alpha)*ones(npulse,1),phi,TR,T1x,T2x,fx,kx,'kmax',inf,'delta',delta);
mxyx = ifftshift(fft(ifftshift(fnx,1),[],1),1);
mxyx2 = sum(mxyx,3); % sum compartments

% This part is analytical result?
dw=0;% on resonance 
ss0 = ssSSFP(d2r(alpha),TR,dw,T1,T2);
ssx1 = ssSSFP_BM(d2r(alpha),TR,dw,T1x,T2x,fx,kx,0);
ssx2 = ssSSFP_BM(d2r(alpha),TR,dw,T1x,T2x,fx,kx,delta);
ssmt= ssSSFP_MT(d2r(alpha),b1sqrdtau,TR,dw,T1_MT,T2_MT,f_MT,k_MT,G);


%%% steady state
phiTR = linspace(-pi/TR,pi/TR,npulse);
ssfpx1SS = zeros(size(phiTR));
ssfpx2SS = zeros(size(phiTR));
ssfpmtSS= zeros(size(phiTR));
ssfpmtSS2= zeros(size(phiTR));
ssfpss = zeros(size(phiTR));

for ii=1:npulse
    ssfpss(ii) = abs(ssSSFP(d2r(alpha),TR,phiTR(ii),T1,T2));
    ssfpx1SS(ii) = abs(ssSSFP_BM(d2r(alpha),TR,phiTR(ii),T1x,T2x,fx,kx,0));
    ssfpx2SS(ii) = abs(ssSSFP_BM(d2r(alpha),TR,phiTR(ii),T1x,T2x,fx,kx,delta));
    ssfpmtSS(ii) = abs(ssSSFP_MT(d2r(alpha),b1sqrdtau,TR,phiTR(ii),T1_MT,T2_MT,f_MT,k_MT,G));
    ssfpmtSS2(ii) = abs(ssSSFP_MT(d2r(alpha),b1sqrdtau,TRmin,phiTR(ii),T1_MT,T2_MT,f_MT,k_MT,G));
end
% analytic form
mxyGloor = ssSSFP_Gloor(d2r(alpha),b1sqrdtau, TR, T1_MT,T2_MT,f_MT,k_MT,G);
%% Figure just with plots
figure();
clf % Clear current figure window

subplot(221)
pp=plot(linspace(-pi,pi,npulse),ssfpss,'-');
pp.LineWidth = 2;
hold
plot(linspace(-pi,pi,size(mxys,1)),squeeze(abs(mxys(:,end))),'--','linewidth',2)
grid on
xlim([-pi pi])
ylim([0 0.2])
ylabel('Signal / M_0')
xlabel('\psi / rad')
legend( 'Direct Steady-State','EPG','location','south')
title('Single compartment')
set(gca,'fontsize',14)

subplot(222)
pp=plot(linspace(-pi,pi,npulse),ssfpmtSS,'-');
pp.LineWidth = 2;
hold
pp=plot(linspace(-pi,pi,npulse),ssfpss,'-');
% pp=plot(linspace(-pi,pi,npulse),ssfpmtSS2,'.');
plot(linspace(-pi,pi,size(mxymt,1)),squeeze(abs(mxymt(:,end))),'--','linewidth',2)
plot(linspace(-pi,pi,size(mxymt2,1)),squeeze(abs(mxymt2(:,end))),'-.','linewidth',2)
plot(0,mxyGloor,'*r','MarkerSize',12,'Linewidth',2)
grid on;
xlim([-pi pi])
ylim([0 0.2])
ylabel('Signal / M_0')
xlabel('\psi / rad')
legend( 'Direct Steady-State', 'Single Pool','EPG-X (MT)','EPG-X (MT) minTR','Analytic (Gloor 2008)','location','north')
title('White matter MT model')
set(gca,'fontsize',14)


subplot(223)
pp=plot(linspace(-pi,pi,npulse),ssfpx1SS,'-');
pp.LineWidth = 2;
hold
plot(linspace(-pi,pi,size(mxyx1,1)),squeeze(abs(mxyx1(:,end))),'--','linewidth',2)
grid on
xlim([-pi pi])
ylim([0 0.2])
ylabel('Signal / M_0')
xlabel('\psi / rad')
legend( 'Direct Steady-State','EPG-X (BM)','location','south')
title('Myelin water exchange model, \delta_b=0')
set(gca,'fontsize',14)

subplot(224)
pp=plot(linspace(-pi,pi,npulse),ssfpx2SS,'-');
pp.LineWidth = 2;
hold
plot(linspace(-pi,pi,size(mxyx1,1)),squeeze(abs(mxyx2(:,end))),'--','linewidth',2)
plot(linspace(-pi,pi,size(mxyx1,1)),squeeze(abs(mxyx(:,end,1))),'--','linewidth',2)
plot(linspace(-pi,pi,size(mxyx1,1)),squeeze(abs(mxyx(:,end,2))),'--','linewidth',2)
grid on
xlim([-pi pi])
ylim([0 0.2])
ylabel('Signal / M_0')
xlabel('\psi / rad')
legend( 'Direct Steady-State','EPG-X (BM)','location','south')
title('Myelin water exchange model, \delta_b=0.1ppm')
set(gca,'fontsize',14)

setpospap([300 100 700 550])

gg=get(gcf,'children');
ii=[8 6 4 2];
lbls={'(a)','(b)','(c)','(d)'};
for jj=1:4
    axes(gg(ii(jj)))
    text(-4,-0.02,lbls{jj},'fontsize',20,'fontweight','bold')
end
%
print -dpng -r300 bin/Test1_fig3.png
%% For my own purpose
figure('Position', get(0, 'Screensize'));
subplot(121)
pp=plot(linspace(-pi,pi,npulse),ssfpmtSS,'-');
pp.LineWidth = 2;
hold
pp=plot(linspace(-pi,pi,npulse),ssfpss,'-');
% pp=plot(linspace(-pi,pi,npulse),ssfpmtSS2,'.');
plot(linspace(-pi,pi,size(mxymt,1)),squeeze(abs(mxymt(:,end))),'--','linewidth',2)
plot(linspace(-pi,pi,size(mxymt2,1)),squeeze(abs(mxymt2(:,end))),'-.','linewidth',2)
plot(0,mxyGloor,'*r','MarkerSize',12,'Linewidth',2)
grid on;
xlim([-pi pi])
ylim([0 0.2])
ylabel('Signal / M_0')
xlabel('\psi / rad')
legend( 'Direct Steady-State', 'Single Pool','EPG-X (MT)','EPG-X (MT) minTR','Analytic (Gloor 2008)','location','north')
title('4% Agar MT model')
set(gca,'fontsize',14)

%% Heart
TR = 9.5;
TRmin = 4;
heart_t1 = 1175;
heart_t2 = 54.4;
%%% Relaxation parameters: single pool, copy MT model
T1=heart_t1;
T2=heart_t2;

%%% Relaxation parameters: MT

T1_MT = [heart_t1 1000]; % ms
f_MT = 0.07; % Unitless
k_MT = 4.1e-3; % ms^-1
T2_MT = heart_t2; % ms

%%% Relaxation parameters: exchange
T1x = [heart_t1 1000]; % ms
T2x = [heart_t2 8.5e-3]; % ms
kx = 4.1e-3; % ms
fx = 0.07;

[ff,G] = SuperLorentzian(T2x(2)*1e-3);% us super-Lorentzian absorption lineshape
G = G(idx);
%% bSSFP with same parameters as above

%%% Do four cases: EPG, EPG-X(MT), EPG-X(BM) with delta=0 and 0.1ppm
delta = 1e-6 * 3 * 42.6e3;

% EPG simulations
npulse=floor(5*T1/TR);
%npulse = 779;
phi = RF_phase_cycle(npulse,'balanced'); % alpha, -alpha, alpha, -alpha ...

% single pool
[~,fn0] = EPG_GRE(d2r(alpha)*ones(npulse,1),phi,TR,T1,T2,'kmax',inf);
mxys = size(fn0,1)*ifftshift(ifft(ifftshift(fn0,1),[],1),1);

% MT
[~,fnmt] = EPGX_GRE_MT(d2r(alpha)*ones(npulse,1),phi,b1sqrdtau*ones(npulse,1),...
    TR,T1_MT,T2_MT,f_MT,k_MT,G,'kmax',inf);
mxymt = size(fnmt,1)*ifftshift(ifft(ifftshift(fnmt,1),[],1),1);

[~,fnmt2] = EPGX_GRE_MT(d2r(alpha)*ones(npulse,1),phi,b1sqrdtau*ones(npulse,1),...
    TRmin,T1_MT,T2_MT,f_MT,k_MT,G,'kmax',inf);
mxymt2 = size(fnmt2,1)*ifftshift(ifft(ifftshift(fnmt2,1),[],1),1);

% BM
[~,fnx] = EPGX_GRE_BM(d2r(alpha)*ones(npulse,1),phi,TR,T1x,T2x,fx,kx,'kmax',inf);
mxyx = ifftshift(fft(ifftshift(fnx,1),[],1),1);
mxyx1 = sum(mxyx,3); % sum compartments

% delta = 0.1ppm
[~,fnx] = EPGX_GRE_BM(d2r(alpha)*ones(npulse,1),phi,TR,T1x,T2x,fx,kx,'kmax',inf,'delta',delta);
mxyx = ifftshift(fft(ifftshift(fnx,1),[],1),1);
mxyx2 = sum(mxyx,3); % sum compartments

% This part is analytical result?
dw=0;% on resonance 
ss0 = ssSSFP(d2r(alpha),TR,dw,T1,T2);
ssx1 = ssSSFP_BM(d2r(alpha),TR,dw,T1x,T2x,fx,kx,0);
ssx2 = ssSSFP_BM(d2r(alpha),TR,dw,T1x,T2x,fx,kx,delta);
ssmt= ssSSFP_MT(d2r(alpha),b1sqrdtau,TR,dw,T1_MT,T2_MT,f_MT,k_MT,G);


%%% steady state
phiTR = linspace(-pi/TR,pi/TR,npulse);
ssfpx1SS = zeros(size(phiTR));
ssfpx2SS = zeros(size(phiTR));
ssfpmtSS= zeros(size(phiTR));
ssfpmtSS2= zeros(size(phiTR));
ssfpss = zeros(size(phiTR));

for ii=1:npulse
    ssfpss(ii) = abs(ssSSFP(d2r(alpha),TR,phiTR(ii),T1,T2));
    ssfpx1SS(ii) = abs(ssSSFP_BM(d2r(alpha),TR,phiTR(ii),T1x,T2x,fx,kx,0));
    ssfpx2SS(ii) = abs(ssSSFP_BM(d2r(alpha),TR,phiTR(ii),T1x,T2x,fx,kx,delta));
    ssfpmtSS(ii) = abs(ssSSFP_MT(d2r(alpha),b1sqrdtau,TR,phiTR(ii),T1_MT,T2_MT,f_MT,k_MT,G));
    ssfpmtSS2(ii) = abs(ssSSFP_MT(d2r(alpha),b1sqrdtau,TRmin,phiTR(ii),T1_MT,T2_MT,f_MT,k_MT,G));
end

% analytic form
mxyGloor = ssSSFP_Gloor(d2r(alpha),b1sqrdtau, TR, T1_MT,T2_MT,f_MT,k_MT,G);
subplot(122)
pp=plot(linspace(-pi,pi,npulse),ssfpmtSS,'-');
pp.LineWidth = 2;
hold
pp=plot(linspace(-pi,pi,npulse),ssfpss,'-');
% pp=plot(linspace(-pi,pi,npulse),ssfpmtSS2,'.');
plot(linspace(-pi,pi,size(mxymt,1)),squeeze(abs(mxymt(:,end))),'--','linewidth',2)
plot(linspace(-pi,pi,size(mxymt2,1)),squeeze(abs(mxymt2(:,end))),'-.','linewidth',2)
plot(0,mxyGloor,'*r','MarkerSize',12,'Linewidth',2)
grid on;
xlim([-pi pi])
ylim([0 0.2])
ylabel('Signal / M_0')
xlabel('\psi / rad')
legend( 'Direct Steady-State', 'Single Pool','EPG-X (MT)','EPG-X (MT) minTR','Analytic (Gloor 2008)','location','north')
title('Porcine Heart MT model')
set(gca,'fontsize',14)

print -dpng -r300 bin/Sim_AgarHeart_FA60_2TRs.png