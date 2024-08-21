% BMP229
% Balanced SSFP frequency response
%
% 2014/03/17
% Holden Wu

close all; clear all;

T1ms = 1000; %T2ms = 1000;

TRms = 5; TEms = TRms/2;

% flip angle
RFthdeg = [0:1:180];
RFthrad = RFthdeg * pi/180;
RFthplot = 60;
[y, RFI] = min(abs(RFthdeg-RFthplot));
% off-resonance phase accrual
dfvec = [-400:20:400]; %Hz
dfplot = 1000/TRms/2;
[y, dfI] = min(abs(dfvec-dfplot));

T2vec = [100, 200, 500, 1000];
Mxyss = zeros(numel(T2vec), numel(RFthdeg), numel(dfvec));
for t2idx=1:numel(T2vec),
T2ms = T2vec(t2idx);
%-------------------------------------------------------------------------
% Signal as a function of flip angle
Mss = zeros(3, numel(RFthrad));
%beta = betarad;
for rfidx = 1:numel(RFthrad),
  RFth = RFthrad(rfidx);
  for dfidx = 1:numel(dfvec),
    df = dfvec(dfidx);
    [Msig0, Mss0] = sssignal(RFth,T1ms,T2ms,TEms,TRms,df);
    Mxyss(t2idx,rfidx,dfidx) = Msig0;
  end  
end
%-------------------------------------------------------------------------
end

figure(101);
% magnitude
subplot(2,1,1);
myplot=plot(dfvec, squeeze(abs(Mxyss(:,RFI,:)))); mylgnd=legend(num2str(T2vec.'/T1ms));
xlabel('off-resonance (Hz)'); ylabel('Mxy,ss (normalized mag)');
title(sprintf('Different T2/T1 ratios, \\theta = %d deg', RFthdeg(RFI)));
xlim([min(dfvec), max(dfvec)]);
% modify plot appearance
set(gcf,'Color','k');
set(gca,'Color','k','XColor','w','YColor','w');
set(myplot,'LineWidth',1.5);
set(mylgnd,'color','none','TextColor','w');
set(get(gca,'Title'),'Color','w','FontSize',16);
set(get(gca,'Xlabel'),'FontSize',16);
set(get(gca,'Ylabel'),'FontSize',16);
% phase
subplot(2,1,2);
myplot2=plot(dfvec, squeeze(angle(Mxyss(:,RFI,:)))); mylgnd=legend(num2str(T2vec.'/T1ms));
xlabel('off-resonance (Hz)'); ylabel('Mxy,ss (phase)');
xlim([min(dfvec), max(dfvec)]);
% modify plot appearance
set(gcf,'Color','k');
set(gca,'Color','k','XColor','w','YColor','w');
set(myplot2,'LineWidth',1.5);
set(mylgnd,'color','none','TextColor','w');
set(get(gca,'Title'),'Color','w','FontSize',16);
set(get(gca,'Xlabel'),'FontSize',16);
set(get(gca,'Ylabel'),'FontSize',16);

RFset = [11 51 91]; % for now, hard-coded indices to select RFth = 10, 50, 90 deg
profiles = squeeze( Mxyss(4,RFset,:) ).';
figure(102);
myplot=plot(dfvec, abs(profiles)); mylgnd=legend(num2str(RFthdeg(RFset).'));
xlabel('off-resonance (Hz)'); ylabel('Mxy,ss (normalized mag)');
title(sprintf('T2/T1 = %.2f, different thetas', T2vec(4)/T1ms));
xlim([min(dfvec), max(dfvec)]);
% modify plot appearance
set(gcf,'Color','k');
set(gca,'Color','k','XColor','w','YColor','w');
set(myplot,'LineWidth',1.5);
set(mylgnd,'color','none','TextColor','w');
set(get(gca,'Title'),'Color','w','FontSize',16);
set(get(gca,'Xlabel'),'FontSize',16);
set(get(gca,'Ylabel'),'FontSize',16);

figure(103); 
myplot=plot(RFthdeg, abs(Mxyss(:,:,dfI))); mylgnd=legend(num2str(T2vec.'/T1ms));
xlabel('RF flip angle \theta'); ylabel('Mxy,ss (normalized)');
title(sprintf('Different T2/T1 ratios, df = %d Hz', dfvec(dfI)));
xlim([min(RFthdeg), max(RFthdeg)]);
% modify plot appearance
set(gcf,'Color','k');
set(gca,'Color','k','XColor','w','YColor','w');
set(myplot,'LineWidth',1.5);
set(mylgnd,'color','none','TextColor','w');
set(get(gca,'Title'),'Color','w','FontSize',16);
set(get(gca,'Xlabel'),'FontSize',16);
set(get(gca,'Ylabel'),'FontSize',16);

