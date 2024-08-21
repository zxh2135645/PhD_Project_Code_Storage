% BMP229, Bloch Equation Simulation
% Steady-state signal for bSSFP, SSFP-FID/Echo
%
% Modified: 2014/04/14
% Created: 2014/03/17
% Holden Wu

%close all; clear all;
clear all;

% set path
addpath('/Users/jameszhang/Documents/MATLAB/PBM229_Project/');

% plot settings
PLOT_NOIR = 1;

% Tissue params
VARY_T1 = 1;
if( ~VARY_T1 )
  % Fix T1, vary T2
  T1vec = [1000]; %ms
  T2vec = [100, 200, 500, 1000]; %ms
else
  % Fix T2, vary T1
  T1vec = [100, 200, 500, 1000]; %ms
  T2vec = [40]; %ms
end

% Sequence sim params
% Note: a coarser grid is used for SSFP-FID/Echo to reduce run time
MODE = 1;
if( 1==MODE )
  % bSSFP sequence params
  TRms = 5; TEms = TRms/2;
  % flip angle
  RFthdeg = [0:1:180];
  %RFset = [11 51 91]; % for now, hard-coded indices to select RFth = 10, 50, 90 deg
  % off-resonance phase accrual
  dfvec = [-400:20:400]; %Hz
elseif( 2==MODE )
  % SSFP-FID
  TRms = 10; TEms = 2; FISP = 'FID';
  % flip angle
  RFthdeg = [0:10:180];
  %RFset = [2 6 10];
  % off-resonance phase accrual
  dfvec = [-400:50:400]; %Hz
elseif( 3==MODE )    
  % SSFP-Echo
  TRms = 10; TEms = 8; FISP = 'ECHO';
  % flip angle
  RFthdeg = [0:10:180];
  %RFset = [2 6 10];
  % off-resonance phase accrual
  dfvec = [-400:50:400]; %Hz
else
  fprintf(1,'MODE = %d not supported\n', MODE);
  return;
end

% flip angle
%RFthdeg = [0:1:180];
RFthrad = RFthdeg * pi/180;
RFthplot = 60;
[y, RFI] = min(abs(RFthdeg-RFthplot)); % for plotting
% off-resonance phase accrual
%dfvec = [-400:20:400]; %Hz
dfplot = 1000/TRms/2;
[y, dfI] = min(abs(dfvec-dfplot)); % for plotting

Mxyss = zeros(numel(T1vec), numel(T2vec), numel(RFthdeg), numel(dfvec));
for t1idx=1:numel(T1vec),
  T1ms = T1vec(t1idx);
  for t2idx=1:numel(T2vec),
    T2ms = T2vec(t2idx);
    %-------------------------------------------------------------------------
    for rfidx = 1:numel(RFthrad),
      RFth = RFthrad(rfidx);
      for dfidx = 1:numel(dfvec),
        df = dfvec(dfidx);
        if( 1==MODE )
          [Msig0, Mss0] = sssignal(RFth,T1ms,T2ms,TEms,TRms,df);
        else  
          [Msig0, Mss0] = gresignal(RFth,T1ms,T2ms,TEms,TRms,df, FISP);
        end  
        Mxyss(t1idx,t2idx,rfidx,dfidx) = Msig0;
      end  % df loop  
    end % FA loop
    %-------------------------------------------------------------------------
  end % T2 loop
end % T1 loop

figure(100*MODE + 10*VARY_T1 + 1);
% magnitude
subplot(2,1,1);
myplot=plot(dfvec, squeeze(abs(Mxyss(:,:,RFI,:)))); 
if( 1==numel(T1vec) )
  mylgnd=legend(num2str(T2vec.'/T1ms));
else
  mylgnd=legend(num2str(T2ms./T1vec.'));
end  
xlabel('off-resonance (Hz)'); ylabel('Mxy,ss (normalized mag)');
title(sprintf('Different T2/T1 ratios, \\theta = %d deg', RFthdeg(RFI)));
xlim([min(dfvec), max(dfvec)]); ylim([0 0.5]);
% modify plot appearance
myplotmod(gcf, gca, myplot, mylgnd, PLOT_NOIR);
% phase
subplot(2,1,2);
myplot2=plot(dfvec, squeeze(angle(Mxyss(:,:,RFI,:))));
if( 1==numel(T1vec) )
  mylgnd=legend(num2str(T2vec.'/T1ms));
else
  mylgnd=legend(num2str(T2ms./T1vec.'));
end  
xlabel('off-resonance (Hz)'); ylabel('Mxy,ss (phase)');
xlim([min(dfvec), max(dfvec)]);
% modify plot appearance
myplotmod(gcf, gca, myplot2, mylgnd, PLOT_NOIR);

% if( 1==numel(T1vec) )
%   profiles = squeeze( Mxyss(1,end,RFset,:) ).';
% else
%   profiles = squeeze( Mxyss(end,1,RFset,:) ).';
% end  
% figure(100*MODE + 10*VARY_T1 + 2);
% myplot=plot(dfvec, abs(profiles)); mylgnd=legend(num2str(RFthdeg(RFset).'));
% xlabel('off-resonance (Hz)'); ylabel('Mxy,ss (normalized mag)');
% title(sprintf('T2/T1 = %.2f, different thetas', T2vec(end)/T1vec(end)));
% xlim([min(dfvec), max(dfvec)]);
% % modify plot appearance
% myplotmod(gcf, gca, myplot, mylgnd, PLOT_NOIR);

figure(100*MODE + 10*VARY_T1 +3); 
% dfI = 40
myplot=plot(RFthdeg, squeeze(abs(Mxyss(:,:,:,dfI))));
if( 1==numel(T1vec) )
  mylgnd=legend(num2str(T2vec.'/T1ms));
else
  mylgnd=legend(num2str(T2ms./T1vec.'));
end  
xlabel('RF flip angle \theta'); ylabel('Mxy,ss (normalized)');
title(sprintf('Different T2/T1 ratios, df = %d Hz', dfvec(dfI)));
xlim([min(RFthdeg), max(RFthdeg)]); ylim([0 0.5]);
% modify plot appearance
myplotmod(gcf, gca, myplot, mylgnd, PLOT_NOIR);

