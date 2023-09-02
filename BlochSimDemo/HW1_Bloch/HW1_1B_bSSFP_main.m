% BMP229, Bloch Equation Simulation
% Balanced SSFP transition to steady state
%
% Modified: 2014/04/14
% created: 2014/03/05 
% Holden H. Wu

%close all; clear all;
clear all;

% set path
addpath('/Users/jameszhang/Documents/MATLAB/PBM229_Project/');

% Debug flags
PLOT_EACHSPIN = 0;
PLOT_SS = 0;


% Sim settings
% ========================================================================
% spin properties
M0 = [0 0 1].'; % equilibrium
BSIM.M = M0;
%BSIM.M = [0 0 0].'; % assume a perfect sat pulse to start
BSIM.T1 = 600;  %ms
BSIM.T2 = 100;  %ms
BSIM.df = 0;    %Hz
% sequence properties
BSIM.TR = 5;    %ms
BSIM.TE = BSIM.TR/2; %ms
BSIM.flipy = 70 *(pi/180); %rad
%   for RF phase cycling
BSIM.RFCYC = 1; % 0:none, 1:SSFP phase cycling, 2:SPGR spoiling
BSIM.RFphi0= 0;
BSIM.RFdph = 2*0.5*pi; 
% sim properties
BSIM.NTR = 200;  %number of TRs to sim
BSIM.dt  = 0.5;  %ms, sim step size in time
BSIM.NstepTR = BSIM.TR/BSIM.dt;  % sim steps per TR
BSIM.NstepTot= BSIM.NTR*BSIM.NstepTR; % total sim steps

% for bSSFP catalyzation
% - no prep, use catmode=0 and NcatTR=0
% - theta/2-TR/2 prep, use catmode=1 and NcatTR=1
% - linear ramp, use catmode=2 and NcatTR=desired
% - cosine ramp, use catmode=3 and NcatTR=desired 
bSSFPcat.catmode  = 2;
bSSFPcat.NcatTR   = 20;
bSSFPcat.NstepTR  = BSIM.TR/BSIM.dt;
if( bSSFPcat.catmode==0 )
    bSSFPcat.NstepTot = 0;
elseif( bSSFPcat.catmode==1 )
    bSSFPcat.NstepTot = bSSFPcat.NcatTR*bSSFPcat.NstepTR/2;
elseif( bSSFPcat.catmode==2 )
    bSSFPcat.NstepTot = bSSFPcat.NcatTR*bSSFPcat.NstepTR;
elseif( bSSFPcat.catmode==3 )
    bSSFPcat.NstepTot = bSSFPcat.NcatTR*bSSFPcat.NstepTR;    
else
    fprintf(1,'bSSFPcat.catmode = %d not supported\n',bSSFPcat.catmode);
    return;
end    

%time = [1:(bSSFP.NstepTot+BSIM.NstepTot)]*BSIM.dt;
time = [1:BSIM.NstepTot]*BSIM.dt;
% ========================================================================


% an ensemble of spins w/range of df
f_vec = [-400:10:400]; %Hz
Nf = numel(f_vec);

% keep track of temporal evolution
Msimf = zeros(3, BSIM.NstepTot + bSSFPcat.NstepTot, Nf);

% Simulate!
for f=1:Nf,
    BSIM.df = f_vec(f);
    
    % bSSFP catalyzation cycles
    if( bSSFPcat.catmode>0 )
        [MsimCat, RF_phiCat] = bSSFPprepfunc(BSIM, 0, bSSFPcat);
        
        % update magnetization history
        BSIM.M = MsimCat(:,end);
        BSIM.RFphi0 = RF_phiCat(end);
    else
        MsimCat = []; RF_phiCat = [];
    end
        
    % Sim TR by TR
    [Msim, RF_phi] = SRsimfunc(BSIM);
    % plot
    if( PLOT_EACHSPIN )
        figure(101); plot( time,Msim(1,:), time,Msim(2,:), time,Msim(3,:) );
        legend('Mx', 'My', 'Mz'); ylim([-1 1]); grid on;
        title(sprintf('GRE with gradient + RF spoiling, spin %d/%d',p,Nphi));
        xlabel('ms'); ylabel('Normalized Magnetization');
        pause;
    end
    
    % keep track of each spin
    Msimf(:, :, f) = cat(2, MsimCat,Msim);
    
    % clear history for each df
    BSIM.RFphi0 = 0;
    BSIM.M = M0;
    %BSIM.M = [0 0 0].'; % assume a perfect sat pulse to start
end

% extract Mxy
Mxy = squeeze( Msimf(1, :, :) + 1i*Msimf(2, :, :) );

% plot freq profile at TE, over all TRs
if( bSSFPcat.catmode==0 )
    plotTEidx = [0:(BSIM.NTR-1)]*BSIM.NstepTR + BSIM.TE/BSIM.dt;    
elseif( bSSFPcat.catmode==1 )
    plotTEidx = bSSFPcat.NstepTot+[0:(BSIM.NTR-1)]*BSIM.NstepTR + BSIM.TE/BSIM.dt;
    plotTEidx = [1, plotTEidx];
elseif( bSSFPcat.catmode==2 )
    plotTEidx = [0:(bSSFPcat.NcatTR+BSIM.NTR-1)]*BSIM.NstepTR + BSIM.TE/BSIM.dt;    
elseif( bSSFPcat.catmode==3 )
    plotTEidx = [0:(bSSFPcat.NcatTR+BSIM.NTR-1)]*BSIM.NstepTR + BSIM.TE/BSIM.dt;        
end
MxyTE = Mxy(plotTEidx, :);

figure; imagesc(f_vec,1:numel(plotTEidx), abs(MxyTE)); colormap gray;
xlabel('Hz'); ylabel('TR #');

% demod RF phase from MxyTE
RF_phi_all = [RF_phiCat; RF_phi];
RF_phi_all = exp(-1i*RF_phi_all) * ones(1,Nf);
MxyTE = MxyTE .* RF_phi_all;

% plot signal at TE, over all TRs, for specific df
[Y, fidx1] = min( abs(f_vec-0) ); % locate the pass-band idx
[Y, fidx2] = min( abs(f_vec-1000/2/BSIM.TR) ); % locate the stop-band idx
fidx3 = round( mean([fidx1, fidx2]) );

figure; 
subplot(2,1,1);
plot( 1:numel(plotTEidx),abs(MxyTE(:,fidx1)), 1:numel(plotTEidx),abs(MxyTE(:,fidx3)), 1:numel(plotTEidx),abs(MxyTE(:,fidx2)) );
hold on; line(1+[bSSFPcat.NcatTR bSSFPcat.NcatTR], [0 1], 'Color','k');
xlim([1 numel(plotTEidx)]); legend(num2str(f_vec([fidx1 fidx3 fidx2]).')); ylim([0 1]);
xlabel('TR #'); ylabel('Mxy (normalized)');
title('evolution of pass-band and stop-band over time');
subplot(2,1,2);
plot( 1:numel(plotTEidx),angle(MxyTE(:,fidx1)), 1:numel(plotTEidx),angle(MxyTE(:,fidx3)), 1:numel(plotTEidx),angle(MxyTE(:,fidx2)) );
hold on; line(1+[bSSFPcat.NcatTR bSSFPcat.NcatTR], [-pi pi], 'Color','k');
xlim([1 numel(plotTEidx)]); legend(num2str(f_vec([fidx1 fidx3 fidx2]).')); ylim([-pi pi]);
xlabel('TR #'); ylabel('Mxy (phase)');
title('evolution of pass-band and stop-band over time');

if( PLOT_SS )
    % plot steady-state freq profile at TE (for last TR)
    MxyTElast = squeeze( MxyTE(end, :) );
    % % demod RF phase
    % MxyTElast = MxyTElast .* exp(-1i*RF_phi(end));
    figure;
    subplot(2,1,1);
    plot(f_vec, abs(MxyTElast));
    title('Balanced SSFP frequency profile: mag');
    xlabel('Hz'); ylabel('Normalized |Mxy|');
    subplot(2,1,2);
    plot( f_vec, angle( MxyTElast ) );
    title('Balanced SSFP frequency profile: phase');
    xlabel('Hz'); ylabel('Phase of Mxy');
end

