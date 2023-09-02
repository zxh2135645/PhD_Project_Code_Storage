function [MxyTE, BSIM, plotTEidx, bSSFPcat, Msimf, bSSFPend] = bSSFP_engine(M0, T1, T2, df, TR, FA, NTR, dt, f_vec, PLOT_EACHSPIN, PLOT_SS, varargin)
% RAMPDOWN has an effect in Mxy in last time point
% extract if RAMP_DOWN is input
if isempty(varargin)
    RAMP_DOWN = 0;
else
    for ii=1:length(varargin)
        if varargin{ii} == 1
            RAMP_DOWN = 1;
        elseif varargin{ii} == 0
            RAMP_DOWN = 0;
        end
    end
end


BSIM.M = M0;
BSIM.T1 = T1;  %ms
BSIM.T2 = T2;  %ms
BSIM.df = df;    %Hz

% sequence properties
BSIM.TR = TR;    %ms
BSIM.TE = BSIM.TR/2; %ms
BSIM.flipy = FA *(pi/180); %rad

%   for RF phase cycling
BSIM.RFCYC = 1; % 0:none, 1:SSFP phase cycling, 2:SPGR spoiling
BSIM.RFphi0= 0;
BSIM.RFdph = 2*0.5*pi;

% sim properties
BSIM.NTR = NTR;  %number of TRs to sim
BSIM.dt  = dt;  %ms, sim step size in time
BSIM.NstepTR = int32(BSIM.TR/BSIM.dt);  % sim steps per TR
BSIM.NstepTot= int32(BSIM.NTR*BSIM.NstepTR); % total sim steps

% for bSSFP catalyzation
% - no prep, use catmode=0 and NcatTR=0
% - theta/2-TR/2 prep, use catmode=1 and NcatTR=1
% - linear ramp, use catmode=2 and NcatTR=desired
% - cosine ramp, use catmode=3 and NcatTR=desired 
bSSFPcat.catmode  = 2;
bSSFPcat.NcatTR   = 5; % for T1 MOLLI, 5 linear ramp-up
bSSFPcat.NstepTR  = int32(BSIM.TR/BSIM.dt);
if( bSSFPcat.catmode==0 )
    bSSFPcat.NstepTot = 0;
elseif( bSSFPcat.catmode==1 )
    bSSFPcat.NstepTot = int32(bSSFPcat.NcatTR*bSSFPcat.NstepTR/2);
elseif( bSSFPcat.catmode==2 )
    bSSFPcat.NstepTot = int32(bSSFPcat.NcatTR*bSSFPcat.NstepTR);
elseif( bSSFPcat.catmode==3 )
    bSSFPcat.NstepTot = int32(bSSFPcat.NcatTR*bSSFPcat.NstepTR);    
else
    fprintf(1,'bSSFPcat.catmode = %d not supported\n',bSSFPcat.catmode);
    return;
end   

% Ramp-down at the end
if RAMP_DOWN == 0
    bSSFPend.RAMP_DOWN = 0;
    bSSFPend.NendTR = 0; 
elseif RAMP_DOWN == 1
    bSSFPend.RAMP_DOWN = 1;
    bSSFPend.NendTR = 1; % Not sure if ramp-down is TR or TR/2
end
bSSFPend.NstepTR = int32(BSIM.TR/BSIM.dt); % Not sure if ramp-down is TR or TR/2
bSSFPend.NstepTot = int32(bSSFPend.NendTR*bSSFPend.NstepTR);



time = [1:BSIM.NstepTot]*BSIM.dt;
% ========================================================================


% an ensemble of spins w/range of df
% f_vec = [-400:10:400]; %Hz
Nf = numel(f_vec);

% keep track of temporal evolution
Msimf = zeros(3, BSIM.NstepTot + bSSFPcat.NstepTot + bSSFPend.NstepTot, Nf);


% Simulate!
for f=1:Nf
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
    % update magnetization history
    BSIM.M = Msim(:,end);
    BSIM.RFphi0 = RF_phi(end);
        
    % bSSFP half-angle restore pulse
    if (bSSFPend.RAMP_DOWN == 1)
        [MsimEnd, RF_phiEnd] = bSSFPrestorefunc(BSIM, 0, bSSFPend);
    else
        MsimEnd = []; RF_phiEnd = [];
    end
    
    % plot
    if( PLOT_EACHSPIN )
        figure(101); plot( time,Msim(1,:), time,Msim(2,:), time,Msim(3,:) );
        legend('Mx', 'My', 'Mz'); ylim([-1 1]); grid on;
        %title(sprintf('GRE with gradient + RF spoiling, spin %d/%d',p,Nphi));
        xlabel('ms'); ylabel('Normalized Magnetization');
        pause;
    end
    
    % keep track of each spin
    Msimf(:, :, f) = cat(2, MsimCat, Msim, MsimEnd);
    
    % clear history for each df
    BSIM.RFphi0 = 0;
    BSIM.M = M0;
    %BSIM.M = [0 0 0].'; % assume a perfect sat pulse to start
end

% extract Mxy
Mxy = squeeze( Msimf(1, :, :) + 1i*Msimf(2, :, :) );
% plot freq profile at TE, over all TRs

if( bSSFPcat.catmode==0 )
    plotTEidx = [0:(BSIM.NTR+bSSFPend.NendTR-1)]*BSIM.NstepTR + BSIM.TE/BSIM.dt;
elseif( bSSFPcat.catmode==1 )
    plotTEidx = bSSFPcat.NstepTot+[0:(BSIM.NTR+bSSFPend.NendTR-1)]*BSIM.NstepTR + BSIM.TE/BSIM.dt;
    plotTEidx = [1, plotTEidx];
elseif( bSSFPcat.catmode==2 )
    % If TR is TR/2 for RAMP_DOWN, this line needs to be changed
    plotTEidx = [0:(bSSFPcat.NcatTR+BSIM.NTR+bSSFPend.NendTR-1)]*double(BSIM.NstepTR) + BSIM.TE/BSIM.dt;
elseif( bSSFPcat.catmode==3 )
    plotTEidx = [0:(bSSFPcat.NcatTR+BSIM.NTR+bSSFPend.NendTR-1)]*BSIM.NstepTR + BSIM.TE/BSIM.dt;
end


MxyTE = Mxy(int32(plotTEidx), :);

%figure; imagesc(f_vec,1:numel(plotTEidx), abs(MxyTE)); colormap gray;
%xlabel('Hz'); ylabel('TR #');

% demod RF phase from MxyTE
RF_phi_all = [RF_phiCat; RF_phi; RF_phiEnd];
RF_phi_all = exp(-1i*RF_phi_all) * ones(1,Nf);
MxyTE = MxyTE .* RF_phi_all;


if( PLOT_SS )
    % RAMPDOWN is taking effect
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
    
    figure();
    pos_idx = imag(MxyTElast) >= 0;
    neg_idx = imag(MxyTElast) < 0;
    %pos_idx = real(MxyTElast) >= 0;
    %neg_idx = real(MxyTElast) < 0;
    
    MxyTElast_phaseenc = -abs(MxyTElast .* pos_idx ) + abs(MxyTElast .* neg_idx);
    % MxyTElast_phaseenc = -abs(MxyTElast .* pos_idx) + abs(MxyTElast .* neg_idx);
    plot(f_vec, MxyTElast_phaseenc);
    title('Balanced SSFP frequency profile: phased-enc');
    xlabel('Hz'); ylabel('Normalized Mxy');
    
end

end