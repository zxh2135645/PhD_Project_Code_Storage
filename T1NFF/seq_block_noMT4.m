function [t, M_total, Mxy_readout, t_readout] = seq_block_noMT4(TD, npulse, T1, T2, alpha, TR, prep, M0, df, RAMP_DOWN, adiabatic, ddt)

% This version Incorperated adiabatic pulse
% adiabatic is struct storing adiabatic pulse
% Add ddt as temporal resolution 0.1 ms
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments:
%       TD:         trigger delay (ms)        e.g. 561.8
%       npulse:     number of pulses          e.g. 60
%       T1:         T1 (ms)                   e.g. 1000
%       T2:         T2 (ms)                   e.g. 45
%       alpha:      flip angle (degree)       e.g. 35
%       TR:         Repetition time (ms)      e.g. 2.4
%       prep:       1x1 struct
%                   flip, t_delay
%       M0:         Initial magnetization     e.g. [0;0;1]
%       df:         Off-resonance (Hz)        e.g. 0 (On-resonance)
%       RAMP_DOWN:  If ramp-down pulse is applied  e.g. 0 or 1
%       adiabatic:  1x1 struct
%                   mu, beta1, pulseWidth, A0
%       ddt:        temporal resolution (ms)       e.g. 0.1 
%
%
% Outputs:
%       t:           Time array (ms)           e.g. 1x8174 double (0 - 817.3)
%       M_total:     temporal evolution of M   e.g. 3x8174 double [5e-05;3e-04;-0.3745]
%       Mxy_readout: Mxy at the readout        e.g. 0.0003 + 0.1538i  complex double 
%       t_readout:   at the readout time point e.g. 745.3 (ms) (t still doesn't match)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = ddt;
t_delay = prep.t_delay;
total_time = TD + TR * npulse + t_delay; % should be half of the duration strictly - Nah
one_rep_t = 0:dt:(total_time-dt); % Start with 1 or 0 
one_rep_M = repmat(M0, [1, length(one_rep_t)]);
num_rampup = 5; % linear ramp-up of 5 is built-in

[Adt,Bdt] = freeprecess(dt,T1,T2,df);

% Need to relax before prep pulse
recovery_timepoint = round((TD - adiabatic.pulseWidth/2)/dt);
M_recov_array = zeros(3, recovery_timepoint);
for i = 1:recovery_timepoint
    if i == 1
        M_recov_array(:,i) = M0;
    else
        M_recov_array(:,i) = Adt*M_recov_array(:,i-1) + Bdt;
    end
end
one_rep_M(:, 1:recovery_timepoint) = M_recov_array;


% Apply inversion pulse and bSSFP readout
initM = M_recov_array(:,end);

% Add adiabatic pulse
ad_tp = round(adiabatic.pulseWidth / dt);
if ad_tp ~= 0
    [mx0, my0, mz0] = AdiabaticPulse(adiabatic, T1, T2, df, M_recov_array(:,end));
    pos = length(mx0);
    div = round(pos / ad_tp);
    ad_m = zeros(3, ad_tp);
    
    for i = 1:ad_tp
        ad_m(1, ad_tp-i+1) = mx0(pos);
        ad_m(2, ad_tp-i+1) = my0(pos);
        ad_m(3, ad_tp-i+1) = mz0(pos);
        pos = pos - div;
    end
    one_rep_M(:,(recovery_timepoint+1):(recovery_timepoint+ad_tp)) = ad_m;
else
    ad_m = initM;
end


invtM = [ad_m(1, end), ad_m(2, end), ad_m(3, end)]'; 
% between prep pulse and readout
ti_tp = round(t_delay/dt) - ad_tp / 2;
TI_timepoint = length(one_rep_M(1,recovery_timepoint+ad_tp+1:recovery_timepoint+ad_tp+ti_tp));
if TI_timepoint ~= 0
    M_TI_array = zeros(3, TI_timepoint);
    for i = 1:TI_timepoint
        if i == 1
            M_TI_array(:,i) = invtM;
        else
            M_TI_array(:,i) =  Adt*M_TI_array(:,i-1) + Bdt;
        end
    end
    % Add inversion pulse
    one_rep_M(:,recovery_timepoint+ad_tp+1:recovery_timepoint+ad_tp+ti_tp) = M_TI_array;
    initM0 = M_TI_array(:,end);
else
    initM0 = initM;
end

% Debug flags
PLOT_EACHSPIN = 0;
PLOT_SS = 0; % PLOT steady state

% Sim settings
% ========================================================================
% spin properties
M02 = initM0;
FA = alpha;
NTR = npulse-num_rampup-RAMP_DOWN;

% ddt = 0.1; %ms
f_vec = [-1000:10:1000];
df_idx = find(df == f_vec);

%%%%%%%%%
%% TODO Here (debug)
[MxyTE, BSIM, plotTEidx, bSSFPcat, Msimf, bSSFPend] = bSSFP_engine(M02, T1, T2, df, TR, FA, NTR, ddt, f_vec, PLOT_EACHSPIN, PLOT_SS, RAMP_DOWN);

% Update one_rep_M
one_rep_M(:,end-npulse*TR/dt+1:end) = Msimf(:,:,df_idx);

% extract M
MTE = Msimf(:, int32(plotTEidx), :);

t = one_rep_t;
M_total = one_rep_M;

t_center_kspace_idx = num_rampup + round((npulse - num_rampup - RAMP_DOWN) / 2);
Mxy_readout = MxyTE(t_center_kspace_idx, df_idx);
fprintf([num2str(abs(Mxy_readout)), '\n'])
MxyTE2 = MTE(1, t_center_kspace_idx, df_idx)+1i*MTE(2, t_center_kspace_idx, df_idx);
fprintf([num2str(abs(MxyTE2)), '\n'])

t_center_kspace_idx_ddt = round(t_center_kspace_idx * (TR / ddt));
t_readout = t(recovery_timepoint + ad_tp + TI_timepoint + t_center_kspace_idx_ddt);

end