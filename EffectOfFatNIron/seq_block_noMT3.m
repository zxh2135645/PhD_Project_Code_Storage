function [t, M_total, Mxy_readout, t_readout] = seq_block_noMT3(TD, npulse, T1, T2, alpha, TR, prep, M0, df, RAMP_DOWN, adiabatic)

% This version Incorperated adiabatic pulse
% adiabatic is struct storing adiabatic pulse
% M0 is 3x1 array

dt = TR;
t_delay = prep.t_delay;
total_time = TD + TR * npulse + t_delay;
one_rep_t = 0:dt:total_time;
one_rep_M = repmat([0;0;1], [1, length(one_rep_t)]);
num_rampup = 5; % linear ramp-up of 5 is built-in

[Adt,Bdt] = freeprecess(dt,T1,T2,df);

% Need to relax before prep pulse
recovery_timepoint = round(TD/dt);
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

% invtM = cos(prep.flip) * initM
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
    one_rep_M(:,end-npulse-round(t_delay/dt)+1:(end-npulse-round(t_delay/dt)+ad_tp)) = ad_m;
else
    ad_m = initM;
end


invtM = [ad_m(1, end), ad_m(2, end), ad_m(3, end)]'; 
% between prep pulse and readout
TI_timepoint = length(one_rep_M(1,end-npulse-round(t_delay/dt)+ad_tp+1:end-npulse));
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
    one_rep_M(:,end-npulse-round(t_delay/dt)+ad_tp+1:end-npulse) = M_TI_array;
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

ddt = 0.1; %ms
f_vec = [-1000:10:1000];
df_idx = find(df == f_vec);

%%%%%%%%%
%% TODO Here (debug)
[MxyTE, BSIM, plotTEidx, bSSFPcat, Msimf, bSSFPend] = bSSFP_engine(M02, T1, T2, df, TR, FA, NTR, ddt, f_vec, PLOT_EACHSPIN, PLOT_SS, RAMP_DOWN);

% extract M
MTE = Msimf(:, int32(plotTEidx), :);

one_rep_M(:,end-npulse+1:end) = MTE(:,:,df_idx);

t = one_rep_t;
M_total = one_rep_M;

t_center_kspace_idx = num_rampup + round((npulse - num_rampup - RAMP_DOWN) / 2);
Mxy_readout = MxyTE(t_center_kspace_idx, df_idx);
fprintf([num2str(abs(Mxy_readout)), '\n'])
MxyTE2 = MTE(1, t_center_kspace_idx, df_idx)+1i*MTE(2, t_center_kspace_idx, df_idx);

fprintf([num2str(abs(MxyTE2)), '\n'])
t_readout = t(recovery_timepoint+TI_timepoint + t_center_kspace_idx);

end