function [t, Mz_total, Mxy_readout, t_readout] = seq_block_noMT2(TD, npulse, T1, T2, alpha, TR, prep, Mz0, df, RAMP_DOWN)

dt = TR;
t_delay = prep.t_delay;
total_time = TD + TR * npulse + t_delay;
one_rep_t = 0:dt:total_time;
one_rep_Mz = ones(1, length(one_rep_t));
num_rampup = 5; % linear ramp-up of 5 is built-in

% Need to relax before prep pulse
recovery_timepoint = round(TD/dt);
E1 = exp(-dt / T1);
Mz_recov_array = zeros(1, recovery_timepoint);
for i = 1:recovery_timepoint
    if i == 1
        Mz_recov_array(i) = Mz0;
    else
        Mz_recov_array(i) = 1 - E1 + Mz_recov_array(i-1) * E1;
    end
end
one_rep_Mz(1:recovery_timepoint) = Mz_recov_array;


% Apply inversion pulse and bSSFP readout
initM = [0 0 Mz_recov_array(end)]';

invtM = cos(prep.flip) * initM;
% TODO here 
% Add adiabatic pulse

% between prep pulse and readout
TI_timepoint = length(one_rep_Mz(end-npulse-round(t_delay/dt)+1:end-npulse));
if TI_timepoint ~= 0
    Mz_TI_array = zeros(1, TI_timepoint);
    for i = 1:TI_timepoint
        if i == 1
            Mz_TI_array(i) = invtM(3);
        else
            Mz_TI_array(i) = 1 - E1 + Mz_TI_array(i-1) * E1;
        end
    end
    % Add inversion pulse
    one_rep_Mz(end-npulse-round(t_delay/dt)+1:end-npulse) = Mz_TI_array;
    initM0 = [0 0 Mz_TI_array(end)]';
else
    initM0 = initM;
end

% Debug flags
PLOT_EACHSPIN = 0;
PLOT_SS = 0; % PLOT steady state


% Sim settings
% ========================================================================
% spin properties
M0 = initM0;
FA = alpha;
NTR = npulse-num_rampup-RAMP_DOWN;

dt = 0.1; %ms
f_vec = [-1000:10:1000];
df_idx = find(df == f_vec);

[MxyTE, BSIM, plotTEidx, bSSFPcat, Msimf, bSSFPend] = bSSFP_engine(M0, T1, T2, df, TR, FA, NTR, dt, f_vec, PLOT_EACHSPIN, PLOT_SS, RAMP_DOWN);
% extract Mz
Mz = squeeze( Msimf(3, :, :));
MzTE = Mz(int32(plotTEidx), :);

one_rep_Mz(end-npulse+1:end) = MzTE(:,df_idx)';

t = one_rep_t;
Mz_total = one_rep_Mz;

t_center_kspace_idx = num_rampup + round((npulse - num_rampup - RAMP_DOWN) / 2);
Mxy_readout = MxyTE(t_center_kspace_idx, df_idx);
t_readout = t(recovery_timepoint+TI_timepoint + t_center_kspace_idx);

end