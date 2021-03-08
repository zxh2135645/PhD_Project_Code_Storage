function [t, Mz_total, mxys_array, mzs_array, Mxy_readout, corespond_t] = seq_block_noMT(TD, npulse, T1, T2, alpha, TR, prep, num_rampup, Mz0, varargin)
%% Extra variables
for ii=1:length(varargin)
    
    if varargin{ii} == 1
        npulse = npulse + 1;
    end
end

t_delay = prep.t_delay;
phi = RF_phase_cycle(npulse,'balanced');

% Initialization
fa_flow = d2r(alpha)*ones(npulse,1); % Change to linear ramp-up
rampup = linspace(0, alpha, num_rampup+1);
rampup = rampup(end-num_rampup+1:end);
fa_flow(1:num_rampup) = d2r(rampup);
for ii=1:length(varargin)
    if varargin{ii} == 1
        fa_flow(end) = d2r(alpha/2);
    end
end

dt = TR;
total_time = TD + TR * npulse + t_delay;
one_rep_t = 0:dt:total_time;
one_rep_Mz = ones(1, length(one_rep_t));

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
[~,fn0, Zn0] = EPG_GRE(fa_flow,phi,TR,T1,T2,'kmax',inf, 'prep', prep, 'initM', initM);

mxys = size(fn0,1)*ifftshift(ifft(ifftshift(fn0,1),[],1),1);
mzs = size(Zn0,1)*ifftshift(ifft(ifftshift(Zn0,1),[],1),1);
mxys_array = abs(mxys(floor(size(mxys, 1)/2),:));
if mod(npulse, 2) == 0
    mzs_array = real(mzs(floor(size(mzs, 1)/2+1),:));
else
    mzs_array = real(mzs(floor(size(mzs, 1)/2),:));
end
one_rep_Mz(end-npulse+1:end) = mzs_array;


% To retrospectively get mag between prep pulse and readout
% Don't worry, EXP simulation is doing the same thing.
TI_timepoint = length(one_rep_Mz(end-npulse-round(t_delay/dt)+1:end-npulse));
Mz_TI_array = zeros(1, TI_timepoint);
for i = 1:TI_timepoint
    if i == 1
        Mz_TI_array(i) = -Mz_recov_array(end);
    else
        Mz_TI_array(i) = 1 - E1 + Mz_TI_array(i-1) * E1;
    end
end

% Add inversion pulse
one_rep_Mz(end-npulse-round(t_delay/dt)+1:end-npulse) = Mz_TI_array;

t = one_rep_t;
% disp(length(t))
Mz_total = one_rep_Mz;

Mxy_readout = mxys_array(round((npulse - num_rampup - 1) / 2));
corespond_t = t(recovery_timepoint+TI_timepoint+round((npulse - num_rampup - 1) / 2));
end