function [t, Mz_total, mxys_array, mxys_array2, mzs_array, mzs_array2] = seq_plot_noMT(TD, npulse, T1, T2, alpha, TR, prep)
t_delay = prep.t_delay;
phi = RF_phase_cycle(npulse,'balanced');

num_rampup = 10;
fa_flow = d2r(alpha)*ones(npulse,1); % Change to linear ramp-up
rampup = linspace(0, alpha, num_rampup+1);
rampup = rampup(end-num_rampup+1:end);
fa_flow(1:num_rampup) = d2r(rampup);
[~,fn0, Zn0] = EPG_GRE(fa_flow,phi,TR,T1,T2,'kmax',inf, 'prep', prep);

mxys = size(fn0,1)*ifftshift(ifft(ifftshift(fn0,1),[],1),1);
mzs = size(Zn0,1)*ifftshift(ifft(ifftshift(Zn0,1),[],1),1);
mxys_array = abs(mxys(floor(size(mxys, 1)/2),:));
mzs_array = real(mzs(floor(size(mzs, 1)/2+1),:));

Mz = mzs_array(end);
Mz2 = 1 - (1 - Mz).*exp(-TD/T1);
initM = [0 0 Mz2]';
[~,fn02, Zn02] = EPG_GRE(fa_flow,phi,TR,T1,T2,'kmax',inf, 'prep', prep, 'initM', initM);
mxys2 = size(fn02,1)*ifftshift(ifft(ifftshift(fn02,1),[],1),1);
mzs2 = size(Zn02,1)*ifftshift(ifft(ifftshift(Zn02,1),[],1),1);
mxys_array2 = abs(mxys2(floor(size(mxys2, 1)/2),:));
mzs_array2 = real(mzs2(floor(size(mzs2, 1)/2+1),:));

% For better accuracy
dt = TR;
total_time = TD + TR * npulse + t_delay;
one_rep_t = 0:dt:total_time;
one_rep_Mz = ones(1, length(one_rep_t));

one_rep_Mz(end-npulse+1:end) = mzs_array;

TI_timepoint = length(one_rep_Mz(end-npulse-round(t_delay/dt)+1:end-npulse));
E1 = exp(-dt / T1);
Mz_TI_array = zeros(1, TI_timepoint);
for i = 1:TI_timepoint
    if i == 1
        Mz_TI_array(i) = -1;
    else
        Mz_TI_array(i) = 1 - E1 + Mz_TI_array(i-1) * E1;
    end
end

% Add inversion pulse
one_rep_Mz(end-npulse-round(t_delay/dt)+1:end-npulse) = Mz_TI_array;


% Second rep
two_rep_t = 0:dt:total_time;
% disp(length(two_rep_t))
two_rep_Mz = ones(1, length(two_rep_t));

% Trigger Delay
TD_timepoint = length(two_rep_Mz(1:end-npulse-round(t_delay/dt)));
Mz_TD_array2 = zeros(1, TD_timepoint);
for i = 1:TD_timepoint
    if i == 1
        Mz_TD_array2(i) = 1 - E1 + one_rep_Mz(end) * E1;
    else
        Mz_TD_array2(i) = 1 - E1 + Mz_TD_array2(i-1) * E1;
    end
end

% Inversion time
TI_timepoint = length(two_rep_Mz(end-npulse-round(t_delay/dt)+1:end-npulse));
Mz_TI_array2 = zeros(1, TI_timepoint);
for i = 1:TI_timepoint
    if i == 1
        Mz_TI_array2(i) = 1 - E1 + (-Mz_TD_array2(end)) * E1;
    else
        Mz_TI_array2(i) = 1 - E1 + Mz_TI_array2(i-1) * E1;
    end
end

% To sum up
two_rep_Mz = [Mz_TD_array2, Mz_TI_array2, mzs_array2];
t = [one_rep_t, one_rep_t(end)+two_rep_t];
% disp(length(t))
Mz_total = [one_rep_Mz, two_rep_Mz];

end