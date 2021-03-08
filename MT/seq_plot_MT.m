function [t, Mzmt_total, mxymt_array, mxymt_array2, mzmt_array, mzmt_array2, Mzmt_bound_total] ...
    = seq_plot_MT(TD, npulse, alpha, TR, MT_para, MT_prep)

% linear ramp-up
num_rampup = 10;
phi = RF_phase_cycle(npulse,'balanced');
fa_flow = d2r(alpha)*ones(npulse,1);
rampup = linspace(0, alpha, num_rampup+1);
rampup = rampup(end-num_rampup+1:end);
fa_flow(1:num_rampup) = d2r(rampup);
t_delay = MT_prep.t_delay;

gam = 267.5221 *1e-3; % rad /ms /uT
trf = MT_para.trf; % ms
b1 = d2r(alpha)./(trf.*gam); % uT
b1sqrdtau = 2^2 * b1.^2.*trf;
b1sqrdtau0 = 2^2 * (d2r(rampup)./(trf.*gam)).^2.*trf;
b1sqrdtau_array = [b1sqrdtau0'; b1sqrdtau * ones(npulse-1,1)];

[ff,G] = SuperLorentzian(MT_para.T2x(2)*1e-3);
idx = find(ff == min(abs(ff)));
MT_para.G = G(idx);
MT_para.b1sqrdtau_array = b1sqrdtau_array;

% Add recovery with MT effect here

% Do bSSFP readout
[~,fnmt,Znmt] = EPGX_GRE_MT(fa_flow,phi,MT_para.b1sqrdtau_array,...
    TR,MT_para.T1x,MT_para.T2x(1),MT_para.F,MT_para.Kf,MT_para.G,...
    'kmax',inf, 'prep', MT_prep);
mxymt = size(fnmt,1)*ifftshift(ifft(ifftshift(fnmt,1),[],1),1);
mzmt = size(Znmt,1)*ifftshift(ifft(ifftshift(Znmt,1),[],1),1);

mxymt_array = abs(mxymt(floor(size(mxymt, 1)/2),:));
mzmt_array = real(mzmt(floor(size(mzmt, 1)/2+1),:)); % size is 1x(2*npulse)
% Thus, first half is free pool Mz
% Second half is bound pool Mz
% Need to use mzmt_array(end) as input for next EPG simulation


dt = TR;
total_time = TD + TR * npulse + t_delay;
one_rep_t = 0:dt:total_time;
one_rep_Mzmt = ones(1, length(one_rep_t));
TI_timepoint = length(one_rep_Mzmt(end-npulse-round(t_delay/dt)+1:end-npulse));

ti_flow = zeros(TI_timepoint, 1);
phi_ti = zeros(TI_timepoint, 1);
b1sqrdtau_array = zeros(TI_timepoint,1);
MT_prep2 = MT_prep;
MT_prep2.t_delay = 0;
[~,~,Znmt_ti] = EPGX_GRE_MT(ti_flow,phi_ti,b1sqrdtau_array,...
    TR,MT_para.T1x,MT_para.T2x(1),MT_para.F,MT_para.Kf,MT_para.G,...
    'kmax',inf, 'prep', MT_prep2);
mzmt_ti = size(Znmt_ti,1)*ifftshift(ifft(ifftshift(Znmt_ti,1),[],1),1);
mzmt_array_ti = real(mzmt_ti(floor(size(mzmt_ti, 1)/2),:));

% Sum up first shot
one_rep_Mzmt(end-npulse+1:end) = mzmt_array(1:npulse);
one_rep_Mzmt(end-npulse-TI_timepoint+1:end-npulse) = mzmt_array_ti(1:TI_timepoint);
% Sum up bound pool
one_rep_Mzmt_bound = MT_para.F * ones(1, length(one_rep_t));
one_rep_Mzmt_bound(end-npulse+1:end) = mzmt_array(1+npulse:end);
one_rep_Mzmt_bound(end-npulse-TI_timepoint+1:end-npulse) = mzmt_array_ti(TI_timepoint+1:end);


% For rep two
two_rep_t = 0:dt:total_time;
two_rep_Mzmt = ones(1, length(two_rep_t));
TD_timepoint = length(two_rep_Mzmt(1:end-npulse-round(t_delay/dt)));
% For better visualization
td_flow = zeros(TD_timepoint, 1);
phi_td = zeros(TD_timepoint, 1);
b1sqrdtau_array = zeros(TD_timepoint,1);
initMmt2 = [0 0 one_rep_Mzmt(end) one_rep_Mzmt_bound(end)]';
[~,~,Znmt_td] = EPGX_GRE_MT(td_flow,phi_td,b1sqrdtau_array,TR,...
    MT_para.T1x,MT_para.T2x(1),MT_para.F,MT_para.Kf,MT_para.G,...
    'kmax',inf, 'initMmt', initMmt2);
mzmt_td = size(Znmt_td,1)*ifftshift(ifft(ifftshift(Znmt_td,1),[],1),1);
mzmt_array_td = real(mzmt_td(floor(size(mzmt_td, 1)/2+1),:));


% The second TI
ti_flow2 = zeros(TI_timepoint, 1);
phi_ti2 = zeros(TI_timepoint, 1);
b1sqrdtau_array2 = zeros(TI_timepoint,1);
MT_prep3 = MT_prep;
MT_prep3.t_delay = 0;
initMmt3 = [0 0 mzmt_array_td(TD_timepoint) mzmt_array_td(end)]';
[~,~,Znmt_ti2] = EPGX_GRE_MT(ti_flow2,phi_ti2,b1sqrdtau_array2,TR,...
    MT_para.T1x,MT_para.T2x(1),MT_para.F,MT_para.Kf,MT_para.G,...
    'kmax',inf, 'prep', MT_prep3, 'initMmt', initMmt3);
mzmt_ti2 = size(Znmt_ti2,1)*ifftshift(ifft(ifftshift(Znmt_ti2,1),[],1),1);
mzmt_array_ti2 = real(mzmt_ti2(floor(size(mzmt_ti2, 1)/2),:));

% second shot readout
initMmt = [0 0 mzmt_array_ti2(TI_timepoint) mzmt_array_ti2(end)]';
[~,fnmt2,Znmt2] = EPGX_GRE_MT(fa_flow,phi,MT_para.b1sqrdtau_array,...
    TR,MT_para.T1x,MT_para.T2x(1),MT_para.F,MT_para.Kf,MT_para.G,...
    'kmax',inf, 'initMmt', initMmt);

mxymt2 = size(fnmt2,1)*ifftshift(ifft(ifftshift(fnmt2,1),[],1),1);
mzmt2 = size(Znmt2,1)*ifftshift(ifft(ifftshift(Znmt2,1),[],1),1);
mxymt_array2 = abs(mxymt2(floor(size(mxymt2, 1)/2),:));
mzmt_array2 = real(mzmt2(floor(size(mzmt2, 1)/2+1),:));

% Sum up 2nd shot
two_rep_Mzmt = [mzmt_array_td(1:TD_timepoint), mzmt_array_ti2(1:TI_timepoint), mzmt_array2(1:npulse)];
two_rep_Mzmt_bound = [mzmt_array_td(1+TD_timepoint:end), mzmt_array_ti2(1+TI_timepoint:end), mzmt_array2(1+npulse:end)];
% combine two shots
t = [one_rep_t, one_rep_t(end)+two_rep_t];
Mzmt_total = [one_rep_Mzmt, two_rep_Mzmt];
Mzmt_bound_total = [one_rep_Mzmt_bound, two_rep_Mzmt_bound];
end