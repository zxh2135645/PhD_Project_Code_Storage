function [t, Mzmt_total, Mxymt_total, mxymt_array, mzmt_array, Mzmt_bound_total, Mxy_readout, corespond_t] ...
    = seq_block_MT2(TD, npulse, alpha, TR, MT_para, MT_prep, M0, num_rampup, RAMP_DOWN, ddt)

% MT energy deposition (b1sqrdtau) was 4 times of the theory (maybe for 
% inline with agar simulation?)
% In this version it is now switched back
% Add ddt as temporal resolution 0.1 ms
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% arguments:
%       TD:         trigger delay (ms)        e.g. 554.6
%       npulse:     number of pulses          e.g. 60
%       alpha:      flip angle (degree)       e.g. 35
%       TR:         Repetition time (ms)      e.g. 2.4
%       MT_para:    1x1 struct
%                   T1x, T2x, F, Kf, trf, G
%       MT_prep:    1x1 struct
%                   flip, t_delay, B1SqrdTau
%       M0:         Initial magnetization     e.g. [0;0;1]
%       num_rampup: number of linear ramp-up  e.g. 5
%       RAMP_DOWN:  if apply alpha/2 restore pulse; Boolean
%       ddt:        temporal resolution (ms)  e.g. 0.1
%
%
% Outputs:
%       t:            Time array (ms)         e.g. 1x8222 double (0 - 822.1)
%       Mzmt_total:   Mz temporal evolution   e.g. 1x8222  complex double
%       Mxymt_total:  Mxy temporal evolution  e.g. 1x8222  complex double
%       mxymt_array:  Mxy of bSSFP            e.g. 1x1584  double 
%       mzmt_array:   Mz of bSSFP (first half free pool/ second half bound pool)
%                     e.g. 1x3168 double
%       Mzmt_bound_total:    Mz of complete of bound pool    e.g. 1x8222 double
%       Mxy_readout:  Mxy of the readout   e.g. 0.0012 - 0.1107i   complex double
%       correspond_t: at the readout time point e.g. 748.9 (ms)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initiation
npulse = npulse + num_rampup + RAMP_DOWN;
dt = ddt;
t_delay = MT_prep.t_delay;
phi = RF_phase_cycle(npulse,'balanced');

% linear ramp-up
fa_flow = d2r(alpha)*ones(npulse,1);
rampup = linspace(0, alpha, num_rampup+1);
rampup = rampup(end-num_rampup+1:end);
fa_flow(1:num_rampup) = d2r(rampup);

% MT Energy deposition
gam = 267.5221 *1e-3; % rad /ms /uT
trf = MT_para.trf; % ms
b1 = d2r(alpha)./(trf.*gam); % uT
b1sqrdtau = b1.^2.*trf; 
b1sqrdtau0 = (d2r(rampup)./(trf.*gam)).^2.*trf;

if RAMP_DOWN == 1
    fa_flow(end) = d2r(alpha/2);
    b1sqrdtau_end = (d2r(alpha/2)./(trf.*gam)).^2.*trf;
    b1sqrdtau_array = [b1sqrdtau0'; b1sqrdtau * ones(npulse-1-num_rampup,1); b1sqrdtau_end];
    % fprintf('%d\n', length(b1sqrdtau_array))
else
    b1sqrdtau_array = [b1sqrdtau0'; b1sqrdtau * ones(npulse-num_rampup,1)];
end

% Super Lorentzian absorption lineshape
[ff,G] = SuperLorentzian(MT_para.T2x(2)*1e-3);
idx = find(ff == min(abs(ff)));
MT_para.G = G(idx);

% For better accuracy
dt = ddt;
total_time = TD + t_delay + TR * npulse ;
one_rep_t = 0:dt:(total_time-dt);
one_rep_Mzmt = ones(1, length(one_rep_t));
one_rep_Mxymt = ones(1, length(one_rep_t));

% 1) Add recovery here
recovery_timepoint = round(TD/dt);
fa_flow0 = zeros(1, recovery_timepoint);
phi0 = zeros(1, recovery_timepoint);
b1sqrdtau_array0 = zeros(1, recovery_timepoint);

[~,fnmt0,Znmt0] = EPGX_GRE_MT(fa_flow0,phi0,b1sqrdtau_array0,...
    dt,MT_para.T1x,MT_para.T2x(1),MT_para.F,MT_para.Kf,MT_para.G,...
    'kmax',inf, 'initMmt', M0);
mxymt0 = size(fnmt0,1)*ifftshift(ifft(ifftshift(fnmt0,1),[],1),1);
mzmt0 = size(Znmt0,1)*ifftshift(ifft(ifftshift(Znmt0,1),[],1),1);

mxymt_array0 = mxymt0(floor(size(mxymt0, 1)/2),:);
if mod(recovery_timepoint, 2) == 0
    mzmt_array0 = real(mzmt0(floor(size(mzmt0, 1)/2+1),:));
else
    mzmt_array0 = real(mzmt0(floor(size(mzmt0, 1)/2),:));
end
    

% MT_prep 
initM1 = [real(mxymt_array0(end)), imag(mxymt_array0(end)), mzmt_array0(recovery_timepoint), mzmt_array0(end)]';
% Sparsify fa_flow, phi, b1sqrdtau_array and TR
fa_flow_reshape = reshape([fa_flow'; zeros(round(TR/dt)-1, size(fa_flow',2))], [], 1);
phi_reshape = reshape([phi; zeros(round(TR/dt)-1, size(phi,2))], [], 1);
b1sqrdtau_array_reshape = reshape([b1sqrdtau_array'; zeros(round(TR/dt)-1, size(b1sqrdtau_array',2))], [], 1);
% IR-bSSFP simulation
[~,fnmt,Znmt] = EPGX_GRE_MT(fa_flow_reshape,phi_reshape, b1sqrdtau_array_reshape,...
    dt,MT_para.T1x,MT_para.T2x(1),MT_para.F,MT_para.Kf,MT_para.G,...
    'kmax',inf, 'prep', MT_prep, 'initMmt', initM1);
mxymt = size(fnmt,1)*ifftshift(ifft(ifftshift(fnmt,1),[],1),1);
mzmt = size(Znmt,1)*ifftshift(ifft(ifftshift(Znmt,1),[],1),1);

mxymt_array = mxymt(floor(size(mxymt, 1)/2),:);
mzmt_array = real(mzmt(floor(size(mzmt, 1)/2+1),:)); % size is 1x(2*npulse)
% Thus, first half is free pool Mz
% Second half is bound pool Mz
% Need to use mzmt_array(end) as input for next EPG simulation

TI_timepoint = length(one_rep_Mzmt(end-npulse*(TR/dt)-round(t_delay/dt)+1:end-npulse*(TR/dt)));
one_rep_Mzmt_bound = MT_para.F * ones(1, length(one_rep_t));

% 
if TI_timepoint ~= 0
    ti_flow = zeros(TI_timepoint, 1);
    phi_ti = zeros(TI_timepoint, 1);
    b1sqrdtau_array = zeros(TI_timepoint,1);
    MT_prep2 = MT_prep;
    MT_prep2.t_delay = 0;
    initM2 = [real(mxymt_array0(end)), imag(mxymt_array0(end)), mzmt_array0(recovery_timepoint), mzmt_array0(end)]';
    [~,fnmt_ti,Znmt_ti] = EPGX_GRE_MT(ti_flow,phi_ti,b1sqrdtau_array,...
        dt,MT_para.T1x,MT_para.T2x(1),MT_para.F,MT_para.Kf,MT_para.G,...
        'kmax',inf, 'prep', MT_prep2, 'initMmt', initM2);
    mxymt_ti = size(fnmt_ti,1)*ifftshift(ifft(ifftshift(fnmt_ti,1),[],1),1);
    mzmt_ti = size(Znmt_ti,1)*ifftshift(ifft(ifftshift(Znmt_ti,1),[],1),1);
    
    mxymt_array_ti = mxymt_ti(floor(size(mxymt_ti, 1)/2),:);
    if mod(TI_timepoint, 2) == 0
        mzmt_array_ti = real(mzmt_ti(floor(size(mzmt_ti, 1)/2+1),:));
    else
        mzmt_array_ti = real(mzmt_ti(floor(size(mzmt_ti, 1)/2),:));
    end
    % To sum up
    one_rep_Mzmt(end-npulse*round(TR/dt)-TI_timepoint+1:end-npulse*round(TR/dt)) = mzmt_array_ti(1:TI_timepoint);
    one_rep_Mxymt(end-npulse*round(TR/dt)-TI_timepoint+1:end-npulse*round(TR/dt)) = mxymt_array_ti(1:TI_timepoint);
    % Sum up bound pool
    
    one_rep_Mzmt_bound(end-npulse*round(TR/dt)-TI_timepoint+1:end-npulse*round(TR/dt)) = mzmt_array_ti(TI_timepoint+1:end);
    % one_rep_Mxymt(end-npulse*round(TR/dt)-TI_timepoint+1:end-npulse*round(TR/dt)) = mxymt_array_ti(TI_timepoint+1:end);

end

% Sum up first shot
one_rep_Mzmt(1:recovery_timepoint) = mzmt_array0(1:recovery_timepoint);
one_rep_Mzmt(end-npulse*round(TR/dt)+1:end) = mzmt_array(1:npulse*round(TR/dt));
one_rep_Mxymt(1:recovery_timepoint) = mxymt_array0(1:recovery_timepoint);
one_rep_Mxymt(end-npulse*round(TR/dt)+1:end) = mxymt_array(1:npulse*round(TR/dt));

one_rep_Mzmt_bound(1:recovery_timepoint) = mzmt_array0(1+recovery_timepoint:end);
one_rep_Mzmt_bound(end-npulse*round(TR/dt)+1:end) = mzmt_array(1+npulse*round(TR/dt):end);

% combine two shots
t = one_rep_t;
Mzmt_total = one_rep_Mzmt;
Mzmt_bound_total = one_rep_Mzmt_bound;

TEidx = round(((num_rampup + round((npulse - num_rampup - RAMP_DOWN) / 2)) * TR + TR/2) / dt);
Mxy_readout = mxymt_array(TEidx);
corespond_t = t(recovery_timepoint+TI_timepoint+TEidx); % Readout t

Mxymt_total = one_rep_Mxymt;

end