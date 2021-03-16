function [t, Mzmt_total, mxymt_array, mzmt_array, Mzmt_bound_total, Mxy_readout, corespond_t] ...
    = seq_block_MT(TD, npulse, alpha, TR, MT_para, MT_prep, M0, num_rampup, varargin)

% MT energy deposition (b1sqrdtau) was 4 times of the theory (maybe for 
% inline with agar simulation?)
% In this version it is now switched back
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
%
%
%   optional arguments (use string then value as next argument)
%       
%       restore_pulse: if apply alpha/2 restore pulse; Boolean
%
% Outputs:
%       t:           Time array (ms)        e.g. 1x340 double (0 - 813.6)
%       Mzmt_total:  Mz temporal evolution  e.g. 1x340 double
%       mxymt_array: Mxy of bSSFP           e.g. 1x66  double 
%       mzmt_array:  Mz of bSSFP (first half free pool/ second half bound pool)
%                    e.g. 1x132 double
%       Mzmt_bound_total:    Mz of complete of bound pool    e.g. 1x340 double
%       Mxy_readout: abs(Mxy) of the readout   e.g. 0.1221
%       correspond_t:at the readout time point e.g. 727.2 (ms)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for ii=1:length(varargin)
    if varargin{ii} == 1
        npulse = npulse + 1;
    end
end

% linear ramp-up
phi = RF_phase_cycle(npulse,'balanced');
fa_flow = d2r(alpha)*ones(npulse,1);
rampup = linspace(0, alpha, num_rampup+1);
rampup = rampup(end-num_rampup+1:end);
fa_flow(1:num_rampup) = d2r(rampup);
t_delay = MT_prep.t_delay;

% MT Energy deposition
gam = 267.5221 *1e-3; % rad /ms /uT
trf = MT_para.trf; % ms
b1 = d2r(alpha)./(trf.*gam); % uT
b1sqrdtau = b1.^2.*trf; 
%b1sqrdtau = 2*2 * b1.^2.*trf;
b1sqrdtau0 = (d2r(rampup)./(trf.*gam)).^2.*trf;
%b1sqrdtau0 = 2*2 * (d2r(rampup)./(trf.*gam)).^2.*trf;
for ii=1:length(varargin)
    if varargin{ii} == 1
        fa_flow(end) = d2r(alpha/2);
        b1sqrdtau_end = (d2r(alpha/2)./(trf.*gam)).^2.*trf;
        b1sqrdtau_array = [b1sqrdtau0'; b1sqrdtau * ones(npulse-1-num_rampup,1); b1sqrdtau_end];
        % fprintf('%d\n', length(b1sqrdtau_array))
    else
        b1sqrdtau_array = [b1sqrdtau0'; b1sqrdtau * ones(npulse-num_rampup,1)];
    end
end

[ff,G] = SuperLorentzian(MT_para.T2x(2)*1e-3);
idx = find(ff == min(abs(ff)));
MT_para.G = G(idx);

% For better accuracy
dt = TR;
total_time = TD + t_delay + TR * npulse ;
one_rep_t = 0:dt:total_time;
one_rep_Mzmt = ones(1, length(one_rep_t));
TI_timepoint = length(one_rep_Mzmt(end-npulse-round(t_delay/dt)+1:end-npulse));

% Add recovery here
recovery_timepoint = round(TD/dt);
fa_flow0 = zeros(1, recovery_timepoint);
phi0 = zeros(1, recovery_timepoint);
b1sqrdtau_array0 = zeros(1, recovery_timepoint);

[~,fnmt0,Znmt0] = EPGX_GRE_MT(fa_flow0,phi0,b1sqrdtau_array0,...
    TR,MT_para.T1x,MT_para.T2x(1),MT_para.F,MT_para.Kf,MT_para.G,...
    'kmax',inf, 'initMmt', M0);
mxymt0 = size(fnmt0,1)*ifftshift(ifft(ifftshift(fnmt0,1),[],1),1);
mzmt0 = size(Znmt0,1)*ifftshift(ifft(ifftshift(Znmt0,1),[],1),1);

if mod(recovery_timepoint, 2) == 0
    mzmt_array0 = real(mzmt0(floor(size(mzmt0, 1)/2+1),:));
else
    mzmt_array0 = real(mzmt0(floor(size(mzmt0, 1)/2),:));
end
    
one_rep_Mzmt_bound = MT_para.F * ones(1, length(one_rep_t));

% MT_prep 
initM1 = [0 0 mzmt_array0(recovery_timepoint) mzmt_array0(end)]';
[~,fnmt,Znmt] = EPGX_GRE_MT(fa_flow,phi, b1sqrdtau_array,...
    TR,MT_para.T1x,MT_para.T2x(1),MT_para.F,MT_para.Kf,MT_para.G,...
    'kmax',inf, 'prep', MT_prep, 'initMmt', initM1);
mxymt = size(fnmt,1)*ifftshift(ifft(ifftshift(fnmt,1),[],1),1);
mzmt = size(Znmt,1)*ifftshift(ifft(ifftshift(Znmt,1),[],1),1);

mxymt_array = abs(mxymt(floor(size(mxymt, 1)/2),:));
mzmt_array = real(mzmt(floor(size(mzmt, 1)/2+1),:)); % size is 1x(2*npulse)
% Thus, first half is free pool Mz
% Second half is bound pool Mz
% Need to use mzmt_array(end) as input for next EPG simulation

% 
if TI_timepoint ~= 0
    ti_flow = zeros(TI_timepoint, 1);
    phi_ti = zeros(TI_timepoint, 1);
    b1sqrdtau_array = zeros(TI_timepoint,1);
    MT_prep2 = MT_prep;
    MT_prep2.t_delay = 0;
    initM2 = [0 0 mzmt_array0(recovery_timepoint) mzmt_array0(end)]';
    [~,~,Znmt_ti] = EPGX_GRE_MT(ti_flow,phi_ti,b1sqrdtau_array,...
        TR,MT_para.T1x,MT_para.T2x(1),MT_para.F,MT_para.Kf,MT_para.G,...
        'kmax',inf, 'prep', MT_prep2, 'initMmt', initM2);
    mzmt_ti = size(Znmt_ti,1)*ifftshift(ifft(ifftshift(Znmt_ti,1),[],1),1);
    
    if mod(TI_timepoint, 2) == 0
        mzmt_array_ti = real(mzmt_ti(floor(size(mzmt_ti, 1)/2+1),:));
    else
        mzmt_array_ti = real(mzmt_ti(floor(size(mzmt_ti, 1)/2),:));
    end
    % To sum up
    one_rep_Mzmt(end-npulse-TI_timepoint+1:end-npulse) = mzmt_array_ti(1:TI_timepoint);
    % Sum up bound pool
   
    one_rep_Mzmt_bound(end-npulse-TI_timepoint+1:end-npulse) = mzmt_array_ti(TI_timepoint+1:end);
end

% Sum up first shot
one_rep_Mzmt(1:recovery_timepoint) = mzmt_array0(1:recovery_timepoint);
one_rep_Mzmt(end-npulse+1:end) = mzmt_array(1:npulse);
% one_rep_Mzmt(end-npulse-TI_timepoint+1:end-npulse) = mzmt_array_ti(1:TI_timepoint);

one_rep_Mzmt_bound(1:recovery_timepoint) = mzmt_array0(1+recovery_timepoint:end);
one_rep_Mzmt_bound(end-npulse+1:end) = mzmt_array(1+npulse:end);

% combine two shots
t = one_rep_t;
Mzmt_total = one_rep_Mzmt;
Mzmt_bound_total = one_rep_Mzmt_bound;

Mxy_readout = mxymt_array(round((npulse - num_rampup - 1) / 2));
corespond_t = t(recovery_timepoint+TI_timepoint+round((npulse - num_rampup - 1) / 2)); 
% Readout t

end