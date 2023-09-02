function bSSFP_engine_EPGX(M0, T1, T2, TR, FA, NTR, dt, f_vec, num_rampup, varargin)
% Try to create a bSSFP simulation with EPG-X engine
% But it seems unnecessary
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
%BSIM.df = df;    %Hz

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
bSSFPcat.NcatTR   = num_rampup; % for T1 MOLLI, 5 linear ramp-up
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

fa_flow = zeros(BSIM.NstepTot + bSSFPcat.NstepTot + bSSFPend.NstepTot, 1);
% for RF phase cycling
phi = RF_phase_cycle(BSIM.NstepTot + bSSFPcat.NstepTot + bSSFPend.NstepTot,'balanced');

% MT Energy deposition
gam = 267.5221 *1e-3; % rad /ms /uT
trf = MT_para.trf; % ms
b1 = d2r(alpha)./(trf.*gam); % uT
b1sqrdtau = b1.^2.*trf; 
b1sqrdtau0 = (d2r(rampup)./(trf.*gam)).^2.*trf;
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

[~,~,Znmt_ti] = EPGX_GRE_MT(fa_flow,phi,b1sqrdtau_array,...
        TR,MT_para.T1x,MT_para.T2x(1),MT_para.F,MT_para.Kf,MT_para.G,...
        'kmax',inf, 'prep', MT_prep2, 'initMmt', initM2);
    
end