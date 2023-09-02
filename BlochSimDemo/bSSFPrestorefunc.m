function [Msim, RF_phi] = bSSFPrestorefunc(BSIM, spoil_phi, bSSFPend)

if(nargin<2)
    % phase from gradient spoiling
    spoil_phi = 0;
end
if(nargin<3)
    bSSFPend.RAMP_DOWN  = 0; %0:none, 1:half prep
    bSSFPend.NendTR   = 0;
    bSSFPend.NstepTR  = int32(BSIM.TR/BSIM.dt);
    bSSFPend.NstepTot = int32(bSSFPend.NendTR*bSSFPend.NstepTR/2);
end

RF_phi0 = BSIM.RFphi0;
RF_phi  = zeros(bSSFPend.NendTR, 1);
if( BSIM.RFCYC==1 )    
    RF_incrphi = BSIM.RFdph;
elseif( BSIM.RFCYC==2 )
    RF_incrphi = BSIM.RFdph; %-BSIM.RFdph;
else
    RF_incrphi = 0;
end

Msim = zeros(3, bSSFPend.NstepTot); %keep track of temporal evolution
step = 1;

% init free precession matrix
[Adt, Bdt] = freeprecess( BSIM.dt, BSIM.T1, BSIM.T2, BSIM.df );

% Sim restore pulse step by step
% assuming instantaneous RF excitation and grad spoiling 
% (Assumption puts here)

%  Not sure if TR needs to be halved
nn = 1;
% RF phase increment
if( BSIM.RFCYC==2 )
    % quadratic for RF spoiling
    RF_incrphi = RF_incrphi + BSIM.RFdph;
end
% excitation and free precession
RF_phi(nn) = RF_phi0 + RF_incrphi;
Mnew = throt(BSIM.flipy/2, RF_phi(nn))*BSIM.M;

Mnew = Adt*Mnew + Bdt;
Msim(:,step) = Mnew; step = step+1;

% subsequent steps include free precession
for(ss=2:bSSFPend.NstepTR)
    Mnew = Adt*Msim(:,step-1) + Bdt;
    % final step includes gradient spoiling
    % There is not gradient spoiling for this case
    if(ss==bSSFPend.NstepTR)
        Mnew = zrot(spoil_phi)*Mnew;
    end
    
    Msim(:,step) = Mnew; step = step+1;
end



end