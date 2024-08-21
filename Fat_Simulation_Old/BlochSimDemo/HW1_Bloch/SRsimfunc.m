function [Msim, RF_phi] = SRsimfunc(BSIM, spoil_phi)
% Bloch Equation Simulation
%   Gradient echo with options for gradient spoiling and/or RF phase incr
%
% Modified:
%   2014/03/05
% Created: 06/28/2005 
% Holden H. Wu

if(nargin<2)
    % phase from gradient spoiling
    spoil_phi = 0;
end

%RF_phi0 = 0;
RF_phi0 = BSIM.RFphi0;
RF_phi  = zeros(BSIM.NTR, 1);
if( BSIM.RFCYC==1 )    
    RF_incrphi = BSIM.RFdph;
elseif( BSIM.RFCYC==2 )
    RF_incrphi = BSIM.RFdph; %-BSIM.RFdph;
else
    RF_incrphi = 0;
end

% % equilibrium 
% M0z = 1;
% BSIM.M = [0 0 M0z]';
% % sim settings
% BSIM.T1 = 600;  %ms
% BSIM.T2 = 100;  %ms
% BSIM.df = 0;    %Hz
% BSIM.TR = 500;  %ms
% BSIM.TE = 1;    %ms
% BSIM.flipy = 60 *(pi/180); %pi/3; %rad
% 
% BSIM.NTR = 10;  %number of TRs to sim
% BSIM.dt  = 0.5;   %ms, sim step size in time
% BSIM.NstepTR = BSIM.TR/BSIM.dt;  % sim steps per TR
% BSIM.NstepTot= BSIM.NTR*BSIM.NstepTR; % total sim steps
% % for RF phase cycling
% BSIM.RFCYC = 1; % 0:none, 1:SSFP phase cycling, 2:SPGR spoiling
% BSIM.RFdph = pi; 

Msim = zeros(3, BSIM.NstepTot); %keep track of temporal evolution
step = 1;

% init free precession matrix
[Adt, Bdt] = freeprecess( BSIM.dt, BSIM.T1, BSIM.T2, BSIM.df );

% Sim TR by TR
% assuming instantaneous RF excitation and grad spoiling
for nn=1:BSIM.NTR,
    % RF phase increment
    if( BSIM.RFCYC==2 )
        % quadratic for RF spoiling
        RF_incrphi = RF_incrphi + BSIM.RFdph;
    end
    % excitation and free precession
    if(nn==1)
        %Mnew = yrot(BSIM.flipy)*BSIM.M;
        RF_phi(nn) = RF_phi0 + RF_incrphi;
        Mnew = throt(BSIM.flipy, RF_phi(nn))*BSIM.M;
    else
        %Mnew = yrot(BSIM.flipy)*Msim(:,step-1);
        RF_phi(nn) = RF_phi(nn-1) + RF_incrphi;
        Mnew = throt(BSIM.flipy, RF_phi(nn))*Msim(:,step-1);
    end    
    Mnew = Adt*Mnew + Bdt;
    Msim(:,step) = Mnew; step = step+1;
    
    % subsequent steps include free precession
    for(ss=2:BSIM.NstepTR)
        Mnew = Adt*Msim(:,step-1) + Bdt;
        % final step includes gradient spoiling
        if(ss==BSIM.NstepTR)
            Mnew = zrot(spoil_phi)*Mnew;
        end
        
        Msim(:,step) = Mnew; step = step+1;
    end
end

% Mxy = Msim(1,:) + 1i*Msim(2,:);
% 
% time = [1:BSIM.NstepTot]*BSIM.dt;
% figure; plot( time,Msim(1,:), time,Msim(2,:), time,abs(Mxy), time,Msim(3,:) );
% legend('Mx', 'My', 'Mxy', 'Mz'); ylim([-1 1]); grid on;
% title('Saturation Recovery');
% xlabel('ms'); ylabel('Normalized Magnetization');

return;