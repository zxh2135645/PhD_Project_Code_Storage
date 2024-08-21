function [Msim, RF_phi] = bSSFPprepfunc(BSIM, spoil_phi, bSSFPcat)
% Bloch Equation Simulation
%   Gradient echo with options for gradient spoiling and/or RF phase incr
%
% Modified:
%   2014/03/24 add options for bSSFP catalyzation
%   2014/03/05
% Created: 06/28/2005 
% Holden H. Wu

if(nargin<2)
    % phase from gradient spoiling
    spoil_phi = 0;
end
if(nargin<3)
    bSSFPcat.catmode  = 1; %0:none, 1:half prep
    bSSFPcat.NcatTR   = 1;
    bSSFPcat.NstepTR  = BSIM.TR/BSIM.dt;
    bSSFPcat.NstepTot = bSSFPcat.NcatTR*bSSFPcat.NstepTR/2;
end

RF_phi0 = BSIM.RFphi0;
RF_phi  = zeros(bSSFPcat.NcatTR, 1);
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

Msim = zeros(3, bSSFPcat.NstepTot); %keep track of temporal evolution
step = 1;

% init free precession matrix
[Adt, Bdt] = freeprecess( BSIM.dt, BSIM.T1, BSIM.T2, BSIM.df );

% Sim prep cycles TR by TR
% assuming instantaneous RF excitation and grad spoiling

% special treatment for th/2-TR/2 prep (TR is halved)
if( bSSFPcat.catmode==1 )
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
    for(ss=2:bSSFPcat.NstepTR/2)
        Mnew = Adt*Msim(:,step-1) + Bdt;
        % final step includes gradient spoiling
        if(ss==bSSFPcat.NstepTR/2)
            Mnew = zrot(spoil_phi)*Mnew;
        end
        
        Msim(:,step) = Mnew; step = step+1;
    end    

end


if( bSSFPcat.catmode>1 )
    
    % linear ramp
    if( bSSFPcat.catmode==2 )
        RF_thrad = BSIM.flipy*[1:bSSFPcat.NcatTR]/bSSFPcat.NcatTR;
        %figure; plot(RF_thrad);
    end
    
    % cosine-weighted ramp
    if( bSSFPcat.catmode==3 )
        RF_thrad = BSIM.flipy* 0.5*(1 - cos(pi*[1:bSSFPcat.NcatTR]/bSSFPcat.NcatTR));
        %figure; plot(RF_thrad);
    end
    
    for nn=1:bSSFPcat.NcatTR,
        % RF phase increment
        if( BSIM.RFCYC==2 )
            % quadratic for RF spoiling
            RF_incrphi = RF_incrphi + BSIM.RFdph;
        end
        % excitation and free precession
        if(nn==1)
            RF_phi(nn) = RF_phi0 + RF_incrphi;
            Mnew = throt(RF_thrad(nn), RF_phi(nn))*BSIM.M;
        else
            RF_phi(nn) = RF_phi(nn-1) + RF_incrphi;
            Mnew = throt(RF_thrad(nn), RF_phi(nn))*Msim(:,step-1);
        end    
        Mnew = Adt*Mnew + Bdt;
        Msim(:,step) = Mnew; step = step+1;

        % subsequent steps include free precession
        for(ss=2:bSSFPcat.NstepTR)
            Mnew = Adt*Msim(:,step-1) + Bdt;
            % final step includes gradient spoiling
            if(ss==bSSFPcat.NstepTR)
                Mnew = zrot(spoil_phi)*Mnew;
            end

            Msim(:,step) = Mnew; step = step+1;
        end    

    end % TR loop
end


return;
