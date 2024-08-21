% Bloch Equation Simulation
% B -- Basic Sequence Simulations, Part 1
% B-1. Saturation-Recovery
%
% Modified: 2014/03/05
% created: 06/28/2005 
% Holden H. Wu

close all; clear all;

% equilibrium 
M0z = 1;
BSIM.M = [0 0 M0z]';
% sim settings
BSIM.T1 = 600;  %ms
BSIM.T2 = 100;  %ms
BSIM.df = 0;    %Hz
BSIM.TR = 500;  %ms
BSIM.TE = 1;    %ms
BSIM.flipy = 60 *(pi/180); %pi/3; %rad

BSIM.NTR = 10;  %number of TRs to sim
BSIM.dt  = 0.5;   %ms, sim step size in time
BSIM.NstepTR = BSIM.TR/BSIM.dt;  % sim steps per TR
BSIM.NstepTot= BSIM.NTR*BSIM.NstepTR; % total sim steps

Msim = zeros(3, BSIM.NstepTot); %keep track of temporal evolution
step = 1;

% init free precession matrix
% note depedence on df
[Adt, Bdt] = freeprecess( BSIM.dt, BSIM.T1, BSIM.T2, BSIM.df );

% Sim TR by TR
% assuming instantaneous RF excitation and grad spoiling
for nn=1:BSIM.NTR,
    % excitation
    if(nn==1)
        Mnew = yrot(BSIM.flipy)*BSIM.M;
    else
        Mnew = yrot(BSIM.flipy)*Msim(:,step-1);
    end
    % first free precession step in this TR
    % included in the same time interval as excitation
    Mnew = Adt*Mnew + Bdt;
    Msim(:,step) = Mnew; step = step+1;
    
    % subsequent steps include free precession
    for(ss=2:BSIM.NstepTR)
        Mnew = Adt*Msim(:,step-1) + Bdt;
        % final step includes gradient spoiling
        % if(ss==NstepTR)
        % ...
        % end
        
        Msim(:,step) = Mnew; step = step+1;
    end
end

Mxy = Msim(1,:) + 1i*Msim(2,:);

time = [1:BSIM.NstepTot]*BSIM.dt;
figure; plot( time,Msim(1,:), time,Msim(2,:), time,abs(Mxy), time,Msim(3,:) );
legend('Mx', 'My', 'Mxy', 'Mz'); ylim([-1 1]); grid on;
title('Saturation Recovery');
xlabel('ms'); ylabel('Normalized Magnetization');
