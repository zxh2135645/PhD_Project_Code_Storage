%
%	function [Msig,Mss] = gresignal(flip,T1,T2,TE,TR,dfreq)
% 
%	Calculate the steady state gradient-spoiled signal at TE for repeated
%	excitations given T1,T2,TR,TE in ms.  dfreq is the resonant
%	frequency in Hz.  flip is in radians.


function [Msig,Mss] = gresignal(flip,T1,T2,TE,TR,dfreq,label)


N = 100;
M = zeros(3,N);
phi = ([1:N]/N-0.5 ) * 4*pi;

if strcmp(label, 'FID')
    for k = 1:N    
        [M1sig,M1] = gssignal(flip,T1,T2,TE,TR,dfreq,phi(k));
        M(:,k)=M1;
    end
elseif strcmp(label, 'Echo')
    for k = 1:N
        [M1sig,M1] = ssfpechosignal(flip,T1,T2,TE,TR,dfreq,phi(k));
        M(:,k)=M1;
    end
else
    error('Please use FID or Echo as labels for spoiled-GRE simulation');
end
            

Mss = mean(M, 2);
Msig = Mss(1)+1i*Mss(2);