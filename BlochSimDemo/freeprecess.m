function [Afp,Bfp]=freeprecess(T,T1,T2,df) 
% T, T1, T2 in ms 
% df in Hz 

% Relaxation
M0 = 1;
A = [exp(-T/T2) 0 0; 0 exp(-T/T2) 0; 0 0 exp(-T/T1)];
B = M0*[0 0 1-exp(-T/T1)]';

% df in Hz  
phi = 2*pi * df*T*10^-3; %omega = 2pi * f, in radians
Rz = zrot( phi );

Afp = A*Rz;
% Bfp = B*Rz;
% same as: 
Bfp = B;
