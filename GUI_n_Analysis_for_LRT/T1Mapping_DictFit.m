clear all;
close all;

%
% Inputs: Nseg,TR,alpha_deg,alpha0_deg,TI
alpha = alpha_deg*pi/180;

e = @(R1)exp(-TR*R1);
Mss = @(e,alpha)(1-e) / (1-cos(alpha)*e);
n = 1:Nseg;
Sint = @(A,e,alpha,B)A * Mss(e,alpha) * (1 + (B-1)*(e*cos(alpha)).^(n-1)) * sin(alpha);
S = @(A,R1,alpha,B)Sint(A,e(R1),alpha,B);

R1s = 1./logspace(log10(.1),log10(3),101);
alphas = (.5:.5:(alpha_deg*1.5))*pi/180;

%   minB = (1 - e(1/3))/Mss(e(1/3),alpha) - e(1/3)/(Mss(e(1/3),alpha) * sin(alpha))
minB = -1
%   B = S(1,1/3,alpha,minB);
%   maxB = (1 - e(1/3)^(TI/TR))/Mss(e(1/3),alpha) - B(end) * e(1/3)/(Mss(e(1/3),alpha) * sin(alpha))
maxB = -0.5
Bs = linspace(minB,maxB,21);

curves=zeros(Nseg,numel(R1s),numel(alphas),numel(Bs));

for j=1:numel(R1s)
    for k=1:numel(alphas)
        for l=1:numel(Bs)
            curves(:,j,k,l) = S(1,R1s(j),alphas(k),Bs(l));
        end
    end
end