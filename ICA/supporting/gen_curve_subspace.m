function [curveU,curveS]=gen_curve_subspace(Nseg,TR,alpha_deg,alpha0_deg,TI)

if nargin < 5
  TI=2*TR;
end

alpha = alpha_deg*pi/180;

e = @(R1)exp(-TR*R1);
Mss = @(e,alpha)(1-e) / (1-cos(alpha)*e);
n = 1:Nseg;
Sint = @(A,e,alpha,B)A * Mss(e,alpha) * (1 + (B-1)*(e*cos(alpha)).^(n-1)) * sin(alpha);
S = @(A,R1,alpha,B)Sint(A,e(R1),alpha,B);
% S = @(A,R1,alpha,B)Sint(A,e(R1),alpha,B)+Sint(A,e(R1),alpha/2,B); %slice prof (unnecessary! already linear combination of multiple flip angles)

R1s = 1./logspace(log10(.1),log10(3),101);
alphas = (.5:.5:(alpha_deg*1.5))*pi/180;
if strcmp(alpha0_deg,'T2IR')
  minB=-1
  maxB=0
  Bs = linspace(minB,maxB,21);
elseif strcmp(alpha0_deg,'T2')
  minB=0
  maxB=1
  Bs = linspace(minB,maxB,21);
elseif alpha0_deg == 90
  minB = -.1; %(1 - e(1/3)^(TI/TR))/(Mss(e(1/3),alpha))
  maxB = .25; %(1 - e(1/.1)^(TI/TR))/(Mss(e(1/.1),alpha)) %assume 2 TR to first imaging pulse, shortest T1
  Bs = linspace(minB,maxB,21);
elseif alpha0_deg == 180
  %   minB = (1 - e(1/3))/Mss(e(1/3),alpha) - e(1/3)/(Mss(e(1/3),alpha) * sin(alpha))
  minB = -1
  %   B = S(1,1/3,alpha,minB);
  %   maxB = (1 - e(1/3)^(TI/TR))/Mss(e(1/3),alpha) - B(end) * e(1/3)/(Mss(e(1/3),alpha) * sin(alpha))
  maxB = -0.5
  Bs = linspace(minB,maxB,21);
  %   Bs = unique([-2:.25:-1, linspace(-1,1,21)  1:.25:maxB]); %IR
else
  disp('Error: unrecognized preparation pulse')
  return
end
curves=zeros(Nseg,numel(R1s),numel(alphas),numel(Bs));

for j=1:numel(R1s)
  for k=1:numel(alphas)
    for l=1:numel(Bs)
      curves(:,j,k,l) = S(1,R1s(j),alphas(k),Bs(l));
    end
  end
end

[~,curveS,curveU] = svd(curves(:,:)','econ');
curveS=diag(curveS);

return