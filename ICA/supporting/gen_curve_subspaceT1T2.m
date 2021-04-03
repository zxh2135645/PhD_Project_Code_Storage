function [curveU,curveS]=gen_curve_subspaceT1T2(Nseg,TR,alpha_deg,alpha0_deg,TI)

if nargin < 5
  TI=2*TR;
end

alpha = alpha_deg*pi/180;

vec = @(x)x(:);
e1 = @(R1)exp(-TR*R1);
e2 = @(R2)exp(-[12 20 30 40 50]*1e-3*R2);
Mss = @(e1,alpha)(1-e1) / (1-cos(alpha)*e1);
n = 1:Nseg;
Sint = @(A,e1,e2,alpha,B)vec(A * Mss(e1,alpha) * (1 + ((e1*cos(alpha)).^(n-1)).'*(B*e2-1)) * sin(alpha)).';
S = @(A,R1,R2,alpha,B)Sint(A,e1(R1),e2(R2),alpha,B);
% S = @(A,R1,alpha,B)Sint(A,e(R1),alpha,B)+Sint(A,e(R1),alpha/2,B); %slice prof (unnecessary! already linear combination of multiple flip angles)

R1s = 1./logspace(log10(.1),log10(3),21);
R2s = 1./logspace(log10(.01),log10(3),21);
alphas = (.5:.5:alpha_deg)*pi/180;
if alpha0_deg == 90
  minB = -.1; %(1 - e(1/3)^(TI/TR))/(Mss(e(1/3),alpha))
  maxB = .25; %(1 - e(1/.1)^(TI/TR))/(Mss(e(1/.1),alpha)) %assume 2 TR to first imaging pulse, shortest T1
  Bs = linspace(minB,maxB,6);
elseif alpha0_deg == 180
%   minB = (1 - e(1/3))/Mss(e(1/3),alpha) - e(1/3)/(Mss(e(1/3),alpha) * sin(alpha))
  minB = -1
%   B = S(1,1/3,alpha,minB);
%   maxB = (1 - e(1/3)^(TI/TR))/Mss(e(1/3),alpha) - B(end) * e(1/3)/(Mss(e(1/3),alpha) * sin(alpha))
maxB = -.5
  Bs = linspace(minB,maxB,6);
%   Bs = unique([-2:.25:-1, linspace(-1,1,21)  1:.25:maxB]); %IR
else
  disp('Error: unrecognized preparation pulse')
  return
end
curves=zeros(Nseg*5,numel(R1s),numel(alphas),numel(Bs),numel(R2s));

for j=1:numel(R1s)
  for k=1:numel(alphas)
    for l=1:numel(Bs)
      for m=1:numel(R2s)
        curves(:,j,k,l,m) = S(1,R1s(j),R2s(m),alphas(k),Bs(l));
      end
    end
  end
end

[~,curveS,curveU] = svd(curves(:,:)','econ');
curveS=diag(curveS);

return