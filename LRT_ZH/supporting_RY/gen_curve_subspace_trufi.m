function [curveU,curveS]=gen_curve_subspace_trufi(Nseg,TR,alpha_deg,alpha0_deg,TI)

if nargin < 5
  TI=2*TR;
end

alpha = alpha_deg*pi/180;

e1 = @(R1)exp(-TR*R1);
e2 = @(R2)exp(-TR*R2);
e1star=@(e1,e2,alpha)e1*cos(alpha/2)^2+e2*sin(alpha/2)^2;
Mss = @(e1,e2,alpha)(1-e1)/(1-(e1-e2)*cos(alpha)-e1*e2);
n = 1:Nseg;
Sint = @(A,e1,e2,alpha,B)A * Mss(e1,e2,alpha) * (1 + (B-1)*e1star(e1,e2,alpha).^n) * sin(alpha);
S = @(A,R1,R2,alpha,B)Sint(A,e1(R1),e2(R2),alpha,B);

R1s = 1./logspace(log10(.1),log10(3),51);
R2s = 1./logspace(log10(.05),log10(2),51);
alphas = (0:5:(alpha_deg*1.5))*pi/180;
if alpha0_deg == 90
  minB = -.1; %(1 - e(1/3)^(TI/TR))/(Mss(e(1/3),alpha))
  maxB = .25; %(1 - e(1/.1)^(TI/TR))/(Mss(e(1/.1),alpha)) %assume 2 TR to first imaging pulse, shortest T1
  Bs = linspace(minB,maxB,21);
elseif alpha0_deg == 180
%   minB = (1 - e(1/3))/Mss(e(1/3),alpha) - e(1/3)/(Mss(e(1/3),alpha) * sin(alpha))
  minB = -1
%   B = S(1,1/3,alpha,minB);
%   maxB = (1 - e(1/3)^(TI/TR))/Mss(e(1/3),alpha) - B(end) * e(1/3)/(Mss(e(1/3),alpha) * sin(alpha))
maxB = 0
  Bs = linspace(minB,maxB,21);
%   Bs = unique([-2:.25:-1, linspace(-1,1,21)  1:.25:maxB]); %IR
else
  disp('Error: unrecognized preparation pulse')
  return
end
curves=zeros(Nseg,numel(R1s),numel(R2s),numel(alphas),numel(Bs));

for j=1:numel(R1s)
    for k=1:numel(alphas)
        for l=1:numel(Bs)
            temp=zeros(numel(R2s),Nseg);
            R1=R1s(j);
            alpha=alphas(k);
            B=Bs(l);
            for m=1:numel(R2s) %parfor
                temp(m,:)= S(1,R1,R2s(m),alpha,B);
            end
            curves(:,j,:,k,l)=temp.';
        end
    end
end

[~,curveS,curveU] = svd(curves(:,:)','econ');
curveS=diag(curveS);

return