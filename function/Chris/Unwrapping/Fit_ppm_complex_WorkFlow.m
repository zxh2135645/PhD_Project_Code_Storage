% %Modidied Chris Huang 2019 July
% Projection onto Dipole Fields (PDF)
%   [p1, dp1, relres, p0]=Fit_ppm_complex(M)
%    
%   output
%   p1 - field map, may need further unwrapping
%   dp1 - a priori error estimate
%   relres - relative residual
%   p0 - initial phase
%
%   input
%   M - a multi-echo and could be a multi-channel dataset
%       echo needs to be the 4th dimension
%       channel needs to be the 5th dimension
%
%   When using the code, please cite 
%   T. Liu et al. MRM 2013;69(2):467-76
%   B. Kressler et al. IEEE TMI 2010;29(2):273-81
%   de Rochefort et al. MRM 2008;60(4):1003-1009
%
%   The coil combination method is similar to
%   MA. Bernstein et al. MRM 1994;32:330-334
%
%   Adapted from a linear fitting created by Ludovic de Rochefort
%   Modified by Tian Liu on 2011.06.01
%   Last modified by Alexey Dimov on 2016.05.12


function [p1, dp1, iFreq_ref, centerF_ref, Y, relres, p0]=Fit_ppm_complex_WorkFlow(M,TE,Mask)
%Modification to handle one echo datasets - assuming zero phase at TE = 0;
%- AD, 05/12/2016
if size(M,4) == 1
    M = cat(4,abs(M),M);
end
%---------- Coil Combination ----------%
if size(M,5)>1
% combine multiple coils together, assuming the coil is the fifth dimension
    M = sum(M.*conj(repmat(M(:,:,:,1,:),[1 1 1 size(M,4) 1])),5);  
    M = sqrt(abs(M)).*exp(1i*angle(M));
end

%M= conj(M);
s0=size(M);
L_s0=length(s0);
nechos=size(M,L_s0);
TE = TE(1:nechos)*1e-6;
TE = TE-TE(1);
M=reshape(M,[prod(s0(1:L_s0-1)),s0(L_s0)]);
s=size(M);
for i = 1:size(TE,2)-1
    deltaTE(i) = TE(i+1)-TE(i);
end

%% for 3 TE object
%format long
% -----Step 1 ------%
Y(:,:) = angle((M(:,1:min(3,nechos))));
d_deltaTE = deltaTE(2) - deltaTE(1);
% ------ Step2 PD images --------%
for i = 1:(min(3,nechos)-1)
    c=(Y(:,i+1)-Y(:,i)); %initial phase difference between echo PD, concern about whether we need 1i here or not
    c(c<-pi) = c(c<-pi)+2*pi;
    c(c>pi) = c(c>pi)-2*pi;
    PD(:,i) = c;
end
% ------ Step 3 DPD images -------%
for i = 1:(min(3,nechos)-2)
    d=(Y(:,i+2)-2*Y(:,i+1)+Y(:,i)); 
    d(d<-pi) = d(d<-pi)+2*pi;
    d(d>pi) = d(d>pi)-2*pi;
    DPD(:,i) = d;
end
% Filter for DPD

nDPD = size(DPD,2);
DPD_tempt = reshape(DPD,[s0(1:L_s0-1),nDPD]);
centerF_ref = DPD_tempt;

% DPD_mf = medfilt3(DPD_tempt,[1 1 2]);
% XZ
DPD_mf = medfilt3(DPD_tempt,[1 1 1]);
DPD = reshape(DPD_mf,[prod(s0(1:L_s0-1)),1]);
% ---- Step 4/Step 5 identify wraps in PD images ------- % 
W1 = DPD/d_deltaTE;
for i = 1:(size(PD,2))
    N(:,i) = fix((PD(:,i)-deltaTE(i)*W1)/2/pi);
end

% -------Step 6/Step 7 picewise unwraping for PD images and normalize for higher SNR ------%
PD_uw = PD-2*pi*N;
W2 = (PD_uw(:,1).*deltaTE(1)+PD_uw(:,2).*deltaTE(2))/(deltaTE(1)+deltaTE(2))/(deltaTE(1));% The weight here need to be optimized
% ------ Step 8/Step9  unwrap for TE images ------%
for i = 1:min(3,nechos)
    n(:,i) = fix((Y(:,i)-TE(i)*W2)/2/pi);
end
Y_uw_u = Y-2*pi*n;
Y_UW = reshape(Y_uw_u,[prod(s0(1:L_s0-1)),3]);
%%  caculate initial p0&p1 from Y
%for i = 1:size(Y_UW,2)
%   Y_UW(:,i) = filloutliers(Y_UW(:,i),'linear','quartiles');
%end
a1=ones(1,min(4,nechos))';
%a1=ones(1,nechos)';
a2=TE(1:min(4,nechos))';
%a2=TE(1:nechos)';
A = [a1,a2];
ip = A(1:min(4,nechos),:)\Y_UW'; % A*ip = Y
%ip = A(1:nechos,:)\Y_uw';
p0 = ip(1,:)';
p1 = ip(2,:)';
% parametr initialization
dp1 = p1;
dp1_Mask = dp1(Mask(:));
tol = norm(p1(Mask(:)))*1*1e-9;
iter = 0;
max_iter = 200;
%% --------  weigthed least square -------- %%
% calculation of WA'*WA
%v1=ones(1,nechos-1);%[1111], 
%v2=reshape(TE(2:nechos),size(v1));
v1=ones(1,nechos);%[1111], 
v2=reshape(TE(1:nechos),size(v1));
a11=sum(abs(M).^2.*(ones(s(1),1)*(v1.^2)),2); %different combination of magnitude map 1+2+3+4
a12=sum(abs(M).^2.*(ones(s(1),1)*(v1.*v2)),2);%2+2*3+3*4
a22=sum(abs(M).^2.*(ones(s(1),1)*(v2.^2)),2);%1*2+4*3+9*4
% inversion
d=a11.*a22-a12.^2;
ai11=a22./d;
ai12=-a12./d;
ai22=a11./d;

while ((norm(dp1_Mask)>tol) &&(iter<max_iter))
    iter = iter+1;
    W2 = abs(M).*exp(1i*(p0*v1 + p1*v2)); % reconstruct W for 4 echoes from abs(M),p0&p1

    % projection
    pr1=sum(conj(1i*W2).*(ones(s(1),1)*v1).*(M-W2),2);
    pr2=sum(conj(1i*W2).*(ones(s(1),1)*v2).*(M-W2),2);

    dp0=real(ai11.*pr1+ai12.*pr2); % error estimation
    dp1=real(ai12.*pr1+ai22.*pr2);
    dp0(isnan(dp0))=0;
    dp1(isnan(dp1))=0;
    dp1_Mask = dp1(Mask(:));
    %update
    p1 = p1+dp1;
    p0 = p0+dp0;
   % a(iter) = norm(dp1_Mask);

end
%figure;plot(a);

% error propagation
dp1=sqrt(ai22);
dp1(isnan(dp1)) = 0;
dp1(isinf(dp1)) = 0;

% relative residual
res = M - abs(M).*exp(1i*(p0*v1 + p1*v2) );
relres = sum(abs(res).^2,2)./sum(abs(M).^2,2);
relres(isnan(relres)) = 0;



p0=reshape(p0,s0(1:L_s0-1)).*deltaTE(1);
p1=reshape(p1,s0(1:L_s0-1)).*deltaTE(1);
dp1=reshape(dp1,s0(1:L_s0-1)).*deltaTE(1);
relres = reshape(relres,s0(1:L_s0-1)).*deltaTE(1);
iFreq_ref = reshape(p1,[(s0(1:L_s0-1)),1]);
% wrap p1 within range -pi to pi
p1(p1>pi)=mod(p1(p1>pi)+pi,2*pi)-pi; 
p1(p1<-pi)=mod(p1(p1<-pi)+pi,2*pi)-pi;
end