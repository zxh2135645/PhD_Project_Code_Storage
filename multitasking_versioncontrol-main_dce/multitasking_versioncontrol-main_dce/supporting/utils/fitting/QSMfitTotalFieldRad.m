function [p1, dp1, relres, p0, iter] = QSMfitTotalFieldRad(reconME,TEarray)

if size(reconME,5)>1
% combine multiple coils together, assuming the coil is the fifth dimension
    reconME = sum(reconME.*conj(reconME(:,:,:,1,:)),5);  
    reconME = sqrt(abs(reconME)).*exp(1i*angle(reconME));
end

reconME = conj(reconME);
[Ny,Nx,Nz,Necho] = size(reconME);

reconME = reshape(reconME,[],Necho);

phaseME = angle(reconME(:,1:min(3,Necho)));
TEarray = TEarray(1:min(3,Necho));
Necho = min(3,Necho);

% % estimate the slope
% c = ((Y(:,2)-Y(:,1)));
% [m ind]=min([abs(c-2*pi),abs(c),abs(c+2*pi)],[],2);
% c(ind==1)=c(ind==1)-2*pi;
% c(ind==3)=c(ind==3)+2*pi;

% unwrap the second echo
% cd = Y(:,2)-Y(:,1)-c;
% Y(cd<-pi,2)=Y(cd<-pi,2)+2*pi;
% Y(cd>pi,2)=Y(cd>pi,2)-2*pi;

phaseIncE2 = phaseME(:,2) - phaseME(:,1);
phaseIncE2(phaseIncE2> pi) = phaseIncE2(phaseIncE2> pi) - 2*pi;
phaseIncE2(phaseIncE2<-pi) = phaseIncE2(phaseIncE2<-pi) + 2*pi;

phaseME(:,2) = phaseME(:,1) + phaseIncE2;

% unwrap the third echo
% cd=(Y(:,3)-Y(:,2))-(TE(3)-TE(2))/(TE(2)-TE(1))*c;
% cd_minus=(cd<-pi);
% Y(:,3)=Y(:,3)+cd_minus.*abs(fix((cd-pi)./(2*pi)))*2*pi;
% cd_plus=(cd>pi);
% Y(:,3)=Y(:,3)-cd_plus.*fix((cd+pi)./(2*pi))*2*pi;

if Necho == 3
    phaseIncE3 = phaseME(:,3) - phaseME(:,1);
    phaseIncE3(phaseIncE3> pi) = phaseIncE2(phaseIncE3> pi) - 2*pi;
    phaseIncE3(phaseIncE3<-pi) = phaseIncE2(phaseIncE3<-pi) + 2*pi;
    
    phaseTemp  = phaseIncE2*(TEarray(3)-TEarray(1))/(TEarray(2)-TEarray(1));
    phaseShift = cos(phaseTemp - phaseIncE3 - pi) > cos(phaseTemp - phaseIncE3);
    phaseIncE3(phaseShift) = phaseIncE3(phaseShift) + sign(phaseIncE2(phaseShift))*2*pi;

    phaseME(:,3) = phaseME(:,1) + phaseIncE3;
end

A = [ones(Necho,1) TEarray(:)-TEarray(1)];
ip = A\(phaseME(:,1:3).'-phaseME(:,1).');
p0 = ip(1,:)';
p1 = ip(2,:)';

dp1 = p1;
tol = norm(p1(:))*1e-4;
iter = 0;
max_iter = 30;

reconME = reconME.*exp(-1j*phaseME(:,1));

% weigthed least square
% calculation of WA'*WA
v1 = ones(1,Necho);
v2 = reshape(TEarray-TEarray(1),size(v1));
a11 = sum(abs(reconME).^2.*(ones(Ny*Nx*Nz,1)*(v1.^2)),2);
a12 = sum(abs(reconME).^2.*(ones(Ny*Nx*Nz,1)*(v1.*v2)),2);
a22 = sum(abs(reconME).^2.*(ones(Ny*Nx*Nz,1)*(v2.^2)),2);
% inversion
d = a11.*a22-a12.^2;
ai11 =  a22./d;
ai12 = -a12./d;
ai22 =  a11./d;

while ((norm(dp1)>tol) &&(iter<max_iter))
    iter = iter+1;
    W = abs(reconME).*exp(1i*(p0*v1 + p1*v2) );

    % projection
    pr1 = sum(conj(1i*W).*(ones(Ny*Nx*Nz,1)*v1).*(reconME-W),2);
    pr2 = sum(conj(1i*W).*(ones(Ny*Nx*Nz,1)*v2).*(reconME-W),2);

    dp0 = real(ai11.*pr1+ai12.*pr2);
    dp1 = real(ai12.*pr1+ai22.*pr2);
    dp1(isnan(dp1))=0;
    dp0(isnan(dp0))=0;
    
    % update
    p1 = p1+dp1;
    p0 = p0+dp0;
end

% error propagation
dp1 = sqrt(ai22);
dp1(isnan(dp1)) = 0;
dp1(isinf(dp1)) = 0;

% relative residual
res = reconME - abs(reconME).*exp(1i*(p0*v1 + p1*v2) );
relres = sum(abs(res).^2,2)./sum(abs(reconME).^2,2);
relres(isnan(relres)) = 0;

p0  = reshape(p0,[Ny Nx Nz])*(TEarray(2)-TEarray(1));
p1  = reshape(p1,[Ny Nx Nz])*(TEarray(2)-TEarray(1));
dp1 = reshape(dp1,[Ny Nx Nz])*(TEarray(2)-TEarray(1));
relres = reshape(relres,[Ny Nx Nz])*(TEarray(2)-TEarray(1));
p1(p1>pi)  = mod(p1(p1> pi)+pi,2*pi)-pi;
p1(p1<-pi) = mod(p1(p1<-pi)+pi,2*pi)-pi;
    

