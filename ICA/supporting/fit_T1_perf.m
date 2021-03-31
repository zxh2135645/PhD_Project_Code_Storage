cphase=1; %5;

recon = Gr\reshape(Phi(:,:,cphase,1,:),L,[]);
recon=dispim(reshape(reshape(U,Ny*Nx,[])*recon,Ny,Nx,[]));
recon=reshape(recon,N-ovs,N-ovs,Nseg/2,[]);

h=figure;
imshow(abs(mean(recon(:,:,end,:),4)),[])
title('Draw LV')
LV_roi = impoly(gca,'Closed',1);
LV_mask = createMask(LV_roi);
close(h)

h=figure;
imshow(abs(mean(recon(:,:,end,:),4)),[])
title('Draw Myo')
myo_roi = impoly(gca,'Closed',1);
myo_mask = createMask(myo_roi);
close(h)

%%
Nt=sizes(5);
recon=reshape(recon,[],Nseg/2,Nt);
% LV=squeeze(mean(recon(LV_mask,:,:)));
% myo=squeeze(mean(recon(myo_mask,:,:)));
LV=squeeze(mean(recon(LV_mask,:,:)./repmat(sign(recon(LV_mask,end,end)),[1 Nseg/2 Nt])));
myo=squeeze(mean(recon(myo_mask,:,:)./repmat(sign(recon(myo_mask,end,end)),[1 Nseg/2 Nt])));

alpha = params.adFlipAngleDegree*pi/180;

e = @(R1)exp(-params.lEchoSpacing*R1);
Mss = @(e,alpha)(1-e) ./ (1-cos(alpha)*e);

n = 1:2:Nseg;
Sint = @(A,e,alpha,B)A * repmat(Mss(e,alpha),[1 Nseg/2]) .* (1 + repmat(B-1,[1 Nseg/2]).*(repmat(e*cos(alpha),[1 Nseg/2]).^repmat(n-1,[Nt 1]))) * sin(alpha);
S = @(A,R1,alpha,B)Sint(A,e(R1),alpha,B);

ppinv=@(x,y)(y(:)'*x(:))/norm(x(:))^2; %fast left-sided pseudoinverse function

minB = -.5;
maxB = 1;


opts=[];
opts.MaxFunEvals = 10000;

%%
curve=LV.';

normcurve = max(curve(end,:));

Avp = @(R1,alpha,B)ppinv(S(1,R1,alpha,B),curve/normcurve); %parameterize solution to A as function of R1,alpha,B

x0=[ones(1,Nt), alpha/10, zeros(1,Nt)];
xmin=[ones(1,Nt)/3, 0, minB*ones(1,Nt)];
xmax=[ones(1,Nt)*100, alpha, maxB*ones(1,Nt)];
cost=@(x)abs(vec(S(Avp(vec(x(1:Nt)),x(Nt+1),vec(x(Nt+1+(1:Nt)))),vec(x(1:Nt)),x(Nt+1),vec(x(Nt+1+(1:Nt))))-curve/normcurve));
x=lsqnonlin(cost, x0, xmin, xmax, opts);

LVfit.A=Avp(vec(x(1:Nt)),x(Nt+1),vec(x(Nt+1+(1:Nt))))*normcurve;
LVfit.R1=x(1:Nt);
LVfit.alpha=x(Nt+1);
LVfit.B=x(Nt+1+(1:Nt));

%%
curve=myo.';

normcurve = max(curve(end,:));

Avp = @(R1,alpha,B)ppinv(S(1,R1,alpha,B),curve/normcurve); %parameterize solution to A as function of R1,alpha,B

x0=[ones(1,Nt), alpha, zeros(1,Nt)];
xmin=[ones(1,Nt)/3, 0, minB*ones(1,Nt)];
xmax=[ones(1,Nt)*100, alpha, maxB*ones(1,Nt)];
cost=@(x)abs(vec(S(Avp(vec(x(1:Nt)),x(Nt+1),vec(x(Nt+1+(1:Nt)))),vec(x(1:Nt)),x(Nt+1),vec(x(Nt+1+(1:Nt))))-curve/normcurve));
x=lsqnonlin(cost, x0, xmin, xmax, opts);

myofit.A=Avp(vec(x(1:Nt)),x(Nt+1),vec(x(Nt+1+(1:Nt))))*normcurve;
myofit.R1=x(1:Nt);
myofit.alpha=x(Nt+1);
myofit.B=x(Nt+1+(1:Nt));

%%
dt=find(wall_clock==2,1)*params.alTR_seconds/(Nseg/2)/60 % in minutes

avgLV=9;
len=40; %Nt;

l=LVfit.R1(1:len);
LVbase=mean(l(1:avgLV));
m=myofit.R1(1:len);

maxm=max(m);
l=(l-LVbase)/maxm;
m=m/maxm;

fermi=@(f)[zeros(1,f(4)), f(1)./(exp(((1:(2*len))-f(2))/f(3))+1)];
resmax=inf;
for j=len+(1:(2*len))
    [f,res]=lsqnonlin(@(f)conv(l,fermi([f(1:3) j])*dt,'same')-m+f(4),[1 1 1 1/1.3],[0 0 0 .2],[inf inf inf 2]);
    if res < resmax
        fbest=[f(1:3) j];
        myobase=f(4);
        resmax = res;
    end
end
figure,plot(fermi(fbest))
figure,plot(conv(l,fermi(fbest)*dt*maxm,'same'))
hold all
plot(m*maxm-myobase*maxm)
plot(l*maxm)

% max(diff(m)/dt)/max(l)/1.05
ambf=max(fermi(fbest)/1.05)