ds = round(1/(40*params.lEchoSpacing)); %downsampling to get to 20fps

Phi_rt_temp=sgolayfilt(double(Phi_rt.'),0,5).';

recon=abs(dispim(reshape(reshape(U,Ny*Nx,[])*Phi_rt_temp(:,ceil(ds/2):ds:end),Ny,Nx,[])));
[histn,histx]=hist(recon(:),1000);
cw=histx(find(cumsum(histn)/sum(histn)>.99,1));
implay(abs(recon)/cw);
clear Phi_rt_temp;