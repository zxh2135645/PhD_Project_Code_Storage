%% Noise estimation
% temp=log(sum(sum(abs(ifft(nav_data(:,:,1,:),[],2)).^2,1),4));
% [ind, cs]=kmeans(temp.',3);
% [~,cs]=min(cs);
% figure,plot(temp)
% hold all
% plot(max(temp)+(min(temp)-max(temp)).*(ind==cs))
% Psi=ifft(nav_data(:,:,1,:)*sqrt(size(nav_data,2)),[],2);
% Psi=reshape(fft(Psi(:,ind==cs,:),[],2)/sqrt(sum(ind==cs)),[],Ncoils);


Psi=reshape(kspace_data(:,[1 end],:),[],Ncoils);
Psi=cov(Psi);
msdev=sqrt(trace(Psi)/Ncoils);
Psi=Psi/msdev^2;

%% Pre-whiten and transform data
nav_data=fftshift(ifft(ifftshift(nav_data,2),[],2),2)*sqrt(size(nav_data,2));
nav_data=nav_data(:,DC+((-Norig/2):(Norig/2-1)),:,:);

kspace_data(:,:,1,:)=reshape(reshape(kspace_data(:,:,1,:),[],Ncoils)/sqrtm(Psi),size(kspace_data(:,:,1,:)));
nav_data=reshape(reshape(nav_data,[],Ncoils)/sqrtm(Psi),size(nav_data));

msdev=std(vec(kspace_data(:,[1 end],:)));

