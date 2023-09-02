%% Noise estimation
Psi=reshape(kspace_data(abs(ParOrder-DC_kz)>Nz/4,[1 end],:),[],Ncoils); %from readout ends and outer half of partitions
Psi=cov(Psi);
msdev=sqrt(trace(Psi)/Ncoils);
Psi=Psi/msdev^2;

%% Pre-whiten and transform data
nav_data=fftshift(ifft(ifftshift(nav_data,2),[],2),2)*sqrt(size(nav_data,2));
nav_data=nav_data(:,DC+((-Norig/2):(Norig/2-1)),:,:);

try
    kspace_data(:,:,1,:)=reshape(reshape(kspace_data(:,:,1,:),[],Ncoils)/sqrtm(Psi),size(kspace_data(:,:,1,:)));
catch
    kspace_data(:,:,1,:)=reshape(reshape(kspace_data(:,:,1,:),[],Ncoils)*inv(sqrtm(Psi)),size(kspace_data(:,:,1,:)));
end

try
    nav_data=reshape(reshape(nav_data,[],Ncoils)/sqrtm(Psi),size(nav_data));
catch
    nav_data=reshape(reshape(nav_data,[],Ncoils)*inv(sqrtm(Psi)),size(nav_data));
end

msdev=std(vec(kspace_data(abs(ParOrder-DC_kz)>Nz/4,[1 end],:)));
if iscartesian %Correct for readout oversampling
  kspace_data=fftshift(ifft(ifftshift(kspace_data,2),[],2),2)*sqrt(size(kspace_data,2));
  kspace_data=kspace_data(:,DC+((-Norig/2):(Norig/2-1)),:,:);
  
  
end
