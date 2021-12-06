%% Noise estimation
Psi=reshape(kspace_data(abs(ParOrder-DC_kz)>Nz/4,[1 end],:),[],Ncoils); %from readout ends and outer half of partitions
Psi=cov(Psi);
msdev=sqrt(trace(Psi)/Ncoils);
Psi=Psi/msdev^2;

%% Pre-whiten and transform data (truncation and symmetrization)

sizeNav = size(nav_data);
evenEcho_idx = (2:2:sizeNav(1));
oddEcho_idx = (1:2:sizeNav(1));


nav_data_temp(oddEcho_idx,:,:,:) = nav_data(oddEcho_idx,1:(2*(sizeNav(2)-params.BaseResolution)),:,:);
nav_data_temp(evenEcho_idx,:,:,:) = nav_data(evenEcho_idx,2*params.BaseResolution-sizeNav(2)+1:sizeNav(2),:,:);

nav_data = nav_data_temp;
nav_data=fftshift(ifft(ifftshift(nav_data,2),[],2),2)*sqrt(size(nav_data,2));


  
% nav_data_temp(oddEcho_idx,:,:,:) = nav_data(oddEcho_idx,DC+((-Norig/2):(Norig/2-1)),:,:);
% nav_data_temp(evenEcho_idx,:,:,:) = nav_data(evenEcho_idx,sizeNav(2)-DC+((-Norig/2):(Norig/2-1)),:,:);
% nav_data = nav_data_temp;

DC = sizeNav(2)-params.BaseResolution+1
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
    
  sizekspace = size(kspace_data);
%   
  oddEcho_idx_1 = (1:4:sizekspace(1));
  oddEcho_idx_2 = (2:4:sizekspace(1));
  evenEcho_idx_1 = (3:4:sizekspace(1));
  evenEcho_idx_2 = (4:4:sizekspace(1));
%   
  kspace_data_temp(oddEcho_idx_1,:,:,:) = kspace_data(oddEcho_idx_1,1:(2*(sizeNav(2)-params.BaseResolution)),:,:);
  kspace_data_temp(oddEcho_idx_2,:,:,:) = kspace_data(oddEcho_idx_2,1:(2*(sizeNav(2)-params.BaseResolution)),:,:);
  kspace_data_temp(evenEcho_idx_1,:,:,:) = kspace_data(evenEcho_idx_1,2*params.BaseResolution-sizeNav(2)+1:sizeNav(2),:,:);
  kspace_data_temp(evenEcho_idx_2,:,:,:) = kspace_data(evenEcho_idx_2,2*params.BaseResolution-sizeNav(2)+1:sizeNav(2),:,:);
  
  kspace_data = kspace_data_temp;

%   kspace_data = kspace_data(:,1:(2*(sizeNav(2)-params.BaseResolution)),:,:);

  kspace_data=fftshift(ifft(ifftshift(kspace_data,2),[],2),2)*sqrt(size(kspace_data,2));
  kspace_data=kspace_data(:,DC+((-Norig/2):(Norig/2-1)),:,:);
end