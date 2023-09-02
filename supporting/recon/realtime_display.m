
ds = round(1/(params.lEchoSpacing)); %downsampling to get to 20fps
recon=reshape(U,Ny,Nx,Nz,L);
recon=dispim(recon); %dispim uses 1st partition
recon=reshape(recon,[],L);
%recon=(reshape(recon*Phi_rt_full(:,1:8:4000),Ny,Nx,[]));

if strcmp(ScanType,'T2prep')
  recon=abs(reshape(recon*Phi_rt_full(:,ceil(ds/2):ds:min(Nread,2000)),N,N,[]));
elseif iscartesian
%   recon=(reshape(recon*Phi_rt_full(:,ceil(ds/2):ds:min(Nread,2000)),Ny,Nx,[]));
  recon=(reshape(recon*Phi_rt_full(:,1:1:400),Ny,Nx,[])); % XZ
else
  recon=abs(reshape(recon*Phi_rt_full(:,ceil(ds/2):ds:min(Nread,2000)),Norig,Norig,[]));
end

%recon=reshape(U_init,Ny,Nx,Nz,L_init);
%recon=dispim(recon); %dispim uses 1st partition
%recon=reshape(recon,[],L_init);
%recon=(reshape(recon*Phi_rt_full_init(:,3:10:4000),Ny,Nx,[]));

cw=max(recon(:));
implay(abs(recon)/abs(cw)*2);

% implay(abs(recon(:,:,1:5:end))/cw);
% implay(abs(recon(:,:,2:5:21475))/cw);
% implay(abs(recon(:,:,3:5:21475))/cw);
% implay(abs(recon(:,:,4:5:21475))/cw);
% implay(abs(recon(:,:,5:5:21475))/cw);

% precon=mean(recon(60:80,55:95,601:5:1600),3);
% postcon=mean(recon(60:80,55:95,5001:5:10000),3);
% figure,imshow([abs(precon), abs(postcon), (abs(postcon)-abs(precon))]/cw);
