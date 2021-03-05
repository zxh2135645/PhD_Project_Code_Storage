
ds = round(1/(params.lEchoSpacing)); %downsampling to get to 20fps
recon=reshape(U,Ny,Nx,Nz,L);
recon=dispim(recon); %dispim uses 1st partition
recon=reshape(recon,[],L);
if strcmp(ScanType,'T2prep')
  recon=abs(reshape(recon*Phi_rt_full(:,ceil(ds/2):ds:min(Nread,2000)),N,N,[]));
elseif iscartesian
%   recon=(reshape(recon*Phi_rt_full(:,ceil(ds/2):ds:min(Nread,2000)),Ny,Nx,[]));
  recon=(reshape(recon*Phi_rt_full(:,3:16:4000),Ny,Nx,[]));
else
  recon=abs(reshape(recon*Phi_rt_full(:,ceil(ds/2):ds:min(Nread,2000)),Norig,Norig,[]));
end
       

cw=max(recon(:));

implay(abs(recon)/abs(cw));

% implay(abs(recon(:,:,1:5:end))/cw);
% implay(abs(recon(:,:,2:5:21475))/cw);
% implay(abs(recon(:,:,3:5:21475))/cw);
% implay(abs(recon(:,:,4:5:21475))/cw);
% implay(abs(recon(:,:,5:5:21475))/cw);

% precon=mean(recon(60:80,55:95,601:5:1600),3);
% postcon=mean(recon(60:80,55:95,5001:5:10000),3);
% figure,imshow([abs(precon), abs(postcon), (abs(postcon)-abs(precon))]/cw);
