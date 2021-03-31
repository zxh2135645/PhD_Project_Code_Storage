resp_time_scale = params.alTR_seconds/phi;
resp_time_scale = resp_time_scale/floor(resp_sr*resp_time_scale);
Hidx=ones(size(nav_data,1),1);
Ridx=Hidx;
wall_clock=ceil((1:numel(Hidx))*params.lEchoSpacing*2/resp_time_scale);
Segidx=repmat((1:(Nseg/2)).',[numel(Hidx)*2/Nseg 1]);

%select region
h=figure;
imshow(dispim(abs(composite)),[]);
title('Select Respiratory ROI')
resp_roi=imrect;
resp_roi=padarray(createMask(resp_roi),[ovs/2 ovs/2]);
close(h)
drawnow

%construct compressed respiratory navdata_temp
sigcov=reshape(SEs.*repmat(composite,[1 1 Ncoils]),[],Ncoils);
sigcov=sigcov(resp_roi(:),:);
sigcov=cov(sigcov);
sigcov=sigcov/trace(sigcov)*Ncoils;

[respmixer,~]=eig(sigcov);
navdata_temp=reshape(reshape(nav_data,[],Ncoils)*sqrtm(sigcov)*respmixer(:,end),size(nav_data,1),[]);

tensor_subspace;

ls_recon;
L_init = L;
Phi_init=Phi;
Phi_rt_init = Phi_rt;
U_init = U;

Phi=reshape(Phi,[L sizes(2:end)]);
respim=Gr\reshape(Phi(:,end,1,1,:),L,[]);
respim=reshape(reshape(U,N^2,[])*respim,N,N,[]);
implay(abs(dispim(respim))/max(abs(vec(dispim(respim)))));