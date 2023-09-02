switch Trajectory
  case 'Radial'
    if useGPU
      delete(gcp('nocreate'));
      parpool('local',8);
      gpuWorkerReset();
    end
    
    %backproject
    [~,~,c]=svd([vec(repmat(phantom(N),1,1,st.Nd(3))) vec(Ainv(A(repmat(phantom(N),1,1,st.Nd(3)),st),st,1))],'econ');
    c=real(c(1)/c(2));
    [Phi_rt,~,~]=svde(nav_data(:,:));
    Phi_rt1=Phi_rt(:,1)/mean(Phi_rt(:,1));
    Phi_rt1=vec(repmat(Phi_rt1.',[SGblock 1]));
    
    fbp_data=zeros(trajs,N,Nz,Ncoils);
    for j=1:trajs
      for np=1:Nz
        t_ind = (thetas==thetas(j)) & (ParOrder.'==np);
        fbp_data(j,:,np,:)=reshape((Phi_rt1(t_ind,:)'*kspace_data(t_ind,:)/norm(Phi_rt1(t_ind))^2),[],Ncoils);
      end
    end
    fbp_data(~isfinite(fbp_data))=0;
    
    fbp=Ainv(fbp_data,st,c);
    
    %approx. sampling mask in Cartesian coords
    fbp_window=false(N^2,Nz);
    for j=1:Nz
      t_ind=ParOrder==j;
      kx=round(sin(thetas(t_ind)*pi/180).'*(-N/2:N/2-1))+N/2;
      ky=round(cos(thetas(t_ind)*pi/180).'*(-N/2:N/2-1))+N/2;
      fbp_window(sub2ind([N N],ky((ky~=0) & (kx~=0)),kx((ky~=0) & (kx~=0))),j)=1;
    end
    fbp_window=fftshift(reshape(fbp_window,N,N,Nz));
    
    composite_fbp=sqrt(sum(abs(fbp).^2,4));
    implay(composite_fbp/max(abs(composite_fbp(:))));
    
    % Coil compression for faster computation
    roi_mask=false(Ny,Nx,Nz); %ROI is center region of center slices
    roi_mask(ovs/2+(1:Norig),ovs/2+(1:Norig),[1:Nzorig/2, Nzorig/2+sl_ovs+1:end])=true;
    
    roi=reshape(fbp,[],Ncoils); %region of interest
    roi=roi(roi_mask,:);
    roi=roi'*roi;
    
    roni=reshape(fbp,[],Ncoils); %region of NO interest
    roni=roni(~roi_mask,:);
    roni=roni'*roni;
    
    [cardmixer,cardl]=svde(roi/roni);
    cards=sqrt(cardl);
    
    newcoils=max(find(sqrt(cumsum(diag(cards).^2))/norm(diag(cards))>.99,1),12); %keep 99% of RMS or 12 coils, whichever is larger...
    newcoils=min(newcoils,Ncoils) %...unless there are already <12 coils!
    
    kspace_data=reshape(reshape(kspace_data,[],Ncoils)*cardmixer(:,1:newcoils),size(kspace_data(:,:,:,1:newcoils)));
    nav_data=reshape(reshape(nav_data,[],Ncoils)*cardmixer(:,1:newcoils),size(nav_data(:,:,:,1:newcoils)));
    fbp=reshape(reshape(fbp,[],Ncoils)*cardmixer(:,1:newcoils),size(fbp(:,:,:,1:newcoils)));
    fbp_data=reshape(reshape(fbp_data,[],Ncoils)*cardmixer(:,1:newcoils),size(fbp_data(:,:,:,1:newcoils)));
    Psi=cardmixer(:,1:newcoils)'*Psi*cardmixer(:,1:newcoils);
    Ncoils=newcoils;
    setup_functions; %redo setup
    mixer=eye(Ncoils);
  case 'Cartesian'
    % Spatially-varying coil compression for faster computation
    newcoils=12;
    %%% TODO
    nonavs=true(Nread,1); nonavs(nav_indices)=0;
    
    % Start with Nx/2+1 (center of readout dimension)
    [~,~,CoilV]=svde(squeeze(kspace_data(nonavs,Nx/2+1,:)));
    CoilV=CoilV(:,1:newcoils); %Compression matrix for center of readout dimension
    kspace_data(:,Nx/2+1,1,1:newcoils)=reshape(reshape(kspace_data(:,Nx/2+1,:,:),[],Ncoils)*CoilV,[],1,newcoils); %Replace first NCha coils with compressed version
    nav_data(:,Nx/2+1,1,1:newcoils)=reshape(reshape(nav_data(:,Nx/2+1,:,:),[],Ncoils)*CoilV,[],1,newcoils); %Replace first NCha coils with compressed version
    Psi_new(:,:,Nx/2+1)=CoilV(:,1:newcoils)'*Psi*CoilV(:,1:newcoils); %update noise matrix
    CoilV_orig=CoilV;
    
    %Forward from center
    CoilV_last=CoilV_orig;
    for j=Nx/2+2:Nx
      [~,~,CoilV]=svde(squeeze(kspace_data(nonavs,j,:)));
      CoilV=CoilV(:,1:newcoils);
      [U,~,V]=svd(CoilV_last'*CoilV,'econ');
      CoilV=CoilV*V*U'; %Compression matrix for index j (as close as possible to preceding compressor)
      kspace_data(:,j,1,1:newcoils)=reshape(reshape(kspace_data(:,j,:,:),[],Ncoils)*CoilV,[],1,newcoils); %Replace first NCha coils with compressed version
      nav_data(:,j,1,1:newcoils)=reshape(reshape(nav_data(:,j,:,:),[],Ncoils)*CoilV,[],1,newcoils); %Replace first NCha coils with compressed version
      Psi_new(:,:,j)=CoilV(:,1:newcoils)'*Psi*CoilV(:,1:newcoils); %update noise matrix
      CoilV_last=CoilV;
    end
    
    %Backward from center
    CoilV_last=CoilV_orig;
    for j=Nx/2:-1:1
      [~,~,CoilV]=svde(squeeze(kspace_data(nonavs,j,:)));
      CoilV=CoilV(:,1:newcoils);
      [U,~,V]=svd(CoilV_last'*CoilV,'econ');
      CoilV=CoilV*V*U'; %Compression matrix for index j (as close as possible to preceding compressor)
      kspace_data(:,j,1,1:newcoils)=reshape(reshape(kspace_data(:,j,:,:),[],Ncoils)*CoilV,[],1,newcoils); %Replace first NCha coils with compressed version
      nav_data(:,j,1,1:newcoils)=reshape(reshape(nav_data(:,j,:,:),[],Ncoils)*CoilV,[],1,newcoils); %Replace first NCha coils with compressed version
      Psi_new(:,:,j)=CoilV(:,1:newcoils)'*Psi*CoilV(:,1:newcoils); %update noise matrix
      CoilV_last=CoilV;
    end
    
    %Replace data with compressed version
    kspace_data=kspace_data(:,:,:,1:newcoils); %Keep only compressed coils;
    nav_data=nav_data(:,:,:,1:newcoils);
    Psi=mean(Psi_new,3); %average noise matrix
    Ncoils=newcoils;
    setup_functions; %redo setup
    mixer=eye(Ncoils);
    
    %Calculate first basis image (i.e., a high-SNR static image)
%     % ------ZH: fbp based on all echoes.start
%     [Phi_rt,~,~]=svde(nav_data(:,:));
%     Phi_rt1=Phi_rt(:,1)/mean(Phi_rt(:,1));
% %     Phi_rt1(:)=1; disp('using direct averaging. Change me when you have 4 navs/recovery period');
% %     Phi_rt1=reshape(Phi_rt1,params.NEco,[]);
% %     Phi_rt1(2:end,:)=0;
% %     Phi_rt1=vec(repmat(Phi_rt1(:).',[SGblock 1]));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     Phi_rt1=vec(repmat(Phi_rt1.',[SGblock 1]));
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     fbp_data=zeros(Ny*Nz,Nx,Ncoils);
%     for j=1:Ny*Nz
%       t_ind = (st.pes == j);
%       fbp_data(j,:,:)=reshape((Phi_rt1(t_ind,:)'*kspace_data(t_ind,:)/(norm(Phi_rt1(t_ind))^2+1e-3)),[],Ncoils);
%     end
%     fbp_data(~isfinite(fbp_data))=0;
%     
%     fbp_data=permute(reshape(fbp_data,Ny,Nz,Nx,Ncoils),[1 3 2 4]);
%     fbp=Ainv(fbp_data,st);
%     % ------ZH: fbp based on all echoes.end
    
    % ------ZH: fbp based on 1st echo.start
%     %%%%%%%%%%%%%%%
    [Phi_rt,~,~]=svde(nav_data(1:params.NEco:end,:));
    Phi_rt1=Phi_rt(:,1)/mean(Phi_rt(:,1));
    Phi_rt1=vec(interp1(1:SGblock:Nread/params.NEco,Phi_rt1,1:Nread/params.NEco,'pchip','extrap')); %XG
%     Phi_rt1=vec(interp1(1:SGblock:Nread/params.NEco,Phi_rt1,1:Nread/params.NEco,'linear'));
    Phi_rt1(1:SGblock:end)=0;
    
    fbp_data=zeros(Ny*Nz,Nx,Ncoils);
    s_index=repmat([1:SGblock:size(kspace_data,1)],[SGblock 1])+repmat([0:SGblock-1].',[1 Nread/SGblock]);
    ss_index=vec(s_index(:,1:params.NEco:end));
    pes=st.pes(ss_index,:);
    kspace_data_small=kspace_data(ss_index,:,:,:);
%     pes = st.pes(1:8:end,:);
%     kspace_data_small=kspace_data(1:8:end,:,:,:);
    for j=1:Ny*Nz
        t_ind=(pes==j);
        fbp_data(j,:,:)=reshape((Phi_rt1(t_ind,:)'*kspace_data_small(t_ind,:)/(norm(Phi_rt1(t_ind))^2+1e-3)),[],Ncoils);
    end
    fbp_data(~isfinite(fbp_data))=0;
    
    fbp_data=permute(reshape(fbp_data,Ny,Nz,Nx,Ncoils),[1 3 2 4]);
    fbp=Ainv(fbp_data,st);
    clear kspace_data_small
%     %%%%%%%%%%%%%%%%%
    % ------ZH: fbp based on 1st echo.end
   
%     %SNR-weighted sampling mask in Cartesian coords
%     fbp_window=zeros(Ny*Nz,Nx);
%     for j=1:Ny*Nz
%          t_ind = (pes == j);
%          fbp_window(j,:)=sum(abs(Phi_rt1(t_ind,:)).^2)/(norm(Phi_rt1(t_ind))^2+1e-3)^2;
%     end
%     fbp_window=1./fbp_window;
%     fbp_window(~isfinite(fbp_window))=0;
%     fbp_window=permute(reshape(fbp_window,Ny,Nz,Nx),[1 3 2]);
%     fbp_window=fftshift(sqrt(abs(fbp_window)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     [Phi_rt,~,~] = svde(nav_data(1:params.NEco:end,:));
%     Phi_rt1 = Phi_rt(:,1)/mean(Phi_rt(:,1));
%     Phi_rt1 = interp1(1:SGblock:Nread/params.NEco,Phi_rt1,1:Nread/params.NEco,'pchip',0).';
% %     Phi_rt1(end) = Phi_rt1(end-1);
%     Phi_rt1(1:SGblock:end) = 0;
%     
%     kspace_data_fbp = kspace_data(1:params.NEco:end,:,:,:);
% %     s_index_1 = 1:SGblock:size(kspace_data_fbp,1);
% %     s_index_2 = s_index_1+1;
% %     s_index = [s_index_1;s_index_2];
% %     s_index = s_index(:)';
% %     kspace_data_fbp = kspace_data_fbp(s_index,:,:,:);
%     pes_fbp = st.pes(1:params.NEco:end);
% %     pes_fbp = pes_fbp(s_index);
%     
%     fbp_data=zeros(Ny*Nz,Nx,Ncoils);
%     for j=1:Ny*Nz
%       t_ind = (pes_fbp == j);
%       fbp_data(j,:,:)=reshape((Phi_rt1(t_ind,:)'*kspace_data_fbp(t_ind,:)/(norm(Phi_rt1(t_ind))^2+1e-3)),[],Ncoils);
%     end
%     clear kspace_data_fbp;
%     fbp_data(~isfinite(fbp_data))=0;
%     fbp_data=permute(reshape(fbp_data,Ny,Nz,Nx,Ncoils),[1 3 2 4]);
%     fbp=Ainv(fbp_data,st);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    fbp_window=fftshift(logical(repmat(reshape(st.w,Ny,1,Nz),[1 Nx 1])));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    composite_fbp=sqrt(sum(abs(fbp).^2,4));
    implay(composite_fbp/max(abs(composite_fbp(:))));
end

%% estimate SEs
if iscartesian
  SEs=squeeze(ifftshift(ifftshift(sensemaps_cbd(bsxfun(@times,bsxfun(@times,fftshift(fftshift(fbp,3),1),reshape((-1).^(1:Nz),1,1,[])),((-1).^(1:Ny).')*((-1).^(1:Nx))),[7 7 5],'figures','acs',fbp_window),3),1)); %slightly slower, but more accurate when fbp has significant streaking
else
  SEs=squeeze(ifftshift(sensemaps_cbd(bsxfun(@times,fftshift(fbp,3),reshape((-1).^(1:Nz),1,1,[])),'figures','acs',fbp_window),3)); %slightly slower, but more accurate when fbp has significant streaking
end

SE_corr=reshape(sqrt(sum(abs(reshape(SEs,[],Ncoils)*sqrtm(Psi)).^2,2)),Ny,Nx,Nz);
SEs=bsxfun(@rdivide,SEs,SE_corr);
% 
% % %Bias correction
% % pd=abs(sum(conj(SEs).*fbp,4)./sum(abs(SEs).^2,4));
% % pd(pd==0)=min(pd(pd~=0));
% % if iscartesian
% %   [x,y,z]=ndgrid(ifftshift(-Ny/2:Ny/2-1),-Nx/2:Nx/2-1,ifftshift(-Nz/2:Nz/2-1));
% % else
% %   [x,y,z]=ndgrid(-Ny/2:Ny/2-1,-Nx/2:Nx/2-1,ifftshift(-Nz/2:Nz/2-1));
% % end
% % polyterm=[x(:).^2, y(:).^2, z(:).^2, x(:).*y(:), x(:).*z(:), y(:).*z(:), x(:), y(:), z(:), ones(size(x(:)))];
% % [polyterm,~,~]=svd(bsxfun(@times,polyterm,pd(:).^2),'econ');
% % w=polyterm'*pd(:);
% % %   norm(polyterm*w-pd(:),1)
% % y=zeros(size(pd(:)));
% % rho=1/max(abs(polyterm*w-pd(:)));
% % for it=1:10
% %   res=polyterm*w-pd(:);
% %   z=res+y/rho;
% %   z=sign(z).*max(abs(z)-1/rho,0);
% %   y=y+rho*(res-z);
% %   rho=rho*1.25;
% %   w=polyterm'*(z+pd(:)-y/rho);
% %   %     norm(polyterm*w-pd(:),1)
% % end
% % [~,~,v]=svd([(polyterm*w)./pd(:), pd(:)],'econ');
% % w=w*v(2)/v(1);
% % bias=reshape(polyterm*w,Ny,Nx,Nz).\pd.^2;
% % bias(1./bias<1/2)=2;
% % bias(bias<1/4)=1/4;
% % implay(abs([pd pd./bias])/max(abs(pd(:)./bias(:))))
% % SEs=bsxfun(@times,SEs,bias);
%% bart SEs

% bart_path = '/home/guanx1/bart';
% addpath(fullfile(bart_path,'matlab'));
% setenv('TOOLBOX_PATH',bart_path);
% 
% %%
% % fbps = bsxfun(@times,bsxfun(@times,fftshift(fftshift(fbp,3),1),reshape((-1).^(1:Nz),1,1,[])),((-1).^(1:Ny).')*((-1).^(1:Nx)));
% fbps = bsxfun(@times,bsxfun(@times,fftshift(fftshift(fbp,3),1),reshape((-1).^(1:Nz),1,1,[])),((-1).^(1:Ny).')*((-1).^(1:Nx)));
% fbp_k = fft(fft(fft(fbps,[],1),[],2),[],3)/sqrt(Ny*Nz*Nx);
% tic;
% [SEs,emaps] = bart('ecalib -m 1 -r 30,30,10 -k 5,5,3 -c 0.8',fbp_k);
% % [SEs,emaps] = bart('ecalib -m 1 -c 0.8',fbp_k);
% % [SEs,emaps] = bart('ecalib',fbp_k);
% toc;
% % SEs = ifftshift(ifftshift(SEs,3),1);
% SEs = ifftshift(ifftshift(SEs,3),2);
% % emaps = ifftshift(ifftshift(emaps,3),1);
% clear fbp_k fbps
% 
% % SE_corr=reshape(sqrt(sum(abs(reshape(SEs,[],Ncoils)*sqrtm(Psi)).^2,2)),Ny,Nx,Nz);
% % SEs=bsxfun(@rdivide,SEs,SE_corr);
% % SEs(isnan(SEs))=0;


%% SENSE on composite
if ~iscartesian
  SEs(:,:,2:2:end,:)=-SEs(:,:,2:2:end,:); %fix alternating sign in partition direction
end
composite_init=sum(conj(SEs).*fbp,4)./sum(abs(SEs).^2,4);
composite_init(isnan(composite_init))=0;
composite=prep(pcg(@(x)vec(AhA_sense(x,st,SEs,GhG,1)),...   
  vec(Ah_sense(fbp_data,st,SEs,1)),[],20,M,[],composite_init(:)),st);
figure,imshow(abs([dispim(composite_fbp) dispim(composite_init) dispim(prep(composite,st))]),[])
title('Composite Reconstructions'),drawnow;

switch Trajectory
  case 'Radial'
    clear st GhG
  case 'Cartesian'
    st.Ahmat=[];
end
