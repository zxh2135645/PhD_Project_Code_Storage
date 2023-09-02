switch Trajectory
  case 'Radial'
    
    delete(gcp('nocreate'));
    parpool('local',8);
    gpuWorkerReset();
    
    [thetast,IA,IC]=unique(thetas(:),'first');
    
    sp.IA = IA;
    sp.IC = IC;
    sp.ParOrder = ParOrder;
    
    om1=cos(thetast*pi/180)*r;
    om2=sin(thetast*pi/180)*r;
    om=[om1(:), om2(:)];
    
    st.useGPU=useGPU;
    
    st.Nd=[N N Nz];
    st.M=size(om,1);
    st.om=om;
    
    clear om1 om2 om
    
    FSU = zeros(trajs,numel(r),Nz,Ncoils,L);
    for traj=1:max(IC) %should this be max(IA)??
      spoke_ind = (IC==traj); %use only current spoke
      Phitemp=Phi_rt(:,spoke_ind)';
      kstemp=kspace_data(spoke_ind,:).';
      parfor np=1:Nz
        %replace (j:trajs:end) with trajs and Nzs
        t_ind = ParOrder(spoke_ind)==np;
        FSU(traj,:,np,:,:)=reshape(kstemp(:,t_ind)*Phitemp(t_ind,:),numel(r),Ncoils,L);
      end
    end
    clear Phitemp kstemp;
    
    FSU=reshape(FSU,st.M,Nz,Ncoils,L);
    Ahb=0;
    for coil=1:Ncoils
      Ahb=Ahb+bsxfun(@times,NUFFT_SoS_adj(FSU(:,:,coil,:),st),conj(SEs(:,:,:,coil)));
    end
    Ahb=Ahb(:);
    clear FSU
    
    U0=zeros(size(Ahb));
    fprintf('Preconditioned conjugate gradient\n')
    tic;
    U0=pcg(@(x)AhA_ps3D(x,st,SEs,Phi_rt,sp),Ahb(:),[],9,M,[],U0(:));
    U=U0;
    toc;
    
  case 'Cartesian'
    FSU = zeros(Ny*Nz,N*Ncoils*L);
    for traj=1:Ny*Nz
      t_ind=find(st.pes==traj);
      FSU(traj,:) = reshape((conj(Phi_rt(:,t_ind))*kspace_data(t_ind,:)).',[],1);
    end
    FSU=reshape(FSU,Ny,Nz,Nx,Ncoils,L);
    Ahb=vec(sum(bsxfun(@times,conj(SEs),reshape(Ah(permute(FSU,[1 3 2 4 5]),st),Ny,Nx,Nz,size(SEs,4),L)),4));
    filtAhb=vec(bsxfun(@rdivide,sum(bsxfun(@times,conj(SEs),reshape(Ah(permute(bsxfun(@times,reshape(st.winv,Ny,Nz),FSU),[1 3 2 4 5]),st),Ny,Nx,Nz,size(SEs,4),L)),4),sqrt(sum(abs(SEs).^2,4))));
    clear FSU
    
    
    %xg 2201117
    filtAhb(isnan(filtAhb)) = 0;
    filtAhb(~isfinite(filtAhb)) = 0;
    
    U0=AhA_ps_cart(filtAhb,st,SEs,Phi_rt,Ny,Nx,Nz);
    U0(isnan(U0))=0;
    U0(~isfinite(U0)) = 0;
    
    %xg
    c2=real(pinv(U0(:))*Ahb(:));
    U0=filtAhb*c2;
    
   
    M = [];
    
    
    fprintf('Preconditioned conjugate gradient\n')
    tic;
    U0=pcg(@(x)AhA_ps_cart(x,st,SEs,Phi_rt,Ny,Nx,Nz),Ahb(:),[],9,M,[],U0(:));
    U=U0;
    toc;
end