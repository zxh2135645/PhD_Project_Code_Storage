
[thetast,IA,IC]=unique(thetas(:),'first');

sp.IA = IA;
sp.IC = IC;

om1=cos(thetast*pi/180)*r;
om2=sin(thetast*pi/180)*r;
om=[om1(:), om2(:)];

if useGPU
  st.Nd=[N N];
  st.M=size(om,1);
  st.om=om;
  %           st.F=gpuNUFFT(om.'/(2*pi),ones(1,st.M),2,6,8,st.Nd,[],true);
  %           st.FS=gpuNUFFT(om.'/(2*pi),ones(1,st.M),2,6,8,st.Nd,SEs,true);
else
  st=nufft_init(om,[N N],[6 6],2*[N N],[N N]/2);
end
st.useGPU=useGPU;
clear om1 om2 om

FSU = zeros(trajs,numel(r)*Ncoils*L);
for traj=1:max(IA)
  t_ind=(IC==traj);
  FSU(traj,:) = reshape((conj(Phi_rt(:,t_ind))*kspace_data(t_ind,:)).',[],1);
end

if useGPU
  FSU=reshape(FSU,st.M,Ncoils,L);
  Ahb=complex(zeros(prod(st.Nd),L));
  filtAhb=Ahb;
  parfor (j=1:L, 6)
    FS=gpuNUFFT(st.om.'/(2*pi),ones(1,st.M),2,6,8,st.Nd,SEs,true);
    Ahb(:,j)=vec(FS'*FSU(:,:,j));
    filtAhb(:,j)=vec(FS'*(diag_sp(dcf(st,c))*FSU(:,:,j)))./vec(sum(abs(SEs).^2,3));
  end
  Ahb=Ahb(:);
  filtAhb=filtAhb(:);
  filtAhb(~isfinite(filtAhb))=0;
else
  Ahb=vec(sum(repmat(conj(SEs),[1 1 1 L]).*reshape(nufft_adj(prep_adj(FSU,st),st),st.Nd(1),st.Nd(2),size(SEs,3),L),3));
  filtAhb=vec(bsxfun(@ldivide,sum(abs(SEs).^2,3),sum(repmat(conj(SEs),[1 1 1 L]).*reshape(nufft_adj(diag_sp(dcf(st,c))*prep_adj(FSU,st),st),st.Nd(1),st.Nd(2),size(SEs,3),L),3)));
  filtAhb(~isfinite(filtAhb))=0;
end
clear FSU

U0=AhA_ps(filtAhb,st,SEs,Phi_rt,sp);
c2=real(pinv(U0(:))*Ahb(:));
U0=filtAhb*c2;

fprintf('Preconditioned conjugate gradient\n')
tic;
U0=pcg(@(x)AhA_ps(x,st,SEs,Phi_rt,sp),Ahb(:),[],9,M,[],U0(:));
U=U0;
toc;

