function AhAx = AhA_ps(x,st,SEs,Phi,sampling_pattern)

L = size(Phi,1);
Ncoils = size(SEs,3);

IA=sampling_pattern.IA;
IC=sampling_pattern.IC;

vec=@(x)x(:);

Phi2=zeros(L,L,max(IA));
for traj=1:max(IA)
  t_ind=(IC==traj);
  Phi2(:,:,traj) = Phi(:,t_ind)*Phi(:,t_ind)';
end

if st.useGPU
  [~,d]=unix('nvidia-smi | grep MiB | head -1 | awk ''{print $9, $11}''');
  d=textscan(d,'%d%s%d%s');
  am=d{3}-d{1};
  if am < 2e3
    delete(gcp('nocreate'));
    parpool('local',6);
  end
  x=reshape(x,[],L);
  AhAx = complex(zeros(size(x)));
  FSU=zeros(st.M*Ncoils,L);
  parfor (j=1:L, 6)
    FS=gpuNUFFT(st.om.'/(2*pi),ones(1,st.M),2,6,8,st.Nd,SEs,true);
    FSU(:,j)=vec(FS*x(:,j));
  end
  FSU=reshape(FSU,numel(IA),[]);
  for traj=1:max(IA)
    FSU(traj,:) = reshape(reshape(FSU(traj,:),[],L)*Phi2(:,:,traj),1,[]);
  end
  FSU=reshape(FSU,st.M,Ncoils,L);
  parfor (j=1:L, 6)
    FS=gpuNUFFT(st.om.'/(2*pi),ones(1,st.M),2,6,8,st.Nd,SEs,true);
    AhAx(:,j)=vec(FS'*FSU(:,:,j));
  end
  AhAx=AhAx(:);
else
  AhAx = zeros(numel(x),1);
  for coil=1:Ncoils
    FU=reshape(nufft(reshape(x,st.Nd(1),st.Nd(2),[]).*repmat(SEs(:,:,coil),[1 1 L]),st),numel(IA),[]);
    for traj=1:max(IA)
      FU(traj,:) = reshape(reshape(FU(traj,:),[],L)*Phi2(:,:,traj),1,[]);
    end
    AhAx=AhAx + reshape(repmat(conj(SEs(:,:,coil)),[1 1 L]) ...
      .*reshape(nufft_adj(reshape(FU,st.M,[]),st),st.Nd(1),st.Nd(2),[]),[],1);
  end
end

fprintf('.')

end