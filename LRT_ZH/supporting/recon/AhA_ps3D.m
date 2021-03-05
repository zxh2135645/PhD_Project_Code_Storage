function AhAx = AhA_ps3D(x,st,SEs,Phi,sampling_pattern)
try
  L = size(Phi,1);
  Ncoils = size(SEs,4);
  
  IA=sampling_pattern.IA;
  IC=sampling_pattern.IC;
  ParOrder = sampling_pattern.ParOrder;
  
  prep=@(x,st)reshape(x,st.Nd(1),st.Nd(2),st.Nd(3),[],L);
  
  % Phi2=zeros(L,L,max(IA));
  % for traj=1:max(IA)
  %   t_ind=(IC==traj);
  %   Phi2(:,:,traj) = Phi(:,t_ind)*Phi(:,t_ind)';
  % end
  
%   try
%     FSU=NUFFT_SoS(bsxfun(@times,prep(x,st),SEs),st);
%     FSU=reshape(FSU,numel(IA),st.M/numel(IA),st.Nd(3),[]);
%     for traj=1:max(IC)
%       for np=1:st.Nd(3)
%         t_ind = (IC==traj) & (ParOrder==np);
%         FSU(traj,:,np,:) = reshape(reshape(FSU(traj,:,np,:),[],L)*(Phi(:,t_ind)*Phi(:,t_ind)'),[],Ncoils*L);
%       end
%     end
%     FSU=reshape(FSU,st.M,st.Nd(3),Ncoils,L);
%     AhAx=sum(bsxfun(@times,NUFFT_SoS_adj(FSU,st),conj(SEs)),4);
%     AhAx=AhAx(:);
%   catch
    AhAx = zeros(numel(x),1);
    for coil=1:Ncoils
      FU=NUFFT_SoS(bsxfun(@times,prep(x,st),SEs(:,:,:,coil)),st);
      FU=reshape(FU,numel(IA),st.M/numel(IA),st.Nd(3),[]);
      for traj=1:max(IC)
        for np=1:st.Nd(3)
          t_ind = (IC==traj) & (ParOrder==np);
          FU(traj,:,np,:) = reshape(FU(traj,:,np,:),[],L)*(Phi(:,t_ind)*Phi(:,t_ind)');
        end
      end
      FU=reshape(FU,st.M,st.Nd(3),1,L);
      AhAx=AhAx + reshape(bsxfun(@times,NUFFT_SoS_adj(FU,st),conj(SEs(:,:,:,coil))),[],1);
    end
%   end
  
  fprintf('.')
catch errormsg
  errormsg
  pause(rand*30);
  delete(gcp('nocreate')); %delete current parpool
  parpool('local',8); %open parpool with 8 workers
  gpuWorkerReset(); %reassign workers
  AhAx = AhA_ps3D(x,st,SEs,Phi,sampling_pattern);
end
end
