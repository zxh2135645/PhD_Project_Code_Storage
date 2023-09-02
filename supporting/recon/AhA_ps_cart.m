function AhAx = AhA_ps_cart(x,st,SEs,Phi,Ny,Nx,Nz)

L = size(Phi,1);
Ncoils = size(SEs,4);

% try
FSU=fft(fft(bsxfun(@times,reshape(x,Ny,Nx,Nz,1,[]),SEs),[],1),[],3);
%     SFU=diag_sp(shifter)*fft(repmat(reshape(x,Ny,Nx,[]),[1 1 1 Ncoils]).*repmat(reshape(SEs,size(SEs,1),size(SEs,2),1,[]),[1 1 L 1]),Nfft);
FSU=reshape(permute(FSU,[1 3 2 4 5]),Ny*Nz,[]);
for traj=1:Ny*Nz
    t_ind=(st.pes==traj);
    FSU(traj,:) = reshape(reshape(FSU(traj,:),[],L)*(Phi(:,t_ind)*Phi(:,t_ind)'),1,[]);
end
AhAx=reshape(sum(bsxfun(@times,ifft(ifft(permute(reshape(FSU,Ny,Nz,Nx,Ncoils,[]),[1 3 2 4 5]),[],1),[],3),conj(SEs)),4),[],1);

%     AhAx=reshape(sum(repmat(reshape(conj(SEs),size(SEs,1),size(SEs,2),1,[]),[1 1 L 1]) ...
%         .*reshape(crop(ifft(diag_sp(conj(shifter))*SFU)),Ny,Nx,L,[]),4),[],1);
% catch
%     Phi2=zeros(L,L,max(IA));
%     for traj=1:max(IA)
%         t_ind=(IC==traj);
%         Phi2(:,:,traj) = Phi(:,t_ind)*Phi(:,t_ind)';
%     end
%
%     AhAx = zeros(numel(x),1);
%     for coil=1:Ncoils
%         FU=diag_sp(shifter)*fft(reshape(x,Ny,Nx,[]).*repmat(SEs(:,:,coil),[1 1 L]),Nfft);
%         for traj=1:Nfft
%             t_ind=(st.pes==traj);
%             FU(traj,:) = reshape(reshape(FU(traj,:),[],L)*Phi2(:,:,traj),1,[]);
%         end
%         AhAx=AhAx + reshape(repmat(conj(SEs(:,:,coil)),[1 1 L]) ...
%             .*reshape(crop(ifft(diag_sp(conj(shifter))*FU)),Ny,Nx,[]),[],1);
%     end
% end

fprintf('.')

end