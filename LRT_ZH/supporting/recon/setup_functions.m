switch Trajectory
  case 'Radial'
    vec=@(x)x(:);
    prep=@(x,st)reshape(x,st.Nd(1),st.Nd(2),st.Nd(3),[]);%Y x X x Z
    prep_adj=@(x,st)reshape(x,st.M,Nz,[]);
    dcff=@(k)(k+(k==0)*min(k(k~=0))/4);
    dcf=@(st,c)c*dcff(sqrt(sum(abs(st.om).^2,2)));
    crop=@(x,st)x(1:st.Nd(1),1:st.Nd(2),:,:,:);
    dispim=@(x)x(ovs/2 + (1:Norig),ovs/2 + (1:Norig),1,:);
    if strcmp(ScanType,'T2prep')
      dispim=@(x)x(:,:,1,:); %no crop
    end
    if useGPU
      A=@(x,st)NUFFT_SoS(x,st);
      Ah=@(x,st)NUFFT_SoS_adj(prep_adj(x,st),st);
      Ainv=@(x,st,c)Ah(bsxfun(@times,dcf(st,c),prep_adj(x,st)),st);
      AhA=@(x,st)Ah(A(x,st),st);%@(x,st,GhG)vec(st.F'*(st.F*reshape(x,prod(st.Nd),[])));
      A_sense=@(x,st,SEs,Nt)NUFFT_SoS(bsxfun(@times,prep(x,st),SEs),st);
      Ah_sense=@(x,st,SEs,Nt)sum( (NUFFT_SoS_adj(prep_adj(x,st),st).*conj(SEs)) ,4); %sum over coils
      AhA_sense=@(x,st,SEs,GhG,Nt)Ah_sense(A_sense(x,st,SEs,Nt),st,SEs,Nt);
    else
      disp('If you''re trying to do this without the GPU, you''re gonna have a bad time.')
    end
    
    % set up preconditioner
    [kx, ky]=ndgrid(-(N/2):(N/2-1),-(N/2):(N/2-1));
    k = sqrt(kx.^2+ky.^2);
    window = k;
    window(window==0)=1/8*N/Nori.g;
    window=ifftshift(window);
    window(window>(N/2-1))=min(window(:));
    window=window/mean(window(:));
    
    Mf=@(x)ifft2(fft2(x).*repmat(window,[1 1 Nz size(x,4)]));
    M=@(x)vec(Mf(prep(x,st)));
  case 'Cartesian'
    vec=@(x)x(:);
    prep=@(x,st)reshape(x,Ny,Nx,Nz,[]);%Y x X x Z
    prep_adj=@(x,st)reshape(x,Ny,Nx,Nz,[]);
    % ZH
    % dispim=@(x,st)fftshift(x(:,:,1,:),1); %par 1 for now
    dispim=@(x,st)fftshift(x(:,:,4,:),1);
    A=@(x,st)fft(fft(prep(x),[],1),[],3)/sqrt(Ny*Nz);
    Ah=@(x,st)ifft(ifft(prep_adj(x),[],1),[],3)*sqrt(Ny*Nz);
    Ainv=@(x,st,c)Ah(x,st);
    AhA=@(x,st)prep(x,st);
    A_sense=@(x,st,SEs,Nt)fft(fft(bsxfun(@times,prep(x),SEs),[],1),[],3)/sqrt(Ny*Nz);
    Ah_sense=@(x,st,SEs,Nt)sum( ifft(ifft(prep_adj(x),[],1),[],3).*conj(SEs)/sqrt(Ny*Nz) ,4);
    AhA_sense=@(x,st,SEs,GhG,Nt)Ah_sense(A_sense(x,st,SEs,Nt),st,SEs,Nt);
    
    M=@(x)vec(ifft(ifft(bsxfun(@times,fft(fft(prep(x),[],1),[],3),reshape(1./(st.w+1),Ny,1,Nz)),[],1),[],3));
end