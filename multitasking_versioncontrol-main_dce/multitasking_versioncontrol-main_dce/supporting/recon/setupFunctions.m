vec = @(x) x(:);
    
switch Trajectory
  case {'Radial','Spiral'}
    prep     = @(x,st) reshape(x, st.Nd(1), st.Nd(2), Nz, []);  
    prep_adj = @(x,st) reshape(x, st.M, Nz, []);
    dcff     = @(k)    (k+(k==0)*min(k(k~=0))/4);
    dcf      = @(st,c) c*dcff(sqrt(sum(abs(st.om).^2,2)));
    crop     = @(x,st) x(1:st.Nd(1),1:st.Nd(2),:,:,:);
    dispim   = @(x)    x(floor(st.Nd(1)/2)-floor(Nydisp/2) + (1:Nydisp), floor(st.Nd(2)/2)-floor(Nxdisp/2) + (1:Nxdisp), 1, :);
    if strcmp(params.ScanType,'T2prep')
        dispim = @(x) x(:,:,ceil(Nz/2),:); %no crop
    end
    if flagUseBart && isPreparedbartnufft
        A         = @(x,st) NUFFT_SoS_bart(x,st);
        Ah        = @(x,st) NUFFT_SoS_adj_bart(x,st);
        Ainv      = @(x,st,c) NUFFT_SoS_adj_bart(x,st,1);
        AhA       = @(x,st)   Ah(A(x,st),st);  
        A_sense   = @(x,st,SEs,Nt)     NUFFT_SoS_bart(bsxfun(@times,prep(x,st),SEs),st);
        Ah_sense  = @(x,st,SEs,Nt)     sum((NUFFT_SoS_adj_bart(prep_adj(x,st),st).*conj(SEs)) ,4); %sum over coils
        AhA_sense = @(x,st,SEs,GhG,Nt) Ah_sense(A_sense(x,st,SEs,Nt),st,SEs,Nt);
    elseif flagUsefinufft && flagUseGPU && isPreparedcufinufft
        A         = @(x,st)   finufft_SoS(prep(x,st), st);
        Ah        = @(x,st)   cufinufft_SoS_adj(prep_adj(x,st), st);
        Ainv      = @(x,st,c) Ah(bsxfun(@times,dcf(st,c),prep_adj(x,st)), st);
        AhA       = @(x,st)   Ah(A(x,st),st);  
        A_sense   = @(x,st,SEs,Nt)     finufft_SoS(bsxfun(@times,prep(x,st),SEs),st);
        Ah_sense  = @(x,st,SEs,Nt)     sum(prep(finufft_SoS_adj(prep_adj(x,st),st),st).*conj(SEs), 4); %sum over coils
        AhA_sense = @(x,st,SEs,GhG,Nt) Ah_sense(A_sense(x,st,SEs,Nt),st,SEs,Nt);
    elseif flagUsefinufft && isPreparedfinufft
        reconOptions.flagUseGPU = 0;
        A         = @(x,st)   finufft_SoS(prep(x,st), st);
        Ah        = @(x,st)   finufft_SoS_adj(prep_adj(x,st), st);
        Ainv      = @(x,st,c) Ah(bsxfun(@times,dcf(st,c),prep_adj(x,st)), st);
        AhA       = @(x,st)   Ah(A(x,st),st);  
        A_sense   = @(x,st,SEs,Nt)     finufft_SoS(bsxfun(@times,prep(x,st),SEs),st);
        Ah_sense  = @(x,st,SEs,Nt)     sum(prep(finufft_SoS_adj(prep_adj(x,st),st),st).*conj(SEs), 4); %sum over coils
        AhA_sense = @(x,st,SEs,GhG,Nt) Ah_sense(A_sense(x,st,SEs,Nt),st,SEs,Nt);
    elseif flagUseGPU && isPreparedgpuNUFFT
        A         = @(x,st) NUFFT_SoS(x,st);
        Ah        = @(x,st) NUFFT_SoS_adj(prep_adj(x,st),st);
        Ainv      = @(x,st,c) Ah(bsxfun(@times,dcf(st,c),prep_adj(x,st)),st);
        AhA       = @(x,st)   Ah(A(x,st),st);  
        A_sense   = @(x,st,SEs,Nt)     NUFFT_SoS(bsxfun(@times,prep(x,st),SEs),st);
        Ah_sense  = @(x,st,SEs,Nt)     sum((NUFFT_SoS_adj(prep_adj(x,st),st).*conj(SEs)) ,4); %sum over coils
        AhA_sense = @(x,st,SEs,GhG,Nt) Ah_sense(A_sense(x,st,SEs,Nt),st,SEs,Nt);
    elseif isPreparedirtnufft
        reconOptions.flagUsefinufft = 0;
        A         = @(x,st)   irt_nufft_SoS(prep(x,st), st);
        Ah        = @(x,st)   irt_nufft_SoS_adj(prep_adj(x,st), st);
        Ainv      = @(x,st,c) Ah(bsxfun(@times,dcf(st,c),prep_adj(x,st)), st);
        AhA       = @(x,st)   Ah(A(x,st),st);  
        A_sense   = @(x,st,SEs,Nt)     irt_nufft_SoS(bsxfun(@times,prep(x,st),SEs),st);
        Ah_sense  = @(x,st,SEs,Nt)     sum(prep(irt_nufft_SoS_adj(prep_adj(x,st),st),st).*conj(SEs), 4); %sum over coils
        AhA_sense = @(x,st,SEs,GhG,Nt) Ah_sense(A_sense(x,st,SEs,Nt),st,SEs,Nt);
    else
        error('NUFFT package not available!'); 
    end
    
    % set up preconditioner
    [ky, kx] = ndgrid(-(st.Nd(1)/2):(st.Nd(1)/2-1),-(st.Nd(2)/2):(st.Nd(2)/2-1));
    k        = sqrt(ky.^2+kx.^2);
    k        = k .* (0.54 - 0.46*cos(2*pi*(k+st.Nd(2)/sqrt(2))/(sqrt(2)*st.Nd(2)))); % multiplied by a hamming window
    window                 = k;
    window(window==0)      = 1/8 * N/Norig;
    window                 = ifftshift(window);
    %window(window>(N/2-1)) = min(window(:));
    window                 = window/mean(window(:));
    
    Mf = @(x) ifft2(fft2(x).*repmat(window, [1 1 Nz size(x,4) size(x,5)]));
    M  = @(x) vec(Mf(prep(x,st)));
    
  case 'Cartesian'
    prep      = @(x,st) reshape(x, st.Nd(1), st.Nd(2), st.Nd(3), []);  
    prep_adj  = @(x,st) reshape(x, st.Nd(1), st.Nd(2), st.Nd(3), []);  
    dispim    = @(x,st) x((floor(st.Nd(1)/2)-floor(Nydisp/2))+(1:Nydisp),(floor(st.Nd(2)/2)-floor(Nxdisp/2))+(1:Nxdisp),floor(st.Nd(3)/2)+1,:); 
    
    fftyz     = @(x,st) fft(fft(x,[],1),[],3)/sqrt(size(x,1))/sqrt(size(x,3));
    ifftyz    = @(x,st) ifft(ifft(x,[],1),[],3)*sqrt(size(x,1))*sqrt(size(x,3));
    A         = @(x,st) fftyz(prep(x,st));
    Ah        = @(x,st) ifftyz(prep_adj(x,st));
    Ainv      = @(x,st) Ah(x,st);
    AhA       = @(x,st) Ah(A(prep(x,st)));
    A_sense   = @(x,st,SEs,Nt)     fftyz(bsxfun(@times,prep(x,st),SEs));
    Ah_sense  = @(x,st,SEs,Nt)     sum(Ah(x,st).*conj(SEs),4);
    AhA_sense = @(x,st,SEs,GhG,Nt) Ah_sense(A_sense(x,st,SEs,Nt),st,SEs,Nt);
    
    M         = @(x) vec(ifftyz(bsxfun(@times,fftyz(prep(x,st)), reshape(st.winv,st.Nd(1),[],st.Nd(3)))));
    %M         = @(x) vec(ifftyz(bsxfun(@times,fftyz(prep(x,st)), reshape(st.winv,Ny,[],Nz))));
%   case 'Spiral'
%     prep     = @(x,st) reshape(x, st.Nd(1), st.Nd(2), Nz, []);  
%     prep_adj = @(x,st) reshape(x, st.M, Nz, []);
%     dcff     = @(k)    (k+(k==0)*min(k(k~=0))/4);
%     dcf      = @(st,c) c*dcff(sqrt(sum(abs(st.om).^2,2)));
%     crop     = @(x,st) x(1:st.Nd(1),1:st.Nd(2),:,:,:);
%     dispim   = @(x)    x(floor(st.Nd(1)/2)-floor(Nydisp/2) + (1:Nydisp), floor(st.Nd(2)/2)-floor(Nxdisp/2) + (1:Nxdisp), 1, :);
%     
%     if flagUseBart && isPreparedbartnufft
%         A         = @(x,st) NUFFT_SoS_bart(x,st);
%         Ah        = @(x,st) NUFFT_SoS_adj_bart(x,st);
%         Ainv      = @(x,st,c) NUFFT_SoS_adj_bart(x,st,1);
%         AhA       = @(x,st)   Ah(A(x,st),st);  
%         A_sense   = @(x,st,SEs,Nt)     NUFFT_SoS_bart(bsxfun(@times,prep(x,st),SEs),st);
%         Ah_sense  = @(x,st,SEs,Nt)     sum((NUFFT_SoS_adj_bart(prep_adj(x,st),st).*conj(SEs)) ,4); %sum over coils
%         AhA_sense = @(x,st,SEs,GhG,Nt) Ah_sense(A_sense(x,st,SEs,Nt),st,SEs,Nt);
%     elseif flagUsefinufft && isPreparedfinufft
%         reconOptions.flagUseGPU = 0;
%         A         = @(x,st)   finufft_SoS(prep(x,st), st);
%         Ah        = @(x,st)   finufft_SoS_adj(prep_adj(x,st), st);
%         Ainv      = @(x,st,c) Ah(bsxfun(@times,dcf(st,c),prep_adj(x,st)), st);
%         AhA       = @(x,st)   Ah(A(x,st),st);  
%         A_sense   = @(x,st,SEs,Nt)     finufft_SoS(bsxfun(@times,prep(x,st),SEs),st);
%         Ah_sense  = @(x,st,SEs,Nt)     sum(prep(finufft_SoS_adj(prep_adj(x,st),st),st).*conj(SEs), 4); %sum over coils
%         AhA_sense = @(x,st,SEs,GhG,Nt) Ah_sense(A_sense(x,st,SEs,Nt),st,SEs,Nt);
%     elseif flagUseGPU && isPreparedgpuNUFFT
%         A         = @(x,st) NUFFT_SoS(x,st);
%         Ah        = @(x,st) NUFFT_SoS_adj(prep_adj(x,st),st);
%         Ainv      = @(x,st,c) Ah(bsxfun(@times,dcf(st,c),prep_adj(x,st)),st);
%         AhA       = @(x,st)   Ah(A(x,st),st);  
%         A_sense   = @(x,st,SEs,Nt)     NUFFT_SoS(bsxfun(@times,prep(x,st),SEs),st);
%         Ah_sense  = @(x,st,SEs,Nt)     sum((NUFFT_SoS_adj(prep_adj(x,st),st).*conj(SEs)) ,4); %sum over coils
%         AhA_sense = @(x,st,SEs,GhG,Nt) Ah_sense(A_sense(x,st,SEs,Nt),st,SEs,Nt);
%     elseif isPreparedirtnufft
%         reconOptions.flagUsefinufft = 0;
%         A         = @(x,st)   irt_nufft_SoS(prep(x,st), st);
%         Ah        = @(x,st)   irt_nufft_SoS_adj(prep_adj(x,st), st);
%         Ainv      = @(x,st,c) Ah(bsxfun(@times,dcf(st,c),prep_adj(x,st)), st);
%         AhA       = @(x,st)   Ah(A(x,st),st);  
%         A_sense   = @(x,st,SEs,Nt)     irt_nufft_SoS(bsxfun(@times,prep(x,st),SEs),st);
%         Ah_sense  = @(x,st,SEs,Nt)     sum(prep(irt_nufft_SoS_adj(prep_adj(x,st),st),st).*conj(SEs), 4); %sum over coils
%         AhA_sense = @(x,st,SEs,GhG,Nt) Ah_sense(A_sense(x,st,SEs,Nt),st,SEs,Nt);
%     else
%         error('NUFFT package not available!'); 
%     end
%     
%     % set up preconditioner
%     [ky, kx] = ndgrid(-(st.Nd(1)/2):(st.Nd(1)/2-1),-(st.Nd(2)/2):(st.Nd(2)/2-1));
%     k        = sqrt(ky.^2+kx.^2);
%     k        = k .* (0.54 - 0.46*cos(2*pi*(k+st.Nd(2)/sqrt(2))/(sqrt(2)*st.Nd(2)))); % multiplied by a hamming window
%     window                 = k;
%     window(window==0)      = 1/8 * N/Norig;
%     window                 = ifftshift(window);
%     %window(window>(N/2-1)) = min(window(:));
%     window                 = window/mean(window(:));
%     
%     Mf = @(x) ifft2(fft2(x).*repmat(window, [1 1 Nz size(x,4) size(x,5)]));
%     M  = @(x) vec(Mf(prep(x,st)));

end
