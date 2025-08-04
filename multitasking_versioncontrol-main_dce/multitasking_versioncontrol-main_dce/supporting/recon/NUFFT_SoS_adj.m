function im = NUFFT_SoS_adj(x, st)

try
    if ~isfield(st, 'gpunufft_adj_scalar')
        st.gpunufft_adj_scalar = 1;
    end
    
    am = getGPUmem(); % get available memory (am)
    if st.flagUseGPU
        if min(am) < 4e3
            gpuWorkerReset(); %reset if < 4e3 MB
        end
    end

    sz   = size(x); sz = [sz 1];
    x    = x(:,:,:);
    Nims = size(x,3);

    st.F = gpuNUFFT(st.om.'/(2*pi), ones(1,st.M), st.osf, 6, 8, st.Nd(1:2), [], true);    
    im = zeros([st.Nd st.Nz Nims]);
    for np = 1:st.Nz
        for nc = 1:Nims  % Ncoils
            im(:,:,np,nc) = st.F'*x(:,np,nc);
        end
    end
  
    if st.Nz > 1
        im = ifft(im,[],3) * sqrt(st.Nz);
    end
    im = reshape(im,[st.Nd,st.Nz,sz(3:end)]) * st.gpunufft_adj_scalar;
  
catch errormsg
    fprintf([errormsg '\n']);
    delete(gcp('nocreate'));
    parpool('local',8)
    im = NUFFT_SoS_adj(x,st);
end

