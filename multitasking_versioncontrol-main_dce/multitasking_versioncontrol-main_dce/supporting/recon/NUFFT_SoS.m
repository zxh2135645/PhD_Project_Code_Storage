function out = NUFFT_SoS(x, st)

try
    if ~isfield(st, 'gpunufft_fwd_scalar')
        st.gpunufft_fwd_scalar = 1;
    end
   
    am = getGPUmem(); %get available memory (am)
    if st.flagUseGPU
       if min(am) < 4e3
            gpuWorkerReset(); %reset if < 4e3 MB
        end
    end

    sz = size(x); sz = [sz 1 1];
    x  = x(:,:,:,:);
    Nims = size(x,4);
        
    if st.Nz > 1   
        x = fft(x,[],3) / sqrt(st.Nz);
    end
    
    out  = zeros(st.M, st.Nz, Nims);
    st.F = gpuNUFFT(st.om.'/(2*pi), ones(1,st.M), st.osf, 6, 8, st.Nd(1:2), [], true);    
    for np = 1:st.Nz
        for nc = 1:Nims
            out(:,np,nc) = st.F*x(:,:,np,nc);
        end
    end

    out = reshape(out,[st.M, st.Nz, sz(4:end)]) * st.gpunufft_fwd_scalar;
    
catch errormsg
    disp(errormsg);
    delete(gcp('nocreate'));
    parpool('local',8)
    out = NUFFT_SoS(x, st);
end

