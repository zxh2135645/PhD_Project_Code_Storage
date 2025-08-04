function AhAx = AhA_ps_Toeplitz(x,st,SEs,Phi2,Om,flagUseGPUToeplitz)

if nargin < 6
    flagUseGPUToeplitz = true;
end

L  = size(Phi2,1);
N  = st.Nd(1);
Nz = st.Nz;
Om_size = size(Om);

cropToeplitz = @(x)x(1:N,1:N,:,:,:);
vec = @(x) x(:);

x = reshape(x, N, N, Nz, 1, L);
x = bsxfun(@times, x, SEs);

try
    if flagUseGPUToeplitz   % try Toeplitz on GPU
        try
            for z = 1:Nz
                AhAx_Temp = fft2(gpuArray(x(:,:,z,:,:)),Om_size(1),Om_size(2));
                temp = zeros(size(AhAx_Temp),'like',AhAx_Temp);
                for j = 1:L
                    temp(:,:,:,j) = sum(bsxfun(@times,squeeze(AhAx_Temp),Om(:,:,:,:,j)),4);
                end
                AhAx(:,:,z,:,:) = sum(bsxfun(@times,cropToeplitz(ifft2(temp)),conj(squeeze(SEs(:,:,z,:)))),3);
            end
            AhAx = vec(AhAx);
        catch
            disp('Toeplitz GPU failed, switching to CPU')
            flagUseGPUToeplitz = false;
            Om = gather(Om);
        end
    end
    
    if ~flagUseGPUToeplitz      % try Toeplitz on CPU
        AhAx = fft2(x,Om_size(1),Om_size(2));
        temp = zeros(size(AhAx),'like',AhAx);
        Om = reshape(Om,Om_size(1),Om_size(2),1,1,L,L);
        for j = 1:L
            temp(:,:,:,:,j) = sum(bsxfun(@times,AhAx,Om(:,:,:,:,:,j)),5);
        end
        AhAx = vec(sum(bsxfun(@times,cropToeplitz(ifft2(temp)),conj(SEs)),4));
    end
catch       %if it fails, try nufft
    disp('Toeplitz failed')
    if st.flagUseGPU
        AhAx = AhA_ps3D(x, st, SEs, Phi2);
    elseif flagUsefinufft
        AhAx = AhA_ps3D_finufft(gather(x), st, SEs, Phi2);
    else
        AhAx = AhA_ps3D_irt_nufft(gather(x), st, SEs, Phi2);
    end
end


