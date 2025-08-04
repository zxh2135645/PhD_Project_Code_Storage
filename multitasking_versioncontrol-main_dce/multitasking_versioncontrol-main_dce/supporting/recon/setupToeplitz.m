function [reconOptions,dataArray] = setupToeplitz(params,reconOptions,dataArray,temporalBasis)

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(temporalBasis);

vec = @(x) x(:);
prep_adj = @(x,st) reshape(x, st.M, 1, []);

try
    disp('Setting up Toeplitz operator...')
    
    if flagUseGPU
        Ah = @(x,st) NUFFT_SoS_adj(prep_adj(x,st),st);
    elseif flagUsefinufft
        Ah = @(x,st) finufft_SoS_adj(prep_adj(x,st), st);
    else
        Ah = @(x,st) irt_nufft_SoS_adj(prep_adj(x,st), st);
    end
    
    % Calculate scaling factor for Om (varies with NUFFT method)
    delta_im = zeros(N,N);
    delta_im(ceil(N/2+1),ceil(N/2+1)) = 1;
    if flagUseGPU
        %AhAd = vec(st.F'*(st.F*delta_im));
        AhAd = vec(NUFFT_SoS_adj(prep_adj(NUFFT_SoS(delta_im,st),st),st));
    elseif flagUsefinufft
        disp('Warning: CPU NUFFT Toeplitz scaling is untested!')
        AhAd = vec(finufft_SoS_adj(finufft_SoS(delta_im,st),st));
    else
        disp('Warning: CPU NUFFT Toeplitz scaling is untested!')
        AhAd = vec(irt_nufft_SoS_adj(irt_nufft_SoS(delta_im,st),st));
    end
    Ah1 = vec(Ah(ones(st.M,1),st));
    toeplitzScalar = real(pinv(AhAd(:))*Ah1(:));

    [reconOptions,dataArray] = calcToeplitzSize(params,reconOptions,dataArray);
    Om_size = reconOptions.Om_size;
    
    L = size(Phi_rt,1);
    Om = zeros(Ntrajs,1,L,L);
    for traj = 1:Ntrajs
        t_ind = (st.linOrder_shift == st.linOrder_shift(traj));
        Om(traj,1,:,:) = Phi_rt(:,t_ind) * Phi_rt(:,t_ind)';
    end

    if reconOptions.flagUseGPUToeplitz
        try
            Om = reshape(fft2(gpuArray(Ah(repmat(Om,[1 Nkx 1 1]),st)),Om_size(1),Om_size(2)),Om_size(1),Om_size(2),1,L,L)/(toeplitzScalar);
        catch                
            reconOptions.flagUseGPUToeplitz = false;
        end
    end
    if ~reconOptions.flagUseGPUToeplitz
        Om = reshape(fft2(Ah(repmat(Om,[1 Nkx 1 1]),st),Om_size(1),Om_size(2)),Om_size(1),Om_size(2),1,L,L)/(toeplitzScalar);
    end
    
    % correct phase
    Om = bsxfun(@rdivide,Om,sign(Om(:,:,1)));
  
    %ensure positive semi-definite
    fprintf('Ensure positive semi-definite Toeplitz matrices.')
    for j = 1:size(Om,1)
        if mod(j,10)==0
            fprintf('.')
        end
        if flagUseGPU
            Om_j = gather(Om(j,:,:,:,:));
        else
            Om_j = Om(j,:,:,:,:);
        end
        parfor k = 1:size(Om,2)
            [V,D] = eig(squeeze(Om_j(1,k,:,:,:)));
            D = real(diag(D));
            D(D<0) = 0;
            D = diag(D);
            Om_jk = V*sqrt(D);
            Om_jk = Om_jk*Om_jk'; %(V*sqrt(D))*(V*sqrt(D))' is numerically better than V*D*V'
            Om_j(1,k,:,:,:) = Om_jk;
        end
        Om(j,:,:,:,:) = Om_j;
    end
    reconOptions.Om = Om;
    fprintf('.\n')
catch errmsg
    disp('Toeplitz canceled')
    errmsg
    reconOptions.flagUseToeplitz = false;
end
