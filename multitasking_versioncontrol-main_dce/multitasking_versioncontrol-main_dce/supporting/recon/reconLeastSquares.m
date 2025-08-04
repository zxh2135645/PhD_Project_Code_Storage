function [reconOptions,dataArray,temporalBasis,spatialCoeff] = reconLeastSquares(params,reconOptions,dataArray,temporalBasis,spatialCoeff,isTensor)

if nargin < 6
    isTensor = false;
end

if nargin < 5
    spatialCoeff = [];
end

vec = @(x) x(:);
row = @(x) x(:).';

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(temporalBasis);

if MBfactor == 1
    SEs = fftshift(dataArray.SEs,3);
else
    SEs = dataArray.SEs;
end
if isCartesian
    SEs = fftshift(SEs,1);
end

if nargin > 4 && flagUseInitialGuess && isfield(spatialCoeff,'U_init')
    try
        U_guess = vec(reshape(spatialCoeff.U_init,[],L_init)*(bsxfun(@times,Phi_rt_init,dataArray.lrw)*pinv(Phi_rt)));
    catch
        flagUseInitialGuess = false;
    end
else
    flagUseInitialGuess = false;
end

kspaceData = dataArray.kspaceData;

Nkx    = size(kspaceData,2);
Ncoils = size(kspaceData,4);
Necho  = size(kspaceData,5);

if ~isTensor
    kspaceData = kspaceData(:,:,:,:,1);
end
kspaceData = permute(kspaceData,[5 1 2 3 4]);
tempNecho  = size(kspaceData,1);
if isfield(dataArray,'lrw') && isTensor
    if numel(dataArray.lrw) == size(kspaceData,2)
        kspaceData = kspaceData.*reshape(dataArray.lrw,1,[]);
    elseif numel(dataArray.lrw) == size(kspaceData,2)*tempNecho
        kspaceData = kspaceData.*reshape(dataArray.lrw,tempNecho,[]);
    end
end

L = size(Phi_rt,1);
if size(Phi_rt,2) == size(kspaceData,2)
    Phi_rt = permute(Phi_rt,[1 3 2]);
elseif size(Phi_rt,2) == size(kspaceData,2)*size(kspaceData,1)
    Phi_rt = reshape(Phi_rt,L,tempNecho,[]);
end

linOrder = st.linOrder_shift;
parOrder = st.parOrder_shift;

setupFunctions;

if ~flagCommandLine
    % Progress bar
    progress = waitbar(0,'','Name', sprintf('Progress'),...
                       'CreateCancelBtn',...
                       'setappdata(gcbf,''canceling'',1)');
    setappdata(progress,'canceling',0)
end

tic; 

temppath = pwd;

try
    
switch Trajectory
  case {'Radial','Spiral'}
      
    st.flagUseGPU = flagUseGPU;
%     if flagUseGPU
%         delete(gcp('nocreate'));
%         parpool('local',8);
%         gpuWorkerReset();
%     end
    
    if reconOptions.flagUseToeplitz
        fprintf('Setup Toeplitz...')
        if ~flagCommandLine
            waitbar(0.1, progress, sprintf('Setup Toeplitz...'));
        end
        [reconOptions,dataArray] = setupToeplitz(params,reconOptions,dataArray,temporalBasis);
        
        if MBfactor == 1
            SEs = fftshift(dataArray.SEs,3);
        else
            SEs = dataArray.SEs;
        end
    end
    
    fprintf('Aligning k-space... ');
    if ~flagCommandLine
        waitbar(0.3, progress, sprintf('Aligning k-space... '));
    end
    
    FSU = zeros(Ntrajs,Nkx,Nz,Ncoils,L, 'single');
   
    for traj = 1:Ntrajs 
        for np = 1:Nz
            t_ind = (linOrder==traj) & (parOrder==np);
            if sum(t_ind) > 0
                tempPhi = reshape(Phi_rt(:,:,t_ind),L,sum(t_ind)*tempNecho);
                tempk   = reshape(kspaceData(:,t_ind,:,:,:),sum(t_ind)*tempNecho,[]);
                FSU(traj,:,np,:,:) = reshape(tempk.'*tempPhi',Nkx,Ncoils,L);
                %FSU(traj,:,np,:,:) = reshape(kspaceData(t_ind,:).'*Phi_rt(:,t_ind)',Nkx,Ncoils,L);
            end
        end
    end
    
    fprintf('Calculating Ahb and filtAhb... ')
    if ~flagCommandLine
        waitbar(0.5, progress, sprintf('Calculating Ahb and filtAhb...'));
    end
    
    if flagUseBart
        % calculte Ahb/filtAhb using BART
        try
            Ahb = sum(bsxfun(@times, reshape(NUFFT_SoS_adj_bart(FSU, st),[size(SEs) L]), conj(SEs)),4);
            if ~flagUseInitialGuess
                filtAhb = sum(bsxfun(@times, reshape(NUFFT_SoS_adj_bart(FSU,st,1),[size(SEs) L]), conj(SEs)),4);
            end
        catch
            Ahb     = 0;
            filtAhb = 0;
            fprintf('%d coils: coil ',Ncoils);
            for coil = 1:Ncoils
                fprintf('#%d... ',coil);
                Ahb = Ahb + bsxfun(@times, NUFFT_SoS_adj_bart(FSU(:,:,:,coil,:), st), conj(SEs(:,:,:,coil)));
                if ~flagUseInitialGuess
                    filtAhb = filtAhb + bsxfun(@times, NUFFT_SoS_adj_bart(FSU(:,:,:,coil,:),st,1), conj(SEs(:,:,:,coil)));
                end
            end
        end
        Ahb = Ahb(:);
        Ahb(~isfinite(Ahb)) = 0;
        
        if ~flagUseInitialGuess
            filtAhb = vec(bsxfun(@rdivide,filtAhb, sum(abs(SEs).^2,4)));
            filtAhb(~isfinite(filtAhb)) = 0;
        end
    elseif flagUsefinufft
        FSU = reshape(FSU, st.M, Nz*Ncoils*L);
        if flagUseGPU
            cd([reconOptions.mainpath '/supporting/nufft/finufft/build']);
            FU = cufinufft_SoS_adj(FSU, st);
        else
            FU = finufft_SoS_adj(FSU, st);
        end
        FU  = reshape(FU, Ny, Nx, Nz, Ncoils, L);
        Ahb = sum(bsxfun(@times, FU, conj(SEs)),4);
        Ahb(~isfinite(Ahb)) = 0;

        % calculate filtAhb
        if ~flagUseInitialGuess
            if flagUseGPU
                FU = cufinufft_SoS_adj(bsxfun(@times,dcf(st,st.c),FSU),st);
            else
                FU = finufft_SoS_adj(bsxfun(@times,dcf(st,st.c),FSU),st);
            end
            FU = reshape(FU, Ny, Nx, Nz, Ncoils, L);
            filtAhb = sum(bsxfun(@times, FU, conj(SEs)),4);
            filtAhb = vec(bsxfun(@rdivide,filtAhb, sum(abs(SEs).^2,4)));
            filtAhb(~isfinite(filtAhb)) = 0;
        end 
    elseif flagUseGPU
        % calculte Ahb/filtAhb on GPU using gpuNUFFT
        FSU = reshape(FSU, st.M, Nz, Ncoils,[]);
  
        Ahb     = 0;
        filtAhb = 0;
        for coil = 1:Ncoils
            Ahb = Ahb + bsxfun(@times, NUFFT_SoS_adj(FSU(:,:,coil,:), st), conj(SEs(:,:,:,coil)));
            if ~flagUseInitialGuess
                filtAhb = filtAhb + bsxfun(@times, NUFFT_SoS_adj(bsxfun(@times, dcf(st,st.c), FSU(:,:,coil,:)),st), conj(SEs(:,:,:,coil)));
            end
        end
        Ahb = Ahb(:);
        Ahb(~isfinite(Ahb)) = 0;
        if ~flagUseInitialGuess
            filtAhb = vec(bsxfun(@rdivide,filtAhb, sum(abs(SEs).^2,4)));
            filtAhb(~isfinite(filtAhb)) = 0;
        end
    else        % calculte Ahb/filtAhb usng irt
        FSU = reshape(FSU, st.M, Nz*Ncoils*L);

        % calculate Ahb 
        FU = irt_nufft_SoS_adj(FSU, st);
        FU  = reshape(FU, Ny, Nx, Nz, Ncoils, L);
        Ahb = sum(bsxfun(@times, FU, conj(SEs)),4);
        Ahb(~isfinite(Ahb)) = 0;

        % calculate filtAhb
        if ~flagUseInitialGuess
            FU = irt_nufft_SoS_adj(bsxfun(@times, dcf(st,st.c), FSU),st);
            FU = reshape(FU, Ny, Nx, Nz, Ncoils, L);
            filtAhb = sum(bsxfun(@times, FU, conj(SEs)),4);
            filtAhb = vec(bsxfun(@rdivide,filtAhb, sum(abs(SEs).^2,4)));
            filtAhb(~isfinite(filtAhb)) = 0;
        end 
    end
    clear FSU FU
    
    %% pre-calculate Phi^2
    fprintf('Calculating Phi^2 ... \n')
    if ~flagCommandLine
        waitbar(0.7, progress, sprintf('Calculating Phi^2 ...'));
    end
    
    Phi2 = zeros(L, L, Ntrajs, Nz);  
    for traj = 1:Ntrajs
        for np = 1:Nz
            t_ind = (linOrder==traj) & (parOrder==np);
            if sum(t_ind) > 0
                tempPhi = reshape(Phi_rt(:,:,t_ind),L,[]);
                Phi2(:,:,traj,np) = tempPhi * tempPhi';
            end
        end
    end
    
    if reconOptions.flagUseToeplitz
        fprintf('Toeplitz recon...')
        if ~flagCommandLine
            waitbar(0.8, progress, sprintf('Toeplitz recon...'));
        end
        [reconOptions,dataArray] = setupToeplitz(params,reconOptions,dataArray,temporalBasis);
    end
    
    if reconOptions.flagUseToeplitz
        if ~flagCommandLine
            waitbar(0.8, progress, sprintf('Toeplitz recon...'));
        end
        if flagUseInitialGuess
            U0 = U_guess;
        else 
            U0 = AhA_ps_Toeplitz(filtAhb, st, SEs, Phi2, reconOptions.Om, flagUseGPUToeplitz);
            c2 = real(pinv(U0(:),eps)*Ahb(:));
            U0 = filtAhb*c2;
        end
        
        if ~flagCommandLine
            delete(progress);
        end
        
        fprintf('Preconditioned conjugate gradient.')
        U = pcg1(@(x) AhA_ps_Toeplitz(x, st, SEs, Phi2, reconOptions.Om, flagUseGPU), Ahb(:), [], 10, M, [], U0(:));
    elseif flagUseBart
        % ls-recon using BART
        fprintf('BART recon...')
        if ~flagCommandLine
            waitbar(0.8, progress, sprintf('BART recon...'));
        end
        if flagUseInitialGuess
            U0 = U_guess;
        else 
            if ~isfield(reconOptions,'ls_lowmem')
                reconOptions.ls_lowmem = 0;
            end

            [U0,reconOptions.ls_lowmem] = AhA_ps3D_bart(filtAhb,st,SEs,Phi2,reconOptions.ls_lowmem);
            c2 = real(pinv(U0(:),eps)*Ahb(:));
            U0 = filtAhb*c2;
        end

        if ~flagCommandLine
            delete(progress);
        end

        fprintf('Preconditioned conjugate gradient.')

        [U,~] = pcg1(@(x) AhA_ps3D_bart(x, st, SEs, Phi2,reconOptions.ls_lowmem), Ahb(:), [], 10, M, [], U0(:));
    elseif flagUsefinufft
        % ls-recon using finufft/cufinufft
        if flagUseInitialGuess
            U0 = U_guess;
        else  
            if flagUseGPU
                if ~flagCommandLine
                    waitbar(0.8, progress, sprintf('cufinufft recon...'));
                end
                fprintf('cufinufft recon...')
%                 U0 = AhA_ps_cufinufft(single(filtAhb), st, single(SEs), single(Phi2));
                U0 = AhA_ps3D_cufinufft(filtAhb, st, SEs, Phi2);
                reconOptions.st.nufftlowmem = 0;
            else
                if ~flagCommandLine
                    waitbar(0.8, progress, sprintf('finufft recon...'));
                end
                fprintf('finufft recon...')
                [U0,st.nufftlowmem] = AhA_ps3D_finufft(filtAhb, st, SEs, Phi2);
                reconOptions.st.nufftlowmem = st.nufftlowmem;
            end
            c2 = real(pinv(U0(:),eps)*Ahb(:));
            U0 = filtAhb*c2;
        end

        if ~flagCommandLine
            delete(progress);
        end

        fprintf('Preconditioned conjugate gradient.')

        if flagUseGPU
%             [U,~] = pcg1(@(x) AhA_ps_cufinufft(single(x), st, single(SEs), single(Phi2)), single(Ahb(:)), [], 10, M, [], single(U0(:)));
            [U,~] = pcg1(@(x) AhA_ps3D_cufinufft(x, st, SEs, Phi2), single(Ahb(:)), [], 10, M, [], U0(:));
        else
            [U,~] = pcg1(@(x) AhA_ps3D_finufft(x, st, SEs, Phi2), single(Ahb(:)), [], 10, M, [], U0(:));
        end
    elseif flagUseGPU
        % ls-recon using gpuNUFFT
        fprintf('gpuNUFFT recon...')
        if ~flagCommandLine
            waitbar(0.8, progress, sprintf('gpuNUFFT recon...'));
        end
        if flagUseInitialGuess
            U0 = U_guess;
        else 
            U0 = AhA_ps3D(filtAhb, st, SEs, Phi2);
            c2 = real(pinv(U0(:),eps)*Ahb(:));
            U0 = filtAhb*c2;
        end

        if ~flagCommandLine
            delete(progress);
        end

        fprintf('Preconditioned conjugate gradient.')

        [U,~] = pcg1(@(x) AhA_ps3D(x, st, SEs, Phi2), Ahb(:), [], 10, M, [], U0(:));
    else
        % ls-recon using irt
        if flagUseInitialGuess
            U0 = U_guess;
        else  
            if ~flagCommandLine
                waitbar(0.8, progress, sprintf('irt_nufft recon...'));
            end
            fprintf('irt_nufft recon...')
            U0 = AhA_ps3D_irt(filtAhb, st, SEs, Phi2);
            c2 = real(pinv(U0(:),eps)*Ahb(:));
            U0 = filtAhb*c2;
        end

        if ~flagCommandLine
            delete(progress);
        end

        fprintf('Preconditioned conjugate gradient.')
        [U,~] = pcg1(@(x) AhA_ps3D_irt(x, st, SEs, Phi2), Ahb(:), [], 10, M, [], U0(:));
    end
    
  case 'Cartesian'
    
    fprintf('Aligning k-space... ');
    if ~flagCommandLine
        waitbar(0.2, progress, sprintf('Aligning k-space... '));
    end 

    FSU = zeros(Ny, Nx, Nz, Ncoils, L, 'single');

    for npy = 1:Ny
        for npz = 1:Nz
            t_ind = (linOrder==npy) & (parOrder==npz);
            if sum(t_ind) > 0
                tempPhi = reshape(Phi_rt(:,:,t_ind),L,sum(t_ind)*tempNecho);
                tempk   = reshape(kspaceData(:,t_ind,:,:,:),sum(t_ind)*tempNecho,[]);
                FSU(npy,:,npz,:,:) = reshape(tempk.'*tempPhi',Nx,Ncoils,L);
                %FSU(npy,Nx-size(kspaceData,2)+1:Nx,npz,:,:) = reshape(kspaceData(t_ind,:).'*Phi_rt(:,t_ind)',Nx,Ncoils,L);
            end
        end
    end    
    clear kspaceData;
    
    fprintf('Calculating Ahb and filtAhb... ');
    if ~flagCommandLine
        waitbar(0.5, progress, sprintf('Calculating Ahb and filtAhb...'));
    end
    Ahb = vec(sum(bsxfun(@times,conj(SEs),reshape(Ah(FSU,st),Ny,Nx,Nz,Ncoils,L)),4));
    if ~flagUseInitialGuess
        filtAhb = vec(bsxfun(@rdivide,sum(bsxfun(@times,conj(SEs),reshape(Ah(bsxfun(@times,st.winv,FSU),st),Ny,Nx,Nz,Ncoils,L)),4),sqrt(sum(abs(SEs).^2,4))));
    end
    clear FSU; 
    
    %% pre-calculate Phi^2
    fprintf('Calculating Phi^2 ... ')
    if ~flagCommandLine
        waitbar(0.7, progress, sprintf('Calculating Phi^2 ...'));
    end
    Phi2 = zeros(L, L, Ny, Nz);  
    for npy = 1:Ny
        for npz = 1:Nz
            t_ind = (linOrder==npy) & (parOrder==npz);
            if sum(t_ind) > 0
                tempPhi = reshape(Phi_rt(:,:,t_ind),L,[]);
                Phi2(:,:,npy,npz) = tempPhi * tempPhi';
            end
        end
    end
     
    if flagUseInitialGuess
        U0 = U_guess;
    else 
        if ~isfield(reconOptions,'ls_lowmem')
            reconOptions.ls_lowmem = 0;
        end

        [U0,reconOptions.ls_lowmem] = AhA_ps_cart(filtAhb,SEs,Phi2,reconOptions.ls_lowmem);

        c2 = real(pinv(U0(:),eps)*Ahb(:));
        U0 = filtAhb*c2;
    end
 
    if ~flagCommandLine
        delete(progress);
    end

    fprintf('Preconditioned conjugate gradient.'); 
    [U,~] = pcg1(@(x)AhA_ps_cart(x,SEs,Phi2,reconOptions.ls_lowmem),Ahb(:),[],10,M,[],U0(:));
end

catch errormsg
    fprintf(2, '%s\n', errormsg.message);
    fprintf(2, 'Least-squares recon failed.\n');
end
toc;

cd(temppath);

temporalBasis.Phi2   = Phi2;
spatialCoeff.Ahb     = Ahb;
spatialCoeff.filtAhb = reshape(U0,Ny,Nx,Nz,[]);
spatialCoeff.U       = U;

if isTensor
    spatialCoeff.U_tensor = U;
else
    spatialCoeff.U_init = U;
end
