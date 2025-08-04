function [params, reconOptions, dataArray] = setupTrajectories(params, reconOptions, dataArray)

flagUsefinufft = 1;
flagUseBart    = 0;

Ntrajs  = 0;
trajInc = 0;

cropFOVx = 1;
cropFOVy = 1;

vec = @(x) x(:);

% load parameters
extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(dataArray);

switch Trajectory
  %% Radial, stack of stars
  case 'Radial'
    if flagCommandLine
        disp('Trajectory: Radial');
    end
    if ~flagIsTrajectorySet

        parOrder = parOrder(:);
        parOrder = parOrder(1:size(kspaceData,1));
        
        % readout trajectory angle
        if sum(diff(linOrder)) == 0
            if Ntrajs == 0 || trajInc == 0
                if isTinyGoldenAngle
                    % Tiny Gloden Angle, by Vincent X. Mao 1/06/20
                    fib = [1 SGBlock-1];
                    tempNtrajs = fib(2) + fib(1)*( SGBlock - 1.0 - 1.0);
                    if Ntrajs == 0
                        bound = params.lBaseResolution * pi / 2;
                    else
                        bound = Ntrajs;
                    end
                    while (tempNtrajs < bound)
                        fib(2)  = fib(1) + fib(2);
                        fib(1)  = fib(2) - fib(1);
                        tempNtrajs  = fib(2) + fib(1)*(SGBlock - 1.0 - 1.0);
                        trajInc = fib(1);
                    end  
                    Ntrajs = tempNtrajs;
                    params.Ntrajs  = Ntrajs;
                    params.trajInc = trajInc;        
                else
                    % Gloden Angle
                    fib = [1 1];

                    while (fib(2) < params.lBaseResolution * pi / 2)
                        fib(2) = fib(1) + fib(2);
                        fib(1) = fib(2) - fib(1);
                    end
                    Ntrajs  = fib(2);
                    trajInc = fib(1);
                    params.Ntrajs  = Ntrajs;
                    params.trajInc = trajInc;  
                end
            end

            % radial angle increment in degree
            theta = 180 * trajInc / Ntrajs;
            % raw linOrder
            linOrder = (mod(trajInc * ((1:Nread) - 1 + (cutoff_shot*linesPerShot - ceil(cutoff_shot*linesPerShot/SGBlock))), 2*Ntrajs) + 1);
        end
        
        linOrder = linOrder(:);
        thetas   = 180 * (linOrder-1) / Ntrajs;

        %% Gradient delay correction
        if flagGradientDelayCorr > 0
            fprintf('Gradient delay correction... ');
            [kspaceData,navData] = gradientDelayCorr(kspaceData,navData,Ntrajs,Norig,linOrder,parOrder,flagGradientDelayCorr,isVE, DC_kz);
            fprintf('done.\n');
        end

        %% Eddy current correction
        if flagEddyCurrentCorr
            fprintf('Eddy current compensation... ');
            [kspaceData, navData] = eddyCompensation(params,kspaceData,navData,thetas,parOrder);
            fprintf('done.\n');
        end

        %%
        % force all angles into [0, 180)-degree range
        kspaceData(linOrder>Ntrajs,:,:,:,:) = kspaceData(linOrder>Ntrajs,[1 end:-1:2],:,:,:);

        % radial trajectory
        dataArray.parOrder   = parOrder;
        dataArray.linOrder   = mod(linOrder-1, Ntrajs) + 1;
        params.thetas        = 180 * (dataArray.linOrder-1) / Ntrajs;
        
        if ~exist('cropFOVy','var')
            cropFOVy = 1;
        end
        if ~exist('cropFOVx','var')
            cropFOVx = 1;
        end
        if cropFOVy > 1
            Ny = floor(Ny/cropFOVy/2)*2;
            if Nydisp>Ny
                Nydisp = Ny;
            end
        end
        if cropFOVx > 1
            Nx = floor(Ny/cropFOVx/2)*2;
            if Nxdisp>Nx
                Nxdisp = Nx;
            end
        end
        N     = max(Ny,Nx);
        Norig = min(min(Ny,Nx),Norig);
        
        params.N = N;
        params.Norig = Norig;
        params.Nx = Nx;
        params.Ny = Ny;
        params.Nxdisp = Nxdisp;
        params.Nydisp = Nydisp;        
    end
 
    % setup NUFFT trajectory
    r   = linspace(-pi, pi, Nx+1); r(end)=[]; r(1:Nx-Nkx) = [];
    om1 = sin(vec(0:Ntrajs-1)*pi/Ntrajs)*r;
    om2 = cos(vec(0:Ntrajs-1)*pi/Ntrajs)*r;
    om  = single([om1(:), om2(:)]);

    traj_bart = zeros([3 size(om1)]);
    traj_bart(1,:,:) = om1/cropFOVy;
    traj_bart(2,:,:) = om2/cropFOVx;
    traj_bart(3,:,:) = 0;
    traj_bart = traj_bart*Nx/2/pi;

    % setup NUFFT object
    fprintf('Configuring and testing NUFFT packages... \n');

    if ~exist('nufftlowmem','var') %|| nufftlowmem == 0
        nufftlowmem = 0;
    end
    if nufftlowmem
        osf = 1.25;
        Kd  = ceil(osf*[N N]/2)*2;
    else
        osf = 2.0;
        Kd  = ceil(osf*[N N]/2)*2;
    end

    try
        % if MIRT
        st = nufft_init(om,[Ny Nx],[6 6],Kd,floor([N N]/2));
    catch
        st.om = om;
        st.M  = size(om,1);
        st.Nd = [Ny Nx];
    end
    st.traj_bart = traj_bart;
    st.Nz     = Nz;
    st.Ndisp  = [Nydisp Nxdisp Nz];
    st.osf    = osf;
    st.nufftlowmem = nufftlowmem;
    st.flagUseGPU  = flagUseGPU;
%         try 
%             % if irt_nufft
%             st.GhG = st.p.G'*st.p.G;
%         catch
%             st.GhG = [];
%         end

    % ifftshift parOrder
    st.parOrder_shift = mod(parOrder-1 - floor(Nz/2), Nz) + 1;
    st.linOrder_shift = dataArray.linOrder;

%% Cartesian
case 'Cartesian'
    disp('Trajectory: Cartesian')
    st.Nd = [Ny Nx Nz];
    st.Ndisp = [Nydisp Nxdisp Nz];
    
    % ifftshift linOrder and parOrder
    st.linOrder_shift = mod(linOrder - DC_ky, Ny) + 1;
    st.parOrder_shift = mod(parOrder - DC_kz, Nz) + 1;
    
    % k-space weighting functions with ifftshift
    st.w = zeros(Ny,1,Nz);
    for npy = 1:Ny
        for npz = 1:Nz
            t_ind = (st.linOrder_shift==npy) & (st.parOrder_shift==npz);
            st.w(npy,:,npz) = sum(t_ind(:));
        end
    end
    %st.w(1,:,1) = st.w(1,:,1) - size(navData,1);
    
    st.winv = 1./st.w;
    st.winv(st.w==0) = 0;
    
    st.flagUseGPU = false;
    
    reconOptions.st = st;
    reconOptions.flagUseGPU = false;

%% Spiral
case 'Spiral'
    disp('Trajectory: Spiral')
    % setup NUFFT trajectory
    kTraj  = loadSpiralTrajectory(params,1);
    Ntrajs = size(kTraj,1);
    Nkx    = size(kTraj,2);
    
    if Ntrajs ~= max(linOrder)
        error('Trajectory does not match kspaceData!');
    end

    om = single(reshape(kTraj,[],2)*dReadoutFOV_mm*1e-3*2*pi/Nx);

    traj_bart = permute(kTraj,[3 1 2])*dReadoutFOV_mm*1e-3;
    traj_bart(1,:,:) = traj_bart(1,:,:)/cropFOVy;
    traj_bart(2,:,:) = traj_bart(2,:,:)/cropFOVx;
    traj_bart(3,:,:) = 0;

    % setup NUFFT object
    fprintf('Configuring and testing NUFFT packages... \n');

    if ~exist('nufftlowmem','var') %|| nufftlowmem == 0
        nufftlowmem = 0;
    end
    if nufftlowmem
        osf = 1.25;
        Kd  = ceil(osf*[N N]/2)*2;
    else
        osf = 2.0;
        Kd  = ceil(osf*[N N]/2)*2;
    end

    try
        % if MIRT
        st = nufft_init(om,[Ny Nx],[6 6],Kd,floor([N N]/2));
    catch
        st.om = om;
        st.M  = size(om,1);
        st.Nd = [Ny Nx];
    end
    st.traj_bart = traj_bart;
    st.Nz     = Nz;
    st.Ndisp  = [Nydisp Nxdisp Nz];
    st.osf    = osf;
    st.nufftlowmem = nufftlowmem;
    st.flagUseGPU  = flagUseGPU;

    if size(kspaceData,2) > size(kTraj,2)
        if ~isMultitasking
            navData = kspaceData(:,1:(size(kspaceData,2)-Nkx),:,:);
        else
            navData(:,1:(size(kspaceData,2)-Nkx),:,:) = [];
        end
        kspaceData(:,1:(size(kspaceData,2)-Nkx),:,:) = [];
    end
        
    % ifftshift parOrder
    st.parOrder_shift = mod(parOrder-1 - floor(Nz/2), Nz) + 1;
    st.linOrder_shift = dataArray.linOrder;

    params.Nkx    = Nkx;
    params.Ntrajs = Ntrajs;
end

% If non-Cartesian trajectory
% Test nufft packages and estimate their forward & adjoint scalars
if ~isCartesian
    delta = zeros([st.Nd]);
    delta(floor(st.Nd(1)/2)+1,floor(st.Nd(2)/2)+1) = 1;

    % BART
    try
        st.Nz = 1;
        st.bartnufft_fwd_scalar = 1;
        st.bartnufft_adj_scalar = 1;
        Adelta = NUFFT_SoS_bart(delta,st);
        st.bartnufft_fwd_scalar = sum(abs(Adelta(:)))/norm(Adelta(:))^2/sqrt(prod(st.Nd));
        Ah1 = NUFFT_SoS_adj_bart(ones(size(st.traj_bart,2),size(st.traj_bart,3)),st);
        st.bartnufft_adj_scalar =  st.M/sqrt(prod(st.Nd))/max(abs(Ah1(:)));
        st.Nz = Nz;
        if isfinite(st.bartnufft_adj_scalar) && isfinite(st.bartnufft_fwd_scalar)
            fprintf('-- BART      passed. bartnufft_fwd_scalar = %.4f, bartnufft_adj_scalar = %.4f\n', st.bartnufft_fwd_scalar, st.bartnufft_adj_scalar);
            isPreparedbartnufft = 1;
        else
            fprintf('-- BART      failed.\n');
            isPreparedbartnufft = 0;
        end
    catch
        fprintf('-- BART      failed.\n');
        isPreparedbartnufft = 0;
    end

    if ~isPreparedbartnufft
        reconOptions.flagUseBart = 0;
    else
        reconOptions.flagUseBart = flagUseBart;
    end

    % gpuNUFFT
    try
        st.F = gpuNUFFT(st.om.'/(2*pi), ones(1,st.M), st.osf, 6, 8, st.Nd(1:2), [], true);  
        st.Nz = 1;
        st.gpunufft_fwd_scalar = 1;
        st.gpunufft_adj_scalar = 1;
        Adelta = NUFFT_SoS(delta,st);
        st.gpunufft_fwd_scalar = sum(abs(Adelta(:)))/norm(Adelta(:))^2/sqrt(prod(st.Nd));
        Ah1 = NUFFT_SoS_adj(ones(st.M,1),st);
        st.gpunufft_adj_scalar =  st.M/sqrt(prod(st.Nd))/max(abs(Ah1(:)));
        st.Nz = Nz;
        if isfinite(st.gpunufft_fwd_scalar) && isfinite(st.gpunufft_adj_scalar)
            fprintf('-- gpuNUFFT  passed. gpunufft_fwd_scalar  = %.4f, gpunufft_adj_scalar  = %.4f\n', st.gpunufft_fwd_scalar, st.gpunufft_adj_scalar);
            isPreparedgpuNUFFT = 1;
        else
            fprintf('-- gpuNUFFT  failed. Switching off gpuNUFFT.\n');
            isPreparedgpuNUFFT = 0;
        end
    catch
        fprintf('-- gpuNUFFT  failed. Switching off gpuNUFFT.\n');
        isPreparedgpuNUFFT = 0;
    end

    % finufft
    try
        st.Nz = 1;
        st.finufft_fwd_scalar = 1;
        st.finufft_adj_scalar = 1;
        Adelta = finufft_SoS(delta,st);
        st.finufft_fwd_scalar = sum(abs(Adelta(:)))/norm(Adelta(:))^2/sqrt(prod(st.Nd));
        Ah1 = finufft_SoS_adj(ones(st.M,1),st);
        st.finufft_adj_scalar = st.M/sqrt(prod(st.Nd))/max(abs(Ah1(:))); 
        st.Nz = Nz;
        if isfinite(st.finufft_fwd_scalar) && isfinite(st.finufft_adj_scalar)
            fprintf('-- finufft   passed. finufft_fwd_scalar   = %.4f, finufft_adj_scalar   = %.4f\n', st.finufft_fwd_scalar, st.finufft_adj_scalar);
            isPreparedfinufft = 1;
        else
            fprintf('-- finufft   failed. Switching off finufft.\n');
            reconOptions.flagUsefinufft = 0;
            isPreparedfinufft = 0;
        end
    catch
        fprintf('-- finufft   failed. Switching off finufft.\n');
        reconOptions.flagUsefinufft = 0;
        isPreparedfinufft = 0;
    end

    if ~isPreparedfinufft
        reconOptions.flagUsefinufft = 0;
    else
        reconOptions.flagUsefinufft = flagUsefinufft;
    end
    
    if ~isPreparedfinufft
        reconOptions.flagUsefinufft = 0;
    else
        reconOptions.flagUsefinufft = flagUsefinufft;
    end

    % cufinufft
    temppath = pwd;
    try
        cd([reconOptions.mainpath '/supporting/nufft/finufft/build']);
        st.Nz = 1;
        st.cufinufft_fwd_scalar = st.finufft_fwd_scalar;
        st.cufinufft_adj_scalar = 1;
%         Adelta = finufft_SoS(delta,st);
%         st.cufinufft_fwd_scalar = sum(abs(Adelta(:)))/norm(Adelta(:))^2/sqrt(prod(st.Nd));
        Ah1 = cufinufft_SoS_adj(ones(st.M,1)*1j,st);
        st.cufinufft_adj_scalar = st.M/sqrt(prod(st.Nd))/max(abs(Ah1(:))); 
        st.Nz = Nz;
        st.cufinufft_AhA_scalar = st.finufft_fwd_scalar * st.cufinufft_adj_scalar / st.Nz;
        if isfinite(st.finufft_fwd_scalar) && isfinite(st.cufinufft_adj_scalar)
            fprintf('-- cufinufft passed. cufinufft_fwd_scalar = %.4f, cufinufft_adj_scalar = %.4f\n', st.cufinufft_fwd_scalar, st.cufinufft_adj_scalar);
            isPreparedcufinufft = 1;
        else
            fprintf('-- cufinufft failed. \n');
            isPreparedcufinufft = 0;
        end
    catch
        fprintf('-- cufinufft failed. \n');
%         reconOptions.flagUsefinufft = 0;
        isPreparedcufinufft = 0;
    end
    cd(temppath);

    % MIRT
    try
        st.Nz = 1;
        st.irtnufft_fwd_scalar = 1;
        st.irtnufft_adj_scalar = 1;
        Adelta = irt_nufft_SoS(delta,st);
        st.irtnufft_fwd_scalar = sum(abs(Adelta(:)))/norm(Adelta(:))^2/sqrt(prod(st.Nd));
        Ah1 = irt_nufft_SoS_adj(ones(st.M,1),st);
        st.irtnufft_adj_scalar =  st.M/sqrt(prod(st.Nd))/max(abs(Ah1(:)));
        st.Nz = Nz;
        if isfinite(st.irtnufft_fwd_scalar) && isfinite(st.irtnufft_adj_scalar)
            fprintf('-- MIRT      passed. irtnufft_fwd_scalar  = %.4f, irtnufft_adj_scalar  = %.4f\n', st.irtnufft_fwd_scalar, st.irtnufft_adj_scalar);
            isPreparedirtnufft = 1;
        else
            fprintf('-- MIRT      failed.\n');
            isPreparedirtnufft = 0;
        end
    catch
        fprintf('-- MIRT      failed.\n');
        isPreparedirtnufft = 0;
    end
    
    st.Nz = Nz;
       
    if isPreparedfinufft && ~reconOptions.flagUsefinufft && ~isPreparedirtnufft
        reconOptions.flagUsefinufft = 1;
    elseif isPreparedgpuNUFFT && ~(isPreparedfinufft || isPreparedirtnufft)
        reconOptions.flagUseGPU = 1;
    end

    if ~(isPreparedbartnufft || isPreparedgpuNUFFT || isPreparedfinufft || isPreparedirtnufft)
        error('All NUFFT package failed!');
    elseif reconOptions.flagUseBart
        fprintf('Use package: BART\n');
        reconOptions.flagUsefinufft = false;
    elseif reconOptions.flagUsefinufft && reconOptions.flagUseGPU && isPreparedcufinufft
        fprintf('Use package: cufinufft\n');
    elseif reconOptions.flagUsefinufft
        fprintf('Use package: finufft\n');
    elseif reconOptions.flagUseGPU && isPreparedgpuNUFFT
        fprintf('Use package: gpuNUFFT\n');
    else
        fprintf('Use package: MIRT\n');
    end
    st.om = single(st.om);
    reconOptions.st = st;
    reconOptions.isPreparedgpuNUFFT  = isPreparedgpuNUFFT;
    reconOptions.isPreparedfinufft   = isPreparedfinufft;
    reconOptions.isPreparedcufinufft = isPreparedcufinufft;
    reconOptions.isPreparedirtnufft  = isPreparedirtnufft;
    reconOptions.isPreparedbartnufft = isPreparedbartnufft;
    dataArray.flagIsTrajectorySet = true;
    dataArray.parOrder = parOrder;
end

dataArray.kspaceData = single(kspaceData);        
dataArray.navData    = single(navData);
        
