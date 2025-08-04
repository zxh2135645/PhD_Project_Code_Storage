function [params,reconOptions,dataArray,temporalBasis] = realtimeSubspace(params, reconOptions, dataArray, temporalBasis)

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(dataArray);

vec = @(x) x(:);

if flagCommandLine
    fprintf('Estimating real-time subspace ... ');
end

if ~flagIsDriftCorrected   
    if flagCommandLine
        fprintf('drift correction...')
    end
    for n = 1:1
        navDataTemp = dataArray.navData(:,:,:,:,n);
        [Phi_rt,S_rt,~] = svde(navDataTemp(:,:));
        L = min(find(diag(S_rt)/S_rt(1)>0.001,1,'last'), L_rt);
        Phi_rt_small = Phi_rt(:,1:L).';
    
        ph_corrs = linspace(-pi,pi,101)/Ntpoint;
        cost     = @(x) std(exp(-1i*x*navTiming).*sign(Phi_rt_small(1,:)));
        costs    = zeros(numel(ph_corrs),1);
        for j = 1:numel(ph_corrs)
            costs(j) = cost(ph_corrs(j));
        end
        [~,j]   = min(costs);
        ph_corr = fminsearch(cost,ph_corrs(j));
        ph_corr = exp(-1i*ph_corr*ACQTiming);
    
        ph_corr_k = ph_corr;
        ph_corr_k(navIndices) = [];
        
        dataArray.navData(:,:,:,:,n) = bsxfun(@times,dataArray.navData(:,:,:,:,n),vec(ph_corr(navIndices)));
    %     if isCartesian
    %         dataArray.kspaceData = bsxfun(@times,dataArray.kspaceData,ph_corr(:));
    %     else
              dataArray.kspaceData(:,:,:,:,n) = bsxfun(@times,dataArray.kspaceData(:,:,:,:,n),ph_corr_k(:));
    %     end
    end
    if flagCommandLine
        fprintf('done\n');
    end
    dataArray.flagIsDriftCorrected = true;
end

navData = dataArray.navData(:,:,:,:,1);

%%
[Phi_rt,S_rt,~] = svde(navData(:,:));
L = L_rt; %min(find(diag(S_rt)/S_rt(1)>0.01,1,'last'), L_rt);
fprintf('real-time recon rank = %d. Calculating temporal basis function...', L);
% if flagCommandLine
%     figure, subplot(2,1,1);plot(abs(fftshift(sqrt(sum(abs(fft(Phi_rt(:,3:L), [], 1)).^2, 2)))));title('unfiltered')
% end

% Apply low pass filter to Phi_rt(:, 3:end)
if ~isCartesian
    L_lp = 3; %input('Which rank to start low-pass filter [3]? ');
    if isempty(L_lp)
        L_lp = 3;
    end
    %fs = 1/(params.SGBlockLength*params.lEchoSpacing); 
    for j = L_lp:L
        Phi_rt(:,j) = sgolayfilt(double(Phi_rt(:,j)), 0, 3);
    end
%     if flagCommandLine
%         subplot(2,1,2);plot(abs(fftshift(sqrt(sum(abs(fft(Phi_rt(:,3:L), [], 1)).^2, 2)))));title('filtered')
%     end
end

%% generate temporal basis functions by interpolation
Phi_rt_small = Phi_rt(:,1:L).';

if strcmp(ScanType,'Cine')
    Phi_rt_full = interp1(navIndices,Phi_rt_small.',1:Ntpoint,'pchip','extrap').';
elseif strcmp(ScanType,'CEST')
    Phi_rt_full = interp1Segmented(Phi_rt_small,navIndices,linesPerShot*moduleLength,'rows','pchip',navTiming,ACQTiming);
else
    Phi_rt_full = interp1Segmented(Phi_rt_small,navIndices,linesPerShot,'rows','pchip',navTiming,ACQTiming);
end
Phi_rt = Phi_rt_full;

if isCartesian
    if strcmp(ScanType,'Cine')
        % k-space weighting functions with ifftshift
        st.w = zeros(Ny,1,Nz);
        for npy = 1:Ny
            for npz = 1:Nz
                t_ind = (st.linOrder_shift==npy) & (st.parOrder_shift==npz);
                st.w(npy,:,npz) = sum(t_ind(:));
            end
        end

        st.winv = 1./st.w;
        st.winv(st.w==0) = 0;

        params.st = st;  
    else    % not Cine
        Phi_rt(:,navIndices) = [];
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

        params.st = st;     
    end
else
    Phi_rt(:,navIndices) = [];
end

%% initial values
reconOptions.L_init = L;
temporalBasis.Phi_rt       = Phi_rt;
temporalBasis.Phi_rt_full  = Phi_rt_full;
temporalBasis.Phi_rt_small = Phi_rt_small;
temporalBasis.Phi_rt_init       = Phi_rt;
temporalBasis.Phi_rt_full_init  = Phi_rt_full;
temporalBasis.Phi_rt_small_init = Phi_rt_small;

if flagCommandLine
    plotPhi(Phi_rt_small,lEchoSpacing*SGBlock,1:4,'Check for respiration and cardiac frequencies.');
    
end

fprintf('done\n');

