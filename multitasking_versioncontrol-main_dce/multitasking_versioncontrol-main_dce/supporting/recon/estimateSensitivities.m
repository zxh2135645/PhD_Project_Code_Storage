function [params, reconOptions, dataArray] = estimateSensitivities(params, reconOptions, dataArray)

dispSlice = floor(params.Nz/2) + 1;

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(dataArray);

reconOptions.dispSlice = dispSlice;

linOrder = st.linOrder_shift;
parOrder = st.parOrder_shift;

setupFunctions;

fprintf('Estimating sensitivity maps... ');

newCoils = Ncoils;

switch Trajectory
  case 'Radial'
%     if flagUseGPU
%         delete(gcp('nocreate'));
%         parpool('local',8);
%         gpuWorkerReset();
%     end
    if ~isfield(reconOptions.st,'c')
        [~,~,c] = svd([vec(repmat(phantom(N),1,1,st.Nz)) vec(Ainv(A(repmat(phantom(N),1,1,st.Nz),st),st,1))],'econ');
        c       = real(c(1)/c(2));
        if ~(c > 0)
            c = 1;
        end
        reconOptions.st.c = c;
    end
    if ~isfield(dataArray,'fbp')
        idx = 1;
        if Necho > 1
            temp1 = fft(dataArray.navData(:,:,1,1,1),[],2);
            temp2 = fft(dataArray.navData(:,:,1,1,2),[],2);
            if mean(abs(temp2(:,1))) > mean(abs(temp1(:,1)))
                idx = 2;
            end
        end
        kspaceData = dataArray.kspaceData(:,:,:,:,idx);
        navData    = dataArray.navData(:,:,:,:,idx);

        [Phi_rt,~,~] = svde(navData(:,:));
        Phi_rt1      = Phi_rt(:,1)/mean(abs(Phi_rt(:,1)));
        Xq           = 1:Ntpoint; Xq(navIndices)=[];
        Phi_rt1      = interp1(navIndices,Phi_rt1,Xq,'pchip',0).';

        t_ind = (linOrder==1) & (parOrder==1);
        tol = msdev*1e3; %norm(Phi_rt1(t_ind,:));

        Ncoils = size(kspaceData,4);
        fbp_data = zeros(Ntrajs, N, Nz, Ncoils);
        for traj = 1:Ntrajs
            for np = 1:Nz
                t_ind = (linOrder==linOrder(traj)) & (parOrder==np);
                fbp_data(traj,:,np,:) = reshape((Phi_rt1(t_ind,:)'*kspaceData(t_ind,:)/(norm(Phi_rt1(t_ind))^2+tol)),[],Ncoils);
            end
        end
        fbp_data(~isfinite(fbp_data)) = 0;

        tic;
        fbp = Ainv(fbp_data,st,st.c); 
        if MBfactor == 1
            fbp = fftshift(fbp,3);
        end
        toc;
        
        fbpComposite = sqrt(sum(abs(fbp).^2,4));
        fbpComposite = fbpComposite/max(abs(fbpComposite(:)));
        dataArray.fbpComposite = fbpComposite;
        dataArray.fbp = fbp;
    end
  case 'Cartesian'
    if ~isfield(dataArray,'fbp')
        kspaceData = dataArray.kspaceData(:,:,:,:,1);
        navData    = dataArray.navData(:,:,:,:,1);

        %Calculate first basis image (i.e., a high-SNR static image)
        [Phi_rt_temp,~,~] = svde(navData(:,:));
        Phi_rt1      = Phi_rt_temp(:,1)/mean(abs(Phi_rt_temp(:,1)));
        Xq           = 1:Nread;
        Phi_rt1      = interp1(navIndices,Phi_rt1,Xq,'pchip',0).';
        Phi_rt1(navIndices) = [];   % set nav lines to 0
        
        t_ind = (linOrder==1) & (parOrder==1);
        tol = msdev*1e3; %norm(Phi_rt1(t_ind,:))*5;
        
        Ncoils = size(kspaceData,4);
        fbp_data = zeros(Ny, Nx, Nz, Ncoils);
        for npy = 1:Ny
            for npz = 1:Nz
                t_ind = (linOrder==npy) & (parOrder==npz);
                fbp_data(npy,(Nx-size(kspaceData,2)+1):Nx,npz,:) = reshape((Phi_rt1(t_ind,:)'*kspaceData(t_ind,:)/(norm(Phi_rt1(t_ind))^2+tol)),[],Ncoils);
            end
        end
        fbp_data(~isfinite(fbp_data)) = 0;

        fbp = Ainv(fbp_data,st);
        fbp = fftshift(fftshift(fbp,1),3);
        
%         wy = sgolayfilt(fftshift(squeeze(st.w(1,1,:))),3,floor(Ny/16)*2+1);
%         winy = wy>0.2*max(abs(wy));
%         wz = sgolayfilt(fftshift(squeeze(st.w(:,1,1))),3,floor(Nz/16)*2+1);
%         winz = wz>0.2*max(abs(wz));
%         fbp_window = repmat(reshape(double(winy)*double(winz)',Ny,1,Nz),[1 Nx 1]);
%         fbp_window(:,1:DC-floor(Nx/16)-1,:) = 0;
%         fbp_window(:,DC+floor(Nx/16):Nx,:) = 0;
        fbp_window = logical(repmat(reshape(st.w,Ny,1,Nz),[1 Nx 1]));
        
        fbpComposite = sqrt(sum(abs(fbp).^2,4));
        fbpComposite = fbpComposite/max(abs(fbpComposite(:)));
        dataArray.fbpComposite = fbpComposite;
        dataArray.fbp = fbp;
    end
end

%% estimate SEs
fprintf('Method: ');

if strcmp(SEmethod,'ESPIRiT') && ~flagUseBart && Nz > 1
    fprintf('3D image, switch to ');
    if exist('bart','file') && Nz > ceil(Ny/64)+2
        flagUseBart = true;
    else
        SEmethod = 'CBD';
    end
end

switch SEmethod
    case 'Walsh'
        fprintf('Walsh\n');
        window_length = 25;
        [SEs, img] = sensemaps_walsh(fbp,window_length,'acs',fbp_window,'figures',0,'voxelSpacing',voxelSpacing);
        %SEs = squeeze(ifftshift(sensemaps_walsh(bsxfun(@times,fftshift(fbp,3),reshape((-1).^(1:Nz),1,1,[])),'figures','acs',fbp_window),3)); % Use walsh method to estimate the sensitivities
    case 'ESPIRiT'
        if flagUseBart
            fprintf('BART ESPIRiT\n');
            calib_size  = [floor(Ny/8)-4, floor(Nx/8)-4, max(floor(Nz/8)-2,ceil(Ny/64)+2)];
            kernel_size = 8;%ceil(Ny/64)+2;
            newNy = 2^ceil(log2(Ny)); shiftNy = floor((newNy-Ny)/2);
            newNx = 2^ceil(log2(Nx)); shiftNx = floor((newNx-Nx)/2);
            newNz = 2^ceil(log2(Nz)); shiftNz = floor((newNz-Nz)/2);
            fbptemp = zeros(Ny,newNx,Nz,size(fbp,4));
            fbptemp(:,shiftNx+(1:Nx),:,:) = fbp;
            fbp_k = single(ft1d(ft1d(ft1d(fbptemp,1),2),3));
            fbp_k = permute(fbp_k, [2 1 3 4]);
            kernel_size = permute(kernel_size,[2 1 3]);
            calib_size = permute(calib_size, [2 1 3]);
            evalstring = ['[SEs, emaps] = bart(''ecalib -m 1 -r ' num2str(calib_size(1)),':',num2str(calib_size(2)),':',num2str(calib_size(3)),' -k ' num2str(kernel_size(1)) ''', fbp_k);'];
            %evalstring = ['[SEs, emaps] = bart(''ecalib -m 1 -r 24:24:12 -k 6'', fbp_k);'];
            tic;eval(evalstring);toc;
            SEs = permute(SEs,[2 1 3 4 5]);
            SEs = crop(SEs(:,:,:,:,1),st);
        else
            fprintf('ESPIRiT\n');
            calib_size = [floor(Ny/8)-4, floor(Nx/8)-4];
            kernel_size = [ceil(Ny/64)+2, ceil(Nx/64)+2];
            SEs = sensemaps_ESPIRiT(squeeze(fbp), calib_size, kernel_size);
            SEs = reshape(SEs,size(fbp));
        end
    case 'CBD'
        fprintf('CBD\n');
%         kernel_size = [3 3 3];
%        [SEs, ~] = sensemaps_cbd(fbp,kernel_size,'acs',fbp_window); 
        kernel_size = [5 5 5];
        [SEs, ~] = sensemaps_cbd_old(fbp,kernel_size,'acs',fbp_window); 
end

SE_corr = reshape(sqrt(sum(abs(reshape(SEs,[],size(SEs,4))*sqrtm(Psi)).^2,2)),Ny,Nx,Nz);
SEs     = bsxfun(@rdivide,SEs,SE_corr);
SEs(~isfinite(SEs)) = 0;
SEs(SEs==0) = min(SEs(SEs~=0));

if ~strcmp(ScanType,'Cine')
    %Bias correction
    pd        = abs(sum(conj(SEs).*fbp,4)./sum(abs(SEs).^2,4));
    pd(~isfinite(pd)) = 0;
    pd(pd==0) = min(pd(pd~=0));
    
    try 
        % SimpleITK N4 bias field correction
        % https://github.com/mikecjz/N4BiasCorrection-Matlab
        fprintf('N4 bias field estimation... ');
        [~,bias] = N4BiasCorrection(abs(pd));
        bias = exp(bias);
        mask2 = estimateMask(fbp,SEs,rawVoxelSpacing);
        bias = bias + mask2/2;
        fprintf('done \n');
    catch errormsg
        fprintf('%s\n', errormsg.message);
        fprintf('N4 Correction failed. Switching to polynomial fitting bias field estimation... \n');

        [x,y,z]        = ndgrid(-Ny/2:Ny/2-1,-Nx/2:Nx/2-1,-Nz/2:Nz/2-1);
        polyterm       = [x(:).^2, y(:).^2, z(:).^2, x(:).*y(:), x(:).*z(:), y(:).*z(:), x(:), y(:), z(:), ones(size(x(:)))]; 
        scalar = norm(polyterm(:));
        polyterm = polyterm/scalar;
        [polyterm,~,~] = svd(bsxfun(@times,polyterm,pd(:).^2),'econ');
        w   = polyterm'*pd(:);
        y   = zeros(size(pd(:)));
        rho = 1/max(abs(polyterm*w-pd(:)));
        for it = 1:10
            res = polyterm*w-pd(:);
            z   = res+y/rho;
            z   = sign(z).*max(abs(z)-1/rho,0);
            y   = y+rho*(res-z);
            rho = rho*1.05;
            w   = polyterm'*(z+pd(:)-y/rho);
            %     norm(polyterm*w-pd(:),1)
        end
        [~,~,v] = svd([(polyterm*w)./pd(:), pd(:)],'econ');
        w       = w*v(2)/v(1);

        bias    = reshape(polyterm*w,Ny,Nx,Nz)./pd.^2;
        bias(1./bias<1/2) = 2;
        bias(bias<1/4)   = 1/4;
        mask2 = estimateMask(fbp,SEs,rawVoxelSpacing);
        bias = bias + mask2/2;
    %     sigma = [3 3 3];
    %     sigma(1) = floor(Ny/16)*2 + 1;
    %     sigma(2) = floor(sigma(1)*voxelSpacing(1)/voxelSpacing(2)/2)*2 + 1;
    %     sigma(3) = floor(sigma(1)*voxelSpacing(1)/voxelSpacing(3)/2)*2 + 1;
    %     filterSize = sigma*4 + 1;
    %     bias = imgaussfilt3(bias,sigma/4,'Padding','symmetric','FilterSize',floor(filterSize/4));

    %      bias = estimateBiasPoly(fbp,SEs,voxelSpacing);
    end
    if flagCommandLine
        cw1 = prctile(abs(pd(:)),99.9);
        cw2 = prctile(abs(pd(:)./bias(:)),99.9);
        implayZoom(1.2*abs([imageOrientLPS(pd,params)/cw1 imageOrientLPS(pd./bias,params)/cw2]),5,[],'[Combined image / Bias corrected image]');
        figure;imagesc(imageOrientLPS(bias(:,:,dispSlice),params),[0 2]);title('bias field');axis equal tight off;
    end
    figure;montage(imageOrientLPS(abs(SEs(:,:,dispSlice,:)),params)/max(abs(SEs(:))));title('Sensitivity Maps (Magnitude)');drawnow;
    SEs = bsxfun(@times,SEs,bias);
    dataArray.bias = bias;
end

SEs(~isfinite(SEs)) = 0;
SEs(SEs==0) = min(SEs(SEs~=0));

figure;montage(imageOrientLPS(abs(SEs(:,:,dispSlice,:)),params)/max(abs(SEs(:))));title('Sensitivity Maps (Magnitude, bias corrected)');drawnow;
figure;montage(imageOrientLPS(angle(SEs(:,:,dispSlice,:)),params)/2/pi+0.5);colormap(hsv);title('Sensitivity Maps (Phase)');drawnow;
 
pd = sum(conj(SEs).*dataArray.fbp,4)./sum(abs(SEs).^2,4);
pd = pd./max(abs(pd(:)));
figure,subplot(1,2,1),imshow(1.2*imageOrientLPS(abs(pd(:,:,dispSlice)),params),[]),title('Coil combined image (Magnitude)');
subplot(1,2,2),imshow(imageOrientLPS(angle(pd(:,:,dispSlice)),params),[-pi pi]),title('Coil combined image (Phase)');drawnow;
   
dataArray.SEs = SEs;
dataArray.fbpComposite = pd;

% dataArray.kspaceData_orig = kspaceData;
% dataArray.navData_orig    = navData;

% if isfield(dataArray,'fbp')
%     dataArray = rmfield(dataArray,{'fbp'});
% end



