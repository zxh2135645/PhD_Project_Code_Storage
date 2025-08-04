function [params, reconOptions, dataArray] = calcfbp(params, reconOptions, dataArray)

reconOptions.flagUsefinufft = 1;

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(dataArray);

idx = 1;
% if Necho > 1
%     temp1 = fft(dataArray.navData(:,:,1,1,1),[],2);
%     temp2 = fft(dataArray.navData(:,:,1,1,2),[],2);
%     if mean(abs(temp2(:,1))) > mean(abs(temp1(:,1)))
%         idx = 2;
%     end
% end
kspaceData = dataArray.kspaceData(:,:,:,:,idx);
navData    = dataArray.navData(:,:,:,:,idx);
    
linOrder = st.linOrder_shift;
parOrder = st.parOrder_shift;

setupFunctions;

fprintf('Calculating back projection... ')

switch Trajectory
  case {'Radial','Spiral'}
%     if flagUseGPU
%         delete(gcp('nocreate'));
%         parpool('local',8);
%         gpuWorkerReset();
%     end
    
%     temp = imresize(phantom(N),[Ny Nx]);
%     [~,~,c] = svd([vec(repmat(temp,1,1,st.Nz)) vec(Ainv(A(repmat(temp,1,1,st.Nz),st),st,1))],'econ');
%     c = real(c(1)/c(2));
%     if ~(c > 0)
%         c = 1;
%     end
%     reconOptions.st.c = c; 
    
    reconOptions.st.c = 1;

    fprintf('Aligning k-space... ')
    [Phi_rt,~,~] = svde(navData(:,:));
    Phi_rt1      = Phi_rt(:,1)/mean(abs(Phi_rt(:,1)));
    if strcmp(ScanType,'Cine')
        Phi_rt1 = interp1(navIndices,Phi_rt1.',1:Ntpoint,'pchip','extrap').';
    elseif strcmp(ScanType,'CEST')
        Phi_rt1 = interp1Segmented(Phi_rt1,navIndices,linesPerShot*moduleLength,'cols');
    else
        Phi_rt1 = interp1Segmented(Phi_rt1,navIndices,linesPerShot,'cols');
    end
    Phi_rt1(navIndices) = []; 
    Phi_rt1 = single(Phi_rt1); 
    
    %%
    Ncoils = size(kspaceData,4);
    fbp_data   = zeros(Ntrajs, Nkx, Nz, Ncoils,'single');
    fbp_window = zeros(Ntrajs, Nkx, Nz);
%     t_ind = (linOrder==1) & (parOrder==1);
    tol = msdev*1e3; %1e-3; %norm(Phi_rt1(t_ind,:));
    for traj = 1:Ntrajs
        for np = 1:Nz
            t_ind = (linOrder==traj) & (parOrder==np);
            if sum(t_ind)
                fbp_data(traj,:,np,:) = reshape((Phi_rt1(t_ind,:)'*kspaceData(t_ind,:)/(norm(Phi_rt1(t_ind))^2+tol)),[],Ncoils);
                fbp_window(traj,:,np) = sum(abs(Phi_rt1(t_ind,:)).^2)/(norm(Phi_rt1(t_ind))^2+tol)^2;
            end
        end
    end
    fbp_data(~isfinite(fbp_data)) = 0;
   
    fprintf('inverse NUFFT... ')
    fbp = Ainv(fbp_data,st,reconOptions.st.c); 
    if MBfactor == 1
        fbp = fftshift(fbp,3);
    end
       
    fbpComposite = sqrt(sum(abs(fbp).^2,4));
    cw = prctile(fbpComposite(:),99.5);
    fbpComposite = fbpComposite/cw;
    if flagCommandLine
        implayZoom(imageOrientLPS(fbpComposite,params));
    end
    
    % SNR-weighted sampling mask in Cartesian coords
    fbp_window = 1./fbp_window;
    fbp_window(~isfinite(fbp_window)) = 0;
    fbp_window = Ah(fbp_window*1j, st);
    fbp_window = sqrt(abs(fftn(fftshift(fftshift(fbp_window,1),2))));
    
  %%
  case 'Cartesian' 
    %Calculate first basis image (i.e., a high-SNR static image)
    [Phi_rt,~,~] = svde(navData(:,:));
    Phi_rt1      = Phi_rt(:,1)/mean(abs(Phi_rt(:,1)));
    if strcmp(ScanType,'Cine')
        Phi_rt1 = interp1(navIndices,Phi_rt1.',1:Ntpoint,'pchip','extrap').';
    elseif strcmp(ScanType,'CEST')
        Phi_rt1 = interp1Segmented(Phi_rt1,navIndices,linesPerShot*moduleLength,'cols');
    else
        Phi_rt1 = interp1Segmented(Phi_rt1,navIndices,linesPerShot,'cols');
    end
    Phi_rt1(navIndices) = []; 
    Phi_rt1 = single(Phi_rt1);  
    t_ind = (linOrder==1) & (parOrder==1);
    tol = msdev*1e3; %1e-3; %norm(Phi_rt1)/numel(Phi_rt1);
    
    fbp_data = zeros(Ny, Nx, Nz, Ncoils,'single');
    for npy = 1:Ny
        for npz = 1:Nz
            t_ind = (linOrder==npy) & (parOrder==npz);
            if ~isempty(t_ind)
                fbp_data(npy,(Nx-size(kspaceData,2)+1):Nx,npz,:) = reshape((Phi_rt1(t_ind,:)'*kspaceData(t_ind,:)/(norm(Phi_rt1(t_ind))^2+tol)),[],Ncoils);
%                 fbp_data(npy,(Nx-size(kspaceData,2)+1):Nx,npz,:) = reshape(pinv(Phi_rt1(t_ind,:))*kspaceData(t_ind,:),[],Ncoils);
            end
        end
    end
    fbp_data(~isfinite(fbp_data)) = 0;
    
    fbp = Ainv(fbp_data,st);
    fbp = fftshift(fftshift(fbp,1),3);
    
%     wy = sgolayfilt(fftshift(squeeze(st.w(1,1,:))),3,floor(Ny/16)*2+1);
%     winy = wy>0.2*max(abs(wy));
%     wz = sgolayfilt(fftshift(squeeze(st.w(:,1,1))),3,floor(Nz/16)*2+1);
%     winz = wz>0.2*max(abs(wz));
%     fbp_window = repmat(reshape(double(winy)*double(winz)',Ny,1,Nz),[1 Nx 1]);
%     fbp_window(:,1:DC-floor(Nx/16)-1,:) = 0;
%     fbp_window(:,DC+floor(Nx/16):Nx,:) = 0;
    fbp_window = logical(repmat(reshape(st.w,Ny,1,Nz),[1 Nx 1]));
    
    fbpComposite = sqrt(sum(abs(fbp).^2,4));
    cw = prctile(fbpComposite(:),99.9);
    fbpComposite = fbpComposite/cw;

    if flagCommandLine
        implayZoom(imageOrientLPS(fbpComposite,params));
    end
    
end

Nzshift = floor(Nz/2)- floor(Nzorig/2);
roi_weighting = zeros(Ny,Nx,Nz);
roi_weighting(:,:,Nzshift+(1:Nzorig)) = repmat(hanning(Ny),[1 Nx Nzorig]); %ROI is center region of center slices

roi = reshape(bsxfun(@times,roi_weighting,fbp),[],Ncoils); %region of interest
roi = roi'*roi;

roni = reshape(bsxfun(@times,max(roi_weighting(:))-roi_weighting,fbp),[],Ncoils); %region of NO interest is outer region
roni = roni'*roni;

[~,cardl] = svde(roi/roni);
cards = sqrt(cardl);
dataArray.coilEnergycumsum = sqrt(cumsum(diag(cards).^2))/norm(diag(cards));

dataArray.fbp        = fbp;
dataArray.fbp_orig   = fbp;
dataArray.fbp_window = fbp_window;
dataArray.fbpComposite = fbpComposite;

reconOptions.interpFactor = [1 1];

fprintf(' done.\n')
