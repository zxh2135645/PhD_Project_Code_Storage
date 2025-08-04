function [params, reconOptions, dataArray] = coilCompression(params, reconOptions, dataArray, selectROI, roi_weighting)

if nargin < 4 || exist('roi_weighting','var')
    selectROI = false;
end

newCoils = 12;

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(dataArray);

if ~isfield(reconOptions,'dispSlice') || (dispSlice < 1 || dispSlice > Nz)
    dispSlice = floor(Nz/2) + 1;
end

if ~isfield(reconOptions,'useNewCoilsAuto') 
    useNewCoilsAuto = false;
end

kspaceData = dataArray.kspaceData;
navData    = dataArray.navData;
Ncoils     = size(navData,4);

linOrder = st.linOrder_shift;
parOrder = st.parOrder_shift;

setupFunctions;

newCoils = min(Ncoils,newCoils);

Nydisp = Ny;
Nxdisp = Nx;
% Define window
if selectROI
    temp = sqrt(sum(abs(fbp(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp),:,:)).^2,4));
    h = figure;imagesc(sqrt(abs(temp(:,:,floor(Nz/2)+1,1))));axis equal tight;colormap('gray');title('Draw ROI')
    roiResp = imrect;
    roiPosition = roiResp.getPosition();
    roiPosition = round(roiPosition);
    close(h);
    Nxshift = max(roiPosition(1),0);
    Nyshift = max(roiPosition(2),0);
    ROINx   = min(roiPosition(3),(Nxdisp-Nxshift));
    ROINy   = min(roiPosition(4),(Nydisp-Nyshift));
    Nxend   = Nxshift + ROINx;
    Nyend   = Nyshift + ROINy;

    if Nz > 1
        ROIy = ones(Nydisp,1); ROIy(Nyshift:(Nyshift+ROINy),1) = 3;
        h = figure;imagesc(sqrt(squeeze(abs(temp(:,Nxshift+floor(ROINx/2),:,1).*ROIy))'));axis equal tight;colormap('gray');title('Draw heart ROI')
        roiResp = imrect;
        roiPosition = roiResp.getPosition();
        roiPosition = round(roiPosition);
        close(h);
        Nyshift = max(min(roiPosition(1),Nyshift),0);
        Nzshift = max(roiPosition(2),0);
        Nyend   = min(max(Nyshift+roiPosition(3),Nyend),Nydisp); 
        ROINy   = Nyend - Nyshift;    ROINy = floor(ROINy/2)*2;
        ROINz   = min(roiPosition(4),(Nz-Nzshift));
        Nzend   = Nzshift + ROINz;

        ROIxz = ones(1,Nxdisp,Nz); ROIxz(1,Nxshift+(1:ROINx),Nzshift+(1:ROINz)) = 3;
        h = figure;imagesc(sqrt(squeeze(abs(temp(Nyshift+floor(ROINy/2),:,:,1).*ROIxz))));axis equal tight;colormap('gray');title('Draw heart ROI')
        roiResp = imrect;
        roiPosition = roiResp.getPosition();
        roiPosition = round(roiPosition);
        close(h);
        Nzshift = max(min(roiPosition(1),Nzshift),0);
        Nxshift = max(min(roiPosition(2),Nxshift),0);
        Nzend   = min(max(Nzshift+roiPosition(3),Nzend),Nz);
        Nxend   = min(max(Nxshift+roiPosition(4),Nxend),Nxdisp); 
        ROINz   = Nzend - Nzshift;    ROINz = floor(ROINz/2)*2;
        ROINx   = Nxend - Nxshift;    ROINx = floor(ROINx/2)*2;
    else
        Nzshift = 0;
        ROINz   = 1;
    end

    windowfun = zeros(Nydisp,Nxdisp,Nz);
    windowfun(Nyshift+(1:ROINy),Nxshift+(1:ROINx),:) = repmat(hanning(ROINy)*hanning(ROINx).',[1 1 Nz]);
    windowZ = zeros(Nz,1);
    windowZ(Nzshift+(1:ROINz)) = hanning(ROINz);
    windowfun = sqrt(sqrt(bsxfun(@times,windowfun,reshape(windowZ,1,1,Nz))));
    windowfun = windowfun/max(windowfun(:));
    roi_weighting = zeros(Ny,Nx,Nz);
    roi_weighting(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp),:) = windowfun;
else
    Nydisp = floor(Nydisp/2)*2;
    Nxdisp = floor(Nxdisp/2)*2;
    Nzshift = floor(Nz/2) - floor(Nzorig/2);
    Nyshift = floor(Ny/2) - floor(Nydisp/2);
    Nxshift = floor(Nx/2) - floor(Nxdisp/2);
    roi_weighting = zeros(Ny,Nx,Nz);
    roi_weighting(Nyshift+(1:Nydisp),Nxshift+(1:Nxdisp),Nzshift+(1:Nzorig)) = repmat(hanning(Nydisp)*hanning(Nxdisp).',[1 1 Nzorig]); %ROI is center region of center slices 
%     roi_weighting = ones(Ny,Nx,Nz);
end
    
roi_mask  = double(roi_weighting>0.05);
roni_mask = double(roi_weighting<0.02);

dataArray.roi_weighting = roi_weighting;

switch Trajectory
  case {'Radial','Spiral'}
   
    %% Coil compression for faster computation
    
    %roi = reshape(bsxfun(@times,roi_weighting,fbp),[],Ncoils); %region of interest
    roi = reshape(bsxfun(@times,roi_mask,fbp),[],Ncoils); %region of interest
    roi = roi'*roi;

    %roni = reshape(bsxfun(@times,1-roi_weighting,fbp),[],Ncoils);   % region of NO interest is outer region
    roni = reshape(bsxfun(@times,roni_mask,fbp),[],Ncoils); 
    roni = roni'*roni;

    if selectROI
        [cardmixer,cardl,V] = svde(roi/roni);
    else
        [cardmixer,cardl,V] = svde(roi);
    end
    cards = sqrt(cardl); 

    fbptemp  = reshape(reshape(fbp,[],Ncoils)*cardmixer,size(fbp(:,:,:,1:Ncoils)));
    sig_roi  = squeeze(sqrt(sum(sum(sum(abs(fbptemp.*(roi_weighting>0.05)).^2,1),2),3)));
    sig_roni = squeeze(sqrt(sum(sum(sum(abs(fbptemp.*(roi_weighting<0.02)).^2,1),2),3)));

    coilEnergycumsum = sqrt(cumsum(sig_roi.^2));
    coilEnergycumsum = coilEnergycumsum/coilEnergycumsum(end);
    %coilEnergycumsum = sqrt(cumsum(diag(cards).^2))/norm(diag(cards));
    newCoilsAuto = find(coilEnergycumsum>=mincoilEnergy,1);
   
    RONIEnergycumsum = sqrt(cumsum(sig_roni.^2));
    RONIEnergycumsum = RONIEnergycumsum/RONIEnergycumsum(end);
    
%     [~,temp2] = kmeans(sig_roi,2);
%     newCoils = find(sig_roi<mean(temp2)/(lambertw(-2*exp(-2)) + 2),1);
   
    if isempty(newCoils) || newCoils < minNewCoils
        newCoils = minNewCoils;
    end
    
    if useNewCoilsAuto
        newCoils = max(min(newCoils,Ncoils),newCoilsAuto);
    else
        newCoils = min(newCoils,Ncoils);
    end

    fprintf('First %d virtual coils retain %.2f %% ROI signal energy.\n', newCoils, 100*coilEnergycumsum(newCoils));

    %cards = sqrt(cardl);
    dataArray.roi_weighting = roi_weighting;
    dataArray.cardmixer = cardmixer;
    dataArray.coilEnergycumsum = coilEnergycumsum;
    
    if flagCommandLine
        cw = prctile(abs(fbptemp(:)),99.5);
        titleText = sprintf('All virtual coil images. Fist %d coils retain %.2f %% of ROI energy', newCoils,100*coilEnergycumsum(newCoils));
        figure;montage(squeeze(imageOrientLPS(abs(fbptemp(:,:,dispSlice,:)/cw),params)));colormap(jet);title(titleText);drawnow;
        temp = sqrt(cumsum(abs(fbptemp).^2,4));
        temp = temp(:,:,dispSlice,newCoils);
        cw = prctile(abs(temp(:)),99.5);
        titleText = sprintf('Combined image of first %d virtual coils, %.2f %% ROI signal energy retained.', newCoils, 100*coilEnergycumsum(newCoils));
        figure;montage(squeeze(imageOrientLPS(abs(temp)/cw,params)));colormap(jet);title(titleText);drawnow;
    end
    
    fprintf('Coil compression: from %d to %d coils... ', Ncoils,newCoils);

    %fbp = fbptemp;
	fbp = fbptemp(:,:,:,1:newCoils);
    kspaceData = reshape(permute(kspaceData,[5 2 1 3 4]),[],Ncoils);
    kspaceData = permute(reshape(kspaceData*cardmixer(1:Ncoils,1:newCoils),Necho,Nkx,[],1,newCoils),[3 2 4 5 1]);
    NechoNav   = size(navData,5);
    Ntemp      = size(navData,2);
    navData    = reshape(permute(navData,[5 2 1 3 4]),[],Ncoils);
    navData    = permute(reshape(navData*cardmixer(1:Ncoils,1:newCoils),NechoNav,Ntemp,[],1,newCoils),[3 2 4 5 1]);
    dataArray.Psi = cardmixer(1:Ncoils,1:newCoils)'*dataArray.Psi*cardmixer(1:Ncoils,1:newCoils);
    %dataArray.Psi = cardmixer'*dataArray.Psi*cardmixer;
    if exist('SEs','var') && size(SEs,4) == Ncoils
        SEs = reshape(reshape(SEs,[],Ncoils)*cardmixer(1:Ncoils,1:newCoils),size(SEs,1),size(SEs,2),size(SEs,3),newCoils);
        dataArray.SEs = SEs;
        fbpComposite = abs(sum(fbp.*conj(SEs),4));
    else
        fbpComposite = sqrt(sum(abs(fbp).^2,4));
    end   
    cw = prctile(fbpComposite(:),98);
    fbpComposite = fbpComposite/cw;
    if flagCommandLine
        titleText = sprintf('New coil number = %d, %.2f %% ROI signal energy retained.', newCoils, 100*coilEnergycumsum(newCoils));
        figure;montage(squeeze(imageOrientLPS(abs(fbp(:,:,dispSlice,:)/cw),params)));title(titleText);drawnow;
    end
    
    Ncoils     = newCoils;
    params.mixer = eye(Ncoils);
    
    dataArray.navData          = navData;
    dataArray.kspaceData       = kspaceData;
    dataArray.flagIsCompressed = true;
    
  %%
  case 'Cartesian'

    fprintf('Coil compression: from %d to %d coils... ', Ncoils,newCoils);

    % Spatially-varying coil compression for faster computation
    % Align coil compression matrices based on nearest spanning vectors in
    % subspaces based on Zhang et. al MRM 2013;69(2):571-82.
    
    NechoNav   = size(navData,5);
    
    kspaceData = permute(kspaceData,[5 2 1 3 4]);
    navData    = permute(navData,[5 2 1 3 4]);
     
    roi_weighting  = roi_weighting.*roi_mask;
    roni_weighting = (1-roi_weighting).*roni_mask;

    % k-space center
    DC = floor(Nx/2) + 1;
    [~,temp,CoilV] = svde(squeeze(kspaceData(1,DC,:,:,:)));
    CoilV = CoilV(:,1:newCoils);      % Compression matrix for DC of readout dimension
    
%     roi  = reshape(bsxfun(@times,roi_weighting(:,DC,:),fbp(:,DC,:,:)),[],Ncoils); 
%     roi  = roi'*roi;
%     roni = reshape(bsxfun(@times,roni_weighting(:,DC,:),fbp(:,DC,:,:)),[],Ncoils); 
%     roni = roni'*roni;
%     if sum(abs(roi(:))~=0) && sum(abs(roni(:))~=0)
%         [CoilV,cardl,V] = svde(roi/roni);
%     elseif sum(roi(:)~=0)
%         [CoilV,cardl,V] = svde(roi);
%     else
%         [~,~,CoilV] = svde(squeeze(kspaceData(1,DC,:,:,:)));
%     end
%     CoilV = CoilV(:,1:newCoils);

    kspaceData(:,DC,:,1,1:newCoils) = reshape(reshape(kspaceData(:,DC,:,:,:),[],Ncoils)*CoilV,Necho,1,[],1,newCoils);   % Replace first NCha coils with compressed version
    navData(:,DC,:,1,1:newCoils)    = reshape(reshape(navData(:,DC,:,:,:),[],Ncoils)*CoilV,NechoNav,1,[],1,newCoils);   % Replace first NCha coils with compressed version
    Psi_new(:,:,DC)                 = CoilV'*Psi*CoilV;  % update noise matrix
    CoilV_orig = CoilV;

    % Forward from center
    CoilV_last = CoilV_orig;
    for j = DC+1:Nx
        [~,~,CoilV] = svde(squeeze(kspaceData(1,j,:,:,:)));
        CoilV = CoilV(:,1:newCoils);
        
%         roi  = reshape(bsxfun(@times,roi_weighting(:,j,:),fbp(:,j,:,:)),[],Ncoils); 
%         roi  = roi'*roi;
%         roni = reshape(bsxfun(@times,roni_weighting(:,j,:),fbp(:,j,:,:)),[],Ncoils); 
%         roni = roni'*roni;
%         if sum(abs(roi(:))~=0) && sum(abs(roni(:))~=0)
%             [CoilV,cardl,V] = svde(roi/roni);
%         elseif sum(roi(:)~=0)
%             [CoilV,cardl,V] = svde(roi);
%         else
%             [~,~,CoilV] = svde(squeeze(kspaceData(1,j,:,:,:)));
%         end
%         CoilV = CoilV(:,1:newCoils);

        [U,~,V] = svd(CoilV_last'*CoilV,'econ');
        CoilV   = CoilV*V*U'; %Compression matrix for index j (as close as possible to preceding compressor)
        kspaceData(:,j,:,1,1:newCoils) = reshape(reshape(kspaceData(:,j,:,:,:),[],Ncoils)*CoilV,Necho,1,[],1,newCoils); % Replace first NCha coils with compressed version
        navData(:,j,:,1,1:newCoils)    = reshape(reshape(navData(:,j,:,:,:),[],Ncoils)*CoilV,NechoNav,1,[],1,newCoils); % Replace first NCha coils with compressed version
        Psi_new(:,:,j) = CoilV(:,1:newCoils)'*Psi*CoilV(:,1:newCoils); %update noise matrix
        CoilV_last     = CoilV;
    end

    % Backward from center
    CoilV_last = CoilV_orig;
    for j = DC-1:-1:1
        [~,~,CoilV] = svde(squeeze(kspaceData(1,j,:,:,:)));
        CoilV = CoilV(:,1:newCoils);

%         roi  = reshape(bsxfun(@times,roi_weighting(:,j,:),fbp(:,j,:,:)),[],Ncoils); 
%         roi  = roi'*roi;
%         roni = reshape(bsxfun(@times,roni_weighting(:,j,:),fbp(:,j,:,:)),[],Ncoils); 
%         roni = roni'*roni;
%         if sum(abs(roi(:))~=0) && sum(abs(roni(:))~=0)
%             [CoilV,cardl,V] = svde(roi/roni);
%         elseif sum(roi(:)~=0)
%             [CoilV,cardl,V] = svde(roi);
%         else
%             [~,~,CoilV] = svde(squeeze(kspaceData(1,j,:,:,:)));
%         end
%         CoilV = CoilV(:,1:newCoils);

        [U,~,V] = svd(CoilV_last'*CoilV,'econ');
        CoilV   = CoilV*V*U'; % Compression matrix for index j (as close as possible to preceding compressor)
        kspaceData(:,j,:,1,1:newCoils) = reshape(reshape(kspaceData(:,j,:,:,:),[],Ncoils)*CoilV,Necho,1,[],1,newCoils); % Replace first NCha coils with compressed version
        navData(:,j,:,1,1:newCoils)    = reshape(reshape(navData(:,j,:,:,:),[],Ncoils)*CoilV,NechoNav,1,[],1,newCoils); % Replace first NCha coils with compressed version
        Psi_new(:,:,j) = CoilV(:,1:newCoils)'*Psi*CoilV(:,1:newCoils); %update noise matrix
        CoilV_last     = CoilV;
    end

    % Replace data with compressed version
    %kspaceData = ift1d(kspaceData(:,:,:,:,1:newCoils),2); % keep only compressed coils;
    %navData    = ift1d(navData(:,:,:,:,1:newCoils),2);
    kspaceData = kspaceData(:,:,:,:,1:newCoils); % keep only compressed coils;
    navData    = navData(:,:,:,:,1:newCoils);
    dataArray.Psi = mean(Psi_new,3);    % average noise matrix
    Ncoils     = newCoils;

    kspaceData = permute(kspaceData,[3 2 4 5 1]);
    navData    = permute(navData,[3 2 4 5 1]);
    
    dataArray.navData    = navData;
    dataArray.kspaceData = kspaceData;
    dataArray.flagIsCompressed = true;
    
    setupFunctions;     % redo setup

    kspaceData = kspaceData(:,:,:,:,1);
    navData    = navData(:,:,:,:,1);
    
    % Calculate first basis image (i.e., a high-SNR static image)
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
    tol = msdev*1e3; %norm(Phi_rt1)/numel(Phi_rt1);%sqrt(sum(abs(Phi_rt1(t_ind,:)).*2)/numel(t_ind));
    
    fbp_data = zeros(Ny, Nx, Nz, Ncoils);
    for npy = 1:Ny
        for npz = 1:Nz
            t_ind = (linOrder==npy) & (parOrder==npz);
            fbp_data(npy,(Nx-size(kspaceData,2)+1):Nx,npz,:) = reshape((Phi_rt1(t_ind,:)'*kspaceData(t_ind,:)/(norm(Phi_rt1(t_ind))^2+tol)),[],Ncoils);
%             fbp_data(npy,(Nx-size(kspaceData,2)+1):Nx,npz,:) = reshape(pinv(Phi_rt1(t_ind,:))*kspaceData(t_ind,:),[],Ncoils);
        end
    end
    fbp_data(~isfinite(fbp_data)) = 0;
    
    fbp = Ainv(fbp_data,st);
    fbp = fftshift(fftshift(fbp,1),3);
    
    fbpComposite = sqrt(sum(abs(fbp).^2,4));
    cw = prctile(fbpComposite(:),99);
    fbpComposite = fbpComposite/cw;
    if flagCommandLine
        titleText = sprintf('Magnitude image: New virtual coil number = %d',Ncoils);
        montage(squeeze(imageOrientLPS(abs(fbp(:,:,dispSlice,:)/cw),params)));title(titleText);drawnow;
    end
end

params.Ncoils = Ncoils;
params.mixer  = eye(Ncoils);

dataArray.fbp = fbp;
dataArray.fbp_window   = fbp_window;
dataArray.fbpComposite = fbpComposite;

fprintf(' done.\n')
