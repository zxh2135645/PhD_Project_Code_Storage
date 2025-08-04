function [params, reconOptions, dataArray] = coilCompressionROVir(params, reconOptions, dataArray, selectROI, ROIweighting, issvROVir)

mincoilEnergy = 0.97;
slwin = 1;

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(dataArray);

if nargin < 6
    issvROVir = 0;
end

if nargin >= 5 
    if size(ROIweighting,1) == Ny && size(ROIweighting,2) == Nx && size(ROIweighting,3) == Nz
        roi_weighting = ROIweighting;
    else
        roi_weighting = [];
    end
else
    roi_weighting = [];
end

if nargin < 4
    selectROI = false;
end

if ~isfield(reconOptions,'dispSlice') || (dispSlice < 1 || dispSlice > Nz)
    dispSlice = floor(Nz/2) + 1;
end

kspaceData = dataArray.kspaceData;
navData    = dataArray.navData;
Ncoils     = size(navData,4);

linOrder = st.linOrder_shift;
parOrder = st.parOrder_shift;

setupFunctions;

newCoils = min(Ncoils,newCoils);

if strcmp(Trajectory,'Cartesian')
%     kspaceData = ft1d(kspaceData,2);
%     navData    = ft1d(navData,2);
    NechoNav   = size(navData,5);
end

% Define window
if isempty(roi_weighting)
    if selectROI
        temp = sqrt(sum(abs(fbp(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp),:,:)).^2,4));
        h = figure;imagesc(sqrt(abs(temp(:,:,floor(Nz/2)+1,1))));axis equal tight;colormap('gray');title('Draw ROI')
        roiResp = imrect;
        roiPosition = roiResp.getPosition();
        roiPosition = round(roiPosition);
        close(h);
        Nxshift = roiPosition(1);
        Nyshift = roiPosition(2);
        ROINx   = roiPosition(3);
        ROINy   = roiPosition(4);

        if Nz > 1
            ROIy = ones(Nydisp,1); ROIy(Nyshift:(Nyshift+ROINy),1) = 3;
            h = figure;imagesc(sqrt(squeeze(abs(temp(:,Nxshift+floor(ROINx/2),:,1).*ROIy))'));axis equal tight;colormap('gray');title('Draw heart ROI')
            roiResp = imrect;
            roiPosition = roiResp.getPosition();
            roiPosition = round(roiPosition);
            close(h);
            Nyshift = min(roiPosition(1),Nyshift);
            Nzshift = roiPosition(2);
            ROINy   = max(roiPosition(3),ROINy);  ROINy = floor(ROINy/2)*2;
            ROINz   = roiPosition(4);

            ROIxz = ones(1,Nxdisp,Nz); ROIxz(1,Nxshift:(Nxshift+ROINx),Nzshift:(Nzshift+ROINz)) = 3;
            h = figure;imagesc(sqrt(squeeze(abs(temp(Nyshift+floor(ROINy/2),:,:,1).*ROIxz))));axis equal tight;colormap('gray');title('Draw heart ROI')
            roiResp = imrect;
            roiPosition = roiResp.getPosition();
            roiPosition = round(roiPosition);
            close(h);
            Nzshift = min(roiPosition(1),Nzshift);
            Nxshift = min(roiPosition(2),Nxshift);
            ROINz   = max(roiPosition(3),ROINz); ROINz = floor(ROINz/2)*2;
            ROINx   = max(roiPosition(4),ROINx); ROINx = floor(ROINx/2)*2;
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
    end
end

roi_mask  = double(roi_weighting>0.05);
roni_mask = double(roi_weighting<0.02);

dataArray.roi_weighting = roi_weighting;

% if flagCommandLine
%     cw = prctile(abs(fbp(:)),99.5);
%     titleText = 'fbp';
%     figure;montage(squeeze(imageOrientLPS(abs(fbp(:,:,dispSlice,:)/cw),params)));colormap(jet);title(titleText);drawnow;
%     figure;montage(squeeze(imageOrientLPS(angle(fbp(:,:,dispSlice,:))/2/pi+0.5,params)));colormap(jet);title(titleText);drawnow;
%     temp = sqrt(sum(abs(fbp(:,:,dispSlice,:)).^2,4));
%     cw = prctile(abs(temp(:)),99.5);
%     titleText = sprintf('Combined image of %d coils.', Ncoils);
%     figure;imagesc(squeeze(imageOrientLPS(abs(temp)/cw,params)));axis equal tight;colormap(jet);title(titleText);drawnow;
% end


switch Trajectory
  case 'Radial'
    % region of interest
    %roi = reshape(bsxfun(@times,roi_weighting,fbp),[],Ncoils); 
    roi = reshape(bsxfun(@times,roi_mask,fbp),[],Ncoils); 
    roi = roi'*roi;

    % region of NO interest
    %roni = reshape(bsxfun(@times,1-roi_weighting,fbp),[],Ncoils);
    roni = reshape(bsxfun(@times,roni_mask,fbp),[],Ncoils); 
    roni = roni'*roni;

    % original mixer
    % [cardmixer,cardl,V] = svde(roi/roni);
    % cards = sqrt(cardl);
    % coilEnergycumsum = sqrt(cumsum(diag(cards).^2))/norm(diag(cards));

    % ROVir
    [cardmixer,Lambda] = genCCmixer(roi,roni);
    % make cardmixer orthonormal
    [u,~,v] = svd(cardmixer,'econ');
    cardmixer = u*v';
                
    fbptemp  = reshape(reshape(fbp,[],Ncoils)*cardmixer,size(fbp(:,:,:,1:Ncoils)));
    sig_roi  = squeeze(sqrt(sum(sum(sum(abs(fbptemp.*roi_mask).^2,1),2),3)));
    sig_roni = squeeze(sqrt(sum(sum(sum(abs(fbptemp.*roni_mask).^2,1),2),3)));

    % sig_roi  = squeeze(sqrt(sum(sum(sum(abs(fbptemp.*(roi_weighting)).^2,1),2),3)));
    % sig_roni = squeeze(sqrt(sum(sum(sum(abs(fbptemp.*(1-roi_weighting)).^2,1),2),3)));

    coilEnergycumsum = sqrt(cumsum(sig_roi.^2));
    coilEnergycumsum = coilEnergycumsum/coilEnergycumsum(end);
    newCoilsAuto = find(coilEnergycumsum>mincoilEnergy,1);

%     RONIEnergycumsum = sqrt(cumsum(sig_roni.^2));
%     RONIEnergycumsum = RONIEnergycumsum/RONIEnergycumsum(end);

    % sig_roi  = smooth(sig_roi);
    % [~,temp2] = kmeans(sig_roi,2);
    % newCoils = find(sig_roi<mean(temp2)/(lambertw(-2*exp(-2)) + 2),1);
    newCoils = newCoilsAuto;
    
    if isempty(newCoils) || newCoils < minNewCoils
        newCoils = minNewCoils;
    end

    fprintf('First %d virtual coils retain %.2f %% ROI signal energy', newCoilsAuto, 100*coilEnergycumsum(newCoilsAuto));

    % newCoils = max(newCoils,4);

    dataArray.mixer_all = cardmixer;
    dataArray.coilEnergycumsum = coilEnergycumsum;

    if flagCommandLine
        cw = prctile(abs(fbptemp(:)),99.5);
        titleText = sprintf('All virtual coil images. # of coils selected = %d', newCoils);
        figure;montage(squeeze(imageOrientLPS(abs(fbptemp(:,:,dispSlice,:)/cw),params)));colormap(jet);title(titleText);drawnow;
        temp = sqrt(sum(abs(fbptemp(:,:,dispSlice,1:newCoils)).^2,4));
        cw = prctile(abs(temp(:)),99.5);
        titleText = sprintf('Combined image of first %d virtual coils, %.2f %% ROI signal energy retained.', newCoils, 100*coilEnergycumsum(newCoils));
        figure;montage(squeeze(imageOrientLPS(abs(temp)/cw,params)));colormap(jet);title(titleText);drawnow;
    end

    fprintf('Coil compression: from %d to %d coils... ', Ncoils,newCoils);

    fbp = fbptemp(:,:,:,1:newCoils);
    kspaceData = reshape(permute(kspaceData,[5 2 1 3 4]),[],Ncoils);
    kspaceData = permute(reshape(kspaceData*cardmixer(1:Ncoils,1:newCoils),Necho,Nx,[],1,newCoils),[3 2 4 5 1]);
    NechoNav   = size(navData,5);
    Ntemp      = size(navData,2);
    navData    = reshape(permute(navData,[5 2 1 3 4]),[],Ncoils);
    navData    = permute(reshape(navData*cardmixer(1:Ncoils,1:newCoils),NechoNav,Ntemp,[],1,newCoils),[3 2 4 5 1]);
    dataArray.Psi = cardmixer(1:Ncoils,1:newCoils)'*dataArray.Psi*cardmixer(1:Ncoils,1:newCoils);
    Ncoils     = newCoils;

    fbpComposite = sqrt(sum(abs(fbp).^2,4));
    cw = prctile(fbpComposite(:),99.5);
    fbpComposite = fbpComposite/cw;
    if flagCommandLine
        titleText = sprintf('New coil number = %d, %.2f %% ROI signal energy retained.', Ncoils, 100*coilEnergycumsum(newCoils));
        figure;montage(squeeze(imageOrientLPS(abs(fbp(:,:,dispSlice,:)/cw),params)));title(titleText);drawnow;
    end

    params.mixer = eye(Ncoils);
  %%
  case 'Cartesian'
   
    NechoNav   = size(dataArray.navData,5);  
    kspaceData = permute(dataArray.kspaceData,[5 2 1 3 4]);
    navData    = permute(dataArray.navData,[5 2 1 3 4]);

    temp = sum(sum(roi_mask,1),3);
    ROstart = find(temp>0,1,'first');
    ROend   = find(temp>0,1,'last');
    DC = ROstart + ceil((ROend-ROstart+1)/2);
    roni_mask([1 end],ROstart:ROend) = 1;
    
    % region of interest
    roi = reshape(bsxfun(@times,roi_mask,fbp),[],Ncoils); 
    roi = roi'*roi;
    % region of NO interest
    roni = reshape(bsxfun(@times,roni_mask,fbp),[],Ncoils); 
    roni = roni'*roni;

    % ROVir
    [cardmixer,Lambda] = genCCmixer(roi,roni);
    % make cardmixer orthonormal
    [u,~,v] = svd(cardmixer,'econ');
    cardmixer = u*v';
    
    fbptemp  = reshape(reshape(fbp,[],Ncoils)*cardmixer,size(fbp(:,:,:,1:Ncoils)));
    sig_roi  = squeeze(sqrt(sum(sum(sum(abs(fbptemp.*roi_mask).^2,1),2),3)));
    sig_roni = squeeze(sqrt(sum(sum(sum(abs(fbptemp.*roni_mask).^2,1),2),3)/sum(roni_mask(:))));
    SIR = sig_roi/sqrt(sum(roi_mask(:)))./sig_roni;    
    coilEnergycumsum = sqrt(cumsum(sig_roi.^2));
    coilEnergycumsum = coilEnergycumsum/coilEnergycumsum(end);
   
    %newCoils = find(SIR>2,1,'last');
    newCoils = find(coilEnergycumsum>mincoilEnergy,1,'first');
    
    if flagCommandLine
        cw = prctile(abs(fbptemp(:)),99.5);
        titleText = sprintf('%d of coils selected', newCoils);
        figure;montage(squeeze(imageOrientLPS(abs(fbptemp(:,:,dispSlice,:)/cw),params)));colormap(jet);title(titleText);drawnow;
        figure;montage(squeeze(imageOrientLPS(angle(fbptemp(:,:,dispSlice,:))/2/pi+0.5,params)));colormap(hsv);title(titleText);drawnow;
        temp = sqrt(sum(abs(fbptemp(:,:,dispSlice,1:newCoils)).^2,4));
        cw = prctile(abs(temp(:)),99.5);
        titleText = sprintf('Combined image of %d virtual coils.', newCoils);
        figure;imagesc(squeeze(imageOrientLPS(abs(temp)/cw,params)));axis equal tight;colormap(jet);title(titleText);drawnow;
    end
    
    % ROI signal energy
    sig_roi  = squeeze(sqrt(sum(sum(sum(abs(fbp.*roi_mask).^2,1),2),3)));
    coilEnergysum = sqrt(sum(sig_roi.^2));
  
    % k-space center
    %[~, DC]     = max(sum(sum(abs(navData(:,:,:)),1),3));
  
    slwin = floor(slwin/2)*2 + 1;
    slwin_half = floor(slwin/2);
    roi_mask(:,1:slwin_half,:) = 0;
    roi_mask(:,end-slwin_half+1:end,:) = 0;
    
    mixer_all = zeros(Ncoils,newCoils,Nx);
    
    if issvROVir
        winidx = DC-slwin_half:DC+slwin_half;
        
        % region of interest
        roi = reshape(bsxfun(@times,roi_mask(:,winidx,:),fbp(:,winidx,:,:)),[],Ncoils); 
        roi = roi'*roi;

        % region of NO interest
        roni = reshape(bsxfun(@times,roni_mask(:,winidx,:),fbp(:,winidx,:,:)),[],Ncoils); 
        roni = roni'*roni;

        % compression matrix
        [cardmixer,Lambda] = genCCmixer(roi,roni);
        %cardmixer = cardmixer(:,1:newCoils);
        % make cardmixer orthonormal
        [u,~,v] = svd(cardmixer,'econ');
        cardmixer = u*v';
                
        cardmixer = cardmixer(:,1:newCoils);
        
        kspaceData(:,DC,:,1,(1:newCoils)) = reshape(reshape(kspaceData(:,DC,:,:,:),[],Ncoils)*cardmixer,Necho,1,[],1,newCoils);  
        navData(:,DC,:,1,(1:newCoils))    = reshape(reshape(navData(:,DC,:,:,:),[],Ncoils)*cardmixer,NechoNav,1,[],1,newCoils); 
        Psi_new(:,:,DC)                   = cardmixer'*Psi_orig*cardmixer;  % update noise matrix
        fbp(:,DC,:,(1:newCoils))          = reshape(reshape(fbp(:,DC,:,:),[],Ncoils)*cardmixer,size(fbp(:,1,:,1:newCoils)));
        
        mixer_all(:,:,DC) = cardmixer;
        cardmixer_orig = cardmixer;

        % Forward from center
        cardmixer_last = cardmixer_orig;
        for j = DC+1:Nx
            winidx = j-slwin_half:j+slwin_half;
            
            % region of interest
            roi = reshape(bsxfun(@times,roi_mask(:,winidx,:),fbp(:,winidx,:,:)),[],Ncoils); 
            roi = roi'*roi;

            if sum(abs(roi(:))) > 0
                % region of NO interest
                roni = reshape(bsxfun(@times,roni_mask(:,winidx,:),fbp(:,winidx,:,:)),[],Ncoils); 
                roni = roni'*roni;

                % compression matrix
                [cardmixer,Lambda] = genCCmixer(roi,roni);
                %cardmixer = cardmixer(:,1:newCoils);
                % make cardmixer orthonormal
                [u,~,v] = svd(cardmixer,'econ');
                cardmixer = u*v';
                
                cardmixer = cardmixer(:,1:newCoils);
                
                % phase alignment
                [U,~,V] = svde(cardmixer_last'*cardmixer);
                cardmixer   = cardmixer*V*U';

                kspaceData(:,j,:,1,1:newCoils) = reshape(reshape(kspaceData(:,j,:,:,:),[],Ncoils)*cardmixer,Necho,1,[],1,newCoils); 
                navData(:,j,:,1,1:newCoils)    = reshape(reshape(navData(:,j,:,:,:),[],Ncoils)*cardmixer,NechoNav,1,[],1,newCoils); 
                Psi_new(:,:,j)                 = cardmixer'*Psi_orig*cardmixer;     % update noise matrix
                fbp(:,j,:,1:newCoils)          = reshape(reshape(fbp(:,j,:,:),[],Ncoils)*cardmixer,size(fbp(:,1,:,1:newCoils)));

                mixer_all(:,:,j) = cardmixer;
                cardmixer_last   = cardmixer;
            else
                ROend = j-1;
                kspaceData(:,j:end,:,1,:) = 0; 
                navData(:,j:end,:,1,:)    = 0; 
                fbp(:,j:end,:,:)          = 0;
                mixer_all(:,:,j:end)      = 0;
                break;
            end
        end

        % Backward from center
        cardmixer_last = cardmixer_orig;
        for j = DC-1:-1:1
            winidx = j-slwin_half:j+slwin_half;
            
            % region of interest
            roi = reshape(bsxfun(@times,roi_mask(:,winidx,:),fbp(:,winidx,:,:)),[],Ncoils); 
            roi = roi'*roi;

            if sum(abs(roi(:))) > 0
                % region of NO interest
                roni = reshape(bsxfun(@times,roni_mask(:,winidx,:),fbp(:,winidx,:,:)),[],Ncoils); 
                roni = roni'*roni;

                % compression matrix
                [cardmixer,Lambda] = genCCmixer(roi,roni);
                %cardmixer = cardmixer(:,1:newCoils);
                % make cardmixer orthonormal
                [u,~,v] = svd(cardmixer,'econ');
                cardmixer = u*v';
                
                cardmixer = cardmixer(:,1:newCoils);
                % phase alignment
                %[U,~,V] = svd(cardmixer_last'*cardmixer,'econ');
                [U,~,V] = svde(cardmixer_last'*cardmixer);
                cardmixer   = cardmixer*V*U';

                
                kspaceData(:,j,:,1,1:newCoils) = reshape(reshape(kspaceData(:,j,:,:,:),[],Ncoils)*cardmixer,Necho,1,[],1,newCoils); 
                navData(:,j,:,1,1:newCoils)    = reshape(reshape(navData(:,j,:,:,:),[],Ncoils)*cardmixer,NechoNav,1,[],1,newCoils); 
                Psi_new(:,:,j)                 = cardmixer'*Psi_orig*cardmixer; %update noise matrix
                fbp(:,j,:,1:newCoils)          = reshape(reshape(fbp(:,j,:,:),[],Ncoils)*cardmixer,size(fbp(:,1,:,1:newCoils)));

                mixer_all(:,:,j) = cardmixer;
                cardmixer_last = cardmixer;
            else
                ROstart = j + 1;
                kspaceData(:,1:j,:,1,:) = 0; 
                navData(:,1:j,:,1,:)    = 0; 
                Psi_new(:,:,1:j)        = [];
                fbp(:,1:j,:,:)          = 0;
                mixer_all(:,:,1:j)      = 0;
                break;
            end
        end
        
        kspaceData = kspaceData(:,ROstart:ROend,:,:,:);
        navData    = navData(:,ROstart:ROend,:,:,:);
        fbp        = fbp(:,ROstart:ROend,:,:);
        fbp_window = fbp_window(:,ROstart:ROend,:);
        roi_mask   = roi_mask(:,ROstart:ROend,:);
        roni_mask  = roni_mask(:,ROstart:ROend,:);
    else
        % region of interest
        roi = reshape(bsxfun(@times,roi_mask,fbp),[],Ncoils); 
        roi = roi'*roi;
        % region of NO interest
        roni = reshape(bsxfun(@times,roni_mask,fbp),[],Ncoils); 
        roni = roni'*roni;

        % compression matrix
        [cardmixer,Lambda] = genCCmixer(roi,roni);
        %cardmixer = cardmixer(:,1:newCoils);
        % make cardmixer orthonormal
        [u,~,v] = svd(cardmixer,'econ');
        cardmixer = u*v';
        cardmixer = cardmixer(:,1:newCoils);
        
        kspaceData = reshape(reshape(kspaceData,[],Ncoils)*cardmixer,Necho,Nx,[],1,newCoils);  
        navData    = reshape(reshape(navData,[],Ncoils)*cardmixer,NechoNav,Nx,[],1,newCoils); 
        Psi_new    = cardmixer'*Psi_orig*cardmixer;  % update noise matrix
        fbp        = reshape(reshape(fbp,[],Ncoils)*cardmixer,size(fbp(:,:,:,1:newCoils)));
        mixer_all  = cardmixer;
        
        kspaceData = kspaceData(:,ROstart:ROend,:,:,:);
        navData    = navData(:,ROstart:ROend,:,:,:);
        fbp        = fbp(:,ROstart:ROend,:,:);
        fbp_window = fbp_window(:,ROstart:ROend,:);
        roi_mask   = roi_mask(:,ROstart:ROend,:);
        roni_mask  = roni_mask(:,ROstart:ROend,:);
    end
    
    fbptemp = fbp(:,:,:,1:newCoils);
    
    % new ROI signal energy
    sig_roi  = squeeze(sqrt(sum(sum(sum(abs(fbptemp.*roi_mask).^2,1),2),3)));
    sig_roni = squeeze(sqrt(sum(sum(sum(abs(fbptemp.*roni_mask).^2,1),2),3)/sum(roni_mask(:))));
    SIR = sig_roi/sqrt(sum(roi_mask(:)))./sig_roni;
    
    coilEnergysumNew = sqrt(sum(sig_roi(1:newCoils).^2));
    coilEnergycumsum = sqrt(cumsum(sig_roi.^2));
    coilEnergycumsum = coilEnergycumsum/coilEnergycumsum(end);

    RONIEnergycumsum = sqrt(cumsum(sig_roni.^2));
    RONIEnergycumsum = RONIEnergycumsum/RONIEnergycumsum(end);

    % sig_roi  = smooth(sig_roi);
    % [~,temp2] = kmeans(sig_roi,2);
    % newCoils = find(sig_roi<mean(temp2)/(lambertw(-2*exp(-2)) + 2),1);
%     if isempty(newCoils) || newCoils < minNewCoils
%         newCoils = minNewCoils;
%     end
    % if newCoils < minNewCoils
    %     newCoils = minNewCoils;
    % end

    fprintf('First %d virtual coils retain %.2f %% ROI signal energy\n', newCoils, 100*coilEnergysumNew/coilEnergysum);

    dataArray.mixer_all        = mixer_all;
    dataArray.coilEnergycumsum = coilEnergycumsum;
    
    if flagCommandLine
        cw = prctile(abs(fbptemp(:)),99.5);
        titleText = sprintf('%d of coils selected', newCoils);
        figure;montage(squeeze(imageOrientLPS(abs(fbptemp(:,:,dispSlice,:)/cw),params)));colormap(jet);title(titleText);drawnow;
        figure;montage(squeeze(imageOrientLPS(angle(fbptemp(:,:,dispSlice,:))/2/pi+0.5,params)));colormap(hsv);title(titleText);drawnow;
        temp = sqrt(sum(abs(fbptemp(:,:,dispSlice,1:newCoils)).^2,4));
        cw = prctile(abs(temp(:)),99.5);
        titleText = sprintf('Combined image of %d virtual coils.', newCoils);
        figure;imagesc(squeeze(imageOrientLPS(abs(temp)/cw,params)));axis equal tight;colormap(jet);title(titleText);drawnow;
    end

    fprintf('Coil compression: from %d to %d coils... \n', Ncoils,newCoils);
    
    % Replace data with compressed version
    kspaceData = permute(kspaceData(:,:,:,:,1:newCoils),[3 2 4 5 1]);
    navData    = permute(navData(:,:,:,:,1:newCoils),[3 2 4 5 1]);
        
    fbp = fbp(:,:,:,1:newCoils);
    fbpComposite = sqrt(sum(abs(fbp).^2,4));
    cw  = prctile(fbpComposite(:),99.5);
    fbpComposite = fbpComposite/cw;
%     if flagCommandLine
%         titleText = sprintf('New coil number = %d, %.2f %% ROI signal energy retained.', Ncoils, 100*coilEnergycumsum(newCoils));
%         figure;montage(squeeze(imageOrientLPS(abs(fbp(:,:,dispSlice,:)/cw),params)));title(titleText);drawnow;
%     end

    dataArray.ROstart = ROstart;
    dataArray.ROend   = ROend;
    dataArray.Psi = mean(Psi_new(1:newCoils,1:newCoils),3);    % average noise matrix
    Ncoils        = newCoils;
    
    newNx = size(navData,2);
    params.Nx     = newNx;
    params.Nxdisp = newNx;
    params.N     = newNx;
    params.Norig = newNx;
    reconOptions.st.Nd(2)    = newNx;
    reconOptions.st.Ndisp(2) = newNx;
    st.Nd(2)    = newNx;
    setupFunctions;     % redo setup
end

params.Ncoils = Ncoils;
params.mixer  = eye(Ncoils);

dataArray.navData    = navData;
dataArray.kspaceData = kspaceData;
dataArray.flagIsCompressed = true;

dataArray.fbp = fbp;
dataArray.fbp_window   = fbp_window;
dataArray.fbpComposite = fbpComposite;

fprintf(' done.\n')

%%
function [mixer, Lambda] = genCCmixer(roi,roni)
    %[eigVector, eigValues] = eig(roni\roi);
    [eigVector, eigValues] = eig(roi,roni);
    Lambda = sqrt(abs(diag(eigValues)));
    [~, idx]  = sort(Lambda,'descend');
    mixer = eigVector(:,idx);
    Lambda = Lambda(idx);
  
    