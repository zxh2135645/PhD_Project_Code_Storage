function [params, reconOptions, dataArray] = cropRO(params, reconOptions, dataArray)

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);

if isCartesian
    h = figure;imagesc(abs(dataArray.fbpComposite(:,:,dispSlice,1)));axis equal tight;colormap('gray');title('Draw ROI')
    roiResp = imrect;
    roiPosition = roiResp.getPosition();
    roiPosition = round(roiPosition);
    close(h);

    Nxshift = roiPosition(1);
    ROINx   = roiPosition(3);
    ROINx   = ROINx + mod(ROINx,2);
    newXIdx = Nxshift+(1:ROINx);

    ROItitle = ['New readout ROI: ' num2str(ROINx)];
    h = figure;imagesc(abs(dataArray.fbpComposite(:,newXIdx,dispSlice,1)));axis equal tight;colormap('gray');title(ROItitle);

    answer = questdlg(  'Crop to the new ROI?', ...
                        'Confirm new ROI', ...
                        'Crop','Cancel','Crop');
    close(h);

    if strcmp(answer,'Crop')
        params.Nx       = ROINx;
        params.Nkx      = ROINx;
        params.Nxdisp   = ROINx;
        params.N        = ROINx;
        params.Norig    = ROINx;
        params.DC       = floor(ROINx/2) + 1; 
        params.DC_kx    = params.DC;
        reconOptions.st.Nd(2)    = ROINx;
        reconOptions.st.Ndisp(2) = ROINx;
        dataArray.kspaceData   = dataArray.kspaceData(:,newXIdx,:,:,:);
        dataArray.navData      = dataArray.navData(:,newXIdx,:,:,:); 
        dataArray.fbpComposite = dataArray.fbpComposite(:,newXIdx,:,:,:);
        dataArray.SEs          = dataArray.SEs(:,newXIdx,:,:,:);
        fprintf('Data cropped in readout dimension. New image size: %d x %d x %d.\n',Ny,params.Nx,Nz);
    else
        fprintf('Cancelled');
    end
else
    fprintf('Cropping is not supported for non-Cartesian data!\n');
end
