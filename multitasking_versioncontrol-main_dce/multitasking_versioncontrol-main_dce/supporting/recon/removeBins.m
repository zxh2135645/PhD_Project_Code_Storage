function dataArray = removeBins(binsToRemove,params,dataArray,flagCommandLine)

if nargin < 4
    flagCommandLine = true;
end

extractVarFromStruct(params);
extractVarFromStruct(dataArray);

rbins = max(Ridx);
if ~isempty(binsToRemove)
    binsToRemove = sort(binsToRemove,'descend');
    binsToRemove(binsToRemove>rbins) = [];
    binsToRemove(binsToRemove<1) = [];
    for n = 1:numel(binsToRemove)
        dataArray.Ridx(dataArray.Ridx == binsToRemove(n)) = 0;
        dataArray.Ridx(dataArray.Ridx > binsToRemove(n)) = dataArray.Ridx(dataArray.Ridx > binsToRemove(n)) - 1;
        dataArray.binsRespMean(:,:,:,n) = [];
        rbins = rbins - 1;
        if isfield(dataArray,'bins_slice')
            dataArray.bins_slice(:,:,n) = [];
        end
    end
end

% for 3D volume acquired in transverse orientation, also display resp bin images in coronal view
if isfield(dataArray,'bins_slice') && flagCommandLine && rbins > 1
    implayZoom(dataArray.bins_slice,2);
end

dt = lEchoSpacing*SGBlock;

if flagCommandLine && rbins > 1
    h = findall(groot,'Type','figure','Name','Binning Results');
    if isempty(h)
        figure('Name','Binning Results','units','normalized','OuterPosition',[0.1 0.4 0.8 0.5]);
    else
        figure(h);
    end
    subplot(2,2,3),plot(dt:dt:dt*numel(Ridx),Ridx,'.-');axis([-inf inf 0 rbins+1+numel(binsToRemove)]);title(sprintf('Ridx, mean cycle %.2f sec',dataArray.meanRPeriod));
    subplot(2,2,4),plot(Segidx(:), Ridx(:),'.');axis([0 linesPerShot*moduleLength 0 rbins+1]);title('Ridx / Segment Index');
    
    implayZoom(imageOrientLPS(dataArray.binsRespMean(:,:,1,:),params),2);
end
