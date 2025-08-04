function dataArray = flipRidx(params,dataArray,flagCommandLine)

if nargin < 3
    flagCommandLine = true;
end

dataArray.Ridx = max(dataArray.Ridx) - dataArray.Ridx + 1;
dataArray.bestresResp = flip(dataArray.bestresResp);
if flagCommandLine
    h = findall(groot,'Type','figure','Name','Binning Results');
    if isempty(h)
        figure('Name','Binning Results','units','normalized','OuterPosition',[0.1 0.4 0.8 0.5]);
    else
        figure(h);
    end
    subplot(2,2,3),plot(dataArray.Ridx,'.-');axis([-inf inf 0 max(dataArray.Ridx)+1]);title(sprintf('Ridx, mean cycle %.2f sec',dataArray.meanRPeriod));
    subplot(2,2,4),plot(dataArray.Segidx(:), dataArray.Ridx(:),'.');axis([0 max(dataArray.Segidx(:)) 0 max(dataArray.Ridx)+1]);title('Ridx / Segment Index')
end

if isfield(dataArray,'binsRespMean')
    dataArray.binsRespMean = flip(dataArray.binsRespMean,4);
    if flagCommandLine
        implayZoom(imageOrientLPS(dataArray.binsRespMean(:,:,1,:),params),2);
    end
end

if isfield(dataArray,'binsResp')
    rbins = length(dataArray.binsResp);
    temp = dataArray.binsResp;
    for j = 1:rbins
        dataArray.binsResp{j} = temp{rbins-j+1};
    end
end

