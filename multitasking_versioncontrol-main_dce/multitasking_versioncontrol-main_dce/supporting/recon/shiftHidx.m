function dataArray = shiftHidx(params,reconOptions,dataArray,temporalBasis,spatialCoeff)

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(dataArray);
extractVarFromStruct(temporalBasis);
extractVarFromStruct(spatialCoeff);

dataArray.Hidx = mod(Hidx - diastoleIdx, max(Hidx)) + 1;
dataArray.diastoleIdx = 1;
dataArray.systoleIdx  = floor(max(Hidx)/2) + 1;

dataArray.binsCard = circshift(dataArray.binsCard,1-diastoleIdx,4);

dt = lEchoSpacing*SGBlock;

if flagCommandLine && cbins > 1
    h = findall(groot,'Type','figure','Name','Binning Results');
    if isempty(h)
        figure('Name','Binning Results','units','normalized','OuterPosition',[0.1 0.4 0.8 0.5]);
    else
        figure(h);
    end
    subplot(2,2,3),plot(dt:dt:dt*numel(Ridx),Ridx,'.-');axis([-inf inf 0 max(Ridx)+1]);title(sprintf('Ridx, mean cycle %.2f sec',dataArray.meanRPeriod));
    subplot(2,2,4),plot(Segidx(:), dataArray.Ridx(:),'.');axis([0 max(Segidx(:)) 0 max(Ridx)+1]);title('Ridx / Segment Index');
    subplot(2,2,1),plot(dt:dt:dt*numel(Hidx),Hidx,'.-');axis([-inf inf 0 max(Hidx)+1]);title(sprintf('Hidx, mean BPM %d',dataArray.meanHBPM));
    subplot(2,2,2);
    for rphase = 1:rbins
        plot(Segidx((Ridx==rphase)), Hidx((Ridx==rphase)),'.');hold on;
    end
    hold off; axis([0 linesPerShot*moduleLength 0 max(Hidx)+1]);title('Hidx / Segment Index');
    
    temp = imageOrientLPS(dataArray.binsCard(:,:,1,:,1),params);
    implayZoom(squeeze(abs(temp)),max(Hidx));
end