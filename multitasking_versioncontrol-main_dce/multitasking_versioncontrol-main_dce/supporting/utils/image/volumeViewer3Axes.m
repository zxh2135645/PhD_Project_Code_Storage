function selectedXYZ = volumeViewer3Axes(img,titleString)

if nargin < 2
    titleString = 'Select voxel';
end

try
    f = figure(201); 
    s = orthosliceViewer(img,'Parent',f);  
%     [hXYAxes, hYZAxes, hXZAxes] = getAxesHandles(s);
%     title(hXZAxes,titleText,'FontSize',11,'FontWeight','bold','Color',[1.0 0.3250 0.0980]);
    f.Units = 'normalized';
    f.Position = [0.1 0.2 0.8 0.7];
    titleText = uicontrol('Parent',f,'Style','text','String',titleString,'FontSize',16,'FontWeight','bold','Units','normalized','Position',[0.75 0.55 0.2 0.1],'Visible','on');
    tempValue = img(s.SliceNumbers(2),s.SliceNumbers(1),s.SliceNumbers(3));
    valueString = sprintf('Selected voxel value = %.2f',tempValue);
    valueText = uicontrol('Parent',f,'Style','text','String',valueString,'Units','normalized','Position',[0.75 0.3 0.2 0.1],'Visible','on');
    addlistener(s,'CrosshairMoving',@(src,evt)allevents(src,evt,s,valueText,img));
    addlistener(s,'CrosshairMoved',@(src,evt)allevents(src,evt,s,valueText,img));
    ButtonH = uicontrol('Parent',f,'Style','togglebutton','String','Confirm voxel selection','FontSize',11,'FontWeight','bold','Units','normalized','Position',[0.75 0.45 0.2 0.1],'Visible','on');
    waitfor(ButtonH,'Value');
    selectedXYZ(1) = s.SliceNumbers(2);
    selectedXYZ(2) = s.SliceNumbers(1);
    selectedXYZ(3) = s.SliceNumbers(3);
    close(f);
catch errormsg
    fprintf(2, '%s\n', errormsg.message);
    fprintf(2, 'Voxel selection failed.\n');
end

%fprintf('image value = %g at (%d, %d, %d).\n',img(selectedXYZ(1),selectedXYZ(2),selectedXYZ(3)),selectedXYZ(1),selectedXYZ(2),selectedXYZ(3));

function allevents(src,evt,s,valueText,img)
tempValue = img(s.SliceNumbers(2),s.SliceNumbers(1),s.SliceNumbers(3));
valueString = sprintf('Selected voxel value = %.2f',tempValue);
valueText.String = valueString;

% [Ny,Nx,Nz] = size(img);
% currentX = floor(Nx/2);
% currentY = floor(Ny/2);
% currentZ = floor(Nz/2);
% 
% volumeViewer3Axes = figure;
% % x-y slice
% ax1 = axes('Parent',volumeViewer3Axes,'position',[0.07 0.55  0.4 0.4]);
% imshow(img(:,:,currentZ), 'Parent', ax1);
% 
% % z-x slice
% ax2 = axes('Parent',volumeViewer3Axes,'position',[0.52 0.55  0.4 0.4]);
% imshow(squeeze(img(currentY,:,:)), 'Parent', ax2);
% 
% % z-y slice
% ax3 = axes('Parent',volumeViewer3Axes,'position',[0.07 0.07  0.4 0.4]);
% imshow(squeeze(img(:,currentX,:)), 'Parent', ax3);
% 
% roiPOC1 = drawpoint(ax1);    %Use Mouse To Select a point ROI
% PosPOC1 = round(get(roiPOC1,'Position'));
% currentX = PosPOC1(1);
% currentY = PosPOC1(2);
% 
% roiPOC2 = drawpoint(ax2);    %Use Mouse To Select a point ROI
% PosPOC2 = round(get(roiPOC2,'Position'));
% 
% roiPOC3 = drawpoint(ax3);    %Use Mouse To Select a point ROI
% PosPOC3 = round(get(roiPOC3,'Position'));
% 
% img = repmat(temp1,[1 1 32]);
% currentXYZ = [1 1 1];
% addlistener(roiPOC1,'MovingROI',@(src,evt)allevents1(src,evt,img,currentXYZ));
% addlistener(roiPOC1,'ROIMoved', @(src,evt)allevents1(src,evt,img,currentXYZ));
% 
% end
% 
% function allevents1(src,evt,img,currentXYZ)
% evname = evt.EventName;
% round(evt.CurrentPosition)
% size(round(evt.CurrentPosition))
% newXYZ = ([round(evt.CurrentPosition),1]);
% newXYZ = newXYZ(:)'
% size(img)
% temp1 = img(newXYZ(2),newXYZ(1),newXYZ(3))
%     switch(evname)
%         case{'MovingROI'}
%             disp(['ROI moving Previous Position: ' mat2str(evt.PreviousPosition)]);
%             disp(['ROI moving Current Position: ' mat2str(evt.CurrentPosition)]);
%         case{'ROIMoved'}
%             disp(['ROI moved Previous Position: ' mat2str(evt.PreviousPosition)]);
%             disp(['ROI moved Current Position: ' mat2str(evt.CurrentPosition)]);
%     end
% end

