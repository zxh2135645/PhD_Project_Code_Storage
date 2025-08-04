function selectedFrame = montageSliceViewer(img,titleString)

if nargin < 2
    titleString = 'Select slice';
end

thumbnailSize = "auto";
cmap          = [];
montageSize   = [];
parent        = [];
borderSize    = [0 0];
backgroundColor = [];
interpolation = 'nearest';
indices = [];
waitbarEnabled = true;

[Ny,Nx,Nz,Nv] = size(img);

for v = 1:size(img,4)
[bigImage(:,:,v), cmap] = images.internal.createMontage(reshape(img(:,:,:,v),Ny,Nx,1,Nz), thumbnailSize,...
    montageSize, borderSize, backgroundColor, indices, cmap, ...
    waitbarEnabled);
end

try
    f = figure(201); 
    s = sliceViewer(bigImage,'Parent',f);  
    hXYAxe = getAxesHandle(s);
%     title(hXZAxes,titleString,'FontSize',11,'FontWeight','bold','Color',[1.0 0.3250 0.0980]);
    f.Units = 'normalized';
    f.Position = [0.1 0.2 0.9 0.9];
    titleText = uicontrol('Parent',f,'Style','text','String',titleString,'Units','normalized','Position',[0.8 0.01 0.15 0.05],'FontSize',16,'FontWeight','bold','Visible','on');
    ButtonH = uicontrol('Parent',f,'Style','togglebutton','String','Confirm','FontSize',16,'FontWeight','bold','Units','normalized','Position',[0.8 0.45 0.15 0.05],'Visible','on');
    waitfor(ButtonH,'Value');
    selectedFrame = s.SliceNumber;
    close(f);
catch errormsg
    fprintf(2, '%s\n', errormsg.message);
    fprintf(2, 'Frame selection failed.\n');
end
