function htemp = implayZoom(image,fps,zoom,title)

if nargin < 3 || isempty(zoom)
    zoom = 2;
end
if nargin < 2
    fps = 20;
end

if size(image,3) ~= 3
    image = image(:,:,:);
end

htemp = implay(abs(image),fps);
%htemp.Visual.Axes.Position(3:4) = htemp.Visual.Axes.Position(3:4)*zoom;
htemp.Parent.Position(3:4) = htemp.Visual.Axes.Position(3:4)*zoom;
htemp.DataSource.Controls.Repeat = 1;
htemp.Parent.findobj('TooltipString','Maintain fit to window').ClickedCallback();
htemp.Parent.findobj('Tag','uimgr.spctoggletool_Repeat').State = 'on';
if nargin > 3
    set(htemp.Parent, 'Name', title);
end
