% Sliders for 2 spatial dimensions and 4 time dimensions. Fast(ish), but high memory usage.

%% Select recon MAT file

if ispc
  sep = '\';
else
  sep = '/';
end

folders=dir;
found=[];
for j=1:numel(folders);
  if folders(j).isdir
    file=dir(sprintf('%s%smultitasking.mat',folders(j).name,sep));
    if numel(file) > 0;
      found{end+1}=folders(j).name;
    end
  end
end

folder=listdlg('ListString',found,'SelectionMode','single','ListSize',[500 150]);
load(sprintf('%s%smultitasking.mat',found{folder},sep),'U','Phi','L','Gr',...
  'sizes','Ny','Nx','N','Nseg','ovs','params','dispim','vec','cw')

if sizes(2)>100
%   inds=sizes(2):-floor(sizes(2)/100):1;
%   inds=inds(100:-1:1);
  inds=floor(sizes(2)/100):floor(sizes(2)/100):sizes(2);
  Phi=reshape(Phi,[L sizes(2:end)]);
  Phi=Phi(:,inds,:,:,:);
  sizes(2)=numel(inds);
end

UV=dispim(reshape(U,Ny,Nx,L));
[rows, cols, ~]=size(UV);
UV=reshape(UV,[],L);

UV=uint8(255*abs(reshape(UV*(Gr\Phi(:,:)),[rows,cols,sizes(2:end)]))/cw);

% UV=reshape(UV*(Gr\Phi(:,:)),[rows,cols,sizes(2:end)]);
% UV=real(UV.*exp(-1i*angle(repmat(UV(:,:,end,:,:),[1 1 sizes(2) 1 1]))));
% UV=UV/(2*cw)+.5;
% UV=uint8(255*UV);

%%

ti=1;
tc=1;
tr=1;
te=1;

f = figure;
h = imshow(UV(:,:,ti,tc,tr,te));

width=get(f, 'Position');
width=width(3);

jSlider1 = javax.swing.JSlider;
javacomponent(jSlider1,[width/4-50,20,100,20]);
jSlider1.setValue(ceil(100/sizes(2)));

jSlider2 = javax.swing.JSlider;
javacomponent(jSlider2,[3*width/4-50,20,100,20]);
jSlider2.setValue(ceil(100/sizes(3)));

jSlider3 = javax.swing.JSlider;
javacomponent(jSlider3,[width/4-50,0,100,20]);
jSlider3.setValue(ceil(100/sizes(4)));

jSlider4 = javax.swing.JSlider;
javacomponent(jSlider4,[3*width/4-50,0,100,20]);
jSlider4.setValue(ceil(100/sizes(5)));

indround=@(val,dim)max(1,round(val/100*sizes(dim)));

hjSlider1 = handle(jSlider1, 'CallbackProperties');
hjSlider1.StateChangedCallback = @(hjSlider,eventData)...
  imshow(UV(:,:,indround(hjSlider1.getValue,2),indround(jSlider2.getValue,3),indround(jSlider3.getValue,4),indround(jSlider4.getValue,5)));

hjSlider2 = handle(jSlider2, 'CallbackProperties');
hjSlider2.StateChangedCallback = @(hjSlider,eventData)...
  imshow(UV(:,:,indround(jSlider1.getValue,2),indround(hjSlider2.getValue,3),indround(jSlider3.getValue,4),indround(jSlider4.getValue,5)));

hjSlider3 = handle(jSlider3, 'CallbackProperties');
hjSlider3.StateChangedCallback = @(hjSlider,eventData)...
  imshow(UV(:,:,indround(jSlider1.getValue,2),indround(jSlider2.getValue,3),indround(hjSlider3.getValue,4),indround(jSlider4.getValue,5)));

hjSlider4 = handle(jSlider4, 'CallbackProperties');
hjSlider4.StateChangedCallback = @(hjSlider,eventData)...
  imshow(UV(:,:,indround(jSlider1.getValue,2),indround(jSlider2.getValue,3),indround(jSlider3.getValue,4),indround(hjSlider4.getValue,5)));