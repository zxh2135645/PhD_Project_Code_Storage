%% Select IR-FLASH recon

if ispc
    sep = '\';
else
    sep = '/';
end

folders=dir;
found=[];
for j=1:numel(folders);
    if folders(j).isdir
        file=dir(sprintf('%s%sAC_recon.mat',folders(j).name,sep));
        if numel(file) > 0;
        found{end+1}=folders(j).name;
        end
    end
end

folder=listdlg('ListString',found,'SelectionMode','single','ListSize',[500 150]);
load(sprintf('%s%sAC_recon.mat',found{folder},sep),'U','Phi','L','Gr',...
    'sizes','Ny','Nx','N','Nseg','ovs','params','dispim','vec','cw')

V=reshape(Gr\Phi(:,:),[L sizes(2:end)]);
U=dispim(reshape(U,Ny,Nx,L));
[Ny, Nx, ~]=size(U);
U=reshape(U,[],L)/cw;

%%

ti=1;
tc=1;
tr=1;

f = figure;
h = imshow(reshape(abs(U*V(:,ti,tc,tr)),Ny,Nx));

jSlider1 = javax.swing.JSlider;
javacomponent(jSlider1,[0,40,200,20]);
jSlider1.setValue(ceil(100/sizes(2)));

jSlider2 = javax.swing.JSlider;
javacomponent(jSlider2,[0,20,200,20]);
jSlider2.setValue(ceil(100/sizes(3)));

jSlider3 = javax.swing.JSlider;
javacomponent(jSlider3,[0,0,200,20]);
jSlider3.setValue(ceil(100/sizes(4)));

indround=@(val,dim)max(1,round(val/100*sizes(dim)));

hjSlider1 = handle(jSlider1, 'CallbackProperties');
hjSlider1.StateChangedCallback = @(hjSlider,eventData)...
  imshow(reshape(abs(U*V(:,indround(hjSlider1.getValue,2),indround(jSlider2.getValue,3),indround(jSlider3.getValue,4))),Ny,Nx));

hjSlider2 = handle(jSlider2, 'CallbackProperties');
hjSlider2.StateChangedCallback = @(hjSlider,eventData)...
  imshow(reshape(abs(U*V(:,indround(jSlider1.getValue,2),indround(hjSlider2.getValue,3),indround(jSlider3.getValue,4))),Ny,Nx));

hjSlider3 = handle(jSlider3, 'CallbackProperties');
hjSlider3.StateChangedCallback = @(hjSlider,eventData)...
  imshow(reshape(abs(U*V(:,indround(jSlider1.getValue,2),indround(jSlider2.getValue,3),indround(hjSlider3.getValue,4))),Ny,Nx));