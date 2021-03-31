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

UV=dispim(reshape(U,Ny,Nx,L));
[rows, cols, ~]=size(UV);
UV=reshape(UV,[],L);

Phi_sm=reshape(Phi(:,:,1,1,:),L,[]); %change cardiac/resp phase if desired
UV=UV.*repmat(exp(-1i*angle(UV*(Gr\Phi_sm(:,end)))),[1 L]);

%%
I=abs(reshape(UV*(Gr\Phi_sm(:,end)),rows,cols));
figure, imshow(I/cw);
title('Select voxel to suppress')
h = impoint(gca,[]);
vox_null = round(getPosition(h));
vox_null=sub2ind([rows,cols],vox_null(2),vox_null(1));
sig_null=UV(vox_null,:)*(Gr\Phi_sm);
[~,null_time]=min(sig_null);


I=abs(reshape(UV*(Gr\Phi_sm(:,null_time)),rows,cols));
imshow(I/cw);
title('Select voxel to emphasize')
h = impoint(gca,[]);
vox_emph = round(getPosition(h));
vox_emph=sub2ind([rows,cols],vox_emph(2),vox_emph(1));
sig_emph=UV(vox_emph,:)*(Gr\reshape(Phi(:,:,1,1,:),L,[]));

% best_time=null_time;
% I=abs(reshape(UV*(Gr\Phi_sm(:,best_time)),rows,cols));
% imshow((I-abs(sig_null(best_time)))/(abs(sig_emph(best_time))-abs(sig_null(best_time))));

[~,best_time]=max(abs(real(sig_emph)).*(sign(real(sig_emph))~=sign(real(sig_null))));
I=real(reshape(UV*(Gr\Phi_sm(:,best_time)),rows,cols))*(-sign(real(sig_null(best_time))));
figure,imshow(I/min(cw,max(I(:))))
