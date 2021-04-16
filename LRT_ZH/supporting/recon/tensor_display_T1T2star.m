% XZ
Phi=reshape(Phi,[L sizes(2:end)]);
temp=Gr\reshape(Phi(:,:,1,1,:),L,[]);

temp=(reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp,Ny,Nx,[],params.NEco));
cw=max(vec(abs(temp)));

h = implay(abs(temp(:,:,:,1))/cw);
set(h.Parent,'Name','old_basal_echo1');

%% resphape it into chronologic order
tempp = reshape(permute(temp, [1 2 4 3]), Ny, Nx, []);
h = implay(abs(tempp)/cw);
set(h.Parent,'Name','old_basal_echo1');