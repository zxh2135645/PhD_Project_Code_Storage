% XZ
Phi=reshape(Phi,[L sizes(2:end)]);
temp=Gr\reshape(Phi(:,:,1,1,:),L,[]);

temp=(reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp,Ny,Nx,[],params.NEco));
cw=max(vec(abs(temp)));

h = implay(abs(temp(:,:,:,1))/cw);
set(h.Parent,'Name','old_basal_echo1');

%% Pick remote ROI for SNR analysis
slc = 192;
figure('Position', [100 100 1600 1000]); 
imagesc(abs(temp(:,:,slc,1))); axis image;
roi = drawpolygon;
roi_pos = roi.Position;
%%
x = roi_pos(:,1);
y = roi_pos(:,2);
BW = poly2mask(x,y,192,192);
%figure(); imagesc(BW)
sig_mean = mean(nonzeros(BW.*abs(temp(:,:,slc,1))));
sig_sd = std(nonzeros(BW.*abs(temp(:,:,slc,1))));

SNR = sig_mean/sig_sd;
%% resphape it into chronologic order
tempp = reshape(permute(temp, [1 2 4 3]), Ny, Nx, []);
h = implay(abs(tempp)/cw);
set(h.Parent,'Name','old_basal_echo1');