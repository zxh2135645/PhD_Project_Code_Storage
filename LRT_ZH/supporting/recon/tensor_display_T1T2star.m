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

%% 5D 
Phi=reshape(Phi,[L sizes(2:end)]);
temp=Gr\reshape(Phi(:,:,1,1,:,4),L,[]);



temp=(reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp,Ny,Nx,[],params.NEco));
cw=max(vec(abs(temp(:,:,:,1))));

h = implay(abs(temp(:,:,:,8))/cw);
set(h.Parent,'Name','old_basal_echo1');
%%
img = abs(temp(:,:,96,1));
figure, imagesc(img);
roi = impoly;
roi_mask = createMask(roi);

%% 
t1_array = zeros(params.lSegments, size(Phi,6));

for avg = 1:size(Phi,6)
temp=Gr\reshape(Phi(:,:,1,1,:,avg),L,[]);

temp=(reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp,Ny,Nx,[],params.NEco));
cw=max(vec(abs(temp(:,:,:,1))));


for i = 1:size(temp, 3)
   temp_mask = roi_mask .* abs(temp(:,:,i,1));
   t1_array(i,avg) = mean(nonzeros(temp_mask(:)));
end
end

figure,
hold on;
for avg = 1:size(Phi,6)
    plot(t1_array(:,avg))
end
legend({'Seg1', 'Seg2', 'Seg3', 'Seg4'}, 'Location', 'SouthEast');
xlabel('Segment #');
grid on;
