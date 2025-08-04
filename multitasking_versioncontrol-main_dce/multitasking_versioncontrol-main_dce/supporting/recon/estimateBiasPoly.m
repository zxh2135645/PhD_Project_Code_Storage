function bias = estimateBiasPoly(fbp,SEs,voxelSpacing)

fprintf('Estimating bias field... ');

if nargin < 3
    voxelSpacing = [1 1 1];
end
if nargin == 2 && numel(fbp) == numel(SEs)
    fbpComposite = abs(sum(fbp.*conj(SEs),4));
else
    fbpComposite = sqrt(sum(abs(fbp).^2,4));
end
fbpComposite = fbpComposite/max(abs(fbpComposite(:)));

[Ny,Nx,Nz] = size(fbpComposite);

[a, b] = kmeans(fbpComposite(:),4);
[~,idx] = min(b);
a(a==idx) = 0;
a(a>0) = 1;
mask = reshape(a,Ny,Nx,Nz);
se = strel("sphere",3);
mask1 = imdilate(mask,se);
for n = 1:Nz
    mask1(:,:,n) = imfill(mask1(:,:,n),"holes");
end
% mask1 = imopen(mask1,se);
% mask1 = imerode(mask1,se);
mask1 = imdilate(mask1,se);

sigma = [3 3 3];
sigma(1) = floor(Ny/64)*2 + 1;
sigma(2) = floor(sigma(1)*voxelSpacing(1)/voxelSpacing(2)/2)*2 + 1;
sigma(3) = floor(sigma(1)*voxelSpacing(1)/voxelSpacing(3)/2)*2 + 1;
filterSize = sigma*4 + 1;

[x,y,z] = ndgrid(-Ny/2:Ny/2-1,-Nx/2:Nx/2-1,-Nz/2:Nz/2-1); c = ones(size(x(:)));

% x(mask1==0) = 0;
% y(mask1==0) = 0;
% z(mask1==0) = 0;
% c(mask1(:)==0) = 0; %mean(fbpComposite(mask(:)==0))/mean(fbpComposite(mask(:)==1));

polyterm       = [x(:).^2, y(:).^2, z(:).^2, x(:).*y(:), x(:).*z(:), y(:).*z(:), x(:), y(:), z(:), c(:)]; 
[polyterm1,~,~] = svd(bsxfun(@times,polyterm,fbpComposite(:).^2),'econ');
polyterm = polyterm1;

polyterm(mask1(:)==0,:) = 0;
%polyterm(mask1(:)==0,end) = fbpComposite(mask1(:)==0).^2.*mean(fbpComposite(mask1(:)==0).^2)./mean(fbpComposite(mask1(:)==1).^2);

w   = polyterm'*fbpComposite(:);
yy   = zeros(size(fbpComposite(:)));
rho = 1/max(abs(polyterm*w-fbpComposite(:)));
for it = 1:10
    res = polyterm*w-fbpComposite(:);
    zz   = res+yy/rho;
    zz   = sign(zz).*max(abs(zz)-1/rho,0);
    yy   = yy+rho*(res-zz);
    rho = rho*1.05;
    w   = polyterm'*(zz+fbpComposite(:)-yy/rho);
    %     norm(polyterm*w-pd(:),1)
end
[~,~,v] = svd([(polyterm*w)./fbpComposite(:), fbpComposite(:)],'econ');
w       = w*v(2)/v(1);

bias = reshape(polyterm1*w,Ny,Nx,Nz).\fbpComposite.^2;
bias = bias/mean(bias(mask(:)>0));
bias(1./bias<1/2) = 2;
bias(bias<1/5)    = 1/5;
bias = bias/mean(bias(mask(:)>0));
%bias = imgaussfilt3(bias,sigma/2,'Padding','symmetric','FilterSize',floor(filterSize));

mask1 = imdilate(mask1,se);
mask1 = imgaussfilt3(mask1,sigma,'Padding','symmetric','FilterSize',floor(filterSize));
mask2 = 1 - mask1;

bias = bias + mask2;
bias = bias/mean(bias(mask(:)>0));
%bias = imgaussfilt3(bias,sigma/2,'Padding','symmetric','FilterSize',floor(filterSize));

bias(1./bias<1/2) = 2;
bias(bias<1/5)    = 1/5;

fprintf('done.\n');


