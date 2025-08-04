function bias = estimateBiasGaussian(fbp,SEs,voxelSpacing)

if nargin < 3
    voxelSpacing = [1 1 1];
end

if nargin == 2 && size(fbp) == size(SEs)
    fbpComposite = abs(sum(fbp.*conj(SEs),4));
else
    fbpComposite = sqrt(sum(abs(fbp).^2,4));
end
fbpComposite(~isfinite(fbpComposite)) = 0;
fbpComposite(fbpComposite==0) = min(fbpComposite(fbpComposite~=0));
fbpComposite = fbpComposite/prctile(abs(fbpComposite(:)),98);

[Ny,Nx,Nz] = size(fbpComposite);

sigma = [3 3 3];
sigma(1) = floor(Ny/16)*2 + 1;
sigma(2) = floor(sigma(1)*voxelSpacing(1)/voxelSpacing(2)/2)*2 + 1;
sigma(3) = floor(sigma(1)*voxelSpacing(1)/voxelSpacing(3)/2)*2 + 1;
filterSize = sigma*4 + 1;

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
mask1 = imopen(mask1,se);

se = strel("sphere",3);
mask1 = imdilate(mask1,se);

mask1 = imgaussfilt3(mask1,sigma/4,'Padding','symmetric','FilterSize',floor(filterSize/4));
mask2 = 1 - mask1;

bias = imgaussfilt3(fbpComposite,sigma,'Padding','symmetric','FilterSize',filterSize);
bias = bias / mean(bias(mask(:)==1));
bias = bias.^2;
bias = bias + mask2;
bias(bias>2) = 2;
bias(bias<0.2) = 0.2;



