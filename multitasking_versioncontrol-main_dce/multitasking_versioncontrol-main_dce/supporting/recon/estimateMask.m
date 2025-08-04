function mask2 = estimateMask(fbp,SEs,voxelSpacing)

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

% [a, b] = kmeans(fbpComposite(:),4);
% [~,idx] = min(b);
% a(a==idx) = 0;
% a(a>0) = 1;
% mask = reshape(a,Ny,Nx,Nz);
% se = strel("disk",2);
% mask1 = imdilate(mask,se);
% for n = 1:Nz
%     mask1(:,:,n) = imfill(mask1(:,:,n),"holes");
% end
% % mask1 = imopen(mask1,se);
% % mask1 = imerode(mask1,se);
% mask1 = imdilate(mask1,se);
mask1 = genMask3D(fbpComposite);

sigma = [3 3 3];
sigma(1) = floor(Ny/64)*2 + 1;
sigma(2) = floor(sigma(1)*voxelSpacing(1)/voxelSpacing(2)/2)*2 + 1;
sigma(3) = floor(sigma(1)*voxelSpacing(1)/voxelSpacing(3)/2)*2 + 1;
filterSize = sigma*4 + 1;

mask1 = imdilate(mask1,strel("disk",2));
mask1 = imgaussfilt3(mask1,sigma,'Padding','symmetric','FilterSize',floor(filterSize));
mask2 = 1 - mask1;



