function bias = estimateBiasPoly(fbp,SEs)

if nargin == 2 && size(fbp) == size(SEs)
    fbpComposite = sqrt(sum(fbp.*conj(SEs),4));
else
    fbpComposite = sqrt(sum(abs(fbp).^2,4));
end
[Ny,Nx,Nz] = size(fbpComposite);

[x,y,z]        = ndgrid(-Ny/2:Ny/2-1,-Nx/2:Nx/2-1,-Nz/2:Nz/2-1);
polyterm       = [x(:).^2, y(:).^2, z(:).^2, x(:).*y(:), x(:).*z(:), y(:).*z(:), x(:), y(:), z(:), ones(size(x(:)))]; 
[polyterm,~,~] = svd(bsxfun(@times,polyterm,fbpComposite(:).^2),'econ');
w   = polyterm'*fbpComposite(:);
y   = zeros(size(fbpComposite(:)));
rho = 1/max(abs(polyterm*w-fbpComposite(:)));
for it = 1:10
    res = polyterm*w-fbpComposite(:);
    z   = res+y/rho;
    z   = sign(z).*max(abs(z)-1/rho,0);
    y   = y+rho*(res-z);
    rho = rho*1.05;
    w   = polyterm'*(z+fbpComposite(:)-y/rho);
    %     norm(polyterm*w-pd(:),1)
end
[~,~,v] = svd([(polyterm*w)./fbpComposite(:), fbpComposite(:)],'econ');
w       = w*v(2)/v(1);

bias    = reshape(polyterm*w,Ny,Nx,Nz).\fbpComposite.^2;
bias(1./bias<1/2) = 2;
bias(bias<1/4)    = 1/4;