function data = ift1d(data, dim, n)

if nargin > 2
    data = fftshift(ifft(ifftshift(data,dim),n,dim),dim)*sqrt(size(data,dim));
elseif nargin == 2
    data = fftshift(ifft(ifftshift(data,dim),[],dim),dim)*sqrt(size(data,dim));
else
    dim = find(size(data) > 1, 1);
    data = fftshift(ifft(ifftshift(data,dim),[],dim),dim)*sqrt(size(data,dim));
end