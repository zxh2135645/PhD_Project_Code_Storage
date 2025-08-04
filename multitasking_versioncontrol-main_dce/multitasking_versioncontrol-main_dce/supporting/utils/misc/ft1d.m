function data = ft1d(data, dim, n)

if nargin > 2
    data = fftshift(fft(ifftshift(data,dim),n,dim),dim)/sqrt(size(data,dim));
elseif nargin == 2
    data = fftshift(fft(ifftshift(data,dim),[],dim),dim)/sqrt(size(data,dim));
else
    dim = find(size(data) > 1, 1);
    data = fftshift(fft(ifftshift(data,dim),[],dim),dim)/sqrt(size(data,dim));
end