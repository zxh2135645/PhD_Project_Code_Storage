function x = iwavelet2D(c,s,varargin)

narginchk(3,4)
x = appcoef2(conj(c(1,:)),s,varargin{:},0);

Nim = size(c,1);
if Nim > 1
    x(end,end,Nim) = 0; %preallocate
end

for j = 2:Nim
    x(:,:,j) = appcoef2(conj(c(j,:)),s,varargin{:},0);
end

return