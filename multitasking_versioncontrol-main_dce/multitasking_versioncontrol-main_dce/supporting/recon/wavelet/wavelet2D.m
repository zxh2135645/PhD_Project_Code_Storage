function [c,s] = wavelet2D(x,n,IN3,IN4)

% Check arguments.
if nargin == 3
    [Lo_D,Hi_D] = wfilters(IN3,'d');
else
    Lo_D = IN3;   Hi_D = IN4;
end

[c,s] = wavedec2(x(:,:,1),n,Lo_D,Hi_D);

Nim = size(x(:,:,:),3);
if Nim > 1
    c(Nim,end) = 0; % preallocate
end

for j = 2:Nim
    c(j,:) = wavedec2(x(:,:,j),n,Lo_D,Hi_D);
end

return