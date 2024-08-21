
function mask = vert2mask(vertices,x,y,z)
% VERT2MASK - convert a set of 3D convex hull vertices into a 3D volume mask

[A,b] = vert2con(vertices);
p=[x(:),y(:),z(:)]';
p=bsxfun(@le,A*p,b); %applies the element-wise binary operation specified by the function handle @le to arrays A*p and b.
p=min(p); %Minimum elements of an array
mask=reshape(p,size(x));

end