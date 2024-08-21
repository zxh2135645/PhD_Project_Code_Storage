function [mask, fvc]=VObjEllipsoid(p,x,y,z,flag)
%create 3D cylinder virtual object

% Initialize parameters
xr = p.a;
yr = p.b;
zr = p.c;

FaceNum=p.FaceNum;
Length=p.Length;
CenterX=p.CenterX;
CenterY=p.CenterY;
CenterZ=p.CenterZ;

% Generate cylinder coordinate
[X,Y,Z] = ellipsoid(CenterX, CenterY, CenterZ, xr, yr, zr, FaceNum); 


fvc = surf2patch(X,Y,Z);

if flag~=1 % Render object only
    mask=0;
    return;
end

mask = vert2mask(fvc.vertices,x,y,z);
   
end