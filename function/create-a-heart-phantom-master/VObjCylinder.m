
function [mask, fvc]=VObjCylinder(p,x,y,z, flag, theta, Groove)
%create 3D cylinder virtual object
if nargin < 6
    theta = 2*pi;
    Groove = 180;
end

% Initialize parameters
Radius=p.Radius;
FaceNum=p.FaceNum;
Length=p.Length;
CenterX=p.CenterX;
CenterY=p.CenterY;
CenterZ=p.CenterZ;

% Generate cylinder coordinate
[X,Y,Z] = cylinder(Radius,FaceNum); 
X=X+CenterX;
Y=Y+CenterY;
Z=(Z-0.5)*Length/2+CenterZ;

fvc = surf2patch(X,Y,Z);

if flag~=1 % Render object only
    mask=0;
    return;
end

mask = vert2mask(fvc.vertices,x,y,z);
   
for slc = 1:size(mask,3)
    Mask = mask(:,:,slc);
    [X,Y] = meshgrid(1:size(Mask,2),1:size(Mask,1));
    center=[sum(reshape(Mask.*X,1,[]))/sum(Mask(:)),sum(reshape(Mask.*Y,1,[]))/sum(Mask(:))];

    X=X-center(1);
    Y=Y-center(2);
    AngleMask=angle(X+1i*Y);
    
    nmax=theta+Groove/180*pi;
    nmin= 0 + Groove/180*pi;
    %project to pi and -pi
    nmaxproj=nmax-2*pi*floor((nmax+pi)/2/pi);
    nminproj=nmin-2*pi*floor((nmin+pi)/2/pi);

    if (nminproj>nmaxproj)
        nMask=Mask.*((AngleMask<=nmaxproj)+(AngleMask>=nminproj));
    else
        nMask=Mask.*(AngleMask<=nmaxproj).*(AngleMask>=nminproj);
    end

    mask(:,:,slc) = nMask;
end

end