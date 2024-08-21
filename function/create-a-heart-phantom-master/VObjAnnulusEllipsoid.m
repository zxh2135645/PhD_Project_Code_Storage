
function [mask, fvc1, fvc2]=VObjAnnulusEllipsoid(p,x,y,z,theta, Groove, flag)
%create 3D cylinder virtual object

% Initialize parameters
% Initialize parameters
xr1 = p.a1;
yr1 = p.b1;
zr1 = p.c1;
xr2 = p.a2;
yr2 = p.b2;
zr2 = p.c2;

FaceNum=p.FaceNum;
Length=p.Length;
CenterX=p.CenterX;
CenterY=p.CenterY;
CenterZ=p.CenterZ;

% Generate Ellipsoid coordinate
[X1,Y1,Z1] = ellipsoid(CenterX, CenterY, CenterZ, xr1, yr1, zr1, FaceNum); 


% Generate Ellipsoid coordinate
[X2,Y2,Z2] = ellipsoid(CenterX, CenterY, CenterZ, xr2, yr2, zr2, FaceNum); 


fvc1 = surf2patch(X1,Y1,Z1);
fvc2 = surf2patch(X2,Y2,Z2);


mask1 = vert2mask(fvc1.vertices,x,y,z);
mask2 = vert2mask(fvc2.vertices,x,y,z);

mask = mask2 - mask1;
%theta = pi/2;
%Groove = 0;

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

if flag~=1 % Render object only
    mask=0;
    return;
end

end