function Rz=zrot(phi)

% use right-handed convention 
Rz = [ cos(phi) -sin(phi) 0; sin(phi) cos(phi) 0; 0 0 1 ];
