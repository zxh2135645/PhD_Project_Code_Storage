clear BzT Bzresize Bzforshim
%% ResolutionMatch

if UNIC_B0Map.Parameters.PhaseDir==BzMap.Parameters.PhaseDir % BZ is in COL
    %do nothing
else
    % flip readout and phase direction
    UNIC_B0Map.Parameters.Mask=permute(UNIC_B0Map.Parameters.Mask, [2 1 3]);
    UNIC_B0Map.Parameters.MatrixSize=(UNIC_B0Map.Parameters.MatrixSize([2 1 3]));
    UNIC_B0Map.Parameters.Voxelsize=(UNIC_B0Map.Parameters.Voxelsize([2 1 3]));
    UNIC_B0Map.Parameters.CornerCoordinate=UNIC_B0Map.Parameters.CornerCoordinate([2 1 3]);
    UNIC_B0Map.Map=permute(UNIC_B0Map.Map, [2 1 3]);
end
dim_B0=size(UNIC_B0Map.Map); 
B0voxel=single(UNIC_B0Map.Parameters.Voxelsize);
B0Center=single(UNIC_B0Map.Parameters.CornerCoordinate+ ...
    UNIC_B0Map.Parameters.Voxelsize.*(UNIC_B0Map.Parameters.MatrixSize'...
    ./[2;1;1]/2 ... % correct for Readout oversampling
    +[0;0;0])); % Correct for 0.5 pixel shift in slicer direction
%B0Center=[0;0;-140]
B0CenterPix=dim_B0(1:3)/2;

Bzin=BzMap.Maps.Bz2_5m1_5;
Bzin(isnan(Bzin))=0;
dim_Bz=size(Bzin);    
%Bzvoxel=[3.57; 3.57; 5.79];
%Bzsize=
%Bzcenter=[0;0;0];
Bzvoxel=single(BzMap.Parameters.Voxelsize);
BzCenter=single(BzMap.Parameters.CornerCoordinate+...
    BzMap.Parameters.Voxelsize.*(BzMap.Parameters.MatrixSize'...
    ./[2;1;1]/2 ... % correct for Readout oversampling
    +[0;0;0])); % Correct for 0.5 pixel shift in slicer direction

BzShift=(rem((B0Center-BzCenter),B0voxel)*-1)./B0voxel;
TranslateMatrix=eye(4);
TranslateMatrix(4,1:3)=BzShift';
Tform=affine3d(TranslateMatrix);
BzTCenter=BzCenter-BzShift;
BzCenterPix=dim_Bz(1:3)/2;

Bz2B0scale=[B0voxel(1)/Bzvoxel(1) B0voxel(2)/Bzvoxel(2) B0voxel(3)/Bzvoxel(3)];

B0FOVStart=BzCenterPix+single(((UNIC_B0Map.Parameters.CornerCoordinate-B0Center).*[2;1;1]-BzCenter)...% correct for readoutfov
    ./UNIC_B0Map.Parameters.Voxelsize)'+1;
B0FOVEnd=B0FOVStart+dim_B0(1:3)-1;
%BzFOV=[B0FOVStart(1):B0FOVEnd(1),B0FOVStart(2):B0FOVEnd(2),B0FOVStart(3):B0FOVEnd(3)];
for j=1:dim_Bz(4)
    BzT(:,:,:,j)=imwarp(Bzin(:,:,:,j),Tform);

    BZresize(:,:,:,j)=imresize3(BzT(:,:,:,j),...
        [dim_Bz(1)/Bz2B0scale(1) dim_Bz(2)/Bz2B0scale(2) dim_Bz(3)/Bz2B0scale(3)]);
    Bztemp=BZresize(:,:,:,j);
    Bzforshim(:,:,:,j)=Bztemp(B0FOVStart(1):B0FOVEnd(1), ...
        B0FOVStart(2):B0FOVEnd(2),B0FOVStart(3):B0FOVEnd(3));
end
Bz=Bzforshim;

%% Location match



