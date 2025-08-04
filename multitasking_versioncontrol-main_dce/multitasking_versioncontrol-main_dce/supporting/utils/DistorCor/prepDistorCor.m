function prepDistorCor(prefix,img,twix_obj,hdr)
% Create image binary and header file required for distortion correction
% Files will be saved as [prefix '_DistorCorHdr']
% Tested on data acquired with software version XA20A and VE11P
%
% Input:
%
%   prefix      - file name prefix for image and header
%   img         - 2D/3D/4D image array (Npe x Nro x Nz x Nt)
%   twix_obj    - scan information extracted with mapVBVD
%   hdr         - header struct containing the info for the first slice
%   (optional)    generated in genDicomfromSiemens()
%
%
%   To have the correct image orientation, the data should be reconstructed 
% 	as described below:
%
%   1)  Apply complex conjugation to kspaceData before recon if patient 
% 		position is head-first. Do kspaceData  = conj(kspaceData) right 
%       after loading data from twix object.
%       
%   2)  For radial data, the trajectory om for NUFFT should be set as
%
%           r   = linspace(-pi, pi, Nkx+1); r(end)=[];
%           om1 = sin(thetas(1:Ntrajs)*pi/180)*r;
%           om2 = cos(thetas(1:Ntrajs)*pi/180)*r;
%           om  = [om1(:), om2(:)];
%
%       where thetas is the dAzimuthalAngle during data acquisition
%
%   -- Hsu-Lei Lee, last edited on 10/20/2020


%% Extract header info
hdr.filenameHdr = [prefix '_DistorCorHdr'];

asSlice = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1};
if ~isfield(hdr, 'SliceThickness')
    if strcmp(twix_obj.hdr.Dicom.tMRAcquisitionType, '2D')
        hdr.SliceThickness = asSlice.dThickness;
    else
        hdr.SliceThickness = asSlice.dThickness/twix_obj.hdr.Config.NoImagesPerSlab;
    end
end

if ~isfield(hdr, 'vecNorm')
    % 3D rotation of column and row vectors
    RotMatQ = quat2RotMat(twix_obj.image.slicePos(4:7,1));
    hdr.vecCol  = RotMatQ*[1;0;0];
    hdr.vecRow  = RotMatQ*[0;1;0];
    hdr.vecNorm = RotMatQ*[0;0;1];
    
    % slice normal vector length = slice thickness + gap            (2D multi-slice)
    %                            = slice thickness                  (3D)
    if strcmp(twix_obj.hdr.Dicom.tMRAcquisitionType, '2D') && isfield(twix_obj.hdr.MeasYaps.sGroupArray.asGroup{1},'dDistFact')
        vecNormScale = (1 + twix_obj.hdr.MeasYaps.sGroupArray.asGroup{1}.dDistFact) * hdr.SliceThickness;
    else
        vecNormScale = hdr.SliceThickness;
    end
    % scaled normal vector (shift between consecutive slices)
    hdr.vecNorm = hdr.vecNorm*vecNormScale; 
end

if ~isfield(hdr, 'vecPos')
    % slice array
    if isfield(asSlice,'sPosition')
        sPosition = twix_obj.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition;
    else
        sPosition.dSag = 0;
        sPosition.dCor = 0;
        sPosition.dTra = 0;
    end

    % position vector of the first slice/partition
    if isfield(sPosition,'dSag')
        hdr.vecPos(1) = sPosition.dSag;
    else
        hdr.vecPos(1) = 0;
    end
    if isfield(sPosition,'dCor')
        hdr.vecPos(2) = sPosition.dCor;
    else
        hdr.vecPos(2) = 0;
    end
    if isfield(sPosition,'dTra')
        hdr.vecPos(3) = sPosition.dTra;
    else
        hdr.vecPos(3) = 0;
    end
    if strcmp(twix_obj.hdr.Dicom.tMRAcquisitionType, '3D')
        hdr.vecPos = hdr.vecPos - hdr.vecNorm*(size(img,3)-1)/2;
    end
end

if ~isfield(hdr, 'pixelSpacing')
    hdr.pixelSpacing = [asSlice.dPhaseFOV/twix_obj.hdr.Config.NImageLins; asSlice.dReadoutFOV/twix_obj.hdr.Dicom.lBaseResolution];    % mm
end

dReadoutFOV = hdr.pixelSpacing(2)*size(img,2);
dPhaseFOV   = hdr.pixelSpacing(1)*size(img,1);

temp = hdr.vecRow;
hdr.vecRow = hdr.vecCol;
hdr.vecCol = temp;

% write to header file
fileID      = fopen(hdr.filenameHdr,'w');

fprintf(fileID,'iWidth=%d\n' , size(img,1));
fprintf(fileID,'iHeight=%d\n', size(img,2));
fprintf(fileID,'iSlices=%d\n', size(img,3));
fprintf(fileID,'NR=%d\n'     , size(img,4));
fprintf(fileID,'fovRead=%f\n' , dReadoutFOV);
fprintf(fileID,'fovPhase=%f\n', dPhaseFOV);
fprintf(fileID,'fovPar=%f\n'  , hdr.SliceThickness);
fprintf(fileID,'vecPosX=%.20f\n', hdr.vecPos(1));
fprintf(fileID,'vecPosY=%.20f\n', hdr.vecPos(2));
fprintf(fileID,'vecPosZ=%.20f\n', hdr.vecPos(3));
fprintf(fileID,'vecRowX=%.20f\n', hdr.vecRow(1));
fprintf(fileID,'vecRowY=%.20f\n', hdr.vecRow(2));
fprintf(fileID,'vecRowZ=%.20f\n', hdr.vecRow(3));
fprintf(fileID,'vecColX=%.20f\n', hdr.vecCol(1));
fprintf(fileID,'vecColY=%.20f\n', hdr.vecCol(2));
fprintf(fileID,'vecColZ=%.20f\n', hdr.vecCol(3));
fprintf(fileID,'vecNormX=%.20f\n', hdr.vecNorm(1));
fprintf(fileID,'vecNormY=%.20f\n', hdr.vecNorm(2));
fprintf(fileID,'vecNormZ=%.20f\n', hdr.vecNorm(3));
tempStr = ['fnameGradCoilCfg=coeff_' twix_obj.hdr.Dicom.tGradientCoil '.grad\n'];   fprintf(fileID,tempStr);
tempStr = ['pcPatientPos=' twix_obj.hdr.Dicom.tPatientPosition '\n'];               fprintf(fileID,tempStr);
fprintf(fileID,'##END');
                    
fclose(fileID);


%% Write image binary
%img_float   = single(abs(img));
filenameImg = [prefix '_Img'];
writeImgBinary(filenameImg,single(abs(img)),'single');

end

%% convert quaternion to 3x3 rotation matrix
function RotMatQ = quat2RotMat(quat)
w = quat(1);
x = quat(2);
y = quat(3);
z = quat(4);
Rxx = 1 - 2*(y^2 + z^2);
Rxy = 2*(x*y - z*w);
Rxz = 2*(x*z + y*w);
Ryx = 2*(x*y + z*w);
Ryy = 1 - 2*(x^2 + z^2);
Ryz = 2*(y*z - x*w );
Rzx = 2*(x*z - y*w );
Rzy = 2*(y*z + x*w );
Rzz = 1 - 2 *(x^2 + y^2);
RotMatQ = [Rxx,    Rxy,    Rxz;
           Ryx,    Ryy,    Ryz;
           Rzx,    Rzy,    Rzz];
end