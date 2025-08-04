function dcminfo = genDicomInfo(twixObj, dcminfo)
% ========================================================================
%
%   Create dicom files from the image array and Siemens data header
%   Files will be saved in the same folder as the raw data file
%
%   You should be able to import the generated dicom files into dicom 
%   viewers like OsiriX, MicroDicom. If you intend to further process the 
%   images with tools like ITK, it's recommanded that you convert them 
%   to NIfTI format.
%
%   Tested with SoS data acquired on XA20A/VE11P systems
%   with patient position 'Head Fist Supine' ('HFS'), where
%
%           -x <--> +x = R <--> L
%           -y <--> +y = A <--> P
%           -z <--> +z = I <--> S
%
%   Output image top/bottom/left/right will be oriented to patient
%
%           S/I/A/P  (sagittal plane)
%           S/I/R/L  (coronal plane)
%           A/P/R/L  (transverse plane)
%
%   Input:
%
%       twixObj    - raw data struct extracted by mapVBVD()
%
%   Output:
%
%       dcminfo - dicom header information 
%
%
%  - last updated 2023-04-27
%  - Hsu-Lei Lee @ BIRI, Cedars-Sinai Medical Center
%
% ========================================================================

warning('off','all');

if nargin < 2
    dcminfo = [];
end

try

%% Calculate image position and dimensions

% Slice info
asSlice = twixObj.hdr.MeasYaps.sSliceArray.asSlice{1};

% for radial/spiral image, lPhaseResolution = lBaseResolution
lBaseResolution  = twixObj.hdr.Dicom.lBaseResolution;
if twixObj.hdr.Config.x2DInterpolation == 1     
    % 2DInterpolation was on
    lPhaseResolution = twixObj.hdr.Config.NImageLins / 2;
else
    lPhaseResolution = twixObj.hdr.Config.NImageLins;
end
if strcmp(twixObj.hdr.Dicom.tMRAcquisitionType, '2D')
    Nzorig = twixObj.hdr.MeasYaps.sSliceArray.lSize;
else
    try
        Nzorig = twixObj.hdr.Config.NImagePar;
    catch
        Nzorig = twixObj.hdr.Config.NoImagesPerSlab;
    end
end

% Spatial resolution
% pixelSpacing = [distance between rows, distance between columns] in mm
pixelSpacing = [asSlice.dPhaseFOV/lPhaseResolution;asSlice.dReadoutFOV/lBaseResolution]; 

% Pixel bandwidth
pixelBandwidth = round(1e9 / (twixObj.hdr.MeasYaps.sRXSPEC.alDwellTime{1}*lBaseResolution*twixObj.hdr.Meas.flReadoutOSFactor));

% Slice position: coordinate of the center of the slice/slab
if isfield(asSlice,'sPosition')
    sPosition = twixObj.hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition;
    if ~isfield(sPosition, 'dSag')
        sPosition.dSag = 0;
    end
    if ~isfield(sPosition, 'dCor')
        sPosition.dCor = 0;
    end    
    if ~isfield(sPosition, 'dTra')
        sPosition.dTra = 0;
    end
else
    sPosition.dSag = 0;
    sPosition.dCor = 0;
    sPosition.dTra = 0;
end
centerPosition = [sPosition.dSag; sPosition.dCor; sPosition.dTra];

% Slice normal vector length = slice thickness + gap            (2D)
%                            = slab thickness / # of partitions (3D)
if strcmp(twixObj.hdr.Dicom.tMRAcquisitionType, '2D')
    sliceThickness = asSlice.dThickness;
else
    sliceThickness = asSlice.dThickness/Nzorig;
end

% Column/row/slice vectors
RotMatQ = quat2RotMat(twixObj.image.slicePos(4:7,1));
vecCol  = RotMatQ*[1;0;0];
vecRow  = RotMatQ*[0;1;0];

%% Populate dicom info struct

% Enhanced MR Image Storage UID
dcminfo.SOPClassUID             = '1.2.840.10008.5.1.4.1.1.4';
dcminfo.MediaStorageSOPClassUID = dcminfo.SOPClassUID;

% Protocol info
dcminfo.SliceThickness          = sliceThickness;
dcminfo.PixelSpacing            = pixelSpacing;
dcminfo.PixelBandwidth          = pixelBandwidth;
dcminfo.ImageOrientationPatient = round([vecRow; vecCol], 6);
dcminfo.RepetitionTime          = twixObj.hdr.MeasYaps.alTR{1}/1000;  % ms
dcminfo.EchoTime                = twixObj.hdr.MeasYaps.alTE{1}/1000;  % ms
dcminfo.FlipAngle               = twixObj.hdr.Dicom.adFlipAngleDegree;
dcminfo.ProtocolName            = twixObj.hdr.Dicom.tProtocolName;
dcminfo.SequenceName            = twixObj.hdr.Config.SequenceString;
dcminfo.MRAcquisitionType       = twixObj.hdr.Dicom.tMRAcquisitionType;

dcminfo.centerPosition          = centerPosition;

% Study info
temp = split(twixObj.hdr.Config.FrameOfReference,'.');
dcminfo.StudyDate               = temp{end-3}(1:8);
dcminfo.StudyTime               = [temp{end-3}(9:14) '.000000'];
if isfield(twixObj.hdr.Meas, 'StudyInstanceUID')
    dcminfo.StudyInstanceUID    = twixObj.hdr.Meas.StudyInstanceUID;
else
    dcminfo.StudyInstanceUID    = twixObj.hdr.Config.FrameOfReference;
end
dcminfo.StudyDescription        = twixObj.hdr.Config.tStudyDescription;

% Series/Acquisition info
try
    if isfield(twixObj.hdr.MeasYaps,'tReferenceImage0')
        temp = split(twixObj.hdr.MeasYaps.tReferenceImage0,';');
    elseif isfield(twixObj.hdr.MeasYaps,'tReferenceImage1')
        temp = split(twixObj.hdr.MeasYaps.tReferenceImage1,';');
    end
    if length(temp) == 1    % VE11
        SeriesInstanceUID = temp{1};
    else                    % XA20
        SeriesInstanceUID = temp{2};
    end
    temp = split(SeriesInstanceUID,'.');
    dcminfo.SeriesTime      = [temp{end}(9:14) '.000000'];
    dcminfo.AcquisitionTime = [temp{end}(9:14) '.000000'];
catch
    SeriesInstanceUID       = twixObj.hdr.Config.FrameOfReference;
    dcminfo.SeriesTime      = dcminfo.StudyTime;
    dcminfo.AcquisitionTime = dcminfo.StudyTime;
end
dcminfo.SeriesInstanceUID       = SeriesInstanceUID;
dcminfo.SeriesDescription       = [twixObj.hdr.Config.SequenceDescription '_phase' sprintf('%02d',1)];
dcminfo.InstanceCreationTime    = dcminfo.AcquisitionTime;

% Equipment info
dcminfo.InstitutionName         = twixObj.hdr.Dicom.InstitutionName;
dcminfo.InstitutionAddress      = twixObj.hdr.Dicom.InstitutionAddress;
dcminfo.DeviceSerialNumber      = num2str(twixObj.hdr.Dicom.DeviceSerialNumber);
dcminfo.Manufacturer            = twixObj.hdr.Dicom.Manufacturer;
dcminfo.ManufacturerModelName   = twixObj.hdr.Dicom.ManufacturersModelName;
dcminfo.Modality                = twixObj.hdr.Dicom.Modality;
dcminfo.SoftwareVersion         = twixObj.hdr.Dicom.SoftwareVersions;
dcminfo.MagneticFieldStrength   = twixObj.hdr.Dicom.flMagneticFieldStrength;
if dcminfo.MagneticFieldStrength > 2.5 && dcminfo.MagneticFieldStrength < 3.5
    dcminfo.MagneticFieldStrength = 3;
end
dcminfo.ImagingFrequency        = twixObj.hdr.Dicom.lFrequency/1e6; % MHz
try
    dcminfo.ImagedNucleus       = twixObj.hdr.MeasYaps.sTXSPEC.asNucleusInfo{1}.tNucleus{1};
catch
    dcminfo.ImagedNucleus       = twixObj.hdr.MeasYaps.sRXSPEC.asNucleusInfo{1}.tNucleus;
end    

% Patient info
dcminfo.PatientPosition         = twixObj.hdr.Dicom.tPatientPosition;
dcminfo.PatientWeight           = round(twixObj.hdr.Dicom.flUsedPatientWeight,4);
dcminfo.PatientAge              = sprintf('%03dY',twixObj.hdr.Dicom.flPatientAge);
dcminfo.PatientSize             = twixObj.hdr.Meas.flPatientHeight/1000;
dcminfo.BodyPartExamined        = twixObj.hdr.Dicom.tBodyPartExamined;

dcminfo.PatientName = ['MID' sprintf('%04d',twixObj.hdr.Config.MeasUID)];
dcminfo.PatientID   = [dcminfo.StudyDate(1:8) sprintf('%04d',twixObj.hdr.Config.MeasUID)];


catch
fprintf('Extract Dicom info failed.\n');
end

warning('on','all');

end


%% ========================================================================
%  create 3D rotation matrix from quaternions
function RotMatQ = quat2RotMat(quat)

w = quat(1);
x = quat(2);
y = quat(3);
z = quat(4);
Rxx = 1 - 2 * (y^2 + z^2);
Rxy = 2 * (x*y - z*w);
Rxz = 2 * (x*z + y*w);
Ryx = 2 * (x*y + z*w);
Ryy = 1 - 2 * (x^2 + z^2);
Ryz = 2 * (y*z - x*w);
Rzx = 2 * (x*z - y*w);
Rzy = 2 * (y*z + x*w);
Rzz = 1 - 2 * (x^2 + y^2);
RotMatQ = [Rxx,    Rxy,    Rxz;
           Ryx,    Ryy,    Ryz;
           Rzx,    Rzy,    Rzz];
end

