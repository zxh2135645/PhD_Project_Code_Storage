function [image,dcminfo,hdr] = genDicomfromSiemens(image, twixObj, interpFactor, flagNormalize, filePrefix, flagDummyName, seriesNumber, volumeSeries, tempResol, dcminfo, mainpath)
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
%       image      - 3D or 4D image (Npe x Nro x Nz x Nt)
%
%       twixObj    - raw data struct extracted by mapVBVD()
%
%       interpFactor  (optional 1x2 double array, default is [1 1]) 
%                  - if the input image has been interpolated, the image size 
% 					 should be interpFactor*[Nydisp x Nxdisp]
%                    PixelSpacing in the dicom header will be scaled accordingly
%
%       flagNormalize (optional bool, default is true) 
%                  - true : rescale image pixel value to [0,4095]
%                    false: don't rescale image
%
%       filePrefix    (optional char array, default is [MID + protocol name]) 
%                  - dicom file name prefix
%
%       flagDummyName (optional bool, default is true)
%                  - true : use MID to generate dummy patient name and ID
%                    false: use patient name and ID in twixObj
%
%       seriesNumber  (optional positive int, default is 1)
%                  - change SeriesInstanceUID to manually create new series
%                    for dicom viewers
%
%       volumeSeries  (optional int, default is 1)
%                  - 0: save Nz dicoms, each has Nt frames
%                    1: save Nt dicoms, each contains a 3D volume
%                    2: save NzxNt dicoms
%
%       tempResol     (optional positive double, default is 50.0)
%                  - temporal resolution (ms) between consecutive frames
%
%       mainpath      (optional char array)   
%                  - path to the multitasking recon tool, specified at 
%                    the beginning of multitasking_recon.m. If provided, 
%                    2D distortion correction will be applied.
%                    Make sure you have the DistorCor package in the utils 
%                    folder and the binary is executable.
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
%       For Cartesian data, the image array should be oriented as
%       phase-encoding x readout x slice x time
%
%  Output:
%
%       image   - flipped/cropped/rescaled 3D/4D image 
%                 (interpFactor*(lPhaseResolution x lBaseResolution) x Nzorig x Nframe)
%       dcminfo - dicom header information corresponding to the first
%                 slice of image
%       hdr     - slice information of image to be used in run_DistorCor()
%
%
%  - last updated 2022-09-23
%  - Hsu-Lei Lee @ BIRI, Cedars-Sinai Medical Center
%
% ========================================================================

warning('off','all');

% Temporal resolution, used if there are multiple time frames
if nargin < 9
    tempResol = 1000/size(image,4);
end

% Generate series of 3D volumes or multi-frame single-slice files?
% 1: 3D volumes; 0: single-slice images
if nargin < 8
    volumeSeries = 1;
end

% Series number
if nargin < 7
    seriesNumber = 1;
elseif seriesNumber < 1
    seriesNumber = 1;
end

% flagDummyName
if nargin < 6
    flagDummyName = 1;
elseif seriesNumber < 1
    flagDummyName = 0;
end

% Dicom file name and path
if nargin < 5 
    filePrefix = sprintf('MID%d_%s',twixObj.hdr.Meas.MeasUID,twixObj.hdr.Dicom.tProtocolName);
elseif isempty(filePrefix) 
    filePrefix = sprintf('MID%d_%s',twixObj.hdr.Meas.MeasUID,twixObj.hdr.Dicom.tProtocolName);
else
    filePrefix = sprintf('MID%d_%s_%s',twixObj.hdr.Meas.MeasUID,twixObj.hdr.Dicom.tProtocolName,filePrefix);
end

try
    %filePath = [twixObj.image.filename(1:end-4) '_dicom'];
    filePath = [fileparts(twixObj.image.filename) '/dicom'];
    if ~exist(filePath,'dir')
        mkdir(filePath)
    end
catch
    filePath = 'dicom';
    if ~exist(filePath,'dir')
        mkdir(filePath)
    end
end
filePrefix = [filePath '/' filePrefix];

% Rescale image?
if nargin < 4
    flagNormalize = true;
elseif isempty(flagNormalize)
    flagNormalize = true;
end
    
if nargin < 3
    interpFactor = [1,1];
elseif isempty(interpFactor)
    interpFactor = [1,1];
end

fprintf('Processing images... ');

% apply 2D distortion correction if mainpath is provided as a funtion input
if nargin >= 11
    try
        image = runDistorCor(image, mainpath, twixObj);
    catch
        disp('2D distortion correction failed. No correction was performed.');
    end
end

% Crop the image if necessary
% for radial/spiral image, lPhaseResolution = lBaseResolution
lBaseResolution  = twixObj.hdr.Dicom.lBaseResolution;
if twixObj.hdr.Config.x2DInterpolation == 1     
    % 2DInterpolation was on
    lPhaseResolution = twixObj.hdr.Config.NImageLins / 2;
else
    lPhaseResolution = twixObj.hdr.Config.NImageLins;
end
if strcmp(twixObj.hdr.Dicom.tMRAcquisitionType, '2D')
    Nzorig = size(image,3);
else
    try
        Nzorig = twixObj.hdr.Config.NImagePar;
    catch
        Nzorig = twixObj.hdr.Config.NoImagesPerSlab;
    end
end

image = padcrop(image, 1, lPhaseResolution*interpFactor(1));
image = padcrop(image, 2, lBaseResolution *interpFactor(2));
image = padcrop(image, 3, Nzorig);

% Normalize pixel values to [0,4095] if flagNormalize is on
if flagNormalize
    image = 4095 * image/max(abs(image(:)));
end


%% Calculate image position and dimensions

% Slice info
asSlice = twixObj.hdr.MeasYaps.sSliceArray.asSlice{1};

% Spatial resolution
% pixelSpacing = [distance between rows, distance between columns] in mm
pixelSpacing = [asSlice.dPhaseFOV/lPhaseResolution/interpFactor(1);asSlice.dReadoutFOV/lBaseResolution/interpFactor(2)]; 

% Pixel bandwidth
pixelBandwidth = round(1e9 / (twixObj.hdr.MeasYaps.sRXSPEC.alDwellTime{1}*lBaseResolution*twixObj.hdr.Meas.flReadoutOSFactor) / interpFactor(2));


% Slice position: coordinate of the center of the slice/slab
if isfield(asSlice,'sPosition')
    sPosition = asSlice.sPosition;
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
    if isfield(twixObj.hdr.MeasYaps.sGroupArray.asGroup{1},'dDistFact')
        vecNormScale = (1 + twixObj.hdr.MeasYaps.sGroupArray.asGroup{1}.dDistFact) * sliceThickness;
    else
        vecNormScale = sliceThickness;
    end
else
    sliceThickness = asSlice.dThickness/Nzorig;
    vecNormScale   = sliceThickness;
end

% Column/row/slice vectors
RotMatQ = quat2RotMat(twixObj.image.slicePos(4:7,1));
vecCol  = RotMatQ*[1;0;0];
vecRow  = RotMatQ*[0;1;0];
vecNorm = RotMatQ*[0;0;1]*vecNormScale;

% Rotate and flip the image if necessary
[~, Mdir] = max(abs(vecNorm));
[~, Midx] = max(abs([vecRow vecCol]), [], 2);

if Mdir == 1        % sagittal plane
    if Midx(3) == 1     % align the column vector to S/I direction
        vecCol  = RotMatQ*[0;1;0];
        vecRow  = RotMatQ*[1;0;0];
        image = permute(image,[2 1 3 4]);
        pixelSpacing = flip(pixelSpacing);
    end
    if vecRow(2) < 0    % align image left/right to A/P
        vecRow = -vecRow;
        image = flip(image,2);
    end
    if vecCol(3) > 0    % align image top/bottom to S/I
        vecCol = -vecCol;
        image = flip(image,1);
    end 
elseif Mdir == 2    % coronal plane
    if Midx(3) == 1     % align the column vector to S/I direction
        vecCol  = RotMatQ*[0;1;0];
        vecRow  = RotMatQ*[1;0;0];
        image = permute(image,[2 1 3 4]);
        pixelSpacing = flip(pixelSpacing);
    end
    if vecRow(1) < 0    % align image left/right to R/L
        vecRow = -vecRow;
        image = flip(image,2);
    end
    if vecCol(3) > 0    % align image top/bottom to S/I
        vecCol = -vecCol;
        image = flip(image,1);
    end 
elseif Mdir == 3    % transverse plane
    if Midx(2) == 1     % align the column vector to A/P direction
        vecCol  = RotMatQ*[0;1;0];
        vecRow  = RotMatQ*[1;0;0];
        image = permute(image,[2 1 3 4]);
        pixelSpacing = flip(pixelSpacing);
    end
    if vecRow(1) < 0    % align image left/right to R/L
        vecRow = -vecRow;
        image = flip(image,2);
    end
    if vecCol(2) < 0    % align image top/bottom to A/P
        vecCol = -vecCol;
        image = flip(image,1);
    end 
end

% The location information of the first slice
% ImagePositionPatient1: coordinates of the first (top-left) pixel
% SliceLocation1: perpendicular distance from isocenter to the slice plane
% vecPos1: coordinate of the center of the slice
[Ny, Nx, Nz, Nt] = size(image);
if strcmp(twixObj.hdr.Dicom.tMRAcquisitionType, '2D')
    ImagePositionPatient1 = centerPosition - vecCol*pixelSpacing(1)*Ny/2 - vecRow*pixelSpacing(2)*Nx/2;
    SliceLocation1 = dot(centerPosition, vecNorm/vecNormScale);
    vecPos1 = centerPosition;
else
    ImagePositionPatient1 = centerPosition - vecCol*pixelSpacing(1)*Ny/2 - vecRow*pixelSpacing(2)*Nx/2 - vecNorm*(Nz-1)/2;
    SliceLocation1 = dot(centerPosition - vecNorm*(Nz-1)/2, vecNorm/vecNormScale);
    vecPos1 = centerPosition - vecNorm*(Nz-1)/2;
end

% Store info of the new slice orientation 
% in case distortion correction is needed later
hdr.vecCol  = vecCol;
hdr.vecRow  = vecRow;
hdr.vecNorm = vecNorm;
hdr.vecPos  = vecPos1;
hdr.pixelSpacing   = pixelSpacing;
hdr.sliceThickness = sliceThickness;


%% Populate dicom info struct

% Enhanced MR Image Storage UID
dcminfo.SOPClassUID             = '1.2.840.10008.5.1.4.1.1.4';
dcminfo.MediaStorageSOPClassUID = dcminfo.SOPClassUID;

% Protocol info
dcminfo.SliceThickness          = sliceThickness;
dcminfo.PixelSpacing            = pixelSpacing;
dcminfo.PixelBandwidth          = pixelBandwidth;
dcminfo.ImageOrientationPatient = round([vecRow; vecCol], 6);
dcminfo.ImagePositionPatient    = round(ImagePositionPatient1, 4);
dcminfo.SliceLocation           = round(SliceLocation1, 5);
dcminfo.RepetitionTime          = twixObj.hdr.MeasYaps.alTR{1}/1000;  % ms
dcminfo.EchoTime                = twixObj.hdr.MeasYaps.alTE{1}/1000;  % ms
dcminfo.FlipAngle               = twixObj.hdr.Dicom.adFlipAngleDegree;
dcminfo.ProtocolName            = twixObj.hdr.Dicom.tProtocolName;
dcminfo.SequenceName            = twixObj.hdr.Config.SequenceString;
dcminfo.MRAcquisitionType       = twixObj.hdr.Dicom.tMRAcquisitionType;

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
dcminfo.SeriesDescription       = [twixObj.hdr.Config.SequenceDescription '_phase' sprintf('%02d',seriesNumber)];

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

if flagDummyName    % use MID as dummy patient name/ID
    dcminfo.PatientName = ['MID' sprintf('%04d',twixObj.hdr.Config.MeasUID)];
    dcminfo.PatientID   = [dcminfo.StudyDate(1:8) sprintf('%04d',twixObj.hdr.Config.MeasUID)];
else
    dcminfo.PatientName = twixObj.hdr.Dicom.tPatientName;
    dcminfo.PatientID   = twixObj.hdr.Config.PatientID;
end

%% Export dicom files

fprintf('Saving dicom to %s/ ', filePath);
tempinfo = dcminfo;
if volumeSeries == 0            % export slice-by-slice multi-frame images
    for slice = 1:Nz
        % current slice
        tempinfo.ImagePositionPatient = round(ImagePositionPatient1 + vecNorm*(slice-1), 4);
        tempinfo.SliceLocation        = round(SliceLocation1 + sliceThickness*(slice-1), 5);
        tempinfo.NumberOfFrames       = Nt;
        tempinfo.NumberOfTemporalPositions = Nt;
        for frame = 1:Nt        % time stamp for each frame
            tempstr = ['tempinfo.PerFrameFunctionalGroupsSequence.Item_' num2str(frame) '.CardiacTriggerSequence.Item_1.CardiacTriggerDelayTime = (' num2str(frame) '-1)*' num2str(round(tempResol)) ';'];
            eval(tempstr);
        end
        % file name for current slice
        fileName = sprintf('%s_sl%03d.dcm', filePrefix, slice);
        % write multi-frame image to dicom file
        dicomwrite(uint16(round(abs(image(:,:,slice,:)))), fileName, tempinfo, 'CreateMode', 'Copy');        
        if mod(slice,round(Nz/10)+1) == 0
            fprintf('.');
        end
    end
elseif volumeSeries == 1      % export frame-by-frame 3D volumes
    hh = str2double(dcminfo.AcquisitionTime(1:2));
    mm = str2double(dcminfo.AcquisitionTime(3:4));
    ss = str2double(dcminfo.AcquisitionTime(5:6));
    timeStart = hh*60*60 + mm*60 + ss;
    for frame = 1:Nt
        % current frame
        temp = num2str(1000 + uint16(Nt*(seriesNumber-1)));
        SeriesInstanceUID = [SeriesInstanceUID(1:end-4) temp];
        tempinfo.SeriesInstanceUID = [SeriesInstanceUID '.0.0.0'];
        tempinfo.SeriesNumber      = seriesNumber;
        timeCurrent = timeStart + tempResol * 1e-3 * (frame-1);
        hh = floor(timeCurrent/3600);
        mm = mod(floor(timeCurrent/60),60);
        ss = mod(floor(timeCurrent),60);
        timeDigit = floor((timeCurrent - floor(timeCurrent)) * 1e6);
        tempinfo.AcquisitionTime   = sprintf('%02d%02d%02d.%06d',hh,mm,ss,timeDigit);
        tempinfo.AcquisitionNumber = frame;
        tempinfo.InstanceCreationTime = tempinfo.AcquisitionTime;
        tempinfo.InstanceNumber       = frame;
        for slice = 1:Nz        % location info for each slice
            tempstr = ['tempinfo.PerFrameFunctionalGroupsSequence.Item_' num2str(slice) '.PlanePositionSequence.Item_1.ImagePositionPatient = round(ImagePositionPatient1 + vecNorm*(' num2str(slice) '-1), 4);'];
            eval(tempstr);
        end

        % file name for current volume
        fileName = sprintf('%s_v%03d.dcm', filePrefix, frame);
        % write 3D volume to dicom file
        dicomwrite(uint16(round(abs(reshape(image(:,:,:,frame),Ny,Nx,1,[])))), fileName, tempinfo, 'CreateMode', 'Copy');        
        if mod(frame,round(frame/10)+1) == 0
            fprintf('.');
        end
    end
else                            % export single-slice, single-frame images
    hh = str2double(dcminfo.AcquisitionTime(1:2));
    mm = str2double(dcminfo.AcquisitionTime(3:4));
    ss = str2double(dcminfo.AcquisitionTime(5:6));
    timeStart = hh*60*60 + mm*60 + ss;
    for slice = 1:Nz
        % current slice
        tempinfo.ImagePositionPatient = round(ImagePositionPatient1 + vecNorm*(slice-1), 4);
        tempinfo.SliceLocation        = round(SliceLocation1 + sliceThickness*(slice-1), 5);
        for frame = 1:Nt        % time stamp for each frame
            % current frame
            temp = num2str(1000 + uint16(Nt*(seriesNumber-1)));
            SeriesInstanceUID = [SeriesInstanceUID(1:end-4) temp];
            tempinfo.SeriesInstanceUID = [SeriesInstanceUID '.0.0.0'];
            tempinfo.SeriesNumber      = seriesNumber;
            timeCurrent = timeStart + tempResol * 1e-3 * (frame-1);
            hh = floor(timeCurrent/3600);
            mm = mod(floor(timeCurrent/60),60);
            ss = mod(floor(timeCurrent),60);
            timeDigit = floor((timeCurrent - floor(timeCurrent)) * 1e6);
            tempinfo.AcquisitionTime   = sprintf('%02d%02d%02d.%06d',hh,mm,ss,timeDigit);
            tempinfo.AcquisitionNumber = frame;
            tempinfo.InstanceCreationTime = tempinfo.AcquisitionTime;
            tempinfo.InstanceNumber       = frame;
            % file name for current slice
            fileName = sprintf('%s_s%03d_ph%03d.dcm', filePrefix, slice, frame);
            % write multi-frame image to dicom file
            dicomwrite(uint16(round(abs(reshape(permute(image,[1 2 4 3]),Ny,Nx,1,[])))), fileName, tempinfo, 'CreateMode', 'Copy'); 
        end       
    end
% else                            % export 4D images
%     for slice = 1:Nz
%         % current slice
%         tempinfo.ImagePositionPatient = round(ImagePositionPatient1 + vecNorm*(slice-1), 4);
%         tempinfo.SliceLocation        = round(SliceLocation1 + sliceThickness*(slice-1), 5);
%         tempinfo.NumberOfFrames       = Nt*Nz;
%         tempinfo.NumberOfTemporalPositions = Nt;
%         for frame = 1:Nt        % time stamp for each frame
%             itemIdx = Nt*(slice-1) + frame;
%             tempstr = ['tempinfo.PerFrameFunctionalGroupsSequence.Item_' num2str(itemIdx) '.PlanePositionSequence.Item_1.ImagePositionPatient = round(ImagePositionPatient1 + vecNorm*(' num2str(slice) '-1), 4);'];
%             eval(tempstr);
%             tempstr = ['tempinfo.PerFrameFunctionalGroupsSequence.Item_' num2str(itemIdx) '.CardiacTriggerSequence.Item_1.CardiacTriggerDelayTime = (' num2str(frame) '-1)*' num2str(round(tempResol)) ';'];
%             eval(tempstr);
%         end
%         % file name for current slice
%         fileName = sprintf('%s_4D.dcm', filePrefix);
%         % write multi-frame image to dicom file
%         dicomwrite(uint16(round(abs(reshape(permute(image,[1 2 4 3]),Ny,Nx,1,[])))), fileName, tempinfo, 'CreateMode', 'Copy');        
%     end
end
warning('on','all');
disp('done.')

end


%% ========================================================================
%  pad or crop array to make size(img,dim) = N
function img = padcrop(img, dim, N)
if length(size(img)) >= dim
    img_dim_order = 1:length(size(img)); 
    img_dim_order(dim) = []; 
    img_dim_order = [dim img_dim_order];
    img = permute(img, img_dim_order);

    Norig = size(img, 1);
% skip padding
%     if Norig < N    % pad array
%         Npad_pre = floor(N/2) - floor(Norig/2);
%         Npad_post = N - Norig - Npad_pre;
%         img = padarray(img, Npad_pre, 'pre');
%         img = padarray(img, Npad_post,'post');
%     end
    if Norig > N    % crop array
        crop_start = floor(Norig/2) - floor(N/2);
        img = img(crop_start+(1:N),:,:,:);
    end

    img = ipermute(img, img_dim_order);
end
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

