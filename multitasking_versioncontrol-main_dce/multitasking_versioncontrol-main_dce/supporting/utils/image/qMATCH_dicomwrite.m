function exportFolder = qMATCH_dicomwrite(params,reconOptions,dataArray,temporalBasis,spatialCoeff,fitResult)

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(temporalBasis);
extractVarFromStruct(fitResult.fitParams);

L = size(Gr,1);
Necho = size(Phi,6);
tempNecho = size(Phi,1)/L;
slice = 1:Nzdisp;
useL = 1:L;
eIdx = 1;
rIdx = 1;
cIdx = 1;
tIdx = linesPerShot;

echo = 1;
Phi   = Phi(min(eIdx,tempNecho):tempNecho:end,:,:,:,:,eIdx);
if size(Phi,2) ~= Nseg*moduleLength
    Phi = reshape(permute(Phi,[1 2 5 3 4]),L,Nseg*moduleLength,size(Phi,3),[]);
end

Nzdisp = numel(fitSlice);
[Ny,Nx,Nz,~] = size(U);
dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), :, :, :);
Utemp = reshape(dispim(U),Nydisp,Nxdisp,Nzdisp,tempNecho,L);
Utemp = reshape(Utemp(:,:,:,min(eIdx,tempNecho),:),[],L);

reconVolume = reshape(Utemp*Phi(:,Nseg),Nydisp,Nxdisp,Nzdisp,1);
reconVolume = real(reconVolume.*exp(-1j*angle(reconVolume(:,:,:,end))));
reconVolume = reconVolume/prctile(abs(reconVolume(:)),99.7);
reconVolume(reconVolume>1) = 1;

maskBackgoundRough = genMask3D(reconVolume,5);
maskBackgoundRough = imgaussfilt3(maskBackgoundRough+0.1,2);

% Select blood/muscle voxels and extract signal curves
voxelBlood  = volumeViewer3Axes(reconVolume,'Select BLOOD voxel near the bifurcation');
voxelMuscle = volumeViewer3Axes(reconVolume,'Select MUSCLE voxel');

maskBlood = zeros(size(reconVolume));
maskBlood(voxelBlood(1),voxelBlood(2),voxelBlood(3)) = 1;
maskBlood = imdilate(maskBlood,strel("sphere",1));

maskMuscle = zeros(size(reconVolume));
maskMuscle(voxelMuscle(1),voxelMuscle(2),voxelMuscle(3)) = 1;
maskMuscle = imdilate(maskMuscle,strel("sphere",1));

PhiRecovery = Gr\reshape(Phi(eIdx:Necho:end,:,cIdx,rIdx,1,eIdx),L,[]);

Ublood = reshape(Utemp(logical(maskBlood(:)),:),[],L);
signalBlood = mean(realify(Ublood*PhiRecovery,'rows'),1);
Umuscle   = reshape(Utemp(logical(maskMuscle(:)),:),[],L);
signalMuscle = mean(realify(Umuscle*PhiRecovery,'rows'),1);

% Generate dark blood image
[~,dbIdx] = min(abs(signalBlood)./abs(signalMuscle));

frameRange = 15;
dbIdxrange = dbIdx + (-frameRange:frameRange);
dbIdxrange(dbIdxrange<1) = [];
dbIdxrange(dbIdxrange>size(Phi,2)) = [];
reconDBtemp = abs(reshape(Utemp*(Gr\reshape(Phi(eIdx:Necho:end,dbIdxrange,cIdx,rIdx,1,eIdx),L,[])),Nydisp,Nxdisp,Nzdisp,[]));
reconDBtemp2 = zeros(Nydisp,Nxdisp,Nzdisp,numel(dbIdxrange)-4);
for frame = 1:size(reconDBtemp2,4)
    reconDBtemp2(:,:,:,frame) = min(reconDBtemp(:,:,:,frame:frame+4),[],4);
end
reconDBtemp2 = reconDBtemp2/prctile(reconDBtemp2(:),99.9);
reconDBtemp2(reconDBtemp2>1) = 1;

sliceIdx = voxelBlood(3)+(-4:4);
sliceIdx(sliceIdx<1) = [];
sliceIdx(sliceIdx>Nzdisp) = [];
selectedFrame = montageSliceViewer(imageOrientLPS(reconDBtemp2(:,:,sliceIdx,:),params),'Use slider to select dark-blood frame');
dbIdx = dbIdx - frameRange + selectedFrame;

reconDB = reconDBtemp2(:,:,:,selectedFrame);
reconDB = reconDB/prctile(reconDB(:),99.9);
reconDB(reconDB>1) = 1;

% Generate bright blood MRA image
[~,bbIdx] = min(abs(signalMuscle)./abs(signalBlood));

% if sign(signalMuscle(bbIdx)*signalBlood(bbIdx)) < 0
%     if signalBlood(bbIdx) < 0
%         bbIdx = bbIdx - 1;
%     else
%         bbIdx = bbIdx + 1;
%     end
% end
% 
% temp1 = mean(signalMuscle(bbIdx:bbIdx+1));
% temp2 = mean(signalMuscle(bbIdx-1:bbIdx));
% if abs(temp1) < abs(signalMuscle(bbIdx)) && abs(temp1) < abs(temp2)
%     bbIdx = bbIdx:bbIdx+1;
% elseif abs(temp2) < abs(signalMuscle(bbIdx)) && abs(temp2) < abs(temp1)
%     bbIdx = bbIdx-1:bbIdx;
% end

bbIdx = bbIdx + (-3:3);

reconMRA = abs(reshape(Utemp*(Gr\reshape(Phi(eIdx:Necho:end,bbIdx,cIdx,rIdx,1,eIdx),L,[])),Nydisp,Nxdisp,Nzdisp,[]));
reconMRA = min(reconMRA,[],4);
reconMRA = reconMRA/prctile(reconMRA(:),99.9);
reconMRA(reconMRA>1) = 1;
reconMRA = reconMRA.^2.*(1-reconDB).^2.*maskBackgoundRough;
reconMRA = reconMRA/prctile(reconMRA(:),99.9);
reconMRA(reconMRA>1) = 1;

maskVessel = reconMRA.^2;
maskVessel(maskVessel<0.2) = 0.2;
maskVessel = imfill(maskVessel,4);
reconDB = reconDB./(maskVessel*5).*maskBackgoundRough;
reconDB = reconDB/prctile(reconDB(:),99.9);
reconDB(reconDB>1) = 1;

% T1/T2 maps
T1map = fitResult.T1map;
T2map = fitResult.T2map;
BIRmap = fitResult.BIRmap;
B1map  = fitResult.B1map;

% Generate weighted images
[~,tIdx] = min(fitResult.fitParams.TEs);
tIdx = linesPerShot*tIdx;
reconMixed = abs(mean(reshape(Utemp*(Gr\reshape(mean(Phi(eIdx:Necho:end,tIdx,cIdx,rIdx,1,eIdx),2),L,[])),Nydisp,Nxdisp,Nzdisp,[]),4));
reconMixed = reconMixed./prctile(reconMixed(:),99.9);

% Steady-state image weighting
alpha = flipAngleArray(1)*pi/180;
e1 = exp(-lEchoSpacing./T1map);
Mss = (1-e1) ./ (1-cos(B1map).*e1);
ssWeight = abs(Mss.*(1 - (BIRmap+1).*(e1.*cos(B1map)).^(linesPerShot-1)) .* sin(B1map));
ssWeight = abs((1-e1) ./ (1-cos(alpha).*e1).*(1 - 2.*(e1.*cos(B1map)).^(linesPerShot-1)));
ssWeight = ssWeight./prctile(ssWeight(:),99.9);
ssWeight(ssWeight==0) = 1;

% T1maptemp = T1map;
% T1maptemp(T1map==0) = 0.1;
% T1weighting = abs(1-2*exp(-0.5./T1maptemp));
% reconT1w = reconMixed.*imgaussfilt(T1weighting./ssWeight,2);
% reconT1w = reconT1w/prctile(reconT1w(:),99.9);
reconT1w = reconMixed;

% T2-weighted image = steady-state image .* T2-weighting ./ steady-state weighting
T2maptemp = T2map;
T2maptemp(T2map==0) = 0.01;
T2weighting = exp(-0.04./T2maptemp);   % 40ms T2 relaxation
reconT2w = reconMixed.*imgaussfilt(T2weighting./ssWeight,2);
reconT2w = reconT2w/prctile(reconT2w(:),99.9);

% % Generate T1-weighted image
% reconT1w = abs(reshape(Utemp*(Gr\reshape(mean(Phi(eIdx:Necho:end,end,cIdx,rIdx,1,eIdx),2),L,[])),Nydisp,Nxdisp,Nzdisp));
% 
% % Generate T2-weighted image
% [~,tIdx] = max(fitResult.fitParams.TEs);
% tIdx = linesPerShot * (tIdx-1) + 1;
% reconT2w = abs(reshape(Utemp*(Gr\reshape(mean(Phi(eIdx:Necho:end,tIdx,cIdx,rIdx,1,eIdx),2),L,[])),Nydisp,Nxdisp,Nzdisp));


% Sequence-specific info
% Spatial resolution
% pixelSpacing = [distance between rows, distance between columns] in mm
pixelSpacing = voxelSpacing(1:2);

% Slice normal vector length = slice thickness + gap            (2D)
%                            = slab thickness / # of partitions (3D)
vecNorm = vecNorm*vecNormScale;

% Rotate and flip the image if necessary
[~, Mdir] = max(abs(vecNorm));
[~, Midx] = max(abs([vecRow vecCol]), [], 2);

if Mdir == 1        % sagittal plane
    if Midx(3) == 1     % align the column vector to S/I direction
        vecCol = RotMatQ*[0;1;0];
        vecRow = RotMatQ*[1;0;0];
        pixelSpacing = flip(pixelSpacing);
        reconDB  = permute(reconDB,[2 1 3 4]);
        reconMRA = permute(reconMRA,[2 1 3 4]);
        reconT1w = permute(reconT1w,[2 1 3 4]);
        reconT2w = permute(reconT2w,[2 1 3 4]);
        T1map = permute(T1map,[2 1 3 4]);
        T2map = permute(T2map,[2 1 3 4]);
    end
    if vecRow(2) < 0    % align image left/right to A/P
        vecRow = -vecRow;
        reconDB  = flip(reconDB,2);
        reconMRA = flip(reconMRA,2);
        reconT1w = flip(reconT1w,2);
        reconT2w = flip(reconT2w,2);
        T1map = flip(T1map,2);
        T2map = flip(T2map,2);
    end
    if vecCol(3) > 0    % align image top/bottom to S/I
        vecCol = -vecCol;
        reconDB  = flip(reconDB,1);
        reconMRA = flip(reconMRA,1);
        reconT1w = flip(reconT1w,1);
        reconT2w = flip(reconT2w,1);
        T1map = flip(T1map,1);
        T2map = flip(T2map,1);
    end 
elseif Mdir == 2    % coronal plane
    if Midx(3) == 1     % align the column vector to S/I direction
        vecCol  = RotMatQ*[0;1;0];
        vecRow  = RotMatQ*[1;0;0];
        pixelSpacing = flip(pixelSpacing);
        reconDB  = permute(reconDB,[2 1 3 4]);
        reconMRA = permute(reconMRA,[2 1 3 4]);
        reconT1w = permute(reconT1w,[2 1 3 4]);
        reconT2w = permute(reconT2w,[2 1 3 4]);
        T1map = permute(T1map,[2 1 3 4]);
        T2map = permute(T2map,[2 1 3 4]);
    end
    if vecRow(1) < 0    % align image left/right to R/L
        vecRow = -vecRow;
        reconDB  = flip(reconDB,2);
        reconMRA = flip(reconMRA,2);
        reconT1w = flip(reconT1w,2);
        reconT2w = flip(reconT2w,2);
        T1map = flip(T1map,2);
        T2map = flip(T2map,2);
    end
    if vecCol(3) > 0    % align image top/bottom to S/I
        vecCol = -vecCol;
        reconDB  = flip(reconDB,1);
        reconMRA = flip(reconMRA,1);
        reconT1w = flip(reconT1w,1);
        reconT2w = flip(reconT2w,1);
        T1map = flip(T1map,1);
        T2map = flip(T2map,1);
    end 
elseif Mdir == 3    % transverse plane
    if Midx(2) == 1     % align the column vector to A/P direction
        vecCol  = RotMatQ*[0;1;0];
        vecRow  = RotMatQ*[1;0;0];
        pixelSpacing = flip(pixelSpacing);
        reconDB  = permute(reconDB,[2 1 3 4]);
        reconMRA = permute(reconMRA,[2 1 3 4]);
        reconT1w = permute(reconT1w,[2 1 3 4]);
        reconT2w = permute(reconT2w,[2 1 3 4]);
        T1map = permute(T1map,[2 1 3 4]);
        T2map = permute(T2map,[2 1 3 4]); 
    end
    if vecRow(1) < 0    % align image left/right to R/L
        vecRow = -vecRow;
        reconDB  = flip(reconDB,2);
        reconMRA = flip(reconMRA,2);
        reconT1w = flip(reconT1w,2);
        reconT2w = flip(reconT2w,2);
        T1map = flip(T1map,2);
        T2map = flip(T2map,2);
    end
    if vecCol(2) < 0    % align image top/bottom to A/P
        vecCol = -vecCol;
        reconDB  = flip(reconDB,1);
        reconMRA = flip(reconMRA,1);
        reconT1w = flip(reconT1w,1);
        reconT2w = flip(reconT2w,1);
        T1map = flip(T1map,1);
        T2map = flip(T2map,1);
    end 
end

% Convert image format
reconDB  = uint16(round(4095*abs(reconDB)));    % normalize to uint16
reconMRA = uint16(round(4095*abs(reconMRA)));
reconT1w = uint16(round(4095*abs(reconT1w)));
reconT2w = uint16(round(4095*abs(reconT2w)));
T1map    = uint16(round(1000*abs(T1map)));      % convert sec to msec
T2map    = uint16(round(10*1000*abs(T2map)));   % convert sec to msec*10


%% dicom info
if isunix
    sep = '/';
else
    sep = '\';
end

% Get scan orientation
[~,normIdx] = max(abs(vecNorm));
if normIdx == 3
    orientation = 'TRA';
elseif normIdx == 2
    orientation = 'COR';
else
    orientation = 'SAG';
end

% Get case name
% caseNameIdx = strfind(lower(fidString),'dtect');
% if isempty(caseNameIdx)
%     temp = fileparts(lower(fidString));
%     tempIdx = strfind(temp,sep);
%     if isempty(tempIdx) 
%         caseName = 'DTECT-XX-YYY-ZZZ';
%     else
%         tempIdx = tempIdx(end) + 1;
%         caseName = ['DTECT-' upper(temp(tempIdx:end))];
%     end
% else
%     caseNameEndIdx = regexpi(fidString(caseNameIdx:end),'[_/\\]');
%     if isempty(caseNameEndIdx)
%         caseName = upper(fidString(caseNameIdx:end));
%     else
%         caseName = upper(fidString(caseNameIdx+(0:caseNameEndIdx(1)-2)));
%     end
% end
caseName = 'DTECT-';
prompt = {'Enter case name: (format: DTECT-XX-YYY-ZZZ)'};
fieldsize = [1 45];
definput = {caseName};
caseName = inputdlg(prompt,'Input',fieldsize,definput);
caseName = upper(caseName{1});

% Choose dicom template and read dicom info
currentPath = pwd;
dicomPath = uigetdir(path,'Choose dicom template folder'); 
cd(dicomPath);
dcm_dir = dir(strcat('*'));
if dcm_dir(1).name == '.'
    dcm_dir = dcm_dir(3:end);
end
dcm_info_template = [];
for i = 1:numel(dcm_dir)
    if isdicom(dcm_dir(i).name) 
        dcm_info_template = dicominfo(dcm_dir(i).name);
        break;
    end
end
cd(currentPath);
if isempty(dcm_info_template)
    fprintf(2,'No dicom template detected!\n');
end

% Create dicom export folder
exportPath = fileparts(params.fidString);
exportFolder = strcat('T1-T2-MAP-3D-reconed','_',caseName,'-',orientation);
mkdir([exportPath sep exportFolder]);

% Generate file prefix and create folders for each contrast image
datetimeStr = char(datetime('now','Format','yyyyMMdd''T''HHmmss'));

filePrefixDB  = [exportPath sep exportFolder sep 'dcm_vw_'  datetimeStr];
filePrefixMRA = [exportPath sep exportFolder sep 'dcm_mra_' datetimeStr];
filePrefixT1w = [exportPath sep exportFolder sep 'dcm_T1w_' datetimeStr];
filePrefixT2w = [exportPath sep exportFolder sep 'dcm_T2w_' datetimeStr];
filePrefixT1map = [exportPath sep exportFolder sep 'dcm_T1Map_' datetimeStr];
filePrefixT2map = [exportPath sep exportFolder sep 'dcm_T2Map_' datetimeStr];

mkdir(filePrefixDB);
mkdir(filePrefixMRA);
mkdir(filePrefixT1w);
mkdir(filePrefixT2w);
mkdir(filePrefixT1map);
mkdir(filePrefixT2map);

% The location information of the first slice
% ImagePositionPatient1: coordinates of the first (top-left) pixel
% SliceLocation1: perpendicular distance from isocenter to the slice plane
% vecPos1: coordinate of the center of the slice
[Nydisp, Nxdisp, Nzdisp, ~] = size(reconDB);
if strcmp(MRAcquisitionType, '2D')
    ImagePositionPatient1 = centerPosition - vecCol*pixelSpacing(1)*Nydisp/2 - vecRow*pixelSpacing(2)*Nxdisp/2;
    SliceLocation1 = dot(centerPosition, vecNorm/vecNormScale);
else
    ImagePositionPatient1 = centerPosition - vecCol*pixelSpacing(1)*Nydisp/2 - vecRow*pixelSpacing(2)*Nxdisp/2 - vecNorm*(Nzdisp-1)/2;
    SliceLocation1 = dot(centerPosition - vecNorm*(Nzdisp-1)/2, vecNorm/vecNormScale);
end

% Loop through slices
fprintf('Exporting dicom images to %s\n --- %d slices in total: ',exportFolder,Nzdisp)
for slice = 1:Nzdisp
    fprintf('%d. ',slice);
    dcminfo = dcm_info_template;

    % Protocol info
    dcminfo.SliceThickness          = dThickness_mm;
    dcminfo.PixelSpacing            = pixelSpacing;
    dcminfo.PixelBandwidth          = pixelBandwidth;
    dcminfo.ImageOrientationPatient = round([vecRow; vecCol], 6);
    dcminfo.ImagePositionPatient    = round(ImagePositionPatient1, 4);
    dcminfo.RepetitionTime          = alTR_seconds*1000;  % convert from sec to msec
    dcminfo.EchoTime                = alTE_seconds*1000;  % convert from sec to msec
    dcminfo.FlipAngle               = flipAngleArray(1);
    dcminfo.SequenceName            = SequenceName;
    dcminfo.MRAcquisitionType       = MRAcquisitionType;
    dcminfo.InstanceNumber          = slice;
    
    % Current slice location
    dcminfo.ImagePositionPatient = round(ImagePositionPatient1 + vecNorm*(slice-1), 4);
    dcminfo.SliceLocation        = round(SliceLocation1 + dThickness_mm*(slice-1), 5);

    % Save dark-blood images
    dcminfo.SeriesNumber = 2006; 
    dcminfo.ProtocolName = 'qMATCH_vw';
    dcminfo.SeriesDescription = 'qMATCH_vw';
    dcminfo.PatientID = caseName;
    dcminfo.StudyID = '1';
    dcminfo.WindowCenter = 2048;
    dcminfo.WindowWidth  = 4095;
    fileName = sprintf('%s%s%s%s%d%s%04d.dcm',filePrefixDB,sep,caseName,'.MR.qMATCH.vw.',dcminfo.SeriesNumber,'.',slice);
    dicomwrite(reconDB(:,:,slice), fileName, dcminfo); 

    % save bright-blood MRA images
    dcminfo.SeriesNumber = 2007;
    dcminfo.ProtocolName = 'qMATCH_MRA';
    dcminfo.SeriesDescription = 'qMATCH_MRA';
    dcminfo.PatientID = caseName;
    dcminfo.StudyID = '1';
    dcminfo.WindowCenter = 2048;
    dcminfo.WindowWidth  = 4095;
    fileName = sprintf('%s%s%s%s%d%s%04d.dcm',filePrefixMRA,sep,caseName,'.MR.qMATCH.mra.',dcminfo.SeriesNumber,'.',slice);
    dicomwrite(reconMRA(:,:,slice), fileName, dcminfo); 

    % Save T1-weighted images
    dcminfo.SeriesNumber = 2008;
    dcminfo.ProtocolName = 'qMATCH_T1W';
    dcminfo.SeriesDescription = 'qMATCH_T1W';
    dcminfo.PatientID = caseName;
    dcminfo.StudyID = '1';
    dcminfo.WindowCenter = 2048;
    dcminfo.WindowWidth  = 4095;
    fileName = sprintf('%s%s%s%s%d%s%04d.dcm',filePrefixT1w,sep,caseName,'.MR.qMATCH.T1w.',dcminfo.SeriesNumber,'.',slice);
    dicomwrite(reconT1w(:,:,slice), fileName, dcminfo); 

    % Save T2-weighted images
    dcminfo.SeriesNumber = 2009;
    dcminfo.ProtocolName = 'qMATCH_T2W';
    dcminfo.SeriesDescription = 'qMATCH_T2W';
    dcminfo.PatientID = caseName;
    dcminfo.StudyID = '1';
    dcminfo.WindowCenter = 2048;
    dcminfo.WindowWidth  = 4095;
    fileName = sprintf('%s%s%s%s%d%s%04d.dcm',filePrefixT2w,sep,caseName,'.MR.qMATCH.T2w.',dcminfo.SeriesNumber,'.',slice);
    dicomwrite(reconT2w(:,:,slice), fileName, dcminfo); 

    % Save T1 map
    dcminfo.SeriesNumber = 2010;
    dcminfo.ProtocolName = 'qMATCH_T1map';
    dcminfo.SeriesDescription = 'qMATCH_T1map';
    dcminfo.PatientID = caseName;
    dcminfo.StudyID = '1';
    dcminfo.WindowCenter = 1500; %floor(min(1000*fitResult.fitParams.maxT1,4096)/2);
    dcminfo.WindowWidth  = 3000; %min(1000*fitResult.fitParams.maxT1,4095);
    fileName = sprintf('%s%s%s%s%d%s%04d.dcm',filePrefixT1map,sep,caseName,'.MR.qMATCH.T1Map.',dcminfo.SeriesNumber,'.',slice);
    dicomwrite(T1map(:,:,slice), fileName, dcminfo); 

    % Save T2 map
    dcminfo.SeriesNumber = 2011;
    dcminfo.ProtocolName = 'qMATCH_T2map';
    dcminfo.SeriesDescription = 'qMATCH_T2map';
    dcminfo.PatientID = caseName;
    dcminfo.StudyID = '1';
    dcminfo.WindowCenter = 50;   %floor(min(1000*fitResult.fitParams.maxT2,4096)/2);
    dcminfo.WindowWidth  = 100;  %min(1000*fitResult.fitParams.maxT2,4095);
    fileName = sprintf('%s%s%s%s%d%s%04d.dcm',filePrefixT2map,sep,caseName,'.MR.qMATCH.T2Map.',dcminfo.SeriesNumber,'.',slice);
    dicomwrite(T2map(:,:,slice), fileName, dcminfo); 
end

% Compress all images
fprintf('\nCompressing all dicom files to %s ... (this might take a while) ',[exportPath sep exportFolder '.zip']);
zip([exportPath sep exportFolder '.zip'],exportFolder,exportPath);
fprintf('\nDone.\n ',slice);

