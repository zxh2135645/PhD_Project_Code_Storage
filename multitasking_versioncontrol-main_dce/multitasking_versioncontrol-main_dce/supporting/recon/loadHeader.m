function [dataParams, reconOptions, twixObj, dataArray, temporalBasis, spatialCoeff, fitParams, status, msg] = loadHeader(dataParams, reconOptions)

if ~isfield(reconOptions,'flagCommandLine')
    reconOptions.flagCommandLine = true;
end
status = 0;
msg = [];
dataArray = [];

fprintf('Loading %s\n', dataParams.fidString);
if strcmp(dataParams.fidString((end-3):end),'.dat')    % Load Siemens rawdata
    try
        %% Load raw data header
        twixObj = mapVBVD(dataParams.fidString);
        if iscell(twixObj)
            if numel(twixObj) > 1
                if isfield(twixObj{1},'noise')
                    dataArray.noiseData = twixObj{1}.noise(); 
                end
                twixObj = twixObj{end};
            else
                twixObj = twixObj{1};
            end
        end
        
        if ~isfield(twixObj,'image') && isfield(twixObj,'RTfeedback')
            twixObj.image = twixObj.RTfeedback;
            twixObj = rmfield(twixObj,'RTfeedback');
        end
        
        % Manual anonymization (rawdata should already be anonymized, but just in case)
        twixObj.hdr.Config.tPatientName = 'xxxxxx';
        twixObj.hdr.Dicom.tPatientName  = 'xxxxxx';
        twixObj.hdr.Meas.tPatientName   = 'xxxxxx';
        twixObj.hdr.Config.PatientID = 'xxxxxxxxxxxxxxxx';
        twixObj.hdr.Meas.PatientID   = 'xxxxxxxxxxxxxxxx';
        twixObj.hdr.Config.PatientBirthDay = 'xxxxxxxx';
        
        % Parse header
        dataParams.MID = twixObj.hdr.Meas.MeasUID;
        dataParams.VerString = twixObj.image.softwareVersion;
        dataParams.lResonanceFrequency = twixObj.hdr.MeasYaps.sTXSPEC.asNucleusInfo{1}.lFrequency;  % Hz
        dataParams.alLarmorConstant = twixObj.hdr.Meas.alLarmorConstant(1);                         % Hz/T
        dataParams.MagneticFieldStrength = twixObj.hdr.Dicom.flMagneticFieldStrength;
        dataParams.SequenceName = twixObj.hdr.Config.SequenceString;
        dataParams.MRAcquisitionType = twixObj.hdr.Dicom.tMRAcquisitionType;
        dataParams.BodyPartExamined = twixObj.hdr.Dicom.tBodyPartExamined;
        switch dataParams.VerString
            case 'vd'   % also VE/XA
                paramsMeas = twixObj.hdr.Config;
                dataParams.isVE              = true;
                dataParams.Necho             = paramsMeas.NEco;
                dataParams.alTR_seconds      = paramsMeas.TR(1)*1e-6;
                dataParams.alTE_seconds      = [twixObj.hdr.MeasYaps.alTE{1:dataParams.Necho}]*1e-6;
                dataParams.adFlipAngleDegree = twixObj.hdr.MeasYaps.adFlipAngleDegree{1};
                try
                    dataParams.dPhaseFOV_mm = paramsMeas.PhaseFOV;
                catch
                    dataParams.dPhaseFOV_mm = paramsMeas.PhaseFoV; 
                end
                try
                    dataParams.dReadoutFOV_mm = paramsMeas.ReadoutFOV;
                catch
                    dataParams.dReadoutFOV_mm = paramsMeas.ReadFoV;
                end
                dataParams.lSegments       = twixObj.hdr.MeasYaps.sFastImaging.lSegments;
                dataParams.lEchoSpacing    = twixObj.hdr.Meas.lEchoSpacing*1e-6;   
                dataParams.lBaseResolution = paramsMeas.BaseResolution;
                dataParams.lPartitions     = paramsMeas.NPar;
                
                % Trajectory type
                if twixObj.hdr.MeasYaps.sKSpace.ucTrajectory == 1
                    dataParams.Trajectory  = 'Cartesian';
                    dataParams.isCartesian = true;
                elseif twixObj.hdr.MeasYaps.sKSpace.ucTrajectory == 2
                    dataParams.Trajectory  = 'Radial';
                    dataParams.isCartesian = false;  
                else%if twixObj.hdr.MeasYaps.sKSpace.ucTrajectory == 4
                    dataParams.Trajectory  = 'Spiral';
                    dataParams.isCartesian = false;  
                end   
            case 'vb'
                paramsMeas = twixObj.hdr.Meas;
                dataParams.isVE            = false;
                dataParams.Necho           = paramsMeas.NEco;
                dataParams.alTR_seconds    = paramsMeas.alTR(1)*1e-6;
                dataParams.alTE_seconds    = paramsMeas.alTE(1)*1e-6;
                dataParams.dPhaseFOV_mm    = paramsMeas.PhaseFoV;
                dataParams.dReadoutFOV_mm  = paramsMeas.ReadFoV;
                dataParams.lEchoSpacing    = paramsMeas.lEchoSpacing*1e-6;
                dataParams.lPartitions     = paramsMeas.NPar;
                dataParams.lSegments       = paramsMeas.lSegments;
                dataParams.lBaseResolution = paramsMeas.BaseResolution;
                 
                % Trajectory type
                if ischar(twixObj.hdr.MeasYaps.sKSpace.ucTrajectory)
                    if strcmp(twixObj.hdr.MeasYaps.sKSpace.ucTrajectory, '0x1')
                        dataParams.Trajectory  = 'Cartesian';
                        dataParams.isCartesian = true;
                    elseif strcmp(twixObj.hdr.MeasYaps.sKSpace.ucTrajectory, '0x2')
                        dataParams.Trajectory  = 'Radial';
                        dataParams.isCartesian = false;  
                    else %if strcmp(twixObj.hdr.MeasYaps.sKSpace.ucTrajectory, '0x3')
                        dataParams.Trajectory  = 'Spiral';
                        dataParams.isCartesian = false;  
                    end
                else
                    if twixObj.hdr.MeasYaps.sKSpace.ucTrajectory == 1
                        dataParams.Trajectory  = 'Cartesian';
                        dataParams.isCartesian = true;
                    elseif twixObj.hdr.MeasYaps.sKSpace.ucTrajectory == 2
                        dataParams.Trajectory  = 'Radial';
                        dataParams.isCartesian = false;  
                    elseif twixObj.hdr.MeasYaps.sKSpace.ucTrajectory == 3
                        dataParams.Trajectory  = 'Spiral';
                        dataParams.isCartesian = false;  
                    end   
                end
            otherwise
                fprintf(2,'Version error!\n');
        end
        dataParams.scanParams = paramsMeas;
           
        % Dicom info struct
        dataParams.dcminfo = genDicomInfo(twixObj);
        
        % Multitasking parameters
        dataParams = getMTfreeparams(twixObj, dataParams);
        
        if ~isfield(dataParams,'linesPerShot') || dataParams.linesPerShot == 0
            dataParams.linesPerShot = dataParams.lSegments;
        else
            dataParams.lSegments = dataParams.linesPerShot;
        end
        
        % Readout mode, only valid in multi-echo scans. 1: monopolar; 2: bipolar
        try
            if twixObj.image.NEco > 1
                dataParams.ReadoutMode = twixObj.hdr.Meas.ucReadOutMode;
            else
                dataParams.ReadoutMode = 1;
            end
        catch
            dataParams.ReadoutMode = 1;
        end
        
        % Geometry parameters
        try 
            dataParams.Nzorig = paramsMeas.NImagePar;
        catch
            try
                dataParams.Nzorig = paramsMeas.NoImagesPerSlab;
            catch
                dataParams.Nzorig = paramsMeas.NPar;
            end
        end
        dataParams.Nz = max(twixObj.image.Par);

        % if SMS, reset partition encoding array and Nz
        if dataParams.MBfactor > 1 || dataParams.Nz < dataParams.Nzorig
            dataParams.Nzdisp = dataParams.Nz;
        else
            dataParams.Nzdisp = dataParams.Nzorig;
        end
      
        % Slice position: coordinate of the center of the slice/slab
        dataParams.centerPosition = twixObj.image.slicePos(1:3,1);
        if ~isempty(twixObj.hdr.Meas.lGlobalTablePosTra)
            dataParams.centerPosition(3) = dataParams.centerPosition(3) + twixObj.hdr.Meas.lGlobalTablePosTra;
        end

        if strcmp(twixObj.hdr.Dicom.tMRAcquisitionType, '2D')
            dataParams.dThickness_mm = twixObj.hdr.MeasYaps.sSliceArray.asSlice{1,1}.dThickness;
            dataParams.dSlabFOV_mm   = dataParams.dThickness_mm * dataParams.Nzorig;
            if isfield(twixObj.hdr.MeasYaps.sGroupArray.asGroup{1},'dDistFact')
                dataParams.vecNormScale = (1 + twixObj.hdr.MeasYaps.sGroupArray.asGroup{1}.dDistFact) * dataParams.dThickness_mm;
            else
                dataParams.vecNormScale = dataParams.dThickness_mm;
            end
        else
            dataParams.dSlabFOV_mm   = twixObj.hdr.MeasYaps.sSliceArray.asSlice{1,1}.dThickness;
            dataParams.dThickness_mm = dataParams.dSlabFOV_mm/dataParams.Nzorig;
            dataParams.vecNormScale  = dataParams.dThickness_mm;
        end
        dataParams.volumeFOV    = [dataParams.dPhaseFOV_mm; dataParams.dReadoutFOV_mm; dataParams.dSlabFOV_mm];
        try
            NyImg = twixObj.hdr.Config.NImageLins;
        catch
            try 
                NyImg = twixObj.hdr.Config.imageNy;
            catch
                NyImg = twixObj.hdr.Config.N0FImageLines;
            end
        end
        try
            NyMeas = twixObj.hdr.Config.NoOfFourierLines;
        catch
            try
                NyMeas = twixObj.hdr.Config.NLinMeas;
            catch
                NyMeas = round(NyImg*twixObj.hdr.Meas.dPhaseResolution);
            end
        end
        if strcmp(dataParams.Trajectory,'Cartesian')
            dataParams.NImageLines = NyImg;
            dataParams.NyMeas = NyMeas;
            dataParams.voxelSpacing = [dataParams.dPhaseFOV_mm/NyImg;dataParams.dReadoutFOV_mm/dataParams.lBaseResolution;dataParams.dThickness_mm]; 
            dataParams.rawVoxelSpacing = [dataParams.dPhaseFOV_mm/NyMeas;dataParams.dReadoutFOV_mm/dataParams.lBaseResolution;dataParams.dThickness_mm]; 
        else
            dataParams.NImageLines = dataParams.lBaseResolution;
            dataParams.NyMeas = dataParams.lBaseResolution;
            dataParams.voxelSpacing = [dataParams.dReadoutFOV_mm/dataParams.lBaseResolution;dataParams.dReadoutFOV_mm/dataParams.lBaseResolution;dataParams.dThickness_mm]; 
            dataParams.rawVoxelSpacing = dataParams.voxelSpacing;
        end
        rawVoxelSpacing = dataParams.rawVoxelSpacing;
        [~,minSpacing] = min(rawVoxelSpacing);
        temp = 1:3; temp(minSpacing) = [];
        tempIdx1 = temp(1);
        tempIdx2 = temp(2);
        newRatio1 = rawVoxelSpacing(temp(1))/rawVoxelSpacing(minSpacing);
        newRatio2 = rawVoxelSpacing(temp(2))/rawVoxelSpacing(minSpacing);
        [~,Idx] = sort([tempIdx1 tempIdx2 minSpacing]);
        newRatio = [newRatio1 newRatio2 1];
        dataParams.dispRatio = newRatio(Idx(1:2));

        for n = 1:dataParams.Necho
            dataParams.pixelBandwidth(n) = round(1e9 / (twixObj.hdr.MeasYaps.sRXSPEC.alDwellTime{n}*dataParams.lBaseResolution*twixObj.hdr.Meas.flReadoutOSFactor));
        end
        
        % Find slice orientation for image display: 
        % (x,y,z) is in LPS patient coordinate
        RotMatQ = quat2RotMat(twixObj.image.slicePos(4:7,1));
        dataParams.vecCol  = RotMatQ*[1;0;0];
        dataParams.vecRow  = RotMatQ*[0;1;0];
        dataParams.vecNorm = RotMatQ*[0;0;1];
        vec3 = [dataParams.vecCol dataParams.vecRow dataParams.vecNorm];
        [~,zDir] = max(abs(vec3(3,:)));
        vec2 = vec3; vec2(:,zDir) = 0;
        [~,xDir] = max(abs(vec2(1,:)));
        vec2(:,xDir) = 0;
        [~,yDir] = max(abs(vec2(2,:)));
        dataParams.xDir = xDir*sign(vec3(1,xDir));
        dataParams.yDir = yDir*sign(vec3(2,yDir));
        dataParams.zDir = zDir*sign(vec3(3,zDir));
        dataParams.RotMatQ = RotMatQ;

        % Data processing flags
        dataArray.flagIsTrajectorySet  = false;
        dataArray.flagIsPrewhitened    = false;
        dataArray.flagIsCompressed     = false;
        dataArray.flagIsDriftCorrected = false;
        dataArray.flagIsBinned         = false;
        dataArray.flagTensor           = false;
        dataArray.flagWavelet          = false;

        % Initialize other return struct variables
        temporalBasis = [];
        spatialCoeff  = [];
        fitParams     = [];
    
        if reconOptions.flagCommandLine
            fprintf(' Header loaded. ScanType: %s. Check scan parameters.\n',dataParams.ScanType);
        else
            msg = sprintf(' Header loaded. ScanType: %s. Check scan parameters.\n',dataParams.ScanType);
        end
        status = 0;
    catch errormsg
        if reconOptions.flagCommandLine
            fprintf(2,' Load Siemens rawdata failed.\n');
            fprintf(2,'%s\n', errormsg.message);
        end
        msg = errormsg.message;
        status = 2;
    end
elseif strcmp(dataParams.fidString((end-3):end),'.mat')    % Load saved Matlab workspace
    try
        load(dataParams.fidString);
        if exist('Params','var')
            dataParams = Params;
        else
            dataParams = [];
        end
        if exist('ReconOptions','var')
            reconOptions = ReconOptions;
        else
            reconOptions = [];
        end
        if exist('TwixObj','var')
            twixObj = TwixObj;
        else
            twixObj = [];
        end
        if exist('DataArray','var')
            dataArray = DataArray;
        else
            dataArray = [];
        end
        if exist('TemporalBasis','var')
            temporalBasis = TemporalBasis;
        else
            temporalBasis = [];
        end
        if exist('SpatialCoeff','var')
            spatialCoeff = SpatialCoeff;
        else
            spatialCoeff = [];
        end
        if exist('FitParams','var')
            fitParams = FitParams;
        else
            fitParams = [];
        end
        if reconOptions.flagCommandLine
            fprintf(' Data loaded.\n');
        else
            msg = ' Data loaded.';
        end
    catch errormsg
        if reconOptions.flagCommandLine
            fprintf(2,' Not a valid multitasking recon workspace!\n');
            fprintf(2,'%s\n', errormsg.message);
        end
        msg = errormsg.message;
        status = 2;
    end
else  
    if reconOptions.flagCommandLine
        fprintf(2, ' Wrong file format!\n');
    end
    msg = ' Wrong file format';
    status = 2;
end


%% convert quaternions to a 3x3 rotation matrix
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
