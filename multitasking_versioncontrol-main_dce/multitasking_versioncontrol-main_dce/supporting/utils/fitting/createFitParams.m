function fitParams = createFitParams(params,reconOptions,dataArray,temporalBasis,spatialCoeff,fitParams)

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(dataArray);
extractVarFromStruct(temporalBasis);
extractVarFromStruct(spatialCoeff);

try
    fitParams.license.mainpath  = reconOptions.mainpath;
    fitParams.license.licenseID = reconOptions.licenseID;
catch
    fitParams.license.mainpath  = '';
    fitParams.license.licenseID = 'NONE';
end

fitParams.MID      = MID;
fitParams.filePath = filePath;
fitParams.ScanType = ScanType;
fitParams.isCartesian  = isCartesian;
fitParams.interpFactor = interpFactor;
fitParams.Nx = Nx;
fitParams.Ny = Ny;
fitParams.Nz = Nz;
fitParams.Nxdisp = Nxdisp;
fitParams.Nydisp = Nydisp;
fitParams.Nzdisp = Nzorig;
fitParams.Norig  = Norig;
fitParams.MBfactor = MBfactor;
fitParams.moduleLength = moduleLength;
fitParams.Nseg = linesPerShot;
fitParams.TR   = lEchoSpacing;
fitParams.TEarray = alTE_seconds;
fitParams.numFA          = numFA;
fitParams.flipAngleArray = flipAngleArray;
fitParams.numSRns        = numSRns;
fitParams.numIRns        = numIRns;
fitParams.numIRsel       = numIRsel;
fitParams.numT2prep      = numT2prep;
fitParams.T2prepDuration = T2prepDuration;
fitParams.numT1rhoPrep   = numT1rhoPrep;
fitParams.T1rhoDuration  = T1rhoDuration;
fitParams.lEchoSpacing   = lEchoSpacing;
fitParams.rbins = rbins;
fitParams.cbins = cbins;
rep_VFA    = moduleLength/numFA;
fitParams.alphaArray = repmat(flipAngleArray, 1, rep_VFA)*pi/180;
fitParams.totalTime = reconOptions.totalTime;

if strcmp(ScanType,'CEST')
    fitParams.CESTnPools = fitParams.ParamsCEST.nPools;
    fitParams.CESTSatFreqOffsetppmList = CESTSatFreqOffsetppmListUnique;
    fitParams.CESTNumRep = CESTNumRep;
    fitParams.CESTMetabolite = CESTMetabolite;
    fitParams.CESTSatFA = CESTSatFA;
    fitParams.CESTSatDuration = CESTSatDuration;
    fitParams.CESTSatSpoilerTotalDuration = CESTSatSpoilerTotalDuration;
end

if exist('fitw_full','var')
    fitParams.fitw_full = fitw_full;
end
if exist('Phi','var')
    fitParams.Phi = Phi;
    L = size(Phi,1);
end
if exist('Gr','var')
    fitParams.Gr  = Gr;
    L = size(Gr,1);
end

dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), floor(Nz/2)-floor(Nzorig/2) + (1:Nzorig), :, :);

if exist('U','var')
    Utemp = reshape(U,Ny,Nx,Nz,[]);
    if MBfactor == 1
        Utemp = fftshift(Utemp,3);
    end
    if params.isCartesian
        Utemp = fftshift(Utemp,1);
    end
    fitParams.U = dispim(Utemp);
    
    if exist('roi_weighting','var')
        fitParams.roi_weighting = dispim(roi_weighting);
        if ~exist('ROstart','var')
            ROstart = 1;
            ROend   = size(fitParams.roi_weighting,2);
        end
        %fitParams.U_ROI = roi_weighting(:,ROstart:ROend,floor(Nz/2)-floor(Nzorig/2) + (1:Nzorig)).*fitParams.U;
    end
end

fitParams.flagLargeROI = 0;
fitParams.params = params;
fitParams.reconOptions = reconOptions;