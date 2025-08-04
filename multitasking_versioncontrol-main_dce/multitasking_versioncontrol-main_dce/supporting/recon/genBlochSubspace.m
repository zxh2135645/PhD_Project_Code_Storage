function [temporalBasis,fitParams] = genBlochSubspace(params, fitParams, reconOptions, temporalBasis, cLin, dataArray)

extractVarFromStruct(params);

if ~isfield(reconOptions,'flagDataDriven')
    reconOptions.flagDataDriven = false;
end

fitParams.params = params;
fitParams.filePath = filePath;
fitParams.ScanType = ScanType;
fitParams.isCartesian = isCartesian;
fitParams.Nx = Nx;
fitParams.Ny = Ny;
fitParams.Nz = Nz;
fitParams.Nxdisp = Nxdisp;
fitParams.Nydisp = Nydisp;
fitParams.Norig  = Norig;
fitParams.MBfactor = MBfactor;
fitParams.moduleLength = moduleLength;
fitParams.Nseg = linesPerShot;
fitParams.TR   = lEchoSpacing;
fitParams.numFA          = numFA;
fitParams.flipAngleArray = flipAngleArray;
fitParams.numIRns    = numIRns;
fitParams.numT2prep  = numT2prep;
fitParams.T2prepDuration = T2prepDuration;
fitParams.numT1rhoPrep   = numT1rhoPrep;
fitParams.T1rhoDuration  = T1rhoDuration;
fitParams.lEchoSpacing   = lEchoSpacing;

reconOptions.flagUseBT2 = true;

if strcmp(ScanType,'CEST')
    cL = 9;
    [curvePhi, curvePhi_binning, temporalBasis, fitParams, curvePhi_alt] = genCurveSubspace_CEST(params, fitParams, reconOptions, dataArray, temporalBasis);
    curvePhi = curvePhi(:,1:cL);
    temporalBasis.curvePhi_binning = curvePhi_binning(:,1:cL);
    temporalBasis.curvePhi_alt = curvePhi_alt(:,1:cL);
elseif ~strcmp(ScanType,'Cine')
    switch ScanType
        case {'SR','IR'}
            cL = 10;
        case 'T1rhoT2IR'
            cL = 10;
        otherwise
            cL = 8;
    end
  
    if nargin > 4
        cL = cLin;
    end
  
    switch ScanType
        case {'SR','IR'}
            [curveU,~,fitParams] = genCurveSubspace_T1_1FAn(params, fitParams, reconOptions);
        case {'SR_VFA','IR_VFA'}
            [curveU,~,fitParams] = genCurveSubspace_T1_VFA(params, fitParams, reconOptions);
        otherwise
            if reconOptions.flagDataDriven
                [curveU,~,fitParams] = genCurveSubspace(params, fitParams, reconOptions);
                temporalBasis.curvePhi_binning = curveU(:,1:cL);
                [curveU,~,fitParams] = genCurveSubspace_T1_1FA(params, fitParams, reconOptions);
                
            else
                [curveU,~,fitParams] = genCurveSubspace(params, fitParams, reconOptions);
            end
    end
    
%     if strcmp(ScanType,'IR')
%         [curveU,~,fitParams] = genCurveSubspace_T1_1FA(params, fitParams, reconOptions);
%     elseif strcmp(ScanType,'T2IR_VFA') || strcmp(ScanType,'T2IR')
%         [curveU,~,fitParams] = genCurveSubspace_symT2(params, fitParams, reconOptions, reconOptions.flagUseBT2);
%     elseif strcmp(ScanType,'T1rho_VFA') || strcmp(ScanType,'T1rho')
%         [curveU,~,fitParams] = genCurveSubspace_symT1rho(params, fitParams, reconOptions, reconOptions.flagUseBT2);
%     else
%         [curveU,~,fitParams] = genCurveSubspace(params, fitParams, reconOptions, reconOptions.flagUseBT2);
%     end

    curvePhi = curveU(:,1:cL);
    clear curveU
else
    curvePhi = ones(linesPerShot,1)/sqrt(linesPerShot);
    cL = 1;
end

if reconOptions.flagCommandLine
    figure;plot(realify(curvePhi(:,1:min(size(curvePhi,2),5)),'cols'));axis([1 size(curvePhi,1) -max(abs(curvePhi(:)))*1.05 max(abs(curvePhi(:)))*1.05]);title('curvePhi');
end

temporalBasis.curvePhi = curvePhi;
temporalBasis.cL = cL;
