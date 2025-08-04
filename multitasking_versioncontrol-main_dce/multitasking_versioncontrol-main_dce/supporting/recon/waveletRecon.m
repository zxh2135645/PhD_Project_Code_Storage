function [reconOptions,spatialCoeff] = waveletRecon(params,reconOptions,dataArray,temporalBasis,spatialCoeff)

fprintf('Wavelet recon... \n');

if params.Nz == 1
    [reconOptions,spatialCoeff] = waveletRecon2D(params,reconOptions,dataArray,temporalBasis,spatialCoeff);
else
    if isfield(reconOptions.wavelet,'flagForce2D') && reconOptions.wavelet.flagForce2D
        [reconOptions,spatialCoeff] = waveletRecon2DSlices(params,reconOptions,dataArray,temporalBasis,spatialCoeff);
    else
        [reconOptions,spatialCoeff] = waveletRecon3D_aniso(params,reconOptions,dataArray,temporalBasis,spatialCoeff);
    end
end
