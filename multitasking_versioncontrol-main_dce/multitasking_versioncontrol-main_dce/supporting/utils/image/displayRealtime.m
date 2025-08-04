function reconstructedImages = displayRealtime(params,reconOptions,dataArray,temporalBasis,spatialCoeff,reconstructedImages,slices,useL,fps)

extractVarFromStruct(params);
extractVarFromStruct(reconOptions);
extractVarFromStruct(dataArray);
extractVarFromStruct(temporalBasis);
extractVarFromStruct(spatialCoeff);

if nargin < 9
    fps = 20;
end

if nargin < 8 || isempty(useL)
    useL = 1:L_init;
end

if nargin < 7 || isempty(slices)
    slices = 1 + [floor(Nz/4) floor(Nz/2) floor(Nz/4)*3];
end
if (slices(1) < 1 || slices(end) > Nz) && Nz > 3
    slices = [floor(Nz/4) floor(Nz/2) floor(Nz/4)*3];
elseif Nz <= 3
    slices = 1:Nz;
end    

ds  = round(1/(fps*lEchoSpacing));    % downsampling to get to 20fps
Nim = linesPerShot*ceil(15/(linesPerShot*lEchoSpacing)); % display the first 15 sec
Nim = floor(min(Ntpoint,Nim)/ds);
frameTime = ds*lEchoSpacing;

dispim = @(x) x(floor(Ny/2)-floor(Nydisp/2) + (1:Nydisp), floor(Nx/2)-floor(Nxdisp/2) + (1:Nxdisp), slices, :);

Utemp = reshape(U_init,Ny,Nx,Nz,[]);
if MBfactor == 1
    Utemp = fftshift(Utemp,3);
end
if params.isCartesian
    Utemp = fftshift(Utemp,1);
end
Utemp = dispim(Utemp);

Nd = [Nydisp Nxdisp numel(slices)];
recon = genImageTimeSeries(Utemp, Nd, Phi_rt_full_init, ds, Nim, useL);

recon_slice = recon(:,:,ceil(Nd(3)/2),:);
cw = prctile(abs(recon_slice(:)),99.9);
recon_slice = recon_slice/cw;

if reconOptions.flagCommandLine
    implayZoom(imageOrientLPS(abs(recon_slice),params),fps);
end

reconstructedImages.reconRealtime = recon/cw;
reconstructedImages.reconRealtimeFrameTime = frameTime;

try
    [~,filePrefix] = fileparts(params.filePath);
    saveGif(abs(recon_slice),filePath,[filePrefix '_reconRealtime.gif'],1/frameTime);
catch
    fprintf("Saving images to %s failed. Saving them in current folder.", filePath);
    saveGif(abs(recon_slice),'.','reconRealtime.gif',1/frameTime);
end

