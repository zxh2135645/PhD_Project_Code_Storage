function image = imageOrientLPS(image, params, flagResize, flagOrient)
% ========================================================================
%
%   Flip and rotate image to match the LPS orientation used in dicom where
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
%       image   - flipped/rotated 3D/4D image 
%
%  - last updated 2023-04-27
%  - Hsu-Lei Lee @ BIRI, Cedars-Sinai Medical Center
%
% ========================================================================

warning('off','all');

if nargin < 4
    flagOrient = 1;
end
if nargin < 3
    flagResize = 1;
end

%% Resize image to isotropic pixels if newRatio is available
if flagResize
    try
        rawVoxelSpacing = params.rawVoxelSpacing;
        [~,minSpacing] = min(rawVoxelSpacing);
        temp = 1:3; temp(minSpacing) = [];
        tempIdx1 = temp(1);
        tempIdx2 = temp(2);
        newRatio1 = rawVoxelSpacing(temp(1))/rawVoxelSpacing(minSpacing);
        newRatio2 = rawVoxelSpacing(temp(2))/rawVoxelSpacing(minSpacing);
        [~,Idx] = sort([tempIdx1 tempIdx2 minSpacing]);
        newRatio = [newRatio1 newRatio2 1];
        newRatio = newRatio(Idx(1:2));
        size2 = [size(image,1),size(image,2)];
        image = imresize(image,newRatio.*size2);
    catch
        fprintf('Resizing information not available.\n');
    end
end
%% Calculate image position and dimensions

if flagOrient
    % Column/row/slice vectors
    try
        vecCol  = params.vecCol;
        vecRow  = params.vecRow;
        vecNorm = params.vecNorm;
    catch
        vecCol  = [1;0;0];
        vecRow  = [0;-1;0];
        vecNorm = [0;0;1];
        fprintf('Image orientation information not available.\n')
    end
    % Rotate and flip the image if necessary
    [~, Mdir] = max(abs(vecNorm));
    [~, Midx] = max(abs([vecRow vecCol]), [], 2);
    
    if Mdir == 1        % sagittal plane
        if Midx(3) == 1     % align the column vector to S/I direction
            temp   = vecCol;
            vecCol = vecRow;
            vecRow = temp;
            image  = permute(image,[2 1 3 4]);
        end
        if vecRow(2) < 0    % align image left/right to A/P
            image = flip(image,2);
        end
        if vecCol(3) > 0    % align image top/bottom to S/I
            image = flip(image,1);
        end 
    elseif Mdir == 2    % coronal plane
        if Midx(3) == 1     % align the column vector to S/I direction
            temp   = vecCol;
            vecCol = vecRow;
            vecRow = temp;
            image  = permute(image,[2 1 3 4]);
        end
        if vecRow(1) < 0    % align image left/right to R/L
            image = flip(image,2);
        end
        if vecCol(3) > 0    % align image top/bottom to S/I
            image = flip(image,1);
        end 
    elseif Mdir == 3    % transverse plane
        if Midx(2) == 1     % align the column vector to A/P direction
            temp   = vecCol;
            vecCol = vecRow;
            vecRow = temp;
            image  = permute(image,[2 1 3 4]);
        end
        if vecRow(1) < 0    % align image left/right to R/L
            image = flip(image,2);
        end
        if vecCol(2) < 0    % align image top/bottom to A/P
            image = flip(image,1);
        end 
    end
end

