function contourData = draw_freehand_contours_on_cine(cineData, seriesIdx, frameIdx)
% draw_freehand_contours_on_cine
%
% Freehand-draw contours on cine images from cineData.
%
% INPUTS:
%   cineData   : struct returned by read_all_cine_dicoms
%   seriesIdx  : index of cineData entry to use
%   frameIdx   : cardiac phase to display for drawing (default = 1)
%
% OUTPUT:
%   contourData : struct with fields
%       .seriesIdx
%       .frameIdx
%       .baseName
%       .masks         [Ny Nx Nz] logical
%       .contours      {Nz x 1} cell, each cell is [x y] contour points
%       .sliceDrawn    [Nz x 1] logical
%
% USAGE:
%   contourData = draw_freehand_contours_on_cine(cineData, 1, 1);

    if nargin < 3 || isempty(frameIdx)
        frameIdx = 1;
    end

    % ----------------------------
    % Basic checks
    % ----------------------------
    assert(seriesIdx >= 1 && seriesIdx <= numel(cineData), ...
        'seriesIdx out of range.');

    V = cineData(seriesIdx).vol;   % [Ny Nx Nz Nt]
    baseName = cineData(seriesIdx).baseName;

    assert(ndims(V) == 4, 'cineData(seriesIdx).vol must be [Ny Nx Nz Nt].');

    [Ny, Nx, Nz, Nt] = size(V);

    assert(frameIdx >= 1 && frameIdx <= Nt, ...
        'frameIdx must be between 1 and %d.', Nt);

    masks = false(Ny, Nx, Nz);
    contours = cell(Nz, 1);
    sliceDrawn = false(Nz, 1);

    fprintf('Drawing contours for series: %s\n', baseName);
    fprintf('Using cardiac phase/frame: %d\n', frameIdx);
    fprintf('Number of slices: %d\n\n', Nz);

    % ----------------------------
    % Loop over slices
    % ----------------------------
    for z = 1:Nz
        I = V(:, :, z, frameIdx);

        hFig = figure('Name', sprintf('%s | Slice %d/%d | Frame %d', ...
            baseName, z, Nz, frameIdx), ...
            'NumberTitle', 'off', 'Color', 'k');

        imagesc(I);
        axis image off;
        colormap gray;
        title(sprintf('%s\nSlice %d/%d, Frame %d', baseName, z, Nz, frameIdx), ...
            'Interpreter', 'none', 'Color', 'w', 'FontWeight', 'bold');

        % contrast helper
        clim = prctile(I(:), [1 99]);
        if all(isfinite(clim)) && clim(2) > clim(1)
            caxis(clim);
        end

        drawnow;

        fprintf('Slice %d/%d:\n', z, Nz);
        fprintf('  Draw contour with mouse.\n');
        fprintf('  Double-click to finish.\n');
        fprintf('  Press Enter without drawing to skip this slice.\n');

        try
            h = drawfreehand('Color', 'r', 'LineWidth', 1.5);
        catch
            % older MATLAB fallback
            h = imfreehand(gca);
        end

        % If user closed the figure or did not draw
        if isempty(h) || ~isvalid_handle_like(h)
            close_if_valid(hFig);
            fprintf('  Skipped slice %d.\n\n', z);
            continue;
        end

        % Wait for completion
        try
            pos = wait(h);
        catch
            try
                pos = h.Position;
            catch
                pos = [];
            end
        end

        if isempty(pos)
            close_if_valid(hFig);
            fprintf('  Skipped slice %d.\n\n', z);
            continue;
        end

        % Create mask
        try
            BW = createMask(h);
        catch
            % fallback for old imfreehand object
            BW = h.createMask();
        end

        masks(:, :, z) = logical(BW);
        contours{z} = pos;
        sliceDrawn(z) = true;

        hold on;
        plot(pos(:,1), pos(:,2), 'y-', 'LineWidth', 1.5);
        drawnow;

        fprintf('  Saved contour for slice %d.\n\n', z);

        uiwait(msgbox(sprintf('Slice %d done. Click OK for next slice.', z), ...
            'Continue', 'modal'));

        close_if_valid(hFig);
    end

    % ----------------------------
    % Package output
    % ----------------------------
    contourData = struct();
    contourData.seriesIdx   = seriesIdx;
    contourData.frameIdx    = frameIdx;
    contourData.baseName    = baseName;
    contourData.masks       = masks;
    contourData.contours    = contours;
    contourData.sliceDrawn  = sliceDrawn;

    fprintf('Finished drawing.\n');
    fprintf('Drawn slices: %d / %d\n', nnz(sliceDrawn), Nz);
end

% =========================
% Helpers
% =========================
function tf = isvalid_handle_like(h)
    tf = false;
    if isempty(h)
        return;
    end
    try
        tf = isvalid(h);
    catch
        try
            tf = ishandle(h);
        catch
            tf = false;
        end
    end
end

function close_if_valid(h)
    try
        if ~isempty(h) && isvalid(h)
            close(h);
        end
    catch
        try
            if ishandle(h)
                close(h);
            end
        catch
        end
    end
end