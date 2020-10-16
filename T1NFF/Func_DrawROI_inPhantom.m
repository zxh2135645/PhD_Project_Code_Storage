function [roi, row, col, N] = Func_DrawROI_inPhantom(img, mask_epi, mask_save, caxis_rg, dim)

if ~exist(mask_save, 'file')
    roi = struct;
    if dim == 1
        N = input('Number of vials: ');
        vial_coords_cell = cell(N, 1);
        vial_mask_cell = cell(N, 1);
        for i = 1:N
            figure(); imagesc(img .* mask_epi); axis image; caxis(caxis_rg); colorbar;
            temp = drawpolygon(gca); % You may need to use roipoly and its equivalent if MATLAB version is earlier than 2020a
            vial_coords_cell{i} = temp.Position;
            vial_mask_cell{i} = createMask(temp);
        end
    elseif dim == 2
        row = input('Number of rows: ');
        col = input('Number of cols: ');
        N = row * col;
        vial_coords_cell = cell(row, col);
        vial_mask_cell = cell(row, col);
        for k = 1:row
            for j = 1:col
                figure(); imagesc(img .* mask_epi); axis image; caxis(caxis_rg); colorbar;
                temp = drawpolygon(gca); % You may need to use roipoly and its equivalent if MATLAB version is earlier than 2020a
                vial_coords_cell{k, j} = temp.Position;
                vial_mask_cell{k, j} = createMask(temp);
            end
        end
    end
    close all;
    roi.vial_coords_cell = vial_coords_cell;
    roi.vial_mask_cell = vial_mask_cell;
else
    roo = load(mask_save);
    fn = fieldnames(roo);    
    roi = roo.(fn{1});
    roi_more = roi.vial_mask_cell;
    row = size(roi_more, 1);
    col = size(roi_more, 2);
    N = row * col;
end
end