function center_mask = Func_DrawCenterLine(img, caxis_rg)
figure('Position', [100 0 1600 1600]); imagesc(img); axis image;
caxis(caxis_rg);
center_line = drawpolygon(gca);
center_coords = center_line.Position;
center_line = drawpolygon(gca,'Position', [center_coords(:,1), center_coords(:,2)]);
center_mask = createMask(center_line);
end