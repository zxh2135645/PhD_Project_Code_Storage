function [Mask_Segn_3D] = SpreadLabels(Mask_Segn)
% convert 2D labels to 3D binary labels
    num_labels = length(unique(Mask_Segn));
    Mask_Segn_3D = zeros(size(Mask_Segn, 1), size(Mask_Segn, 2), num_labels-1);
    for i = 1:(num_labels-1)
        Mask_Segn_2D = zeros(size(Mask_Segn));
        Mask_Segn_2D(Mask_Segn == i) = 1;
        Mask_Segn_3D(:,:,i) = Mask_Segn_2D;
    end
end