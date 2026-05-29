function [idx_reordered] = Func_AlignSliceLoc(target_array, obj_array)
    idx_reordered = zeros(1, length(target_array));
    for i = 1:length(target_array)
        slc = target_array(i);
        [M, I] = min(round(abs(obj_array - slc), 4));
        idx_reordered(i) = I;
    end
end