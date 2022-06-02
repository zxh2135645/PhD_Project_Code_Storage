function [exvivo_match, cine_match] = matching_slices(slice_data_cine, sep_ratio_cine, sep_ratio_exvivo_med)

% be used in Exvivo_Invivo_Matching
% Match slices CINE vs exvivo and save as gif
% Assuming exvivo always has larger coverage

thickness_cine = slice_data_cine(1).SliceThickness;

sep_ratio_cine(isnan(sep_ratio_cine)) = [];
n_cine = length(sep_ratio_cine);

sep_ratio_cine_mat = repmat(sep_ratio_cine, [n_cine-3+1, 1]);
temp = zeros(size(sep_ratio_cine_mat, 1), n_cine - 3);
for i = 1:size(sep_ratio_cine_mat, 1)
    temp(i,:) = cat(2, ones(1, n_cine - 3 - i + 1), zeros(1, i-1));
end
multiplier_mat = cat(2, ones(size(sep_ratio_cine_mat, 1), 3), temp);
multiplier_mat(multiplier_mat == 0) = nan;
sep_ratio_cine_mat = sep_ratio_cine_mat .* multiplier_mat;

a_array = zeros(size(sep_ratio_cine_mat, 1), 1);
b_array = zeros(size(sep_ratio_cine_mat, 1), 1);
sep_ratio_cine_interp_norm_cell = cell(size(sep_ratio_cine_mat, 1), 1);
exvivo_isnan_cell = cell(size(sep_ratio_cine_mat, 1), 1);

for k = 1:size(sep_ratio_cine_mat, 1)
    sep_ratio_cine_new = sep_ratio_cine_mat(k,:);
    sep_ratio_cine_new(isnan(sep_ratio_cine_new)) = [];
    n_cine_k = length(sep_ratio_cine_new);
    
    interp_slice_cine = ((n_cine_k-1)*thickness_cine+1);
    
    sep_ratio_cine_interp = interp1(1:thickness_cine:interp_slice_cine, sep_ratio_cine_new, 1:interp_slice_cine);
    %figure(); plot(sep_ratio_cine_interp, 'LineWidth', 1.5); hold on; plot(sep_ratio_exvivo_med, 'LineWidth', 1.5);
    sep_ratio_cine_interp_norm = sep_ratio_cine_interp ./ max(sep_ratio_cine_interp);
    
    sep_ratio_exvivo_med_norm = sep_ratio_exvivo_med ./ max(sep_ratio_exvivo_med);
    %figure(); plot(sep_ratio_cine_interp_norm, 'LineWidth', 1.5); hold on; plot(sep_ratio_exvivo_med_norm, 'LineWidth', 1.5);
    
    exvivo_isnan_array = isnan(sep_ratio_exvivo_med_norm);
    sep_ratio_exvivo_med_norm(exvivo_isnan_array) = [];
        
    n = length(sep_ratio_exvivo_med_norm)-length(sep_ratio_cine_interp_norm)+1;
    
    RMS = zeros(1, n);
    
    for i = 1:n
        RMS(i) = sqrt(sum((sep_ratio_exvivo_med_norm(i:(length(sep_ratio_cine_interp_norm)+i-1)) - sep_ratio_cine_interp_norm).^2)) / length(sep_ratio_cine_interp_norm);
    end
    
    [a_array(k),b_array(k)] = min(RMS);
    
    sep_ratio_cine_interp_norm_cell{k} = sep_ratio_cine_interp_norm;
    exvivo_isnan_cell{k} = exvivo_isnan_array;
end

    [aa, bb] = min(a_array);
    
    x_shift = b_array(bb):(b_array(bb)+length(sep_ratio_cine_interp_norm_cell{bb})-1);
    figure(); plot(x_shift, sep_ratio_cine_interp_norm_cell{bb}, 'LineWidth', 1.5); hold on; plot(sep_ratio_exvivo_med_norm, 'LineWidth', 1.5);
    
    % Map back to 64 slices (nan was ignored and now should add back)
    first_nan = find(diff(exvivo_isnan_cell{bb}));
    first_nan = first_nan(1);
    
    sep_ratio_cine_new = sep_ratio_cine_mat(bb,:);
    sep_ratio_cine_new(isnan(sep_ratio_cine_new)) = [];
    n_cine_new = length(sep_ratio_cine_new);
    
    cine_match = 1:n_cine_new;
    exvivo_match = linspace(first_nan+b_array(bb), first_nan+b_array(bb)+(n_cine_new-1)*thickness_cine, n_cine_new);
end