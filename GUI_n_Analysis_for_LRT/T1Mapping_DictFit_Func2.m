
t1_map_1d = zeros(1, size(ipt_2d, 1));

for ii = (1:size(ipt_2d, 1))
    if mask_1d(ii)
       ipt_1d = ipt_2d(ii,:);
       ipt_1d_norm = ipt_1d ./ max(ipt_1d);
       
%        ipt_2d_norm = repmat(ipt_1d_norm', [1 size(Mz_dict_norm_abs, 1)]);
%        temp2 = (Mz_dict_norm_abs' - ipt_2d_norm).^2;
       ipt_2d_norm = repmat(ipt_1d_norm', [1 size(Mz_dict_norm_abs_truc, 1)]);
       temp2 = (Mz_dict_norm_abs_truc' - ipt_2d_norm).^2;
       
       rms = squeeze(sum(temp2, 1));
       
       [minimum, I] = min(rms(:));
       
       t1_map_1d(ii) = 1 ./ R1s(I);
    end
end

t1_map_2d = reshape(t1_map_1d, Ny, Nx) * 1000;