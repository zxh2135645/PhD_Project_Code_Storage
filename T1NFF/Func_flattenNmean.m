function [Mipix_flat, Mipix_flat2, Mipix_flat3, Mipix_mean, Mipix_mean2, Mipix_mean3] = ...
    Func_flattenNmean(i, Mi_idx_intersect, Mipix, Mipix2, Mipix3, Mipix_mean, Mipix_mean2, Mipix_mean3)
    Mipix_flat = {};
    Mipix_flat2 = {};
    Mipix_flat3 = {};
    
    if ~isempty(Mi_idx_intersect)
        idx = 1;
        for j = 1:length(Mi_idx_intersect)
            ind = Mi_idx_intersect(j);
            Mipix_flat{idx} = Mipix{ind,i};
            Mipix_flat2{idx} = Mipix2{ind,i};
            Mipix_flat2{idx}(Mipix_flat2{idx} > 100) = 100;
            Mipix_flat2{idx}(Mipix_flat2{idx} < 0) = 0;
            
            Mipix_flat3{idx} = Mipix3{ind,i};
            Mipix_flat3{idx}(Mipix_flat3{idx} < 0) = 0;
            
            Mipix_mean(idx, i) = mean(Mipix_flat{idx});
            Mipix_mean2(idx, i) = mean(Mipix_flat2{idx});
            Mipix_mean3(idx, i) = mean(Mipix_flat3{idx});
            idx = idx + 1;
        end
        
    else
        Mipix_mean(:,i) = zeros(size(Mipix_mean, 1), 1);
        Mipix_mean2(:,i) = zeros(size(Mipix_mean2, 1), 1);
        Mipix_mean3(:,i) = zeros(size(Mipix_mean3, 1), 1);
    end

end