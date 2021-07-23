sizez = size(curves);
R1s = 1./logspace(log10(.1),log10(2),401);
t1_map_1d = zeros(1, size(ipt_2d, 1));


parfor ii = (1:size(ipt_2d, 1))
    if mask_1d(ii)
        ipt_1d = ipt_2d(ii, :);
        ipt_1d_norm = ipt_1d ./ max(ipt_1d);
        rms = zeros(sizez(2:end));
        for jj=1:sizez(2)
            for kk=1:sizez(3)
                for ll=1:sizez(4)
                    curves_norm = curves(:,jj,kk,ll)'./max(curves(:,jj,kk,ll));
                    temp = ipt_1d_norm - curves_norm;
                    rms(jj,kk,ll) = sum(temp.^2);
                end
            end
        end
        
        [minimum, I] = min(rms(:));
        sz = sizez(2:end);
        [I1, I2, I3] = ind2sub(sz,I);
        t1_map_1d(ii) = 1./R1s(I1);
    else
        
    end
end

t1_map_2d = reshape(t1_map_1d, Ny, Nx).*1000;