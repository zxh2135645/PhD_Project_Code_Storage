sizez = size(curves);
R1s = 1./logspace(log10(.1),log10(2),401);
t1_map_1d = zeros(1, size(ipt_2d, 1));
sz = sizez(2:end);

curves_norm = zeros(size(curves));
for jj=1:sizez(2)
    for kk=1:sizez(3)
        for ll=1:sizez(4)
            curves_norm(:,jj,kk,ll) = curves(:,jj,kk,ll)'./max(curves(:,jj,kk,ll));
        end
    end
end


parfor ii = (1:size(ipt_2d, 1))
    if mask_1d(ii)
        % tic;
        ipt_1d = ipt_2d(ii, :);
        ipt_1d_no1rm = ipt_1d ./ max(ipt_1d);
        % rms = zeros(sizez(2:end));

        
        ipt_4d_norm = repmat(ipt_1d_norm', [1, sz]);
        temp = (curves_norm - ipt_4d_norm).^2;
        rms = squeeze(sum(temp, 1));
        
        [minimum, I] = min(rms(:));
        [I1, I2, I3] = ind2sub(sz,I);
        t1_map_1d(ii) = 1./R1s(I1);
        % toc;
    else
        
    end
end

t1_map_2d = reshape(t1_map_1d, Ny, Nx).*1000;