function AhAx = AhA_ps3D_cufinufft(x, st, SEs, Phi2)

if ~isfield(st,'cufinufft_AhA_scalar')
    st.cufinufft_AhA_scalar = 1;
end

Ncoils = size(SEs,4);

try
    % Call mex function and apply scalar
    if numel(x)*Ncoils < 2e8
        AhAx = AhA_ps_cufinufft(single(x), st, single(SEs), single(Phi2)) * single(st.cufinufft_AhA_scalar);
    else
        AhAx = single(zeros(size(x)));
        % AhA_ps_cufinufft needs SEs to have at least 2 coils , so we do it
        % in 2-coil blocks
        nBlocks = floor(Ncoils/2);
        for nc = 1:nBlocks-1
            AhAx = AhAx + AhA_ps_cufinufft(single(x), st, single(SEs(:,:,:,(nc-1)*2+(1:2))), single(Phi2)) * single(st.cufinufft_AhA_scalar);
        end
        AhAx = AhAx + AhA_ps_cufinufft(single(x), st, single(SEs(:,:,:,(nBlocks-1)*2+(1:2+mod(Ncoils,2)))), single(Phi2)) * single(st.cufinufft_AhA_scalar);
    end
catch errormsg
    fprintf(2,'AhA_ps3D_cufinufft failed. \n');
    fprintf(2,'%s\n', errormsg.message);
end

AhAx = AhAx(:);
