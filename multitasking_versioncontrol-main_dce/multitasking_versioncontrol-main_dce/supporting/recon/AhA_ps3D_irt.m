function AhAx = AhA_ps3D_irt(x, st, SEs, Phi2, loopcoil)

if nargin < 5
    loopcoil = 0;
end

L      = size(Phi2,1);
Ntrajs = size(Phi2,3);
Nx = st.Nd(1);
Ny = st.Nd(2); 
Nz = st.Nz;
NM = st.M;
Nkx    = NM/Ntrajs;
Ncoils = size(SEs,4);

if loopcoil == 0
    x = reshape(x, Nx, Ny, Nz, 1, L);
    x = bsxfun(@times, x, SEs);
    FU = irt_nufft_SoS((x), st);
    FU = reshape(FU, Ntrajs, Nkx, Nz, Ncoils, L);
    for traj = 1:Ntrajs
        for np = 1:Nz
            FU(traj,:,np,:,:) = reshape(reshape(FU(traj,:,np,:,:),[],L) * Phi2(:,:,traj,np),Nkx,Ncoils,L);
        end
    end
    FU   = reshape(FU, NM, Nz, []);
    AhAx = irt_nufft_SoS_adj((FU), st);
    AhAx = reshape(AhAx, Nx, Ny, Nz, Ncoils, L);
    AhAx = sum(bsxfun(@times, (AhAx), conj(SEs)),4);
else
    x = reshape(x, Nx, Ny, Nz, L);
    AhAx = zeros(size(x));
    for coil = 1:Ncoils
        xc = bsxfun(@times, x, SEs(:,:,:, coil));
        FU = irt_nufft_SoS(conj(xc), st);
        FU = reshape(conj(FU), Ntrajs, Nkx, Nz, L);
        for traj = 1:Ntrajs
            for np = 1:Nz
                FU(traj,:,np,:) = reshape(FU(traj,:,np,:),[],L) * Phi2(:,:,traj,np);
            end
        end
        FU  = reshape(FU, NM, Nz, []);
        FSU = irt_nufft_SoS_adj(FU, st);
        FSU = reshape(FSU, Nx, Ny, Nz, L);
        AhAx = AhAx + bsxfun(@times, FSU, conj(SEs(:,:,:,coil)));
    end
end

AhAx = AhAx(:);

    


