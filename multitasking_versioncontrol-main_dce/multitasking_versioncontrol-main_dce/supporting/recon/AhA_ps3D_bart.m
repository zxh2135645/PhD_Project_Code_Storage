function [AhAx, lowmem] = AhA_ps3D_bart(x, st, SEs, Phi2, lowmem)

if nargin < 4
    lowmem = 0;
end

L                 = size(Phi2,1);
[Ny,Nx,Nz,Ncoils] = size(SEs);
[~,Ntrajs,Nkx]    = size(st.traj_bart);

AhAx = complex(zeros([st.Nd Nz L]),0);

try
    if lowmem == 0
        x = reshape(x, st.Nd(1), st.Nd(2), st.Nz, 1, L);
        FU = NUFFT_SoS_bart(bsxfun(@times, x, SEs), st);
        FU = reshape(FU, Ntrajs, Nkx, Nz, Ncoils, L);
        tempPhi2 = permute(Phi2,[3 5 4 6 1 2]);
        FU = reshape(sum(FU.*tempPhi2,5),Ntrajs,Nkx,Nz,Ncoils,L);
        AhAx = sum(bsxfun(@times,reshape(NUFFT_SoS_adj_bart(FU,st),[size(SEs) L]),conj(SEs)),4);
    end
catch errormsg
    fprintf('%s\n ', errormsg.message);
    fprintf('\n Try looping through trajectories... \n');
    lowmem = 1;
end

try
    if lowmem == 1
        x = reshape(x, st.Nd(1), st.Nd(2), st.Nz, 1, L);
        FU = NUFFT_SoS_bart(bsxfun(@times, x, SEs), st);
        FU = reshape(FU, Ntrajs, Nkx, Nz, Ncoils, L);
        for traj = 1:Ntrajs
            for np = 1:Nz
                FU(traj,:,np,:,:) = reshape(reshape(FU(traj,:,np,:,:),[],L) * Phi2(:,:,traj,np),1,Nkx,1,Ncoils,L);
            end
        end
        AhAx = sum(bsxfun(@times,reshape(NUFFT_SoS_adj_bart(FU,st),[size(SEs) L]),conj(SEs)),4);
    end
catch errormsg
    fprintf('%s\n ', errormsg.message);
    fprintf('\n Try looping through coils and trajectories... \n');
    lowmem = 2;
end

try
    if lowmem == 2
        x = reshape(x, st.Nd(1), st.Nd(2), st.Nz, L);
        for coil = 1:Ncoils
            %fprintf('AhA_ps3D_bart coil %d.\n',coil);
            FU = NUFFT_SoS_bart(bsxfun(@times, x, SEs(:,:,:,coil)), st);
            FU = reshape(FU, Ntrajs, Nkx, Nz, L);
            for traj = 1:Ntrajs
                for np = 1:Nz
                    FU(traj,:,np,:) = reshape(FU(traj,:,np,:),[],L) * Phi2(:,:,traj,np);
                end
            end
            AhAx = AhAx + bsxfun(@times,NUFFT_SoS_adj_bart(FU,st),conj(SEs(:,:,:,coil)));
        end
    end
catch errormsg
    fprintf(2,'\n AhA_ps3D_bart failed. \n');
    fprintf(2,'%s\n ', errormsg.message);
end

AhAx = AhAx(:);
    


