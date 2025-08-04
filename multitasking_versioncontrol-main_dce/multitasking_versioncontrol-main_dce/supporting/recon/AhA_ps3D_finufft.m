function [AhAx,lowmem] = AhA_ps3D_finufft(x, st, SEs, Phi2)

if isfield(st,'nufftlowmem')
    lowmem = st.nufftlowmem;
else
    lowmem = 0;
end

L      = size(Phi2,1);
Ntrajs = size(Phi2,3);
Ny = st.Nd(1);
Nx = st.Nd(2); 
Nz = st.Nz;
NM = st.M;
Nkx    = NM/Ntrajs;
Ncoils = size(SEs,4);

try 
    if lowmem == 0
        x = reshape(x, Ny, Nx, Nz, 1, L);
        FU = finufft_SoS(bsxfun(@times, x, SEs), st);
        FU = reshape(FU, Ntrajs, Nkx, Nz, Ncoils, L);
        Phi2temp = permute(Phi2,[3 5 4 6 1 2]);
        FU = reshape(sum(FU.*Phi2temp,5),NM,Nz,[]);
        AhAx = finufft_SoS_adj(FU, st);
        AhAx = reshape(AhAx, Ny, Nx, Nz, Ncoils, L);
        AhAx = sum(bsxfun(@times, (AhAx), conj(SEs)), 4);
    end
catch errormsg
    fprintf('%s\n', errormsg.message);
    fprintf('Try reducing memory usage by looping through trajectories... \n');
    lowmem = 1;
end

try
    if lowmem == 1
        x = reshape(x, Ny, Nx, Nz, 1, L);
        FU = finufft_SoS(bsxfun(@times, x, SEs), st);
        FU = reshape(FU, Ntrajs, Nkx, Nz, Ncoils, L);
        for traj = 1:Ntrajs
            for np = 1:Nz
                FU(traj,:,np,:,:) = reshape(reshape(FU(traj,:,np,:,:),[],L) * Phi2(:,:,traj,np),Nkx,Ncoils,L);
            end
        end
        FU   = reshape(FU, NM, Nz, []);
        AhAx = finufft_SoS_adj((FU), st);
        AhAx = reshape(AhAx, Ny, Nx, Nz, Ncoils, L);
        AhAx = sum(bsxfun(@times, (AhAx), conj(SEs)),4);
    end
catch errormsg
    fprintf('%s\n', errormsg.message);
    fprintf('Try reducing memory usage by looping through coils... \n');
    lowmem = 2;
end

try
    if lowmem == 2
        x = reshape(x, Ny, Nx, Nz, L);
        AhAx = zeros(size(x));
        for coil = 1:Ncoils
            xc = bsxfun(@times, x, SEs(:,:,:, coil));
            FU = finufft_SoS(xc, st);
            FU = reshape(FU, Ntrajs, Nkx, Nz, L);
            for traj = 1:Ntrajs
                for np = 1:Nz
                    FU(traj,:,np,:) = reshape(FU(traj,:,np,:),[],L) * Phi2(:,:,traj,np);
                end
            end
            FU  = reshape(FU, NM, Nz, []);
            FSU = finufft_SoS_adj(FU, st);
            FSU = reshape(FSU, Ny, Nx, Nz, L);
            AhAx = AhAx + bsxfun(@times, FSU, conj(SEs(:,:,:,coil)));
        end
    end
catch errormsg
    fprintf(2,'AhA_ps3D_finufft failed. \n');
    fprintf(2,'%s\n', errormsg.message);
end

AhAx = AhAx(:);


% function AhAx = AhA_ps3D_finufft(x, st, SEs, Phi2)
% 
% L      = size(Phi2,1);
% Ntrajs = size(Phi2,3);
% Nx = st.Nd(1);
% Ny = st.Nd(2); 
% Nz = st.Nz;
% NM = st.M;
% Nkx    = NM/Ntrajs;
% Ncoils = size(SEs,4);
% 
% if st.nufftlowmem == 0
%     x = reshape(x, Nx, Ny, Nz, 1, L);
%     x = bsxfun(@times, x, SEs);
%     FU = finufft_SoS(x, st);
%     FU = reshape(FU, Ntrajs, Nkx, Nz, Ncoils, L);
%     for traj = 1:Ntrajs
%         for np = 1:Nz
%             FU(traj,:,np,:,:) = reshape(reshape(FU(traj,:,np,:,:),[],L) * Phi2(:,:,traj,np),Nkx,Ncoils,L);
%         end
%     end
%     FU   = reshape(FU, NM, Nz, []);
%     AhAx = finufft_SoS_adj((FU), st);
%     AhAx = reshape(AhAx, Nx, Ny, Nz, Ncoils, L);
%     AhAx = sum(bsxfun(@times, (AhAx), conj(SEs)),4);
% else
%     x = reshape(x, Nx, Ny, Nz, L);
%     AhAx = zeros(size(x));
%     for coil = 1:Ncoils
%         xc = bsxfun(@times, x, SEs(:,:,:, coil));
%         FU = finufft_SoS(xc, st);
%         FU = reshape(FU, Ntrajs, Nkx, Nz, L);
%         for traj = 1:Ntrajs
%             for np = 1:Nz
%                 FU(traj,:,np,:) = reshape(FU(traj,:,np,:),[],L) * Phi2(:,:,traj,np);
%             end
%         end
%         FU  = reshape(FU, NM, Nz, []);
%         FSU = finufft_SoS_adj(FU, st);
%         FSU = reshape(FSU, Nx, Ny, Nz, L);
%         AhAx = AhAx + bsxfun(@times, FSU, conj(SEs(:,:,:,coil)));
%     end
% end
% 
% AhAx = AhAx(:);
% 
%     
% 
% 
