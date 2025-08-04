function AhAx = AhA_ps3D(x, st, SEs, Phi2)

try
    L      = size(Phi2,1);
    Ntrajs = size(Phi2,3);
    Nz     = size(Phi2,4);
    Ncoils = size(SEs,4);
    Nkx    = st.M/Ntrajs;
    
    prep = @(x,st) reshape(x, st.Nd(1),st.Nd(2), st.Nz, [], L);
    %F = gpuNUFFT(st.om.'/(2*pi), ones(1,st.M), st.osf, 6, 8, st.Nd(1:2), [], true);    
    
    AhAx = zeros(numel(x),1);
    for coil = 1:Ncoils
        FU = NUFFT_SoS(bsxfun(@times, prep(x,st), SEs(:,:,:,coil)), st);
        FU = reshape(FU, Ntrajs, Nkx, Nz, L);
        for traj = 1:Ntrajs
            for np = 1:Nz
                FU(traj,:,np,:) = reshape(FU(traj,:,np,:),[],L) * Phi2(:,:,traj,np);
            end
        end
        FU   = reshape(FU, st.M, Nz, 1, L);
        AhAx = AhAx + reshape(bsxfun(@times, NUFFT_SoS_adj(FU, st),conj(SEs(:,:,:,coil))),[],1);
    end
catch errormsg
    disp(errormsg);
    pause(rand*30);
    delete(gcp('nocreate'));  % delete current parpool
    parpool('local',3);       % open parpool with 8 workers
    gpuWorkerReset();         % reassign workers
    AhAx = AhA_ps3D(x,st,SEs,Phi2,flagDisplay);
end

