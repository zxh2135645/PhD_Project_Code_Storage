function recon = genImageTimeSeries(U, Nd, Phi, ds, Nim, useL)
    if nargin < 6
        useL = 1:size(Phi,1);
    end
    for l = numel(useL):-1:1
        if useL(l)>size(Phi,1)
            useL(l) = [];
        end
    end
    if nargin < 5
        Nim = min(1000,size(Phi,2));
    end
    
    if nargin < 4
        ds = 1;
    end
    
    Nim = min(Nim, floor(size(Phi,2)/ds));
    Phi_display = Phi(useL, (floor(ds/2)+1):ds:ds*Nim);
    
    recon = reshape(U,[],size(Phi,1));
    recon = recon(:,useL);
    recon = (reshape(recon*Phi_display,Nd(1), Nd(2), Nd(3),[]));
    