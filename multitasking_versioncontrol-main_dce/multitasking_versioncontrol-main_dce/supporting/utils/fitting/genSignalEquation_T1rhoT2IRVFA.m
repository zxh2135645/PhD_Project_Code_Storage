function fitParams = genSignalEquation_T1rhoT2IRVFA(fitParams,ns)

%% load fitting parameters

extractVarFromStruct(fitParams);

if nargin < 2
    ns = 1:Nseg;
end

%% Signal equation

e1  = @(T1) exp(-TR /T1);
e2  = @(T2) exp(-TEs/T2);
e1rho = @(T1rho)exp(-TSLs/T1rho);

Mss      = @(e1,alphas) (1-e1) ./ (1-cos(alphas)*e1);
step     = @(e1,alphas) (bsxfun(@power, e1*cos(alphas).', (ns-1))).';
sin_step = @(alphas)    sin(alphas);

Sint = @(A,e1,e2,e1rho,BalphaArray,BIR,BT2,Eff) A .* Mss(e1,BalphaArray) .* (1 + (step(e1,BalphaArray)) .* ((cos(BIR*pi).*IRs + (invSign.*sin(BT2*pi/2)^2).*(T2s+T1rhos).*e2.*e1rho + cos(BT2*pi/2)^2.*(T2s+T1rhos)).*Eff-1)) .* sin_step(BalphaArray);

if Nz > 1 && MBfactor == 1
    S = @(A,T1,T2,T1rho,Beta,BIR,BT2) row(Sint(  A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray, BIR, BT2, Eff_cellfunc(T1,T2,T1rho,Beta*alphaArray,BIR,BT2,TR,TEs,TSLs,Nseg,invSign)));
else
    S = @(A,T1,T2,T1rho,Beta,BIR,BT2) row(Sint(  A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.9387, BIR, BT2, Eff_cellfunc(T1,T2,T1rho,Beta*alphaArray*.9387,BIR,BT2,TR,TEs,TSLs,Nseg,invSign))...
                                        + Sint(2*A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.5049, BIR, BT2, Eff_cellfunc(T1,T2,T1rho,Beta*alphaArray*.5049,BIR,BT2,TR,TEs,TSLs,Nseg,invSign))...
                                        + Sint(2*A, e1(T1), e2(T2), e1rho(T1rho), Beta*alphaArray*.0525, BIR, BT2, Eff_cellfunc(T1,T2,T1rho,Beta*alphaArray*.0525,BIR,BT2,TR,TEs,TSLs,Nseg,invSign)));
end


%% store signal equation

fitParams.handles.e1 = e1;
fitParams.handles.e2 = e2;
fitParams.handles.Mss       = Mss;
fitParams.handles.Mss_scale = Mss_scale;
fitParams.handles.Eff  = Eff;
fitParams.handles.Sint = Sint;
fitParams.handles.S    = S;

