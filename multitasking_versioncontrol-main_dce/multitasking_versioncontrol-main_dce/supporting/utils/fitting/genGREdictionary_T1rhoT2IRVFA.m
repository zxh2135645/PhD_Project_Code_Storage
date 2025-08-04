function [curves_compress,basis_compress,fitParams] = genGREdictionary_T1rhoT2IRVFA(fitParams)

%% initialize fitting parameters

% load fitting parameters, overwrite initial values
extractVarFromStruct(fitParams);

% % slice profile
% halfSlice = numel(B1_coeff)/Nz/2;
% B1_coeff  = interp1((1:numel(B1_coeff))-0.5,B1_coeff,linspace(halfSlice,numel(B1_coeff)-halfSlice,Nz),'pchip');


%% Signal equation

e1  = @(T1) exp(-TR /T1);
e2  = @(T2) exp(-TEs/T2);
e1rho = @(T1rho)exp(-TSLs/T1rho);

Mss      = @(e1,alphas) (1-e1) ./ (1-cos(alphas).*e1);
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


%% Generate dictionaries

tic;

T1s = [300:20:1200,1240:40:1800,1900:100:2200,2400:200:3000]*1e-3;
T2s = [20:2:60,65:5:100,110:10:150]*1e-3;
BIRs = [0.5:0.05:0.75 0.76:0.02:1];
Betas = 0.05:0.05:1.0;

curves = zeros(Nseg*moduleLength,numel(T1s),numel(T2s),numel(Betas),numel(BIRs));
for j=1:numel(T1s)
    for k=1:numel(T2s)
        for l = 1:numel(Betas)
            for m = 1:numel(BIRs)
                curves(:,j,k,l,m) = S(1,T1s(j),T2s(k),BIRs(m),Betas(l),BIRs(m));
            end
        end
    end
end
curves = curves(:,:).';
[~,curveS,curveU] = svde(curves);
curveS = diag(curveS);
basis_compress  = curveU(:,1:100);
curves_compress = curves*basis_compress;

toc;

fitParams.T1s = T1s;
fitParams.T2s = T2s;
fitParams.Betas = Betas;
fitParams.BIRs = BIRs;
fitParams.basis_compress = basis_compress;
fitParams.curves_compress = curves_compress;

return