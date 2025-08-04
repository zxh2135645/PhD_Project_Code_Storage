function simResult = simCurve_T1rhoT2IRVFA(fitResult,position,cIdx,curve)

if nargin < 3
    cIdx = 1;
end

row = @(x) x(:).';

y = position(1); x = position(2); z = position(3);

if z > size(fitResult.T1map,3)
    z = 1;
end

extractVarFromStruct(fitResult);
extractVarFromStruct(fitParams);

if ~exist('fitw','var')
    fitw = curve*0 + 1;
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

ppinv = @(x,y)(y*x')/norm(x)^2; 

% sim signal
pvalue  = [T1map(y,x,z,cIdx),T2map(y,x,z,cIdx),T1rhomap(y,x,z,cIdx),B1map(y,x,z,cIdx),BIRmap(y,x,z,cIdx),BT2map(y,x,z,cIdx)];  
        
if nargin > 3 
    if exist('fitw','var')
        normcurve = norm(curve(fitw>0));
    else
        normcurve = norm(curve);
    end
else
    curve = ones(1,numel(ns)*moduleLength);
end

Avp  = @(T1,T2,T1rho,Beta,BIR,BT2) ppinv(S(1,T1,T2,T1rho,Beta,BIR,BT2).*fitw,curve.*fitw); %parameterize solution to A as function of T1, T2, B1_alpha, B1_IR, B1_T2IR, assume BT2=1 for now
simResult = S(Avp(pvalue(1),pvalue(2),pvalue(3),pvalue(4),pvalue(5),pvalue(6)),pvalue(1),pvalue(2),pvalue(3),pvalue(4),pvalue(5),pvalue(6));            

normsimResult = norm(simResult(fitw>0));
simResult = simResult(:)*normcurve/normsimResult;

