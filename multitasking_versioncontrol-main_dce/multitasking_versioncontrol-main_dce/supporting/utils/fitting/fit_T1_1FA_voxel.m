function fitResult = fit_T1_1FA_voxel(curve, fitParams)

% load fitting parameters, overwrite initial values
extractVarFromStruct(fitParams);

% inline functions
vec   = @(x) x(:);
row   = @(x) x(:).';
ppinv = @(x,y)(y*x')/norm(x)^2; %fast right-sided pseudoinverse function (for later)


e1  = @(T1) exp(-TR./T1);
Mss = @(e1,alpha)(1-e1)./(1-cos(alpha)*e1);

Sint = @(A,e1,B,Balpha) vec(A * Mss(e1,Balpha) * (1 - (B+1)*(e1*cos(Balpha)).^(vec(ns)-1)) * sin(Balpha));
S    = @(A,T1,B,Beta)   row(Sint(A,e1(T1),B,Beta*alphaArray(1)));


%% Fitting for curve

normcurve = curve(end);
curve = curve./normcurve;

Avp = @(T1,BIR,Beta) ppinv(S(1,T1,BIR,Beta).*fitw,curve.*fitw);      % parameterize solution to A as function of R1,B
[tempfit, res] = lsqnonlin(@(x)abs(S(Avp(x(1:wallClock_count),x(wallClock_count+1),x(wallClock_count+2)),x(1:wallClock_count),x(wallClock_count+1),x(wallClock_count+2))-curve).*fitw,...
                           xinit, xlb, xub, opts);               
fitResult = [Avp(tempfit(1:wallClock_count),tempfit(wallClock_count+1),tempfit(wallClock_count+2))*normcurve tempfit res];
        

