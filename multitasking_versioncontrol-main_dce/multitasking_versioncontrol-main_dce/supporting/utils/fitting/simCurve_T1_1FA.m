function simResult = simCurve_T1_1FA(fitResult,position,cIdx,curve)

if nargin < 4
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

%% Setup functions

%ns   = 1:size(Phi,2);

e1  = @(T1) exp(-TR /T1);
Mss = @(e1,alpha)(1-e1)/(1-cos(alpha)*e1);

Sint = @(A,e1,BIR,Balpha) A * Mss(e1,Balpha) * (1 - (BIR+1)*(e1*cos(Balpha)).^(ns-1)) * sin(Balpha);
S    = @(A,T1,BIR,Beta)   row(Sint(A,e1(T1),BIR,Beta*alphaArray(1)));

ppinv = @(x,y)(y*x')/norm(x)^2; %fast right-sided pseudoinverse function (for later)


%% sim signal
     
if nargin > 3 
    if exist('fitw','var')
        normcurve = norm(curve(fitw>0));
    else
        normcurve = norm(curve);
    end
else
    curve = ones(1,numel(ns)*moduleLength);
end

pvalue = [T1map(y,x,z,cIdx),BIRmap(y,x,z,cIdx),B1map(y,x,z,cIdx)];

Avp  = @(T1,BIR,Beta) ppinv(S(1,T1,BIR,Beta).*fitw,curve.*fitw);
simResult = S(Avp(pvalue(1),pvalue(2),pvalue(3)),pvalue(1),pvalue(2),pvalue(3));            

normsimResult = norm(simResult(fitw>0));
simResult = simResult(:)*normcurve/normsimResult;

%% Fitting with 20 cardiac motions

% alpha0 = initBeta * flipAngleArray(1)*pi/180;
% 
% if nargin > 3 
%     normcurve = curve(end);
%     curve = curve/normcurve;
% else
%     normcurve = 1;
%     curve = ones(1,Nseg);
% end
% 
% if ChooseConst == 0
%     Avp = @(T1,B,alpha) ppinv(S(1,T1,B,alpha),curve);      % parameterize solution to A as function of R1,B
%     simResult = S(Avp(x(1),x(2),x(3)),x(1),x(2),x(3));
% else
%     Avp = @(T1,B) ppinv(S(1,T1,B,alpha0),curve);           % parameterize solution to A as function of R1,B
%     simResult = S(Avp(x(1),x(2)),x(1),x(2),alpha0);
% end
% 
% simResult = simResult(:)*normcurve;


