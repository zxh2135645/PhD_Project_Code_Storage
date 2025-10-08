function Eff_result = Eff_cellfunc_VTR(T1,T2,T1rho,alphaArray,BIR,BT2,TR,TEs,TSLs,Nseg,Preps)

if nargin < 11
    Preps = -(TEs==0) + (TEs>0)*2 + (TSLs>0)*3;
end
invSign = Preps./abs(Preps);
invSign(Preps==0.5) = 0.5;

moduleLength = numel(TEs);
Eff_result   = zeros(1,moduleLength);

% approach 2, explicitly writing out the linear equation
Eff_result(1) = solve_Eff(T1,T2,T1rho,alphaArray,BIR,BT2,TR,TEs,TSLs,Nseg,Preps);

if numel(TR) > 1
    TRnav = TR(2);
    TR    = TR(1); 
else
    TRnav = TR;
end

e1    = @(T1) exp(-TR /T1);
e1nav = @(T1) exp(-TRnav/T1);
e2    = @(T2,shot)    exp(-TEs(shot-1)/T2);
e1rho = @(T1rho,shot) exp(-TSLs(shot-1)/T1rho);

Mss1 = @(e1, e1nav, alpha, alpha_prev) ((1-e1nav)+(1-e1)*e1nav*cos(alpha))*(1-e1*e1nav*cos(alpha_prev).^2) / ((1-e1*e1nav*cos(alpha).^2)*((1-e1nav)+(1-e1)*e1nav*cos(alpha_prev)));
Mss2 = @(e1, e1nav, alpha, alpha_prev) ((1-e1)+(1-e1nav)*e1*cos(alpha))*(1-e1*e1nav*cos(alpha_prev).^2) / ((1-e1*e1nav*cos(alpha).^2)*((1-e1)+(1-e1nav)*e1*cos(alpha_prev)));
Mss_scale = @(e1, e1nav, alpha, alpha_prev) Mss1(e1, e1nav, alpha, alpha_prev);

step_end  = @(T1,alpha) (bsxfun(@power, e1(T1)*e1nav(T1)*cos(alpha).^2, (ceil((Nseg+1)/2))));
    
for shot = 2:moduleLength
    FA1 = alphaArray(mod(shot-2,moduleLength) + 1);
    FA2 = alphaArray(mod(shot-3,moduleLength) + 1);

    % signal model 1
    if abs(Preps(shot-1)) <= 1                  % IR/SR prep
        Eff_result(shot) = Mss_scale(exp(-TR/T1),exp(-TRnav/T1),FA1,FA2) * (1 + (step_end(T1,FA1))*(cos(invSign(shot-1)*BIR*pi)*Eff_result(shot-1)-1)); %B is relative B1 field, not inversion efficiency
    else                                        % T2/T1rho prep
        Eff_result(shot) = Mss_scale(exp(-TR/T1),exp(-TRnav/T1),FA1,FA2) * (1 + (step_end(T1,FA1))*((cos(BT2*pi/2)^2 + invSign(shot-1)*sin(BT2*pi/2)^2*e2(T2,shot)*e1rho(T1rho,shot))*Eff_result(shot-1)-1)); 
    end
  
end

end     % end of Eff_cellfunc_VFA()


function res = solve_Eff(T1,T2,T1rho,alphaArray,BIR,BT2,TR,TEs,TSLs,Nseg,Preps)
invSign = Preps./abs(Preps);
invSign(Preps==0.5) = 0.5;

if numel(TR) > 1
    TRnav = TR(2);
    TR    = TR(1); 
else
    TRnav = TR;
end

e1    = @(T1) exp(-TR /T1);
e1nav = @(T1) exp(-TRnav/T1);

Mss1 = @(e1, e1nav, alpha, alpha_prev) ((1-e1nav)+(1-e1)*e1nav*cos(alpha))*(1-e1*e1nav*cos(alpha_prev).^2) / ((1-e1*e1nav*cos(alpha).^2)*((1-e1nav)+(1-e1)*e1nav*cos(alpha_prev)));
Mss2 = @(e1, e1nav, alpha, alpha_prev) ((1-e1)+(1-e1nav)*e1*cos(alpha))*(1-e1*e1nav*cos(alpha_prev).^2) / ((1-e1*e1nav*cos(alpha).^2)*((1-e1)+(1-e1nav)*e1*cos(alpha_prev)));
Mss_scale = @(e1, e1nav, alpha, alpha_prev) Mss1(e1, e1nav, alpha, alpha_prev);

step_end  = @(T1,alpha) (bsxfun(@power, e1(T1)*e1nav(T1)*cos(alpha).^2, (ceil((Nseg+1)/2))));

moduleLength = numel(TEs);
a = 1;
b = 0;
for shot = 1:moduleLength
    FA1 = alphaArray(mod(shot-1,moduleLength) + 1);
    FA2 = alphaArray(mod(shot-2,moduleLength) + 1);

    % signal model 1
    if abs(Preps(shot)) <= 1                % IR/SR prep
        c = Mss_scale(exp(-TR/T1),exp(-TRnav/T1),FA1,FA2) * (step_end(T1,FA1))*cos(invSign(shot).*BIR*pi);
        d = Mss_scale(exp(-TR/T1),exp(-TRnav/T1),FA1,FA2) * (1 + (step_end(T1,FA1))*(-1));
    else                                    % T2T1rho prep
        c = Mss_scale(exp(-TR/T1),exp(-TRnav/T1),FA1,FA2) * step_end(T1,FA1)*((cos(BT2*pi/2)^2 + invSign(shot)*sin(BT2*pi/2)^2*exp(-TEs(shot)/T2)*exp(-TSLs(shot)/T1rho))); 
        d = Mss_scale(exp(-TR/T1),exp(-TRnav/T1),FA1,FA2) * (1 + step_end(T1,FA1)*(-1)); 
    end

%     % signal model 2
%     c = Mss_scale(exp(-TR/T1),FA1,FA2) * ((exp(-TR/T1)*cos(FA1))^Nseg).'*(cos(BIR*pi)*(cos(BT2*pi/2)^2-sin(BT2*pi/2)^2*exp(-TEs(shot)/T2)*exp(-TSLs(shot)/T1rho))); 
%     d = Mss_scale(exp(-TR/T1),FA1,FA2) * (1 + ((exp(-TR/T1)*cos(FA1))^Nseg).'*(-1)); 
       
    a = c*a;
    b = c*b + d;
end
a = a - 1;
res = -b/a;

end     % end of solve_Eff()

