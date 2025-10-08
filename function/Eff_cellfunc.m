function Eff_result = Eff_cellfunc(T1,T2,T1rho,alphaArray,BIR,BT2,TR,TEs,TSLs,Nseg,Preps,TRfill)

if nargin < 12
    TRfill = 0;
end

if nargin < 11
    Preps = -(TEs==0) + (TEs>0)*2 + (TSLs>0)*3;
end
invSign = Preps./abs(Preps);
invSign(Preps==0.5) = 0.5;

% if nargin < 11
%     invSign = -(sign(TEs) + sign(TSLs));
% end

Mss_scale = @(e1, alpha, alpha_prev) (1-cos(alpha_prev)*e1) / (1-cos(alpha)*e1);

moduleLength = numel(TEs);
Eff_result   = zeros(1,moduleLength);

% approach 2, explicitly writing out the linear equation
Eff_result(1) = solve_Eff(T1,T2,T1rho,alphaArray,BIR,BT2,TR,TEs,TSLs,Nseg,Preps,TRfill);
for shot = 2:moduleLength
    FA1 = alphaArray(mod(shot-2,moduleLength) + 1);
    FA2 = alphaArray(mod(shot-3,moduleLength) + 1);

    % signal model 1
    if abs(Preps(shot-1)) <= 1                  % IR/SR prep
        Eff_result(shot) = Mss_scale(exp(-TR/T1),FA1,FA2) * (1 + ((exp(-TR/T1)*cos(FA1))^Nseg).'*(cos(invSign(shot-1)*BIR*pi)*Eff_result(shot-1)-1)) * exp(-TRfill/T1) + (1 - exp(-TRfill/T1)); %B is relative B1 field, not inversion efficiency
    else                                        % T2/T1rho prep
        Eff_result(shot) = Mss_scale(exp(-TR/T1),FA1,FA2) * (1 + ((exp(-TR/T1)*cos(FA1))^Nseg).'*((cos(BT2*pi/2)^2 + invSign(shot-1)*sin(BT2*pi/2)^2*exp(-TEs(shot-1)/T2)*exp(-TSLs(shot-1)/T1rho))*Eff_result(shot-1)-1)) * exp(-TRfill/T1) + (1 - exp(-TRfill/T1)); 
    end
    
%     % signal model 2
%     Eff_result(shot) = Mss_scale(exp(-TR/T1),FA1,FA2) * (1 + ((exp(-TR/T1)*cos(FA1))^Nseg).'*(cos(BIR*pi)*(cos(BT2*pi/2)^2+sin(BT2*pi/2)^2*exp(-TEs(shot-1)/T2)*exp(-TSLs(shot-1)/T1rho))*Eff_result(shot-1)-1)); 
    
end

end     % end of Eff_cellfunc_VFA()


function res = solve_Eff(T1,T2,T1rho,alphaArray,BIR,BT2,TR,TEs,TSLs,Nseg,Preps,TRfill)
invSign = Preps./abs(Preps);
invSign(Preps==0.5) = 0.5;

Mss_scale = @(e1, alpha, alpha_prev) (1-cos(alpha_prev)*e1) / (1-cos(alpha)*e1);
moduleLength = numel(TEs);
a = 1;
b = 0;
for shot = 1:moduleLength
    FA1 = alphaArray(mod(shot-1,moduleLength) + 1);
    FA2 = alphaArray(mod(shot-2,moduleLength) + 1);

    % signal model 1
    if abs(Preps(shot)) <= 1                % IR/SR prep
%     if  TEs(shot) == 0 && TSLs(shot) == 0   % IR/SR prep
        c = Mss_scale(exp(-TR/T1),FA1,FA2) * ((exp(-TR/T1)*cos(FA1))^Nseg).'*cos(invSign(shot).*BIR*pi) * exp(-TRfill/T1);
        d = Mss_scale(exp(-TR/T1),FA1,FA2) * (1 + ((exp(-TR/T1)*cos(FA1))^Nseg).'*(-1)) * exp(-TRfill/T1) + (1 - exp(-TRfill/T1));
    else                                    % T2T1rho prep
%         if TEs(shot) > 0
            c = Mss_scale(exp(-TR/T1),FA1,FA2) * ((exp(-TR/T1)*cos(FA1))^Nseg).'*((cos(BT2*pi/2)^2 + invSign(shot)*sin(BT2*pi/2)^2*exp(-TEs(shot)/T2)*exp(-TSLs(shot)/T1rho))) * exp(-TRfill/T1); 
%         else
%             c = Mss_scale(exp(-TR/T1),FA1,FA2) * ((exp(-TR/T1)*cos(FA1))^Nseg).'*((cos(BT2*pi/2)^2 + invSign(mod(shot-1,moduleLength)+1)*sin(BT2*pi/2)^2*exp(-TEs(shot)/T2)*exp(-TSLs(shot)/T1rho))); 
%         end
        d = Mss_scale(exp(-TR/T1),FA1,FA2) * (1 + ((exp(-TR/T1)*cos(FA1))^Nseg).'*(-1)) * exp(-TRfill/T1) + (1 - exp(-TRfill/T1));
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

