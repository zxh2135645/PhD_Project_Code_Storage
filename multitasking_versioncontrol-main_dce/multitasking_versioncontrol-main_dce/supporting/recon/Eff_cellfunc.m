function Eff_result = Eff_cellfunc(T1,T2,T1rho,alphaArray,BIR,BT2,TR,TEs,TSLs,Nseg,invSign)

if nargin < 11
    invSign = -(sign(TEs) + sign(TSLs));
end

Mss_scale = @(e1, alpha, alpha_prev) (1-cos(alpha_prev)*e1) / (1-cos(alpha)*e1);

moduleLength = numel(TEs);
Eff_result   = zeros(1,moduleLength);

% approach 2, explicitly writing out the linear equation
Eff_result(1) = solve_Eff(T1,T2,T1rho,alphaArray,BIR,BT2,TR,TEs,TSLs,Nseg,invSign);
for shot = 2:moduleLength
    FA1 = alphaArray(mod(shot-2,moduleLength) + 1);
    FA2 = alphaArray(mod(shot-3,moduleLength) + 1);

    % signal model 1
    if  TEs(shot-1) == 0 && TSLs(shot-1) == 0   % IR pulse
        Eff_result(shot) = Mss_scale(exp(-TR/T1),FA1,FA2) * (1 + ((exp(-TR/T1)*cos(FA1))^Nseg).'*(cos(BIR*pi)*Eff_result(shot-1)-1)); %B is relative B1 field, not inversion efficiency
    else                                        % T2-IR pulse
        Eff_result(shot) = Mss_scale(exp(-TR/T1),FA1,FA2) * (1 + ((exp(-TR/T1)*cos(FA1))^Nseg).'*((cos(BT2*pi/2)^2 + invSign(shot-1)*sin(BT2*pi/2)^2*exp(-TEs(shot-1)/T2)*exp(-TSLs(shot-1)/T1rho))*Eff_result(shot-1)-1)); 
    end
    
%     % signal model 2
%     Eff_result(shot) = Mss_scale(exp(-TR/T1),FA1,FA2) * (1 + ((exp(-TR/T1)*cos(FA1))^Nseg).'*(cos(BIR*pi)*(cos(BT2*pi/2)^2+sin(BT2*pi/2)^2*exp(-TEs(shot-1)/T2)*exp(-TSLs(shot-1)/T1rho))*Eff_result(shot-1)-1)); 
    
end

end     % end of Eff_cellfunc_VFA()


function res = solve_Eff(T1,T2,T1rho,alphaArray,BIR,BT2,TR,TEs,TSLs,Nseg,invSign)

Mss_scale = @(e1, alpha, alpha_prev) (1-cos(alpha_prev)*e1) / (1-cos(alpha)*e1);
moduleLength = numel(TEs);
a = 1;
b = 0;
for shot = 1:moduleLength
    FA1 = alphaArray(mod(shot-1,moduleLength) + 1);
    FA2 = alphaArray(mod(shot-2,moduleLength) + 1);

    % signal model 1
    if  TEs(shot) == 0 && TSLs(shot) == 0   % IR pulse
        c = Mss_scale(exp(-TR/T1),FA1,FA2) * ((exp(-TR/T1)*cos(FA1))^Nseg).'*cos(BIR*pi);
        d = Mss_scale(exp(-TR/T1),FA1,FA2) * (1 + ((exp(-TR/T1)*cos(FA1))^Nseg).'*(-1));
    else                                    % T2-IR pulse
        if TEs(shot) > 0
            c = Mss_scale(exp(-TR/T1),FA1,FA2) * ((exp(-TR/T1)*cos(FA1))^Nseg).'*((cos(BT2*pi/2)^2 + invSign(mod(shot-2,moduleLength)+1)*sin(BT2*pi/2)^2*exp(-TEs(shot)/T2)*exp(-TSLs(shot)/T1rho))); 
        else
            c = Mss_scale(exp(-TR/T1),FA1,FA2) * ((exp(-TR/T1)*cos(FA1))^Nseg).'*((cos(BT2*pi/2)^2 + invSign(mod(shot-2,moduleLength)+1)*sin(BT2*pi/2)^2*exp(-TEs(shot)/T2)*exp(-TSLs(shot)/T1rho))); 
        end
        d = Mss_scale(exp(-TR/T1),FA1,FA2) * (1 + ((exp(-TR/T1)*cos(FA1))^Nseg).'*(-1)); 
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

