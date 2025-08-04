function Eff_result = Eff_cellfunc_T1_2FA(T1,FA1,FA2,B,TR,Nseg)

moduleLength = 2;
Eff_result = zeros(1,moduleLength);

Mss_scale = @(e1, alpha, alpha_prev) (1-cos(alpha_prev)*e1) / (1-cos(alpha)*e1);
Eff_result(1) = solve_Eff(T1,FA1,FA2,B,TR,Nseg);
Eff_result(2) = Mss_scale(exp(-TR/T1),FA2,FA1) * (1 + ((exp(-TR/T1)*cos(FA2))^Nseg).'*(-B*Eff_result(1)-1)); 

end     % end of Eff_cellfunc_T1_2FA()


function res = solve_Eff(T1,FA1,FA2,B,TR,Nseg)

moduleLength = 2;
Mss_scale = @(e1, alpha, alpha_prev) (1-cos(alpha_prev)*e1) / (1-cos(alpha)*e1);
a = 1;
b = 0;
for shot = 1:moduleLength
    if mod(shot,2) == 1
        c = Mss_scale(exp(-TR/T1),FA1,FA2) * ((exp(-TR/T1)*cos(FA1))^Nseg).'*(-B);
        d = Mss_scale(exp(-TR/T1),FA1,FA2) * (1 + ((exp(-TR/T1)*cos(FA1))^Nseg).'*(-1));
    else
        c = Mss_scale(exp(-TR/T1),FA2,FA1) * ((exp(-TR/T1)*cos(FA2))^Nseg).'*(-B);
        d = Mss_scale(exp(-TR/T1),FA2,FA1) * (1 + ((exp(-TR/T1)*cos(FA2))^Nseg).'*(-1));
    end

    a = c*a;
    b = c*b + d;
end
a = a - 1;
res = -b/a;

end     % end of solve_Eff()
