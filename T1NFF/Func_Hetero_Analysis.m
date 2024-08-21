function [Q] = Func_Hetero_Analysis(slc, roi_in_myo_t1, t1)
% Heterogeneity analysis adapted by Frank J Brooks, 2013 paper
% Quantification of heterogeneity observed in medical images
% slc = 4;

[X,Y] = find(roi_in_myo_t1(:,:,slc));
C_array = 1:length(X);
C = nchoosek(C_array, 2);

%Xm = 59; Ym = 52; Im = 1221;
%Xn = 62; Yn = 67; In = 1383;
Delta_I_ens = zeros(size(C,1),2);

% example of i = 7534
for i = 1:size(C,1)
    %for i = 7534:7534
    %for i = 7:7
    Xm = X(C(i,1)); Ym = Y(C(i,1));
    Im = t1(Xm, Ym, slc);
    Xn = X(C(i,2)); Yn = Y(C(i,2));
    In = t1(Xn, Yn, slc);
    %if Xm ~= Xn && Ym ~= Yn
    rmn = sqrt((Xm-Xn).^2 + (Ym-Yn).^2);
    Imn = In - Im;

    [XL, YL] = bresenham(Xm, Ym, Xn, Yn);
    %figure();
    %plot(XL, YL, 'or');

    % I_rmL = zeros(length(XL), 1);
    Delta_I_array = zeros(length(XL), 1);
    L = length(XL);
    count = 0;
    for l = 1:L
        Xl = XL(l); Yl = YL(l);
        Il = t1(Xl, Yl, slc) .* (roi_in_myo_t1(Xl,Yl,slc));

        if Il ~= 0 % To make it shape unaware
            rml = sqrt((Xm-Xl).^2 + (Ym-Yl).^2);
            I_rml = Im + (In - Im)/rmn * rml;
            Delta_I_array(l) = abs(I_rml - Il);
        else
            count = count + 1;
        end
    end

    L = L - count; % To make it shape unaware

    Delta_I = sum(Delta_I_array) / L;
    Delta_I_ens(i,1) = L;
    Delta_I_ens(i,2) = Delta_I;
    %end
end

L_array = unique(Delta_I_ens(:,1));
L_max = max(L_array);
L_array_norm = L_array / L_max;
Delta_I_ens_avg = zeros(length(L_array), 1);
for l = 1:length(L_array)
    L = L_array(l);
    [idx] = find(Delta_I_ens(:,1) == L);
    Delta_I_ens_avg(l) = mean(Delta_I_ens(idx,2));
end
%figure();
%plot(L_array_norm, Delta_I_ens_avg);
Q = trapz(L_array_norm,Delta_I_ens_avg);
% Q

end