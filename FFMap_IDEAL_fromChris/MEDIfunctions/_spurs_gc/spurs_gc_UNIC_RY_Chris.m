%Modidied Randy Yang 2019 July modify for shimming field derivation
% Simultanes Phase Unwrapping and Removal of chemical Shift (SPURS) via Graph Cut

% [wwater wfat wfreq wunwph_uf unwphw N_std] = spurs_gc(iField,TE,CF,voxel_size,SUBSAMPLE)
%% HOW TO USE
% SPURS with full-matrix: [wwater wfat wfreq wunwph_uf unwphw N_std] = spurs_gc(iField,TE,CF,voxel_size);
% SPURS with subsample: [wwater wfat wfreq wunwph_uf unwphw N_std] = spurs_gc(iField,TE,CF,voxel_size,2);
% Note that if the water fat map totally swap, try conj(iField) instead of iField as input
%% output: 
%     - wwater: the water map 
%     - wfat:   the fat map 
%     - wfreq:  the field map in rad after running IDEAL as fine tunning,input for QSM
%     - wunwph_uf:  the field map after unwrapping and unfat,initial guess for IDEAL
%     - unwphw: phase unwrapping result

%% input: 
%      - iField : a multi-echo 4 dimentional data (Note that if the water
%      fat map totally swap, try conj(iField) instead of iField as input)
%      - how to choose iField or conj(iField) as input: 
%             if PrecessionIsClockwise = 1, [] = spurs_gc(conj(iField),TE,CF,voxel_size);
%             if PrecessionIsClockwise = -1, [] = purcs_gc(iField,TE,CF,voxel_size);

% written by Jianwu Dong  2014.2.10
% last modified by Jianwu 2014.9.3
% last modified by Jianwu, add voxel_size when computing the gradient

function [wunwph_uf ,unwphw ,N_std ,iFreq_raw, iFreq_ref,Mask0] = spurs_gc_UNIC_RY_Chris(L_Fat,iField,TE,CF,voxel_size,SUBSAMPLE,dfat)
if nargin<6
    SUBSAMPLE = 1;
end

energy = [];
iField0 = iField;
[sx sy sz necho] = size(iField);


iMag = sqrt(sum(abs(iField).^2,4));
[x,y,z,t] = size(iField);
Mask=autoMask(iMag,voxel_size);
Mask0 = Mask;
% use only the first 2 echos to estimate the initial field
if 0 %abs((TE(2)-TE(1))-(TE(3)-TE(2)))< 0.0002
    [iFreq_raw N_std] = Fit_ppm_complex(iField,TE);
else
    %[iFreq_raw N_std,iFreq_ref] = Fit_ppm_complex_TE(iField,TE);
    method = 2;  
    [iFreq_raw,N_std,iFreq_ref] = Fit_ppm_complex_test(iField,TE,Mask,method);
    %[iFreq_raw,N_std] = Fit_ppm_complex_Laplacian(iField,TE,Mask);
end
if ndims(iFreq_raw)<3
    iFreq_raw(:,:,2) = iFreq_raw;
    iMag(:,:,2) = iMag;
end    
%iFreq_raw = iFreq_raw.*Mask;
iFreq_raw(isnan(iFreq_raw))=0;
iFreq_raw(isinf(iFreq_raw))=0;
iFreq_raw0 = iFreq_raw(:,:,:);
iFreq_raw1 = iFreq_raw(:,:,:);

if SUBSAMPLE == 2
    iField = iField(1:SUBSAMPLE:end,1:SUBSAMPLE:end,:,:);
    iFreq_raw1 = iFreq_raw0(1:SUBSAMPLE:end,1:SUBSAMPLE:end,:);
    N_std = N_std(1:SUBSAMPLE:end,1:SUBSAMPLE:end,:);
    Mask = Mask(1:SUBSAMPLE:end,1:SUBSAMPLE:end,:);
    iMag = iMag(1:SUBSAMPLE:end,1:SUBSAMPLE:end,:);
    voxel_size(1) = voxel_size(1)/0.5;
    voxel_size(2) = voxel_size(2)/0.5;
end

%iMag = sqrt(sum(abs(iField).^2,4));

%Mask=autoMask(iMag,voxel_size);

delta_TE = (TE(2) - TE(1))*1e-6;

if nargin < 7
    dfat = -3.5e-6*CF;
end
    dyna_range = 1/delta_TE;
    effect_fat_Hz = dfat + floor( (0.5*dyna_range-dfat)/dyna_range)*dyna_range;
    effect_fat_rad = effect_fat_Hz/dyna_range*2*pi;


p = 2;
w1 = effect_fat_rad/pi;
z = size(iFreq_raw1,3);
w1
if (w1 > 0)
    w = w1;
    [unwphw,iter,erglist] = phase_unwrap_3d_UNIC(iFreq_raw1,p,iMag,voxel_size,Mask);
    %[unwphw,iter,erglist] = phase_unwrap_3d_UNIC_1002(iFreq_raw1,p,iMag,voxel_size,Mask); 
    %[unwphw,iter,erglist] = phase_unwrap_3d_UNIC_Chris(iFreq_raw1,p,iMag,voxel_size,Mask); 
    energy = erglist;
    if L_Fat == 0
        wunwph_uf = unwphw;
    else
        %[wkappa,wm_fat,wunwph_uf,iter,erglist,wkiter] = unwrap_unfat_3dP_UNIC_Chris_1002(voxel_size,iMag,w,unwphw,p); 
        %[wkappa,wm_fat,wunwph_uf,iter,erglist,wkiter] = unwrap_unfat_3dP_UNIC_Chris(voxel_size,iMag,w,unwphw,p); 
        [wkappa,wm_fat,wunwph_uf,iter,erglist,wkiter] = unwrap_unfat_3dP_UNIC_Chris_1017(voxel_size,iMag,w,unwphw,p,Mask); % a small mistake found here. 2014.2.21
    end
    energy = [energy erglist];
end


if (w1 < 0)
    w = -w1;
    [unwphw,iter,erglist] = phase_unwrap_3d_UNIC_Chris(iFreq_raw1,p,iMag,voxel_size,Mask);
    %[unwphw,iter,erglist] = phase_unwrap_3d_UNIC_Chris_1002(iFreq_raw1,p,iMag,voxel_size,Mask);  
    %[unwphw,iter,erglist] = phase_unwrap_3d_UNIC(iFreq_raw1,p,iMag,voxel_size,Mask); 
    energy = erglist;
    if L_Fat == 0
        wunwph_uf = unwphw;
    else
        %[wkappa,wm_fat,wunwph_uf,iter,erglist,wkiter] = unwrap_unfat_3dN(voxel_size,iMag,w,unwphw,p);
        [wkappa,wm_fat,wunwph_uf,iter,erglist,wkiter] = unwrap_unfat_3dN_UNIC_Chris_1017(voxel_size,iMag,w,unwphw,p,Mask);
        %[wkappa,wm_fat,wunwph_uf,iter,erglist,wkiter] = unwrap_unfat_3dN_UNIC_Chris_1002(voxel_size,iMag,w,unwphw,p,Mask);
    end
    energy = [energy erglist];
end

% interpolation to get the field map
%{
if SUBSAMPLE == 2
    allX = 1:sx;
    allY = 1:sy;
    subX = 1:SUBSAMPLE:sx;
    subY = 1:SUBSAMPLE:sy;
    [ALLX,ALLY] = meshgrid(allY(:),allX(:));
    [SUBX,SUBY] = meshgrid(subY(:),subX(:));
    fm = zeros(sx,sy,sz);
    for ind = 1:sz
        fm(:,:,ind) = interp2(SUBX,SUBY,wunwph_uf(:,:,ind),ALLX,ALLY,'*spline');
    end 
    
    k_vals = unique(wkappa);
    k_app = (fm - iFreq_raw0)/pi;

    for i = 1:length(k_vals)
        dk(:,:,:,i) = abs(k_app - k_vals(i));
    end 
    [a b] = min(dk,[],4);

    K = b;
    for i = 1:length(k_vals)
        K(b==i)= k_vals(i);
    end
    wunwph_uf = iFreq_raw + K*pi;
 end
%}
if SUBSAMPLE == 2
wunwph_uf = imresize3(wunwph_uf,[x,y,z],'linear');
unwphw = imresize3(unwphw,[x,y,z], 'linear');
N_std = imresize3(N_std,[x,y,z],'linear');
iFreq_raw = imresize3(iFreq_raw,[x,y,z],'linear');
iFreq_ref = imresize3(iFreq_ref,[x,y,z],'linear');
end
iField = iField0;
% 


