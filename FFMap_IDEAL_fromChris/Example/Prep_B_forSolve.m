%Prep_B_forSolve

mask=Mask;

%clearvars -except DCf stdf slice0
DC_limit0 = [3.5];
%% import
%B0=B0_test;
%Bz=Bz_30slice;
if IS_Slice
    slab=ShimSlice;
else
    slab = 1:size(B0,3); %whole heart
end
mask(isnan(B0))=0;
Bz(isnan(Bz))=0;
%% 
clear B0f Bzf
 [nx,ny,nz,nc] = size(Bz);
 for i=1:size(Bz,4)
 Bz_temp = Bz(:,:,:,i);
 Bzf(:,i)=Bz_temp(mask==1);
 end
 clear Bz_temp;


 B0f = B0(mask==1);
