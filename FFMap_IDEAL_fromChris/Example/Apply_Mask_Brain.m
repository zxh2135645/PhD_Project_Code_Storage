
%% Creat Mask
% Creat mask

AllPhasemap(n).Mask=BET(iMag,matrix_size,voxel_size);
%one more time to trim the FOV
AllPhasemap(n).Mask(:,:,20:end)=0;%!!!!hard code for 1 time use 
[Mi,Mj,Mk]=ind2sub(size(AllPhasemap(n).Mask),find(AllPhasemap(n).Mask));
Maskind=[Mi,Mj,Mk];
tempFOVstart=[min(Mi) min(Mj) min(Mk)];
tempFOVend=[max(Mi) max(Mj) max(Mk)];
tempiMag=iMag(tempFOVstart(1):tempFOVend(1), tempFOVstart(2):tempFOVend(2), tempFOVstart(3):tempFOVend(3));
tempMask=BET(tempiMag,size(tempiMag),voxel_size);
Mask=zeros(size(iMag));
Mask(tempFOVstart(1):tempFOVend(1), tempFOVstart(2):tempFOVend(2), tempFOVstart(3):tempFOVend(3))=tempMask;
AllPhasemap(n).Mask=Mask;
%save([MRdat_path,'Manual_Mask.mat'],'Mask');
%else
%load([MRdat_path,'Manual_Mask.mat']);
%AllPhasemap(n).Mask=Mask;
%end
%clear Mask