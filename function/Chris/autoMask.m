function [Mask]=autoMask(iMag,voxel_size)
% tultithreshold
maskpercent=multithresh(iMag/(max(iMag(:))),4);
    Mask=iMag>(max(iMag(:))*min(maskpercent));%generate a rough mask
    %open and close to avoid holes
       % Mask= bwareaopen(Mask,100,6);
    
%closed
disksize= 15/voxel_size(1);%close with 15mm disk
se = strel('disk',4);
      Mask=imclose(Mask,se);
