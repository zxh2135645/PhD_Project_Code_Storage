function [feedbak] =mat2tif(im_Nonbin,save_fname)%
%input matrix and directory name
if ~exist(save_fname, 'dir')
    mkdir(save_fname);
end
counter = size(im_Nonbin,3);
 im_Nonbin=(im_Nonbin-min(im_Nonbin(:)))/(max(im_Nonbin(:))-min(im_Nonbin(:)));
for n = 1:counter,  
    if(n<10)
        imwrite(squeeze(im_Nonbin(:,:,n)),strcat(save_fname,'00',num2str(n),'.tif'),'tif');
    elseif (n < 99)
        imwrite(squeeze(im_Nonbin(:,:,n)),strcat(save_fname,'0',num2str(n),'.tif'),'tif');
    else
        imwrite(squeeze(im_Nonbin(:,:,n)),strcat(save_fname,num2str(n),'.tif'),'tif');
    end
end
feedbak=true;