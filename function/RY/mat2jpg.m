function [feedbak] =mat2jpg(im_Nonbin,save_fname)%
%input matrix and directory name

mkdir(save_fname);
counter = size(im_Nonbin,3);
 im_Nonbin=(im_Nonbin-min(im_Nonbin(:)))/(max(im_Nonbin(:))-min(im_Nonbin(:)));
for n = 1:counter,  
    if(n<10)
        imwrite(squeeze(im_Nonbin(:,:,n)),strcat(save_fname,'00',num2str(n),'.jpeg'),'jpeg','Quality',100);
    elseif (n < 99)
        imwrite(squeeze(im_Nonbin(:,:,n)),strcat(save_fname,'0',num2str(n),'.jpeg'),'jpeg','Quality',100);
    else
        imwrite(squeeze(im_Nonbin(:,:,n)),strcat(save_fname,num2str(n),'.jpeg'),'jpeg','Quality',100);
    end
end
feedbak=true;