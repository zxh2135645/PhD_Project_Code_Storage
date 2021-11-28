function [feedbak] =mat2dicom(im_Nonbin,save_fname)%,Dicominfo,onlinedirbottom)%
%input matrix and directory name (filesep at the end)
%if write dincominfo then find the Dicom
%sourcefile directory to copy dicominfo

mkdir(save_fname);
counter = size(im_Nonbin,3);
%im_Nonbin=double((im_Nonbin-min(im_Nonbin(:)))/(max(im_Nonbin(:))-min(im_Nonbin(:))));
% im_Nonbin=im_Nonbin*65535;
 im_Nonbin=im_Nonbin;
for n = 1:counter,  
    if(n<10)
        dicomwrite(squeeze(im_Nonbin(:,:,n)),strcat(save_fname,'00',num2str(n),'.ima'),'CompressionMode','None','Quality',100);
    elseif (n < 99)
        dicomwrite(squeeze(im_Nonbin(:,:,n)),strcat(save_fname,'0',num2str(n),'.ima'),'CompressionMode','None','Quality',100);
    else
        dicomwrite(squeeze(im_Nonbin(:,:,n)),strcat(save_fname,num2str(n),'.ima'),'CompressionMode','None','Quality',100);
    end
%     if Dicominfo
%        if(n<10)
%             targetfile=strcat(save_fname,'00',num2str(n),'.ima');
%        elseif (n < 99)
%            targetfile=strcat(save_fname,'0',num2str(n),'.ima');
%        else
%            targetfile=strcat(save_fname,num2str(n),'.ima');
%        end
%             listonlinebottom=dir(onlinedirbottom);
%             oversampleslice=0;
%            if counter>(size(listonlinebottom,1)-2)
%                oversampleslice=counter/2-(size(listonlinebottom,1)-2)/2;
%            end
%             if ((n+2)<=(oversampleslice+2))
%                sourcefile=[onlinedirbottom,'\',listonlinebottom(3).name]; 
%             elseif((n+2)>size(listonlinebottom,1))
%             sourcefile=[onlinedirbottom,'\',listonlinebottom(size(listonlinebottom,1)).name];
%             else
%                sourcefile=[onlinedirbottom,'\',listonlinebottom(n+2-oversampleslice).name]; 
%             end
%             info = dicominfo(sourcefile);
%             
%                 I = dicomread(targetfile);
%                    dicomwrite(I,targetfile,info)
%         
%     end
end
feedbak=true;