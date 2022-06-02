%Please cite the article 
%Prakash, R. M., & Kumari, R. S. (2016). Spatial Fuzzy C Means and 
%Expectation Maximization Algorithms with Bias Correction for Segmentation
%of MR Brain Images. Journal of Medical Systems, 41(1). 

clear all;
tic
global mu  ima original mask r mix pp r gk maskg mug bias smooth bias2 bias3 s biasout
[q,i]=IBSRread('C:\Users\yanghj\Downloads\SFCM_bias.zip\SFCM_bias\20Normals_T1_seg\','C:\Users\yanghj\Downloads\SFCM_bias.zip\SFCM_bias\20Normals_T1_brain\',9);
image=i;
output=zeros(size(image));
cf=zeros(4,4,49);
gk1=[1/9 1/9 1/9;1/9 1/9 1/9;1/9 1/9 1/9;];
%%% Wavelet Transform Parameters.
for y=1:49
%y=60;
        x=image(:,:,y);
        original=x;
       [maskg mug vg pg]=EMSeg(original,7);
        copy=x;
        original=x;
         im=q(:,:,y);
    
        k=3;
       
     
   [cc mask]=fuzzybiassmooth(x,k);
   [m n]=size(x);
  
    
    
  for i=1:m
    for j=1:n
          if mask(i,j)==1
            if maskg(i,j)==2
            mask(i,j)=128;
           end
            if maskg(i,j)~=2
                mask(i,j)=192;
            end
        end
         if mask(i,j)==2
            mask(i,j)=192;
         end
         if mask(i,j)==3
            mask(i,j)=254;
         end
      
    end
end
output(:,:,y)=mask;
output=uint8(output);
    manualsegmented=reshape(q(:,:,y),1,m*n);
actualsegmented=reshape(output(:,:,y),1,m*n);
disp(y);
cf(:,:,y) = cfmatrix(manualsegmented,actualsegmented,[0,128,192,254],0);
end
figure(1);
subplot(3,3,1);imshow(uint8(output(:,:,1)));subplot(3,3,2);imshow(uint8(output(:,:,2)));subplot(3,3,3);imshow(uint8(output(:,:,3)));
subplot(3,3,4);imshow(uint8(output(:,:,4)));subplot(3,3,5);imshow(uint8(output(:,:,5)));subplot(3,3,6);imshow(uint8(output(:,:,6)));
subplot(3,3,7);imshow(uint8(output(:,:,7)));subplot(3,3,8);imshow(uint8(output(:,:,8)));subplot(3,3,9);imshow(uint8(output(:,:,9)));
figure(2);
subplot(3,3,1);imshow(uint8(output(:,:,10)));subplot(3,3,2);imshow(uint8(output(:,:,11)));subplot(3,3,3);imshow(uint8(output(:,:,12)));
subplot(3,3,4);imshow(uint8(output(:,:,13)));subplot(3,3,5);imshow(uint8(output(:,:,14)));subplot(3,3,6);imshow(uint8(output(:,:,15)));
subplot(3,3,7);imshow(uint8(output(:,:,16)));subplot(3,3,8);imshow(uint8(output(:,:,17)));subplot(3,3,9);imshow(uint8(output(:,:,18)));
figure(3);
subplot(3,3,1);imshow(uint8(output(:,:,19)));subplot(3,3,2);imshow(uint8(output(:,:,20)));subplot(3,3,3);imshow(uint8(output(:,:,21)));
subplot(3,3,4);imshow(uint8(output(:,:,22)));subplot(3,3,5);imshow(uint8(output(:,:,23)));subplot(3,3,6);imshow(uint8(output(:,:,24)));
subplot(3,3,7);imshow(uint8(output(:,:,25)));subplot(3,3,8);imshow(uint8(output(:,:,26)));subplot(3,3,9);imshow(uint8(output(:,:,27)));
figure(4);
subplot(3,3,1);imshow(uint8(output(:,:,28)));subplot(3,3,2);imshow(uint8(output(:,:,29)));subplot(3,3,3);imshow(uint8(output(:,:,30)));
subplot(3,3,4);imshow(uint8(output(:,:,31)));subplot(3,3,5);imshow(uint8(output(:,:,32)));subplot(3,3,6);imshow(uint8(output(:,:,33)));
subplot(3,3,7);imshow(uint8(output(:,:,34)));subplot(3,3,8);imshow(uint8(output(:,:,35)));subplot(3,3,9);imshow(uint8(output(:,:,36)));
figure(5);
subplot(3,3,1);imshow(uint8(output(:,:,37)));subplot(3,3,2);imshow(uint8(output(:,:,38)));subplot(3,3,3);imshow(uint8(output(:,:,39)));
subplot(3,3,4);imshow(uint8(output(:,:,40)));subplot(3,3,5);imshow(uint8(output(:,:,41)));subplot(3,3,6);imshow(uint8(output(:,:,42)));
subplot(3,3,7);imshow(uint8(output(:,:,43)));subplot(3,3,8);imshow(uint8(output(:,:,44)));subplot(3,3,9);imshow(uint8(output(:,:,45)));
figure(6);
imshow(biasout,[]);
toc
cfresult=zeros(4,4);
for i=1:49
    cfresult=cfresult+cf(:,:,i);
end
disp(cfresult);
tanimotocsf=cfresult(2,2)/(cfresult(2,2)+cfresult(2,1)+cfresult(2,3)+cfresult(2,4)+cfresult(1,2)+cfresult(3,2)+cfresult(4,2));
tanimotogray=cfresult(3,3)/(cfresult(3,3)+cfresult(3,1)+cfresult(3,2)+cfresult(3,4)+cfresult(1,3)+cfresult(2,3)+cfresult(4,3));
tanimotowhite=cfresult(4,4)/(cfresult(4,4)+cfresult(4,1)+cfresult(4,3)+cfresult(4,2)+cfresult(1,4)+cfresult(3,4)+cfresult(2,4));
disp(tanimotocsf);
disp(tanimotogray);
disp(tanimotowhite);
dicsf=(2*cfresult(2,2))/((2*cfresult(2,2))+cfresult(2,1)+cfresult(2,3)+cfresult(2,4)+cfresult(1,2)+cfresult(3,2)+cfresult(4,2));
digray=(2*cfresult(3,3))/((2*cfresult(3,3))+cfresult(3,1)+cfresult(3,2)+cfresult(3,4)+cfresult(1,3)+cfresult(2,3)+cfresult(4,3));
diwhite=(2*cfresult(4,4))/((2*cfresult(4,4))+cfresult(4,1)+cfresult(4,3)+cfresult(4,2)+cfresult(1,4)+cfresult(3,4)+cfresult(2,4));
disp(dicsf);
disp(digray);
disp(diwhite);






