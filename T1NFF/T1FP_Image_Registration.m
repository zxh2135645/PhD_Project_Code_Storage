clear all;
close all;
%% Try different registration methods
% img2(img2>100) = 100;
% img2(img2<0) = 0;
% I1 = img2; I2 = img;

I1 = moving; I2 = fixed;
% Set static and moving image
S=I2; M=I1;

% [Ireg,Bx,By,Fx,Fy] = register_images(M,S);
% Alpha (noise) constant
alpha=2.5;

% Velocity field smoothing kernel
Hsmooth=fspecial('gaussian',[60 60],10);

% The transformation fields
Tx=zeros(size(M)); Ty=zeros(size(M));

[Sy,Sx] = gradient(S);
for itt=1:200
	    % Difference image between moving and static image
        Idiff=M-S;

        % Default demon force, (Thirion 1998)
        %Ux = -(Idiff.*Sx)./((Sx.^2+Sy.^2)+Idiff.^2);
        %Uy = -(Idiff.*Sy)./((Sx.^2+Sy.^2)+Idiff.^2);

        % Extended demon force. With forces from the gradients from both
        % moving as static image. (Cachier 1999, He Wang 2005)
        [My,Mx] = gradient(M);
        Ux = -Idiff.*  ((Sx./((Sx.^2+Sy.^2)+alpha^2*Idiff.^2))+(Mx./((Mx.^2+My.^2)+alpha^2*Idiff.^2)));
        Uy = -Idiff.*  ((Sy./((Sx.^2+Sy.^2)+alpha^2*Idiff.^2))+(My./((Mx.^2+My.^2)+alpha^2*Idiff.^2)));
 
        % When divided by zero
        Ux(isnan(Ux))=0; Uy(isnan(Uy))=0;

        % Smooth the transformation field
        Uxs=3*imfilter(Ux,Hsmooth);
        Uys=3*imfilter(Uy,Hsmooth);

        % Add the new transformation field to the total transformation field.
        Tx=Tx+Uxs;
        Ty=Ty+Uys;
        M=movepixels(I1,Tx,Ty); 
end

M_bi = M>0.5;

figure();
subplot(1,3,1), imshow(I1,[]); title('image 1');
subplot(1,3,2), imshow(I2,[]); title('image 2');
subplot(1,3,3), imshow(M_bi,[]); title('Registered image 1');

figure();
subplot(1,2,1);
imshowpair(S,M_bi,'Scaling','joint');
subplot(1,2,2); imagesc(S&M_bi); axis image;

%% register_images
I1 = moving; I2 = fixed;
% Set static and moving image
S=I2; M=I1;

[Ireg,Bx,By,Fx,Fy] = register_images(M,S);

M_transform=movepixels(img2,Bx,By);

Ireg_bi = Ireg > 0.5;
M_transform_bi = M_transform;
figure();
subplot(1,4,1), imshow(I1,[]); title('image 1');
subplot(1,4,2), imshow(I2,[]); title('image 2');
subplot(1,4,3), imshow(Ireg_bi,[]); title('Registered image 1');
subplot(1,4,4), imshow(M_transform_bi,[]); title('Registered image 2'); caxis([0 50])

figure();
subplot(2,2,1);
imshowpair(S,Ireg_bi,'Scaling','joint');
subplot(2,2,2); imagesc(S&Ireg_bi); axis image;
subplot(2,2,3);
imshowpair(S,M_transform_bi,'Scaling','joint');
subplot(2,2,4); imagesc(S&M_transform_bi); axis image;