clear all;
close all;
%%
load('./mri.mat');
figure(); 
subplot(2,2,1);
imagesc(I);

I_fft = fftshift(fftshift(fft2(I),1),2);

subplot(2,2,2); imagesc(abs(I_fft));

I_fft(1:2:end,:) = 0;

subplot(2,2,3); imagesc(abs(I_fft));

I_fourier = ifftshift(ifft2(I_fft),1);

subplot(2,2,4); imagesc(abs(I_fourier));

%% Direct interpolation does not work!
x = 1:size(I_fft,1);
x_o = (1:2:length(I_fft)).';
y = size(I_fft,2);

I_fft_interp = zeros(size(I_fft));
for i = 1:y
    I_ffty = I_fft(:, i);
    I_ffty(I_ffty == 0) = [];
    
    I_ffty_interp = interp1(x_o, I_ffty, x);
    I_fft_interp(:,i) = I_ffty_interp;
end

figure(); imagesc(abs(I_fft_interp));
I_fourier_interp = ifft2(I_fft_interp);
figure(); imagesc(abs(I_fourier_interp));

%% This may work
I_fft = fftshift(fftshift(fft2(I),1),2);
I_fft_expand = zeros(size(I_fft,1)*2, size(I_fft,2));
I_fft_expand(1:2:end,:) = I_fft;

I_fourier_expand = ifftshift(ifft2(I_fft_expand),1);

figure(); 
subplot(2,2,1);
imagesc(abs(I_fft_expand)); axis image;
subplot(2,2,2);
imagesc(abs(I_fourier_expand)); axis image;

DC = size(I_fourier_expand,1)/2;
I_fft_shrink = I_fft_expand(DC-DC/2+1:DC+DC/2,:);
I_fourier_shrink = ifftshift(ifft2(I_fft_shrink),1);

subplot(2,2,3);
imagesc(abs(I_fft_shrink)); axis image;
subplot(2,2,4);
imagesc(abs(I_fourier_shrink)); axis image;