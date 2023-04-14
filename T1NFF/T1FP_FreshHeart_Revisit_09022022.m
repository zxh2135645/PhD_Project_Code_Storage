clear all;
close all;

ff = zeros(size(fat));
fat_flag = fat > water;
ff(fat_flag) = abs(fat(fat_flag)) ./ abs(fat(fat_flag) + water(fat_flag));
ff(~fat_flag) = 1 - (abs(water(~fat_flag))) ./ abs(fat(~fat_flag) + water(~fat_flag));

figure(); 
% subplot(2,2,1); imagesc(transpose(ff(:,:,1)));caxis([0 0.2]);
% subplot(2,2,2); imagesc(transpose(ff(:,:,2)));caxis([0 0.2]);
% subplot(2,2,3); imagesc(transpose(ff(:,:,3)));caxis([0 0.2]);
% subplot(2,2,4); imagesc(transpose(ff(:,:,4)));caxis([0 0.2]);
subplot(2,2,1); imagesc(ff(:,:,1));caxis([0 0.5]); axis image;
subplot(2,2,2); imagesc(ff(:,:,2));caxis([0 0.5]); axis image;
subplot(2,2,3); imagesc(ff(:,:,3));caxis([0 0.5]); axis image;
subplot(2,2,4); imagesc(ff(:,:,4));caxis([0 0.5]); axis image;

figure(); 
subplot(2,2,1); imagesc(transpose(R2s(:,:,1)));caxis([0 50]);
subplot(2,2,2); imagesc(transpose(R2s(:,:,2)));caxis([0 50]);
subplot(2,2,3); imagesc(transpose(R2s(:,:,3)));caxis([0 50]);
subplot(2,2,4); imagesc(transpose(R2s(:,:,4)));caxis([0 50]);


