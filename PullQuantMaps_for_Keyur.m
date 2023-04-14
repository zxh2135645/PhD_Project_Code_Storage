close all;
clear all;
%%
addpath('function/')
% Sahara Chronic for Keyur
% Merry  Chronic for Keyur
ff = abs(fat) ./ (abs(fat) + abs(water));
figure();
for i = 1:size(fat, 3)
    subplot(2,2,i); imagesc(abs(ff(:,:,i))); caxis([0 0.5]); axis off; axis image;
    colorbar;
    title(cat(2, 'Slice ', num2str(i)));
end
colormap(brewermap([], '*PuBu'));
%%
figure();
for i = 1:size(R2s, 3)
    subplot(2,2,i); imagesc(abs(R2s(:,:,i))); caxis([0 100]); axis off; axis image;
    colorbar;
    title(cat(2, 'Slice ', num2str(i)));
end
colormap(brewermap([], '*RdBu'));

%%
T2s = 1000 ./ R2s;
figure();
for i = 1:size(R2s, 3)
    subplot(2,2,i); imagesc(abs(T2s(:,:,i))); caxis([0 50]); axis off; axis image;
    colorbar;
end
colormap(brewermap([], '*RdBu'));

%% Fat Fraction
ff = zeros(size(fat));
fat_flag = fat > water;
ff(fat_flag) = abs(fat(fat_flag)) ./ abs(fat(fat_flag) + water(fat_flag));
ff(~fat_flag) = 1 - (abs(water(~fat_flag))) ./ abs(fat(~fat_flag) + water(~fat_flag));
i = 4;
figure();
imagesc(abs(ff(:,:,i))); caxis([0 0.5]); axis off; axis image;
colorbar;
%title(cat(2, 'Slice ', num2str(i)));
colormap(brewermap([], '*PuBu'));

figure();
imagesc(abs(R2s(:,:,i))); caxis([0 100]); axis off; axis image;
colorbar;
%title(cat(2, 'Slice ', num2str(i)));
colormap(brewermap([], '*RdBu'));

T2s = 1000 ./ R2s;
figure();
imagesc(abs(T2s(:,:,i))); caxis([0 50]); axis off; axis image;
colorbar;
%title(cat(2, 'Slice ', num2str(i)));
colormap(brewermap([], '*RdBu'));