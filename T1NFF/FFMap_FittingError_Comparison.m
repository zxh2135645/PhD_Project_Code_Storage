clear all;
close all;
%%
addpath('../function/');
dicom_dir = uigetdir();

resdir=GetFullPath([dicom_dir, '/..','/Result/']);
if ~exist(resdir)
    mkdir(resdir);
end
strings = strsplit(dicom_dir, '/');
subject_name = strings{end};
fib_struct = load([resdir,'AllPhasemap_', subject_name, '_MP_fibrosis.mat']);
sp_struct = load([resdir,'AllPhasemap_', subject_name, '_SinglePeak.mat']);
mp_struct = load([resdir,'AllPhasemap_', subject_name, '_MultiPeak.mat']);

figure(); 
subplot(2,2,1); imagesc(fib_struct.fitting_error(:,:,1));caxis([0 0.2]);
subplot(2,2,2); imagesc(fib_struct.fitting_error(:,:,2));caxis([0 0.2]);
subplot(2,2,3); imagesc(fib_struct.fitting_error(:,:,3));caxis([0 0.2]);
subplot(2,2,4); imagesc(fib_struct.fitting_error(:,:,4));caxis([0 0.2]);

figure(); 
subplot(2,2,1); imagesc(sp_struct.fitting_error(:,:,1));caxis([0 0.2]);
subplot(2,2,2); imagesc(sp_struct.fitting_error(:,:,2));caxis([0 0.2]);
subplot(2,2,3); imagesc(sp_struct.fitting_error(:,:,3));caxis([0 0.2]);
subplot(2,2,4); imagesc(sp_struct.fitting_error(:,:,4));caxis([0 0.2]);

figure(); 
subplot(2,2,1); imagesc(mp_struct.fitting_error(:,:,1));caxis([0 0.2]);
subplot(2,2,2); imagesc(mp_struct.fitting_error(:,:,2));caxis([0 0.2]);
subplot(2,2,3); imagesc(mp_struct.fitting_error(:,:,3));caxis([0 0.2]);
subplot(2,2,4); imagesc(mp_struct.fitting_error(:,:,4));caxis([0 0.2]);

%%
mask_f = cat(2, dicom_dir, '/../mask_roi.mat');
load(mask_f);
fitting_errors_fib = mean(nonzeros(fib_struct.fitting_error .* mask_roi), 'omitnan');
fitting_errors_sp = mean(nonzeros(sp_struct.fitting_error .* mask_roi), 'omitnan');
fitting_errors_mp = mean(nonzeros(mp_struct.fitting_error .* mask_roi), 'omitnan');
fitting_errors_fib_sd = std(nonzeros(fib_struct.fitting_error .* mask_roi), 'omitnan');
fitting_errors_sp_sd = std(nonzeros(sp_struct.fitting_error .* mask_roi), 'omitnan');
fitting_errors_mp_sd = std(nonzeros(mp_struct.fitting_error .* mask_roi), 'omitnan');
x = 1:3;
errhigh = [fitting_errors_fib+fitting_errors_fib_sd/sqrt(length(x)), fitting_errors_sp + fitting_errors_sp_sd/sqrt(length(x)), fitting_errors_mp+fitting_errors_mp_sd/sqrt(length(x))];
errlow = [fitting_errors_fib-fitting_errors_fib_sd/sqrt(length(x)), fitting_errors_sp-fitting_errors_sp_sd/sqrt(length(x)), fitting_errors_mp-fitting_errors_mp_sd/sqrt(length(x))];

figure(); bar(x, [fitting_errors_fib, fitting_errors_sp, fitting_errors_mp]);
hold on;
er = errorbar(x,[fitting_errors_fib, fitting_errors_sp, fitting_errors_mp], errlow, errhigh);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
xticklabels({'Fib','SP','MP'}); ylabel('Fitting Error');
hold off

%% 
fat = fib_struct.fat;
water = fib_struct.water;
fib = fib_struct.fib;
R2s = fib_struct.R2s;

ff = zeros(size(fat));
fat_flag = fat > (water+fib);
ff(fat_flag) = abs(fat(fat_flag)) ./ abs(fat(fat_flag) + water(fat_flag) + fib(fat_flag));
ff(~fat_flag) = 1 - (abs(water(~fat_flag))+abs(fib(~fat_flag))) ./ abs(fat(~fat_flag) + water(~fat_flag) + fib(~fat_flag));

figure(); subplot(2,2,1); imagesc(ff(:,:,1));caxis([0 0.2]);
subplot(2,2,2); imagesc(ff(:,:,2));caxis([0 0.2]);
subplot(2,2,3); imagesc(ff(:,:,3));caxis([0 0.2]);
subplot(2,2,4); imagesc(ff(:,:,4));caxis([0 0.2]);

fibf = zeros(size(fib));
fib_flag = fib > (water+fat);
fibf(fib_flag) = abs(fib(fib_flag)) ./ abs(fat(fib_flag) + water(fib_flag) + fib(fib_flag));
fibf(~fib_flag) = 1 - (abs(water(~fib_flag)+ fib(~fib_flag))) ./ abs(fat(~fib_flag) + water(~fib_flag) + fib(~fib_flag));

figure(); subplot(2,2,1); imagesc(fibf(:,:,1));caxis([0 0.2]);
subplot(2,2,2); imagesc(fibf(:,:,2));caxis([0 0.2]);
subplot(2,2,3); imagesc(fibf(:,:,3));caxis([0 0.2]);
subplot(2,2,4); imagesc(fibf(:,:,4));caxis([0 0.2]);


waterf = zeros(size(water));
water_flag = water > (fib+fat);
waterf(water_flag) = abs(water(water_flag)) ./ abs(fat(water_flag) + water(water_flag) + fib(water_flag));
waterf(~water_flag) = 1 - (abs(fat(~water_flag)+ fib(~water_flag))) ./ abs(fat(~water_flag) + water(~water_flag) + fib(~water_flag));

figure(); 
subplot(2,2,1); imagesc(waterf(:,:,1));caxis([0 0.2]);
subplot(2,2,2); imagesc(waterf(:,:,2));caxis([0 1]);
subplot(2,2,3); imagesc(waterf(:,:,3));caxis([0 1]);
subplot(2,2,4); imagesc(waterf(:,:,4));caxis([0 1]);

figure(); 
subplot(2,2,1); imagesc(R2s(:,:,1));caxis([0 100]);
subplot(2,2,2); imagesc(R2s(:,:,2));caxis([0 100]);
subplot(2,2,3); imagesc(R2s(:,:,3));caxis([0 100]);
subplot(2,2,4); imagesc(R2s(:,:,4));caxis([0 100]);

%%  Single-Peak
fat = sp_struct.fat;
water = sp_struct.water;
R2s = sp_struct.R2s;

ff = zeros(size(fat));
fat_flag = fat > (water);
ff(fat_flag) = abs(fat(fat_flag)) ./ abs(fat(fat_flag) + water(fat_flag));
ff(~fat_flag) = 1 - (abs(water(~fat_flag))) ./ abs(fat(~fat_flag) + water(~fat_flag));

figure(); 
subplot(2,2,1); imagesc(ff(:,:,1));caxis([0 0.2]);
subplot(2,2,2); imagesc(ff(:,:,2));caxis([0 0.2]);
subplot(2,2,3); imagesc(ff(:,:,3));caxis([0 0.2]);
subplot(2,2,4); imagesc(ff(:,:,4));caxis([0 0.2]);

waterf = zeros(size(water));
water_flag = water > (fat);
waterf(water_flag) = abs(water(water_flag)) ./ abs(fat(water_flag) + water(water_flag));
waterf(~water_flag) = 1 - (abs(fat(~water_flag))) ./ abs(fat(~water_flag) + water(~water_flag));

figure(); 
subplot(2,2,1); imagesc(waterf(:,:,1));caxis([0 0.2]);
subplot(2,2,2); imagesc(waterf(:,:,2));caxis([0 1]);
subplot(2,2,3); imagesc(waterf(:,:,3));caxis([0 1]);
subplot(2,2,4); imagesc(waterf(:,:,4));caxis([0 1]);

figure(); 
subplot(2,2,1); imagesc(R2s(:,:,1));caxis([0 100]);
subplot(2,2,2); imagesc(R2s(:,:,2));caxis([0 100]);
subplot(2,2,3); imagesc(R2s(:,:,3));caxis([0 100]);
subplot(2,2,4); imagesc(R2s(:,:,4));caxis([0 100]);

%%  Multi-Peak
fat = mp_struct.fat;
water = mp_struct.water;
R2s = mp_struct.R2s;

ff = zeros(size(fat));
fat_flag = fat > (water);
ff(fat_flag) = abs(fat(fat_flag)) ./ abs(fat(fat_flag) + water(fat_flag));
ff(~fat_flag) = 1 - (abs(water(~fat_flag))) ./ abs(fat(~fat_flag) + water(~fat_flag));

figure(); 
subplot(2,2,1); imagesc(ff(:,:,1));caxis([0 0.2]);
subplot(2,2,2); imagesc(ff(:,:,2));caxis([0 0.2]);
subplot(2,2,3); imagesc(ff(:,:,3));caxis([0 0.2]);
subplot(2,2,4); imagesc(ff(:,:,4));caxis([0 0.2]);

waterf = zeros(size(water));
water_flag = water > (fat);
waterf(water_flag) = abs(water(water_flag)) ./ abs(fat(water_flag) + water(water_flag));
waterf(~water_flag) = 1 - (abs(fat(~water_flag))) ./ abs(fat(~water_flag) + water(~water_flag));

figure(); 
subplot(2,2,1); imagesc(waterf(:,:,1));caxis([0 0.2]);
subplot(2,2,2); imagesc(waterf(:,:,2));caxis([0 1]);
subplot(2,2,3); imagesc(waterf(:,:,3));caxis([0 1]);
subplot(2,2,4); imagesc(waterf(:,:,4));caxis([0 1]);

figure(); 
subplot(2,2,1); imagesc(R2s(:,:,1));caxis([0 100]);
subplot(2,2,2); imagesc(R2s(:,:,2));caxis([0 100]);
subplot(2,2,3); imagesc(R2s(:,:,3));caxis([0 100]);
subplot(2,2,4); imagesc(R2s(:,:,4));caxis([0 100]);