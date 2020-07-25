clear all;
close all;
% Image Pullup automation for exvivo scan which I'm expecting 3D acquistion
% for images

%% 20P10, 20P11 Pig scan

addpath('D:\src\function');

base_dir = uigetdir;
folder_glob = glob(cat(2, base_dir, '\*'));

labels = {'T1', 'T1MAP', 'T2', 'T2MAP', 'T2STAR', '_T2STAR', 'MGRE', 'T1_TSE', 'FATSAT', 'MT'};

label = labels{5};
idx_array = contains(folder_glob, label);

if any(idx_array)
    num = numel(nonzeros(idx_array));
    ind_array = find(idx_array == 1);
    dst_files = cell(num, 1);
    for i = 1:num
        dst_files{i} = folder_glob{ind_array(i)};
    end
    
    dst_name = ExtractNames(dst_files);
    disp(dst_name);
    
    sel_array = input('Please add an array here:  ');
    
    char_array = num2str(sel_array', '%04.f');
    
    ind_array2 = zeros(size(dst_name, 1), 1);
    for i = 1:size(char_array, 1)
        cha = char_array(i, :);
        ind_array2 = ind_array2 + contains(dst_name, cha);
    end
    
    ind_array3 = find(ind_array2 == 1);
    list_to_read = dst_files(ind_array3);
    
    name_to_compare = ExtractNames(list_to_read);
    
    order_to_read = zeros(length(list_to_read), 1);
    for i = 1:length(list_to_read)
        order_to_read(i) = find(contains(name_to_compare, char_array(i, :)) == 1);
    end
    
    save_dir = GetFullPath(cat(2, base_dir, '\..\img\'));
    if ~exist(save_dir, 'dir')
        mkdir(save_dir)
    end
    
    whatsinit = cell(length(list_to_read), 1);
    for i = 1:length(list_to_read)
        f = list_to_read{order_to_read(i)};
        whatsinit{i} = dicom23D(f);
    end
    
end

roi_dir = GetFullPath(cat(2, base_dir, '\..\roi2.mat'));
if exist(roi_dir, 'file')
    roi = load(roi_dir);
    
end
%% 
img = whatsinit{4};
coords_cell = cell(size(img,3), 2);

for i = 6:(size(img, 3)-2)
    disp('Epicardium: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 50])
    epi = drawpolygon(gca);
    epi_coords = epi.Position;
    
    disp('Endocardium: ');
    figure('Position', [100 0 1600 1600]); imagesc(img(:,:,i)); axis image;
    caxis([0 50])
    endo = drawpolygon(gca);
    endo_coords = endo.Position;
    
    coords_cell{i, 1} = epi.Position;
    coords_cell{i, 2} = endo.Position;
    epi_mask = createMask(epi);
    endo_mask = createMask(endo);
    
    myo_mask(:,:,i-5) = epi_mask - endo_mask;
end
%%
center_coords = DivideMyoInHalf(epi_coords, endo_coords);
figure('Position', [100 0 1600 1600]); imagesc(img(:,:,11)); axis image;  caxis([0 50])
centerline = drawpolygon(gca,'Position', [center_coords(:,1), center_coords(:,2)]);
center_mask = createMask(centerline);

figure('Position', [100 0 1600 1600]); imagesc(img2); axis image;  caxis([0 50])
centerline2 = drawpolygon(gca,'Position', [center_coords(:,1)/2+1/2, center_coords(:,2)/2+1/2]);
center_mask2 = createMask(centerline2);

figure('Position', [100 0 1600 1600]); imagesc(img3); axis image;  caxis([0 50])
centerline3 = drawpolygon(gca,'Position', [center_coords(:,1)/3+2/3, center_coords(:,2)/3+2/3]);
center_mask3 = createMask(centerline3);

figure('Position', [100 0 1600 1600]); imagesc(img4); axis image;  caxis([0 50])
centerline4 = drawpolygon(gca,'Position', [center_coords(:,1)/4+3/4, center_coords(:,2)/4+3/4]);
center_mask4 = createMask(centerline4);

myo_mask_endo = myo_mask(:,:,6) .* center_mask;
myo_mask_epi = myo_mask(:,:,6) - myo_mask_endo;

myo_mask_endo2 = myo_mask2 .* center_mask2;
myo_mask_epi2 = myo_mask2 - myo_mask_endo2;

myo_mask_endo3 = myo_mask3 .* center_mask3;
myo_mask_epi3 = myo_mask3 - myo_mask_endo3;

myo_mask_endo4 = myo_mask4 .* center_mask4;
myo_mask_epi4 = myo_mask4 - myo_mask_endo4;
%%
% figure();
% 
% for i = 1:size(myo_mask)
%     subplot(3,3,i)
%     imagesc(myo_mask(:,:,i).*img(:,:,i+5)); caxis([0 50]); 
%     axis image;
% end

img2 = whatsinit{3};
figure();
imagesc(img2(:,:,11)); caxis([0 50]);
epi2 = drawpolygon(gca,'Position', [coords_cell{11,1}(:,1)/2 + 1/2, coords_cell{11,1}(:,2)/2 + 1/2]);
endo2 = drawpolygon(gca,'Position', [coords_cell{11,2}(:,1)/2 + 1/2, coords_cell{11,2}(:,2)/2 + 1/2]);
epi_mask2 = createMask(epi2);
endo_mask2 = createMask(endo2);

img3 = whatsinit{2};
figure();
imagesc(img3(:,:,11)); caxis([0 50]);
epi3 = drawpolygon(gca,'Position', [coords_cell{11,1}(:,1)/3 + 2/3, coords_cell{11,1}(:,2)/3 + 2/3]);
endo3 = drawpolygon(gca,'Position', [coords_cell{11,2}(:,1)/3 + 2/3, coords_cell{11,2}(:,2)/3 + 2/3]);
epi_mask3 = createMask(epi3);
endo_mask3 = createMask(endo3);

img4 = whatsinit{1};
figure();
imagesc(img4(:,:,11)); caxis([0 50]);
epi4 = drawpolygon(gca,'Position', [coords_cell{11,1}(:,1)/4 + 3/4, coords_cell{11,1}(:,2)/4 + 3/4]);
endo4 = drawpolygon(gca,'Position', [coords_cell{11,2}(:,1)/4 + 3/4, coords_cell{11,2}(:,2)/4 + 3/4]);
epi_mask4 = createMask(epi4);
endo_mask4 = createMask(endo4);

myo_mask2 = epi_mask2 - endo_mask2;
myo_mask3 = epi_mask3 - endo_mask3;
myo_mask4 = epi_mask4 - endo_mask4;

%%
figure(); imagesc(myo_mask(:,:,6).*img(:,:,11)); caxis([0 50]);
figure(); imagesc(myo_mask2.*img2(:,:,11)); caxis([0 50]);
figure(); imagesc(myo_mask3.*img3(:,:,11)); caxis([0 50]);
figure(); imagesc(myo_mask4.*img4(:,:,11)); caxis([0 50]);
%%
figure();
imagesc(img(:,:,11)); axis image; caxis([0 50]);
roi_mi1 = drawpolygon(gca);
roi_coords = roi_mi1.Position;
mask_mi1 = createMask(roi_mi1);
figure(); imagesc(mask_mi1.*img(:,:,11)); axis image; caxis([0 50]);

figure();
imagesc(img2); axis image; caxis([0 50])
roi_mi2 = drawpolygon(gca, 'Position', [roi_coords(:,1)/2+1/2 roi_coords(:,2)/2+1/2]);
mask_mi2 = createMask(roi_mi2);
figure(); imagesc(mask_mi2.*img2); axis image; caxis([0 50]);

figure();
imagesc(img3); axis image; caxis([0 50])
roi_mi3 = drawpolygon(gca, 'Position', [roi_coords(:,1)/3+2/3 roi_coords(:,2)/3+2/3]);
mask_mi3 = createMask(roi_mi3);
figure(); imagesc(mask_mi3.*img3); axis image; caxis([0 50]);

figure();
imagesc(img4); axis image; caxis([0 50])
roi_mi4 = drawpolygon(gca, 'Position', [roi_coords(:,1)/4+3/4 roi_coords(:,2)/4+3/4]);
mask_mi4 = createMask(roi_mi4);
figure(); imagesc(mask_mi4.*img4); axis image; caxis([0 50]);

%%
figure(); imagesc(mask_mi1.*myo_mask(:,:,6).*img(:,:,11)); caxis([0 50]);
figure(); imagesc(myo_mask2.*mask_mi2.*img2); caxis([0 50]);
figure(); imagesc(myo_mask3.*img3.*mask_mi3); caxis([0 50]);
figure(); imagesc(myo_mask4.*img4.*mask_mi4); caxis([0 50]);

%%

addpath('D:\src\AHA16Segment');
Segn = 50;
Groove = 0;
[Segmentpix1, stats1, Mask_Segn1] = AHASegmentation(img(:,:,11), myo_mask(:,:,6), Segn, Groove);
[Segmentpix2, stats2, Mask_Segn2] = AHASegmentation(img2, myo_mask2, Segn, Groove);
[Segmentpix3, stats3, Mask_Segn3] = AHASegmentation(img3, myo_mask3, Segn, Groove);
[Segmentpix4, stats4, Mask_Segn4] = AHASegmentation(img4, myo_mask4, Segn, Groove);

%%
figure(); imagesc(Mask_Segn1.* mask_mi1);
uniq = unique(Mask_Segn1.* mask_mi1);
uniq = nonzeros(uniq);

MIpix1 = cell(length(uniq), 1);
MIpix2 = cell(length(uniq), 1);
MIpix3 = cell(length(uniq), 1);
MIpix4 = cell(length(uniq), 1);
Img = img(:,:,11);
for i = 1:length(uniq)
    MIpix1{i} = Img(Mask_Segn1.* mask_mi1 == uniq(i));
    MIpix2{i} = Img(Mask_Segn2.* mask_mi2 == uniq(i));
    MIpix3{i} = Img(Mask_Segn3.* mask_mi3 == uniq(i));
    MIpix4{i} = Img(Mask_Segn4.* mask_mi4 == uniq(i));
end
%Segmentpix{n,m}=Img(nMask==1);
%%
figure(); imagesc(Mask_Segn1.* myo_mask_endo);
figure(); imagesc(Mask_Segn1.* myo_mask_epi);

Segpix1 = cell(Segn, 2);
Segpix2 = cell(Segn, 2);
Segpix3 = cell(Segn, 2);
Segpix4 = cell(Segn, 2);
Img = img(:,:,11);
for i = 1:Segn
    Segpix1{i,1} = Img(Mask_Segn1.* myo_mask_endo == i);
    Segpix2{i,1} = img2(Mask_Segn2.* myo_mask_endo2 == i);
    Segpix3{i,1} = img3(Mask_Segn3.* myo_mask_endo3 == i);
    Segpix4{i,1} = img4(Mask_Segn4.* myo_mask_endo4 == i);
end

for i = 1:Segn
    Segpix1{i,2} = Img(Mask_Segn1.* myo_mask_epi == i);
    Segpix2{i,2} = img2(Mask_Segn2.* myo_mask_epi2 == i);
    Segpix3{i,2} = img3(Mask_Segn3.* myo_mask_epi3 == i);
    Segpix4{i,2} = img4(Mask_Segn4.* myo_mask_epi4 == i);
end

%%
clear res04 res08 res12 res16 perc_array1 perc_array2 perc_array3 perc_array4
perc_array1 = zeros(Segn*2, 1);
perc_array2 = zeros(Segn*2, 1);
perc_array3 = zeros(Segn*2, 1);
perc_array4 = zeros(Segn*2, 1);
for k = 1:2
    for i = 1:Segn
        perc_array1(i+(k-1)*Segn) = sum(Segpix1{i,k}<thresh1) / length(Segpix1{i,k});
        perc_array2(i+(k-1)*Segn) = sum(Segpix2{i,k}<thresh2) / length(Segpix2{i,k});
        perc_array3(i+(k-1)*Segn) = sum(Segpix3{i,k}<thresh3) / length(Segpix3{i,k});
        perc_array4(i+(k-1)*Segn) = sum(Segpix4{i,k}<thresh4) / length(Segpix4{i,k});
    end
end

res04 = perc_array1 > 0.1;
res08 = perc_array2 > 0.1;
res12 = perc_array3 > 0.1;
res16 = perc_array4 > 0.1;
%%
clear res04 res08 res12 res16 perc_array1 perc_array2 perc_array3 perc_array4
perc_array1 = zeros(length(uniq), 1);
perc_array2 = zeros(length(uniq), 1);
perc_array3 = zeros(length(uniq), 1);
perc_array4 = zeros(length(uniq), 1);
for i = 1:length(uniq)
    perc_array1(i) = sum(MIpix1{i}<thresh1) / length(MIpix1{i});
    perc_array2(i) = sum(MIpix2{i}<thresh2) / length(MIpix2{i});
    perc_array3(i) = sum(MIpix3{i}<thresh3) / length(MIpix3{i});
    perc_array4(i) = sum(MIpix4{i}<thresh4) / length(MIpix4{i});
end

res04 = perc_array1 > 0.1;
res08 = perc_array2 > 0.1;
res12 = perc_array3 > 0.1;
res16 = perc_array4 > 0.1;
%%
clear res04 res08 res12 res16 perc_array1 perc_array2 perc_array3 perc_array4
perc_array1 = zeros(Segn, 1);
perc_array2 = zeros(Segn, 1);
perc_array3 = zeros(Segn, 1);
perc_array4 = zeros(Segn, 1);

thresh1 = mean(nonzeros(img(:,:,11).*mask)) - 2*std(nonzeros(img(:,:,11).*mask));
thresh2 = mean(nonzeros(img2.*mask2)) - 2*std(nonzeros(img2.*mask2));
thresh3 = mean(nonzeros(img3.*mask3)) - 2*std(nonzeros(img3.*mask3));
thresh4 = mean(nonzeros(img4.*mask4)) - 2*std(nonzeros(img4.*mask4));
for i = 1:Segn
    perc_array1(i) = sum(Segmentpix1{i}<thresh1) / length(Segmentpix1{i});
    perc_array2(i) = sum(Segmentpix2{i}<thresh2) / length(Segmentpix2{i});
    perc_array3(i) = sum(Segmentpix3{i}<thresh3) / length(Segmentpix3{i});
    perc_array4(i) = sum(Segmentpix4{i}<thresh4) / length(Segmentpix4{i});
end

res04 = perc_array1 > 0.1;
res08 = perc_array2 > 0.1;
res12 = perc_array3 > 0.1;
res16 = perc_array4 > 0.1;
%%
[cm, order] = confusionmat(res04,res16);
figure();confusionchart(cm, order); % sens 0.40, spec 0.93
set(gca, 'FontSize', 18);
xlabel('1.6 x 1.6 mm^2'); ylabel('0.4 x 0.4 mm^2');
cm(2,2) / (cm(2,2) + cm(2,1))
cm(1,1) / (cm(1,1) + cm(1,2))

[cm, order] = confusionmat(res04,res12); % 0.80, 0.96
figure();confusionchart(cm, order);
set(gca, 'FontSize', 18);
xlabel('1.2 x 1.2 mm^2'); ylabel('0.4 x 0.4 mm^2');
cm(2,2) / (cm(2,2) + cm(2,1))
cm(1,1) / (cm(1,1) + cm(1,2))

[cm, order] = confusionmat(res04,res08); % 0.80, 1.0
figure();confusionchart(cm, order);
set(gca, 'FontSize', 18);
xlabel('0.8 x 0.8 mm^2'); ylabel('0.4 x 0.4 mm^2');
cm(2,2) / (cm(2,2) + cm(2,1))
cm(1,1) / (cm(1,1) + cm(1,2))

%% 
img = whatsinit{end};
img_size = size(img);

center = [img_size(1)/2, img_size(2)/2];
figure();
imagesc(img(:,:,11));axis image;caxis([0 50]);
hold on;
plot(center(1), center(2), '*r');
[x, y] = ginput(2);

%x = [70 180];
%y = [50 200];
k = (y(2) - y(1))/(x(2) - x(1));
angle = abs(atand(k));
%angle = 0;
w = sqrt((y(2) - y(1))^2 + (x(2) - x(1))^2);
h = 24;
c = [(x(2)+x(1))/2, (y(2)+y(1))/2]';


%% Remote %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure();
%[hdl, rot_coords] = drawRectangleonImageAtAngle(img(:,:,11),c,w,h,angle);
%roi = drawpolygon(gca,'Position',[rot_coords(1,:);rot_coords(2,:)]');
imagesc(img(:,:,11)); axis image; caxis([0 50])
%roi = drawpolygon(gca);

roi = drawpolygon(gca, 'Position', [roi.roi2.remote_coords(:,1) roi.roi2.remote_coords(:,2)]);
roi_coords = roi.Position;
mask = createMask(roi);
figure(); imagesc(mask.*img(:,:,11)); axis image; caxis([0 50]);
%% read images
img2 = whatsinit{3}(:,:,11);
img3 = whatsinit{2}(:,:,11);
img4 = whatsinit{1}(:,:,11);
%%
figure();
%[hdl2, rot_coords2] = drawRectangleonImageAtAngle(img2,c/2+1/2,w/2,h/2,angle);
imagesc(img2); axis image; caxis([0 50])
roi_coords2 = roi_coords / 2;
roi2 = drawpolygon(gca,'Position',[roi_coords2(:,1)+1/2 roi_coords2(:,2)+1/2]);
mask2 = createMask(roi2);
figure(); imagesc(mask2.*img2); axis image; caxis([0 50])

figure();
imagesc(img3); axis image; caxis([0 50])
% [hdl3, rot_coords3] = drawRectangleonImageAtAngle(img3,c/3+2/3,w/3,h/3,angle);
% roi3 = drawpolygon(gca,'Position',[rot_coords3(1,:);rot_coords3(2,:)]');
roi_coords3 = roi_coords / 3;
roi3 = drawpolygon(gca,'Position',[roi_coords3(:,1)+2/3 roi_coords3(:,2)+2/3]);
mask3 = createMask(roi3);
figure(); imagesc(mask3.*img3); axis image; caxis([0 50])


figure();
imagesc(img4); axis image; caxis([0 50])
%[hdl4, rot_coords4] = drawRectangleonImageAtAngle(img4,c/4+3/4,w/4,h/4,angle);
%roi4 = drawpolygon(gca,'Position',[rot_coords4(1,:);rot_coords4(2,:)]');
roi_coords4 = roi_coords / 4;
roi4 = drawpolygon(gca,'Position',[roi_coords4(:,1)+3/4 roi_coords4(:,2)+3/4]);
mask4 = createMask(roi4);
figure(); imagesc(mask4.*img4); axis image; caxis([0 50])

%% ----------- MI ---------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure();
imagesc(img(:,:,11)); axis image; caxis([0 50])
roi_mi = drawpolygon(gca, 'Position', [roi.roi2.roi_coords(:,1) roi.roi2.roi_coords(:,2)]);
roi_coords = roi_mi.Position;
mask_mi1 = createMask(roi_mi);
figure(); imagesc(mask_mi1.*img(:,:,11)); axis image; caxis([0 50]);

figure();
imagesc(img2); axis image; caxis([0 50])
roi_mi2 = drawpolygon(gca, 'Position', [roi_coords(:,1)/2+1/2 roi_coords(:,2)/2+1/2]);
mask_mi2 = createMask(roi_mi2);
figure(); imagesc(mask_mi2.*img2); axis image; caxis([0 50]);

figure();
imagesc(img3); axis image; caxis([0 50])
roi_mi3 = drawpolygon(gca, 'Position', [roi_coords(:,1)/3+2/3 roi_coords(:,2)/3+2/3]);
mask_mi3 = createMask(roi_mi3);
figure(); imagesc(mask_mi3.*img3); axis image; caxis([0 50]);

figure();
imagesc(img4); axis image; caxis([0 50])
roi_mi4 = drawpolygon(gca, 'Position', [roi_coords(:,1)/4+3/4 roi_coords(:,2)/4+3/4]);
mask_mi4 = createMask(roi_mi4);
figure(); imagesc(mask_mi4.*img4); axis image; caxis([0 50]);
%%
masked1 = img(:,:,11).*mask;
B = imgaussfilt(img(:,:,11));
B = round(B .* mask);
B(B>60) = 0;
masked1(masked1>60) = 0;
forhist1 = nonzeros(masked1);
masked2 = img2.*mask2;
masked2(masked2>60) = 0;
forhist2 = nonzeros(masked2);
masked3 = img3.*mask3;
masked3(masked3>60) = 0;
forhist3 = nonzeros(masked3);
masked4 = img4.*mask4;
masked4(masked4>60) = 0;
forhist4 = nonzeros(masked4);
figure();

subplot(2,2,1)
histogram(forhist1);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('0.4x0.4 mm^2')
subplot(2,2,2)
histogram(forhist2);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('0.8x0.8 mm^2')
subplot(2,2,3)
histogram(forhist3);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('1.2x1.2 mm^2')
subplot(2,2,4) 
histogram(forhist4);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('1.6x1.6 mm^2')
%%
x = [15:65];
GMModel1 = fitgmdist(forhist1, 1);
GMModel2 = fitgmdist(forhist2, 1);
GMModel3 = fitgmdist(forhist3, 1);
GMModel4 = fitgmdist(forhist4, 1);

y_mi1 = normpdf(x,GMModel1.mu(1),GMModel1.Sigma(1));
y_mi2 = normpdf(x,GMModel2.mu(1),GMModel2.Sigma(1));
y_mi3 = normpdf(x,GMModel3.mu(1),GMModel3.Sigma(1));
y_mi4 = normpdf(x,GMModel4.mu(1),GMModel4.Sigma(1));


%%
x = [15:65];
nel = unique(forhist1);
count_mi1 = zeros(1, numel(nel));
for n = 1:length(nel)
count_mi1(n) = sum(forhist1 == nel(n));
end

nel_r = unique(forhist1_remote);
count_remote1 = zeros(1, numel(nel_r));
for n = 1:length(nel_r)
count_remote1(n) = sum(forhist1_remote == nel_r(n));
end
%%
nel2 = unique(forhist2);
count_mi2 = zeros(1, numel(nel2));
for n = 1:length(nel2)
count_mi2(n) = sum(forhist2 == nel2(n));
end
%%
forhist1_B = nonzeros(B);
nel_B = unique(forhist1_B);
count_mi1_B = zeros(1, numel(nel_B));
for n = 1:length(nel_B)
count_mi1_B(n) = sum(forhist1_B == nel_B(n));
end
% B = smoothdata(count_mi1);
figure(); plot(nel, count_mi1); hold on; 
plot(nel2, count_mi2);
plot(nel_B, count_mi1_B);
yyaxis right; 
plot(x, y_mi3_1*GMModel3.ComponentProportion(1)); hold on; plot(x, y_mi3_2*GMModel3.ComponentProportion(2)); plot(x, y_mi3_3*GMModel3.ComponentProportion(3));

GMModel_B3 = fitgmdist( forhist1_B, 3,'RegularizationValue',0.01,'Options',statset('MaxIter',1500));
y_B3_1 = GMModel_B3.ComponentProportion(1) * normpdf(x, GMModel_B3.mu(1), GMModel_B3.Sigma(1));
y_B3_2 = GMModel_B3.ComponentProportion(2) * normpdf(x, GMModel_B3.mu(2), GMModel_B3.Sigma(2));
y_B3_3 = GMModel_B3.ComponentProportion(3) * normpdf(x, GMModel_B3.mu(3), GMModel_B3.Sigma(3));
figure(); plot(x, y_B3_1); hold on;
plot(x, y_B3_2); plot(x, y_B3_3);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% To interpolate image first
masked1 = img(:,:,11).*mask;
masked1(masked1>60) = 0;
forhist1 = nonzeros(masked1);

masked2 = Zq.*mask;
masked2(masked2>60) = 0;
forhist2 = nonzeros(masked2);

masked3 = Zq2.*mask;
masked3(masked3>60) = 0;
forhist3 = nonzeros(masked3);

masked4 = Zq3.*mask;
masked4(masked4>60) = 0;
forhist4 = nonzeros(masked4);

figure();

subplot(2,2,1)
histogram(forhist1);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('0.4x0.4 mm^2')
ylim([0 40])
subplot(2,2,2)
histogram(forhist2,50);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('0.8x0.8 mm^2')
ylim([0 40])
subplot(2,2,3)
histogram(forhist3,50);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('1.2x1.2 mm^2')
ylim([0 40])
subplot(2,2,4) 
histogram(forhist4, 50);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('1.6x1.6 mm^2')
ylim([0 40])
 
%%
%I = zeros(img_size(1), img_size(2));
%J = regionfill(I,rot_coords(1,:),rot_coords(2,:));
%roi = drawpolygon(gca,'Position',[rot_coords(1,:);rot_coords(2,:)]');
%mask = createMask(roi);
%imshow(mask)
%figure(); imagesc(mask); axis image;

%% Image Interpolation

Z = img2;
[X,Y] = meshgrid(1:size(img2, 1));
figure(); imagesc(Z); caxis([0 50])
dd = size(img2,1)/img_size(1);

[Xq, Yq] = meshgrid(dd:dd:size(img2,1));
Zq = interp2(X,Y,Z,Xq,Yq);
figure();imagesc(Zq); caxis([0 50]);

img3 = whatsinit{2}(:,:,11);
Z = img3;
[X,Y] = meshgrid(1:size(img3, 1));
figure(); imagesc(Z); caxis([0 50])
dd = size(img3,1)/img_size(1);

[Xq, Yq] = meshgrid(dd:dd:size(img3,1));
Zq2 = interp2(X,Y,Z,Xq,Yq);
figure();imagesc(Zq2); caxis([0 50])

img4 = whatsinit{1}(:,:,11);
Z = img4;
[X,Y] = meshgrid(1:size(img4, 1));
figure(); imagesc(Z); caxis([0 50])
dd = size(img4,1)/img_size(1);

[Xq, Yq] = meshgrid(dd:dd:size(img4,1));
Zq3 = interp2(X,Y,Z,Xq,Yq);
figure();imagesc(Zq3); caxis([0 50])
%%
h1 = figure(); plot(sort(nonzeros(mask.*img(:,:,11)))); grid on;
hold on;
plot(sort(nonzeros(mask.*Zq)));
plot(sort(nonzeros(mask.*Zq2)));
plot(sort(nonzeros(mask.*Zq3)));
legend({'0.4x0.4 mm^2', '0.8x0.8 mm^2', '1.2x1.2 mm^2', '1.6x1.6 mm^2'})
ylim([0 50])

%%
%inter = 50 - 70*k;
%linx = x(1):1:x(2);
%liny = round(k*linx + inter);

%I = zeros(256,256);
%for i = 1:length(linx)
%    I(linx(i),liny(i)) = 1;
%end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Water/formalin

figure();
imagesc(img(:,:,11)); axis image; caxis([0 50])
roi = drawpolygon(gca);

water_coords = roi.Position;
water_mask = createMask(roi);
figure(); imagesc(water_mask.*img(:,:,11)); axis image; caxis([0 50]);

%%
figure();
%[hdl2, rot_coords2] = drawRectangleonImageAtAngle(img2,c/2+1/2,w/2,h/2,angle);
imagesc(img2); axis image; caxis([0 50])
water_coords2 = water_coords / 2;
roi2 = drawpolygon(gca,'Position',[water_coords2(:,1)+1/2 water_coords2(:,2)+1/2]);
water_mask2 = createMask(roi2);
figure(); imagesc(water_mask2.*img2); axis image; caxis([0 50])

figure();
imagesc(img3); axis image; caxis([0 50])
% [hdl3, rot_coords3] = drawRectangleonImageAtAngle(img3,c/3+2/3,w/3,h/3,angle);
% roi3 = drawpolygon(gca,'Position',[rot_coords3(1,:);rot_coords3(2,:)]');
water_coords3 = water_coords / 3;
roi3 = drawpolygon(gca,'Position',[water_coords3(:,1)+2/3 water_coords3(:,2)+2/3]);
water_mask3 = createMask(roi3);
figure(); imagesc(water_mask3.*img3); axis image; caxis([0 50])


figure();
imagesc(img4); axis image; caxis([0 50])
%[hdl4, rot_coords4] = drawRectangleonImageAtAngle(img4,c/4+3/4,w/4,h/4,angle);
%roi4 = drawpolygon(gca,'Position',[rot_coords4(1,:);rot_coords4(2,:)]');
water_coords4 = water_coords / 4;
roi4 = drawpolygon(gca,'Position',[water_coords4(:,1)+3/4 water_coords4(:,2)+3/4]);
water_mask4 = createMask(roi4);
figure(); imagesc(water_mask4.*img4); axis image; caxis([0 50])

%% 
masked1 = img(:,:,11).*water_mask;
%masked1(masked1>60) = 0;
forhist1 = nonzeros(masked1);
masked2 = img2.*water_mask2;
%masked2(masked2>60) = 0;
forhist2 = nonzeros(masked2);
masked3 = img3.*water_mask3;
%masked3(masked3>60) = 0;
forhist3 = nonzeros(masked3);
masked4 = img4.*water_mask4;
%masked4(masked4>60) = 0;
forhist4 = nonzeros(masked4);
figure();

subplot(2,2,1)
histogram(forhist1,10);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('0.4x0.4 mm^2')
subplot(2,2,2)
histogram(forhist2,10);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('0.8x0.8 mm^2')
subplot(2,2,3)
histogram(forhist3,10);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('1.2x1.2 mm^2')
subplot(2,2,4) 
histogram(forhist4,10);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('1.6x1.6 mm^2')

%% To interpolate image first
masked1 = img(:,:,11).*water_mask;
%masked1(masked1>60) = 0;
forhist1 = nonzeros(masked1);

masked2 = Zq.*water_mask;
%masked2(masked2>60) = 0;
forhist2 = nonzeros(masked2);

masked3 = Zq2.*water_mask;
%masked3(masked3>60) = 0;
forhist3 = nonzeros(masked3);

masked4 = Zq3.*water_mask;
%masked4(masked4>60) = 0;
forhist4 = nonzeros(masked4);

figure();

subplot(2,2,1)
histogram(forhist1,50);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('0.4x0.4 mm^2')
ylim([0 40])
subplot(2,2,2)
histogram(forhist2,50);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('0.8x0.8 mm^2')
ylim([0 40])
subplot(2,2,3)
histogram(forhist3,50);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('1.2x1.2 mm^2')
ylim([0 40])
subplot(2,2,4) 
histogram(forhist4, 50);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('1.6x1.6 mm^2')
ylim([0 40])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remote myocardium

figure();
imagesc(img(:,:,11)); axis image; caxis([0 50])
roi = drawpolygon(gca);

remote_coords = roi.Position;
remote_mask = createMask(roi);
figure(); imagesc(remote_mask.*img(:,:,11)); axis image; caxis([0 50]);


figure();
%[hdl2, rot_coords2] = drawRectangleonImageAtAngle(img2,c/2+1/2,w/2,h/2,angle);
imagesc(img2); axis image; caxis([0 50])
remote_coords2 = remote_coords / 2;
roi2 = drawpolygon(gca,'Position',[remote_coords2(:,1)+1/2 remote_coords2(:,2)+1/2]);
remote_mask2 = createMask(roi2);
figure(); imagesc(remote_mask2.*img2); axis image; caxis([0 50])

figure();
imagesc(img3); axis image; caxis([0 50])
% [hdl3, rot_coords3] = drawRectangleonImageAtAngle(img3,c/3+2/3,w/3,h/3,angle);
% roi3 = drawpolygon(gca,'Position',[rot_coords3(1,:);rot_coords3(2,:)]');
remote_coords3 = remote_coords / 3;
roi3 = drawpolygon(gca,'Position',[remote_coords3(:,1)+2/3 remote_coords3(:,2)+2/3]);
remote_mask3 = createMask(roi3);
figure(); imagesc(remote_mask3.*img3); axis image; caxis([0 50])


figure();
imagesc(img4); axis image; caxis([0 50])
%[hdl4, rot_coords4] = drawRectangleonImageAtAngle(img4,c/4+3/4,w/4,h/4,angle);
%roi4 = drawpolygon(gca,'Position',[rot_coords4(1,:);rot_coords4(2,:)]');
remote_coords4 = remote_coords / 4;
roi4 = drawpolygon(gca,'Position',[remote_coords4(:,1)+3/4 remote_coords4(:,2)+3/4]);
remote_mask4 = createMask(roi4);
figure(); imagesc(remote_mask4.*img4); axis image; caxis([0 50])

%% 
masked1 = img(:,:,11).*remote_mask;
%masked1(masked1>60) = 0;
forhist1_remote = nonzeros(masked1);
masked2 = img2.*remote_mask2;
%masked2(masked2>60) = 0;
forhist2_remote = nonzeros(masked2);
masked3 = img3.*remote_mask3;
%masked3(masked3>60) = 0;
forhist3_remote = nonzeros(masked3);
masked4 = img4.*remote_mask4;
%masked4(masked4>60)= 0;
forhist4_remote = nonzeros(masked4);

GMModel_remote1 = fitgmdist(forhist1_remote, 1);
GMModel_remote2 = fitgmdist(forhist2_remote, 1);
GMModel_remote3 = fitgmdist(forhist3_remote, 1);
GMModel_remote4 = fitgmdist(forhist4_remote, 1);
x_remote = [20:40];
y_remote1 = normpdf(x_remote,GMModel_remote1.mu(1),GMModel_remote1.Sigma(1));
y_remote2 = normpdf(x_remote,GMModel_remote2.mu(1),GMModel_remote2.Sigma(1));
y_remote3 = normpdf(x_remote,GMModel_remote3.mu(1),GMModel_remote3.Sigma(1));
y_remote4 = normpdf(x_remote,GMModel_remote4.mu(1),GMModel_remote4.Sigma(1));

figure();
subplot(2,2,1)
histogram(forhist1_remote,10);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('0.4x0.4 mm^2')
xlim([20 40]);
hold on; yyaxis right; plot(x_remote, y_remote1); 

subplot(2,2,2)
histogram(forhist2_remote,10);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('0.8x0.8 mm^2');
xlim([20 40]);
hold on; yyaxis right; plot(x_remote, y_remote2); 

subplot(2,2,3)
histogram(forhist3_remote,10);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('1.2x1.2 mm^2');
xlim([20 40]);
hold on; yyaxis right; plot(x_remote, y_remote3); 

subplot(2,2,4) 
histogram(forhist4_remote,10);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('1.6x1.6 mm^2');
xlim([20 40]);
hold on; yyaxis right; plot(x_remote, y_remote4); 

%% 
nel = unique(forhist1_remote);
count_remote = zeros(1, numel(nel));
for n = 1:length(nel)
count_remote(n) = sum(forhist1_remote == nel(n));
end
figure(); plot(nel, count_remote); hold on; yyaxis right; plot(x_remote, y_remote1); 
%% To interpolate image first
masked1 = img(:,:,11).*remote_mask;
%masked1(masked1>60) = 0;
forhist1 = nonzeros(masked1);

masked2 = Zq.*remote_mask;
%masked2(masked2>60) = 0;
forhist2 = nonzeros(masked2);

masked3 = Zq2.*remote_mask;
%masked3(masked3>60) = 0;
forhist3 = nonzeros(masked3);

masked4 = Zq3.*remote_mask;
%masked4(masked4>60) = 0;
forhist4 = nonzeros(masked4);

figure();
subplot(2,2,1)
histogram(forhist1,10);xlabel('T2* (ms)'); ylabel('Frequency');
%histfit(forhist1);
%x = [20:40];

set(gca, 'FontSize', 16); title('0.4x0.4 mm^2')
ylim([0 40]); xlim([20 40]);
subplot(2,2,2)
histogram(forhist2,10);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('0.8x0.8 mm^2')
ylim([0 40]); xlim([20 40]);
subplot(2,2,3)
histogram(forhist3,10);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('1.2x1.2 mm^2')
ylim([0 40]); xlim([20 40]);
subplot(2,2,4)
histogram(forhist4, 10);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('1.6x1.6 mm^2')
ylim([0 40]); xlim([20 40]);

%pd = fitdist(forhist1,'Normal');

%% Overlay with guassian curves



x = [10:70];
figure();
histogram(forhist1);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('0.4x0.4 mm^2')
yyaxis right; 
plot(x, y_mi1*GMModel1.ComponentProportion(1)); 


figure();
histogram(forhist1);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('0.4x0.4 mm^2')
yyaxis right; 
plot(x, y_mi2_1*GMModel2.ComponentProportion(1)); hold on; plot(x, y_mi2_2*GMModel2.ComponentProportion(2));

figure();
histogram(forhist1);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('0.4x0.4 mm^2')
yyaxis right; 
plot(x, y_mi3_1*GMModel3.ComponentProportion(1)); hold on; plot(x, y_mi3_2*GMModel3.ComponentProportion(2)); plot(x, y_mi3_3*GMModel3.ComponentProportion(3));

figure();
histogram(forhist1);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('0.4x0.4 mm^2')
yyaxis right; 
plot(x, y_mi4_1*GMModel4.ComponentProportion(1)); hold on; plot(x, y_mi4_2*GMModel4.ComponentProportion(2)); plot(x, y_mi4_3*GMModel4.ComponentProportion(3));plot(x, y_mi4_4*GMModel4.ComponentProportion(4));

figure();
histogram(forhist1);xlabel('T2* (ms)'); ylabel('Frequency');
set(gca, 'FontSize', 16); title('0.4x0.4 mm^2')
yyaxis right; 
plot(x, y_mi5_1*GMModel5.ComponentProportion(1)); hold on; plot(x, y_mi5_2*GMModel5.ComponentProportion(2)); plot(x, y_mi5_3*GMModel5.ComponentProportion(3));plot(x, y_mi5_4*GMModel5.ComponentProportion(4));plot(x, y_mi5_5*GMModel5.ComponentProportion(5));



%% Save coordinates

roi2 = struct;
roi2.remote_coords = remote_coords;
roi2.water_coords = water_coords;
roi2.roi_coords = roi_coords;
f_dst = GetFullPath(cat(2, base_dir, '\..\roi2.mat'));
save(f_dst,'roi2');
