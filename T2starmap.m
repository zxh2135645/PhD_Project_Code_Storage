% T2* fitting

clear all;

%% import dciom

images = [];
files = dir('*.dcm');
for i = 1:length(files)
    images(:,:,i) = double(dicomread(files(i).name));
    info = dicominfo(files(i).name);
    xvector(i) = info.EchoTime;
end

% temp = images(45:150,25:135,:);
% clear images;
% images = temp;

%%
% Y_lower = 39;
% Y_upper = 100;
% X_lower = 18;
% X_upper = 75;

images = abs(temp/cw);
% images = abs(temp(Y_lower:Y_upper,X_lower:X_upper,1:8)/cw);
% sz = size(images);

%twix_obj.hdr.MeasYaps.alTE 



% % xvector = [1.12 2.88 4.64 6.4 8.16 9.92 11.68 13.44];
% xvector = [1.37 3.38 5.39 7.4 9.41 11.42 13.43 15.44];
% xvector = [1.42 3.27 5.13 6.99 8.85 10.71 12.57 14.43];
xvector = [1.42 3.53 5.64 7.75 9.86 11.97 14.08 16.19];
% clear all
% directory = dir('img*.*');
% for i = 1:length(directory)
%     filename = directory(i).name;
%     images(:,:,i) = dicomread(filename);
%     info = dicominfo(filename);
%     xvector(i) = info.EchoTime;
% end

%%
sz = size(images);

ffunction_mono = [];
para_mono = [];
t2star_map_mono = [];
rsquare_map = [];
adjrsquare_map = [];



for j = 1:sz(1)
    parfor i = 1:sz(2)
        yvector = double(squeeze(images(j,i,:)));
%         if max(yvector<10)
%             t2star_map_mono(j,i) = 0;
%             rsquare_map(j,i) = 0;
%             adjrsquare_map(j,i) = 0;
%         else
%             [ffunction_mono, para_mono] = fit(xvector,yvector,'exp1','lower',[0 -Inf], 'upper',[inf -0.005]);
            [ffunction_mono, para_mono] = fit(xvector',yvector,'exp1');
            t2star_map_mono(j,i) = -1/(ffunction_mono.b);
            rsquare_map(j,i) = para_mono.rsquare;
            adjrsquare_map(j,i) = para_mono.adjrsquare;
%         end
    end
    if mod(j,10) == 0;
        j
    end
end
0


figure,
subplot(2,2,1), imshow(images(:,:,1));
subplot(2,2,2), imshow(t2star_map_mono,'DisplayRange',[0 250]);
subplot(2,2,3),imshow(rsquare_map);
subplot(2,2,4),imshow(adjrsquare_map);

% dicomwrite(uint16(adjrsquare_map.*100),'rsquare_map');
% dicomwrite(uint16(t2star_map_mono),'t2star_map');


%%

t2star_mask = [];
t2star_indices = [];


figure,imshow(t2star_map_mono,'DisplayRange',[0 250]);
roi = impoly;
roi_mask = createMask(roi);
t2star_mask = t2star_map_mono .*roi_mask;
t2star_mask = t2star_mask(:);


% t2star_indices = find(t2star_mask);
t2star_indices = t2star_mask>0;

t2star = mean(t2star_mask(t2star_indices))

sd = std(t2star_mask(t2star_indices))

%%
save('T2star');





