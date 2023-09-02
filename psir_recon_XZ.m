clear all;
close all;

%% Load Data
[fid_file, fid_path] = uigetfile('*.mat');

load(cat(2, fid_path, fid_file));

%% Display Image
dispim = @(x,st)fftshift(x(:,:,3,:),1);

Phi=reshape(Phi,[L sizes(2:end)]);
temp=Gr\reshape(Phi(:,:,1,end,end),L,[]);

temp=(reshape(reshape(dispim(reshape(U,Ny,Nx,Nz,[])),[],L)*temp,Ny,Nx,[],params.NEco));
cw=0.5*max(vec(abs(temp)));

h = implay(abs(temp(:,:,:))/cw);
set(h.Parent,'Name','old_basal_echo1');

%% 
phase_temp = angle(temp(:,:,:));

figure(); montage(abs(temp)/cw); axis image;
s_array = [6, 11, 16, 21, 41, 61, 81, 101, 121, 141, 161, 181];

figure(); 
for i = 1:length(s_array)
    subplot(3,4,i);
    imagesc(abs(temp(:,:,s_array(i)))); axis image; colormap gray;
end

figure(); 
for i = 1:length(s_array)
    subplot(3,4,i);
    imagesc(abs(phase_temp(:,:,s_array(i)))); axis image; colormap gray;
end
%%
phase_diff = phase_temp - phase_temp(:,:,end);
cos_phase_diff = cos(phase_diff);
cos_phase_diff(cos_phase_diff>=0) = 1;
cos_phase_diff(cos_phase_diff<0) = -1;

figure();
for i = 1:length(s_array)
    subplot(3,4,i);
    imagesc(cos_phase_diff(:,:,s_array(i))); axis image; colormap gray;
    %caxis([-pi pi]);
end
%%
%phase_diff(phase_diff >= 0) = 1;
%phase_diff(phase_diff < 0) = -1;
%figure();
%for i = 1:length(s_array)
%    subplot(3,4,i);
%    imagesc(phase_diff(:,:,s_array(i))); axis image; colormap gray;
%end
%%
temp_psir = abs(temp) .* cos(phase_diff);
%temp_psir = abs(temp) .* cos_phase_diff;

psir = (temp_psir - min(temp_psir(:))) ./ (max(temp_psir(:)) - min(temp_psir(:)));
figure(); 
for i = 1:length(s_array)
   subplot(3,4,i);
   imagesc(psir(:,:,s_array(i))); axis image; colormap gray;
   caxis([0.2 0.8]);
end

figure(); plot(squeeze(psir(110,100,:))); hold on;
plot(squeeze(psir(104,102,:)));

%%
phase_temp2 = angle(temp./temp(:,:,end));
figure(); 
for i = 1:length(s_array)
    subplot(3,4,i);
    imagesc(phase_temp2(:,:,s_array(i))); axis image; colormap gray;
end

figure(); 
for i = 1:length(s_array)
    subplot(3,4,i);
    imagesc(phase_temp2(:,:,s_array(i)) - phase_temp2(:,:,181)); axis image; colormap gray;
end