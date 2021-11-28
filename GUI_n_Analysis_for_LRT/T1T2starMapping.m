clear all; 
close all;
%% 
[fid_file, fid_path] = uigetfile('*.mat');
load(strcat(fid_path, fid_file), 'dispim', 'Gr', 'Phi', 'L', 'U', 'Ny', 'Nx', 'Nz', 'vec','params', 'sizes');
% dispim = @(x)fftshift(x(:,:,SliceSlider.Value,:),1);
qMRinfo('mono_t2'); % set it up first in qMRLab
% T2star
% sizes(2) -> T1
% sizes(3) -> Cardiac
% sizes(4) -> resp
TE_array = [1.41, 3.38, 5.39, 7.40, 9.41, 11.42]; % 13.43, 15.44 ms
IR_array = repmat(TE_array, [sizes(2), 1]) + repmat((0:1:(sizes(2)-1)) * params.lEchoSpacing, [sizes(5), 1]).'*1000;

mask = zeros(Ny, Nx, Nz);
for i = 1:Nz
    dispim = @(x,st)fftshift(x(:,:,i,:),1);
    for j = 1:1
        for k = 1:1
            temp = Gr\reshape(Phi(:,:,j,k,:), L, []);
            temp = reshape(reshape(dispim(reshape(U, Ny, Nx, Nz, [])),[],L) * temp, Ny, Nx, [], params.NEco);
            cw = 0.5*max(vec(abs(temp)));
            figure();
            imagesc(abs(temp(:,:,end,1))/cw); axis image;
            roi = drawpolygon;
            mask(:,:,i) = createMask(roi);
        end
    end
end

%%
t1_map = zeros(Ny, Nx, Nz, 1, 1, 1);
t2star_map = zeros(Ny, Nx, Nz, 1, 1, sizes(2));
for i = 1:Nz
    dispim = @(x,st)fftshift(x(:,:,i,:),1);
    temp = Gr\reshape(Phi(:,:,j,k,:), L, []);
    temp = reshape(reshape(dispim(reshape(U, Ny, Nx, Nz, [])),[],L) * temp, Ny, Nx, [], params.NEco);
    for j = 1:1
        for k = 1:1
            temp = Gr\reshape(Phi(:,:,j,k,:), L, []);
            temp = reshape(reshape(dispim(reshape(U, Ny, Nx, Nz, [])),[],L) * temp, Ny, Nx, [], params.NEco);
            for l = 1:sizes(2)
                % Reshape matrix as [Width x Height x #Slice x #TE]
                ipt = abs(temp(:,:,l,:));
                Model = mono_t2;  % Create class from model
                %Model = Custom_OptionsGUI(Model);
                Model.Prot.SEdata.Mat = TE_array.'; %
                Model.st = [100 2000];
                Model.lb = [1 2000];
                Model.fx = [0 0];
                Model.voxelwise = 1;
                Model.options.FitType = 'Linear';
                data = struct;  % Create data structure
                data.SEdata = ipt;
                data.Mask = mask(:,:,i);
                FitResults = FitData(data, Model); %fit data
                t2star_map(:,:,i,j,k,l) = FitResults.T2;
            end
            
            for m = 1:1
                % a - create object
                Model = inversion_recovery;
                % Reshape matrix as [Width x Height x #Slice x #TE]
                ipt = abs(reshape(temp(:,:,:,m), Ny, Nx, 1, []));
                
                data = struct;
                data.IRData= double(ipt);
                data.Mask= double(mask(:,:,i));

                Model.Prot.IRData.Mat = IR_array(:,m);
                Model.voxelwise = 1;

                % b- fit dataset
                FitResults = FitData(data,Model,0);
                t1_map(:,:,i,j,k,m) = FitResults.T1 .* (-FitResults.rb ./ FitResults.ra - 1);
            end
        end
    end
end

map_to_save = struct;
map_to_save.mask = mask;
map_to_save.t2star_map = t2star_map;
map_to_save.t1_map = t1_map;
save_f = cat(2, fid_path, fid_file(1:15), 'LRT_Mappings.mat');
save(save_f, 'map_to_save');

%% AHA analysis
