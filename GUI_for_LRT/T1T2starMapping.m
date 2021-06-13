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
t1_map = zeros(Ny, Nx, Nz, sizes(3), sizes(4), sizes(5));
t2star_map = zeros(Ny, Nx, Nz, 1, 1, sizes(2));
for i = 1:Nz
    dispim = @(x,st)fftshift(x(:,:,i,:),1);
    for j = 1:1
        for k = sizes(4):sizes(4)
            temp = Gr\reshape(Phi(:,:,j,k,:), L, []);
            temp = reshape(reshape(dispim(reshape(U, Ny, Nx, Nz, [])),[],L) * temp, Ny, Nx, [], params.NEco);
            cw = 0.5*max(vec(abs(temp)));
            figure();
            imagesc(abs(temp(:,:,end,1))/cw); axis image;
            roi = drawpolygon;
            mask = createMask(roi);
            
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
                data.Mask = mask;
                FitResults = FitData(data, Model); %fit data
                
                t2star_map(:,:,i,1,1,l) = FitResults.T2;
            end
        end
    end
end



ax1 = implay(abs(temp(:,:,:,1)/cw));