%% Import data

% clc
% clear
% close all
% path_datfile='D:\MATLAB\IDEAL\Import_data';
% cd(path_datfile);
% path_temp = 'D:\MATLAB\codes';
% addpath(path_datfile,path_temp);
% 
% Files = dir('*.ima');
%A = dicomread(Files(1).name);

% for i = 1:length(Files)
%     Import_data(:,:,i) = double(dicomread(Files(i).name));
%     info = dicominfo(Files(i).name);
%     TE(i) = info.EchoTime;
% end
% 
% disp('Reading data completed. Data size:');
% disp(size(Import_data));

%% estimate psi0

Import_data = temp(Y_lower:Y_upper,X_lower:X_upper,1:8)/cw;
TE = xvector;

M = 2;% water + 1 fat peak
delta_frequencies = [0, 420]; %Hz
size_of_import_data = size(Import_data);
resolution = [size_of_import_data(1),size_of_import_data(2)];
echoes = size_of_import_data(3);

% psi0 = repmat((42.58e6)*2.89,resolution);%Hz B0 = 2.89T.
% psi0t = repmat((42.58e6)*2.89, [resolution echoes]);

psi0 = repmat(double(0),resolution);%Hz B0 = 2.89T.
psi0t = repmat(double(0), [resolution echoes]);

for i = 1:echoes
    psi0t(:,:,i) = psi0.*TE(i);
end

signal_hat = Import_data.*exp(-1i*2*pi*psi0t);
%%
% matrix A

A = repmat(double(0),[2*M 2*echoes]);
for a = 1:echoes
    for b = 1:M
        A(2*b-1,a) = cos(2*pi*delta_frequencies(b)*TE(a));
        A(2*b,a) = -sin(2*pi*delta_frequencies(b)*TE(a));
        A(2*b-1,a+echoes) = -A(2*b,a);
        A(2*b,a+echoes) = A(2*b-1,a);
    end
end

matrix_A = A';
LS_A = (matrix_A' * matrix_A)\matrix_A';

g_r = [];
g_i = [];
s_hat_hat =[];

rho_hat = repmat(double(0),[resolution 2*M]);
delta_rho_hat = repmat([0],[resolution 2*M+1]);

%%
for iter_num = 1:30
    
    %initial estimation of each chemical species:
    
    for i = 1:size_of_import_data(1)
        for j = 1:size_of_import_data(2)
            s_hat = [squeeze(real(signal_hat(i,j,:)));squeeze(imag(signal_hat(i,j,:)))];
            rho_hat(i,j,:) = LS_A * s_hat; 
        end
    end
    
    disp('Estimation of each chemical species finished.');

    % calculate the the error to the field map.

    
    for i = 1:size_of_import_data(1)
        for j = 1:size_of_import_data(2)
            for m = 1:echoes
                g_r(m) = 2*pi*TE(m)*(-rho_hat(i,j,2)...
                    -(rho_hat(i,j,3)*A(2*M-1,m+echoes))...
                    -(rho_hat(i,j,4)*A(2*M,m+echoes)));
                g_i(m) = 2*pi*TE(m)*(rho_hat(i,j,1)...
                    +(rho_hat(i,j,3)*A(2*M,m+echoes))...
                    -(rho_hat(i,j,4)*A(2*M-1,m+echoes)));

                s_hat_hat(m) = real(signal_hat(i,j,m))...
                    -rho_hat(i,j,1)...
                    -(rho_hat(i,j,3)*A(2*M,m+echoes))...
                    +(rho_hat(i,j,4)*A(2*M-1,m+echoes));
                s_hat_hat(m+echoes) = imag(signal_hat(i,j,m))...
                    -rho_hat(i,j,2)...
                    -(rho_hat(i,j,3)*A(2*M-1,m+echoes))...
                    -(rho_hat(i,j,4)*A(2*M,m+echoes));
            end
            matrix_B = [[g_r, g_i]',matrix_A];

            LS_B = (matrix_B'*matrix_B)\matrix_B';

            delta_rho_hat(i,j,:) = LS_B * s_hat_hat';
        end
    end
%     

    disp('Field map error calculation finished.');

    delta_psi0 = delta_rho_hat(:,:,1);
    
%     if abs(delta_psi0)<40
%         disp('Iteration ends at iter_num=');
%         disp(iter_num);
%         break;
%     end
    
    psi0 = psi0 + delta_psi0;
    
    for i = 1:echoes
        psi0t(:,:,i) = psi0.*TE(i);
    end

%signalhat = repmat(double(0),[resolution coil_num echonum]);

%     signalhat = squeeze(images(:,:,:,slicenumber,:)).*exp(-1i*2*pi*psi0t);

    signal_hat = Import_data.*exp(-1i*2*pi*psi0t);   
    
%     signalhat(:,:,:,1) = images01(:,:,:,slicenumber).*exp(-1i*2*pi*psi0*t(1));
%     signalhat(:,:,:,2) = images02(:,:,:,slicenumber).*exp(-1i*2*pi*psi0*t(2));
%     signalhat(:,:,:,3) = images03(:,:,:,slicenumber).*exp(-1i*2*pi*psi0*t(3));
    
    disp('iter_num');
    disp(iter_num);
%     disp('finished');
%     disp('max_delta_psi');
%     disp(max(max(abs(delta_psi0))));
%     disp('mean_delta_psi')
%     disp(mean(mean(abs(delta_psi0))));
    
    res(iter_num) = max(max(abs(delta_psi0)));
%     iter_num = iter_num +1;
    
%     if (rem(iter_num, 10) == 0)
%         figurename = sprintf('psi0 coil 1 iter num %03d',iter_num);
%         figure, imagesc(sum(psi0(:,:,1),3))
%         axis equal
%         colormap gray
%         title(figurename)
%     end
end

plot(res);
















