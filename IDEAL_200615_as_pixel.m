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

%% initialize psi0

Import_data = temp(Y_lower:Y_upper,X_lower:X_upper,1:8)/cw;
TE = xvector;

M = 2;% water + 1 fat peak
delta_frequencies = [0, 420]; %Hz
sz = size(Import_data);
resolution = [size_of_import_data(1),size_of_import_data(2)];
echoes = size_of_import_data(3);

% psi0 = repmat((42.58e6)*2.89,resolution);%Hz B0 = 2.89T.
% psi0t = repmat((42.58e6)*2.89, [resolution echoes]);

psi0 = repmat(double(0),resolution);%Hz B0 = 2.89T.
psi0t = repmat(double(0), [resolution echoes]);

for i = 1:echoes
    psi0t(:,:,i) = psi0.*TE(i);
end

% signal_hat = Import_data.*exp(-1i*2*pi*psi0t);
%% define matrix A

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
water_signal = [];
fat_signal = [];
PDFF = [];

rho_hat = repmat(double(0),[resolution 2*M]);
delta_rho_hat = repmat([0],[resolution 2*M+1]);

%% Calculate psi

for i = 1 : sz(1)
    for j = 1 : sz(2)
        for iter_num = 1:30
            
            signal_hat(i,j,:) = Import_data(i,j,:).*exp(-1i*2*pi*psi0t(i,j,:));
            s_hat = [squeeze(real(signal_hat(i,j,:)));squeeze(imag(signal_hat(i,j,:)))];
            rho_hat(i,j,:) = LS_A * s_hat; 

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
                    -(rho_hat(i,j,3)*A(2*M-1,m+echoes))...signal_hat(i,j,:) = Import_data(i,j,:).*exp(-1i*2*pi*psi0t(i,j,:));
                    -(rho_hat(i,j,4)*A(2*M,m+echoes));
            end
            
            matrix_B = [[g_r, g_i]',matrix_A];
            LS_B = (matrix_B'*matrix_B)\matrix_B';
            delta_rho_hat(i,j,:) = LS_B * s_hat_hat';
            psi0(i,j) = psi0(i,j) + delta_rho_hat(i,j,1);
            
            for m = 1:echoes
                psi0t(i,j,m) = psi0(i,j).*TE(m);
            end
           
            if abs(delta_rho_hat(i,j,1))<0.05
                break;
            end
        end
        
        %smooth psi
        
        signal_hat(i,j,:) = Import_data(i,j,:).*exp(-1i*2*pi*psi0t(i,j,:));
        s_hat = [squeeze(real(signal_hat(i,j,:)));squeeze(imag(signal_hat(i,j,:)))];
        rho_hat(i,j,:) = LS_A * s_hat; 
        water_signal(i,j) = rho_hat(i,j,1)+rho_hat(i,j,2)*1i;
        fat_signal(i,j) = rho_hat(i,j,3)+rho_hat(i,j,4)*1i;
        PDFF(i,j) = abs(fat_signal(i,j))/(abs(fat_signal(i,j))+abs(water_signal(i,j)));
    end
end
%%

imshow(PDFF)










