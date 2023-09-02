%% Set up paths
mainpath = uigetdir; %Put your path here.

addpath('../function/');
addpath(genpath(fullfile(mainpath,'supporting','recon'))); %Supporting reconstruction scripts
addpath(genpath(fullfile(mainpath,'supporting','utils'))); %Supporting utilities

%Third-party paths are below:

% NYU GPU NUFFT code. Not included. You will have to download, compile, and manage your own copy
% addpath(genpath('~/gpuNUFFT'));

% Siemens import code. Included. If you already have it, you can remove this folder from your branch and update this path to your existing folder
addpath(genpath((fullfile(mainpath,'supporting','mapVBVD'))));

% Michigan irt code. Included.  If you already have it, you can remove this folder from your branch and update this path to your existing folder
% irtdir=fullfile(mainpath,'supporting','irt');
% current_dir=pwd;
% cd(irtdir);
% irtdir=pwd; %reset irtdir to use absolute pathname
% setup;
% cd(current_dir);
addpath(genpath(fullfile(mainpath,'supporting_RY'))); %Supporting utilities
%% Setting up parameters
rbins=4; %Set # of respiratory bins here
cbins=24; %Set # of cardiac bins here

resp = 1;
card = 1;

[fid_file, fid_path] = uigetfile('*.dat');

Sub_ID = input('What is subject ID:  ', 's');
%% Loading and preprocessing
useGPU=(exist('gpuNUFFT', 'file')>1); %true if using gpuNUFFT. false if using irt.

total_time = inf; %amount of data to use, in seconds. Use "inf" to use all.
%rep = 1; %only for repeatability studies
 
%load_data;
num_eco_total = 6;
L_input = 64;

for num_eco = 1:num_eco_total
% for num_eco = 1:1  
    load_data_SingleEcho;
    disp('Data loaded');
    
    setup_trajectories;
    disp('Finished setup trajectories');
    %
    setup_functions;
    disp('Finished setup functions');
    %
    %pre_whiten_asym;
    pre_whiten;
    disp('Finished setup pre whiten');
    %
    estimate_sensitivities;
    disp('Finished SEs.')
    
    % initial reconstruction
    realtime_subspace; %default is probably overconstrained, but enough to identify respiration

    phase_drift_correction; %For long scans, such as 3D, correct for phase drift and re-calculate real-time subspace
    realtime_subspace;

    if iscartesian
        Phi_rt(:,nav_indices)=0; %set nav lines to 0
        
        %     st.w(1)=st.w(1)-size(nav_data,1);
        %     st.winv = 1./st.w;
        %     st.winv(st.w==0)=0;
        
        M=@(x)vec(ifft(ifft(bsxfun(@times,fft(fft(prep(x),[],1),[],3),reshape(1./(st.w+1),Ny,1,Nz)),[],1),[],3));
    end

ls_recon;

%Store real-time reconstruction
L_init = L;
Phi_rt_init = Phi_rt;
Phi_rt_full_init =  Phi_rt_full;
Phi_rt_small_init = Phi_rt_small;
U_init = U;

realtime_display;

if num_eco == 1
    % generate relaxation subspace
    gen_bloch_subspace;
    
    disp(cat(2, 'Number of Averages: ', params.Averages));
    seg = input(sprintf('# of Sections [%d]: ', 3));
    do_binning_ica_SingleEcho;
end

%
% tensor subspace estimation
set_tensor_params; %change smoothness/low-rankness parameters in here
tensor_subspace_SingleEcho;

%VERY approximate preview
%U=vec(reshape(U_init,[],L_init)*(Phi_rt_full_init*pinv(Phi_rt_full)));


%frac_array = [0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1];
frac_array = [0.05, 0.5];

for f = 1:length(frac_array)
%for f = 3:3

frac = frac_array(f);
U=vec(reshape(U_init,[],L_init)*(bsxfun(@times,Phi_rt_init,lrw)*pinv(Phi_rt)));
tensor_display;

disp('Finished section tensor subspace estimation.')
%u recon
ls_recon; %set least-squares iterations in here
tensor_display;

% Whiten
Wt=reshape(U0,[],L);
Wt=Wt'*Wt/size(Wt,1);
Wt = inv(sqrtm(Wt));
WtWh = Wt*Wt';

if false  % for brain DSC perfusion, should use this part
    Wt=mean(diag(Wt))*eye(L);
    WtWh=Wt*Wt';
end

% tv or wavelet
% wavelet_recon_aniso3D;
% tensor_display;


tv_recon_aniso3D;
tensor_display;

TPR = params.dThickness_mm/params.NTruePar;
USR = size(kspace_data, 1)/SGblock/params.NEco/Ny/Nz/rbins/cbins*100
% specialNotes = input('Special Notes?:  ', 's');
specialNotes = cat(2, 'Echo', num2str(num_eco), '_Seg', num2str(seg), '_lambda', num2str(frac));

if ~isempty(specialNotes)
    specialNotes = cat(2, '_', specialNotes);
end
save_results = [fid_file(15:22), '_', Sub_ID, '_', num2str(round(TPR)), 'mm_USR', num2str(round(USR)), ...
    '%_L', num2str(L), '_', 'results_', datestr(now, 'yyyy_mm_dd_HH_MM'), specialNotes, '.mat']
%clear twix.obj;
%cd(fid_path);
save_dir = GetFullPath(cat(2, mainpath, '/../../../Data/Results/', Sub_ID, '/'));
if ~exist(save_dir, 'dir')
   mkdir(save_dir); 
end
save(cat(2, save_dir, save_results), 'save_results', 'lambda', 'rbins', 'cbins', 'fid_file', 'Gr', 'L', 'Nx', 'Ny', 'Nz', 'Phi', 'U', 'U_init', 'SGblock',...
    'dispim', 'vec', 'sizes', 'ScanType', 'Phi_rt_small_init', 'Phi_rt_full', 'Phi_rt_full_init', 'Phi', 'params', 'Ridx', 'Hidx', 'total_time', 'USR', ...
    'temp_Phi_rt_small', 'Ridx_init', 'Hidx_init', 'Norig', 'ccL', 'cL', 'RR_int', 'Phiresp', 'Phicard', 'L_init');

% Saving binning results
specialNotes = cat(2, 'BinningResults');

if ~isempty(specialNotes)
    specialNotes = cat(2, '_', specialNotes);
end
save_results2 = [fid_file(15:22), '_', Sub_ID, '_', num2str(round(TPR)), 'mm_USR', num2str(round(USR)), ...
    '%_L', num2str(L), '_', 'results_', datestr(now, 'yyyy_mm_dd_HH_MM'), specialNotes, '.mat'];

if num_eco == 1 && f == 1
    save(cat(2, save_dir, save_results2), 'L','Nx', 'Ny', 'Nz','U_init','dispim', 'vec','sizes','ScanType','Phi_rt_full_init','params','Ridx','Hidx',...
        'Ridx_init', 'Hidx_init', 'Norig', 'ccL', 'cL', 'RR_int', 'Phiresp', 'Phicard', 'L_init');
end
% save('20P48_4wk_6mm_2meas_0125_modifiedBinning.mat','-v7.3');
disp('Results saved.')
end

end