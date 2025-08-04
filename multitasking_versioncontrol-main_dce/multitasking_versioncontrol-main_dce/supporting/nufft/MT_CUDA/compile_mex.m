%% Set The following paths according to your server environments

% check your local machine cuda path
cuda_lib_path = '/usr/local/cuda/lib64';

% run
% > nvidia-smi --query-gpu=compute_cap --format=csv,noheader
% to find GPUs compute capability
% or check at https://arnon.dk/matching-sm-architectures-arch-and-gencode-for-various-nvidia-cards/
gpu_arch = '75';   

% finufft paths
finufft_include = [mainpath '/supporting/nufft/finufft/include'];
finufft_build   = [mainpath '/supporting/nufft/finufft/build'];


%% Realease compilation

mexcuda('-v','AhA_ps_cufinufft.cu',['NVCCFLAGS=-gencode=arch=compute_' gpu_arch ',code=sm_' gpu_arch],'-R2018a',['-I' finufft_include],['-L' finufft_build],['-L' cuda_lib_path],'-lcufinufft','-lcufft')
mexcuda('-v','cufinufftf2d1.cu',['NVCCFLAGS=-gencode=arch=compute_' gpu_arch ',code=sm_' gpu_arch],'-R2018a',['-I' finufft_include],['-L' finufft_build],['-L' cuda_lib_path],'-lcufinufft','-lcufft')

%% Debug compilation

% mexcuda('-v','-G','AhA_ps_cufinufft.cu',['NVCCFLAGS=-gencode=arch=compute_' gpu_arch ',code=sm_' gpu_arch],'-R2018a',['-I' finufft_include],['-L' finufft_build],['-L' cuda_lib_path],'-lcufinufft','-lcufft')
% mexcuda('-v','-G','cufinufftf2d1.cu',['NVCCFLAGS=-gencode=arch=compute_' gpu_arch ',code=sm_' gpu_arch],'-R2018a',['-I' finufft_include],['-L' finufft_build],['-L' cuda_lib_path],'-lcufinufft','-lcufft')


