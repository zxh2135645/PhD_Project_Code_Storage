# MT_CUDA
## Description
This repository contains c and cuda-based source codes to perform efficient GPU-based AhA operation for stack-of-stars trajectory in MR-Multitasking reconconstruction. Under the hood, it makes api calls to [cuFFT](https://docs.nvidia.com/cuda/cufft/index.html) and [cufinufft](https://finufft.readthedocs.io/en/latest/c_gpu.html) for both uniform and non-uniform FFTs. This repo also contains self-defined CUDA kernels to perform operations such as coil-sensitivity maps multiplication (broadcast) and block-diagnal low-rank kspace update with $\Phi\Phi^H$. The main function is [AhA_ps_cufinufft.cu](AhA_ps_cufinufft.cu). Additionally, this repo contains example simple nufft API calls to [cufinufftf2d1(GPU)](cufinufftf2d1.cu) and [finufftf2d1(CPU)](finufftf2d1.c).

At the moment, the AhA operation takes around 8 seconds for a image matrix size of $320 \times 320 \times 56$ with $L = 16$ and $Ncoils = 12$.

>Hardware environments:
> * NVIDIA GeForce RTX 3090 (24G)
> * NVIDIA-SMI 535.171.04
> * Driver Version: 535.171.04
> * CUDA Version: 12.2

## Compilation and Usage
>CAUTION: This project is in its early phase of developements. Please exercise patientce when setting up necessary environments for compilation and runtime. The original author has undertaken efforts to reduce the chance of unexpected errors but they still might occur and lead to matlab crashing.

#### Dependencies
1. You will need to have [CUDA toolkit](https://developer.nvidia.com/cuda-toolkit) installed. 
2. You will need to download and compile [finufft](https://github.com/flatironinstitute/finufft) with [CUDA enabled](https://finufft.readthedocs.io/en/latest/install_gpu.html).
3. You will need matlab version R2018a or later.

#### Compilation
1. Set the paths in [compile_mex.m](/compile_mex.m).
2. run [compile_mex.m](/compile_mex.m) in matlab.

>To compile with debug flags, uncomment the debug section.

#### Runtime Environments
Because the mex is compiled with references to compiled dynamic libraries from finufft and cuFFT, it is necessary to let matlab know where to find these libraries. This is achieved by setting the environment variable `LD_LIBRARY_PATH` before launching matlab. Additionally, in a multi-GPU server, we want to handle the GPU selection by setting the encironment variable `CUDA_VISIBLE_DEVICES` before launching matlab since we are making low-level API calls to CUDA codes. 

A shell script [matlab_with_env.sh](/matlab_with_env.sh) is provided to launch matlab after selecting the GPU with the largest available memory and setting the `LD_LIBRARY_PATH` variable.

1. Set the `FINUFFT_LIB` and `CUDA_LIB` variables inside the script.
2. Make sure the script is executable. If not, run in shell:
    ```bash
    chmod +x ./matlab_with_env.sh
    ```
3. Launch matlab using the script
    ```bash
    ./matlab_with_env.sh
    ```
    >Note: Copy the script in one of your `PATH` folders (e.g. `$HOME/.local/bin`) to launch without the relative path indicator `./`

#### Running the function  `AhA_ps_cufinufft`

Run the function in matlab with 
```matlab
AhAx = AhA_ps_cufinufft(x, st, SEs, Phi2);
```

>The function assumes all major data arrays `x`, `st.om`, `SEs`, `Phi2` are stored on CPU and in single (float) precision. For detailed description of these input variables, refer to the method description for `mexFunction` in [AhA_ps_cufinufft.cu](/AhA_ps_cufinufft.cu#L408)

A test matlab script [cu_finufft_test.m](/cu_finufft_test.m) is provided. 

You can download the test data `AhAx_test.mat` from this [link](https://junzhou.chen.engineer/MT_TEST/AhAx_test.mat) or
```bash
wget https://junzhou.chen.engineer/MT_TEST/AhAx_test.mat
```

## Development and Debugging in Visual Studio Code
This code was developed in Visual Studio Code remotely with ssh. Set the `includePath` in [c_cpp_properties.json](/.vscode/c_cpp_properties.json) for IntelliSense and auto code completion.

To launch debug session:
1. Compile mex with debug flags enabled.
2. Change `program` in [launch.json](/.vscode/launch.json) to the matlab executable file in your server. 
3. Make sure `environment` for `LD_LIBRARY_PATH` in [launch.json](/.vscode/launch.json) is set to the specification described above.
4. Set breakpoints in [AhA_ps_cufinufft.cu](/AhA_ps_cufinufft.cu) and launch directly in VS Code.
