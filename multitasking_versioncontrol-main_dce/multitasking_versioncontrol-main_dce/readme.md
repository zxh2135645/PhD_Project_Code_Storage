## MR Multitasking reconstruction code v1.0 - VB and VE multitasking sequence support

It currently covers 2D/3D Cartesian, 2D/2D-SMS/3D radial trajectories, T1/T2/T2* contrast

This scipt is majorly for Human Study.

Xinheng Zhang 11/18/2024
This version is working in progress -> to averaging 680 into 170

## Install nufft packages before first use (if you need nufft, choose one from below)

### finufft (recommended)
* Linux 
 
	If you are NOT using recon3:
	
	Before first use, remove the mainpath/supporting/nufft/finufft folder.
	Open a terminal in mainpath/supporting/nufft/ and run
	
	```
	git clone https://github.com/flatironinstitute/finufft.git
	cd finufft && mkdir build && cd build
	cmake -D FINUFFT_USE_CUDA=ON -D FINUFFT_BUILD_MATLAB=ON ..
	cmake --build . -j
	```
	Make sure the 'matlab' folder is under mainpath/supporting/nufft/finufft and the mex files are inside.
	
	Remove the FINUFFT_USE_CUDA flag if not using GPU. 
	If need to build for a specific compute capability, find GPU info with
	```
	nvidia-smi --query-gpu=compute_cap --format=csv,noheader
	```
	and add CMAKE_CUDA_ARCHITECTURES flag to the cmake command:
	```
	cmake -D FINUFFT_USE_CUDA=ON -D FINUFFT_BUILD_MATLAB=ON -D CMAKE_CUDA_ARCHITECTURES=75 ..
	```

* Windows 10 
 
    Download the [finufft package](https://github.com/flatironinstitute/finufft)
 
    and copy the 'matlab' folder to mainpath\supporting\nufft\finufft\

    Download the [mex file](https://users.flatironinstitute.org/~ahb/codes/finufft-binaries/2.0.2/win/finufft.mexw64) 
	and save it to mainpath\supporting\nufft\finufft\matlab\

    If using Windows 10, copy all dll files under mainpath\supporting\nufft\finufft\winlib
    to C:\Windows\system32\

### Matlab cufinufft interface
* Linux Only

	Written by Junzhou Chen (Junzhou.Chen@med.usc.edu)

	If CUDA path is not already in your bash LD_LIBRARY_PATH, run the following lines: (change CUDA_LIB to the proper CUDA library path on yur machine)
	```
	CUDA_LIB='/usr/local/cuda/lib64'
	echo "export LD_LIBRARY_PATH=$CUDA_LIB:\$LD_LIBRARY_PATH" >> ~/.bashrc
	source ~/.bashrc
	```
	If you are NOT using recon3:
	
	Before first use, open mainpath/supporting/nufft/MT_CUDA/compile_mex.m,	change the cuda settings according to your system environment, then run the following lines in Matlab command window:
	```
	compile_mex
	copyfile *.mexa64 ../finufft/matlab
	```
	If cufinufft keeps failing, try open matlab by running (remember to check the paths in the script)
	```
	mainpath/matlab_cufinufft.sh
	```
	
### gpuNUFFT 
* Linux Only 
 
    Open a terminal in mainpath/supporting/nufft/
    and run
    ```
    git clone https://github.com/andyschwarzl/gpuNUFFT.git
    cd gpuNUFFT/CUDA && mkdir -p build && cd build
    cmake ..
    make
    ```
 
 ### MIRT (optional, slow, only use this if other tools are not available)
    
	After adding the MIRT toolbox path, run in MATLAB command window
    ```
    ir_mex_build_table
    ```
	to build the mex files

### BART (optional)
* Linux Only 

	Open a terminal in mainpath/supporting/nufft/
    and run
	```
    sudo apt-get install make gcc libfftw3-dev liblapacke-dev libpng-dev libopenblas-dev
	git clone https://github.com/mrirecon/bart.git
	cd bart
	make
    ```
