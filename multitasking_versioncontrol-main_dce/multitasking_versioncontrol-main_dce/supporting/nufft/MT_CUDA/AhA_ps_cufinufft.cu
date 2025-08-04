/*
 * MATLAB to CUDA interface to run AhA for MR-Multitasking using cufinufft with float precision.
 * Vectorized nufft call possible. This Function assumes all major data variables to be in single (float)
 * precision and on CPU at input. All multidimensional arrays are assumed to have column-major order 
 * (https://en.wikipedia.org/wiki/Row-_and_column-major_order), same as MATLAB.
 * Author: Junzhou Chen (junzhou@chen.engineer)
 *
 */

#include "mex.h"
#include "gpu/mxGPUArray.h"
#include <cufinufft.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <time.h>
#include <string.h>

/**
 * @brief Multiplies image data by sensitivity encoding maps element-wise.
 *
 * This CUDA kernel is designed to handle the multiplication of complex-valued 3D (plus rank) image data with
 * sensitivity encoding maps. The multiplication is broadcasted and performed element-wise and the 
 * results are stored in the output array d_im_SEs.
 *
 * @param d_im Pointer to the input device array containing the image data. The array is a flattened
 *             4D array with dimensions [Nx, Ny, Nz, L], where L is the number of image sets or time points.
 * @param d_SEs Pointer to the input device array containing the sensitivity encoding maps. The array
 *              is a flattened 4D array with dimensions [Nx, Ny, Nz, Ncoils, 1].
 * @param d_im_SEs Pointer to the output device array where the result of the multiplication is stored.
 *                 Dimensions [Nx, Ny, Nz, Ncoils, L].
 * @param Nx The size of the first dimension (width) of the image and sensitivity maps.
 * @param Ny The size of the second dimension (height) of the image and sensitivity maps.
 * @param Nz The size of the third dimension (depth) of the image and sensitivity maps.
 * @param Ncoils The number of coils, corresponding to the fourth dimension in `d_SEs` and `d_im_SEs`.
 * @param L The number of ranks.
 *
 * @note The kernel assumes that all input and output arrays are properly allocated and reside in GPU memory.
 *       It is designed to be executed on a CUDA-capable GPU with a grid size that covers the entire volume
 *       of the image (Nx, Ny, Nz).
 */
__global__ void SEs_times_Im( cuFloatComplex *d_im,  cuFloatComplex *d_SEs,
                             cuFloatComplex *d_im_SEs, int Nx, int Ny, int Nz, int Ncoils, int L)
{
    // Calculate global thread coordinates
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;

    if (x < Nx && y < Ny && z < Nz)
    {
        for (int l = 0; l < L; l++)
        {
            for (int coil = 0; coil < Ncoils; coil++)
            {
                // Calculate the index for d_im (without coil dimension)
                int im_index = ((l * Nz + z) * Ny + y) * Nx + x;

                // Calculate the index for d_SEs(without rank dimension)
                int ses_index = ((coil * Nz + z) * Ny + y) * Nx + x;

                // Calculate the index for d_im_SEs
                int im_ses_index = (((l * Ncoils + coil) * Nz + z) * Ny + y) * Nx + x;

                // Perform the multiplication
                d_im_SEs[im_ses_index] = cuCmulf(d_im[im_index], d_SEs[ses_index]);
            }
        }
    }
}


/**
 * @brief Multiplies multi-coil image data by the complex conjugate of sensitivity encoding maps element-wise.
 * Then sum over the coil dimension.
 *
 * This CUDA kernel is designed to handle the multiplication of complex-valued 3D + coil + rank images data with
 * the conjugate of sensitivity encoding maps across multiple coils. The multiplication is broadcasted performed
 * element-wise. Then the coil dimension is summed and the results are stored in the output array d_im.
 *
 * @param d_im Pointer to the input device array containing the image data. The array is a flattened
 *             4D array with dimensions [Nx, Ny, Nz, L], where L is the number of image sets or time points.
 * @param d_SEs Pointer to the input device array containing the sensitivity encoding maps. The array
 *              is a flattened 4D array with dimensions [Nx, Ny, Nz, Ncoils, 1].
 * @param d_im_SEs Pointer to the output device array where the result of the multiplication is stored.
 *                 Dimensions [Nx, Ny, Nz, Ncoils, L].
 * @param Nx The size of the first dimension (width) of the image and sensitivity maps.
 * @param Ny The size of the second dimension (height) of the image and sensitivity maps.
 * @param Nz The size of the third dimension (depth) of the image and sensitivity maps.
 * @param Ncoils The number of coils, corresponding to the fourth dimension in `d_SEs` and `d_im_SEs`.
 * @param L The number of ranks.
 *
 * @note The kernel assumes that all input and output arrays are properly allocated and reside in GPU memory.
 *       It is designed to be executed on a CUDA-capable GPU with a grid size that covers the entire volume
 *       of the image (Nx, Ny, Nz).
 */
__global__ void SEs_times_Im_adj( cuFloatComplex *d_im,  cuFloatComplex *d_SEs,
                             cuFloatComplex *d_im_SEs, int Nx, int Ny, int Nz, int Ncoils, int L)
{
    // Calculate global thread coordinates
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;

    if (x < Nx && y < Ny && z < Nz)
    {
        for (int l = 0; l < L; l++)
        {
            // Calculate the index for d_im (without coil dimension)
            int im_index = ((l * Nz + z) * Ny + y) * Nx + x;

            cuFloatComplex tempValue = {0.0f, 0.0f};

            for (int coil = 0; coil < Ncoils; coil++)
            {
                
                // Calculate the index for d_SEs(without rank dimension)
                int ses_index = ((coil * Nz + z) * Ny + y) * Nx + x;

                // Calculate the index for d_im_SEs
                int im_ses_index = (((l * Ncoils + coil) * Nz + z) * Ny + y) * Nx + x;

                // Perform the multiplication and summation
                tempValue = cuCaddf(tempValue, cuCmulf(d_im_SEs[im_ses_index], cuConjf(d_SEs[ses_index])));
            }

            d_im[im_index] = tempValue;
            
        }
    }
}

/**
 * @brief CUDA kernel to permute the third dimension to the first dimensnion of a 5D multidimensional 
 * complex data arrays in GPU memory.
 * 
 * This kernel rearranges the dimensions of the input data, specifically targeting
 * a use case involving a 5D array with dimensions corresponding to different parameters
 * such as spatial coordinates (x, y, z), the number of coils (Ncoils), and the number of rank L. 
 * The kernel effectively reorders these dimensions, moving the third dimension
 * (z) to become the first dimension in the output array while preserving the order of the
 * other dimensions. This function exists because we want to make 1D FFT transform in the Z dimension
 * and it is more convineient to call cuFFT when successive signals are arranged next to each other
 * in a pointer array.
 * 
 * @param input Pointer to the input array of complex numbers (float precision)
 *   in GPU memory. The input array is expected to be in a 5D layout flattened into a 1D array,
 *   with dimensions ordered as [Nx, Ny, Nz, Ncoils, L].
 * @param output Pointer to the output array of complex numbers (float precision)
 *   in GPU memory. This array will hold the permuted data, with dimensions reordered to
 *   [Nz, Nx, Ny, Ncoils, L], effectively moving the third dimension (z) of the input to the
 *   first dimension in the output.
 * @param Nx The dimensions of the spatial coordinates (x) of the input data.
 * @param Ny The dimensions of the spatial coordinates (y) of the input data.
 * @param Nz The dimensions of the spatial coordinates (z) of the input data.
 * @param Ncoils The number of coils, representing an additional dimension in the input data.
 * @param L The number of rank in the input data, which, along with Ncoils, succeeds
 *   the spatial dimensions in the input data ordering.
 * 
 */
__global__ void PermuteKernelThirdToFirst(cuFloatComplex *input, cuFloatComplex *output,
                                          int Nx, int Ny, int Nz, int Ncoils, int L)
{
    // Calculate global thread coordinates
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;

    if (x < Nx && y < Ny && z < Nz)
    {
        for (int l = 0; l < L; l++)
        {
            for (int coil = 0; coil < Ncoils; coil++)
            {
                
                // Calculate the index for d_im_SEs
                // int im_ses_index = (((z * Ny + y) * Nx + x) * Ncoils + coil) * L + l;
                int in_index = (((l * Ncoils + coil) * Nz + z) * Ny + y) * Nx + x;

                int out_index = (((l * Ncoils + coil) * Ny + y) * Nx + x) * Nz + z;

                // Perform index exchange
                output[out_index] = input[in_index];
            }
        }
    }
}

/**
 * @brief CUDA kernel to permute the first dimension back to the thrid dimensnion of a 5D multidimensional 
 * complex data arrays in GPU memory.
 * 
 * This kernel rearranges the dimensions of the input data, specifically targeting
 * a use case involving a 5D array with dimensions corresponding to different parameters
 * such as spatial coordinates (x, y, z), the number of coils (Ncoils), and the number of rank L. 
 * The kernel effectively reorders these dimensions, moving the previosuly permuted first dimension
 * (z) to become the third dimension in the output array while preserving the order of the
 * other dimensions. This function exists because we want to make 1D FFT transform in the Z dimension
 * and it is more convineient to call cuFFT when successive signals are arranged next to each other
 * in a pointer array. When the 1D FFT is done, we need to permute back the Z dimension to the original
 * position
 * 
 * @param input Pointer to the input array of complex numbers (float precision)
 *   in GPU memory. The input array is expected to be in a 5D layout flattened into a 1D array,
 *   with dimensions ordered as [Nz, Nx, Ny, Ncoils, L].
 * @param output Pointer to the output array of complex numbers (float precision)
 *   in GPU memory. This array will hold the permuted data, with dimensions reordered to
 *   [Nx, Ny, Ny, Ncoils, L], effectively moving the third dimension (z) of the input to the
 *   first dimension in the output.
 * @param Nx The dimensions of the spatial coordinates (x) of the input data.
 * @param Ny The dimensions of the spatial coordinates (y) of the input data.
 * @param Nz The dimensions of the spatial coordinates (z) of the input data.
 * @param Ncoils The number of coils, representing an additional dimension in the input data.
 * @param L The number of rank in the input data, which, along with Ncoils, succeeds
 *   the spatial dimensions in the input data ordering.
 * 
 */
__global__ void PermuteKernelFirstToThird(cuFloatComplex *input, cuFloatComplex *output,
                                          int Nx, int Ny, int Nz, int Ncoils, int L)
{
    // Calculate global thread coordinates
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;

    if (x < Nx && y < Ny && z < Nz)
    {
        for (int l = 0; l < L; l++)
        {
            for (int coil = 0; coil < Ncoils; coil++)
            {
                
                // Calculate the index for d_im_SEs
                // int im_ses_index = (((z * Ny + y) * Nx + x) * Ncoils + coil) * L + l;
                int out_index = (((l * Ncoils + coil) * Nz + z) * Ny + y) * Nx + x;

                int in_index = (((l * Ncoils + coil) * Ny + y) * Nx + x) * Nz + z;

                // Perform the multiplication
                output[out_index] = input[in_index];
            }
        }
    }
}

/**
 * @brief Executes a 2D Type-2 NUFFT (Non-Uniform Fast Fourier Transform) operation using cuFINUFFT.
 * 
 * This function sets up and executes a 2D Type-2 NUFFT operation, transforming uniformly
 * sampled spatial data into non-uniformly sampled frequency (kspace) data. It involves creating a plan,
 * setting the non-uniform sample points, executing the transform, and then cleaning up the plan.
 * The function assumes that all input arrays and parameters have been properly initialized and
 * that the cuFINUFFT library has been correctly set up in the environment.
 * 
 * @param d_im_SEs Pointer to the input array of complex numbers (float precision) in GPU memory,
 *                 which will hold the multi-coil image data. [Nx, Ny, Nz, Ncoils, L]. Note the thrid dimension is
 *                 already in kspace
 * @param d_k Pointer to the output kspace array of complex numbers (float precision) in GPU memory,
 *            containing the data in the frequenct domain (kspace). Dimension [M, Ntrans], equivalent to [Ntrajs, Nx, Nz, Ncoils, L].
 *            with Ntrajs being the number of unique radial spokes.
 * @param d_x Pointer to the array of x-coordinates (float) in GPU memory for the non-uniform sample points.
 * `          Contrains M points of range [-pi, pi)
 * @param d_y Pointer to the array of y-coordinates (float) in GPU memory for the non-uniform sample points.
 *            Contrains M points of range [-pi, pi)
 * @param modes Pointer to an array of two int64_t values specifying the number of modes in each dimension
 *              (x and y) for the output frequency domain data.
 * @param M The number of non-uniform sample points. i.e. Ntrajs * Nro(readout length)
 * @param Ntrans The number of stacked transformations to perform. This parameter allows for batched execution
 *               of multiple transforms using the same plan but different data.
 * 
 * @note This function supports stacked 2D transform, allowing multiple 2D transforms across Z, coil and rank
 * dimensions. No explicit for-loop needed.
 */
void cufinufftf2d2(cuFloatComplex *d_im_SEs, cuFloatComplex * d_k, float* d_x, float *d_y, 
                   const int64_t *modes, size_t M, size_t Ntrans)
{
    cufinufftf_plan plan; //declare cufinufft plan object

    //Assuing all values properly initialzed.
    //int cufinufftf_makeplan(int type, int dim, const int64_t *n_modes, int iflag, int ntr, float eps, cufinufftf_plan *d_plan_ptr, cufinufft_opts *opts)
    cufinufftf_makeplan(2, 2, modes, -1, Ntrans, 1e-6, &plan, NULL);// make plan. Deafault options.

    //Set Points
    cufinufftf_setpts(plan, M, d_x, d_y, NULL, 0, NULL, NULL, NULL);

    //Execute
    cufinufftf_execute(plan, d_k, d_im_SEs);

    cufinufftf_destroy(plan);

}

/**
 * @brief Executes a 2D Type-1 NUFFT (Non-Uniform Fast Fourier Transform) operation using cuFINUFFT.
 * 
 * This function sets up and executes a 2D Type-1 NUFFT operation, transforming non-uniformly
 * sampled frequency (kspace) data into uniformly sampled image data. It involves creating a plan,
 * setting the non-uniform sample points, executing the transform, and then cleaning up the plan.
 * The function assumes that all input arrays and parameters have been properly initialized and
 * that the cuFINUFFT library has been correctly set up in the environment.
 * 
 * @param d_im_SEs Pointer to the ouput array of complex numbers (float precision) in GPU memory,
 *                 which will hold the multi-coil image data. [Nx, Ny, Nz, Ncoils, L]. Note the thrid dimension is
 *                 already in kspace. 
 * @param d_k Pointer to the input kspace array of complex numbers (float precision) in GPU memory,
 *            containing the data in the frequenct domain (kspace). Dimension [M, Ntrans], equivalent to [Ntrajs, Nx, Nz, Ncoils, L].
 *            with Ntrajs being the number of unique radial spokes.
 * @param d_x Pointer to the array of x-coordinates (float) in GPU memory for the non-uniform sample points.
 * `          Contrains M points of range [-pi, pi)
 * @param d_y Pointer to the array of y-coordinates (float) in GPU memory for the non-uniform sample points.
 *            Contrains M points of range [-pi, pi)
 * @param modes Pointer to an array of two int64_t values specifying the number of modes in each dimension
 *              (x and y) for the output frequency domain data.
 * @param M The number of non-uniform sample points. The number of non-uniform sample points. i.e. Ntrajs * Nro(readout length)
 * @param Ntrans The number of stacked transformations to perform. This parameter allows for batched execution
 *               of multiple transforms using the same plan but different data.
 * 
 * @note This function supports stacked 2D transform, allowing multiple 2D transforms across Z, coil and rank
 * dimensions. No explicit for-loop needed.
 */
void cufinufftf2d1(cuFloatComplex *d_im_SEs, cuFloatComplex * d_k, float* d_x, float *d_y, 
                   const int64_t *modes, size_t M, size_t Ntrans)
{
    cufinufftf_plan plan; //declare cufinufft plan object

    //Assuing all values properly initialzed.

    //int cufinufftf_makeplan(int type, int dim, const int64_t *n_modes, int iflag, int ntr, float eps, cufinufftf_plan *d_plan_ptr, cufinufft_opts *opts)
    cufinufftf_makeplan(1, 2, modes, 1, Ntrans, 1e-6, &plan, NULL);// make plan, deafault options

    //Set Points
    cufinufftf_setpts(plan, M, d_x, d_y, NULL, 0, NULL, NULL, NULL);

    //Execute
    cufinufftf_execute(plan, d_k, d_im_SEs);

    cufinufftf_destroy(plan);

}


/**
 * @brief Updates each k-space [1 x L] line with its corresponding (� * �H) [L x L]. Overall, this amounts to a block-diagnal matrix 
 * multiplication.
 * 
 * 
 * @param d_k Pointer to the input k-space data. In GPU memory. This is a complex-valued array representing the k-space data. 
 *            Dimension [M, Ntrans], equivalent to [Ntrajs, Nx, Nz, Ncoils, L]. with Ntrajs being the number of unique radial spokes.
 * @param d_Phi2 Pointer to the overall (� * �H) transformation matrix. In GPU memory. This is a complex-valued array representing the transformation to be applied.
 *               Dimension [L, L, Ntraj, Nz]
 * @param temp Pointer to the output buffer where the updated k-space data will be stored. This is a complex-valued array.
 *             In GPU memory and has the same dimension as d_k [M, Ntrans], equivalent to [Ntrajs, Nx, Nz, Ncoils, L]
 * @param Ntrajs The number of trajectories or unique radial spokes in the k-space data.
 * @param Nx The size of the k-space data in the x-dimension.
 * @param Nz The size of the k-space data in the z-dimension.
 * @param Ncoils The number of coils. 
 * @param L The number of ranks.
 * 
 * 
 * @note The kernel assumes that the input data (d_k and d_Phi2) are stored in column-major order.
 * 
 * d_k [Ntrajs, Nx, Nz, Ncoils, L]
 *                              ^innerL
 * 
 * d_Phi2 [L,       L,       Ntraj, Nz]
 *         ^innerL  ^l
 */
__global__ void KspaceUpdate(cuFloatComplex *d_k, cuFloatComplex *d_Phi2, cuFloatComplex* temp, 
                            size_t Ntrajs, size_t Nx, size_t Nz, size_t Ncoils, size_t L)
{
    // Calculate global thread coordinates
    int traj = blockIdx.x * blockDim.x + threadIdx.x; //Ntrajs spoke dimension in d_k
    int x = blockIdx.y * blockDim.y + threadIdx.y; //Nx readout dimension in d_k
    int z = blockIdx.z * blockDim.z + threadIdx.z; //Nz partition dimension in d_k

    if (traj < Ntrajs && x < Nx && z < Nz)
    {
        for (int l = 0; l < L; l++)
        {
            for (int coil = 0; coil < Ncoils; coil++)
            {
                size_t out_idx = (((l * Ncoils + coil) * Nz + z) * Nx + x) * Ntrajs + traj;

                cuFloatComplex tempValue = {0.0f,0.0f};
                //Vector Matrix Multiplications
                for (int innerL = 0; innerL < L; innerL++)
                {
                    size_t k_index = (((innerL * Ncoils + coil) * Nz + z) * Nx + x) * Ntrajs + traj;

                    size_t Phi2_index = ((z * Ntrajs + traj) * L + l) * L + innerL;
                    
                    tempValue  =  cuCaddf(tempValue, cuCmulf(d_k[k_index], d_Phi2[Phi2_index]));
                }

                temp[out_idx] = tempValue;

            }
        }
    }

}


/**
 * @brief Main functioin for the MR-Multitasking subspace recon AhA operations for stack-of-stars trajectory
 * 
 * Overall, it does the following:
 * 1. Initialize necessary parameters, data array and does rudimentray data checking.
 * 2. copy data from host (CPU) to GPU device.
 * 3. Perform AhA calculations.
 * 4. Transfer calculated data from GPU device back to host. Then return to matlab workspace.
 * 
 * All major data arrays are assumed to be in single precision and stored in column-major order on CPU.
 * 
 * @param x prhs[0]. Spatial basis function image input of size [Nx, Ny, Nz, L]. 
 * @param st prhs[1]. Matlab structure containing necessary data parameters incluing in-plane matrix size N, Number of NU points M (Ntrajs * Nro)
 *           and actual trajectoris points om [M, 2] for kx and ky.
 * @param SEs prhs[2]. Sensitivity maps of dimension [Nx, Ny, Nz, Ncoils].
 * @param Phi2 prhs[3]. The overall (� * �H) matrix. Dimension [L, L, Ntraj, Nz].
 * 
 * @param AhAx plhs[0]. The flatenned output of the operation AhA(x). Dimension [(Nx * Ny * Nz * L), 1]
 * 
 * @note AhAx = AhA_ps_cufinufft(x, st, SEs, Phi2)
 *
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, mxArray const *prhs[])
{
    clock_t start, end;
    double cpu_time_used;

    /* Declare all variables.*/
    float *x;              // om(1)
    float *y;              // om(2)
    mxComplexSingle *im;   // Pointer to complex image-domain data (x), on host [Nx, Ny, Nz, L]
    mxComplexSingle *SEs;  // Pointer to complex Complex sensitivity maps (SEs), on host [Nx, Ny, Nz, Ncoils]
    mxComplexSingle *Phi2; // Pointer to complex Phi2 data (Phi2, aka, Om), on host [L, L, Ntrajs, Nz]
    int isSign = 1;        // Whether to use + sign in complex exponential. Default to true
    float eps = 1e-6;      // NUFFT precision, default 1e-6;
    int64_t Nx;
    int64_t Ny;
    int64_t Nz;
    int64_t Ncoils;
    int64_t L;          // rank
    int64_t Nkx;        // length of one k-space line // LeeHL
    size_t M;           // number of NU points
    size_t Ntrans;      // Number of stacked transforms
    size_t Ntrajs;      // Number of in-plane tranjectories. e.g. 280
    //int GPUDeviceIndex; // Which GPU to use

    char const *const errId = "parallel:gpu:mexGPUExample:InvalidInput";
    char const *const errMsg = "Invalid input to MEX file.";
    char const *const instruction = "Usage:\n AhAx = AhA_ps_cufinufft(x, st, SEs, Phi2).";

    char errorMessage[1024]; // remember to ensure this buffer is large enough for your message
    // Check input number
    if (nrhs != 4)
    {
        snprintf(errorMessage, sizeof(errorMessage), "Exactly 4 inputs are required. \n%s", instruction);
        mexErrMsgIdAndTxt("AhA_ps_cufinufft:InputError", errorMessage);
    }

    // Check input type
    if (!mxIsSingle(prhs[0]) || !mxIsSingle(mxGetField(prhs[1], 0, "om")) || !mxIsSingle(prhs[2])|| !mxIsSingle(prhs[3]))
    {
        snprintf(errorMessage, sizeof(errorMessage), "Inputs original image (x), NU trajectories (st.om), Sensitivity maps (SEs) and Phi2 must be single. \n%s", instruction);
        mexErrMsgIdAndTxt("AhA_ps_cufinufft:InputError", errorMessage);
    }

    // ==========Assign Pointers and values ========== //

    // Input 1: image domain data
    im = mxGetComplexSingles(prhs[0]); // Get input image domain data



    // Input 2: Trajectories and Image Dimensions
    // Separate om(:,1) and om(:,2)
    size_t numRows = mxGetM(mxGetField(prhs[1], 0, "om"));
    float *om = mxGetSingles(mxGetField(prhs[1], 0, "om"));
    x = (float *)malloc(numRows * sizeof(float));
    y = (float *)malloc(numRows * sizeof(float));
    for (size_t i = 0; i < numRows; i++)
    {
        x[i] = om[i];
        y[i] = om[numRows + i];
    }

    Nx = (int64_t)mxGetDoubles(mxGetField(prhs[1], 0, "Nd"))[0];
    Ny = (int64_t)mxGetDoubles(mxGetField(prhs[1], 0, "Nd"))[1];
    Nz = (int64_t)mxGetScalar(mxGetField(prhs[1], 0, "Nz"));
    M = (size_t)mxGetScalar(mxGetField(prhs[1], 0, "M"));
    //Ntrajs = (size_t)(M/Nx);  // LeeHL
    

    // Input 3: Sensitivity Maps
    const mwSize *coilDims = mxGetDimensions(prhs[2]);
    Ncoils = (int64_t)coilDims[3];
    SEs = mxGetComplexSingles(prhs[2]);

    // Input 4: Phi2 (aka, Om)
    const mwSize *phi2Dims = mxGetDimensions(prhs[3]);
    L = (int64_t)phi2Dims[0];
    Phi2 = mxGetComplexSingles(prhs[3]);

    // Get k-space size from Phi2 and st    // LeeHL
    Ntrajs = (size_t)phi2Dims[2];
    Nkx = (int64_t) (M/Ntrajs);

    Ntrans = Nz * Ncoils * L;

    // ==========Initialize the MathWorks GPU API. ========== //
    mxInitGPU();

    // ==========Broadcast SEs to image input on GPU. (SEs * im)========== //
    cuFloatComplex *d_im, *d_SEs; // coil-combined image and SEs on device
    cuFloatComplex *d_im_SEs;     // im * SEs on devices

    // Allocate device memories
    cudaMalloc(&d_im, Nx * Ny * Nz * L * sizeof(cuFloatComplex));
    cudaMalloc(&d_SEs, Nx * Ny * Nz * Ncoils * sizeof(cuFloatComplex));
    cudaMalloc(&d_im_SEs, Nx * Ny * Nz * Ncoils * L * sizeof(cuFloatComplex));

    // Copy data to device
    cudaMemcpy(d_im, im, Nx * Ny * Nz * L * sizeof(cuFloatComplex), cudaMemcpyHostToDevice);
    cudaMemcpy(d_SEs, SEs, Nx * Ny * Nz * Ncoils * sizeof(cuFloatComplex), cudaMemcpyHostToDevice);

    // Intialze Blocks and Grids and Perform Broadcasting
    dim3 blockSize(16, 16, 1);
    dim3 gridSize((Nx + 15) / 16, (Ny + 15) / 16, Nz);

    SEs_times_Im<<<gridSize, blockSize>>>(d_im, d_SEs, d_im_SEs, Nx, Ny, Nz, Ncoils, L);
    //cudaDeviceSynchronize();

    cudaFree(d_im); // Free subsequently unsed variable d_im on device


    // ========== Forward cuFFT in the 3rd (Z) dimension ========== /

    // Permute 3rd dimension to first dimension
    cuFloatComplex *temp;
    cudaMalloc(&temp, Nx * Ny * Nz * Ncoils * L * sizeof(cuFloatComplex));
    PermuteKernelThirdToFirst<<<gridSize, blockSize>>>(d_im_SEs, temp, Nx, Ny, Nz, Ncoils, L);

    cudaFree(d_im_SEs);
    d_im_SEs = temp; //Assign image to permuted memoiry
    temp = nullptr; //Prevent accidental use of temp

    //Make 1D cuFFT Plan
    cufftHandle fft_plan; //Decalre plan object
    cufftPlan1d(&fft_plan, Nz, CUFFT_C2C, Nx * Ny * Ncoils * L);

    //Execute
    cufftExecC2C(fft_plan, d_im_SEs, d_im_SEs, CUFFT_FORWARD);

    // Permute 1st dimension back to to 3rd dimension
    cudaMalloc(&temp, Nx * Ny * Nz * Ncoils * L * sizeof(cuFloatComplex));
    PermuteKernelFirstToThird<<<gridSize, blockSize>>>(d_im_SEs, temp, Nx, Ny, Nz, Ncoils, L);

    cudaFree(d_im_SEs);
    d_im_SEs = temp; //Assign image to permuted memoiry
    temp = nullptr; //Prevent accidental use of temp

    //Destroy cuFFT Plan
    cufftDestroy(fft_plan);


    // ========== FINUFFT Type 2, U -> NU========== /

    //Declare device variables
    float *d_x, *d_y;
    cuFloatComplex *d_k, *d_Phi2;

    int64_t modes[2] = {Nx, Ny};

     // Allocate device memories
    cudaMalloc(&d_x, M * sizeof(float));
    cudaMalloc(&d_y, M * sizeof(float));
    cudaMalloc(&d_k, M * Ntrans * sizeof(float _Complex));
    cudaMalloc(&d_Phi2, Ntrajs * Nz * L * L * sizeof(float _Complex));

    // Copy values to device
    cudaMemcpy(d_x, x, M * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_y, y, M * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_Phi2, Phi2, Ntrajs * Nz * L * L * sizeof(float _Complex), cudaMemcpyHostToDevice);


    //Do Type 2 transform
    cufinufftf2d2(d_im_SEs, d_k, d_x, d_y, modes, M, Ntrans);

    //values in d_im_SEs no longer used, clear
    cudaFree(d_im_SEs);
    d_im_SEs = nullptr; //Prevent accidental use



    // =========Block-Diagonal kspace update ========/
    cudaMalloc(&temp, Ntrajs * Nkx * Nz * Ncoils * L * sizeof(cuFloatComplex));
    //cudaMemset(temp, 0, Ntrajs * Nx * Nz * Ncoils * L * sizeof(cuFloatComplex));

    dim3 blockSizeBlkdiag(16, 16, 1);
    dim3 gridSizeBlkdiag((Ntrajs + 15) / 16, (Nx + 15) / 16, Nz);
    KspaceUpdate<<<gridSizeBlkdiag, blockSizeBlkdiag>>>(d_k, d_Phi2, temp, Ntrajs, Nkx, Nz, Ncoils, L);

    cudaFree(d_k);
    d_k = temp; //Assign ksspace to updated values memory
    temp = nullptr; //Prevent accidental use of temp




    // ========= FINUFFT Type 1, NU -> U ========/
    cudaMalloc(&d_im_SEs, Nx * Ny * Nz * Ncoils * L * sizeof(cuFloatComplex));// Reallocate device multi-coil image data mem
    cufinufftf2d1(d_im_SEs, d_k, d_x, d_y, modes, M, Ntrans);

    //Clear out device memories
    cudaFree(d_k);
    cudaFree(d_x);
    cudaFree(d_y);
    cudaFree(d_Phi2);




    // ========== Inverse cuFFT in the 3rd (Z) dimension ========== /

    // Permute 3rd dimension to first dimension
    cudaMalloc(&temp, Nx * Ny * Nz * Ncoils * L * sizeof(cuFloatComplex));
    PermuteKernelThirdToFirst<<<gridSize, blockSize>>>(d_im_SEs, temp, Nx, Ny, Nz, Ncoils, L);

    cudaFree(d_im_SEs);
    d_im_SEs = temp; //Assign image to permuted memoiry
    temp = nullptr; //Prevent accidental use of temp

    //Make 1D cuFFT Plan
    cufftHandle ifft_plan; //Decalre plan object
    cufftPlan1d(&ifft_plan, Nz, CUFFT_C2C, Nx * Ny * Ncoils * L);

    //Execute
    cufftExecC2C(ifft_plan, d_im_SEs, d_im_SEs, CUFFT_INVERSE);

    // Permute 1st dimension back to to 3rd dimension
    cudaMalloc(&temp, Nx * Ny * Nz * Ncoils * L * sizeof(cuFloatComplex));
    PermuteKernelFirstToThird<<<gridSize, blockSize>>>(d_im_SEs, temp, Nx, Ny, Nz, Ncoils, L);

    cudaFree(d_im_SEs);
    d_im_SEs = temp; //Assign image to permuted memoiry
    temp = nullptr; //Prevent accidental use of temp

    //Destroy cuFFT Plan
    cufftDestroy(ifft_plan);





    // ========= coil_combination ========/
    cudaMalloc(&d_im, Nx * Ny * Nz * L * sizeof(cuFloatComplex));// Reallocate device single-coil image data mem
    SEs_times_Im_adj<<<gridSize, blockSize>>>(d_im, d_SEs, d_im_SEs, Nx, Ny, Nz, Ncoils, L);



    
    // ========= Post-operation clean-up========/
    //clear out device memories
    cudaFree(d_SEs);
    cudaFree(d_im_SEs);

    // Copy results to host
    mxComplexSingle *AhAx = (mxComplexSingle *)mxMalloc(Ny * Nx * Nz  * L * sizeof(mxComplexSingle));
    cudaMemcpy(AhAx, d_im, Ny * Nx * Nz  * L * sizeof(mxComplexSingle), cudaMemcpyDeviceToHost);

    cudaFree(d_im);
   // cudaFree(d_im_SEs);


    // Create matlab output
    mwSize out_dimensions[1] = {(size_t)Nx * (size_t)Ny * (size_t)Nz  * L};
    const mwSize *ndims = out_dimensions;
    plhs[0] = mxCreateNumericArray(1, ndims, mxSINGLE_CLASS, mxCOMPLEX);

    // mexPrintf("Create output  Complete\n");

    // Set matlab output array Value
    mxSetComplexSingles(plhs[0], AhAx);

    // Free host allocated memory
    free(x);
    free(y);

    mexPrintf(".");
    return;
}