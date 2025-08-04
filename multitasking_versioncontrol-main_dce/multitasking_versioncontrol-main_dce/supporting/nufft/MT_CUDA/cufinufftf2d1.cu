/*
 * MATLAB to CUDA interface to call cufinufft 2D transform type 1 (NU->U) with float precision. 
 * Vectorized nufft call possible. This Function assumes all major data variables to be in single (float) 
 * precision and on CPU at input.
 *
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
#include <time.h>
#include <string.h>

/* 
 * Host code
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, mxArray const *prhs[])
{
    clock_t start, end;
    double cpu_time_used;

    /* Declare all variables.*/
    float *x;          // om(1)
    float *y;          // om(2)
    float _Complex *k; // Pointer to complex kspace-strengths, on host
    mxComplexSingle * im ;// Pointer to complex image-domain data, on host
    int isSign;
    float eps;
    int64_t Nx;
    int64_t Ny;
    size_t M;      // number of NU points
    size_t Ntrans; // Number of stacked transforms
    int GPUDeviceIndex; //Which GPU to use

    char const *const errId = "parallel:gpu:mexGPUExample:InvalidInput";
    char const *const errMsg = "Invalid input to MEX file.";
    char const *const instruction = "Usage:\n im = cufinufft2d1(x,y,k,isign,eps,ms,mt, M, Ntrans).";

    

    char errorMessage[1024]; // remember to ensure this buffer is large enough for your message
    // Check input number
    if (nrhs != 9)
    {
        snprintf(errorMessage, sizeof(errorMessage), "Exactly 9 inputs required. \n%s", instruction);
        mexErrMsgIdAndTxt("cufinufftf2d1:InputError", errorMessage);
    }

    // Check input type
    if (!mxIsSingle(prhs[0]) || !mxIsSingle(prhs[1]) || !mxIsSingle(prhs[2]))
    {   
        snprintf(errorMessage, sizeof(errorMessage), "Inputs 1-3 must be single. \n%s", instruction);
        mexErrMsgIdAndTxt("cufinufftf2d1:InputError", errorMessage);
    }

    // Assign Pointers and value
    x = mxGetSingles(prhs[0]);
    y = mxGetSingles(prhs[1]);
    k = (float _Complex *)mxGetComplexSingles(prhs[2]);
    
    isSign = (int)mxGetScalar(prhs[3]);
    eps = (float)mxGetScalar(prhs[4]);
    Nx = (int64_t)mxGetScalar(prhs[5]);
    Ny = (int64_t)mxGetScalar(prhs[6]);
    M = (size_t)mxGetScalar(prhs[7]);
    Ntrans = (size_t)mxGetScalar(prhs[8]);
    //im = (float _Complex *) malloc((size_t)Nx * (size_t)Ny * Ntrans * sizeof(float _Complex));
    im = (mxComplexSingle*)mxMalloc(Ny * Nx * Ntrans *  sizeof(float _Complex));
    //GPUDeviceIndex = (int)mxGetScalar(prhs[9]) - 1; // k is zero based
    //mexPrintf("Data Assignment Complete\n");

    // Initialize the MathWorks GPU API.
    mxInitGPU();
    //cudaSetDevice(GPUDeviceIndex); //Set which GPU to run. This should strictly matach the one being used by MATLAB

    // Create modes
    int64_t modes[2] = {Nx, Ny};

    // Create default cufinufft plan object
    cufinufftf_plan plan;

    // Create on device variables
    float *d_x, *d_y;
    cuFloatComplex *d_k, *d_im;

    // Allocate device memories
    cudaMalloc(&d_x, M * sizeof(float));
    cudaMalloc(&d_y, M * sizeof(float));
    cudaMalloc(&d_k, M * Ntrans * sizeof(float _Complex));
    cudaMalloc(&d_im, (size_t)Nx * (size_t)Ny * Ntrans * sizeof(float _Complex));

    //mexPrintf("cuda Malloc Complete\n");


    // Copy values to device
    start = clock();
    cudaMemcpy(d_x, x, M * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_y, y, M * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_k, k, M * Ntrans * sizeof(float _Complex), cudaMemcpyHostToDevice);
    end = clock();

    //double timeHostToDevice = ((double)(end-start))/ CLOCKS_PER_SEC;
    //mexPrintf("Host to Device Copying Time (s): %.5f\n", timeHostToDevice);

    //mexPrintf("CUDA MEM COPY Complete\n");

    // make cufinufft plan
    /*
     * int cufinufftf_makeplan(int type, int dim, int64_t* nmodes, int iflag, int ntr, float tol, cufinufftf_plan *plan, cufinufft_opts *opts)
     */
    cufinufftf_makeplan(1, 2, modes, isSign, Ntrans, eps, &plan, NULL);

    //mexPrintf("cufinufftf_makeplan Complete\n");


    // Set points
    /*
     * int cufinufftf_setpts(cufinufftf_plan plan, int M, float* x, float* y,loat* z, int N, float* s, float* t, float *u)
     */
    cufinufftf_setpts(plan, M, d_x, d_y, NULL, 0, NULL, NULL, NULL);

    //mexPrintf("cufinufftf_setpts Complete\n");


    /************ EXECUTE ************/
    start = clock();
    cufinufftf_execute(plan, d_k, d_im);
    end = clock();
    //double timeExecute = ((double)(end-start))/ CLOCKS_PER_SEC;
    //mexPrintf("cufinufftf execution Time (s): %.5f\n", timeExecute);
    //mexPrintf("cufinufftf_execute Complete\n");
    /************ EXECUTE END ************/


    // // transfer the data back onto the host, destroy the plan, and free the device arrays.
    start = clock();
    cudaMemcpy(im, d_im, (size_t)Nx * (size_t)Ny * Ntrans * sizeof(float _Complex), cudaMemcpyDeviceToHost);
    end = clock();
    //double timeDeviceToHost = ((double)(end-start))/ CLOCKS_PER_SEC;
    //mexPrintf("Device to Host Copying Time (s): %.5f\n", timeDeviceToHost);

    cufinufftf_destroy(plan);
    cudaFree(d_x);
    cudaFree(d_y);
    cudaFree(d_k);
    cudaFree(d_im);

    //mexPrintf("cufinufftf_cleancup Complete\n");


    //Create output 
    mwSize out_dimensions[3] = {(size_t)Nx, (size_t)Ny, Ntrans};
    const mwSize * ndims = out_dimensions;
    plhs[0] = mxCreateNumericArray(3, ndims, mxSINGLE_CLASS, mxCOMPLEX);

    //mexPrintf("Create output  Complete\n");

    //Set Array Value
    mxSetComplexSingles(plhs[0], im);


    return;


}
