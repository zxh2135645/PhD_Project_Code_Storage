/*
 * MATLAB to CPU interface to call finufft 2D transform type 1 (NU->U) with float precision. Vectorized nufft call possible
 *
 * Author: Junzhou Chen (junzhou@chen.engineer)
 *
 */

#include "mex.h"
//#include "gpu/mxGPUArray.h"
#include <finufft.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>


/* s
 * Host code
 */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, mxArray const *prhs[])
{
    /* Declare all variables.*/
    float *x;          // om(1)
    float *y;          // om(2)
    float _Complex *k; // Pointer to complex kspace-strengths, on host
    mxComplexSingle  *im; // Pointer to complex image-domain data, on host
    // cuFloatComplex *d_c;
    int isSign;
    float eps;
    int64_t Nx;
    int64_t Ny;
    size_t M;      // number of NU points
    size_t Ntrans; // Number of stacked transforms

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
    
    isSign = (int)mxGetScalar(prhs[3]);
    eps = (float)mxGetScalar(prhs[4]);
    Nx = (int64_t)mxGetScalar(prhs[5]);
    Ny = (int64_t)mxGetScalar(prhs[6]);
    M = (size_t)mxGetScalar(prhs[7]);
    Ntrans = (size_t)mxGetScalar(prhs[8]);
    im = (mxComplexSingle*)mxMalloc(Ny * Nx * Ntrans *  sizeof(float _Complex));

    k = (float _Complex *)mxGetComplexSingles(prhs[2]);
    //mexPrintf("Data Assignment Complete\n");


    // Create modes
    int64_t modes[2] = {Nx, Ny};

    //run CPU Vectorized finufftf2d1
    finufftf2d1many(Ntrans, M, x, y, k, isSign, eps, Nx, Ny, im, NULL);

    //Create output 
    mwSize out_dimensions[3] = {(size_t)Nx, (size_t)Ny, Ntrans};
    const mwSize * ndims = out_dimensions;
    plhs[0] = mxCreateNumericArray(3, ndims, mxSINGLE_CLASS, mxCOMPLEX);

    //mexPrintf("Create output  Complete\n");

    //Set Array Value
    mxSetComplexSingles(plhs[0], im);

    return;


}
