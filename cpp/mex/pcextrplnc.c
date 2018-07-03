#include "mex.h"

#include "inc_test.h"

/* The gateway function that replaces the "main".
 * plhs[] - the array of output values
 * prhs[] - the array of input values
 */

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  mexPrintf("N: %i %i\n", nlhs, nrhs);

  inc_test();

  //double *par, *val;
  //par = mxGetPr(prhs[0]);
  //int a = mxGetScalar(prhs[0]);
  //val = mxGetPr(plhs[0]);

  mwSize ndim = 2;
  mwSize dims[ndim];
  dims[0] = 2;
  dims[1] = 3;
  //plhs[0] = mxCreateDoubleScalar(2.3);
  plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, false);
}
