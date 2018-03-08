#include "mex.h"
#include <math.h>

/*
 * Input: matrices XS (3xP) and YS (3xQ)
 * Output: distance matrix DS (Px3), with DS(i,j) the Euclidian distance 
 * between XS(:,i) and YS(:,j
 *
 * NNO May 2010
 */
void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
     int i, j, m0, n0, m1, n1;
     double *d0, *d1, *o0;
     double k0,k1,k2;
     char msg[80];

     if (nlhs == 0) {
         /* No output argument, so return immediately */
         return;
     }
     if (nlhs != 1) {
         sprintf(msg, "Expected one output argument, but found %i ",nlhs);
         mexErrMsgTxt(msg);
     }

     if (nrhs != 2) {
         sprintf(msg, "Expected two input arguments, but found %i",nrhs);
         mexErrMsgTxt(msg);
     }
     
     if (!(mxIsDouble(prhs[0]) && mxIsDouble(prhs[1]))) {
         mexErrMsgTxt("Only input of type double is supported; use double(X) to convert X to double");
     }


     /* size of first input argument*/
     m0 = mxGetM(prhs[0]);
     n0 = mxGetN(prhs[0]);

     /* size of second input argument*/
     m1 = mxGetM(prhs[1]);
     n1 = mxGetN(prhs[1]);

    if (m0 != 3) {
        mexErrMsgTxt("First input should be 3xP.");
    }

    if (m1 != 3) {
        mexErrMsgTxt("Second argument should be 3xQ");
    }

    /* Input data */
    d0 = mxGetPr(prhs[0]);
    d1 = mxGetPr(prhs[1]);

    /* Space for output */
    plhs[0] = mxCreateDoubleMatrix(n0, n1, mxREAL);

    /* Create a pointer to the output data */
    o0 = mxGetPr(plhs[0]);

    for (i=0; i<n0; i++) {
        for (j=0; j<n1; j++) {
            /* difference along x, y, and z direction */
            k0=d0[i*3]-d1[j*3];
            k1=d0[i*3+1]-d1[j*3+1];
            k2=d0[i*3+2]-d1[j*3+2];

            o0[j*n0+i]=sqrt(k0*k0+k1*k1+k2*k2);
        }
    }  
}
