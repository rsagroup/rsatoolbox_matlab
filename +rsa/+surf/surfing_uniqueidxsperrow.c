#include "mex.h"
#include <string.h>

/*
 * Input:  matrix X (NxP) for N rows and at most P indices per row
 * Output: matrix Y (NxP) where each row contains the unique, positive 
 * indices from the corresponding row in X; if a value V occurs K times
 * in a row, then it is maintained once, and replaced by zero (K-1) times.
 *
 * NNO May 2010
 */
void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
     int m, n, colcount, i, j, k, targj, maxj;
     double *d, *unq, *unqsmall;
     double v;
     mxArray *mxunqsmall; /* output */
     char msg[80];

     if (nlhs == 0) {
         /* No output argument, so return immediately */
         return;
     }
     if (nlhs != 1) {
         sprintf(msg, "Expected one output argument, but found %i ",nlhs);
         mexErrMsgTxt(msg);
     }

     if (nrhs != 1) {
         sprintf(msg, "Expected one input argument, but found %i",nrhs);
         mexErrMsgTxt(msg);
     }
     
     
     if (!(mxIsDouble(prhs[0]))) {
         mexErrMsgTxt("Only input of type double is supported; use double(X) to convert X to double");
     }


     /* size of input argument*/
     m = mxGetM(prhs[0]);
     n = mxGetN(prhs[0]);
     
     /* pointer to input matrix*/
     d = mxGetPr(prhs[0]);

     /* initial "big" matrix (may have columns with zeros) with unique elements and zeros*/
     unq=calloc(m*n,sizeof(double));
     
     maxj=0; /* number of columns needed in output */
     for(i=0; i<m; i++) { /*row index*/
        targj=0; /* target column for the next unique element in the current row*/
        for (j=0; j<n; j++) { /*input matrix row index*/
            v=d[j*m+i]; /* value in m(i,j)*/
            if (v==0) { /* skip zero*/
                continue;
            }
            
            if (j==n-1) { /*last element, always add*/
                unq[targj++*m+i]=v;
                if (targj>maxj) {
                    maxj=targj;
                }   
            } else {
                for (k=j+1; k<n; k++) { /*output matrix row index*/
                    if (v==d[k*m+i]) { /*value occurs later in this row, do not add*/
                        break;
                    }
                    if (k==n-1) { /* reached the last column, good: add this element */
                        unq[targj++*m+i]=v;
                        if (targj>maxj) {
                            maxj=targj;
                        }
                    }
                }
            }
        }
     }
     
     plhs[0]=mxCreateDoubleMatrix(m, maxj, mxREAL);
     unqsmall=mxGetPr(plhs[0]);
     /* as order is row-column, it's safe to copy just the first part of the matrix */
     memcpy(unqsmall,unq,m*maxj*sizeof(double)); 
     free(unq);
}
