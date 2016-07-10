/*==========================================================
 * %  Mex version to compute C*A*C'
 * %  Where C is the pairmatrix of K*(K-1)/2 x K size
 * %  This consitutes the main time sink in RSA-based computations,
 * %  for example for varianceLDC
 * % INPUT:
 * %   A: KxK matrix
 * %
 * % OUTPUT:
 * %   B : K*(K-1)/2 x K matrix
 * %
 * % (c) Joern Diedrichsen
 *========================================================*/

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    // (1) check for proper number of input and output arguments
    if(nrhs!=1)
        mexErrMsgIdAndTxt("rsaToolbox:pairCACt:oneInput","1 input required");
    
    int K,K1,D;
    int d,i1,j1,i2,j2,indxI,indxJ;
    mxClassID myClassID;
    
    // getting dimensions
    K = mxGetM(prhs[0]);
    K1 = mxGetN(prhs[0]);
    
    if(K!=K1)
        mexErrMsgIdAndTxt("rsaToolbox:pairCACt:squareMatrix","Square matrix required");
    
    // Size of distance
    D = K * (K-1) / 2;
    
    // Check class
    myClassID = mxGetClassID(prhs[0]);
    if(myClassID==mxDOUBLE_CLASS){
        double *ptrA, *ptrOut;
        double tmp;
        plhs[0]= mxCreateNumericMatrix(D, D,myClassID, 0);
        ptrA    = mxGetPr(prhs[0]);
        ptrOut  = mxGetPr(plhs[0]);
        
        d=0;
        for (i2=0;i2<K;i2++) {
            indxI=i2*K;
            for (j2=i2+1;j2<K;j2++){
                indxJ=j2*K;
                for (i1=0;i1<K;i1++) {
                    tmp = ptrA[i1+indxI]-ptrA[i1+indxJ];
                    for (j1=i1+1;j1<K;j1++){
                        ptrOut[d]=tmp - ptrA[j1+indxI]+ptrA[j1+indxJ];
                        d++;
                    }
                }
            }
        }
        return;
    } else if (myClassID==mxSINGLE_CLASS){
        float *ptsA, *ptsOut;
        float tmpS;
        plhs[0]= mxCreateNumericMatrix(D, D,myClassID, 0);
        ptsA    = mxGetData(prhs[0]);
        ptsOut  = mxGetData(plhs[0]);
        
        d=0;
        for (i2=0;i2<K;i2++) {
            indxI=i2*K;
            for (j2=i2+1;j2<K;j2++){
                indxJ=j2*K;
                for (i1=0;i1<K;i1++) {
                    tmpS = ptsA[i1+indxI]-ptsA[i1+indxJ];
                    for (j1=i1+1;j1<K;j1++){
                        ptsOut[d]=tmpS - ptsA[j1+indxI]+ptsA[j1+indxJ];
                        d++;
                    }
                }
            }
        }
        return;
    } else { 
        mexErrMsgIdAndTxt("rsaToolbox:pairCACt:numberFormat","Single or Double matrix required");
    } 
}




